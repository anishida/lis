#include "delayed_ldlt.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static int member(const p58_symbolic *p,uint32_t r){
 uint32_t lo=0,hi=p->nrows;while(lo<hi){uint32_t m=lo+(hi-lo)/2;if(p->rows[m]<r)lo=m+1;else hi=m;}
 return lo<p->nrows&&p->rows[lo]==r;
}
static int to_double(long double a,double *v){if(!isfinite(a)||fabsl(a)>DBL_MAX)return 0;*v=(double)a;return 1;}
static void swap(double *a,double *l,uint32_t *perm,size_t n,size_t k,size_t i,size_t j){
 if(i==j)return;
 for(size_t r=k;r<n;r++){double v=a[i*n+r];a[i*n+r]=a[j*n+r];a[j*n+r]=v;}
 for(size_t r=k;r<n;r++){double v=a[r*n+i];a[r*n+i]=a[r*n+j];a[r*n+j]=v;}
 for(size_t r=0;r<k;r++){double v=l[i*n+r];l[i*n+r]=l[j*n+r];l[j*n+r]=v;}
 uint32_t v=perm[i];perm[i]=perm[j];perm[j]=v;
}
int p60_factorize(p56_context **slot,uint32_t n,const double *input,size_t np,const p58_symbolic *s,const p60_limits *lim,p60_stats *stats){
 if(!lim)return P55_INVALID;
 uint32_t max_n=lim->max_n;long double max_growth=lim->max_growth;
 if(!slot||!n||!input||!s||!np||np>n||!max_n||!lim->max_front||n>max_n||!isfinite(max_growth)||(max_growth!=0&&max_growth<1)||SIZE_MAX/n/n<sizeof(double))return P55_INVALID;
 uint32_t end=0;size_t baseline=0;
 for(size_t p=0;p<np;p++){
  if(s[p].first!=end||!s[p].width||s[p].width>n-end)return P55_INVALID;
  end+=s[p].width;
  if(s[p].nrows>n-end||(s[p].nrows&&!s[p].rows))return P55_INVALID;
  for(uint32_t r=0;r<s[p].nrows;r++)if(s[p].rows[r]<end||s[p].rows[r]>=n||(r&&s[p].rows[r]<=s[p].rows[r-1]))return P55_INVALID;
  baseline+=(size_t)s[p].width*(s[p].width-1)/2+(size_t)s[p].width*s[p].nrows;
 }
 if(end!=n)return P55_INVALID;
 size_t nn=(size_t)n*n;int rc=P55_NOMEM;p60_stats st={0};st.baseline_slots=baseline;
 const long double alpha=(1+sqrtl(17))/8;long double initial=0,peak=0;
 double *a=malloc(nn*sizeof *a),*l=calloc(nn,sizeof *l),*diag=calloc(n,sizeof *diag),*off=calloc(n,sizeof *off);
 uint8_t *kind=calloc(n,1);uint32_t *perm=malloc((size_t)n*sizeof *perm);p55_panel_input *out=calloc(np,sizeof *out);uint32_t *ends=calloc(np,sizeof *ends);size_t used=0;
 if(!a||!l||!diag||!off||!kind||!perm||!out||!ends)goto done;
 rc=P55_INVALID;
 for(size_t i=0;i<n;i++){perm[i]=(uint32_t)i;for(size_t j=0;j<=i;j++)if(!isfinite(input[i*n+j])||input[i*n+j]!=input[j*n+i])goto done;}
 for(size_t p=0;p<np;p++)for(uint32_t j=s[p].first;j<s[p].first+s[p].width;j++)for(uint32_t i=s[p].first+s[p].width;i<n;i++)if(input[(size_t)i*n+j]!=0&&!member(s+p,i)){rc=P58_PATTERN;goto done;}
 memcpy(a,input,nn*sizeof *a);
 for(size_t i=0;i<nn;i++)initial=fmaxl(initial,fabsl(a[i]));
 if(initial==0){rc=P58_NO_PIVOT;goto done;}peak=initial;
 size_t k=0;
 for(size_t p=0;p<np;p++){
  size_t last=(size_t)s[p].first+s[p].width, eligible=last, start=k;
  size_t incoming=(size_t)s[p].first-k;
  if(incoming>st.peak_pending)st.peak_pending=incoming;
  if(last-k>st.peak_front)st.peak_front=last-k;
  if(last-k>lim->max_front){rc=P60_LIMIT;goto done;}
  while(k<eligible){
   long double col=0,akk=fabsl(a[k*n+k]);size_t r=k;
   for(size_t i=k+1;i<n;i++)if(fabsl(a[i*n+k])>col){col=fabsl(a[i*n+k]);r=i;}
   size_t step=1,chosen=k;
   if(akk==0&&col==0){rc=P58_NO_PIVOT;goto done;}
   int branch=0;
   if(akk>=alpha*col){branch=0;}
   else{
    long double row=0;for(size_t i=k;i<n;i++)if(i!=r)row=fmaxl(row,fabsl(a[r*n+i]));
    if(akk>=alpha*col*(col/row)){branch=1;}
    else if(fabsl(a[r*n+r])>=alpha*row){chosen=r;branch=2;}
    else{chosen=r;step=2;branch=3;}
   }
   size_t target=k+step-1;
   if(chosen>=eligible||target>=eligible){
    if(p+1==np){rc=P58_NO_PIVOT;goto done;}
    if(st.delays>=lim->max_delays||last-eligible+1>lim->max_pending){rc=P60_LIMIT;goto done;}
    /* Move current candidate to the pending suffix. The suffix is not
       eligible again until the next original panel is reached. */
    --eligible;if(k!=eligible){swap(a,l,perm,n,k,k,eligible);st.swaps++;}
    st.delays++;if(last-eligible>st.peak_pending)st.peak_pending=last-eligible;
    continue;
   }
   if(branch==0)st.direct++;else if(branch==1)st.secondary++;else if(branch==2)st.alternate++;else st.pairs++;
   if(chosen!=target){swap(a,l,perm,n,k,target,chosen);st.swaps++;}
   if(step==1){
    kind[k]=1;diag[k]=a[k*n+k];if(diag[k]==0){rc=P58_NO_PIVOT;goto done;}
    for(size_t i=k+1;i<n;i++)if(!to_double((long double)a[i*n+k]/diag[k],&l[i*n+k])){rc=P55_RANGE;goto done;}
   }else{
    kind[k]=2;diag[k]=a[k*n+k];off[k]=a[(k+1)*n+k];diag[k+1]=a[(k+1)*n+k+1];
    long double aa=diag[k],bb=off[k],cc=diag[k+1],scale=fmaxl(fabsl(aa),fmaxl(fabsl(bb),fabsl(cc)));
    if(scale==0){rc=P58_NO_PIVOT;goto done;}
    aa/=scale;bb/=scale;cc/=scale;long double det=fmal(aa,cc,-bb*bb);
    if(!isfinite(det)||fabsl(det)<=8*LDBL_EPSILON*(fabsl(aa*cc)+bb*bb)){rc=P58_NO_PIVOT;goto done;}
    for(size_t i=k+2;i<n;i++){
     long double u=(long double)a[i*n+k]/scale,v=(long double)a[i*n+k+1]/scale;
     if(!to_double(fmal(cc,u,-bb*v)/det,&l[i*n+k])||!to_double(fmal(aa,v,-bb*u)/det,&l[i*n+k+1])){rc=P55_RANGE;goto done;}
    }
   }
   for(size_t i=k+step;i<n;i++)for(size_t j=k;j<k+step;j++)st.max_multiplier=fmaxl(st.max_multiplier,fabsl(l[i*n+j]));
   for(size_t i=k+step;i<n;i++)for(size_t j=k+step;j<=i;j++){
    long double v=fmal(-(long double)l[i*n+k],a[j*n+k],a[i*n+j]);
    if(step==2)v=fmal(-(long double)l[i*n+k+1],a[j*n+k+1],v);
    if(!to_double(v,&a[i*n+j])){rc=P55_RANGE;goto done;}a[j*n+i]=a[i*n+j];peak=fmaxl(peak,fabsl(a[i*n+j]));
    if(max_growth!=0&&peak/initial>max_growth){rc=P59_GROWTH;goto done;}
   }
   k+=step;
  }
  if(k>start)ends[used++]=(uint32_t)k;
 }
 if(k!=n){rc=P58_NO_PIVOT;goto done;}
 /* Rebuild symbolic panel row lists in final factor order. Count all
    allocated L slots (including panel padding zeros) before allocation. */
 size_t first=0;
 for(size_t p=0;p<used;p++){
  size_t last=ends[p],w=last-first,nr=0;
  for(size_t i=last;i<n;i++){int hit=0;for(size_t j=first;j<last;j++)if(l[i*n+j]!=0){hit=1;break;}nr+=(size_t)hit;}
  size_t nv=w*(w-1)/2+w*nr;
  st.factor_slots+=nv;
  if(st.factor_slots>baseline&&st.factor_slots-baseline>lim->max_extra_slots){rc=P60_LIMIT;goto done;}
  first=last;
 }
 first=0;
 for(size_t p=0;p<used;p++){
  size_t last=ends[p],w=last-first,nr=0,at=0;
  uint32_t *rows=last<n?malloc((n-last)*sizeof *rows):NULL;
  if(last<n&&!rows){rc=P55_NOMEM;goto done;}
  for(size_t i=last;i<n;i++){for(size_t j=first;j<last;j++)if(l[i*n+j]!=0){rows[nr++]=(uint32_t)i;break;}}
  size_t nv=w*(w-1)/2+w*nr;
  double *v=nv?malloc(nv*sizeof *v):NULL;if(nv&&!v){free(rows);rc=P55_NOMEM;goto done;}
  out[p]=(p55_panel_input){(uint32_t)first,(uint32_t)w,(uint32_t)nr,rows,v,nv};
  for(size_t j=first;j<last;j++){
   for(size_t i=j+1;i<last;i++)v[at++]=l[i*n+j];
   for(size_t r=0;r<nr;r++)v[at++]=l[(size_t)rows[r]*n+j];
  }
  first=last;
 }
 st.panels=used;st.growth=peak/initial;rc=p56_rebuild(slot,n,used,out,kind,diag,off,perm);if(!rc&&stats)*stats=st;
 done:if(out)for(size_t p=0;p<np;p++){free((void*)out[p].values);free((void*)out[p].rows);}free(out);free(ends);free(a);free(l);free(diag);free(off);free(kind);free(perm);return rc;
}
