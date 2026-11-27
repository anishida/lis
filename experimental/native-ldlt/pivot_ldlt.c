#include "pivot_ldlt.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static int member(const p58_symbolic *p,uint32_t r){
 uint32_t lo=0,hi=p->nrows;while(lo<hi){uint32_t m=lo+(hi-lo)/2;if(p->rows[m]<r)lo=m+1;else hi=m;}
 return lo<p->nrows&&p->rows[lo]==r;
}
static int compatible(const p58_symbolic *s,size_t q,uint32_t i,uint32_t j){
 for(size_t p=0;p<q;p++)if(member(s+p,i)!=member(s+p,j))return 0;
 return 1;
}
static int to_double(long double a,double *v){if(!isfinite(a)||fabsl(a)>DBL_MAX)return 0;*v=(double)a;return 1;}
static void swap(double *a,double *l,uint32_t *perm,size_t n,size_t k,size_t i,size_t j){
 if(i==j)return;
 for(size_t r=k;r<n;r++){double v=a[i*n+r];a[i*n+r]=a[j*n+r];a[j*n+r]=v;}
 for(size_t r=k;r<n;r++){double v=a[r*n+i];a[r*n+i]=a[r*n+j];a[r*n+j]=v;}
 for(size_t r=0;r<k;r++){double v=l[i*n+r];l[i*n+r]=l[j*n+r];l[j*n+r]=v;}
 uint32_t v=perm[i];perm[i]=perm[j];perm[j]=v;
}
static int pair_ok(const double *a,size_t n,size_t k,size_t j,double tau,long double *scale,long double *aa,long double *bb,long double *cc,long double *det){
 long double x=a[k*n+k],y=a[j*n+k],z=a[j*n+j];
 *scale=fmaxl(fabsl(x),fmaxl(fabsl(y),fabsl(z)));if(*scale==0)return 0;
 x/=*scale;y/=*scale;z/=*scale;long double d=fmal(x,z,-y*y);
 if(!isfinite(d)||fabsl(d)<=8*LDBL_EPSILON*(fabsl(x*z)+y*y))return 0;
 for(size_t r=k+1;r<n;r++)if(r!=j){
  long double u=(long double)a[r*n+k]/ *scale,v=(long double)a[r*n+j]/ *scale;
  long double f=fmal(z,u,-y*v)/d,g=fmal(x,v,-y*u)/d;
  if(!isfinite(f)||!isfinite(g)||fabsl(f)>DBL_MAX||fabsl(g)>DBL_MAX||fabsl(f)*tau>1||fabsl(g)*tau>1)return 0;
 }
 *aa=x;*bb=y;*cc=z;*det=d;return 1;
}
int p58_factorize(p56_context **slot,uint32_t n,const double *input,size_t np,const p58_symbolic *s,double tau,uint32_t max_n,p58_stats *stats){
 if(!slot||!n||!input||!s||!np||np>n||!max_n||n>max_n||!isfinite(tau)||tau<=0||tau>1||SIZE_MAX/n/n<sizeof(double))return P55_INVALID;
 uint32_t end=0;
 for(size_t p=0;p<np;p++){
  if(s[p].first!=end||!s[p].width||s[p].width>n-end)return P55_INVALID;
  end+=s[p].width;
  if(s[p].nrows>n-end||(s[p].nrows&&!s[p].rows))return P55_INVALID;
  for(uint32_t r=0;r<s[p].nrows;r++)if(s[p].rows[r]<end||s[p].rows[r]>=n||(r&&s[p].rows[r]<=s[p].rows[r-1]))return P55_INVALID;
 }
 if(end!=n)return P55_INVALID;
 size_t nn=(size_t)n*n;int rc=P55_NOMEM;p58_stats st={0};
 double *a=malloc(nn*sizeof *a),*l=calloc(nn,sizeof *l),*diag=calloc(n,sizeof *diag),*off=calloc(n,sizeof *off);
 uint8_t *kind=calloc(n,1);uint32_t *perm=malloc((size_t)n*sizeof *perm);p55_panel_input *out=calloc(np,sizeof *out);
 if(!a||!l||!diag||!off||!kind||!perm||!out)goto done;
 rc=P55_INVALID;
 for(size_t i=0;i<n;i++){perm[i]=(uint32_t)i;for(size_t j=0;j<=i;j++)if(!isfinite(input[i*n+j])||input[i*n+j]!=input[j*n+i])goto done;}
 memcpy(a,input,nn*sizeof *a);
 for(size_t p=0;p<np;p++){
  size_t last=(size_t)s[p].first+s[p].width;
  for(size_t k=s[p].first;k<last;){
   long double col=0;for(size_t i=k+1;i<n;i++)col=fmaxl(col,fabsl(a[i*n+k]));
   size_t step=1;long double scale=0,aa=0,bb=0,cc=0,det=0;
   if(a[k*n+k]!=0&&fabsl(a[k*n+k])>=tau*col){
    kind[k]=1;diag[k]=a[k*n+k];st.scalar++;
    for(size_t i=k+1;i<n;i++)if(!to_double((long double)a[i*n+k]/diag[k],&l[i*n+k])){rc=P55_RANGE;goto done;}
   }else{
    size_t best=n;long double strength=0;
    for(size_t j=k+1;j<last;j++)if(fabsl(a[j*n+k])>strength&&compatible(s,p,(uint32_t)(k+1),(uint32_t)j)){
     long double t,x,y,z,d;if(pair_ok(a,n,k,j,tau,&t,&x,&y,&z,&d)){best=j;strength=fabsl(a[j*n+k]);scale=t;aa=x;bb=y;cc=z;det=d;}
    }
    if(best==n){rc=P58_NO_PIVOT;goto done;}
    if(best!=k+1){swap(a,l,perm,n,k,k+1,best);st.swaps++;}
    kind[k]=2;diag[k]=a[k*n+k];off[k]=a[(k+1)*n+k];diag[k+1]=a[(k+1)*n+k+1];st.pairs++;step=2;
    for(size_t i=k+2;i<n;i++){
     long double u=(long double)a[i*n+k]/scale,v=(long double)a[i*n+k+1]/scale;
     if(!to_double(fmal(cc,u,-bb*v)/det,&l[i*n+k])||!to_double(fmal(aa,v,-bb*u)/det,&l[i*n+k+1])){rc=P55_RANGE;goto done;}
    }
   }
   for(size_t i=k+step;i<n;i++)for(size_t j=k+step;j<=i;j++){
    long double v=fmal(-(long double)l[i*n+k],a[j*n+k],a[i*n+j]);
    if(step==2)v=fmal(-(long double)l[i*n+k+1],a[j*n+k+1],v);
    if(!to_double(v,&a[i*n+j])){rc=P55_RANGE;goto done;}a[j*n+i]=a[i*n+j];
   }
   k+=step;
  }
 }
 for(size_t p=0;p<np;p++){
  size_t first=s[p].first,w=s[p].width,last=first+w,nv=w*(w-1)/2+w*s[p].nrows,at=0;
  double *v=nv?malloc(nv*sizeof *v):NULL;if(nv&&!v){rc=P55_NOMEM;goto done;}
  out[p]=(p55_panel_input){s[p].first,s[p].width,s[p].nrows,s[p].rows,v,nv};
  for(size_t k=first;k<last;k++){
   for(size_t i=k+1;i<last;i++)v[at++]=l[i*n+k];
   for(uint32_t i=(uint32_t)last;i<n;i++)if(!member(s+p,i)&&l[(size_t)i*n+k]!=0){rc=P58_PATTERN;goto done;}
   for(uint32_t r=0;r<s[p].nrows;r++)v[at++]=l[(size_t)s[p].rows[r]*n+k];
  }
 }
 rc=p56_rebuild(slot,n,np,out,kind,diag,off,perm);if(!rc&&stats)*stats=st;
 done:if(out)for(size_t p=0;p<np;p++)free((void*)out[p].values);free(out);free(a);free(l);free(diag);free(off);free(kind);free(perm);return rc;
}
