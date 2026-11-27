#include "supernodal_ldlt.h"
#include "schur_simd.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stddef.h>
#define NONE UINT32_MAX
typedef union {max_align_t align;size_t size;} hdr;
typedef struct {size_t live,peak,limit;int error;} budget;
static void *alloc(budget *b,size_t n,size_t s){
 if(!n)return NULL;
 if(s>SIZE_MAX/n||n*s>SIZE_MAX-sizeof(hdr)||n*s>b->limit-b->live){b->error=P60_LIMIT;return NULL;}
 size_t z=n*s;hdr *h=calloc(1,sizeof *h+z);if(!h){b->error=P55_NOMEM;return NULL;}
 h->size=z;b->live+=z;if(b->live>b->peak)b->peak=b->live;return h+1;
}
static void release(budget *b,void *p){if(p){hdr *h=(hdr*)p-1;b->live-=h->size;free(h);}}
typedef struct {uint32_t m,pending;uint32_t *id;double *a;} update;
typedef struct {uint32_t row;double x,y;} entry;
typedef struct {uint32_t first,w,nr;entry *e;} pivot;
static void clear(budget *b,update *u){release(b,u->id);release(b,u->a);*u=(update){0};}
/* All arguments are absolute values of finite quantities, or running
   nonnegative maxima. Avoid a libm call in the hot assembly/pivot loops. */
static inline long double max_finite(long double a,long double b){return a>b?a:b;}
static int cv(long double x,double *v){if(!isfinite(x)||fabsl(x)>DBL_MAX)return 0;*v=(double)x;return 1;}
static void sw(double *a,uint32_t *id,size_t m,size_t k,size_t i,size_t j){
 if(i==j)return;
 for(size_t r=k;r<m;r++){double v=a[i*m+r];a[i*m+r]=a[j*m+r];a[j*m+r]=v;}
 for(size_t r=k;r<m;r++){double v=a[r*m+i];a[r*m+i]=a[r*m+j];a[r*m+j]=v;}
 uint32_t t=id[i];id[i]=id[j];id[j]=t;
}
static int cmp(const void *a,const void *b){uint32_t x=((const entry*)a)->row,y=((const entry*)b)->row;return (x>y)-(x<y);}
static int add(uint32_t id,uint32_t *ids,uint32_t *map,size_t *m,size_t cap){
 if(map[id]!=NONE)return 1;
 if(*m>=cap)return 0;
 map[id]=(uint32_t)*m;ids[(*m)++]=id;return 1;
}
int p65_numeric(p56_context **slot,uint32_t n,const size_t *col,const uint32_t *row,const double *val,const uint32_t *parent,const uint32_t *bounds,size_t ng,const p61_limits *lim,p61_stats *stats){
 if(!slot||!n||n==NONE||!col||!parent||!lim||!lim->max_front||!lim->max_work_bytes||col[0]||(!row&&col[n])||(!val&&col[n])||!isfinite(lim->max_growth)||(lim->max_growth&&lim->max_growth<1))return P55_INVALID;
 budget b={0,0,lim->max_work_bytes,P55_NOMEM};p61_stats st={0};int rc=P55_INVALID;
 uint32_t *head=NULL,*next=NULL,*tin=NULL,*tout=NULL,*stack=NULL,*map=NULL,*ids=NULL,*perm=NULL,*inv=NULL,*group=NULL;uint8_t *kind=NULL;
 double *diag=NULL,*off=NULL,*a=NULL;update *u=NULL;pivot *pv=NULL;p55_panel_input *out=NULL;
 size_t used=0,done=0;long double initial=0,peak=0;size_t cap=lim->max_front<n?lim->max_front:n;
 #define NEW(p,count) do{p=alloc(&b,(count),sizeof *(p));if(!(p)){rc=b.error;goto exit;}}while(0)
 NEW(head,n);NEW(next,n);NEW(tin,n);NEW(tout,n);NEW(stack,n);NEW(map,n);NEW(ids,cap);NEW(perm,n);NEW(inv,n);NEW(kind,n);NEW(diag,n);NEW(off,n);NEW(u,n);NEW(pv,n);NEW(out,n);
 for(uint32_t i=0;i<n;i++){head[i]=map[i]=NONE; if(parent[i]!=NONE&&(parent[i]<=i||parent[i]>=n))goto exit; if(col[i]>col[i+1])goto exit;}
 for(uint32_t i=0;i<n;i++){next[i]=NONE;if(parent[i]!=NONE){next[i]=head[parent[i]];head[parent[i]]=i;}}
 /* DFS entry/exit intervals, without recursion. Reuse inv for child cursor. */
 uint32_t tick=0;
 for(uint32_t root=0;root<n;root++)if(parent[root]==NONE){size_t depth=0;stack[depth++]=root;tin[root]=tick++;inv[root]=head[root];while(depth){uint32_t v=stack[depth-1],c=inv[v];if(c==NONE){tout[v]=tick;depth--;continue;}inv[v]=next[c];tin[c]=tick++;inv[c]=head[c];stack[depth++]=c;}}
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++){
  uint32_t i=row[z];if(i<j||i>=n||(z>col[j]&&i<=row[z-1])||!isfinite(val[z])||!(tin[i]<=tin[j]&&tin[j]<tout[i]))goto exit;
  initial=max_finite(initial,fabsl(val[z]));
 }
 if(bounds){
  if(!ng||ng>n||bounds[0]!=0||bounds[ng]!=n)goto exit;
  NEW(group,n);
  for(size_t g=0;g<ng;g++){
   uint32_t lo=bounds[g],hi=bounds[g+1];if(lo>=hi||hi>n)goto exit;
   for(uint32_t j=lo;j<hi;j++){if(j+1<hi&&parent[j]!=j+1)goto exit;group[j]=hi-1;}
  }
  /* Contract chains into representative (last column) front nodes. */
  for(uint32_t i=0;i<n;i++)head[i]=next[i]=NONE;
  for(size_t g=0;g<ng;g++){uint32_t f=bounds[g+1]-1,p=parent[f];if(p!=NONE){p=group[p];next[f]=head[p];head[p]=f;}}
 }
 if(initial==0){rc=P58_NO_PIVOT;goto exit;}peak=initial;
 release(&b,tin);tin=NULL;release(&b,tout);tout=NULL;release(&b,stack);stack=NULL;
 const long double alpha=(1+sqrtl(17))/8;
 p66_kernel schur=p66_select_kernel();
 for(size_t g=0;g<(bounds?ng:n);g++){
  uint32_t lo=bounds?bounds[g]:(uint32_t)g,f=bounds?bounds[g+1]-1:(uint32_t)g;
  size_t m=0,k=0,eligible;
  /* Candidate prefix: original column and each child's delayed variables. */
  for(uint32_t j=lo;j<=f;j++)if(!add(j,ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  for(uint32_t c=head[f];c!=NONE;c=next[c])for(uint32_t i=0;i<u[c].pending;i++)if(!add(u[c].id[i],ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  eligible=m;
  for(size_t z=col[lo];z<col[f+1];z++)if(!add(row[z],ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  for(uint32_t c=head[f];c!=NONE;c=next[c])for(uint32_t i=u[c].pending;i<u[c].m;i++)if(!add(u[c].id[i],ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  if(m>st.peak_front)st.peak_front=m;
  if(m>SIZE_MAX/m){rc=P60_LIMIT;goto exit;}NEW(a,m*m);
  for(uint32_t j=lo;j<=f;j++)for(size_t z=col[j];z<col[j+1];z++){size_t i=map[row[z]],h=map[j];a[i*m+h]=a[h*m+i]=val[z];}
  for(uint32_t c=head[f];c!=NONE;c=next[c]){
   update *q=&u[c];for(size_t i=0;i<q->m;i++)for(size_t j=0;j<=i;j++){
    size_t x=map[q->id[i]],y=map[q->id[j]];
    if(!cv((long double)a[x*m+y]+q->a[i*q->m+j],&a[x*m+y])){rc=P55_RANGE;goto exit;}a[y*m+x]=a[x*m+y];
    peak=max_finite(peak,fabsl(a[x*m+y]));
   }
   clear(&b,q);
  }
  for(size_t i=0;i<m;i++)map[ids[i]]=NONE;
  if(lim->max_growth&&peak/initial>lim->max_growth){rc=P59_GROWTH;goto exit;}
  size_t candidates=eligible;
  while(k<eligible){
   long double akk=fabsl(a[k*m+k]),colmax=0;size_t r=k,w=1,chosen=k;
   for(size_t i=k+1;i<m;i++)if(fabsl(a[i*m+k])>colmax){colmax=fabsl(a[i*m+k]);r=i;}
   if(akk==0&&colmax==0){rc=P58_NO_PIVOT;goto exit;}
   if(akk<alpha*colmax){long double rm=0;for(size_t i=k;i<m;i++)if(i!=r)rm=max_finite(rm,fabsl(a[r*m+i]));
    if(akk<alpha*colmax*(colmax/rm)){chosen=r;if(fabsl(a[r*m+r])<alpha*rm)w=2;}}
   size_t target=k+w-1;
   if(chosen>=eligible||target>=eligible){
    if(parent[f]==NONE){rc=P58_NO_PIVOT;goto exit;}
    if(st.delays>=lim->max_delays){rc=P60_LIMIT;goto exit;}--eligible;
    if(k!=eligible){sw(a,ids,m,k,k,eligible);st.swaps++;}st.delays++;continue;
   }
   if(chosen!=target){sw(a,ids,m,k,target,chosen);st.swaps++;}
   long double aa=a[k*m+k],bb=0,cc=0,s=1,det=1;
   if(w==2){bb=a[(k+1)*m+k];cc=a[(k+1)*m+k+1];s=max_finite(fabsl(aa),max_finite(fabsl(bb),fabsl(cc)));aa/=s;bb/=s;cc/=s;det=fmal(aa,cc,-bb*bb);if(!isfinite(det)||fabsl(det)<=8*LDBL_EPSILON*(fabsl(aa*cc)+bb*bb)){rc=P58_NO_PIVOT;goto exit;}}
   kind[done]=(uint8_t)w;diag[done]=a[k*m+k];perm[done]=ids[k];
   if(w==2){diag[done+1]=a[(k+1)*m+k+1];off[done]=a[(k+1)*m+k];perm[done+1]=ids[k+1];st.pairs++;}else st.scalars++;
   pivot *p=&pv[used++];p->first=(uint32_t)done;p->w=(uint32_t)w;
   /* Store multipliers in unused upper part of the pivot rows; lower
      columns retain original coupling for the Schur update. */
   for(size_t i=k+w;i<m;i++){
    long double x=a[i*m+k],y=w==2?a[i*m+k+1]:0;
    if(w==1)x/=aa;else{long double z=x/s,v=y/s;x=fmal(cc,z,-bb*v)/det;y=fmal(aa,v,-bb*z)/det;}
    if(!cv(x,&a[k*m+i])||(w==2&&!cv(y,&a[(k+1)*m+i]))){rc=P55_RANGE;goto exit;}
    st.max_multiplier=max_finite(st.max_multiplier,max_finite(fabsl(x),fabsl(y)));
    if(a[k*m+i]!=0||(w==2&&a[(k+1)*m+i]!=0))p->nr++;
   }
   size_t slots=(size_t)p->nr*w+(w-1);
   if(slots>lim->max_factor_slots-st.factor_slots){rc=P60_LIMIT;goto exit;}st.factor_slots+=slots;
   if(p->nr){NEW(p->e,p->nr);size_t at=0;for(size_t i=k+w;i<m;i++)if(a[k*m+i]!=0||(w==2&&a[(k+1)*m+i]!=0))p->e[at++]=(entry){ids[i],a[k*m+i],w==2?a[(k+1)*m+i]:0};}
   if(!schur(a,m,k,w,&peak)){rc=P55_RANGE;goto exit;}
   if(lim->max_growth&&peak/initial>lim->max_growth){rc=P59_GROWTH;goto exit;}
   k+=w;done+=w;
  }
  if(k<m){
   if(parent[f]==NONE){rc=P58_NO_PIVOT;goto exit;}
   update *q=&u[f];q->m=(uint32_t)(m-k);q->pending=(uint32_t)(candidates-k);NEW(q->id,q->m);NEW(q->a,(size_t)q->m*q->m);
   memcpy(q->id,ids+k,(size_t)q->m*sizeof *ids);for(size_t i=0;i<q->m;i++)memcpy(q->a+i*q->m,a+(i+k)*m+k,(size_t)q->m*sizeof *a);
  }
  release(&b,a);a=NULL;st.fronts++;
 }
 if(done!=n){rc=P58_NO_PIVOT;goto exit;}
 for(uint32_t i=0;i<n;i++)inv[perm[i]]=i;
 for(size_t j=0;j<used;j++){
  pivot *p=&pv[j];for(uint32_t i=0;i<p->nr;i++)p->e[i].row=inv[p->e[i].row];
  if(p->nr)qsort(p->e,p->nr,sizeof *p->e,cmp);
  uint32_t *rr=NULL;double *vv=NULL;size_t nv=(size_t)p->nr*p->w+p->w-1;
  out[j]=(p55_panel_input){p->first,p->w,p->nr,NULL,NULL,nv};
  if(p->nr){NEW(rr,p->nr);out[j].rows=rr;}
  if(nv){NEW(vv,nv);out[j].values=vv;}
  for(uint32_t i=0;i<p->nr;i++){rr[i]=p->e[i].row;vv[i+p->w-1]=p->e[i].x;if(p->w==2)vv[1+p->nr+i]=p->e[i].y;}
  release(&b,p->e);p->e=NULL;
 }
 st.growth=peak/initial;st.peak_work_bytes=b.peak;
 rc=p56_rebuild(slot,n,used,out,kind,diag,off,perm);if(!rc&&stats)*stats=st;
 exit:
 release(&b,a);if(u)for(uint32_t i=0;i<n;i++)clear(&b,&u[i]);if(pv)for(size_t i=0;i<used;i++)release(&b,pv[i].e);
 if(out)for(size_t i=0;i<used;i++){release(&b,(void*)out[i].rows);release(&b,(void*)out[i].values);}
 release(&b,group);release(&b,head);release(&b,next);release(&b,tin);release(&b,tout);release(&b,stack);release(&b,map);release(&b,ids);release(&b,perm);release(&b,inv);release(&b,kind);release(&b,diag);release(&b,off);release(&b,u);release(&b,pv);release(&b,out);return rc;
 #undef NEW
}

int p61_factorize(p56_context **slot,uint32_t n,const size_t *col,const uint32_t *row,const double *val,const uint32_t *parent,const p61_limits *lim,p61_stats *stats){return p65_numeric(slot,n,col,row,val,parent,NULL,0,lim,stats);}
