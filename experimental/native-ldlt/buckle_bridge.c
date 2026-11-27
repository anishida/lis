#include "buckle_bridge.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

void p69_defaults(p69_options *o){
 if(!o)return;
 *o=(p69_options){.max_conversion_bytes=256u*1024*1024,
 .ordering={256u*1024*1024,100000000,500000000,8,64,1000000},
 .symbolic={128u*1024*1024,100000000},
 .numeric={1024,256u*1024*1024,100000000,1000000,0},
 .method=P64_AUTO,.width=16,.threads=1};
}
int p69_validate_csr(const p69_csr *a){
 if(!a||!a->n||!a->ptr||a->base>1||a->storage<P69_FULL||a->storage>P69_UPPER||a->nnz>SIZE_MAX-a->base||(a->nnz&&(!a->index||!a->value)))return P55_INVALID;
 if(a->ptr[0]!=a->base||a->ptr[a->n]!=a->nnz+a->base)return P55_INVALID;
 for(uint32_t i=0;i<a->n;i++){
  if(a->ptr[i]<a->base||a->ptr[i]>a->ptr[i+1]||a->ptr[i+1]>a->nnz+a->base)return P55_INVALID;
  uint32_t prev=0;
  for(size_t z=a->ptr[i]-a->base;z<a->ptr[i+1]-a->base;z++){
   uint32_t raw=a->index[z];if(raw<a->base)return P55_INVALID;
   uint32_t j=raw-a->base;
   if(j>=a->n||(z>a->ptr[i]-a->base&&j<=prev)||!isfinite(a->value[z])||(a->storage==P69_LOWER&&j>i)||(a->storage==P69_UPPER&&j<i))return P55_INVALID;
   prev=j;
  }
 }
 if(a->storage==P69_FULL)for(uint32_t i=0;i<a->n;i++)for(size_t z=a->ptr[i]-a->base;z<a->ptr[i+1]-a->base;z++){
  uint32_t j=a->index[z]-a->base;if(i==j)continue;
  size_t lo=a->ptr[j]-a->base,hi=a->ptr[j+1]-a->base,end=hi;uint32_t target=i+a->base;
  while(lo<hi){size_t m=lo+(hi-lo)/2;if(a->index[m]<target)lo=m+1;else hi=m;}
  if(lo==end||a->index[lo]!=target||a->value[lo]!=a->value[z])return P55_INVALID;
 }
 return P55_OK;
}
/* Emit lower CSC entries (c,r). Sorting/merging allows mixed triangle layouts. */
typedef struct {uint32_t c,r;double x;unsigned src;} entry;
static size_t count(const p69_csr *a){size_t n=0;if(a)for(uint32_t i=0;i<a->n;i++)for(size_t z=a->ptr[i]-a->base;z<a->ptr[i+1]-a->base;z++)if(a->storage!=P69_FULL||a->index[z]-a->base>=i)n++;return n;}
static void emit(const p69_csr *a,entry *e,size_t *at,unsigned src){if(!a)return;for(uint32_t i=0;i<a->n;i++)for(size_t z=a->ptr[i]-a->base;z<a->ptr[i+1]-a->base;z++){
 uint32_t j=a->index[z]-a->base;if(a->storage==P69_FULL&&j<i)continue;e[(*at)++]=(entry){i<j?i:j,i<j?j:i,a->value[z],src};}}
static int cmp(const void *aa,const void *bb){const entry *a=aa,*b=bb;if(a->c!=b->c)return a->c<b->c?-1:1;if(a->r!=b->r)return a->r<b->r?-1:1;return (a->src>b->src)-(a->src<b->src);}
static int charge(size_t *bytes,size_t n,size_t s,size_t cap){if(n>(cap-*bytes)/s)return 0;*bytes+=n*s;return 1;}
int p69_factorize(p56_context **slot,const p69_csr *k,const p69_csr *g,double sigma,const p69_options *o,p69_stats *stats){
 if(!slot||!o||!isfinite(sigma)||(!g&&sigma!=0)||!o->width||!o->threads||o->threads>64||o->method<P64_NATURAL||o->method>P64_AUTO)return P55_INVALID;
 int rc=p69_validate_csr(k);if(rc)return rc;if(g){rc=p69_validate_csr(g);if(rc||g->n!=k->n)return P55_INVALID;}
 if((uint64_t)k->n+1>SIZE_MAX/sizeof(size_t))return P60_LIMIT;
 size_t nk=count(k),ng=sigma==0?0:count(g);if(nk>SIZE_MAX-ng)return P60_LIMIT;size_t total=nk+ng,bytes=0;
 if(!charge(&bytes,(size_t)k->n+1,sizeof(size_t),o->max_conversion_bytes)||!charge(&bytes,total,sizeof(entry),o->max_conversion_bytes)||!charge(&bytes,total,sizeof(uint32_t),o->max_conversion_bytes)||!charge(&bytes,total,sizeof(double),o->max_conversion_bytes))return P60_LIMIT;
 entry *e=total?malloc(total*sizeof *e):NULL;size_t *cp=calloc((size_t)k->n+1,sizeof *cp);uint32_t *ri=total?malloc(total*sizeof *ri):NULL;double *v=total?malloc(total*sizeof *v):NULL;
 if(!cp||(total&&(!e||!ri||!v))){rc=P55_NOMEM;goto end;}
 size_t at=0;emit(k,e,&at,0);if(sigma!=0)emit(g,e,&at,1);if(total)qsort(e,total,sizeof *e,cmp);
 size_t m=0;
 for(size_t z=0;z<total;){uint32_t c=e[z].c,r=e[z].r;long double a=0,b=0;
  do{if(e[z].src)b=e[z].x;else a=e[z].x;z++;}while(z<total&&e[z].c==c&&e[z].r==r);
  long double x=a-(long double)sigma*b;if(!isfinite(x)||fabsl(x)>DBL_MAX){rc=P55_RANGE;goto end;}
  /* Retain structural zeros, including exact cancellation. */
  v[m]=(double)x;ri[m++]=r;cp[c+1]++;
 }
 for(uint32_t i=0;i<k->n;i++)cp[i+1]+=cp[i];
 free(e);e=NULL;
 p69_stats st={.input_k=k->nnz,.input_g=g?g->nnz:0,.shifted_nnz=m,.conversion_bytes=bytes};
 rc=p67_ordered_factorize(slot,k->n,m,cp,ri,v,&o->ordering,&o->symbolic,&o->numeric,o->method,o->width,o->threads,&st.ordering,&st.parallel);
 if(!rc&&stats)*stats=st;
 end:free(e);free(cp);free(ri);free(v);return rc;
}
