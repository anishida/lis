#include "supernodal_ldlt.h"
#include "factor_private.h"
#include <stdlib.h>
#include <string.h>
#define NONE UINT32_MAX
typedef struct {uint32_t *v;size_t len,cap;} adj;
typedef struct {size_t used,peak,limit;} budget;
static int charge(budget *b,size_t n,size_t w){if(n>(b->limit-b->used)/w)return 0;b->used+=n*w;if(b->used>b->peak)b->peak=b->used;return 1;}
static int reserve(adj *a,size_t need,budget *b){
 if(need<=a->cap)return 0;
 size_t cap=a->cap?a->cap:2;while(cap<need){if(cap>SIZE_MAX/2)return P60_LIMIT;cap*=2;}
 /* Count old and new buffers simultaneously to conservatively bound resize. */
 if(!charge(b,cap,sizeof(uint32_t)))return P60_LIMIT;
 uint32_t *v=malloc(cap*sizeof *v);if(!v){b->used-=cap*sizeof *v;return P55_NOMEM;}
 if(a->len)memcpy(v,a->v,a->len*sizeof *v);
 free(a->v);b->used-=a->cap*sizeof *v;a->v=v;a->cap=cap;return 0;
}
static int less(const adj *g,uint32_t a,uint32_t b){return g[a].len<g[b].len||(g[a].len==g[b].len&&a<b);}
static void exchange(uint32_t *h,uint32_t *pos,size_t i,size_t j){uint32_t t=h[i];h[i]=h[j];h[j]=t;pos[h[i]]=(uint32_t)i;pos[h[j]]=(uint32_t)j;}
static void down(const adj *g,uint32_t *h,uint32_t *pos,size_t len,size_t i){
 while(i<len/2){size_t j=2*i+1;if(j+1<len&&less(g,h[j+1],h[j]))j++;if(!less(g,h[j],h[i]))break;exchange(h,pos,i,j);i=j;}
}
static void fix(const adj *g,uint32_t *h,uint32_t *pos,size_t len,size_t i){
 while(i&&less(g,h[i],h[(i-1)/2])){exchange(h,pos,i,(i-1)/2);i=(i-1)/2;}
 down(g,h,pos,len,i);
}
static int valid(uint32_t n,size_t nnz,const size_t *col,const uint32_t *row){
 if(!n||n==NONE||!col||col[0]||col[n]!=nnz||(nnz&&!row))return 0;
 for(uint32_t j=0;j<n;j++){if(col[j]>col[j+1]||col[j+1]>nnz)return 0;for(size_t z=col[j];z<col[j+1];z++)if(row[z]<j||row[z]>=n||(z>col[j]&&row[z]<=row[z-1]))return 0;}return 1;
}
void p63_free(p63_order *o){if(o){free(o->perm);free(o->inv);free(o);}}
int p63_analyze(p63_order **out,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const p63_limits *lim){
 if(!out||!lim||!valid(n,nnz,col,row))return P55_INVALID;
 budget b={0,0,lim->max_work_bytes};int rc=P55_NOMEM;
 if(!charge(&b,1,sizeof(p63_order))||!charge(&b,n,sizeof(adj))||!charge(&b,n,5*sizeof(uint32_t)))return P60_LIMIT;
 p63_order *o=calloc(1,sizeof *o);adj *g=calloc(n,sizeof *g);uint32_t *h=malloc((size_t)n*sizeof *h),*pos=malloc((size_t)n*sizeof *pos),*mark=malloc((size_t)n*sizeof *mark);
 if(!o||!g||!h||!pos||!mark)goto done;
 o->perm=malloc((size_t)n*sizeof *o->perm);o->inv=malloc((size_t)n*sizeof *o->inv);if(!o->perm||!o->inv)goto done;o->n=n;
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)if(row[z]!=j){uint32_t i=row[z];g[j].len++;g[i].len++;}
 for(uint32_t j=0;j<n;j++){size_t d=g[j].len;g[j].len=0;if((rc=reserve(&g[j],d,&b)))goto done;h[j]=pos[j]=j;mark[j]=NONE;}
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)if(row[z]!=j){uint32_t i=row[z];g[j].v[g[j].len++]=i;g[i].v[g[i].len++]=j;}
 for(size_t i=n/2;i>0;i--)down(g,h,pos,n,i-1);
 for(uint32_t k=0;k<n;k++){
  size_t len=(size_t)n-k;uint32_t v=h[0];size_t d=g[v].len;
  if(d>lim->max_slots-o->slots){rc=P60_LIMIT;goto done;}o->slots+=d;if(d+1>o->max_front)o->max_front=d+1;
  o->perm[k]=v;o->inv[v]=k;exchange(h,pos,0,len-1);pos[v]=NONE;len--;if(len)fix(g,h,pos,len,0);
  /* Remove v first, then create clique. Maintain heap after EVERY degree change. */
  for(size_t i=0;i<d;i++){uint32_t a=g[v].v[i];adj *ga=&g[a];size_t j=0;while(j<ga->len&&ga->v[j]!=v)j++;if(j==ga->len){rc=P55_INVALID;goto done;}ga->v[j]=ga->v[--ga->len];fix(g,h,pos,len,pos[a]);}
  for(size_t i=0;i<d;i++){
   uint32_t a=g[v].v[i];adj *ga=&g[a];
   /* Mark via the unique vertex a: it can be visited in multiple eliminations,
      so clear touched marks after this row. */
   for(size_t z=0;z<ga->len;z++)mark[ga->v[z]]=a;
   for(size_t j=i+1;j<d;j++){
    if(o->pair_checks==lim->max_pair_checks){rc=P60_LIMIT;goto done;}o->pair_checks++;
    uint32_t c=g[v].v[j];if(mark[c]==a)continue;adj *gc=&g[c];
    if((rc=reserve(ga,ga->len+1,&b))||(rc=reserve(gc,gc->len+1,&b)))goto done;
    ga->v[ga->len++]=c;fix(g,h,pos,len,pos[a]);gc->v[gc->len++]=a;fix(g,h,pos,len,pos[c]);
   }
   for(size_t z=0;z<ga->len;z++)mark[ga->v[z]]=NONE;
  }
  free(g[v].v);b.used-=g[v].cap*sizeof(uint32_t);g[v]=(adj){0};
 }
 o->peak_work_bytes=b.peak;p63_free(*out);*out=o;o=NULL;rc=0;
 done:if(g)for(uint32_t i=0;i<n;i++)free(g[i].v);free(g);free(h);free(pos);free(mark);p63_free(o);return rc;
}
/* Quotient-graph minimum degree: active eliminated elements represent
   cliques implicitly. Heap keys are lower bounds on the exact degree.
   Adapted from the quotient-graph strategy in LIS PR #273, with bounded
   allocations and deterministic smallest-index ties as in p63_analyze. */
static void qfree(adj *a,budget *b){free(a->v);b->used-=a->cap*sizeof(uint32_t);*a=(adj){0};}
static int qless(const size_t *key,uint32_t a,uint32_t b){return key[a]<key[b]||(key[a]==key[b]&&a<b);}
static void qdown(const size_t *key,uint32_t *h,uint32_t *pos,size_t len,size_t i){
 while(i<len/2){size_t j=2*i+1;if(j+1<len&&qless(key,h[j+1],h[j]))j++;if(!qless(key,h[j],h[i]))break;exchange(h,pos,i,j);i=j;}
}
static void qfix(const size_t *key,uint32_t *h,uint32_t *pos,size_t len,size_t i){
 while(i&&qless(key,h[i],h[(i-1)/2])){exchange(h,pos,i,(i-1)/2);i=(i-1)/2;}qdown(key,h,pos,len,i);
}
static int visit(p63_order *o,const p63_limits *l){if(o->pair_checks==l->max_pair_checks)return 0;o->pair_checks++;return 1;}
static int neighbors(uint32_t v,const adj *g,const adj *inc,adj *el,const uint32_t *pos,uint32_t *seen,uint32_t stamp,uint32_t *list,size_t *d,p63_order *o,const p63_limits *l){
 *d=0;seen[v]=stamp;
 for(size_t i=0;i<g[v].len;i++){if(!visit(o,l))return P60_LIMIT;uint32_t u=g[v].v[i];if(pos[u]!=NONE&&seen[u]!=stamp){seen[u]=stamp;list[(*d)++]=u;}}
 for(size_t i=0;i<inc[v].len;i++){if(!visit(o,l))return P60_LIMIT;adj *a=&el[inc[v].v[i]];for(size_t j=0;j<a->len;j++){if(!visit(o,l))return P60_LIMIT;uint32_t u=a->v[j];if(pos[u]!=NONE&&seen[u]!=stamp){seen[u]=stamp;list[(*d)++]=u;}}}
 return 0;
}
int p63q_analyze(p63_order **out,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const p63_limits *lim){
 if(!out||!lim||!valid(n,nnz,col,row))return P55_INVALID;
 budget b={0,0,lim->max_work_bytes};int rc=P55_NOMEM;
 if(!charge(&b,1,sizeof(p63_order))||!charge(&b,n,3*sizeof(adj)+6*sizeof(uint32_t)+sizeof(size_t)))return P60_LIMIT;
 p63_order *o=calloc(1,sizeof *o);adj *g=calloc(n,sizeof *g),*inc=calloc(n,sizeof *inc),*el=calloc(n,sizeof *el);
 uint32_t *h=malloc((size_t)n*sizeof *h),*pos=malloc((size_t)n*sizeof *pos),*seen=calloc(n,sizeof *seen),*list=malloc((size_t)n*sizeof *list);size_t *key=calloc(n,sizeof *key);uint32_t stamp=0;
 if(!o||!g||!inc||!el||!h||!pos||!seen||!list||!key)goto done;
 o->perm=malloc((size_t)n*sizeof *o->perm);o->inv=malloc((size_t)n*sizeof *o->inv);if(!o->perm||!o->inv)goto done;o->n=n;
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)if(row[z]!=j){key[j]++;key[row[z]]++;}
 for(uint32_t j=0;j<n;j++){if((rc=reserve(&g[j],key[j],&b)))goto done;h[j]=pos[j]=j;}
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)if(row[z]!=j){uint32_t i=row[z];g[j].v[g[j].len++]=i;g[i].v[g[i].len++]=j;}
 for(size_t i=n/2;i>0;i--)qdown(key,h,pos,n,i-1);
 for(uint32_t k=0;k<n;k++){
  size_t len=(size_t)n-k,d=0;uint32_t v;
  for(;;){v=h[0];if(++stamp==0){memset(seen,0,(size_t)n*sizeof *seen);stamp=1;}rc=neighbors(v,g,inc,el,pos,seen,stamp,list,&d,o,lim);if(rc)goto done;
   if(d<key[v]){rc=P55_INVALID;goto done;}if(d==key[v])break;key[v]=d;qfix(key,h,pos,len,0);}
  if(d>lim->max_slots-o->slots){rc=P60_LIMIT;goto done;}o->slots+=d;if(d+1>o->max_front)o->max_front=d+1;
  o->perm[k]=v;o->inv[v]=k;exchange(h,pos,0,len-1);pos[v]=NONE;len--;if(len)qfix(key,h,pos,len,0);
  for(size_t j=0;j<inc[v].len;j++){qfree(&el[inc[v].v[j]],&b);}
  qfree(&inc[v],&b);qfree(&g[v],&b);
  if(d){if((rc=reserve(&el[v],d,&b)))goto done;memcpy(el[v].v,list,d*sizeof *list);el[v].len=d;}
  for(size_t j=0;j<d;j++){uint32_t u=list[j];adj *a=&inc[u];if(a->len==a->cap){size_t t=0;for(size_t i=0;i<a->len;i++){if(!visit(o,lim)){rc=P60_LIMIT;goto done;}if(el[a->v[i]].v)a->v[t++]=a->v[i];}a->len=t;}
   if((rc=reserve(a,a->len+1,&b)))goto done;
   a->v[a->len++]=v;
   size_t old=key[u]?key[u]-1:0;key[u]=old>d-1?old:d-1;qfix(key,h,pos,len,pos[u]);}
 }
 o->peak_work_bytes=b.peak;p63_free(*out);*out=o;o=NULL;rc=0;
 done:for(uint32_t i=0;i<n;i++){if(g)free(g[i].v);if(inc)free(inc[i].v);if(el)free(el[i].v);}free(g);free(inc);free(el);free(h);free(pos);free(seen);free(list);free(key);p63_free(o);return rc;
}

typedef struct {uint32_t col,row;double value;} entry;
static int cmp(const void *a,const void *b){const entry *x=a,*y=b;if(x->col!=y->col)return x->col<y->col?-1:1;return (x->row>y->row)-(x->row<y->row);}
static int factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p63_limits *ol,const p62_limits *sl,const p61_limits *nl,p63_stats *stats,int quotient,uint32_t width){
 if(!slot||!ol||!sl||!nl||(nnz&&!val))return P55_INVALID;
 p63_order *o=NULL;int rc=(quotient?p63q_analyze:p63_analyze)(&o,n,nnz,col,row,ol);if(rc)return rc;
 budget b={0,0,ol->max_work_bytes};entry *e=NULL;size_t *cp=NULL;uint32_t *ri=NULL,*perm=NULL;double *a=NULL;p56_context *ctx=NULL;p61_stats st={0};
 if(!charge(&b,1,sizeof *o)||!charge(&b,n,2*sizeof(uint32_t))||!charge(&b,(size_t)n+1,sizeof *cp)||!charge(&b,nnz,sizeof *e+sizeof *ri+sizeof *a)||!charge(&b,n,sizeof *perm)){rc=P60_LIMIT;goto end;}
 e=nnz?malloc(nnz*sizeof *e):NULL;cp=calloc((size_t)n+1,sizeof *cp);ri=nnz?malloc(nnz*sizeof *ri):NULL;a=nnz?malloc(nnz*sizeof *a):NULL;perm=malloc((size_t)n*sizeof *perm);
 if(!cp||!perm||(nnz&&(!e||!ri||!a))){rc=P55_NOMEM;goto end;}
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++){uint32_t x=o->inv[j],y=o->inv[row[z]];e[z]=(entry){x<y?x:y,x<y?y:x,val[z]};}
 if(nnz)qsort(e,nnz,sizeof *e,cmp);
 for(size_t z=0;z<nnz;z++){cp[(size_t)e[z].col+1]++;ri[z]=e[z].row;a[z]=e[z].value;}for(uint32_t j=0;j<n;j++)cp[j+1]+=cp[j];free(e);e=NULL;
 rc=p65_factorize(&ctx,n,nnz,cp,ri,a,sl,nl,width,&st);if(rc)goto end;
 int identity=1;for(uint32_t k=0;k<n;k++){perm[k]=o->perm[ctx->perm?ctx->perm[k]:k];if(perm[k]!=k)identity=0;}
 free(ctx->perm);ctx->perm=identity?NULL:perm;if(!identity)perm=NULL;
 if(stats){p63_stats s={0};s.summary=*o;s.summary.perm=s.summary.inv=NULL;s.permutation_bytes=b.peak;s.numeric=st;*stats=s;}
 p56_destroy(*slot);*slot=ctx;ctx=NULL;
 end:free(e);free(cp);free(ri);free(a);free(perm);p63_free(o);p56_destroy(ctx);return rc;
}

int p63_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p63_limits *ol,const p62_limits *sl,const p61_limits *nl,p63_stats *stats){return factorize(slot,n,nnz,col,row,val,ol,sl,nl,stats,0,1);}
int p63q_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p63_limits *ol,const p62_limits *sl,const p61_limits *nl,uint32_t width,p63_stats *stats){return factorize(slot,n,nnz,col,row,val,ol,sl,nl,stats,1,width);}
