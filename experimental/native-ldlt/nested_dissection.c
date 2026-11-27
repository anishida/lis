#include "supernodal_ldlt.h"
#include "factor_private.h"
#include <stdlib.h>
#include <string.h>
#define NONE UINT32_MAX
static int charge(size_t *b,size_t n,size_t w,size_t limit){if(n>(limit-*b)/w)return 0;*b+=n*w;return 1;}
static int valid(uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri){
 if(!n||n==NONE||!cp||cp[0]||cp[n]!=nnz||(nnz&&!ri))return 0;
 for(uint32_t j=0;j<n;j++){if(cp[j]>cp[j+1]||cp[j+1]>nnz)return 0;for(size_t z=cp[j];z<cp[j+1];z++)if(ri[z]<j||ri[z]>=n||(z>cp[j]&&ri[z]<=ri[z-1]))return 0;}return 1;
}
typedef struct {size_t *rp,*mark,stamp;uint32_t *edges,*nodes,*queue,*dist,*tmp;uint32_t leaf,depth;} graph;
static size_t bfs(graph *g,uint32_t start,size_t tail,size_t stamp){
 size_t head=tail;g->queue[tail++]=start;g->dist[start]=0;
 while(head<tail){uint32_t v=g->queue[head++];for(size_t z=g->rp[v];z<g->rp[v+1];z++){uint32_t w=g->edges[z];if(g->mark[w]==stamp&&g->dist[w]==NONE){g->dist[w]=g->dist[v]+1;g->queue[tail++]=w;}}}return tail;
}
static void split(graph *g,size_t lo,size_t hi,uint32_t depth){
 size_t n=hi-lo;if(n<=g->leaf||depth>=g->depth)return;
 size_t stamp=++g->stamp;for(size_t i=lo;i<hi;i++){uint32_t v=g->nodes[i];g->mark[v]=stamp;g->dist[v]=NONE;}
 size_t tail=0,bound=0,balance=SIZE_MAX;
 for(size_t i=lo;i<hi;i++)if(g->dist[g->nodes[i]]==NONE){tail=bfs(g,g->nodes[i],tail,stamp);if(tail<n){size_t diff=tail>n-tail?tail-(n-tail):(n-tail)-tail;if(diff<balance){balance=diff;bound=tail;}}}
 if(bound){memcpy(g->nodes+lo,g->queue,n*sizeof(uint32_t));split(g,lo,lo+bound,depth+1);split(g,lo+bound,hi,depth+1);return;}
 /* Second sweep from a peripheral vertex. BFS queue is in level order. */
 uint32_t start=g->queue[n-1];for(size_t i=lo;i<hi;i++)g->dist[g->nodes[i]]=NONE;bfs(g,start,0,stamp);
 size_t left=0,right=0,best=SIZE_MAX;
 for(size_t a=0;a<n;){size_t b=a+1;while(b<n&&g->dist[g->queue[b]]==g->dist[g->queue[a]])b++;
  if(a&&b<n){size_t m=a>n-b?a:n-b;size_t score=m+2*(b-a);if(score<best){best=score;left=a;right=b;}}a=b;}
 if(!left)return; /* clique or two BFS levels: natural leaf */
 /* [left][right][separator]. No edges join left and right. */
 memcpy(g->tmp,g->queue,left*sizeof(uint32_t));memcpy(g->tmp+left,g->queue+right,(n-right)*sizeof(uint32_t));memcpy(g->tmp+left+n-right,g->queue+left,(right-left)*sizeof(uint32_t));
 memcpy(g->nodes+lo,g->tmp,n*sizeof(uint32_t));split(g,lo,lo+left,depth+1);split(g,lo+left,lo+left+n-right,depth+1);
}
static int nd(uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,const p64_limits *l,p63_order **out){
 size_t off=0,bytes=0;for(uint32_t j=0;j<n;j++)for(size_t z=cp[j];z<cp[j+1];z++)off+=ri[z]!=j;
 if(off>SIZE_MAX/2||!charge(&bytes,1,sizeof(p63_order),l->max_work_bytes)||!charge(&bytes,(size_t)n+1,sizeof(size_t),l->max_work_bytes)||!charge(&bytes,n,2*sizeof(size_t)+5*sizeof(uint32_t),l->max_work_bytes)||!charge(&bytes,2*off,sizeof(uint32_t),l->max_work_bytes))return P60_LIMIT;
 (void)nnz; graph g={0};g.leaf=l->leaf_size;g.depth=l->max_depth;int rc=P55_NOMEM;p63_order *o=calloc(1,sizeof *o);size_t *cur=malloc((size_t)n*sizeof *cur);
 g.rp=calloc((size_t)n+1,sizeof *g.rp);g.mark=calloc(n,sizeof *g.mark);g.edges=off?malloc(2*off*sizeof *g.edges):NULL;g.nodes=malloc((size_t)n*sizeof *g.nodes);g.queue=malloc((size_t)n*sizeof *g.queue);g.dist=malloc((size_t)n*sizeof *g.dist);g.tmp=malloc((size_t)n*sizeof *g.tmp);
 if(!o||!cur||!g.rp||!g.mark||(off&&!g.edges)||!g.nodes||!g.queue||!g.dist||!g.tmp)goto end;
 o->inv=malloc((size_t)n*sizeof *o->inv);if(!o->inv)goto end;
 for(uint32_t j=0;j<n;j++)for(size_t z=cp[j];z<cp[j+1];z++)if(ri[z]!=j){g.rp[j+1]++;g.rp[ri[z]+1]++;}
 for(uint32_t j=0;j<n;j++){g.rp[j+1]+=g.rp[j];g.nodes[j]=j;}memcpy(cur,g.rp,(size_t)n*sizeof *cur);
 for(uint32_t j=0;j<n;j++)for(size_t z=cp[j];z<cp[j+1];z++)if(ri[z]!=j){g.edges[cur[j]++]=ri[z];g.edges[cur[ri[z]]++]=j;}
 split(&g,0,n,0);o->perm=g.nodes;g.nodes=NULL;o->n=n;o->peak_work_bytes=bytes;for(uint32_t j=0;j<n;j++)o->inv[o->perm[j]]=j;*out=o;o=NULL;rc=0;
 end:free(cur);free(g.rp);free(g.mark);free(g.edges);free(g.nodes);free(g.queue);free(g.dist);free(g.tmp);p63_free(o);return rc;
}
typedef struct {uint32_t c,r;double v;} entry;
static int cmp(const void *a,const void *b){const entry *x=a,*y=b;if(x->c!=y->c)return x->c<y->c?-1:1;return (x->r>y->r)-(x->r<y->r);}
static int reorder(uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,const double *v,const p63_order *o,size_t limit,size_t **pc,uint32_t **pr,double **pv,size_t *peak){
 size_t bytes=0;if(!charge(&bytes,(size_t)n+1,sizeof(size_t),limit)||!charge(&bytes,n,2*sizeof(uint32_t),limit)||!charge(&bytes,nnz,sizeof(entry)+sizeof(uint32_t)+(v?sizeof(double):0),limit))return P60_LIMIT;
 entry *e=nnz?malloc(nnz*sizeof *e):NULL;size_t *c=calloc((size_t)n+1,sizeof *c);uint32_t *r=nnz?malloc(nnz*sizeof *r):NULL;double *a=v&&nnz?malloc(nnz*sizeof *a):NULL;
 if(!c||(nnz&&(!e||!r||(v&&!a)))){free(e);free(c);free(r);free(a);return P55_NOMEM;}
 for(uint32_t j=0;j<n;j++)for(size_t z=cp[j];z<cp[j+1];z++){uint32_t x=o->inv[j],y=o->inv[ri[z]];e[z]=(entry){x<y?x:y,x<y?y:x,v?v[z]:0};}
 if(nnz)qsort(e,nnz,sizeof *e,cmp);
 for(size_t z=0;z<nnz;z++){c[e[z].c+1]++;r[z]=e[z].r;if(a)a[z]=e[z].v;}for(uint32_t j=0;j<n;j++)c[j+1]+=c[j];free(e);*pc=c;*pr=r;*pv=a;*peak=bytes;return 0;
}
static int measure(uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,p63_order *o,const p64_limits *l,const p62_limits *sl,p64_estimate *e){
 size_t *c=NULL,peak=0;uint32_t *r=NULL;double *unused=NULL;p62_symbolic *s=NULL;int rc=0;
 if(o){rc=reorder(n,nnz,cp,ri,NULL,o,l->max_work_bytes,&c,&r,&unused,&peak);if(rc)goto end;}
 rc=p62_analyze(&s,n,nnz,o?c:cp,o?r:ri,sl);if(rc)goto end;
 e->slots=s->slots;e->front=s->max_front;e->work_bytes=peak>s->peak_work_bytes?peak:s->peak_work_bytes;if(o&&o->peak_work_bytes>e->work_bytes)e->work_bytes=o->peak_work_bytes;
 for(uint32_t j=0;j<n;j++)e->work+=(long double)s->count[j]*s->count[j];
 if(s->slots>l->max_slots||s->max_front>l->max_front){rc=P60_LIMIT;goto end;}if(o){o->slots=s->slots;o->max_front=s->max_front;}
 end:free(c);free(r);p62_free(s);return rc;
}
int p64_order(uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,const p64_limits *l,const p62_limits *sl,p64_method method,p63_order **out,p64_stats *stats){
 if(!out||!l||!sl||!l->leaf_size||!l->max_depth||l->max_depth>64||method<P64_NATURAL||method>P64_AUTO||!valid(n,nnz,cp,ri))return P55_INVALID;
 p64_stats st={0};p63_order *best=NULL;int chosen=-1,rc=P60_LIMIT;for(int i=0;i<3;i++)st.candidates[i].status=-1;
 for(int i=0;i<3;i++){if(method!=P64_AUTO&&(int)method!=i)continue;p63_order *o=NULL;int r=0;
  if(i==P64_MD){p63_limits ml={l->max_work_bytes,l->max_slots,l->max_pair_checks};r=p63_analyze(&o,n,nnz,cp,ri,&ml);}else if(i==P64_ND)r=nd(n,nnz,cp,ri,l,&o);
  if(!r)r=measure(n,nnz,cp,ri,o,l,sl,&st.candidates[i]);
  st.candidates[i].status=r;
  if(r){p63_free(o);if(r!=P60_LIMIT){p63_free(best);return r;}continue;}
  p64_estimate *e=&st.candidates[i];if(chosen<0||e->work<st.candidates[chosen].work||(e->work==st.candidates[chosen].work&&e->slots<st.candidates[chosen].slots)){p63_free(best);best=o;chosen=i;}else p63_free(o);
  if(i==P64_NATURAL&&e->front<=2)break; /* natural tree/chain, no fill */
 }
 if(chosen<0)return rc;
 if(!best){size_t b=0;if(!charge(&b,1,sizeof *best,l->max_work_bytes)||!charge(&b,n,2*sizeof(uint32_t),l->max_work_bytes))return P60_LIMIT;best=calloc(1,sizeof *best);if(!best)return P55_NOMEM;best->perm=malloc((size_t)n*sizeof *best->perm);best->inv=malloc((size_t)n*sizeof *best->inv);if(!best->perm||!best->inv){p63_free(best);return P55_NOMEM;}best->n=n;best->slots=st.candidates[chosen].slots;best->max_front=st.candidates[chosen].front;best->peak_work_bytes=b;for(uint32_t j=0;j<n;j++)best->perm[j]=best->inv[j]=j;}
 st.selected=(p64_method)chosen;p63_free(*out);*out=best;if(stats)*stats=st;return 0;
}
int p65_ordered_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,const double *v,const p64_limits *l,const p62_limits *sl,const p61_limits *nl,p64_method method,uint32_t max_width,p64_stats *stats){
 if(!slot||!nl||!max_width||(nnz&&!v))return P55_INVALID;
 p64_stats st={0};p63_order *o=NULL;int rc=p64_order(n,nnz,cp,ri,l,sl,method,&o,&st);if(rc)return rc;
 size_t *c=NULL,peak=0;uint32_t *r=NULL;double *a=NULL;p56_context *ctx=NULL;
 if(st.selected==P64_NATURAL)rc=p65_factorize(&ctx,n,nnz,cp,ri,v,sl,nl,max_width,&st.numeric);
 else{rc=reorder(n,nnz,cp,ri,v,o,l->max_work_bytes,&c,&r,&a,&peak);if(!rc)rc=p65_factorize(&ctx,n,nnz,c,r,a,sl,nl,max_width,&st.numeric);if(!rc){int identity=1;for(uint32_t j=0;j<n;j++){o->inv[j]=o->perm[ctx->perm?ctx->perm[j]:j];if(o->inv[j]!=j)identity=0;}free(ctx->perm);ctx->perm=identity?NULL:o->inv;if(!identity)o->inv=NULL;}}
 if(!rc){p56_destroy(*slot);*slot=ctx;ctx=NULL;if(stats)*stats=st;}p56_destroy(ctx);p63_free(o);free(c);free(r);free(a);return rc;
}

int p64_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,const double *v,const p64_limits *l,const p62_limits *sl,const p61_limits *nl,p64_method method,p64_stats *stats){return p65_ordered_factorize(slot,n,nnz,cp,ri,v,l,sl,nl,method,1,stats);}
