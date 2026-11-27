#include "qmd_plan.h"
#include "supernodal_ldlt.h"
#include "factor_private.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
struct p71_plan {
 uint32_t n;size_t nnz,ng,bytes;
 size_t *col,*pc;uint32_t *row,*pr,*scatter;
 p63_order *order;p62_symbolic *symbolic;
};
typedef struct {uint32_t col,row,source;} record;
static int add(size_t *z,size_t n,size_t s,size_t cap){if(n>(cap-*z)/s)return 0;*z+=n*s;return 1;}
static int cmp(const void *a,const void *b){const record *x=a,*y=b;if(x->col!=y->col)return x->col<y->col?-1:1;return (x->row>y->row)-(x->row<y->row);}
void p71_plan_free(p71_plan *p){if(p){free(p->col);free(p->pc);free(p->row);free(p->pr);free(p->scatter);p63_free(p->order);p62_free(p->symbolic);free(p);}}
int p71_plan_create(p71_plan **out,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const p63_limits *ol,const p62_limits *sl,uint32_t width,size_t cap){
 if(!out||!n||n==UINT32_MAX||nnz>UINT32_MAX||!col||(nnz&&!row)||!ol||!sl||!width||!cap)return P55_INVALID;
 size_t bytes=0,np=(size_t)n+1;
 /* Also count retained ordering and symbolic arrays; their construction
    scratch is bounded by the separate ol/sl limits. */
 if(!add(&bytes,1,sizeof(p71_plan)+sizeof(p63_order)+sizeof(p62_symbolic),cap)||!add(&bytes,np,2*sizeof(size_t)+sizeof(uint32_t),cap)||!add(&bytes,n,4*sizeof(uint32_t),cap)||!add(&bytes,nnz,3*sizeof(uint32_t)+sizeof(record),cap))return P60_LIMIT;
 p71_plan *p=calloc(1,sizeof *p);record *e=NULL;int rc=P55_NOMEM;
 if(!p)return rc;
 rc=p63q_analyze(&p->order,n,nnz,col,row,ol);if(rc)goto done;
 rc=P55_NOMEM;p->n=n;p->nnz=nnz;p->bytes=bytes;
 p->col=malloc(np*sizeof *p->col);p->pc=calloc(np,sizeof *p->pc);
 p->row=nnz?malloc(nnz*sizeof *p->row):NULL;p->pr=nnz?malloc(nnz*sizeof *p->pr):NULL;p->scatter=nnz?malloc(nnz*sizeof *p->scatter):NULL;e=nnz?malloc(nnz*sizeof *e):NULL;
 if(!p->col||!p->pc||(nnz&&(!p->row||!p->pr||!p->scatter||!e)))goto done;
 memcpy(p->col,col,np*sizeof *col);if(nnz)memcpy(p->row,row,nnz*sizeof *row);
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++){uint32_t x=p->order->inv[j],y=p->order->inv[row[z]];e[z]=(record){x<y?x:y,x<y?y:x,(uint32_t)z};}
 if(nnz)qsort(e,nnz,sizeof *e,cmp);
 for(size_t z=0;z<nnz;z++){p->pc[(size_t)e[z].col+1]++;p->pr[z]=e[z].row;p->scatter[z]=e[z].source;}
 for(uint32_t j=0;j<n;j++)p->pc[j+1]+=p->pc[j];
 free(e);e=NULL;
 rc=p62_analyze(&p->symbolic,n,nnz,p->pc,p->pr,sl);if(rc)goto done;
 p62_symbolic *s=p->symbolic;size_t ng=0;
 for(size_t g=0;g<s->nsuper;g++){uint32_t w=s->super[g+1]-s->super[g];ng+=1+(w-1)/width;}
 size_t k=ng;uint32_t hi=n;s->super[k]=n;
 for(size_t g=s->nsuper;g>0;g--){uint32_t lo=s->super[g-1],w=hi-lo;size_t pieces=1+(w-1)/width;
  for(size_t j=pieces;j>0;j--)s->super[--k]=lo+(uint32_t)(j-1)*width;
  hi=lo;
 }
 p->ng=ng;p71_plan_free(*out);*out=p;p=NULL;rc=P55_OK;
 done:free(e);p71_plan_free(p);return rc;
}
int p71_refactor(const p71_plan *p,p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p61_limits *nl,p63_stats *stats){
 if(!p||!slot||!nl||n!=p->n||nnz!=p->nnz||!col||(nnz&&(!row||!val)))return P55_INVALID;
 if(memcmp(col,p->col,((size_t)n+1)*sizeof *col)||(nnz&&memcmp(row,p->row,nnz*sizeof *row)))return P55_INVALID;
 /* Refresh conversion workspace is separately bounded by numeric limit. */
 size_t bytes=0;if(!add(&bytes,nnz,sizeof(double),nl->max_work_bytes)||!add(&bytes,n,sizeof(uint32_t),nl->max_work_bytes))return P60_LIMIT;
 double *a=nnz?malloc(nnz*sizeof *a):NULL;uint32_t *perm=malloc((size_t)n*sizeof *perm);p56_context *ctx=NULL;p61_stats st={0};int rc=P55_NOMEM;
 if(!perm||(nnz&&!a))goto done;
 for(size_t z=0;z<nnz;z++)a[z]=val[p->scatter[z]];
 rc=p65_numeric(&ctx,n,p->pc,p->pr,a,p->symbolic->parent,p->ng==n?NULL:p->symbolic->super,p->ng,nl,&st);if(rc)goto done;
 int identity=1;for(uint32_t k=0;k<n;k++){perm[k]=p->order->perm[ctx->perm?ctx->perm[k]:k];if(perm[k]!=k)identity=0;}
 free(ctx->perm);ctx->perm=identity?NULL:perm;if(!identity)perm=NULL;
 if(stats){p63_stats q={0};q.summary=*p->order;q.summary.perm=q.summary.inv=NULL;q.permutation_bytes=p->bytes;q.numeric=st;*stats=q;}
 p56_destroy(*slot);*slot=ctx;ctx=NULL;
 done:free(a);free(perm);p56_destroy(ctx);return rc;
}
