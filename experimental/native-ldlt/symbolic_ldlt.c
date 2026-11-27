#include "symbolic_ldlt.h"
#include <stdlib.h>
#include <string.h>
#define NONE UINT32_MAX
void p62_free(p62_symbolic *s){if(s){free(s->parent);free(s->count);free(s->super);free(s);}}
static int addbytes(size_t *total,size_t n,size_t width,size_t limit){
 if(n>(limit-*total)/width)return 0;
 *total+=n*width;return 1;
}
int p62_analyze(p62_symbolic **out,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const p62_limits *lim){
 if(!out||!n||n==NONE||!col||!lim||col[0]||col[n]!=nnz||(nnz&&!row))return P55_INVALID;
 size_t off=0,bytes=0;int rc=P55_NOMEM;
 for(uint32_t j=0;j<n;j++){
  if(col[j]>col[j+1]||col[j+1]>nnz)return P55_INVALID;
  for(size_t z=col[j];z<col[j+1];z++){
   uint32_t i=row[z];if(i<j||i>=n||(z>col[j]&&i<=row[z-1]))return P55_INVALID;
   off+=i!=j;
  }
 }
 size_t np=(size_t)n+1;
 if(!addbytes(&bytes,1,sizeof(p62_symbolic),lim->max_work_bytes)||
    !addbytes(&bytes,np,sizeof(size_t),lim->max_work_bytes)||
    !addbytes(&bytes,n,sizeof(size_t),lim->max_work_bytes)||
    !addbytes(&bytes,off,sizeof(uint32_t),lim->max_work_bytes)||
    !addbytes(&bytes,n,3*sizeof(uint32_t),lim->max_work_bytes)||
    !addbytes(&bytes,np,sizeof(uint32_t),lim->max_work_bytes))return P60_LIMIT;
 p62_symbolic *s=calloc(1,sizeof *s);size_t *rp=calloc(np,sizeof *rp),*cursor=malloc((size_t)n*sizeof *cursor);
 uint32_t *edges=off?malloc(off*sizeof *edges):NULL,*mark=malloc((size_t)n*sizeof *mark);
 if(!s||!rp||!cursor||(off&&!edges)||!mark)goto done;
 s->parent=malloc((size_t)n*sizeof *s->parent);s->count=malloc((size_t)n*sizeof *s->count);s->super=malloc(np*sizeof *s->super);
 if(!s->parent||!s->count||!s->super)goto done;
 s->n=n;s->peak_work_bytes=bytes;
 /* Transpose only strict lower structure; row k lists earlier neighbours. */
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)if(row[z]!=j)rp[(size_t)row[z]+1]++;
 for(uint32_t i=0;i<n;i++)rp[i+1]+=rp[i];
 memcpy(cursor,rp,(size_t)n*sizeof *cursor);
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)if(row[z]!=j)edges[cursor[row[z]]++]=j;
 for(uint32_t i=0;i<n;i++){s->parent[i]=mark[i]=NONE;s->count[i]=1;}
 /* Elimination tree by ancestor path compression. */
 for(uint32_t k=0;k<n;k++)for(size_t z=rp[k];z<rp[k+1];z++){
  uint32_t j=edges[z];while(j!=NONE&&j<k){uint32_t next=mark[j];mark[j]=k;if(next==NONE)s->parent[j]=k;j=next;}
 }
 /* Exact structural reach for each row, independent of numerical cancellation. */
 for(uint32_t i=0;i<n;i++)mark[i]=NONE;
 for(uint32_t k=0;k<n;k++)for(size_t z=rp[k];z<rp[k+1];z++){
  uint32_t j=edges[z];while(j!=NONE&&j<k&&mark[j]!=k){
   if(s->slots==lim->max_slots){rc=P60_LIMIT;goto done;}
   mark[j]=k;s->count[j]++;s->slots++;j=s->parent[j];
  }
 }
 /* For consecutive columns with parent[j]=j+1, the later pattern of j
    is a subset of the pattern of j+1 plus row j+1. Equal tail sizes imply equality. */
 s->super[s->nsuper++]=0;
 for(uint32_t j=0;j<n;j++){
  if(s->parent[j]==NONE)s->roots++;
  if(s->count[j]>s->max_front)s->max_front=s->count[j];
  if(j+1<n&&(s->parent[j]!=j+1||s->count[j]!=s->count[j+1]+1u))s->super[s->nsuper++]=j+1;
 }
 s->super[s->nsuper]=n;
 p62_free(*out);*out=s;s=NULL;rc=P55_OK;
 done:free(rp);free(cursor);free(edges);free(mark);p62_free(s);return rc;
}
int p62_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p62_limits *slim,const p61_limits *nlim,p61_stats *stats){
 if(!slot||!nlim||(nnz&&!val))return P55_INVALID;
 p62_symbolic *s=NULL;int rc=p62_analyze(&s,n,nnz,col,row,slim);
 if(!rc)rc=p61_factorize(slot,n,col,row,val,s->parent,nlim,stats);
 p62_free(s);return rc;
}

#include "supernodal_ldlt.h"
int p65_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p62_limits *sl,const p61_limits *nl,uint32_t max_width,p61_stats *stats){
 if(!slot||!nl||!max_width||(nnz&&!val))return P55_INVALID;
 if(max_width==1)return p62_factorize(slot,n,nnz,col,row,val,sl,nl,stats);
 p62_symbolic *s=NULL;int rc=p62_analyze(&s,n,nnz,col,row,sl);if(rc)return rc;
 /* Reuse the already allocated n+1 boundary array. Fill backwards so
    expansion cannot overwrite an unread original boundary. */
 size_t ng=0;
 for(size_t g=0;g<s->nsuper;g++){uint32_t w=s->super[g+1]-s->super[g];ng+=1+(w-1)/max_width;}
 if(ng==n){rc=p61_factorize(slot,n,col,row,val,s->parent,nl,stats);p62_free(s);return rc;}
 size_t k=ng;uint32_t hi=n;s->super[k]=n;
 for(size_t g=s->nsuper;g>0;g--){uint32_t lo=s->super[g-1];uint32_t w=hi-lo;size_t pieces=1+(w-1)/max_width;
  for(size_t j=pieces;j>0;j--)s->super[--k]=lo+(uint32_t)(j-1)*max_width;
  hi=lo;
 }
 rc=p65_numeric(slot,n,col,row,val,s->parent,s->super,ng,nl,stats);p62_free(s);return rc;
}
