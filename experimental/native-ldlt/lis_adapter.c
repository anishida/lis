#include "lis_adapter.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "factor_private.h"
void p56_destroy(p56_context *c) {
    if(c) { p55_free(c->factor); free(c->perm); free(c); }
}
uint32_t p56_size(const p56_context *c) { return c ? c->n : 0; }
size_t p56_block_count(const p56_context *c) { return c ? p55_block_count(c->factor) : 0; }
int p56_rebuild(p56_context **slot,uint32_t n,size_t np,
                const p55_panel_input *panels,const uint8_t *kind,
                const double *diag,const double *off,const uint32_t *perm) {
    if(!slot || !n || SIZE_MAX/n < sizeof(double)) return P55_INVALID;
    p56_context *c=calloc(1,sizeof *c);
    if(!c) return P55_NOMEM;
    c->n=n;
    int rc=P55_OK;
    if(perm) {
        if(SIZE_MAX/n < sizeof(uint32_t)) {rc=P55_INVALID;goto fail;}
        unsigned char *seen=calloc(n,1);
        if(!seen){rc=P55_NOMEM;goto fail;}
        int identity=1;
        for(uint32_t k=0;k<n;k++) {
            uint32_t j=perm[k];
            if(j>=n || seen[j]){rc=P55_INVALID;break;}
            seen[j]=1; if(j!=k)identity=0;
        }
        free(seen);
        if(rc)goto fail;
        if(!identity){
            c->perm=malloc((size_t)n*sizeof *c->perm);
            if(!c->perm){rc=P55_NOMEM;goto fail;}
            memcpy(c->perm,perm,(size_t)n*sizeof *c->perm);
        }
    }
    rc=p55_pack(&c->factor,n,np,panels,kind,diag,off);
    if(rc)goto fail;
    p56_context *old=*slot; *slot=c; p56_destroy(old);
    return P55_OK;
fail:
    p56_destroy(c); return rc;
}
/* uintptr_t is used for practical flat-address C targets; no pointer arithmetic
   on unrelated arrays. Difference comparisons avoid address-end overflow. */
static int overlaps(const void *a,const void *b,size_t bytes) {
    uintptr_t x=(uintptr_t)a,y=(uintptr_t)b;
    return (x>=y ? x-y : y-x)<bytes;
}
int p56_apply(const p56_context *c,const double *rhs,double *out,
              double *scratch,size_t scratch_count) {
    if(!c||!rhs||!out||!scratch||scratch_count<c->n)return P55_INVALID;
    size_t bytes=(size_t)c->n*sizeof(double);
    if(overlaps(scratch,rhs,bytes)||overlaps(scratch,out,bytes))return P55_INVALID;
    if(!c->perm)memcpy(scratch,rhs,bytes);
    else for(uint32_t k=0;k<c->n;k++)scratch[k]=rhs[c->perm[k]];
    int rc=p55_solve(c->factor,scratch);
    if(rc)return rc;
    if(!c->perm)memcpy(out,scratch,bytes);
    else for(uint32_t k=0;k<c->n;k++)out[c->perm[k]]=scratch[k];
    return P55_OK;
}
