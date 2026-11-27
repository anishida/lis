#include "packed_ldlt.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>

#include "factor_private.h"
static int mul(size_t a,size_t b,size_t *v) {
    if(b && a>SIZE_MAX/b) return 0;
    *v=a*b; return 1;
}
static uint32_t row_at(const panel *p,uint32_t r) {
    return p->short_index ? p->first+p->width+((const uint16_t*)p->rows)[r]
                          : ((const uint32_t*)p->rows)[r];
}
void p55_free(p55_factor *f) {
    if(!f) return;
    if(f->p) for(size_t i=0;i<f->np;i++) { free(f->p[i].rows); free(f->p[i].values); }
    free(f->p); free(f->diag); free(f->blocks); free(f);
}
size_t p55_block_count(const p55_factor *f) { return f?f->nb:0; }
size_t p55_index_bytes(const p55_factor *f) { return f?f->index_bytes:0; }
int p55_pack(p55_factor **out,uint32_t n,size_t np,const p55_panel_input *in,
             const uint8_t *kind,const double *diag,const double *off) {
    int rc=P55_INVALID; size_t bytes, nb=0; uint32_t next=0;
    if(!out) return rc;
    *out=NULL;
    if(!n || !np || np>n || !in || !kind || !diag) return rc;
    p55_factor *f=calloc(1,sizeof *f);
    if(!f) return P55_NOMEM;
    f->n=n; f->np=np;
    if(!mul(np,sizeof(panel),&bytes)) goto fail;
    f->p=calloc(1,bytes);
    if(!mul(n,sizeof(double),&bytes)) goto fail;
    f->diag=malloc(bytes);
    if(!f->p || !f->diag) {rc=P55_NOMEM;goto fail;}
    memcpy(f->diag,diag,bytes);
    for(uint32_t k=0;k<n;) {
        if(!isfinite(diag[k])) goto fail;
        if(kind[k]==1) { if(diag[k]==0) {rc=P55_SINGULAR;goto fail;} k++; }
        else if(kind[k]==2 && k+1<n && kind[k+1]==0 && off && isfinite(off[k]) && isfinite(diag[k+1])) { nb++; k+=2; }
        else goto fail;
    }
    if(nb) {
        if(!mul(nb,sizeof(block),&bytes)) goto fail;
        f->blocks=malloc(bytes); if(!f->blocks){rc=P55_NOMEM;goto fail;}
    }
    f->nb=nb;
    for(uint32_t k=0,j=0;k<n;) {
        if(kind[k]==1) {k++;continue;}
        long double a=diag[k],b=off[k],c=diag[k+1];
        long double s=fmaxl(fabsl(a),fmaxl(fabsl(b),fabsl(c)));
        if(s==0){rc=P55_SINGULAR;goto fail;}
        a/=s;b/=s;c/=s;
        long double det=fmal(a,c,-b*b);
        /* Conservative local guard; not a global pivot acceptance criterion. */
        if(!isfinite(det)||fabsl(det)<=8*LDBL_EPSILON*(fabsl(a*c)+b*b)) {rc=P55_SINGULAR;goto fail;}
        f->blocks[j++]=(block){k,a,b,c,det,s,off[k]}; k+=2;
    }
    for(size_t i=0;i<np;i++) {
        const p55_panel_input *q=&in[i]; panel *p=&f->p[i];
        if(q->first!=next || !q->width || q->width>n-next) goto fail;
        uint32_t end=next+q->width;
        if(kind[next]==0 || kind[end-1]==2 || q->nrows>n-end || (q->nrows&&!q->rows)) goto fail;
        p->first=next;p->width=q->width;p->nrows=q->nrows;next=end;
        size_t w=q->width, t,u,nv;
        if(!mul(w%2?w:w/2,w%2?(w-1)/2:w-1,&t)) goto fail;
        if(!mul(w,q->nrows,&u)||t>SIZE_MAX-u) goto fail;
        nv=t+u; p->nvalues=nv;
        if(nv!=q->nvalues || (nv&&!q->values)) goto fail;
        for(uint32_t r=0;r<q->nrows;r++) if(q->rows[r]<end || q->rows[r]>=n || (r&&q->rows[r]<=q->rows[r-1])) goto fail;
        p->short_index=!q->nrows || q->rows[q->nrows-1]-end<=UINT16_MAX;
        if(!mul(q->nrows,p->short_index?sizeof(uint16_t):sizeof(uint32_t),&bytes) || bytes>SIZE_MAX-f->index_bytes) goto fail;
        f->index_bytes+=bytes;
        if(bytes) {
            p->rows=malloc(bytes); if(!p->rows){rc=P55_NOMEM;goto fail;}
            if(p->short_index) for(uint32_t r=0;r<q->nrows;r++) ((uint16_t*)p->rows)[r]=(uint16_t)(q->rows[r]-end);
            else memcpy(p->rows,q->rows,bytes);
        }
        size_t at=0;
        for(uint32_t j=0;j<p->width;j++) {
            if(kind[p->first+j]==2 && q->values[at]!=0) goto fail;
            size_t count=(size_t)(p->width-j-1)+p->nrows;
            for(size_t z=0;z<count;z++) if(!isfinite(q->values[at+z])) goto fail;
            at+=count;
        }
        if(!mul(nv,sizeof(double),&bytes)) goto fail;
        if(bytes) {p->values=malloc(bytes); if(!p->values){rc=P55_NOMEM;goto fail;} memcpy(p->values,q->values,bytes);}
    }
    if(next!=n) goto fail;
    *out=f;return P55_OK;
fail:p55_free(f);return rc;
}
int p55_solve(const p55_factor *f,double *x) {
    if(!f||!x) return P55_INVALID;
    for(uint32_t k=0;k<f->n;k++) if(!isfinite(x[k])) return P55_RANGE;
    for(size_t i=0;i<f->np;i++) {
        const panel *p=&f->p[i]; size_t at=0;
        for(uint32_t j=0;j<p->width;j++) {
            uint32_t k=p->first+j; double v=x[k];
            for(uint32_t l=j+1;l<p->width;l++) x[p->first+l]-=p->values[at++]*v;
            if(p->short_index) for(uint32_t r=0;r<p->nrows;r++) x[row_at(p,r)]-=p->values[at++]*v;
            else for(uint32_t r=0;r<p->nrows;r++) x[((const uint32_t*)p->rows)[r]]-=p->values[at++]*v;
        }
    }
    if(!f->nb) {
        for(uint32_t k=0;k<f->n;k++) x[k]/=f->diag[k];
    } else {
        uint32_t k=0;
        for(size_t j=0;j<f->nb;j++) {
            const block *b=&f->blocks[j];
            for(;k<b->k;k++) x[k]/=f->diag[k];
            long double u=(long double)x[k]/b->scale,v=(long double)x[k+1]/b->scale;
            long double y=fmal(b->c,u,-b->b*v)/b->det,z=fmal(b->a,v,-b->b*u)/b->det;
            if(!isfinite(y)||!isfinite(z)||fabsl(y)>DBL_MAX||fabsl(z)>DBL_MAX) return P55_RANGE;
            x[k]=(double)y;x[k+1]=(double)z;k+=2;
        }
        for(;k<f->n;k++) x[k]/=f->diag[k];
    }
    for(size_t i=f->np;i-->0;) {
        const panel *p=&f->p[i]; size_t at=p->nvalues;
        for(uint32_t j=p->width;j-->0;) {
            at-=(size_t)(p->width-j-1)+p->nrows;size_t z=at;
            double v=x[p->first+j];
            for(uint32_t l=j+1;l<p->width;l++) v-=p->values[z++]*x[p->first+l];
            if(p->short_index) for(uint32_t r=0;r<p->nrows;r++) v-=p->values[z++]*x[row_at(p,r)];
            else for(uint32_t r=0;r<p->nrows;r++) v-=p->values[z++]*x[((const uint32_t*)p->rows)[r]];
            x[p->first+j]=v;
        }
    }
    for(uint32_t k=0;k<f->n;k++) if(!isfinite(x[k])) return P55_RANGE;
    return P55_OK;
}
