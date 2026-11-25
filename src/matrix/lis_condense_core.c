#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_condense_core.h"
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Memory belonging to a plan is kept fully independent of caller-owned CSR. */
struct lis_sc_plan {
    lis_sc_csr a, s;
    size_t nphys, nq, nb;
    size_t *u2full, *full2u, *block_of, *local_of;
    size_t *block_ptr, *block_dof;
    size_t *neighbors_ptr, *neighbors;
    size_t *lu_ptr, *piv_r, *piv_c;
    double *lu;
    double pivot_tol;
};

typedef struct { size_t block, phys; } link_pair;
typedef struct { size_t row, col; } matrix_pair;

static int size_mul(size_t a, size_t b, size_t *out) {
    if (a && b > SIZE_MAX/a) return 0;
    *out = a*b;
    return 1;
}
static int size_add(size_t a, size_t b, size_t *out) {
    if (b > SIZE_MAX-a) return 0;
    *out=a+b;
    return 1;
}
#ifdef LIS_TEST_OOM
static long lis_sc_alloc_fail_index = -1;
static long lis_sc_alloc_counter = 0;
void lis_sc_test_fail_after(long index){
    lis_sc_alloc_fail_index=index;
    lis_sc_alloc_counter=0;
}
#endif
static void *calloc_n(size_t n, size_t item) {
    size_t bytes;
    if (!size_mul(n,item,&bytes)) return NULL;
#ifdef LIS_TEST_OOM
    if(lis_sc_alloc_fail_index>=0 &&
       lis_sc_alloc_counter++==lis_sc_alloc_fail_index)return NULL;
#endif
    return calloc(n ? n : 1, item);
}
static int compare_link(const void *ap, const void *bp) {
    const link_pair *a=ap, *b=bp;
    if(a->block!=b->block) return (a->block>b->block)-(a->block<b->block);
    return (a->phys>b->phys)-(a->phys<b->phys);
}
static int compare_pair(const void *ap, const void *bp) {
    const matrix_pair *a=ap, *b=bp;
    if(a->row!=b->row) return (a->row>b->row)-(a->row<b->row);
    return (a->col>b->col)-(a->col<b->col);
}

static lis_sc_status check_csr(const lis_sc_csr *a) {
    size_t i,k;
    if(!a || !a->n || a->n==SIZE_MAX || !a->rowptr ||
       (a->nnz && (!a->colind || !a->values)))
        return LIS_SC_BAD_ARGUMENT;
    if(a->rowptr[0]!=0 || a->rowptr[a->n]!=a->nnz) return LIS_SC_BAD_CSR;
    for(i=0;i<a->n;i++) {
        if(a->rowptr[i]>a->rowptr[i+1] || a->rowptr[i+1]>a->nnz)
            return LIS_SC_BAD_CSR;
        for(k=a->rowptr[i];k<a->rowptr[i+1];k++) {
            if(a->colind[k]>=a->n || (k>a->rowptr[i] && a->colind[k-1]>=a->colind[k]))
                return LIS_SC_BAD_CSR;
            if(!isfinite(a->values[k])) return LIS_SC_NONFINITE;
        }
    }
    return LIS_SC_OK;
}

static lis_sc_status copy_csr(lis_sc_csr *dst,const lis_sc_csr *src) {
    size_t i;
    size_t *rp=(size_t*)calloc_n(src->n+1,sizeof(size_t));
    size_t *ci=(size_t*)calloc_n(src->nnz,sizeof(size_t));
    double *va=(double*)calloc_n(src->nnz,sizeof(double));
    if(!rp||!ci||!va){free(rp);free(ci);free(va);return LIS_SC_OUT_OF_MEMORY;}
    for(i=0;i<=src->n;i++)rp[i]=src->rowptr[i];
    for(i=0;i<src->nnz;i++){ci[i]=src->colind[i];va[i]=src->values[i];}
    dst->n=src->n;dst->nnz=src->nnz;dst->rowptr=rp;dst->colind=ci;dst->values=va;
    return LIS_SC_OK;
}
static void free_csr(lis_sc_csr *x) {
    free((void*)x->rowptr);free((void*)x->colind);free((void*)x->values);
    memset(x,0,sizeof(*x));
}
void lis_sc_free(lis_sc_plan *p) {
    if(!p)return;
    free_csr(&p->a);free_csr(&p->s);
    free(p->u2full);free(p->full2u);free(p->block_of);free(p->local_of);
    free(p->block_ptr);free(p->block_dof);free(p->neighbors_ptr);free(p->neighbors);
    free(p->lu_ptr);free(p->piv_r);free(p->piv_c);free(p->lu);free(p);
}

static size_t csr_find(const lis_sc_csr *a,size_t row,size_t col) {
    size_t lo=a->rowptr[row],hi=a->rowptr[row+1];
    while(lo<hi){size_t mid=lo+(hi-lo)/2;if(a->colind[mid]<col)lo=mid+1;else hi=mid;}
    return (lo<a->rowptr[row+1] && a->colind[lo]==col)?lo:SIZE_MAX;
}

static lis_sc_status validate_blocks(lis_sc_plan *p,const lis_sc_blocks *b) {
    size_t i,j,next=0;
    if(!b)return LIS_SC_BAD_ARGUMENT;
    p->nb=b->nblocks;
    if(p->nb && (!b->block_ptr || !b->block_dof)) return LIS_SC_BAD_MAP;
    if(p->nb && b->block_ptr[0]!=0) return LIS_SC_BAD_MAP;
    if(p->nb){
        p->nq=b->block_ptr[p->nb];
        if(p->nq>=p->a.n || !p->nq) return LIS_SC_BAD_MAP;
        if(p->nb>p->nq) return LIS_SC_BAD_MAP;
    }
    p->nphys=p->a.n-p->nq;
    p->full2u=(size_t*)calloc_n(p->a.n,sizeof(size_t));
    p->block_of=(size_t*)calloc_n(p->a.n,sizeof(size_t));
    p->local_of=(size_t*)calloc_n(p->a.n,sizeof(size_t));
    p->u2full=(size_t*)calloc_n(p->nphys,sizeof(size_t));
    p->block_ptr=(size_t*)calloc_n(p->nb+1,sizeof(size_t));
    p->block_dof=(size_t*)calloc_n(p->nq,sizeof(size_t));
    if(!p->full2u||!p->block_of||!p->local_of||!p->u2full||!p->block_ptr||!p->block_dof)
        return LIS_SC_OUT_OF_MEMORY;
    for(i=0;i<p->a.n;i++){p->block_of[i]=SIZE_MAX;p->full2u[i]=SIZE_MAX;p->local_of[i]=SIZE_MAX;}
    if(p->nb){
        for(i=0;i<p->nb;i++){
            size_t start=b->block_ptr[i],end=b->block_ptr[i+1];
            if(start!=next||end<=start||end>p->nq) return LIS_SC_BAD_MAP;
            p->block_ptr[i]=start;
            for(j=start;j<end;j++){
                size_t d=b->block_dof[j];
                if(d>=p->a.n||p->block_of[d]!=SIZE_MAX) return LIS_SC_BAD_MAP;
                p->block_of[d]=i;p->local_of[d]=j-start;p->block_dof[j]=d;
            }
            next=end;
        }
        if(next!=p->nq) return LIS_SC_BAD_MAP;
    }
    p->block_ptr[p->nb]=p->nq;
    next=0;
    for(i=0;i<p->a.n;i++) if(p->block_of[i]==SIZE_MAX){p->full2u[i]=next;p->u2full[next++]=i;}
    return next==p->nphys?LIS_SC_OK:LIS_SC_BAD_MAP;
}

static lis_sc_status check_locality(const lis_sc_plan *p) {
    size_t i,k;
    for(i=0;i<p->a.n;i++) if(p->block_of[i]!=SIZE_MAX){
        for(k=p->a.rowptr[i];k<p->a.rowptr[i+1];k++){
            size_t j=p->a.colind[k];
            if(p->block_of[j]!=SIZE_MAX && p->block_of[i]!=p->block_of[j]
                    && p->a.values[k]!=0.0) return LIS_SC_COUPLED_BLOCKS;
        }
    }
    return LIS_SC_OK;
}

static lis_sc_status build_neighbors(lis_sc_plan *p) {
    size_t k,i,count=0,cap=p->a.nnz,unique=0,b;
    link_pair *ls=(link_pair*)calloc_n(cap,sizeof(link_pair));
    if(!ls)return LIS_SC_OUT_OF_MEMORY;
    for(i=0;i<p->a.n;i++){
        size_t bi=p->block_of[i];
        for(k=p->a.rowptr[i];k<p->a.rowptr[i+1];k++){
            size_t j=p->a.colind[k],bj=p->block_of[j];
            /* Neighbors are structural, not numerical: an explicitly stored
             * zero can become nonzero during a Newton update without a
             * sparsity-pattern change.  Keep its potential Schur fill. */
            if(bi!=SIZE_MAX && bj==SIZE_MAX){ls[count].block=bi;ls[count++].phys=p->full2u[j];}
            else if(bi==SIZE_MAX && bj!=SIZE_MAX){ls[count].block=bj;ls[count++].phys=p->full2u[i];}
        }
    }
    qsort(ls,count,sizeof(link_pair),compare_link);
    p->neighbors_ptr=(size_t*)calloc_n(p->nb+1,sizeof(size_t));
    p->neighbors=(size_t*)calloc_n(count,sizeof(size_t));
    if(!p->neighbors_ptr||!p->neighbors){free(ls);return LIS_SC_OUT_OF_MEMORY;}
    k=0;
    for(b=0;b<p->nb;b++){
        p->neighbors_ptr[b]=unique;
        while(k<count&&ls[k].block==b){
            if(unique==p->neighbors_ptr[b]||ls[k].phys!=p->neighbors[unique-1])
                p->neighbors[unique++]=ls[k].phys;
            k++;
        }
    }
    p->neighbors_ptr[p->nb]=unique;
    free(ls);
    return LIS_SC_OK;
}

static lis_sc_status build_symbolic(lis_sc_plan *p){
    size_t i,j,k,b,ntr=0,cap=p->nphys,uniq=0;
    matrix_pair *tr;
    size_t *rowptr=NULL,*colind=NULL;
    double *values=NULL;
    for(i=0;i<p->nphys;i++){
        size_t row=p->u2full[i];
        for(k=p->a.rowptr[row];k<p->a.rowptr[row+1];k++)
            if(p->full2u[p->a.colind[k]]!=SIZE_MAX){if(!size_add(cap,1,&cap))return LIS_SC_SIZE;}
    }
    for(b=0;b<p->nb;b++){
        size_t nu=p->neighbors_ptr[b+1]-p->neighbors_ptr[b],sq;
        if(!size_mul(nu,nu,&sq)||!size_add(cap,sq,&cap))return LIS_SC_SIZE;
    }
    tr=(matrix_pair*)calloc_n(cap,sizeof(matrix_pair));
    if(!tr)return LIS_SC_OUT_OF_MEMORY;
    for(i=0;i<p->nphys;i++){
        size_t row=p->u2full[i];
        tr[ntr++]=(matrix_pair){i,i}; /* structural diagonal */
        for(k=p->a.rowptr[row];k<p->a.rowptr[row+1];k++){
            j=p->full2u[p->a.colind[k]];
            if(j!=SIZE_MAX)tr[ntr++]=(matrix_pair){i,j};
        }
    }
    for(b=0;b<p->nb;b++){
        for(i=p->neighbors_ptr[b];i<p->neighbors_ptr[b+1];i++)
          for(j=p->neighbors_ptr[b];j<p->neighbors_ptr[b+1];j++)
            tr[ntr++]=(matrix_pair){p->neighbors[i],p->neighbors[j]};
    }
    qsort(tr,ntr,sizeof(matrix_pair),compare_pair);
    for(k=0;k<ntr;k++)
        if(!k||compare_pair(&tr[k],&tr[k-1])!=0) uniq++;
    rowptr=(size_t*)calloc_n(p->nphys+1,sizeof(size_t));
    colind=(size_t*)calloc_n(uniq,sizeof(size_t));
    values=(double*)calloc_n(uniq,sizeof(double));
    if(!rowptr||!colind||!values){free(tr);free(rowptr);free(colind);free(values);return LIS_SC_OUT_OF_MEMORY;}
    size_t at=0;
    for(k=0;k<ntr;k++) if(!k||compare_pair(&tr[k],&tr[k-1])!=0){
        size_t r=tr[k].row;
        colind[at++]=tr[k].col;
        rowptr[r+1]++;
    }
    for(i=0;i<p->nphys;i++)rowptr[i+1]+=rowptr[i];
    free(tr);
    p->s=(lis_sc_csr){p->nphys,uniq,rowptr,colind,values};
    return LIS_SC_OK;
}

/* Complete pivoting: PAQ=LU. pr and pc track original row/col indices. */
static lis_sc_status factor_gecp(size_t n,double *a,size_t *pr,size_t *pc,double tol){
    size_t i,j,k,pi,pj;
    double scale=0;
    for(i=0;i<n;i++)for(j=0;j<n;j++){
        double v=fabs(a[i*n+j]);if(!isfinite(v))return LIS_SC_NONFINITE;
        if(v>scale)scale=v;
    }
    if(scale==0.0)return LIS_SC_SINGULAR_BLOCK;
    for(i=0;i<n;i++){pr[i]=i;pc[i]=i;}
    for(k=0;k<n;k++){
        double best=0;pi=k;pj=k;
        for(i=k;i<n;i++)for(j=k;j<n;j++){
            double v=fabs(a[i*n+j]);
            if(v>best){best=v;pi=i;pj=j;}
        }
        if(!isfinite(best))return LIS_SC_NONFINITE;
        if(best==0.0||best/scale<=tol)return LIS_SC_SINGULAR_BLOCK;
        if(pi!=k){size_t t=pr[k];pr[k]=pr[pi];pr[pi]=t;
            for(j=0;j<n;j++){double x=a[k*n+j];a[k*n+j]=a[pi*n+j];a[pi*n+j]=x;}}
        if(pj!=k){size_t t=pc[k];pc[k]=pc[pj];pc[pj]=t;
            for(i=0;i<n;i++){double x=a[i*n+k];a[i*n+k]=a[i*n+pj];a[i*n+pj]=x;}}
        for(i=k+1;i<n;i++){
            double mult=a[i*n+k]/a[k*n+k];
            if(!isfinite(mult))return LIS_SC_NONFINITE;
            a[i*n+k]=mult;
            for(j=k+1;j<n;j++){
                a[i*n+j]-=mult*a[k*n+j];
                if(!isfinite(a[i*n+j]))return LIS_SC_NONFINITE;
            }
        }
    }
    return LIS_SC_OK;
}
static lis_sc_status solve_gecp(size_t n,const double *lu,const size_t *pr,
                                const size_t *pc,const double *rhs,double *out,double *work){
    size_t i,j;
    for(i=0;i<n;i++)work[i]=rhs[pr[i]];
    for(i=0;i<n;i++)for(j=0;j<i;j++)work[i]-=lu[i*n+j]*work[j];
    for(i=n;i>0;){--i;for(j=i+1;j<n;j++)work[i]-=lu[i*n+j]*work[j];
        work[i]/=lu[i*n+i];if(!isfinite(work[i]))return LIS_SC_NONFINITE;}
    for(i=0;i<n;i++)out[pc[i]]=work[i];
    return LIS_SC_OK;
}

static lis_sc_status build_numeric(lis_sc_plan *p){
    size_t b,i,j,k,nlu=0;
    lis_sc_status st;
    /* Copy Kuu into structural Schur. */
    for(i=0;i<p->nphys;i++){
        size_t full=p->u2full[i];
        for(k=p->a.rowptr[full];k<p->a.rowptr[full+1];k++){
            size_t fcol=p->full2u[p->a.colind[k]];
            if(fcol!=SIZE_MAX){size_t pos=csr_find(&p->s,i,fcol);
                if(pos==SIZE_MAX)return LIS_SC_BAD_CSR;
                ((double*)p->s.values)[pos]+=p->a.values[k];}
        }
    }
    p->lu_ptr=(size_t*)calloc_n(p->nb+1,sizeof(size_t));
    p->piv_r=(size_t*)calloc_n(p->nq,sizeof(size_t));
    p->piv_c=(size_t*)calloc_n(p->nq,sizeof(size_t));
    if(!p->lu_ptr||!p->piv_r||!p->piv_c)return LIS_SC_OUT_OF_MEMORY;
    for(b=0;b<p->nb;b++){
        size_t m=p->block_ptr[b+1]-p->block_ptr[b],msq;
        if(!size_mul(m,m,&msq)||!size_add(nlu,msq,&nlu))return LIS_SC_SIZE;
        p->lu_ptr[b+1]=nlu;
    }
    p->lu=(double*)calloc_n(nlu,sizeof(double));
    if(!p->lu)return LIS_SC_OUT_OF_MEMORY;
    for(b=0;b<p->nb;b++){
        size_t start=p->block_ptr[b],m=p->block_ptr[b+1]-start;
        size_t from=p->neighbors_ptr[b],nu=p->neighbors_ptr[b+1]-from;
        double *lu=p->lu+p->lu_ptr[b],*kuq=NULL,*kqu=NULL,*rhs=NULL,*sol=NULL,*work=NULL;
        size_t *pr=p->piv_r+start,*pc=p->piv_c+start;
        size_t entries;
        if(nu){
            if(!size_mul(m,nu,&entries))return LIS_SC_SIZE;
            kuq=(double*)calloc_n(entries,sizeof(double));
            kqu=(double*)calloc_n(entries,sizeof(double));
            rhs=(double*)calloc_n(m,sizeof(double));
            sol=(double*)calloc_n(m,sizeof(double));
            work=(double*)calloc_n(m,sizeof(double));
            if(!kuq||!kqu||!rhs||!sol||!work){st=LIS_SC_OUT_OF_MEMORY;goto fail_block;}
        }
        for(i=0;i<m;i++){
            size_t full=p->block_dof[start+i];
            for(k=p->a.rowptr[full];k<p->a.rowptr[full+1];k++){
                size_t col=p->a.colind[k];
                if(p->block_of[col]==b)lu[i*m+p->local_of[col]]=p->a.values[k];
                else if(p->full2u[col]!=SIZE_MAX){
                    size_t off=from, hi=from+nu, phys=p->full2u[col];
                    while(off<hi){size_t mid=off+(hi-off)/2;
                        if(p->neighbors[mid]<phys)off=mid+1;else hi=mid;}
                    if(off>=from+nu||p->neighbors[off]!=phys){st=LIS_SC_BAD_CSR;goto fail_block;}
                    kqu[i*nu+off-from]=p->a.values[k];
                }
            }
        }
        st=factor_gecp(m,lu,pr,pc,p->pivot_tol);
        if(st!=LIS_SC_OK)goto fail_block;
        /* Build Kuq independently, handling structurally nonsymmetric matrices. */
        for(i=0;i<nu;i++){
            size_t full=p->u2full[p->neighbors[from+i]];
            for(j=0;j<m;j++){
                size_t pos=csr_find(&p->a,full,p->block_dof[start+j]);
                if(pos!=SIZE_MAX)kuq[i*m+j]=p->a.values[pos];
            }
        }
        for(j=0;j<nu;j++){
            for(i=0;i<m;i++)rhs[i]=kqu[i*nu+j];
            st=solve_gecp(m,lu,pr,pc,rhs,sol,work);
            if(st!=LIS_SC_OK)goto fail_block;
            for(i=0;i<nu;i++){
                double dot=0;
                for(k=0;k<m;k++)dot+=kuq[i*m+k]*sol[k];
                if(!isfinite(dot)){st=LIS_SC_NONFINITE;goto fail_block;}
                size_t row=p->neighbors[from+i],col=p->neighbors[from+j];
                size_t pos=csr_find(&p->s,row,col);
                if(pos==SIZE_MAX){st=LIS_SC_BAD_CSR;goto fail_block;}
                ((double*)p->s.values)[pos]-=dot;
                if(!isfinite(p->s.values[pos])){st=LIS_SC_NONFINITE;goto fail_block;}
            }
        }
        st=LIS_SC_OK;
fail_block:
        free(kuq);free(kqu);free(rhs);free(sol);free(work);
        if(st!=LIS_SC_OK)return st;
    }
    return LIS_SC_OK;
}

/* Clone only the symbolic state for a same-pattern numerical update.
 * New numeric values, factorizations and Schur coefficients are built in a
 * fresh plan.  The old plan remains fully usable if any step fails. */
static size_t *copy_index_array(const size_t *src,size_t count){
    size_t bytes;
    size_t *dst;
    if(!size_mul(count,sizeof(size_t),&bytes))return NULL;
    dst=(size_t*)calloc_n(count,sizeof(size_t));
    if(dst && count)memcpy(dst,src,bytes);
    return dst;
}
static lis_sc_status clone_symbolic(lis_sc_plan *p,const lis_sc_plan *src){
    size_t nlinks;
    p->nphys=src->nphys;
    p->nq=src->nq;
    p->nb=src->nb;
    nlinks=src->neighbors_ptr[src->nb];
    p->u2full=copy_index_array(src->u2full,p->nphys);
    p->full2u=copy_index_array(src->full2u,p->a.n);
    p->block_of=copy_index_array(src->block_of,p->a.n);
    p->local_of=copy_index_array(src->local_of,p->a.n);
    p->block_ptr=copy_index_array(src->block_ptr,p->nb+1);
    p->block_dof=copy_index_array(src->block_dof,p->nq);
    p->neighbors_ptr=copy_index_array(src->neighbors_ptr,p->nb+1);
    p->neighbors=copy_index_array(src->neighbors,nlinks);
    p->s.n=src->s.n;
    p->s.nnz=src->s.nnz;
    p->s.rowptr=copy_index_array(src->s.rowptr,p->nphys+1);
    p->s.colind=copy_index_array(src->s.colind,src->s.nnz);
    p->s.values=(double*)calloc_n(src->s.nnz,sizeof(double));
    if(!p->u2full||!p->full2u||!p->block_of||!p->local_of||
       !p->block_ptr||!p->block_dof||!p->neighbors_ptr||!p->neighbors||
       !p->s.rowptr||!p->s.colind||!p->s.values)return LIS_SC_OUT_OF_MEMORY;
    return LIS_SC_OK;
}

static int same_pattern(const lis_sc_csr *a,const lis_sc_csr *b){
    return a->n==b->n && a->nnz==b->nnz &&
           memcmp(a->rowptr,b->rowptr,(a->n+1)*sizeof(size_t))==0 &&
           (!a->nnz||memcmp(a->colind,b->colind,a->nnz*sizeof(size_t))==0);
}
static int same_blocks(const lis_sc_plan *p,const lis_sc_blocks *b){
    if(!b||p->nb!=b->nblocks)return 0;
    if(!p->nb)return 1;
    if(!b->block_ptr||!b->block_dof)return 0;
    return memcmp(p->block_ptr,b->block_ptr,(p->nb+1)*sizeof(size_t))==0 &&
        memcmp(p->block_dof,b->block_dof,p->nq*sizeof(size_t))==0;
}

lis_sc_status lis_sc_stage(lis_sc_plan **candidate, const lis_sc_plan *base,
                           const lis_sc_csr *a, const lis_sc_blocks *blocks,
                           lis_sc_reuse_mode mode, double pivot_tolerance){
    lis_sc_plan *p;
    lis_sc_status st;
    if(candidate)*candidate=NULL;
    if(!candidate||!blocks||!(mode==LIS_SC_FULL||mode==LIS_SC_NUMERIC||mode==LIS_SC_REUSE_ALL))return LIS_SC_BAD_ARGUMENT;
    if(!isfinite(pivot_tolerance)||pivot_tolerance<0||pivot_tolerance>=1)return LIS_SC_BAD_ARGUMENT;
    st=check_csr(a);if(st!=LIS_SC_OK)return st;
    if(mode!=LIS_SC_FULL){
        if(!base||!same_pattern(&base->a,a)||!same_blocks(base,blocks)||
           base->pivot_tol!=pivot_tolerance) return LIS_SC_PATTERN_MISMATCH;
        if(mode==LIS_SC_REUSE_ALL){
            if(!a->nnz || memcmp(base->a.values,a->values,a->nnz*sizeof(double))==0)
                return LIS_SC_OK;
            return LIS_SC_PATTERN_MISMATCH;
        }
    }
    p=(lis_sc_plan*)calloc_n(1,sizeof(*p));
    if(!p)return LIS_SC_OUT_OF_MEMORY;
    p->pivot_tol=pivot_tolerance;
    st=copy_csr(&p->a,a);if(st!=LIS_SC_OK)goto fail;
    if(mode==LIS_SC_NUMERIC){
        /* The earlier symbolic state is valid because both the entire CSR
         * sparsity pattern and block metadata were checked above. */
        st=clone_symbolic(p,base);if(st!=LIS_SC_OK)goto fail;
        st=check_locality(p);if(st!=LIS_SC_OK)goto fail;
    }else{
        st=validate_blocks(p,blocks);if(st!=LIS_SC_OK)goto fail;
        st=check_locality(p);if(st!=LIS_SC_OK)goto fail;
        st=build_neighbors(p);if(st!=LIS_SC_OK)goto fail;
        st=build_symbolic(p);if(st!=LIS_SC_OK)goto fail;
    }
    st=build_numeric(p);if(st!=LIS_SC_OK)goto fail;
    *candidate=p;
    return LIS_SC_OK;
fail:
    lis_sc_free(p);
    return st;
}

/* Historical atomic replacement API, backed by the same transactional path. */
lis_sc_status lis_sc_build(lis_sc_plan **slot, const lis_sc_csr *a,
                           const lis_sc_blocks *blocks, lis_sc_reuse_mode mode,
                           double pivot_tolerance){
    lis_sc_plan *candidate=NULL;
    lis_sc_status st;
    if(!slot)return LIS_SC_BAD_ARGUMENT;
    st=lis_sc_stage(&candidate,*slot,a,blocks,mode,pivot_tolerance);
    if(st!=LIS_SC_OK)return st;
    if(candidate){
        lis_sc_plan *prior=*slot;
        *slot=candidate;
        lis_sc_free(prior);
    }
    return LIS_SC_OK;
}

const lis_sc_csr *lis_sc_matrix(const lis_sc_plan *p){return p?&p->s:NULL;}
size_t lis_sc_physical_dofs(const lis_sc_plan *p){return p?p->nphys:0;}
size_t lis_sc_internal_dofs(const lis_sc_plan *p){return p?p->nq:0;}
const size_t *lis_sc_physical_to_full(const lis_sc_plan *p){return p?p->u2full:NULL;}

lis_sc_status lis_sc_rhs(const lis_sc_plan *p,const double *b,double *g){
    size_t i,k,bl;
    double *zq=NULL,*tmp=NULL;
    if(!p||!b||!g)return LIS_SC_BAD_ARGUMENT;
    for(i=0;i<p->a.n;i++)if(!isfinite(b[i]))return LIS_SC_NONFINITE;
    zq=(double*)calloc_n(p->nq,sizeof(double));
    tmp=(double*)calloc_n(p->nphys,sizeof(double));
    if(!zq||!tmp){free(zq);free(tmp);return LIS_SC_OUT_OF_MEMORY;}
    for(bl=0;bl<p->nb;bl++){
        size_t start=p->block_ptr[bl],m=p->block_ptr[bl+1]-start;
        double *rhs=(double*)calloc_n(m,sizeof(double)),*wrk=(double*)calloc_n(m,sizeof(double));
        lis_sc_status st;
        if(!rhs||!wrk){free(rhs);free(wrk);free(zq);free(tmp);return LIS_SC_OUT_OF_MEMORY;}
        for(i=0;i<m;i++)rhs[i]=b[p->block_dof[start+i]];
        st=solve_gecp(m,p->lu+p->lu_ptr[bl],p->piv_r+start,p->piv_c+start,rhs,zq+start,wrk);
        free(rhs);free(wrk);
        if(st!=LIS_SC_OK){free(zq);free(tmp);return st;}
    }
    for(i=0;i<p->nphys;i++){
        size_t full=p->u2full[i];
        double val=b[full];
        for(k=p->a.rowptr[full];k<p->a.rowptr[full+1];k++){
            size_t col=p->a.colind[k],bl2=p->block_of[col];
            if(bl2!=SIZE_MAX)val-=p->a.values[k]*zq[p->block_ptr[bl2]+p->local_of[col]];
        }
        if(!isfinite(val)){free(zq);free(tmp);return LIS_SC_NONFINITE;}
        tmp[i]=val;
    }
    memcpy(g,tmp,p->nphys*sizeof(double));
    free(zq);free(tmp);
    return LIS_SC_OK;
}

lis_sc_status lis_sc_reconstruct(const lis_sc_plan *p,const double *u,
                                  const double *b,double *out){
    size_t bl,i,k;
    double *tmp;
    if(!p||!u||!b||!out)return LIS_SC_BAD_ARGUMENT;
    for(i=0;i<p->nphys;i++)if(!isfinite(u[i]))return LIS_SC_NONFINITE;
    for(i=0;i<p->a.n;i++)if(!isfinite(b[i]))return LIS_SC_NONFINITE;
    tmp=(double*)calloc_n(p->a.n,sizeof(double));
    if(!tmp)return LIS_SC_OUT_OF_MEMORY;
    for(i=0;i<p->nphys;i++)tmp[p->u2full[i]]=u[i];
    for(bl=0;bl<p->nb;bl++){
        size_t start=p->block_ptr[bl],m=p->block_ptr[bl+1]-start;
        double *rhs=(double*)calloc_n(m,sizeof(double)),*wrk=(double*)calloc_n(m,sizeof(double)),*sol=(double*)calloc_n(m,sizeof(double));
        lis_sc_status st;
        if(!rhs||!wrk||!sol){free(rhs);free(wrk);free(sol);free(tmp);return LIS_SC_OUT_OF_MEMORY;}
        for(i=0;i<m;i++){
            size_t full=p->block_dof[start+i];
            double val=b[full];
            for(k=p->a.rowptr[full];k<p->a.rowptr[full+1];k++){
                size_t ucol=p->full2u[p->a.colind[k]];
                if(ucol!=SIZE_MAX)val-=p->a.values[k]*u[ucol];
            }
            rhs[i]=val;
        }
        st=solve_gecp(m,p->lu+p->lu_ptr[bl],p->piv_r+start,p->piv_c+start,rhs,sol,wrk);
        if(st==LIS_SC_OK)for(i=0;i<m;i++)tmp[p->block_dof[start+i]]=sol[i];
        free(rhs);free(wrk);free(sol);
        if(st!=LIS_SC_OK){free(tmp);return st;}
    }
    memcpy(out,tmp,p->a.n*sizeof(double));
    free(tmp);
    return LIS_SC_OK;
}

lis_sc_status lis_sc_true_residual(const lis_sc_plan *p,const double *b,
                                    const double *x,double *rr){
    size_t i,k;
    double bnorm=0,rnorm=0;
    if(!p||!b||!x||!rr)return LIS_SC_BAD_ARGUMENT;
    for(i=0;i<p->a.n;i++)if(!isfinite(x[i])||!isfinite(b[i]))return LIS_SC_NONFINITE;
    for(i=0;i<p->a.n;i++){
        long double dot=0;
        for(k=p->a.rowptr[i];k<p->a.rowptr[i+1];k++){
            size_t col=p->a.colind[k];
            if(!isfinite(x[col]))return LIS_SC_NONFINITE;
            dot+=(long double)p->a.values[k]*(long double)x[col];
        }
        if(!isfinite(b[i])||!isfinite((double)dot))return LIS_SC_NONFINITE;
        bnorm=hypot(bnorm,b[i]);
        rnorm=hypot(rnorm,(double)(dot-(long double)b[i]));
    }
    *rr=bnorm ? rnorm/bnorm : rnorm;
    return isfinite(*rr)?LIS_SC_OK:LIS_SC_NONFINITE;
}

lis_sc_status lis_sc_project_modes(const lis_sc_plan *p,size_t count,
                                    const double *full,double *phys){
    size_t i,k,total;
    if(!p || (count && (!full||!phys)))return LIS_SC_BAD_ARGUMENT;
    if(!size_mul(count,p->a.n,&total) || !size_mul(count,p->nphys,&total))
        return LIS_SC_SIZE;
    /* Reject invalid input before modifying any of the output modes. */
    for(k=0;k<count;k++)for(i=0;i<p->nphys;i++)
        if(!isfinite(full[k*p->a.n+p->u2full[i]]))return LIS_SC_NONFINITE;
    for(k=0;k<count;k++)for(i=0;i<p->nphys;i++)
        phys[k*p->nphys+i]=full[k*p->a.n+p->u2full[i]];
    return LIS_SC_OK;
}
lis_sc_status lis_sc_project_labels(const lis_sc_plan *p,const long *src,long *dst){
    size_t i;
    if(!p||!src||!dst)return LIS_SC_BAD_ARGUMENT;
    for(i=0;i<p->nphys;i++)dst[i]=src[p->u2full[i]];
    return LIS_SC_OK;
}
const char *lis_sc_status_string(lis_sc_status st){
    switch(st){
    case LIS_SC_OK:return "OK";case LIS_SC_BAD_ARGUMENT:return "BAD_ARGUMENT";
    case LIS_SC_BAD_CSR:return "BAD_CSR";case LIS_SC_BAD_MAP:return "BAD_MAP";
    case LIS_SC_COUPLED_BLOCKS:return "COUPLED_BLOCKS";
    case LIS_SC_SINGULAR_BLOCK:return "SINGULAR_BLOCK";
    case LIS_SC_NONFINITE:return "NONFINITE";case LIS_SC_SIZE:return "SIZE";
    case LIS_SC_OUT_OF_MEMORY:return "OUT_OF_MEMORY";
    case LIS_SC_PATTERN_MISMATCH:return "PATTERN_MISMATCH";
    case LIS_SC_UNSUPPORTED:return "UNSUPPORTED";
    case LIS_SC_SOLVER_ERROR:return "SOLVER_ERROR";
    case LIS_SC_RESIDUAL_REJECTED:return "RESIDUAL_REJECTED";
    case LIS_SC_REBIND_REQUIRED:return "REBIND_REQUIRED";
    default:return "UNKNOWN";
    }
}
