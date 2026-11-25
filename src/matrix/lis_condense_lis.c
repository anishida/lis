#include "lislib.h"
#include "lis_condense_lis.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

struct lis_sc_lis_session {
    lis_sc_plan *plan;
    LIS_MATRIX S;
};

void lis_sc_lis_destroy(lis_sc_lis_session *s){
    if(!s)return;
    if(s->S)lis_matrix_destroy(s->S);
    lis_sc_free(s->plan);
    free(s);
}
LIS_MATRIX lis_sc_lis_matrix(const lis_sc_lis_session *s){return s?s->S:NULL;}

#if defined(USE_COMPLEX) || defined(_COMPLEX) || defined(USE_MPI)
/* Do not silently reinterpret complex scalars / distributed ghost indices. */
lis_sc_status lis_sc_lis_create(LIS_MATRIX A,LIS_INT nblocks,
                                const LIS_INT *block_ptr,const LIS_INT *block_dof,
                                double tol,lis_sc_lis_session **out){
    (void)A;(void)nblocks;(void)block_ptr;(void)block_dof;(void)tol;(void)out;
    return LIS_SC_UNSUPPORTED;
}
lis_sc_status lis_sc_lis_update(lis_sc_lis_session *session, LIS_MATRIX A,
                                 LIS_INT nblocks, const LIS_INT *block_ptr,
                                 const LIS_INT *block_dof, double pivot_tolerance,
                                 lis_sc_reuse_mode mode){
    (void)session;(void)A;(void)nblocks;(void)block_ptr;(void)block_dof;
    (void)pivot_tolerance;(void)mode;return LIS_SC_UNSUPPORTED;
}
lis_sc_status lis_sc_lis_rhs(const lis_sc_lis_session *s,LIS_VECTOR b,LIS_VECTOR g){
    (void)s;(void)b;(void)g;return LIS_SC_UNSUPPORTED;
}
lis_sc_status lis_sc_lis_reconstruct(const lis_sc_lis_session*s,LIS_VECTOR u,LIS_VECTOR b,LIS_VECTOR x){
    (void)s;(void)u;(void)b;(void)x;return LIS_SC_UNSUPPORTED;
}
lis_sc_status lis_sc_lis_residual(const lis_sc_lis_session*s,LIS_VECTOR b,LIS_VECTOR x,double *rr){
    (void)s;(void)b;(void)x;(void)rr;return LIS_SC_UNSUPPORTED;
}
lis_sc_status lis_sc_lis_labels(const lis_sc_lis_session*s,const LIS_INT *a,
                                 size_t an,LIS_INT *b,size_t bn){
    (void)s;(void)a;(void)an;(void)b;(void)bn;return LIS_SC_UNSUPPORTED;
}
lis_sc_status lis_sc_lis_modes(const lis_sc_lis_session*s,LIS_INT nvec,
                                LIS_VECTOR const full_modes[],
                                LIS_VECTOR const reduced_modes[]){
    (void)s;(void)nvec;(void)full_modes;(void)reduced_modes;
    return LIS_SC_UNSUPPORTED;
}
#else

/* This allocator is needed only for the supported serial real-double path. */
static void *sc_lis_calloc(size_t count, size_t elem){
    if (elem && count > SIZE_MAX/elem) return NULL;
    return calloc(count ? count : 1, elem);
}

static int valid_lis_int(size_t v){return (LIS_INT)v >= 0 && (size_t)(LIS_INT)v == v;}

/* Marshal raw LIS CSR to a fully owned portable view. This is required
 * because LIS_INT may have a different width than size_t. No external numeric
 * library is used and the caller's A is never modified. */
/* LIS RCO->CSR assembly preserves insertion order (see
 * lis_matrix_convert_rco2csr); column indices are not necessarily sorted.
 * Normalize a PRIVATE copy of the input, never the caller's matrix.
 * Duplicate columns are combined by summation, matching CSR matvec meaning.
 * The canonical pattern makes FULL/NUMERIC/REUSE_ALL independent of insertion
 * order while preserving every algebraic cross-partition coupling. */
typedef struct { size_t col; double value; } sc_lis_entry;
static int sc_lis_entry_compare(const void *ap,const void *bp){
    const sc_lis_entry *a=(const sc_lis_entry*)ap,*b=(const sc_lis_entry*)bp;
    return (a->col>b->col)-(a->col<b->col);
}
static lis_sc_status read_lis_csr(LIS_MATRIX A, lis_sc_csr *view){
    size_t i,k,n,nnz,maxrow=0,count=0;
    size_t *rp=NULL,*ci=NULL;
    double *v=NULL;
    sc_lis_entry *scratch=NULL;
    lis_sc_status st=LIS_SC_OK;
    if(!A||!view||A->n<=0||A->matrix_type!=LIS_MATRIX_CSR||
       !A->ptr||sizeof(LIS_SCALAR)!=sizeof(double))return LIS_SC_BAD_ARGUMENT;
    /* No assumptions about MPI ownership or split storage. */
    if(A->nprocs!=1 || A->gn!=A->n || A->is_splited)
        return LIS_SC_UNSUPPORTED;
    n=(size_t)A->n;
    if(n==SIZE_MAX||A->ptr[0]!=0||A->ptr[n]<0)return LIS_SC_BAD_CSR;
    nnz=(size_t)A->ptr[n];
    if(nnz&&(!A->index||!A->value))return LIS_SC_BAD_CSR;
    for(i=0;i<n;i++){
        size_t length;
        if(A->ptr[i]<0||A->ptr[i+1]<A->ptr[i]||A->ptr[i+1]>A->ptr[n])
            return LIS_SC_BAD_CSR;
        length=(size_t)(A->ptr[i+1]-A->ptr[i]);
        if(length>maxrow)maxrow=length;
    }
    rp=(size_t*)sc_lis_calloc(n+1,sizeof(*rp));
    ci=(size_t*)sc_lis_calloc(nnz,sizeof(*ci));
    v=(double*)sc_lis_calloc(nnz,sizeof(*v));
    scratch=(sc_lis_entry*)sc_lis_calloc(maxrow,sizeof(*scratch));
    if(!rp||!ci||!v||!scratch){st=LIS_SC_OUT_OF_MEMORY;goto fail;}
    for(i=0;i<n;i++){
        size_t length=(size_t)(A->ptr[i+1]-A->ptr[i]);
        for(k=0;k<length;k++){
            size_t pos=(size_t)A->ptr[i]+k;
            if(A->index[pos]<0 || (size_t)A->index[pos]>=n){
                st=LIS_SC_BAD_CSR;goto fail;
            }
            if(!isfinite((double)A->value[pos])){
                st=LIS_SC_NONFINITE;goto fail;
            }
            scratch[k].col=(size_t)A->index[pos];
            scratch[k].value=(double)A->value[pos];
        }
        qsort(scratch,length,sizeof(*scratch),sc_lis_entry_compare);
        for(k=0;k<length;k++){
            if(count>rp[i] && ci[count-1]==scratch[k].col){
                v[count-1]+=scratch[k].value;
                if(!isfinite(v[count-1])){st=LIS_SC_NONFINITE;goto fail;}
            }else{
                ci[count]=scratch[k].col;
                v[count]=scratch[k].value;
                count++;
            }
        }
        rp[i+1]=count;
    }
    free(scratch);
    *view=(lis_sc_csr){n,count,rp,ci,v};
    return LIS_SC_OK;
fail:
    free(scratch);free(rp);free(ci);free(v);
    return st;
}

static void free_lis_csr(lis_sc_csr *view){
    free((void*)view->rowptr);free((void*)view->colind);
    free((void*)view->values);memset(view,0,sizeof(*view));
}
static lis_sc_status read_lis_blocks(LIS_INT nblocks,
                                     const LIS_INT *block_ptr,
                                     const LIS_INT *block_dof,
                                     size_t n_total,
                                     lis_sc_blocks *blocks){
    size_t nb,nq,i;
    size_t *bp=NULL,*qd=NULL;
    if(!blocks||nblocks<0)return LIS_SC_BAD_ARGUMENT;
    nb=(size_t)nblocks;
    /* Bounds check BEFORE indexing block_ptr[nb]. One internal block needs
     * at least one DOF; at least one physical DOF must survive. */
    if(nb>=n_total || nb==SIZE_MAX ||
       (nb && (!block_ptr||!block_dof)))return LIS_SC_BAD_MAP;
    if(nb && block_ptr[0]!=0)return LIS_SC_BAD_MAP;
    nq=nb?(block_ptr[nb]<0?SIZE_MAX:(size_t)block_ptr[nb]):0;
    if(nq>=n_total || (nb && nq<nb))return LIS_SC_BAD_MAP;
    /* Validate all offsets prior to using the caller's DOF indices. */
    for(i=0;i<nb;i++)
        if(block_ptr[i]<0 || block_ptr[i+1]<=block_ptr[i] ||
           (size_t)block_ptr[i+1]>nq)return LIS_SC_BAD_MAP;
    bp=(size_t*)sc_lis_calloc(nb+1,sizeof(*bp));
    qd=(size_t*)sc_lis_calloc(nq,sizeof(*qd));
    if(!bp||!qd){free(bp);free(qd);return LIS_SC_OUT_OF_MEMORY;}
    for(i=0;i<=nb;i++){
        if(nb && block_ptr[i]<0){free(bp);free(qd);return LIS_SC_BAD_MAP;}
        bp[i]=nb?(size_t)block_ptr[i]:0;
    }
    for(i=0;i<nq;i++){
        if(block_dof[i]<0){free(bp);free(qd);return LIS_SC_BAD_MAP;}
        qd[i]=(size_t)block_dof[i];
    }
    *blocks=(lis_sc_blocks){nb,bp,qd};
    return LIS_SC_OK;
}
static void free_lis_blocks(lis_sc_blocks *b){
    free((void*)b->block_ptr);free((void*)b->block_dof);
    memset(b,0,sizeof(*b));
}

/* Materialize a candidate Schur matrix independently of the live session.
 * On successful LIS set_csr, LIS owns all three arrays. */
static lis_sc_status materialize_schur(const lis_sc_plan *plan,
                                       LIS_MATRIX *matrix){
    const lis_sc_csr *sv=lis_sc_matrix(plan);
    LIS_MATRIX fresh=NULL;
    LIS_INT *rp=NULL,*ci=NULL;
    LIS_SCALAR *values=NULL;
    LIS_INT err;
    size_t i;
    if(matrix)*matrix=NULL;
    if(!matrix||!sv||!valid_lis_int(sv->n)||!valid_lis_int(sv->nnz)||
       sv->n==SIZE_MAX)return LIS_SC_SIZE;
    err=lis_matrix_create(LIS_COMM_WORLD,&fresh);
    if(err)return LIS_SC_OUT_OF_MEMORY;
    err=lis_matrix_set_size(fresh,0,(LIS_INT)sv->n);
    if(err){lis_matrix_destroy(fresh);return LIS_SC_BAD_ARGUMENT;}
    rp=(LIS_INT*)sc_lis_calloc(sv->n+1,sizeof(*rp));
    ci=(LIS_INT*)sc_lis_calloc(sv->nnz,sizeof(*ci));
    values=(LIS_SCALAR*)sc_lis_calloc(sv->nnz,sizeof(*values));
    if(!rp||!ci||!values){free(rp);free(ci);free(values);lis_matrix_destroy(fresh);return LIS_SC_OUT_OF_MEMORY;}
    for(i=0;i<=sv->n;i++){
        if(!valid_lis_int(sv->rowptr[i]))goto size_fail;
        rp[i]=(LIS_INT)sv->rowptr[i];
    }
    for(i=0;i<sv->nnz;i++){
        if(!valid_lis_int(sv->colind[i]))goto size_fail;
        ci[i]=(LIS_INT)sv->colind[i];values[i]=sv->values[i];
    }
    err=lis_matrix_set_csr((LIS_INT)sv->nnz,rp,ci,values,fresh);
    if(err){free(rp);free(ci);free(values);lis_matrix_destroy(fresh);return LIS_SC_BAD_ARGUMENT;}
    /* Ownership transferred to LIS. */
    err=lis_matrix_assemble(fresh);
    if(err){lis_matrix_destroy(fresh);return LIS_SC_BAD_ARGUMENT;}
    *matrix=fresh;
    return LIS_SC_OK;
size_fail:
    free(rp);free(ci);free(values);lis_matrix_destroy(fresh);
    return LIS_SC_SIZE;
}

lis_sc_status lis_sc_lis_update(lis_sc_lis_session *session,LIS_MATRIX A,
                                 LIS_INT nblocks,const LIS_INT *block_ptr,
                                 const LIS_INT *block_dof,double pivot_tolerance,
                                 lis_sc_reuse_mode mode){
    lis_sc_csr av={0};
    lis_sc_blocks bl={0};
    lis_sc_plan *candidate=NULL,*oldplan;
    LIS_MATRIX fresh=NULL,oldmat;
    lis_sc_status st;
    if(!session)return LIS_SC_BAD_ARGUMENT;
    st=read_lis_csr(A,&av);if(st!=LIS_SC_OK)return st;
    st=read_lis_blocks(nblocks,block_ptr,block_dof,av.n,&bl);
    if(st!=LIS_SC_OK)goto cleanup;
    st=lis_sc_stage(&candidate,session->plan,&av,&bl,mode,pivot_tolerance);
    if(st!=LIS_SC_OK||!candidate)goto cleanup; /* REUSE_ALL is a no-op. */
    st=materialize_schur(candidate,&fresh);
    if(st!=LIS_SC_OK)goto cleanup;
    /* Commit point: no allocation or numerical operation can fail beyond here. */
    oldplan=session->plan;
    oldmat=session->S;
    session->plan=candidate;candidate=NULL;
    session->S=fresh;fresh=NULL;
    if(oldmat)lis_matrix_destroy(oldmat);
    lis_sc_free(oldplan);
cleanup:
    if(fresh)lis_matrix_destroy(fresh);
    lis_sc_free(candidate);
    free_lis_blocks(&bl);
    free_lis_csr(&av);
    return st;
}

lis_sc_status lis_sc_lis_create(LIS_MATRIX A,LIS_INT nblocks,
                                 const LIS_INT *block_ptr,const LIS_INT *block_dof,
                                 double pivot_tolerance,lis_sc_lis_session **out){
    lis_sc_lis_session *s;
    lis_sc_status st;
    if(out)*out=NULL;
    if(!out)return LIS_SC_BAD_ARGUMENT;
    s=(lis_sc_lis_session*)sc_lis_calloc(1,sizeof(*s));
    if(!s)return LIS_SC_OUT_OF_MEMORY;
    st=lis_sc_lis_update(s,A,nblocks,block_ptr,block_dof,pivot_tolerance,LIS_SC_FULL);
    if(st!=LIS_SC_OK){lis_sc_lis_destroy(s);return st;}
    *out=s;
    return LIS_SC_OK;
}

static int dims(const lis_sc_lis_session *s,LIS_VECTOR v,size_t want){
    return s && v && v->n>=0 && (size_t)v->n==want && v->value!=NULL;
}
lis_sc_status lis_sc_lis_rhs(const lis_sc_lis_session *s,LIS_VECTOR b,LIS_VECTOR g){
    if(!s||!s->plan||!dims(s,b,lis_sc_physical_dofs(s->plan)+lis_sc_internal_dofs(s->plan))||
       !dims(s,g,lis_sc_physical_dofs(s->plan)))return LIS_SC_BAD_ARGUMENT;
    return lis_sc_rhs(s->plan,b->value,g->value);
}
lis_sc_status lis_sc_lis_reconstruct(const lis_sc_lis_session*s,LIS_VECTOR u,LIS_VECTOR b,LIS_VECTOR x){
    if(!s||!s->plan||!dims(s,u,lis_sc_physical_dofs(s->plan))||
       !dims(s,b,lis_sc_physical_dofs(s->plan)+lis_sc_internal_dofs(s->plan))||
       !dims(s,x,lis_sc_physical_dofs(s->plan)+lis_sc_internal_dofs(s->plan)))return LIS_SC_BAD_ARGUMENT;
    return lis_sc_reconstruct(s->plan,u->value,b->value,x->value);
}
lis_sc_status lis_sc_lis_residual(const lis_sc_lis_session*s,LIS_VECTOR b,LIS_VECTOR x,double *rr){
    if(!s||!s->plan||!rr||!dims(s,b,lis_sc_physical_dofs(s->plan)+lis_sc_internal_dofs(s->plan))||
       !dims(s,x,lis_sc_physical_dofs(s->plan)+lis_sc_internal_dofs(s->plan)))return LIS_SC_BAD_ARGUMENT;
    return lis_sc_true_residual(s->plan,b->value,x->value,rr);
}
lis_sc_status lis_sc_lis_labels(const lis_sc_lis_session*s,const LIS_INT *src,
                                 size_t src_count,LIS_INT *out,size_t out_count){
    size_t i,nphys,ninternal;
    const size_t *u2full;
    if(!s||!s->plan||!src||!out)return LIS_SC_BAD_ARGUMENT;
    nphys=lis_sc_physical_dofs(s->plan);
    ninternal=lis_sc_internal_dofs(s->plan);
    if(nphys>SIZE_MAX-ninternal || src_count<nphys+ninternal ||
       out_count<nphys)return LIS_SC_BAD_ARGUMENT;
    u2full=lis_sc_physical_to_full(s->plan);
    if(!u2full)return LIS_SC_BAD_ARGUMENT;
    /* Validate every index before changing the caller's output. */
    for(i=0;i<nphys;i++)
        if(u2full[i]>=src_count)return LIS_SC_BAD_ARGUMENT;
    for(i=0;i<nphys;i++)out[i]=src[u2full[i]];
    return LIS_SC_OK;
}

lis_sc_status lis_sc_lis_modes(const lis_sc_lis_session *s,LIS_INT nvec,
                                LIS_VECTOR const full_modes[],
                                LIS_VECTOR const reduced_modes[]){
    size_t k,i,nfull,nphys;
    const size_t *u2full;
    if(!s||!s->plan||nvec<0)return LIS_SC_BAD_ARGUMENT;
    if(nvec==0)return LIS_SC_OK;
    if(!full_modes||!reduced_modes)return LIS_SC_BAD_ARGUMENT;
    nphys=lis_sc_physical_dofs(s->plan);
    nfull=nphys+lis_sc_internal_dofs(s->plan);
    u2full=lis_sc_physical_to_full(s->plan);
    /* Two-phase validation prevents partial updates on bad dimensions / NaNs. */
    for(k=0;k<(size_t)nvec;k++){
        if(!dims(s,full_modes[k],nfull)||!dims(s,reduced_modes[k],nphys))
            return LIS_SC_BAD_ARGUMENT;
        for(i=0;i<nphys;i++)
            if(!isfinite((double)full_modes[k]->value[u2full[i]]))
                return LIS_SC_NONFINITE;
    }
    for(k=0;k<(size_t)nvec;k++)for(i=0;i<nphys;i++)
        reduced_modes[k]->value[i]=full_modes[k]->value[u2full[i]];
    return LIS_SC_OK;
}
#endif
