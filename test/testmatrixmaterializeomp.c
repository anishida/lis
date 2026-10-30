/* OpenMP CSR materialization: safe USER/OPERATOR, fallback and failures. */
#include "lis.h"
#include "lis_matrix.h"
#include "lis_matrix_materialize.h"
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define NROWS 768
#define TRACK_THREADS 256

typedef struct {
    LIS_INT n;
    LIS_REAL scale;
    LIS_INT failed_row;
    int seen[TRACK_THREADS];
} CONTEXT;

static void scalar_real(LIS_SCALAR *p,LIS_REAL x)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
    (*p)[0]=x;
    (*p)[1]=(LIS_REAL)0;
#else
    *p=(LIS_SCALAR)x;
#endif
}

static int equal_real(const LIS_SCALAR *p,LIS_REAL x)
{
    LIS_SCALAR expected;
    scalar_real(&expected,x);
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
    return (*p)[0]==expected[0] && (*p)[1]==expected[1];
#else
    return *p==expected;
#endif
}

static LIS_INT row_callback(
    void *data,LIS_INT i,LIS_INT *nnz,LIS_INT *cols,LIS_SCALAR *vals)
{
    CONTEXT *c=(CONTEXT *)data;
    LIS_INT count,k=0;
    int tid=0;
#ifdef _OPENMP
    tid=omp_get_thread_num();
#endif
    if(c==NULL || nnz==NULL || i<0 || i>=c->n)
        return LIS_ERR_ILL_ARG;
    if(tid>=0 && tid<TRACK_THREADS)
        c->seen[tid]=1;
    if(i==c->failed_row)
        return LIS_ERR_ILL_ARG;

    count=1+(i>0)+(i<c->n-1);
    if(cols==NULL && vals==NULL)
    {
        *nnz=count;
        return LIS_SUCCESS;
    }
    if(cols==NULL || vals==NULL || *nnz<count)
    {
        *nnz=count;
        return LIS_ERR_ILL_ARG;
    }
    if(i>0)
    {
        cols[k]=i-1; scalar_real(&vals[k++],-c->scale);
    }
    cols[k]=i; scalar_real(&vals[k++],(LIS_REAL)4*c->scale);
    if(i<c->n-1)
    {
        cols[k]=i+1; scalar_real(&vals[k++],-c->scale);
    }
    *nnz=count;
    return LIS_SUCCESS;
}

static LIS_INT matvec(
    void *data,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CONTEXT *c=(CONTEXT *)data;
    LIS_INT i;
    if(!c || !x || !y) return LIS_ERR_ILL_ARG;
    for(i=0;i<c->n;i++)
    {
        y[i]=(LIS_SCALAR)4*c->scale*x[i];
        if(i>0)y[i]-=c->scale*x[i-1];
        if(i<c->n-1)y[i]-=c->scale*x[i+1];
    }
    return LIS_SUCCESS;
}

static LIS_INT make_user(CONTEXT *c,LIS_INT safe,LIS_MATRIX *out)
{
    LIS_INT err;
    *out=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,out);
    if(err) return err;
    err=lis_matrix_set_size(*out,0,c->n);
    if(!err)err=lis_matrix_set_user(*out,c,matvec,matvec);
    if(!err)err=lis_matrix_set_user_row(*out,row_callback);
    if(!err && safe)err=lis_matrix_set_user_row_thread_safe(*out,LIS_TRUE);
    if(err){lis_matrix_destroy(*out);*out=NULL;}
    return err;
}

static int check_csr(LIS_MATRIX A,LIS_REAL scale)
{
    LIS_INT i,p,k;
    if(A==NULL || A->matrix_type!=LIS_MATRIX_CSR)
        return 1;
    for(i=0;i<NROWS;i++)
    {
        LIS_INT n=1+(i>0)+(i<NROWS-1);
        if(A->ptr[i+1]-A->ptr[i]!=n)return 1;
        p=A->ptr[i]; k=0;
        if(i>0)
        {
            if(A->index[p+k]!=i-1 ||
               !equal_real(&A->value[p+k],-scale))return 1;
            k++;
        }
        if(A->index[p+k]!=i ||
           !equal_real(&A->value[p+k],(LIS_REAL)4*scale))return 1;
        k++;
        if(i<NROWS-1)
        {
            if(A->index[p+k]!=i+1 ||
               !equal_real(&A->value[p+k],-scale))return 1;
        }
    }
    return 0;
}

static int thread_count(const CONTEXT *c)
{
    int i,count=0;
    for(i=0;i<TRACK_THREADS;i++)count+=c->seen[i]!=0;
    return count;
}

int main(int argc,char **argv)
{
    CONTEXT a,b,unsafe;
    LIS_MATRIX A=NULL,B=NULL,U=NULL,op=NULL,csr=NULL,prior=NULL;
    LIS_MATRIX_CSR_CACHE ca=NULL,cb=NULL,cu=NULL;
    LIS_INT action,err;
    int failed=1;
    memset(&a,0,sizeof(a));memset(&b,0,sizeof(b));
    memset(&unsafe,0,sizeof(unsafe));
    a.n=b.n=unsafe.n=NROWS;
    a.scale=(LIS_REAL)1.0;
    b.scale=(LIS_REAL)0.5;
    unsafe.scale=(LIS_REAL)1.0;
    a.failed_row=b.failed_row=unsafe.failed_row=-1;
    err=lis_initialize(&argc,&argv);
    if(err)return 1;
    err=make_user(&a,LIS_TRUE,&A);if(err)goto done;
    err=make_user(&b,LIS_TRUE,&B);if(err)goto done;
    err=make_user(&unsafe,LIS_FALSE,&U);if(err)goto done;

    /* First-pass callback failure must leave no partially built cache. */
    a.failed_row=14;
    err=lis_matrix_materialize_csr_cached(A,&ca,&csr,&action);
    if(err==LIS_SUCCESS || csr!=NULL || ca!=NULL)goto done;
    a.failed_row=-1;

    err=lis_matrix_materialize_csr_cached(A,&ca,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD ||
       check_csr(csr,(LIS_REAL)1.0))goto done;
#ifdef _OPENMP
    if(omp_get_max_threads()>1 && !omp_get_dynamic() &&
       thread_count(&a)<2)
    {
        fprintf(stderr,"safe USER row callback was not parallel\n");
        goto done;
    }
#endif
    prior=csr;
    a.scale=(LIS_REAL)2.0;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=lis_matrix_materialize_csr_cached(A,&ca,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES ||
       csr!=prior || check_csr(csr,(LIS_REAL)2.0))goto done;
    err=lis_matrix_materialize_csr_cached(A,&ca,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REUSE || csr!=prior)
        goto done;

    a.failed_row=14;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    csr=prior;
    err=lis_matrix_materialize_csr_cached(A,&ca,&csr,&action);
    if(err==LIS_SUCCESS || csr!=NULL)goto done;
    a.failed_row=-1;
    err=lis_matrix_materialize_csr_cached(A,&ca,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES ||
       check_csr(csr,(LIS_REAL)2.0))goto done;

    err=lis_matrix_materialize_csr_cached(U,&cu,&csr,&action);
    if(err || check_csr(csr,(LIS_REAL)1.0) ||
       thread_count(&unsafe)!=1)
    {
        fprintf(stderr,"unsafe callback did not use serial fallback\n");
        goto done;
    }
    err=lis_matrix_create_operator((LIS_SCALAR)1.0,A,
                                   (LIS_SCALAR)-1.0,B,&op);
    if(err)goto done;
    memset(a.seen,0,sizeof(a.seen));
    memset(b.seen,0,sizeof(b.seen));
    err=lis_matrix_materialize_csr_cached(op,&cb,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD ||
       check_csr(csr,(LIS_REAL)1.5))goto done;
#ifdef _OPENMP
    if(omp_get_max_threads()>1 && !omp_get_dynamic() &&
       (thread_count(&a)<2 || thread_count(&b)<2))
    {
        fprintf(stderr,"safe nested operator was not parallel\n");
        goto done;
    }
#endif
    b.scale=(LIS_REAL)0.75;
    err=lis_matrix_mark_values_changed(B);if(err)goto done;
    prior=csr;
    err=lis_matrix_materialize_csr_cached(op,&cb,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES ||
       csr!=prior || check_csr(csr,(LIS_REAL)1.25))goto done;

    failed=0;
done:
    if(ca)lis_matrix_csr_cache_destroy(ca);
    if(cb)lis_matrix_csr_cache_destroy(cb);
    if(cu)lis_matrix_csr_cache_destroy(cu);
    if(op)lis_matrix_destroy(op);
    if(A)lis_matrix_destroy(A);
    if(B)lis_matrix_destroy(B);
    if(U)lis_matrix_destroy(U);
    lis_finalize();
    printf("LIS_MATRIX_MATERIALIZE_OMP %s\n",failed?"FAILED":"PASSED");
    return failed;
}
