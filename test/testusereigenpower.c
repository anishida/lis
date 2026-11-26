#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#ifdef LIS_TEST_OOM
#include "lis_system.h"
#endif
#include <math.h>
#include <stdio.h>

#if !defined(USE_MPI) && !defined(_COMPLEX)
typedef struct { LIS_INT n,calls,fail; LIS_SCALAR diag; } CTX;
static LIS_INT apply(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)ptr;
    LIS_INT i;
    ++c->calls;
    if(c->fail && c->calls==c->fail)
    {
        y[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<c->n;i++) y[i]=c->diag*x[i];
    return LIS_SUCCESS;
}
static LIS_INT create(CTX *c,LIS_MATRIX *A)
{
    LIS_INT err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err)return err;
    err=lis_matrix_set_size(*A,0,c->n);
    if(err)return err;
    return lis_matrix_set_user(*A,c,apply,apply);
}
static int run(int generalized,LIS_INT afail,LIS_INT bfail)
{
    CTX ac={4,0,afail,4.0},bc={4,0,bfail,2.0};
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_SCALAR value=0;
    LIS_REAL residual=0;
    LIS_INT err,status,i,iter;
    int failed=1;
    if(create(&ac,&A))goto done;
    if(generalized && create(&bc,&B))goto done;
    if(lis_vector_duplicate(A,&x) || lis_esolver_create(&e))goto done;
    if(lis_esolver_set_option(generalized ?
        "-e gpi -etol 1e-10 -emaxiter 50 -eprint none" :
        "-e pi -etol 1e-10 -emaxiter 50 -eprint none",e))goto done;
    if(lis_vector_set_all(1.0,x))goto done;
    err=generalized ? lis_gesolve(A,B,x,&value,e) : lis_esolve(A,x,&value,e);
    if(lis_esolver_get_status(e,&status) || lis_esolver_get_iter(e,&iter) ||
       lis_esolver_get_residualnorm(e,&residual))goto done;
    if(afail || bfail)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG || iter!=0 ||
           !isnan((double)residual) ||
           (afail && ac.calls!=afail) || (bfail && bc.calls<bfail))
        {
            fprintf(stderr,"fault gen=%d A=%d B=%d api=%d status=%d iter=%d calls=%d/%d resid=%g\n",
                generalized,(int)afail,(int)bfail,(int)err,(int)status,(int)iter,
                (int)ac.calls,(int)bc.calls,(double)residual);
            goto done;
        }
    }
    else if(err!=LIS_SUCCESS || status!=LIS_SUCCESS ||
            !(fabs((double)value-(generalized?2.0:4.0))<=1e-9) ||
            !(residual<=1e-10))goto done;
    for(i=0;i<x->n;i++)if(!isfinite((double)x->value[i]))goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"USER eigen power failure gen=%d A=%d B=%d\n",
                      generalized,(int)afail,(int)bfail);
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* A nontrivial spectrum forces multiple iterations. Sweep every application
 * observed in the healthy run, including final physical certification, and
 * reuse the same eigensolver after every failure. */
typedef struct { LIS_INT calls,fail,code; } SWEEP_CTX;
static LIS_INT sweep_apply(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    SWEEP_CTX *c=(SWEEP_CTX *)raw;
    const double diag[4]={1.0,2.0,3.0,8.0};
    LIS_INT i;
    ++c->calls;
    if(c->fail && c->calls==c->fail){
        y[0]=(LIS_SCALAR)NAN;
        return c->code;
    }
    for(i=0;i<4;i++)y[i]=diag[i]*x[i];
    return LIS_SUCCESS;
}
static int sweep_healthy(LIS_MATRIX A,LIS_VECTOR x,LIS_ESOLVER e,SWEEP_CTX *c)
{
    const double diag[4]={1.0,2.0,3.0,8.0};
    LIS_SCALAR value=0;
    LIS_INT i,status;
    double r2=0.0,x2=0.0;
    c->calls=0; c->fail=0;
    if(lis_vector_set_all(1.0,x) || lis_esolve(A,x,&value,e) ||
       lis_esolver_get_status(e,&status) || status!=LIS_SUCCESS ||
       !(fabs((double)value-8.0)<=1e-9))return 1;
    for(i=0;i<4;i++){
        double v=(double)x->value[i],r=(diag[i]-(double)value)*v;
        r2+=r*r; x2+=v*v;
    }
    return !(x2>0.0 && sqrt(r2/x2)/8.0<=1e-9 && c->calls>4);
}
static int sweep(void)
{
    SWEEP_CTX c={0,0,LIS_ERR_ILL_ARG};
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_INT count=0,k,j,status;
    const LIS_INT codes[]={LIS_ERR_ILL_ARG,LIS_BREAKDOWN,LIS_MAXITER,12345};
    LIS_SCALAR value;
    LIS_REAL residual;
    int failed=1;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,4) ||
       lis_matrix_set_user(A,&c,sweep_apply,sweep_apply) ||
       lis_vector_duplicate(A,&x) || lis_esolver_create(&e) ||
       lis_esolver_set_option("-e pi -etol 1e-10 -emaxiter 200 -eprint none",e) ||
       sweep_healthy(A,x,e,&c))goto done;
    count=c.calls;
    for(j=0;j<(LIS_INT)(sizeof(codes)/sizeof(codes[0]));j++){
        for(k=1;k<=count;k++){
            c.calls=0; c.fail=k; c.code=codes[j]; value=0;
            if(lis_vector_set_all(1.0,x) || lis_esolve(A,x,&value,e) ||
               lis_esolver_get_status(e,&status) || status!=c.code ||
               lis_esolver_get_residualnorm(e,&residual) ||
               !isnan((double)value) || !isnan((double)residual) ||
               c.calls!=k)goto done;
            if(sweep_healthy(A,x,e,&c) || c.calls!=count)goto done;
        }
    }
    failed=0;
done:
    if(failed)fprintf(stderr,"USER power sweep failed: code=%d fail=%d calls=%d baseline=%d\n",
                      (int)c.code,(int)c.fail,(int)c.calls,(int)count);
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}

#ifdef LIS_TEST_OOM
/* Four result arrays, the initial vector descriptor and values, then
 * eigenvector slots and residual history. Fail each staged allocation
 * while a previous successful result is still owned by this solver. */
static int sweep_result_oom(void)
{
    SWEEP_CTX c={0,0,LIS_ERR_ILL_ARG};
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_INT status,err;
    LIS_SCALAR value=0;
    size_t k=0,before,live=lis_test_oom_live();
    int failed=1;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,4) ||
       lis_matrix_set_user(A,&c,sweep_apply,sweep_apply) ||
       lis_vector_duplicate(A,&x) || lis_esolver_create(&e) ||
       lis_esolver_set_option("-e pi -etol 1e-10 -emaxiter 200 -eprint none",e) ||
       sweep_healthy(A,x,e,&c))goto done;
    for(k=1;k<=8;k++){
        before=lis_test_oom_live();
        c.calls=0; value=0;
        lis_test_oom_begin(k);
        err=lis_esolve(A,x,&value,e);
        lis_test_oom_end();
        if(err!=LIS_OUT_OF_MEMORY || lis_esolver_get_status(e,&status) ||
           status!=LIS_OUT_OF_MEMORY || !isnan((double)value) ||
           c.calls!=0 || lis_test_oom_live()!=before ||
           sweep_healthy(A,x,e,&c))goto done;
    }
    failed=0;
done:
    lis_test_oom_end();
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    if(lis_test_oom_live()!=live)failed=1;
    if(failed)fprintf(stderr,"USER eigen result OOM failed at %zu\n",k);
    return failed;
}
#endif
#endif
int main(int argc,char **argv)
{
#if defined(USE_MPI) || defined(_COMPLEX)
    (void)argc;(void)argv;
    puts("USER eigen power SKIP: real serial regression");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(0,1,0);
    failed|=run(0,2,0);
    failed|=run(1,1,0);
    failed|=run(1,2,0);
    failed|=run(1,0,1);
    failed|=run(0,0,0);
    failed|=run(1,0,0);
    failed|=sweep();
#ifdef LIS_TEST_OOM
    failed|=sweep_result_oom();
#endif
    lis_finalize();
    printf("USER_EIGEN_POWER %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
