#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

#if !defined(USE_MPI) && !defined(_COMPLEX)
typedef struct { LIS_INT n,calls,fail; LIS_SCALAR diag; LIS_INT code; } CTX;
static LIS_INT apply(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)ptr;
    LIS_INT i;
    ++c->calls;
    if(c->fail && c->calls==c->fail)
    {
        y[0]=(LIS_SCALAR)NAN;
        return c->code;
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
static int run(const char *method,int generalized,LIS_INT afail,LIS_INT bfail,LIS_INT code)
{
    CTX ac={4,0,afail,4.0,code},bc={4,0,bfail,2.0,code};
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_SCALAR value=0;
    LIS_REAL residual=0;
    LIS_INT err,status,i,iter,attempts;
    int failed=1;
    if(create(&ac,&A))goto done;
    if(generalized && create(&bc,&B))goto done;
    if(lis_vector_duplicate(A,&x) || lis_esolver_create(&e))goto done;
    char options[160];
    snprintf(options,sizeof(options),"-e %s -ss 2 -ie pi -ige gpi -etol 1e-10 -emaxiter 50 -eprint none",method);
    if(lis_esolver_set_option(options,e))goto done;
    if(code==LIS_MAXITER && lis_esolver_set_subspace_retry(e,2,2.0,3))goto done;
    if(lis_vector_set_all(1.0,x))goto done;
    err=generalized ? lis_gesolve(A,B,x,&value,e) : lis_esolve(A,x,&value,e);
    if(lis_esolver_get_status(e,&status) || lis_esolver_get_iter(e,&iter) ||
       lis_esolver_get_residualnorm(e,&residual))goto done;
    if(afail || bfail)
    {
        if(err!=LIS_SUCCESS || status!=code ||
           (afail && ac.calls<afail) || (bfail && bc.calls<bfail))
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
    if(code==LIS_MAXITER &&
       (lis_esolver_get_subspace_retry(e,NULL,NULL,NULL,&attempts,NULL) ||
        attempts!=0 || !isnan((double)value) || !isnan((double)residual) ||
        (afail && ac.calls!=afail) || (bfail && bc.calls!=bfail)))goto done;
    for(i=0;i<x->n;i++)if(!isfinite((double)x->value[i]))goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"USER eigen projection failure gen=%d A=%d B=%d\n",
                      generalized,(int)afail,(int)bfail);
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    return failed;
}
#endif
int main(int argc,char **argv)
{
#if defined(USE_MPI) || defined(_COMPLEX)
    (void)argc;(void)argv;
    puts("USER eigen projection SKIP: real serial regression");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run("ai",0,1,0,LIS_ERR_ILL_ARG);
    failed|=run("gai",1,1,0,LIS_ERR_ILL_ARG);
    failed|=run("gai",1,0,1,LIS_ERR_ILL_ARG);
    failed|=run("li",0,1,0,LIS_ERR_ILL_ARG);
    failed|=run("gli",1,1,0,LIS_ERR_ILL_ARG);
    failed|=run("gli",1,0,1,LIS_ERR_ILL_ARG);
    failed|=run("si",0,1,0,LIS_ERR_ILL_ARG);
    failed|=run("gsi",1,1,0,LIS_ERR_ILL_ARG);
    failed|=run("gsi",1,0,1,LIS_ERR_ILL_ARG);
    failed|=run("ai",0,1,0,LIS_MAXITER);
    failed|=run("li",0,1,0,LIS_MAXITER);
    failed|=run("si",0,1,0,LIS_MAXITER);
    lis_finalize();
    printf("USER_EIGEN_PROJECTION %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
