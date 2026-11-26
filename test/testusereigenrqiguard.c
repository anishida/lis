#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
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
        "-e grqi -eprint none" : "-e rqi -eprint none",e))goto done;
    if(lis_vector_set_all(1.0,x))goto done;
    err=generalized ? lis_gesolve(A,B,x,&value,e) : lis_esolve(A,x,&value,e);
    if(err!=LIS_ERR_NOT_IMPLEMENTED || ac.calls || bc.calls)goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"USER eigen RQI guard failure gen=%d A=%d B=%d\n",
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
    puts("USER eigen RQI guard SKIP: real serial regression");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(0,0,0);
    failed|=run(1,0,0);
    lis_finalize();
    printf("USER_EIGEN_RQI_GUARD %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
