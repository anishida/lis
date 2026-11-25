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
    if(c->fail && c->calls==2)
    {
        y[0]=(LIS_SCALAR)NAN;
        return c->fail;
    }
    for(i=0;i<c->n;i++) y[i]=(c->diag+(LIS_SCALAR)i)*x[i];
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
static int run(LIS_INT callback_error)
{
    CTX ac={4,0,callback_error,4.0};
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_SCALAR value=0;
    LIS_REAL residual=0;
    LIS_INT err,status,i,iter;
    int failed=1;
    if(create(&ac,&A))goto done;
    if(lis_vector_duplicate(A,&x) || lis_esolver_create(&e))goto done;
    if(lis_esolver_set_option("-e pi -etol 1e-10 -emaxiter 250 -eprint mem",e))goto done;
    if(lis_vector_set_all(1.0,x))goto done;
    /* Establish a valid result before testing reuse after an error. */
    ac.fail=0;
    if(lis_esolve(A,x,&value,e) || lis_esolver_get_status(e,&status) ||
       status!=LIS_SUCCESS || !(fabs((double)value-7.0)<=1e-8))goto done;
    ac.calls=0; ac.fail=callback_error;
    if(lis_vector_set_all(1.0,x))goto done;
    err=lis_esolve(A,x,&value,e);
    if(lis_esolver_get_status(e,&status) || lis_esolver_get_iter(e,&iter) ||
       lis_esolver_get_residualnorm(e,&residual))goto done;
    if(err!=LIS_SUCCESS || status!=callback_error || ac.calls!=2 ||
       !isnan((double)value) || !isnan((double)residual))goto done;
    ac.calls=0; ac.fail=0;
    if(lis_vector_set_all(1.0,x) || lis_esolve(A,x,&value,e) ||
       lis_esolver_get_status(e,&status) || status!=LIS_SUCCESS ||
       !(fabs((double)value-7.0)<=1e-8))goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"USER eigen status failure code=%d\n",(int)callback_error);
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
    puts("USER eigen status SKIP: real serial regression");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(LIS_FAILS);
    failed|=run(12345);
    lis_finalize();
    printf("USER_EIGEN_STATUS %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
