#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include "lis.h"

int main(int argc, char *argv[])
{
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER esolver=NULL;
    LIS_SCALAR lambda=0.0;
    LIS_INT err,i,status;
    LIS_INT retries,min_ss,attempts,final_ss;
    LIS_REAL growth;
    const LIS_INT n=12;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A); if(err) goto fail;
    err=lis_matrix_set_size(A,0,n); if(err) goto fail;
    for(i=0;i<n;i++){ err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)(i+1),A); if(err) goto fail; }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR); if(err) goto fail;
    err=lis_matrix_assemble(A); if(err) goto fail;

    err=lis_matrix_create(LIS_COMM_WORLD,&B); if(err) goto fail;
    err=lis_matrix_set_size(B,0,n); if(err) goto fail;
    for(i=0;i<n;i++){ err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,B); if(err) goto fail; }
    err=lis_matrix_set_type(B,LIS_MATRIX_CSR); if(err) goto fail;
    err=lis_matrix_assemble(B); if(err) goto fail;

    err=lis_vector_duplicate(A,&x); if(err) goto fail;
    err=lis_vector_set_all((LIS_SCALAR)1.0,x); if(err) goto fail;
    err=lis_esolver_create(&esolver); if(err) goto fail;

    err=lis_esolver_set_subspace_retry(esolver,2,(LIS_REAL)2.0,6); if(err) goto fail;
    err=lis_esolver_get_subspace_retry(esolver,&retries,&growth,&min_ss,&attempts,&final_ss);
    if(err || retries!=2 || growth!=(LIS_REAL)2.0 || min_ss!=6 || attempts!=0) goto fail;

    err=lis_esolver_set_option("-e gli -nev 2 -ss 3 -emaxiter 1 -etol 1.0e-30 -eprint none",esolver);
    if(err) goto fail;
    err=lis_gesolve(A,B,x,&lambda,esolver); if(err) goto fail;
    err=lis_esolver_get_status(esolver,&status); if(err) goto fail;
    err=lis_esolver_get_subspace_retry(esolver,NULL,NULL,NULL,&attempts,&final_ss); if(err) goto fail;

    if(status==LIS_MAXITER && (attempts<1 || final_ss<6))
    {
        fprintf(stderr,"retry did not grow subspace: attempts=%d final_ss=%d\n",(int)attempts,(int)final_ss);
        goto fail;
    }
    if(status!=LIS_SUCCESS && status!=LIS_MAXITER)
    {
        fprintf(stderr,"unexpected status=%d\n",(int)status);
        goto fail;
    }

    if(lis_esolver_set_subspace_retry(esolver,-1,(LIS_REAL)2.0,0)!=LIS_ERR_ILL_ARG) goto fail;
    if(lis_esolver_set_subspace_retry(esolver,1,(LIS_REAL)1.0,0)!=LIS_ERR_ILL_ARG) goto fail;

    lis_esolver_destroy(esolver);
    lis_vector_destroy(x);
    lis_matrix_destroy(B);
    lis_matrix_destroy(A);
    lis_finalize();
    return 0;
fail:
    if(esolver) lis_esolver_destroy(esolver);
    if(x) lis_vector_destroy(x);
    if(B) lis_matrix_destroy(B);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return 1;
}
