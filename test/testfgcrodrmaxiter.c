/* Reproduce the ordinary FGCRO-DR MAXITER path with bcoef allocated.
 * Under LeakSanitizer this also checks workspace cleanup.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include <math.h>
#include "lis.h"
#include "lis_solver.h"
#include "lislib.h"

#ifdef _COMPLEX
int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    puts("SKIP: FGCRO-DR positive-k real-only regression");
    return 77;
}
#else

#define N 64

int main(int argc, char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR u=NULL,b=NULL,x=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,dim=0,failed=0,status1=-1,status2=-1;
    LIS_UNSIGNED_INT builds1=0,uses1=0,uses2=0;
    LIS_REAL orth=0.0,relation=0.0,residual_norm=0.0;

    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    if((err=lis_matrix_create(LIS_COMM_WORLD,&A))) goto fail;
    if((err=lis_matrix_set_size(A,0,N))) goto fail;
    for(i=0;i<N;i++)
    {
        if((err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A))) goto fail;
        if(i>0 &&
           (err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A))) goto fail;
        if(i<N-1 &&
           (err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.5,A)))
            goto fail;
    }
    if((err=lis_matrix_set_type(A,LIS_MATRIX_CSR))) goto fail;
    if((err=lis_matrix_assemble(A))) goto fail;
    if((err=lis_vector_duplicate(A,&u))) goto fail;
    if((err=lis_vector_duplicate(A,&b))) goto fail;
    if((err=lis_vector_duplicate(A,&x))) goto fail;
    if((err=lis_vector_duplicate(A,&r))) goto fail;
    if((err=lis_solver_create(&solver))) goto fail;
    if((err=lis_solver_set_option(
        "-i fgcrodr -p jacobi -restart 4 -maxiter 200 "
        "-tol 1e-11 -fgcrodr_k 2 -print none",solver))) goto fail;

    /* First solve: build a valid persistent recycle basis. */
    for(i=0;i<N;i++)
        u->value[i]=(LIS_SCALAR)(1.0+(i%11)*0.123);
    if((err=lis_matvec(A,u,b))) goto fail;
    if((err=lis_vector_set_all(0.0,x))) goto fail;
    err=lis_solve_sequence(A,b,x,solver);
    status1=solver->retcode;
    if(err!=LIS_SUCCESS || status1!=LIS_SUCCESS) goto fail;
    if((err=lis_solver_get_fgcrodr_recycle(solver,&dim,NULL,NULL,NULL)))
        goto fail;
    if((err=lis_fgcrodr_get_stats(solver,&builds1,NULL,&uses1,NULL)))
        goto fail;
    if(dim<=0 || builds1==0) {err=LIS_ERR_ILL_ARG; goto fail;}

    /* Second solve: keep same A and recycle object; force ordinary
     * maxiter exit, not the separate no-improvement exit. */
    solver->options[LIS_OPTIONS_MAXITER]=1;
    for(i=0;i<N;i++)
        u->value[i]=(LIS_SCALAR)((i%13)*0.125-0.75);
    if((err=lis_matvec(A,u,b))) goto fail;
    if((err=lis_vector_set_all(0.0,x))) goto fail;
    err=lis_solve_sequence(A,b,x,solver);
    status2=solver->retcode;
    if((err!=LIS_SUCCESS && err!=LIS_MAXITER) || status2!=LIS_MAXITER)
        goto fail;
    if((err=lis_fgcrodr_get_stats(solver,NULL,NULL,&uses2,NULL)))
        goto fail;
    if((err=lis_fgcrodr_validate_recycle(solver,&orth,&relation)))
        goto fail;
    if(uses2<=uses1 || !isfinite(orth) || !isfinite(relation)
       || orth>1e-8 || relation>1e-8)
    {
        err=LIS_ERR_ILL_ARG;
        goto fail;
    }
    if((err=lis_matvec(A,x,r))) goto fail;
    if((err=lis_vector_xpay(b,-1.0,r))) goto fail;
    if((err=lis_vector_nrm2(r,&residual_norm))) goto fail;
    if(!isfinite(residual_norm) || residual_norm<1e-9)
    {
        err=LIS_ERR_ILL_ARG;
        goto fail;
    }
    printf("FGCRODR_ACTIVE_MAXITER_CLEANUP_OK "
           "first=%d second=%d recycle_dim=%d uses=%llu->%llu\n",
           (int)status1,(int)status2,(int)dim,
           (unsigned long long)uses1,(unsigned long long)uses2);
    goto cleanup;
fail:
    fprintf(stderr,"FGCRODR_ACTIVE_MAXITER_CLEANUP_FAIL "
            "api=%d first=%d second=%d dim=%d\n",
            (int)err,(int)status1,(int)status2,(int)dim);
    failed=1;
cleanup:
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(u) lis_vector_destroy(u);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
#endif
