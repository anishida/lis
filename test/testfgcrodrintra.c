/*
 * Standalone (non-MFAMG) regression of within-solve FGCRO-DR recycling.
 * Uses the same matrix and right-hand side with the feature OFF and ON.
 * Verifies a real restart, recycle projection, true residual, and AZ=C.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lis.h"
#include "lis_solver.h"
#include "lislib.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _COMPLEX
int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    puts("SKIP: FGCRO-DR positive-k real-only regression");
    return 77;
}
#else
#define N 120

/* Solver-specific options: do not modify global process environment. */
static int set_intra(LIS_SOLVER solver, int enabled)
{
    return lis_solver_set_option(
        enabled
        ? "-fgcrodr_intra 1 -fgcrodr_projected_relative 1 -fgcrodr_intra_trace 1"
        : "-fgcrodr_intra 0 -fgcrodr_projected_relative 1",
        solver)!=LIS_SUCCESS;
}

static int solve_one(LIS_MATRIX A, LIS_VECTOR b, int intra, int requested_k, int trace)
{
    LIS_VECTOR x=NULL, ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT rc,dim=0,action=0,iter=0,status=0;
    LIS_UNSIGNED_INT builds=0,updates=0,uses=0,rejects=0;
    LIS_REAL orth=0.0,relation=0.0,bnorm=0.0,rnorm=0.0;
    char opts[256];
    int failed=1;
    if((rc=lis_vector_duplicate(A,&x))) goto done;
    if((rc=lis_vector_duplicate(A,&ax))) goto done;
    if((rc=lis_solver_create(&solver))) goto done;
    if(set_intra(solver,intra)) goto done;
    /* Explicitly exercise the safety check without diagnostic trace. */
    if(intra && !trace &&
       (rc=lis_solver_set_option("-fgcrodr_intra_trace 0",solver)))
        goto done;
    snprintf(opts,sizeof(opts),
        "-i fgcrodr -p jacobi -restart 2 -maxiter 300 "
        "-tol 1e-5 -fgcrodr_k %d -print none",requested_k);
    if((rc=lis_solver_set_option(opts,solver))) goto done;
    if((rc=lis_vector_set_all(0.0,x))) goto done;
    rc=lis_solve(A,b,x,solver);
    status=solver->retcode;
    iter=solver->iter;
    if(rc!=LIS_SUCCESS || status!=LIS_SUCCESS || iter<=2)
        goto done;
    if((rc=lis_solver_get_fgcrodr_recycle(
        solver,&dim,&action,NULL,NULL))) goto done;
    if((rc=lis_fgcrodr_get_stats(
        solver,&builds,&updates,&uses,&rejects))) goto done;
    if((rc=lis_fgcrodr_validate_recycle(
        solver,&orth,&relation))) goto done;
    if((rc=lis_matvec(A,x,ax))) goto done;
    if((rc=lis_vector_xpay(b,-1.0,ax))) goto done;
    if((rc=lis_vector_nrm2(ax,&rnorm))) goto done;
    if((rc=lis_vector_nrm2(b,&bnorm))) goto done;
    if(bnorm==0.0 || rnorm/bnorm>1e-5 ||
       !isfinite(orth) || !isfinite(relation) ||
       orth>1e-8 || relation>1e-8 ||
       builds<1 || rejects!=0 || dim!=1 ||
       (intra && (uses!=1 || updates<1)) ||
       (!intra && uses!=0))
        goto done;
    printf("FGCRODR_INTRA_%s_OK iter=%d true_relative=%.12e "
           "k=%d builds=%llu updates=%llu uses=%llu rejects=%llu "
           "orth=%.6e relation=%.6e\n",
           intra?"ON":"OFF",(int)iter,(double)(rnorm/bnorm),
           (int)dim,(unsigned long long)builds,
           (unsigned long long)updates,(unsigned long long)uses,
           (unsigned long long)rejects,(double)orth,(double)relation);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_INTRA_%s_FAIL status=%d iter=%d "
                "dim=%d updates=%llu uses=%llu rejects=%llu "
                "orth=%.6e relation=%.6e rnorm=%.6e\n",
                intra?"ON":"OFF",(int)status,(int)iter,(int)dim,
                (unsigned long long)updates,(unsigned long long)uses,
                (unsigned long long)rejects,(double)orth,
                (double)relation,(double)rnorm);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    return failed;
}


/* A synthetic changed-operator sequence.  Reusing the solver across
 * distinct matrix objects must not project the new residual against
 * a stale A1*Z=C relation for A2.
 */
static int check_changed_operator(LIS_MATRIX A1, LIS_VECTOR b)
{
    LIS_MATRIX A2=NULL,stored=NULL;
    LIS_VECTOR x=NULL,ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,is=0,ie=0,valid=0,failed=1;
    LIS_UNSIGNED_INT sg=0,vg=0;
    LIS_REAL orth=0.0,relation=0.0,rnorm=0.0,bnorm=0.0;

    if((err=lis_matrix_create(LIS_COMM_WORLD,&A2))) goto done;
    if((err=lis_matrix_set_size(A2,0,N))) goto done;
    if((err=lis_matrix_get_range(A2,&is,&ie))) goto done;
    for(i=is;i<ie;i++)
    {
        if((err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.25,A2)))
            goto done;
        if(i>0 &&
           (err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A2)))
            goto done;
        if(i<N-1 &&
           (err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.5,A2)))
            goto done;
    }
    if((err=lis_matrix_set_type(A2,LIS_MATRIX_CSR))) goto done;
    if((err=lis_matrix_assemble(A2))) goto done;
    if((err=lis_vector_duplicate(A1,&x))) goto done;
    if((err=lis_vector_duplicate(A1,&ax))) goto done;
    if((err=lis_solver_create(&solver))) goto done;
    if(set_intra(solver,1)) goto done;
    if((err=lis_solver_set_option(
        "-i fgcrodr -p jacobi -restart 2 -maxiter 300 "
        "-tol 1e-5 -fgcrodr_k 1 -print none",solver))) goto done;

    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A1,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_lifecycle_state(
        solver,&valid,&stored,&sg,&vg))) goto done;
    if(!valid || stored!=A1) goto done;

    /* New assembled matrix object, same layout, modified diagonal. */
    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A2,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_lifecycle_state(
        solver,&valid,&stored,&sg,&vg))) goto done;
    if(!valid || stored!=A2 ||
       sg!=A2->structure_generation ||
       vg!=A2->value_generation)
        goto done;
    if((err=lis_fgcrodr_validate_recycle(
        solver,&orth,&relation))) goto done;
    if((err=lis_matvec(A2,x,ax))) goto done;
    if((err=lis_vector_xpay(b,-1.0,ax))) goto done;
    if((err=lis_vector_nrm2(ax,&rnorm))) goto done;
    if((err=lis_vector_nrm2(b,&bnorm))) goto done;
    if(bnorm==0.0 || !isfinite(rnorm) || rnorm/bnorm>1e-5 ||
       !isfinite(orth) || !isfinite(relation) ||
       orth>1e-8 || relation>1e-8)
        goto done;
    printf("FGCRODR_INTRA_CHANGED_OPERATOR_OK "
           "A1_to_A2_status=%d residual=%.12e relation=%.6e\n",
           (int)solver->retcode,(double)(rnorm/bnorm),(double)relation);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_INTRA_CHANGED_OPERATOR_FAIL "
                "api=%d status=%d valid=%d stored_matches_A2=%d "
                "relation=%.6e\n",
                (int)err,solver?(int)solver->retcode:-1,
                (int)valid,(int)(stored==A2),(double)relation);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    if(A2) lis_matrix_destroy(A2);
    return failed;
}


/* Same matrix POINTER and sparsity, but changed numerical values.
 * A new value_generation must invalidate/rebuild the old AZ=C relation.
 * This tests a different path from the distinct-object A1 -> A2 check.
 */
static int check_same_operator_value_update(LIS_MATRIX A, LIS_VECTOR b)
{
    LIS_VECTOR x=NULL,ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_MATRIX stored=NULL;
    LIS_INT err=0,valid=0,failed=1;
    LIS_UNSIGNED_INT sg_before=0,vg_before=0,sg_after=0,vg_after=0;
    LIS_UNSIGNED_INT saved_sg=0,saved_vg=0;
    LIS_REAL orth=0.0,relation=0.0,rnorm=0.0,bnorm=0.0;

    if((err=lis_vector_duplicate(A,&x))) goto done;
    if((err=lis_vector_duplicate(A,&ax))) goto done;
    if((err=lis_solver_create(&solver))) goto done;
    if(set_intra(solver,1)) goto done;
    if((err=lis_solver_set_option(
        "-i fgcrodr -p jacobi -restart 2 -maxiter 300 "
        "-tol 1e-5 -fgcrodr_k 1 -print none",solver))) goto done;

    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_lifecycle_state(
        solver,&valid,&stored,&saved_sg,&saved_vg))) goto done;
    if(!valid || stored!=A) goto done;
    if((err=lis_matrix_get_generation(A,&sg_before,&vg_before)))
        goto done;
    if(saved_sg!=sg_before || saved_vg!=vg_before) goto done;

    /* shift_diagonal changes numerical storage; the caller must notify
     * the generation tracker so cross-solve recycling is revalidated. */
    if((err=lis_matrix_shift_diagonal(A,(LIS_SCALAR)-0.25)))
        goto done;
    if((err=lis_matrix_mark_values_changed(A)))
        goto done;
    if((err=lis_matrix_get_generation(A,&sg_after,&vg_after)))
        goto done;
    if(sg_after!=sg_before || vg_after<=vg_before)
        goto done;

    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_lifecycle_state(
        solver,&valid,&stored,&saved_sg,&saved_vg))) goto done;
    if(!valid || stored!=A ||
       saved_sg!=sg_after || saved_vg!=vg_after)
        goto done;
    if((err=lis_fgcrodr_validate_recycle(
        solver,&orth,&relation))) goto done;
    if((err=lis_matvec(A,x,ax))) goto done;
    if((err=lis_vector_xpay(b,-1.0,ax))) goto done;
    if((err=lis_vector_nrm2(ax,&rnorm))) goto done;
    if((err=lis_vector_nrm2(b,&bnorm))) goto done;
    if(bnorm==0.0 || !isfinite(rnorm) || rnorm/bnorm>1e-5 ||
       !isfinite(orth) || !isfinite(relation) ||
       orth>1e-8 || relation>1e-8)
        goto done;
    printf("FGCRODR_INTRA_VALUE_GENERATION_OK "
           "generation=%llu->%llu relative=%.12e relation=%.6e\n",
           (unsigned long long)vg_before,
           (unsigned long long)vg_after,
           (double)(rnorm/bnorm),(double)relation);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_INTRA_VALUE_GENERATION_FAIL "
                "api=%d status=%d valid=%d generation=%llu->%llu "
                "saved=%llu relation=%.6e\n",
                (int)err,solver?(int)solver->retcode:-1,(int)valid,
                (unsigned long long)vg_before,
                (unsigned long long)vg_after,
                (unsigned long long)saved_vg,(double)relation);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    return failed;
}


/* A rank-two saved recycle basis built with restart=4 must not be
 * used at restart=2. Instead invalidate it before the first cycle
 * and rebuild at most one direction with the smaller restart.
 */
static int check_restart_shrink(LIS_MATRIX A,LIS_VECTOR b,int intra)
{
    LIS_VECTOR x=NULL,ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT rc=0,dim_before=0,dim_after=0,failed=1;
    LIS_INT action=0,status=0;
    LIS_UNSIGNED_INT builds=0,rejects_before=0,rejects_after=0;
    LIS_REAL orth=0.0,relation=0.0,rnorm=0.0,bnorm=0.0;

    if((rc=lis_vector_duplicate(A,&x))) goto done;
    if((rc=lis_vector_duplicate(A,&ax))) goto done;
    if((rc=lis_solver_create(&solver))) goto done;
    if(set_intra(solver,intra)) goto done;
    if((rc=lis_solver_set_option(
        "-i fgcrodr -p jacobi -restart 4 -maxiter 300 "
        "-tol 1e-5 -fgcrodr_k 2 -print none",solver))) goto done;
    if((rc=lis_vector_set_all(0.0,x))) goto done;
    rc=lis_solve_sequence(A,b,x,solver);
    if(rc!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((rc=lis_solver_get_fgcrodr_recycle(
        solver,&dim_before,&action,NULL,NULL))) goto done;
    if((rc=lis_fgcrodr_get_stats(
        solver,&builds,NULL,NULL,&rejects_before))) goto done;
    if(dim_before!=2 || builds==0) goto done;

    /* Internal sequence test: change the effective restart but retain
     * the solver and its valid cross-solve recycle state. */
    solver->options[LIS_OPTIONS_RESTART]=2;
    if((rc=lis_vector_set_all(0.0,x))) goto done;
    rc=lis_solve_sequence(A,b,x,solver);
    status=solver->retcode;
    if(rc!=LIS_SUCCESS || status!=LIS_SUCCESS) goto done;
    if((rc=lis_solver_get_fgcrodr_recycle(
        solver,&dim_after,&action,NULL,NULL))) goto done;
    if((rc=lis_fgcrodr_get_stats(
        solver,&builds,NULL,NULL,&rejects_after))) goto done;
    if(dim_after!=1 || rejects_after<=rejects_before) goto done;
    if((rc=lis_fgcrodr_validate_recycle(
        solver,&orth,&relation))) goto done;
    if((rc=lis_matvec(A,x,ax))) goto done;
    if((rc=lis_vector_xpay(b,-1.0,ax))) goto done;
    if((rc=lis_vector_nrm2(ax,&rnorm))) goto done;
    if((rc=lis_vector_nrm2(b,&bnorm))) goto done;
    if(bnorm<=0.0 || !isfinite(rnorm) || rnorm/bnorm>1e-5 ||
       !isfinite(orth) || !isfinite(relation) ||
       orth>1e-8 || relation>1e-8) goto done;

    printf("FGCRODR_RESTART_SHRINK_%s_OK "
           "old_k=%d new_k=%d rejects=%llu->%llu "
           "true_relative=%.12e relation=%.6e\n",
           intra?"INTRA":"CROSS",(int)dim_before,(int)dim_after,
           (unsigned long long)rejects_before,
           (unsigned long long)rejects_after,
           (double)(rnorm/bnorm),(double)relation);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_RESTART_SHRINK_%s_FAIL "
                "rc=%d status=%d before=%d after=%d "
                "rejects=%llu->%llu\n",
                intra?"INTRA":"CROSS",
                (int)rc,(int)status,(int)dim_before,(int)dim_after,
                (unsigned long long)rejects_before,
                (unsigned long long)rejects_after);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    return failed;
}


/* A per-rank mode mismatch must fail deterministically before Krylov
 * collectives; otherwise inconsistent paths can hang on MPI.
 */
static int check_mpi_mode_mismatch(LIS_MATRIX A,LIS_VECTOR b,int mode)
{
#ifdef USE_MPI
    int rank=0,nprocs=1;
    LIS_VECTOR x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,failed=1;
    MPI_Comm_rank(A->comm,&rank);
    MPI_Comm_size(A->comm,&nprocs);
    if(nprocs<2) return 0;
    if((err=lis_vector_duplicate(A,&x))) goto done;
    if((err=lis_solver_create(&solver))) goto done;
    if((err=lis_solver_set_option(
        "-i fgcrodr -p jacobi -restart 2 -maxiter 10 "
        "-tol 1e-5 -fgcrodr_k 1 -fgcrodr_projected_relative 1 "
        "-print none",solver))) goto done;
    if(mode==0)
    {
        if((err=lis_solver_set_option(
            rank==0 ? "-fgcrodr_intra 1" : "-fgcrodr_intra 0",
            solver))) goto done;
    }
    else if(mode==1)
    {
        /* Both ranks have k>0; testing a Boolean would miss this. */
        if((err=lis_solver_set_option(
            rank==0 ? "-fgcrodr_k 1" : "-fgcrodr_k 2",
            solver))) goto done;
    }
    else if(mode==2)
    {
        if((err=lis_solver_set_option(
            rank==0 ? "-restart 2" : "-restart 3",
            solver))) goto done;
    }
    else if(mode==3)
    {
        if((err=lis_solver_set_option(
            rank==0 ? "-maxiter 10" : "-maxiter 11",
            solver))) goto done;
    }
    else
        goto done;
    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A,b,x,solver);
    /* LIS solve-sequence API can return success for completed
     * execution while the iterative failure is in retcode. */
    if((err!=LIS_SUCCESS && err!=LIS_ERR_ILL_ARG) ||
       solver->retcode!=LIS_ERR_ILL_ARG)
        goto done;
    if(rank==0)
        printf("FGCRODR_MPI_MODE_MISMATCH_REJECTED_OK mode=%d\n",mode);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_MPI_MODE_MISMATCH_FAIL "
                "mode=%d rank=%d err=%d retcode=%d\n",
                mode,rank,(int)err,solver?(int)solver->retcode:-1);
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    return failed;
#else
    (void)A; (void)b;
    return 0;
#endif
}

int main(int argc, char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR u=NULL,b=NULL;
    LIS_INT rc=0,i,is=0,ie=0,failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,&A))) goto done;
    if((rc=lis_matrix_set_size(A,0,N))) goto done;
    if((rc=lis_matrix_get_range(A,&is,&ie))) goto done;
    for(i=is;i<ie;i++)
    {
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A)))
            goto done;
        if(i>0 &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A)))
            goto done;
        if(i<N-1 &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.5,A)))
            goto done;
    }
    if((rc=lis_matrix_set_type(A,LIS_MATRIX_CSR))) goto done;
    if((rc=lis_matrix_assemble(A))) goto done;
    if((rc=lis_vector_duplicate(A,&u))) goto done;
    if((rc=lis_vector_duplicate(A,&b))) goto done;
    for(i=0;i<u->n;i++)
        u->value[i]=(LIS_SCALAR)((((u->is+i)*37)%23)-11.0);
    if((rc=lis_matvec(A,u,b))) goto done;
    if(solve_one(A,b,0,1,0) ||
       solve_one(A,b,1,1,1) ||
       solve_one(A,b,1,1,0) ||
       solve_one(A,b,1,1000000,0)) goto done;
    if(check_changed_operator(A,b)) goto done;
    /* Run the restart tests on the unmodified reference operator.
     * The value-generation test intentionally changes A in place. */
    if(check_restart_shrink(A,b,0) ||
       check_restart_shrink(A,b,1)) goto done;
    if(check_same_operator_value_update(A,b)) goto done;
    if(check_mpi_mode_mismatch(A,b,0) ||
       check_mpi_mode_mismatch(A,b,1) ||
       check_mpi_mode_mismatch(A,b,2) ||
       check_mpi_mode_mismatch(A,b,3)) goto done;
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_INTRA_REGRESSION_FAIL rc=%d\n",(int)rc);
    if(b) lis_vector_destroy(b);
    if(u) lis_vector_destroy(u);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
#endif
