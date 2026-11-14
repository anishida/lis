#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_solver_options.h"
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* The experimental factorization is deliberately real-only. */
#if defined(_COMPLEX) || defined(USE_COMPLEX)
int main(int argc,char **argv)
{
    (void)argc; (void)argv;
    fprintf(stderr,"testilutp: SKIP (complex LIS build)\n");
    return 77;
}
#else

/* Nonsymmetric, nonsingular matrix with a zero (0,0) diagonal:
   column pivoting is required for an unregularized exact ILUT factor. */
#define N 4
static const LIS_SCALAR a[N][N] = {
    {0.0, 4.0, 1.0, 0.0},
    {3.0, 0.0,-1.0, 0.0},
    {0.0, 2.0, 0.0, 5.0},
    {1.0, 0.0, 2.0, 1.0}
};
static const LIS_SCALAR truth[N] = {1.25,-0.5,2.0,0.75};

static int run_one(LIS_INT rcm, LIS_INT scale)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL, bt=NULL, x=NULL;
    LIS_INT err,i,j;
    LIS_REAL diff;
    int failed=1,swaps=0;
    char options[512];

    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto cleanup;
    err=lis_matrix_set_size(A,0,N);if(err)goto cleanup;
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            if(a[i][j]!=0.0)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,j,a[i][j],A);
                if(err)goto cleanup;
            }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto cleanup;
    err=lis_matrix_assemble(A);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&bt);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    for(i=0;i<N;i++)
    {
        b->value[i]=bt->value[i]=0.0;
        for(j=0;j<N;j++)
        {
            b->value[i]+=a[i][j]*truth[j];
            bt->value[i]+=a[j][i]*truth[j];
        }
    }

    err=lis_solver_create(&solver);if(err)goto cleanup;
    solver->A=A;
    snprintf(options,sizeof(options),
        "-i fgmres -p ilut -ilutp 1 -ilutp_threshold 1 "
        "-ilut_rcm %d -ilut_row_scale %d "
        "-ilut_drop 0 -ilut_rate 100 "
        "-ilu_pivot_tol 0 -scale none -omp_num_threads 1",
        (int)rcm,(int)scale);
    err=lis_solver_set_option(options,solver);if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);if(err){pc=NULL;goto cleanup;}
    if(!pc->ilutp_colperm || !pc->ilut_perm || !pc->D)
        goto cleanup;
    if(!rcm && !scale)
        for(i=0;i<N;i++)
            if(pc->ilutp_colperm[i]!=i) swaps++;
    if(!rcm && !scale && !swaps)
        goto cleanup;

    solver->precon=pc;
    err=lis_psolve_ilut_csr(solver,b,x);if(err)goto cleanup;
    for(i=0;i<N;i++)
    {
        diff=fabs(x->value[i]-truth[i]);
        if(!isfinite((double)diff) || diff>1e-9)
        {
            fprintf(stderr,
                "ILUTP forward mismatch rcm=%d scale=%d i=%d diff=%g\n",
                (int)rcm,(int)scale,(int)i,(double)diff);
            goto cleanup;
        }
    }

    err=lis_psolveh_ilut_csr(solver,bt,x);if(err)goto cleanup;
    for(i=0;i<N;i++)
    {
        diff=fabs(x->value[i]-truth[i]);
        if(!isfinite((double)diff) || diff>1e-9)
        {
            fprintf(stderr,
                "ILUTP adjoint mismatch rcm=%d scale=%d i=%d diff=%g\n",
                (int)rcm,(int)scale,(int)i,(double)diff);
            goto cleanup;
        }
    }
#ifdef _OPENMP
    if(!rcm && !scale)
    {
        /* The number of threads may change after factor construction.
         * Applying a serial ILUTP via an OpenMP-partitioned triangular
         * solve must fail explicitly rather than silently corrupting x. */
        omp_set_num_threads(2);
        err=lis_psolve_ilut_csr(solver,b,x);
        omp_set_num_threads(1);
        if(err!=LIS_ERR_NOT_IMPLEMENTED)goto cleanup;
        err=LIS_SUCCESS;
    }
#endif
    failed=0;

cleanup:
    if(failed)
        fprintf(stderr,"ILUTP test failed (rcm=%d scale=%d err=%d)\n",
            (int)rcm,(int)scale,(int)err);
    if(solver)solver->precon=NULL;
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(bt)lis_vector_destroy(bt);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}


/* A larger nonsymmetric CSR system exercises incomplete fill dropping and
 * different diagonal pivot strengths.  No perfect preconditioner is assumed. */
#define REG_N 32
static LIS_INT make_sparse_matrix(LIS_MATRIX *out,int weak_first)
{
    LIS_INT i,err;
    *out=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,out);if(err)return err;
    err=lis_matrix_set_size(*out,0,REG_N);if(err)return err;
    for(i=0;i<REG_N;i++)
    {
        LIS_SCALAR diag=(LIS_SCALAR)(4.0+0.03*(double)i);
        if(i==0 && weak_first)diag=(LIS_SCALAR)1e-8;
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,diag,*out);
        if(err)return err;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-0.65,*out);
            if(err)return err;
        }
        if(i+1<REG_N)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,0.7,*out);
            if(err)return err;
        }
        if(i+4<REG_N)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+4,0.25,*out);
            if(err)return err;
        }
        if(i>=5)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-5,-0.15,*out);
            if(err)return err;
        }
        if(i%4==0 && i+10<REG_N)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+10,0.38,*out);
            if(err)return err;
        }
    }
    err=lis_matrix_set_type(*out,LIS_MATRIX_CSR);
    if(err)return err;
    return lis_matrix_assemble(*out);
}

static int close_scalar(LIS_SCALAR x,LIS_SCALAR y)
{
    LIS_REAL scale=fabs(x)+fabs(y)+(LIS_REAL)1.0;
    return fabs(x-y)<=(LIS_REAL)2e-13*scale;
}

/* Exact tau=0 deliberately delegates to legacy ILUT. Tiny positive
 * tau exercises the new kernel. For tau>0, compare numerical L/U entries
 * by column with dropping disabled and unlimited fill: legacy sorts
 * retained coefficients in ascending magnitude, native ILUTP in descending
 * magnitude. Truncated factors are NOT required to be continuous at zero.
 * Check the underlying elimination separately from this policy difference. */
static int check_tau_zero(int rcm,int scale,double tau)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER base=NULL,pivot=NULL;
    LIS_PRECON pbase=NULL,ppivot=NULL;
    LIS_INT err,i,k,j;
    char options[512];
    int failed=1;

    err=make_sparse_matrix(&A,0);if(err)goto cleanup;
    err=lis_solver_create(&base);if(err)goto cleanup;
    err=lis_solver_create(&pivot);if(err)goto cleanup;
    base->A=pivot->A=A;
    snprintf(options,sizeof(options),
        "-p ilut -ilut_rcm %d -ilut_row_scale %d "
        "-ilut_drop %.17g -ilut_rate %.17g "
        "-ilu_pivot_tol 1e-9 -scale none -omp_num_threads 1",
        rcm,scale,tau>0 ? 0.0 : 0.02,tau>0 ? 100.0 : 1.1);
    err=lis_solver_set_option(options,base);if(err)goto cleanup;
    err=lis_solver_set_option(options,pivot);if(err)goto cleanup;
    snprintf(options,sizeof(options),"-ilutp 1 -ilutp_threshold %.17g",tau);
    err=lis_solver_set_option(options,pivot);if(err)goto cleanup;
    err=lis_precon_create(base,&pbase);if(err){pbase=NULL;goto cleanup;}
    err=lis_precon_create(pivot,&ppivot);if(err){ppivot=NULL;goto cleanup;}
    if(!ppivot->ilutp_colperm)goto cleanup;

    for(i=0;i<REG_N;i++)
    {
        if(ppivot->ilutp_colperm[i]!=i ||
           pbase->L->nnz[i]!=ppivot->L->nnz[i] ||
           pbase->U->nnz[i]!=ppivot->U->nnz[i] ||
           !close_scalar(pbase->D->value[i],ppivot->D->value[i]))
        {
            fprintf(stderr,"TAU_PARITY tau=%.17g row=%d perm=%d L=%d/%d U=%d/%d D=%.17g/%.17g\\n",
                    tau,(int)i,(int)ppivot->ilutp_colperm[i],
                    (int)pbase->L->nnz[i],(int)ppivot->L->nnz[i],
                    (int)pbase->U->nnz[i],(int)ppivot->U->nnz[i],
                    (double)pbase->D->value[i],(double)ppivot->D->value[i]);
            goto cleanup;
        }
        for(k=0;k<pbase->L->nnz[i];k++)
        {
            int match=0;
            for(j=0;j<ppivot->L->nnz[i];j++)
            {
                if(pbase->L->index[i][k]==ppivot->L->index[i][j])
                {
                    if(!close_scalar(pbase->L->value[i][k],
                                     ppivot->L->value[i][j]))goto cleanup;
                    if(tau==0.0 && k!=j)goto cleanup;
                    match=1;
                    break;
                }
            }
            if(!match)goto cleanup;
        }
        for(k=0;k<pbase->U->nnz[i];k++)
        {
            int match=0;
            for(j=0;j<ppivot->U->nnz[i];j++)
            {
                if(pbase->U->index[i][k]==ppivot->U->index[i][j])
                {
                    if(!close_scalar(pbase->U->value[i][k],
                                     ppivot->U->value[i][j]))goto cleanup;
                    if(tau==0.0 && k!=j)goto cleanup;
                    match=1;
                    break;
                }
            }
            if(!match)goto cleanup;
        }
    }
    failed=0;
cleanup:
    if(failed)fprintf(stderr,
        "ILUTP zero/near-zero factor parity failed (tau=%.17g rcm=%d scale=%d err=%d)\\n",
        tau,rcm,scale,(int)err);
    if(ppivot)lis_precon_destroy(ppivot);
    if(pbase)lis_precon_destroy(pbase);
    if(pivot)lis_solver_destroy(pivot);
    if(base)lis_solver_destroy(base);
    if(A)lis_matrix_destroy(A);
    return failed;
}

static int solution_close(LIS_VECTOR x,LIS_VECTOR truth)
{
    LIS_INT i;
    for(i=0;i<REG_N;i++)
        if(!isfinite((double)x->value[i]) ||
           fabs(x->value[i]-truth->value[i])>1e-6)
            return 0;
    return 1;
}

/* True FGMRES solve, limited L/U fill, tau in (0,1), automatic sequence
 * reuse and option-change invalidation (not just direct psolve). */
static int check_fgmres_sequence(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR truth=NULL,b=NULL,x=NULL;
    LIS_PRECON first=NULL;
    LIS_INT i,err,status;
    int failed=1;

    err=make_sparse_matrix(&A,1);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&truth);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    for(i=0;i<REG_N;i++)
        truth->value[i]=(LIS_SCALAR)(0.75+0.01*(double)i);
    err=lis_matvec(A,truth,b);if(err)goto cleanup;

    err=lis_solver_create(&solver);if(err)goto cleanup;
    err=lis_solver_set_option(
        "-i fgmres -restart 24 -maxiter 300 -tol 1e-10 "
        "-p ilut -ilutp 1 -ilutp_threshold 0.5 "
        "-ilut_rcm 1 -ilut_row_scale 1 "
        "-ilut_drop 0.02 -ilut_rate 1.1 "
        "-ilu_pivot_tol 1e-9 -scale none -omp_num_threads 1",
        solver);
    if(err)goto cleanup;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto cleanup;
    err=lis_solve(A,b,x,solver);if(err)goto cleanup;
    err=lis_solver_get_status(solver,&status);if(err || status!=LIS_SUCCESS)
        goto cleanup;
    if(!solution_close(x,truth))goto cleanup;

    err=lis_solver_set_sequence_policy(solver,LIS_SEQUENCE_POLICY_AUTO);
    if(err)goto cleanup;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto cleanup;
    err=lis_solve_sequence(A,b,x,solver);if(err)goto cleanup;
    if(solver->retcode!=LIS_SUCCESS || !solution_close(x,truth) ||
       solver->sequence_precon==NULL ||
       solver->sequence_precon->ilutp_colperm==NULL)
        goto cleanup;
    first=solver->sequence_precon;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto cleanup;
    err=lis_solve_sequence(A,b,x,solver);if(err)goto cleanup;
    if(solver->retcode!=LIS_SUCCESS || !solution_close(x,truth) ||
       solver->sequence_precon!=first)
        goto cleanup;

    err=lis_solver_set_option("-ilutp_threshold 0.25",solver);if(err)goto cleanup;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto cleanup;
    err=lis_solve_sequence(A,b,x,solver);if(err)goto cleanup;
    if(solver->retcode!=LIS_SUCCESS || !solution_close(x,truth) ||
       solver->sequence_precon==NULL ||
       fabs(solver->sequence_precon_ilutp_threshold-0.25)>1e-12)
        goto cleanup;

    /*
     * Fixed matrix graph, new numerical values, frozen-pattern policy:
     * ILUTP must reject the incompatible in-place frozen refactorization,
     * and the sequence manager must fall back to a fresh pivoted build.
     */
    err=lis_solver_set_sequence_precon_refresh_policy(
        solver,LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN);
    if(err)goto cleanup;
    {
        LIS_INT k,found=0;
        for(k=A->ptr[10];k<A->ptr[11];k++)
            if(A->index[k]==10)
            {
                A->value[k] *= (LIS_SCALAR)1.15;
                found=1;
                break;
            }
        if(!found)goto cleanup;
    }
    err=lis_matrix_mark_values_changed(A);if(err)goto cleanup;
    err=lis_matvec(A,truth,b);if(err)goto cleanup;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto cleanup;
    err=lis_solve_sequence(A,b,x,solver);if(err)goto cleanup;
    if(solver->retcode!=LIS_SUCCESS || !solution_close(x,truth) ||
       solver->sequence_precon==NULL ||
       solver->sequence_precon->ilutp_colperm==NULL)
        goto cleanup;

    failed=0;
cleanup:
    if(failed)fprintf(stderr,
        "ILUTP FGMRES/sequence failed (err=%d status=%d)\n",
        (int)err,(int)(solver?solver->retcode:-1));
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(truth)lis_vector_destroy(truth);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* tau>0 must exercise the new elimination kernel even when no
 * permutation is needed. Retain an explicit structural CSR zero to
 * cover the zero-multiplier / zero-fill short-circuit. */
static int check_no_swap_kernel(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR truth=NULL,b=NULL,x=NULL;
    LIS_INT err=0,i,k;
    int failed=1,found=0;

    err=make_sparse_matrix(&A,0);if(err)goto cleanup;
    for(k=A->ptr[10];k<A->ptr[11];k++)
        if(A->index[k]==9)
        {
            A->value[k]=0.0; /* deliberate stored CSR zero */
            found=1;
            break;
        }
    if(!found)goto cleanup;
    err=lis_vector_duplicate(A,&truth);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    for(i=0;i<REG_N;i++)truth->value[i]=(LIS_SCALAR)(1.0+0.02*i);
    err=lis_matvec(A,truth,b);if(err)goto cleanup;
    err=lis_solver_create(&solver);if(err)goto cleanup;
    solver->A=A;
    err=lis_solver_set_option(
        "-p ilut -ilutp 1 -ilutp_threshold 1e-6 "
        "-ilut_rcm 0 -ilut_row_scale 0 "
        "-ilut_drop 0 -ilut_rate 100 -ilu_pivot_tol 0 "
        "-scale none -omp_num_threads 1",solver);
    if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);if(err){pc=NULL;goto cleanup;}
    if(!pc->ilutp_colperm)goto cleanup;
    for(i=0;i<REG_N;i++)
        if(pc->ilutp_colperm[i]!=i)goto cleanup;
    solver->precon=pc;
    err=lis_psolve_ilut_csr(solver,b,x);if(err)goto cleanup;
    if(!solution_close(x,truth))goto cleanup;
    failed=0;
cleanup:
    if(failed)fprintf(stderr,
        "ILUTP nonswapping positive-tau kernel failed err=%d\n",(int)err);
    if(solver)solver->precon=NULL;
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(truth)lis_vector_destroy(truth);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* The option-profile replacement helpers must distinguish identical
 * profiles from changed preconditioner options, and must release the
 * sequence preconditioner only for a changed profile. */
static int check_option_replace(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR truth=NULL,b=NULL,x=NULL;
    LIS_PRECON saved=NULL;
    LIS_INT err=0,i,changed=LIS_FALSE;
    int failed=1,step=0;
    char options[512];
    const int rcms[]={0,1,1,1};
    const double taus[]={0.5,0.5,0.25,0.25};
    const int limits[]={0,0,0,5000};

    err=make_sparse_matrix(&A,0);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&truth);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    for(i=0;i<REG_N;i++)truth->value[i]=1.0+0.01*(double)i;
    err=lis_matvec(A,truth,b);if(err)goto cleanup;
    err=lis_solver_create(&solver);if(err)goto cleanup;
    err=lis_solver_set_sequence_policy(solver,LIS_SEQUENCE_POLICY_AUTO);
    if(err)goto cleanup;

    for(step=0;step<4;step++)
    {
        snprintf(options,sizeof(options),
            "-i fgmres -restart 24 -maxiter 300 -tol 1e-10 "
            "-p ilut -ilutp 1 -ilutp_threshold %.16g "
            "-ilut_rcm %d -ilut_row_scale 1 "
            "-ilutp_full_limit %d -ilut_drop 0.02 -ilut_rate 1.1 "
            "-ilu_pivot_tol 1e-9 -scale none -omp_num_threads 1",
            taus[step],rcms[step],limits[step]);
        if(step==0)
        {
            err=lis_solver_set_option_replace(options,solver);
            if(err)goto cleanup;
        }
        else
        {
            changed=LIS_FALSE;
            err=lis_solver_set_option_replace_if_changed(options,solver,&changed);
            if(err || changed!=LIS_TRUE || solver->sequence_precon!=NULL)
                goto cleanup;
        }
        if(solver->ilut_rcm!=rcms[step] ||
           solver->ilut_row_scale!=1 ||
           fabs(solver->ilutp_threshold-taus[step])>1e-12 ||
           solver->ilutp_full_limit!=limits[step])
            goto cleanup;
        err=lis_vector_set_all(0.0,x);if(err)goto cleanup;
        err=lis_solve_sequence(A,b,x,solver);if(err)goto cleanup;
        if(solver->retcode!=LIS_SUCCESS || !solution_close(x,truth) ||
           solver->sequence_precon==NULL)goto cleanup;
        saved=solver->sequence_precon;

        /* Same profile must be a no-op, including factor identity. */
        changed=LIS_TRUE;
        err=lis_solver_set_option_replace_if_changed(options,solver,&changed);
        if(err || changed!=LIS_FALSE || solver->sequence_precon!=saved)
            goto cleanup;
    }
    /* The same replace/sequence lifecycle also applies to legacy ILUT.
     * Changing only the pre-existing RCM or row-scale options must
     * invalidate the old factors; identical profiles must reuse them. */
    {
        const int legacy_rcm[]={0,1,1};
        const int legacy_scale[]={0,0,1};
        for(step=0;step<3;step++)
        {
            snprintf(options,sizeof(options),
                "-i fgmres -restart 24 -maxiter 300 -tol 1e-10 "
                "-p ilut -ilutp 0 -ilut_rcm %d -ilut_row_scale %d "
                "-ilut_drop 0.02 -ilut_rate 1.1 -ilu_pivot_tol 1e-9 "
                "-scale none -omp_num_threads 1",
                legacy_rcm[step],legacy_scale[step]);
            changed=LIS_FALSE;
            err=lis_solver_set_option_replace_if_changed(
                options,solver,&changed);
            if(err || changed!=LIS_TRUE || solver->sequence_precon!=NULL)
                goto cleanup;
            if(solver->ilutp_enabled!=0 ||
               solver->ilut_rcm!=legacy_rcm[step] ||
               solver->ilut_row_scale!=legacy_scale[step])
                goto cleanup;
            err=lis_vector_set_all(0.0,x);if(err)goto cleanup;
            err=lis_solve_sequence(A,b,x,solver);if(err)goto cleanup;
            if(solver->retcode!=LIS_SUCCESS ||
               !solution_close(x,truth) ||
               solver->sequence_precon==NULL)goto cleanup;
            saved=solver->sequence_precon;
            changed=LIS_TRUE;
            err=lis_solver_set_option_replace_if_changed(
                options,solver,&changed);
            if(err || changed!=LIS_FALSE ||
               solver->sequence_precon!=saved)goto cleanup;
        }
    }
    failed=0;
cleanup:
    if(failed)fprintf(stderr,
        "ILUTP option_replace lifecycle failed step=%d err=%d changed=%d\n",
        step,(int)err,(int)changed);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(truth)lis_vector_destroy(truth);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* Hybrid must preserve an opt-in ILUTP policy in its inner ILUT
 * preconditioner. The ordinary Hybrid/USER default remains untouched. */
static int check_hybrid_ilutp(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL,inner=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR truth=NULL,b=NULL,x=NULL;
    LIS_INT err=0,i;
    int failed=1;

    err=make_sparse_matrix(&A,0);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&truth);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    for(i=0;i<REG_N;i++)truth->value[i]=(LIS_SCALAR)(1.1+0.02*i);
    err=lis_matvec(A,truth,b);if(err)goto cleanup;
    err=lis_solver_create(&solver);if(err)goto cleanup;
    solver->A=A;
    err=lis_solver_set_option(
        "-i fgmres -p hybrid -hybrid_i gmres -hybrid_p ilut "
        "-hybrid_maxiter 40 -hybrid_restart 20 -hybrid_tol 1e-12 "
        "-ilutp 1 -ilutp_full 1 -ilutp_threshold 0.1 "
        "-ilutp_full_limit 5000 -ilut_rcm 1 -ilut_row_scale 1 "
        "-ilut_drop 0.02 -ilut_rate 1.1 -ilu_pivot_tol 0 "
        "-scale none -omp_num_threads 1",solver);
    if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);if(err){pc=NULL;goto cleanup;}
    inner=pc->solver;
    if(!inner || !inner->precon ||
       !inner->ilutp_enabled || !inner->ilutp_full ||
       !inner->ilut_rcm || !inner->ilut_row_scale ||
       fabs(inner->ilutp_threshold-0.1)>1e-12 ||
       inner->ilutp_full_limit!=5000 ||
       !inner->precon->ilutp_colperm)
        goto cleanup;
    solver->precon=pc;
    err=lis_psolve_hybrid(solver,b,x);
    if(err || !solution_close(x,truth))goto cleanup;

    solver->precon=NULL;
    lis_precon_destroy(pc);pc=NULL;

    /* Half-storage policy must propagate to inner ILUTP and fail closed. */
    err=lis_solver_set_precon_half(solver,LIS_PRECON_HALF_ON);
    if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);
    if(err!=LIS_ERR_NOT_IMPLEMENTED)goto cleanup;
    pc=NULL; /* LIS destroyed the failed nested object */
    err=lis_solver_set_precon_half(solver,LIS_PRECON_HALF_OFF);
    if(err)goto cleanup;

    /* An ILUTP request with a different Hybrid inner preconditioner
     * must fail instead of silently ignoring the option. */
    err=lis_solver_set_option("-hybrid_p none",solver);
    if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);
    if(err!=LIS_ERR_ILL_ARG)goto cleanup;
    pc=NULL;
    failed=0;
cleanup:
    if(failed)fprintf(stderr,"ILUTP Hybrid inner-LU failed err=%d\n",(int)err);
    if(solver)solver->precon=NULL;
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(truth)lis_vector_destroy(truth);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* Half storage is a legacy ILUT optimization: opt-in ILUTP must not
 * silently compress factors or fall back to a different algorithm. */
static int check_ilutp_half_guard(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err=0;
    int failed=1;
    err=make_sparse_matrix(&A,0);if(err)goto cleanup;
    err=lis_solver_create(&solver);if(err)goto cleanup;
    solver->A=A;
    err=lis_solver_set_option(
        "-p ilut -ilutp 1 -ilutp_threshold 0.5 "
        "-scale none -omp_num_threads 1",solver);
    if(err)goto cleanup;
    err=lis_solver_set_precon_half(solver,LIS_PRECON_HALF_ON);
    if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);
    if(err!=LIS_ERR_NOT_IMPLEMENTED)
        goto cleanup;
    pc=NULL; /* failed preconditioner already destroyed by LIS */
    err=lis_solver_set_option("-ilutp 0",solver);
    if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);
    if(err){pc=NULL;goto cleanup;}
    if(pc->L==NULL || pc->U==NULL || pc->D==NULL)goto cleanup;
    failed=0;
cleanup:
    if(failed)fprintf(stderr,"ILUTP half-precision guard failed err=%d\n",
                      (int)err);
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* Complete-fill mode is an in-library alternative to incomplete LU;
 * verify Ax and A^T x, a deliberately weak ILUTP control, and the
 * explicit factor-memory budget without relying on external solvers. */
static int check_full_fill(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR truth=NULL,b=NULL,bt=NULL,x=NULL,r=NULL;
    LIS_INT e,i,nnzl=0,nnzu=0;
    double weak=0.0,bnorm=0.0,rnorm=0.0,maxerr=0.0;
    LIS_INT weak_nnz=0;
    int failed=1;

    e=make_sparse_matrix(&A,0);if(e)goto done;
    e=lis_vector_duplicate(A,&truth);if(e)goto done;
    e=lis_vector_duplicate(A,&b);if(e)goto done;
    e=lis_vector_duplicate(A,&bt);if(e)goto done;
    e=lis_vector_duplicate(A,&x);if(e)goto done;
    e=lis_vector_duplicate(A,&r);if(e)goto done;
    for(i=0;i<REG_N;i++)truth->value[i]=1.0+0.3*sin((double)i*0.13);
    e=lis_matvec(A,truth,b);if(e)goto done;
    e=lis_matvech(A,truth,bt);if(e)goto done;
    e=lis_vector_nrm2(b,&bnorm);if(e || !(bnorm>0))goto done;
    e=lis_solver_create(&solver);if(e)goto done;
    solver->A=A;
    e=lis_solver_set_option(
        "-i fgmres -p ilut -ilutp 1 -ilutp_threshold 1 "
        "-ilut_drop 0.5 -ilut_rate 0.2 -ilu_pivot_tol 0 "
        "-scale none -omp_num_threads 1",solver);if(e)goto done;
    /* Incomplete control has the same pivot policy but tiny fill. */
    e=lis_precon_create(solver,&pc);if(e){pc=NULL;goto done;}
    solver->precon=pc;
    e=lis_psolve_ilut_csr(solver,b,x);if(e)goto done;
    e=lis_matvec(A,x,r);if(e)goto done;
    e=lis_vector_axpy(-1.0,b,r);if(e)goto done;
    e=lis_vector_nrm2(r,&rnorm);if(e)goto done;
    weak=(double)rnorm/(double)bnorm;
    if(!isfinite(weak) || weak<0.0)goto done;
    for(i=0;i<REG_N;i++)weak_nnz+=pc->L->nnz[i]+pc->U->nnz[i];
    solver->precon=NULL;
    lis_precon_destroy(pc);pc=NULL;

    e=lis_solver_set_option("-ilutp_full 1 -ilutp_full_limit 1",solver);
    if(e)goto done;
    e=lis_precon_create(solver,&pc);
    if(e!=LIS_OUT_OF_MEMORY)goto done;
    pc=NULL; /* create destroys failed partial factors */
    e=lis_solver_set_option("-ilutp_full_limit 5000",solver);if(e)goto done;
    e=lis_precon_create(solver,&pc);if(e){pc=NULL;goto done;}
    if(!solver->ilutp_full || pc->ilutp_colperm==NULL)goto done;
    for(i=0;i<REG_N;i++){nnzl+=pc->L->nnz[i];nnzu+=pc->U->nnz[i];}
    if(nnzl+nnzu<=weak_nnz || nnzl+nnzu>5000)goto done;
    solver->precon=pc;
    e=lis_psolve_ilut_csr(solver,b,x);if(e)goto done;
    for(i=0;i<REG_N;i++) {
        double d=fabs(x->value[i]-truth->value[i]);
        if(d>maxerr)maxerr=d;
    }
    if(maxerr>1e-10)goto done;
    /* Check Ax=b while x is still the forward solve result. */
    e=lis_matvec(A,x,r);if(e)goto done;
    e=lis_vector_axpy(-1.0,b,r);if(e)goto done;
    e=lis_vector_nrm2(r,&rnorm);if(e)goto done;
    if((double)rnorm/(double)bnorm>1e-11)goto done;
    e=lis_psolveh_ilut_csr(solver,bt,x);if(e)goto done;
    for(i=0;i<REG_N;i++) {
        double d=fabs(x->value[i]-truth->value[i]);
        if(d>maxerr)maxerr=d;
    }
    if(maxerr>1e-10)goto done;
    /* Transpose solve has its own physical A^T x = bt check. */
    e=lis_matvech(A,x,r);if(e)goto done;
    e=lis_vector_axpy(-1.0,bt,r);if(e)goto done;
    {
        LIS_REAL btnorm,transnorm;
        e=lis_vector_nrm2(bt,&btnorm);if(e || !(btnorm>0.0))goto done;
        e=lis_vector_nrm2(r,&transnorm);if(e)goto done;
        if((double)transnorm/(double)btnorm>1e-11)goto done;
    }
    printf("\nILUTP_FULL n=%d weak_relres=%.9g full_relres=%.9g "
           "maxerr=%.9g weak_nnz=%d factor_nnz=%d\n",(int)REG_N,
           weak,(double)rnorm/(double)bnorm,maxerr,
           (int)weak_nnz,(int)(nnzl+nnzu));
    failed=0;
done:
    if(failed)fprintf(stderr,"ILUTP_FULL FAIL err=%d weak=%g fullmaxerr=%g\n",
                       (int)e,weak,maxerr);
    if(solver)solver->precon=NULL;
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(r)lis_vector_destroy(r);
    if(x)lis_vector_destroy(x);
    if(bt)lis_vector_destroy(bt);
    if(b)lis_vector_destroy(b);
    if(truth)lis_vector_destroy(truth);
    if(A)lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
    int failed=0;
    LIS_INT err;
    LIS_SOLVER solver=NULL;
    err=lis_initialize(&argc,&argv);if(err)return 1;
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    failed|=run_one(0,0);
    failed|=run_one(1,0);
    failed|=run_one(0,1);
    failed|=run_one(1,1);
    failed|=check_tau_zero(0,0,0.0);
    failed|=check_tau_zero(1,0,0.0);
    failed|=check_tau_zero(0,1,0.0);
    failed|=check_tau_zero(1,1,0.0);
    /* Tiny positive tau takes the native kernel without swapping: compare
     * actual factors against legacy rather than convergence alone. */
    failed|=check_tau_zero(0,0,1e-300);
    failed|=check_tau_zero(1,0,1e-300);
    failed|=check_tau_zero(0,1,1e-300);
    failed|=check_tau_zero(1,1,1e-300);
    failed|=check_fgmres_sequence();
    failed|=check_no_swap_kernel();
    failed|=check_option_replace();
    failed|=check_hybrid_ilutp();
    failed|=check_ilutp_half_guard();
    failed|=check_full_fill();
    err=lis_solver_create(&solver);
    if(err)failed=1;
    else
    {
        if(fabs(solver->ilutp_threshold-0.5)>1e-12 || solver->ilutp_full)failed=1;
        if(lis_solver_set_option("-p ilut -ilutp 1 -ilutp_threshold 0.5",
                                  solver)!=LIS_SUCCESS ||
           !solver->ilutp_enabled ||
           fabs(solver->ilutp_threshold-0.5)>1e-12 ||
           lis_solver_set_option("-ilutp_threshold 2.0",solver)==LIS_SUCCESS ||
           lis_solver_set_option("-ilutp_full_limit nope",solver)==LIS_SUCCESS ||
           lis_solver_set_option("-ilutp_full 2",solver)==LIS_SUCCESS)
            failed=1;
        lis_solver_destroy(solver);
    }
    lis_finalize();
    if(!failed)printf("\ntestilutp: PASS\n");
    return failed;
}

#endif /* real-double only */
