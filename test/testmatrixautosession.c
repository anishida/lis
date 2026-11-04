#include "lis.h"
#include "lis_solver.h"
#include "lis_matrix_materialize.h"
#include "lis_solver_options.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _COMPLEX
#include <complex.h>
#endif

typedef struct {
    LIS_REAL scale;
    LIS_INT dense;
    LIS_INT failed_row;
    LIS_SOLVER observed_solver;
    LIS_INT assert_reset_before_row;
    /* Induce failure only after a CSR row has actually been overwritten,
     * not in either of the two pre-refresh pattern-validation passes. */
    LIS_MATRIX watched_csr;
    LIS_SCALAR watched_before;
    LIS_INT fail_after_partial_refresh;
    LIS_INT partial_failure_observed;
} TEST_CONTEXT;

static LIS_INT
test_row(void *data,LIS_INT row,LIS_INT *nnz,LIS_INT *cols,LIS_SCALAR *vals)
{
    const double base[9]={4,-1,0.25,-1,4,-1,0.25,-1,3};
    TEST_CONTEXT *ctx=(TEST_CONTEXT *)data;
    LIS_INT col,n=0,capacity;

    if(ctx==NULL || nnz==NULL || row<0 || row>=3)
        return LIS_ERR_ILL_ARG;
    if(ctx->assert_reset_before_row)
    {
        if(ctx->observed_solver==NULL ||
           ctx->observed_solver->sequence_precon!=NULL)
        {
            fprintf(stderr,"row callback saw stale preconditioner before CSR refresh\n");
            return LIS_ERR_ILL_ARG;
        }
        ctx->assert_reset_before_row=LIS_FALSE;
    }
    if(ctx->failed_row==row)
        return LIS_ERR_ILL_ARG;
    if(ctx->fail_after_partial_refresh && row==1 &&
       ctx->watched_csr!=NULL &&
       memcmp(&ctx->watched_csr->value[0],&ctx->watched_before,
              sizeof(LIS_SCALAR))!=0)
    {
        ctx->partial_failure_observed=LIS_TRUE;
        return LIS_ERR_ILL_ARG;
    }
    if((cols==NULL)!=(vals==NULL))
        return LIS_ERR_ILL_ARG;

    capacity=*nnz;
    for(col=0;col<3;col++)
    {
        if(!ctx->dense &&
           ((row==0 && col==2) || (row==2 && col==0)))
            continue;
        if(base[row*3+col]==0.0)
            continue;
        if(cols!=NULL)
        {
            if(n>=capacity)
                return LIS_ERR_ILL_ARG;
            cols[n]=col;
            vals[n]=(LIS_SCALAR)(base[row*3+col]*ctx->scale);
        }
        n++;
    }
    *nnz=n;
    return LIS_SUCCESS;
}

static LIS_INT
test_matvec(void *data,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    LIS_INT i,j,n,idx[3],err;
    LIS_SCALAR values[3];

    if(data==NULL || x==NULL || y==NULL)
        return LIS_ERR_ILL_ARG;
    for(i=0;i<3;i++)
    {
        n=3;
        err=test_row(data,i,&n,idx,values);
        if(err) return err;
        y[i]=0;
        for(j=0;j<n;j++)
            y[i]+=values[j]*x[idx[j]];
    }
    return LIS_SUCCESS;
}

static LIS_INT
test_make_user(TEST_CONTEXT *ctx,LIS_MATRIX *A)
{
    LIS_INT err;
    *A=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return err;
    err=lis_matrix_set_size(*A,0,3);
    if(err) return err;
    err=lis_matrix_set_user(*A,ctx,test_matvec,test_matvec);
    if(err) return err;
    err=lis_matrix_set_user_row(*A,test_row);
    if(err) return err;
    return lis_matrix_set_user_row_thread_safe(*A,LIS_TRUE);
}

static int
test_solution(LIS_VECTOR x,const LIS_SCALAR expected[3])
{
    LIS_INT i;
    for(i=0;i<3;i++)
    {
#ifdef _COMPLEX
        if(cabs(x->value[i]-expected[i])>1.0e-7)
#else
        if(fabs((double)(x->value[i]-expected[i]))>1.0e-7)
#endif
            return 1;
    }
    return 0;
}

static int
test_solve(LIS_MATRIX_AUTO_SESSION session,LIS_SOLVER solver,
           LIS_MATRIX source,LIS_VECTOR b,LIS_VECTOR x,
           const LIS_SCALAR truth[3],LIS_INT expected_policy,
           LIS_INT expected_action,const char *label)
{
    LIS_INT err,status,i,policy=99,action=99;

    for(i=0;i<3;i++)
        x->value[i]=(LIS_SCALAR)0;
    err=lis_matrix_auto_session_solve(
        session,source,b,x,&policy,&action);
    if(err || policy!=expected_policy || action!=expected_action)
    {
        fprintf(stderr,"%s: error=%lld policy=%lld action=%lld\n",
            label,(long long)err,(long long)policy,(long long)action);
        return 1;
    }
    err=lis_solver_get_status(solver,&status);
    if(err || status!=LIS_SUCCESS || test_solution(x,truth))
    {
        fprintf(stderr,"%s: inaccurate solve/status=%lld error=%lld\n",
            label,(long long)status,(long long)err);
        return 1;
    }
    printf("AUTO_SESSION %s PASS\n",label);
    return 0;
}

int main(int argc,char **argv)
{
    TEST_CONTEXT ctx={1.0,0,-1,NULL,0},ctx2={2.0,0,-1,NULL,0};
    const LIS_SCALAR truth[3]={
        (LIS_SCALAR)1,(LIS_SCALAR)-2,(LIS_SCALAR)3};
    LIS_MATRIX A=NULL,B=NULL;
    LIS_MATRIX_AUTO_SESSION session=NULL,jacobi_session=NULL;
    LIS_SOLVER solver=NULL,jacobi_solver=NULL;
    LIS_VECTOR x=NULL,b=NULL;
    LIS_PRECON saved_precon=NULL,saved_jacobi=NULL;
    LIS_INT err,policy,action,sequence_policy;
    LIS_SCALAR old_jacobi_diagonal;
    int failed=1;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=test_make_user(&ctx,&A);if(err)goto done;
    err=test_make_user(&ctx2,&B);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_matrix_auto_session_create(&session);if(err)goto done;
    err=lis_matrix_auto_session_get_solver(session,&solver);if(err)goto done;
    ctx.observed_solver=solver;
    ctx2.observed_solver=solver;
    err=lis_solver_get_sequence_policy(solver,&sequence_policy);
    if(err || sequence_policy!=LIS_SEQUENCE_POLICY_AUTO)goto done;

    {
        char options[]="-i gmres -p ilut -restart 10 -tol 1e-11 -maxiter 80 -scale none";
        err=lis_solver_set_option(options,solver);if(err)goto done;
    }
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_DIRECT,LIS_MATRIX_POLICY_ACTION_DIRECT,
          "direct-ilut"))goto done;

    /* Changing values in DIRECT mode invalidates the previous sequence. */
    ctx.scale=(LIS_REAL)1.25;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    ctx.assert_reset_before_row=LIS_TRUE;
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_DIRECT,LIS_MATRIX_POLICY_ACTION_DIRECT,
          "direct-value-update-reset"))goto done;
    if(ctx.assert_reset_before_row)goto done;

    {
        char options[]="-i gmres -p ssor -restart 10 -tol 1e-11 -maxiter 80 -scale none";
        err=lis_solver_set_option_replace(options,solver);if(err)goto done;
    }
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "csr-build"))goto done;

    saved_precon=solver->sequence_precon;
    if(saved_precon==NULL)goto done;
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REUSE,
          "csr-reuse"))goto done;
    if(solver->sequence_precon!=saved_precon)
    {
        fprintf(stderr,"unchanged CSR lost persistent preconditioner\n");
        goto done;
    }

    /*
     * A change to solver options alone must switch CACHED_CSR -> DIRECT
     * before the persistent state is reused; no source generation changes.
     * The old cache stays allocated and may be REUSEd when switching back.
     */
    {
        char direct_options[]="-p ilut";
        /* Incremental option change deliberately does NOT explicitly reset. */
        err=lis_solver_set_option(direct_options,solver);if(err)goto done;
    }
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_DIRECT,LIS_MATRIX_POLICY_ACTION_DIRECT,
          "options-csr-to-direct"))goto done;
    {
        char csr_options[]="-p ssor";
        err=lis_solver_set_option(csr_options,solver);if(err)goto done;
    }
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REUSE,
          "options-direct-to-csr"))goto done;

    /* Failure leaves policy/action unchanged, even with invalid inputs. */
    policy=1234;action=5678;
    err=lis_matrix_auto_session_solve(session,NULL,b,x,&policy,&action);
    if(err==LIS_SUCCESS || policy!=1234 || action!=5678)
    {
        fprintf(stderr,"AUTO_SESSION modified outputs on invalid arguments\n");
        goto done;
    }

    ctx.scale=(LIS_REAL)1.5;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    ctx.assert_reset_before_row=LIS_TRUE;
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_VALUES,
          "csr-refresh"))goto done;
    if(ctx.assert_reset_before_row)goto done;

    ctx.dense=1;
    err=lis_matrix_mark_structure_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    ctx.assert_reset_before_row=LIS_TRUE;
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "csr-structural-rebuild"))goto done;
    if(ctx.assert_reset_before_row)goto done;

    /* Switching sources must release old sequence references before
       replacing the cache-owned CSR. */
    err=test_matvec(&ctx2,truth,b->value);if(err)goto done;
    ctx2.assert_reset_before_row=LIS_TRUE;
    if(test_solve(session,solver,B,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "csr-source-switch"))goto done;
    if(ctx2.assert_reset_before_row)goto done;

    err=lis_matrix_auto_session_reset(session);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "explicit-reset"))goto done;

    /* On a failed refresh, subsequent materialization must rebuild and
       must not reuse the partially changed cache or stale preconditioner. */
    ctx.scale=(LIS_REAL)2.0;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    ctx.failed_row=1;
    policy=4321;action=8765;
    err=lis_matrix_auto_session_solve(
        session,A,b,x,&policy,&action);
    if(err==LIS_SUCCESS || policy!=4321 || action!=8765)
    {
        fprintf(stderr,"failed callback or incorrect error outputs\n");
        goto done;
    }
    ctx.failed_row=-1;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "failure-recovery"))goto done;

    /*
     * Repeated structure rebuilds exercise cache identity changes.
     * No pointer equality, including malloc address reuse, may cause REUSE.
     */
    {
        LIS_INT k;
        for(k=0;k<12;k++)
        {
            ctx.dense=k%2;
            err=lis_matrix_mark_structure_changed(A);if(err)goto done;
            err=test_matvec(&ctx,truth,b->value);if(err)goto done;
            if(test_solve(session,solver,A,b,x,truth,
                  LIS_MATRIX_POLICY_CACHED_CSR,
                  LIS_MATRIX_MATERIALIZE_REBUILD,
                  "csr-epoch-rebuild"))goto done;
        }
    }

    {
        char options[]="-i gmres -p ilut -restart 10 -tol 1e-11 -maxiter 80 -scale none";
        err=lis_solver_set_option_replace(options,solver);if(err)goto done;
    }
    if(test_solve(session,solver,A,b,x,truth,
          LIS_MATRIX_POLICY_DIRECT,LIS_MATRIX_POLICY_ACTION_DIRECT,
          "csr-to-direct"))goto done;

    /*
     * Opt-in value-only Jacobi refresh for cached CSR. A USER matrix with a
     * row callback but no diagonal callback selects cached CSR for Jacobi.
     * The preconditioner owns its diagonal; refresh must keep that object,
     * recompute its reciprocal values, and return a correct solution.
     */
    err=lis_matrix_auto_session_create(&jacobi_session);if(err)goto done;
    err=lis_matrix_auto_session_get_solver(jacobi_session,&jacobi_solver);
    if(err)goto done;
    {
        char options[]="-i gmres -p jacobi -restart 10 -tol 1e-11 -maxiter 80 -scale none";
        err=lis_solver_set_option(options,jacobi_solver);if(err)goto done;
    }
    err=lis_solver_set_sequence_precon_refresh_policy(
        jacobi_solver,LIS_SEQUENCE_PRECON_REFRESH_JACOBI_NUMERIC);
    if(err)goto done;

    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    if(test_solve(jacobi_session,jacobi_solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "jacobi-csr-build"))goto done;
    saved_jacobi=jacobi_solver->sequence_precon;
    if(saved_jacobi==NULL || saved_jacobi->D==NULL)goto done;
    old_jacobi_diagonal=saved_jacobi->D->value[0];

    ctx.scale+=(LIS_REAL)0.25;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    if(test_solve(jacobi_session,jacobi_solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_VALUES,
          "jacobi-csr-value-refresh"))goto done;
    if(jacobi_solver->sequence_precon!=saved_jacobi ||
       fabs(saved_jacobi->D->value[0]-old_jacobi_diagonal)<1.0e-12)
    {
        fprintf(stderr,"Jacobi CSR value refresh failed to reuse/update factors\n");
        goto done;
    }

    if(test_solve(jacobi_session,jacobi_solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REUSE,
          "jacobi-csr-unchanged-reuse"))goto done;

    /* A structural change must reset the old Jacobi before row access. */
    ctx.observed_solver=jacobi_solver;
    ctx.dense=!ctx.dense;
    err=lis_matrix_mark_structure_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    ctx.assert_reset_before_row=LIS_TRUE;
    if(test_solve(jacobi_session,jacobi_solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "jacobi-csr-structural-reset"))goto done;
    if(ctx.assert_reset_before_row)goto done;

    /* Even if the application incorrectly marks a structural modification
     * as values-only, the actual pattern check must detect it, reset the
     * persistent Jacobi and rebuild the cache before solving. */
    ctx.dense=!ctx.dense;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    if(test_solve(jacobi_session,jacobi_solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "jacobi-csr-underreported-structure"))goto done;

    /* A failed row callback must drop retained state and invalidate cache. */
    ctx.scale+=(LIS_REAL)0.5;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    ctx.failed_row=1;
    policy=4321;action=8765;
    err=lis_matrix_auto_session_solve(
        jacobi_session,A,b,x,&policy,&action);
    if(err==LIS_SUCCESS || policy!=4321 || action!=8765 ||
       jacobi_solver->sequence_precon!=NULL)
    {
        fprintf(stderr,"Jacobi failed refresh retained unsafe state\n");
        goto done;
    }
    ctx.failed_row=-1;
    if(test_solve(jacobi_session,jacobi_solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "jacobi-csr-failure-recovery"))goto done;

    /* Verify error handling after a *partial* CSR value write.  Row 0
     * changes in the value pass, then row 1 fails. The sole pattern check
     * must finish before this callback failure is injected. */
    ctx.scale+=(LIS_REAL)0.25;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=test_matvec(&ctx,truth,b->value);if(err)goto done;
    /* Inspect the cached CSR via the solver's published matrix pointer;
     * LIS_MATRIX_AUTO_SESSION is intentionally opaque to applications. */
    if(jacobi_solver->sequence_precon_matrix==NULL)goto done;
    ctx.watched_csr=jacobi_solver->sequence_precon_matrix;
    ctx.watched_before=ctx.watched_csr->value[0];
    ctx.partial_failure_observed=LIS_FALSE;
    ctx.fail_after_partial_refresh=LIS_TRUE;
    policy=4321;action=8765;
    err=lis_matrix_auto_session_solve(
        jacobi_session,A,b,x,&policy,&action);
    ctx.fail_after_partial_refresh=LIS_FALSE;
    if(err==LIS_SUCCESS || !ctx.partial_failure_observed ||
       policy!=4321 || action!=8765 ||
       jacobi_solver->sequence_precon!=NULL)
    {
        fprintf(stderr,"Jacobi partially updated CSR was not invalidated\\n");
        goto done;
    }
    ctx.watched_csr=NULL;
    if(test_solve(jacobi_session,jacobi_solver,A,b,x,truth,
          LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD,
          "jacobi-csr-partial-failure-recovery"))goto done;

    failed=0;
done:
    printf("AUTO_SESSION_REGRESSION=%s\n",failed?"FAIL":"PASS");
    lis_matrix_auto_session_destroy(jacobi_session);
    lis_matrix_auto_session_destroy(session);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
