#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_matrix_operator.h"
#include "lis_precon.h"
#include <stdio.h>

static LIS_INT build_matrix(
    LIS_INT n,
    LIS_SCALAR diag,
    LIS_SCALAR off,
    LIS_SCALAR weak,
    LIS_MATRIX *A)
{
    LIS_INT i,err;

    *A=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return err;
    err=lis_matrix_set_size(*A,0,n);
    if(err) return err;

    for(i=0;i<n;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,diag,*A);
        if(err) return err;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,off,*A);
            if(err) return err;
        }
        if(i<n-1)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,off,*A);
            if(err) return err;
        }
        if(i>1)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-2,weak,*A);
            if(err) return err;
        }
        if(i<n-2)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+2,weak,*A);
            if(err) return err;
        }
    }

    err=lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if(err) return err;
    return lis_matrix_assemble(*A);
}

static int solve_and_check(
    LIS_MATRIX C,
    LIS_SOLVER solver,
    LIS_VECTOR xtrue,
    LIS_VECTOR b,
    LIS_VECTOR x,
    const char *label)
{
    LIS_INT err,status=LIS_SUCCESS;
    LIS_REAL norm=0.0;

    err=lis_matvec(C,xtrue,b);
    if(err)
    {
        fprintf(stderr,"%s: matvec failed (%d)\n",label,(int)err);
        return 1;
    }

    err=lis_vector_set_all((LIS_SCALAR)0.0,x);
    if(err) return 1;

    err=lis_solve_sequence(C,b,x,solver);
    if(err)
    {
        fprintf(stderr,"%s: sequence solve API failed (%d)\n",label,(int)err);
        return 1;
    }

    err=lis_solver_get_status(solver,&status);
    if(err || status!=LIS_SUCCESS)
    {
        fprintf(stderr,"%s: solver status=%d query=%d\n",
                label,(int)status,(int)err);
        return 1;
    }

    err=lis_vector_axpy((LIS_SCALAR)-1.0,xtrue,x);
    if(err) return 1;
    err=lis_vector_nrm2(x,&norm);
    if(err) return 1;

    if(!(norm<1.0e-8))
    {
        fprintf(stderr,"%s: solution error %.17e\n",label,(double)norm);
        return 1;
    }

    return 0;
}

int main(int argc,char **argv)
{
    const LIS_INT n=24;
    LIS_MATRIX A=NULL,B=NULL,C=NULL;
    LIS_VECTOR xtrue=NULL,b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon_before=NULL;
    LIS_MATRIX_ILU L_before=NULL,U_before=NULL;
    LIS_VECTOR D_before=NULL;
    LIS_INT i,err,policy=-1;
    int failed=0;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=build_matrix(
        n,(LIS_SCALAR)5.0,(LIS_SCALAR)-1.0,
        (LIS_SCALAR)1.0e-7,&A);
    if(err) { failed=1; goto cleanup; }

    err=build_matrix(
        n,(LIS_SCALAR)1.0,(LIS_SCALAR)0.0,
        (LIS_SCALAR)0.0,&B);
    if(err) { failed=1; goto cleanup; }

    err=lis_matrix_create_operator(
        (LIS_SCALAR)1.0,A,(LIS_SCALAR)-0.20,B,&C);
    if(err) { failed=1; goto cleanup; }

    err=lis_vector_duplicate(C,&xtrue); if(err) { failed=1; goto cleanup; }
    err=lis_vector_duplicate(C,&b);     if(err) { failed=1; goto cleanup; }
    err=lis_vector_duplicate(C,&x);     if(err) { failed=1; goto cleanup; }

    for(i=0;i<n;i++)
        xtrue->value[i]=(LIS_SCALAR)(1.0 + 0.05*(double)i);

    err=lis_solver_create(&solver);
    if(err) { failed=1; goto cleanup; }

    err=lis_solver_get_sequence_precon_refresh_policy(solver,&policy);
    if(err || policy!=LIS_SEQUENCE_PRECON_REFRESH_SAFE)
    {
        fprintf(stderr,"default refresh policy mismatch (%d, %d)\n",
                (int)err,(int)policy);
        failed=1;
        goto cleanup;
    }

    err=lis_solver_set_option(
        "-i gmres -restart 20 -p ilut "
        "-ilut_drop 1.0e-3 -ilut_rate 2.0 "
        "-tol 1.0e-12 -maxiter 1000 "
        "-scale none -initx_zeros true",
        solver);
    if(err) { failed=1; goto cleanup; }

    fprintf(stderr,
        "configured precon=%d solver=%d scale=%d precision=%d\n",
        (int)solver->options[LIS_OPTIONS_PRECON],
        (int)solver->options[LIS_OPTIONS_SOLVER],
        (int)solver->options[LIS_OPTIONS_SCALE],
        (int)solver->options[LIS_OPTIONS_PRECISION]);

    err=lis_solver_set_sequence_policy(
        solver,LIS_SEQUENCE_POLICY_AUTO);
    if(err) { failed=1; goto cleanup; }

    err=lis_solver_set_sequence_precon_refresh_policy(
        solver,LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN);
    if(err) { failed=1; goto cleanup; }

    err=lis_solver_get_sequence_precon_refresh_policy(solver,&policy);
    if(err || policy!=LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN)
    {
        fprintf(stderr,"configured refresh policy mismatch\n");
        failed=1;
        goto cleanup;
    }

    if(solve_and_check(C,solver,xtrue,b,x,"sigma=0.20"))
    {
        failed=1;
        goto cleanup;
    }

    precon_before=solver->sequence_precon;
    fprintf(stderr,
        "after first solve: seq_precon=%p solver_precon=%p seq_type=%d age=%d rebuilds=%d\n",
        (void *)solver->sequence_precon,
        (void *)solver->precon,
        (int)solver->sequence_precon_type,
        (int)solver->sequence_precon_age,
        (int)solver->sequence_diag_precon_rebuilds);
    if(precon_before!=NULL)
    {
        fprintf(stderr,
            "precon details: type=%d L=%p U=%p D=%p half=%p\n",
            (int)precon_before->precon_type,
            (void *)precon_before->L,
            (void *)precon_before->U,
            (void *)precon_before->D,
            (void *)precon_before->half);
    }

    if(precon_before==NULL ||
       precon_before->precon_type!=LIS_PRECON_TYPE_ILUT ||
       precon_before->L==NULL || precon_before->U==NULL ||
       precon_before->D==NULL)
    {
        fprintf(stderr,"initial ILUT state missing\n");
        failed=1;
        goto cleanup;
    }
    L_before=precon_before->L;
    U_before=precon_before->U;
    D_before=precon_before->D;

    err=lis_matrix_operator_set_coefficients(
        C,(LIS_SCALAR)1.0,(LIS_SCALAR)-0.35);
    if(err) { failed=1; goto cleanup; }

    if(solve_and_check(C,solver,xtrue,b,x,"sigma=0.35"))
    {
        failed=1;
        goto cleanup;
    }

    if(solver->sequence_precon!=precon_before ||
       solver->sequence_precon->L!=L_before ||
       solver->sequence_precon->U!=U_before ||
       solver->sequence_precon->D!=D_before)
    {
        fprintf(stderr,"frozen ILUT symbolic state was replaced\n");
        failed=1;
        goto cleanup;
    }

    /* Invalid setter input must be transactional. */
    err=lis_solver_set_sequence_precon_refresh_policy(solver,99);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"invalid refresh policy accepted\n");
        failed=1;
        goto cleanup;
    }
    err=lis_solver_get_sequence_precon_refresh_policy(solver,&policy);
    if(err || policy!=LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN)
    {
        fprintf(stderr,"invalid setter changed active policy\n");
        failed=1;
        goto cleanup;
    }

    /* Unsupported half-storage refresh must not mutate the cache. */
    solver->precon_half_mode=LIS_PRECON_HALF_ON;
    err=lis_precon_update_ilut_frozen(solver,solver->sequence_precon);
    solver->precon_half_mode=LIS_PRECON_HALF_OFF;
    if(err!=LIS_ERR_NOT_IMPLEMENTED ||
       solver->sequence_precon!=precon_before)
    {
        fprintf(stderr,"half-storage frozen refresh was not rejected\n");
        failed=1;
        goto cleanup;
    }

    /* Zero diagonal in the updated operator must reject frozen factors. */
    err=lis_matrix_operator_set_coefficients(
        C,(LIS_SCALAR)1.0,(LIS_SCALAR)-5.0);
    if(err) { failed=1; goto cleanup; }
    err=lis_precon_update_ilut_frozen(solver,solver->sequence_precon);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"invalid frozen ILUT pivot was accepted\n");
        failed=1;
        goto cleanup;
    }

    err=lis_matrix_operator_set_coefficients(
        C,(LIS_SCALAR)1.0,(LIS_SCALAR)-0.35);
    if(err) { failed=1; goto cleanup; }
    /* A structural-generation change must force a complete rebuild. */
    err=lis_matrix_mark_structure_changed(C);
    if(err) { failed=1; goto cleanup; }
    if(solve_and_check(C,solver,xtrue,b,x,"structure-change"))
    {
        failed=1;
        goto cleanup;
    }

    printf("FROZEN_PATTERN_ILUT PASSED\n");

cleanup:
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(xtrue) lis_vector_destroy(xtrue);
    if(C) lis_matrix_destroy(C);
    if(B) lis_matrix_destroy(B);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed ? 1 : 0;
}
