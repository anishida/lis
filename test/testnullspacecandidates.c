/* Focused regression for generic nullspace candidate classification. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include "lis.h"
#include "lis_nullspace_candidates.h"

static int fail(const char *msg, int err)
{
    fprintf(stderr,"FAIL: %s (err=%d)\n",msg,err);
    return 1;
}

int main(int argc, char *argv[])
{
    LIS_MATRIX A = NULL;
    LIS_SOLVER solver = NULL;
    LIS_VECTOR cand[3] = {NULL,NULL,NULL};
    LIS_VECTOR right_only = NULL;
    LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS diag;
    LIS_INT err;
    int i;

    err = lis_initialize(&argc,&argv);
    if(err) return fail("lis_initialize",err);

    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) return fail("matrix_create",err);

    err = lis_matrix_set_size(A,0,3);
    if(err) return fail("matrix_set_size",err);

    err = lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) return fail("matrix_set_type",err);

    err = lis_matrix_set_value(LIS_INS_VALUE,0,0,0.0,A);
    if(err) return fail("set A00",err);
    err = lis_matrix_set_value(LIS_INS_VALUE,1,1,1.0e-7,A);
    if(err) return fail("set A11",err);
    err = lis_matrix_set_value(LIS_INS_VALUE,2,2,1.0,A);
    if(err) return fail("set A22",err);

    err = lis_matrix_assemble(A);
    if(err) return fail("matrix_assemble",err);

    for(i=0;i<3;i++)
    {
        err = lis_vector_duplicate(A,&cand[i]);
        if(err) return fail("vector_duplicate",err);
        err = lis_vector_set_all(0.0,cand[i]);
        if(err) return fail("vector_set_all",err);
        err = lis_vector_set_value(LIS_INS_VALUE,i,1.0,cand[i]);
        if(err) return fail("vector_set_value",err);
    }

    err = lis_solver_create(&solver);
    if(err) return fail("solver_create",err);

    err = lis_solver_classify_nullspace_candidates(
        A,solver,3,cand,NULL,
        1.0,
        1.0e-12,
        1.0e-5,
        LIS_NULLSPACE_CANDIDATE_BIND_NEAR,
        &diag);
    if(err) return fail("classify",err);

    if(diag.input_count!=3 ||
       diag.independent_count!=3 ||
       diag.exact_count!=1 ||
       diag.near_count!=1)
    {
        fprintf(stderr,
            "FAIL: diagnostics input=%d independent=%d exact=%d near=%d\n",
            (int)diag.input_count,
            (int)diag.independent_count,
            (int)diag.exact_count,
            (int)diag.near_count);
        return 1;
    }

    printf(
        "PASS: nullspace candidates input=%d independent=%d exact=%d near=%d\n",
        (int)diag.input_count,
        (int)diag.independent_count,
        (int)diag.exact_count,
        (int)diag.near_count);

    lis_solver_destroy(solver);
    solver = NULL;
    for(i=0;i<3;i++)
    {
        lis_vector_destroy(cand[i]);
        cand[i] = NULL;
    }
    lis_matrix_destroy(A);
    A = NULL;

    /*
     * Jacobi-rotation regression. The candidates e0+e1 and e0-e1 are
     * independent but their images under diag(0,1e-7) are collinear.
     * The classifier must rotate the candidate/image pairs and recover
     * one exact and one near-null direction.
     */
    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) return fail("mixed matrix_create",err);
    err = lis_matrix_set_size(A,0,2);
    if(err) return fail("mixed matrix_set_size",err);
    err = lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) return fail("mixed matrix_set_type",err);
    err = lis_matrix_set_value(LIS_INS_VALUE,0,0,0.0,A);
    if(err) return fail("mixed set A00",err);
    err = lis_matrix_set_value(LIS_INS_VALUE,1,1,1.0e-7,A);
    if(err) return fail("mixed set A11",err);
    err = lis_matrix_assemble(A);
    if(err) return fail("mixed matrix_assemble",err);

    for(i=0;i<2;i++)
    {
        err = lis_vector_duplicate(A,&cand[i]);
        if(err) return fail("mixed vector_duplicate",err);
        err = lis_vector_set_all(0.0,cand[i]);
        if(err) return fail("mixed vector_set_all",err);
        err = lis_vector_set_value(LIS_INS_VALUE,0,1.0,cand[i]);
        if(err) return fail("mixed set e0",err);
        err = lis_vector_set_value(
            LIS_INS_VALUE,1,(i==0)?1.0:-1.0,cand[i]);
        if(err) return fail("mixed set e1",err);
    }

    err = lis_solver_create(&solver);
    if(err) return fail("mixed solver_create",err);

    err = lis_solver_classify_nullspace_candidates(
        A,solver,2,cand,NULL,
        1.0,
        1.0e-12,
        1.0e-5,
        LIS_NULLSPACE_CANDIDATE_BIND_NEAR,
        &diag);
    if(err) return fail("mixed classify",err);

    if(diag.independent_count!=2 ||
       diag.exact_count!=1 ||
       diag.near_count!=1)
        return fail("mixed Jacobi classification",LIS_BREAKDOWN);

    printf("PASS: mixed-candidate Jacobi classification\n");

    lis_solver_destroy(solver);
    solver = NULL;
    for(i=0;i<2;i++)
    {
        lis_vector_destroy(cand[i]);
        cand[i] = NULL;
    }
    lis_matrix_destroy(A);
    A = NULL;

    /*
     * Nonsymmetric ownership regression:
     *
     *     A = [ 0 1 ]
     *         [ 0 1 ]
     *
     * e0 is an exact right-null vector because A*e0=0, but it is not
     * a left-null vector because A^H*e0=[0,1]^T. Requesting left-space
     * classification must not remove the valid right-null vector.
     */
    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) return fail("nonsym matrix_create",err);

    err = lis_matrix_set_size(A,0,2);
    if(err) return fail("nonsym matrix_set_size",err);

    err = lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) return fail("nonsym matrix_set_type",err);

    err = lis_matrix_set_value(LIS_INS_VALUE,0,1,1.0,A);
    if(err) return fail("nonsym set A01",err);
    err = lis_matrix_set_value(LIS_INS_VALUE,1,1,1.0,A);
    if(err) return fail("nonsym set A11",err);

    err = lis_matrix_assemble(A);
    if(err) return fail("nonsym matrix_assemble",err);

    err = lis_vector_duplicate(A,&right_only);
    if(err) return fail("nonsym vector_duplicate",err);

    err = lis_vector_set_all(0.0,right_only);
    if(err) return fail("nonsym vector_set_all",err);

    err = lis_vector_set_value(LIS_INS_VALUE,0,1.0,right_only);
    if(err) return fail("nonsym vector_set_value",err);

    err = lis_solver_create(&solver);
    if(err) return fail("nonsym solver_create",err);

    cand[0] = right_only;

    err = lis_solver_classify_nullspace_candidates(
        A,solver,1,cand,NULL,
        2.0,
        1.0e-12,
        1.0e-6,
        LIS_NULLSPACE_CANDIDATE_BIND_LEFT,
        &diag);
    if(err) return fail("nonsym classify",err);

    if(diag.exact_count!=1 || diag.left_exact_count!=0)
    {
        fprintf(stderr,
            "FAIL: right/left ownership exact=%d left_exact=%d\n",
            (int)diag.exact_count,
            (int)diag.left_exact_count);
        return 1;
    }

    printf(
        "PASS: independent right/left classification exact=%d left_exact=%d\n",
        (int)diag.exact_count,
        (int)diag.left_exact_count);

    lis_solver_destroy(solver);
    lis_vector_destroy(right_only);
    lis_matrix_destroy(A);
    lis_finalize();

    return 0;
}
