#include "lis.h"
#include "lis_solver_capabilities.h"
#include <stdio.h>

#define ALL_CAPS ( \
    LIS_SOLVER_CAP_FLEXIBLE_PRECON | \
    LIS_SOLVER_CAP_RECYCLING | \
    LIS_SOLVER_CAP_NEAR_NULLSPACE | \
    LIS_SOLVER_CAP_NESTED_PRECON )

static int check_query(
    const char *label,
    const char *options,
    LIS_INT expected_precon,
    LIS_INT expected_nested,
    LIS_INT expected_caps)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err,precon,nested,caps;
    int failed = 0;

    err = lis_solver_create(&solver);
    if(err)
    {
        fprintf(stderr,"%s: solver create failed: %d\n",label,(int)err);
        return 1;
    }

    err = lis_solver_set_option((char *)options,solver);
    if(err)
    {
        fprintf(stderr,"%s: option parse failed: %d\n",label,(int)err);
        lis_solver_destroy(solver);
        return 1;
    }

    precon = -1;
    nested = -1;
    caps = -1;

    err = lis_solver_get_effective_preconditioner(
        solver,&precon,&nested);
    if(err || precon!=expected_precon || nested!=expected_nested)
    {
        fprintf(stderr,
            "%s: effective preconditioner mismatch: err=%d precon=%d nested=%d\n",
            label,(int)err,(int)precon,(int)nested);
        failed = 1;
    }

    err = lis_solver_get_effective_capabilities(solver,&caps);
    if(err || caps!=expected_caps)
    {
        fprintf(stderr,
            "%s: capability mismatch: err=%d caps=0x%x expected=0x%x\n",
            label,(int)err,(unsigned)caps,(unsigned)expected_caps);
        failed = 1;
    }

    lis_solver_destroy(solver);
    return failed;
}

static int check_recursive_hybrid_rejected(void)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err,precon,nested,caps;
    int failed = 0;

    err = lis_solver_create(&solver);
    if(err) return 1;

    err = lis_solver_set_option(
        "-i gmres -p hybrid -hybrid_i gmres -hybrid_p hybrid",
        solver);
    if(err)
    {
        lis_solver_destroy(solver);
        return 1;
    }

    precon = -1;
    nested = -1;
    err = lis_solver_get_effective_preconditioner(
        solver,&precon,&nested);
    if(err!=LIS_ERR_NOT_IMPLEMENTED)
    {
        fprintf(stderr,
            "recursive Hybrid effective-preconditioner query was not rejected: err=%d\n",
            (int)err);
        failed = 1;
    }

    caps = -1;
    err = lis_solver_get_effective_capabilities(solver,&caps);
    if(err!=LIS_ERR_NOT_IMPLEMENTED)
    {
        fprintf(stderr,
            "recursive Hybrid capability query was not rejected: err=%d\n",
            (int)err);
        failed = 1;
    }

    lis_solver_destroy(solver);
    return failed;
}

int main(int argc, char **argv)
{
    LIS_INT err;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if(err) return 1;

    /* Ordinary solver + ordinary preconditioner. */
    failed |= check_query(
        "GMRES+ILU",
        "-i gmres -p ilu",
        LIS_PRECON_TYPE_ILU,
        LIS_FALSE,
        0);

    /* Flexible Krylov method without recycling. */
    failed |= check_query(
        "FGMRES+NONE",
        "-i fgmres -p none",
        LIS_PRECON_TYPE_NONE,
        LIS_FALSE,
        LIS_SOLVER_CAP_FLEXIBLE_PRECON);

    /* Recycling without flexible preconditioning. */
    failed |= check_query(
        "GCRODR+NONE",
        "-i gcrodr -p none",
        LIS_PRECON_TYPE_NONE,
        LIS_FALSE,
        LIS_SOLVER_CAP_RECYCLING);

    /* Flexible recycling. */
    failed |= check_query(
        "FGCRODR+ILU",
        "-i fgcrodr -p ilu",
        LIS_PRECON_TYPE_ILU,
        LIS_FALSE,
        LIS_SOLVER_CAP_FLEXIBLE_PRECON |
        LIS_SOLVER_CAP_RECYCLING);

    /* One supported Hybrid layer with an MFAMG leaf. */
    failed |= check_query(
        "FGCRODR+HYBRID(MFAMG)",
        "-i fgcrodr -p hybrid -hybrid_i gmres -hybrid_p mfamg",
        LIS_PRECON_TYPE_MFAMG,
        LIS_TRUE,
        ALL_CAPS);

    /* NONE must remain a real leaf, not a special parsing case. */
    failed |= check_query(
        "CG+NONE",
        "-i cg -p none",
        LIS_PRECON_TYPE_NONE,
        LIS_FALSE,
        0);

    failed |= check_recursive_hybrid_rejected();

    if(lis_solver_get_effective_preconditioner(NULL,NULL,NULL)!=LIS_ERR_ILL_ARG)
        failed = 1;
    if(lis_solver_get_effective_capabilities(NULL,&err)!=LIS_ERR_ILL_ARG)
        failed = 1;

    lis_finalize();

    if(failed)
    {
        fprintf(stderr,"SOLVER_EFFECTIVE_CAPABILITIES FAILED\n");
        return 1;
    }

    printf("SOLVER_EFFECTIVE_CAPABILITIES PASSED\n");
    return 0;
}
