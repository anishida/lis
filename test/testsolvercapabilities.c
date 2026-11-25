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

static int check_mfamg_configuration(void)
{
    LIS_SOLVER s=NULL;
    LIS_VECTOR v=NULL;
    LIS_INT p,n,c,dim,enabled,eq,issues,err;
    LIS_INT labels[3]={0,1,1},different[3]={0,1,2};
    int failed=1;

    if(lis_solver_create(&s)) goto done;
    if(lis_solver_get_near_nullspace_dim(s,&dim) || dim!=0 ||
       lis_solver_get_mfamg_contact_configuration(s,&enabled,&n) || n!=0 ||
       lis_solver_compare_mfamg_contact_partition(s,0,NULL,&eq) || !eq)
        goto done;
    if(lis_solver_set_option("-i fgmres -p hybrid -hybrid_i gmres -hybrid_p mfamg_ilutp",s))
        goto done;
    p=n=c=-1;
    if(lis_solver_get_effective_preconditioner(s,&p,&n)!=LIS_ERR_NOT_IMPLEMENTED ||
       lis_solver_get_effective_capabilities(s,&c)!=LIS_ERR_NOT_IMPLEMENTED ||
       p!=-1 || n!=-1 || c!=-1)
        goto done;
    if(lis_solver_set_option("-i fgmres -p mfamg_ilutp -mfamg_contact_aware 1",s))
        goto done;
    if(lis_vector_create(LIS_COMM_WORLD,&v) ||
       lis_vector_set_size(v,0,3) ||
       lis_vector_set_all((LIS_SCALAR)1,v) ||
       lis_solver_set_near_nullspace(s,1,&v))
        goto done;
    if(lis_solver_get_near_nullspace_dim(s,&dim) || dim!=1 ||
       lis_solver_check_nullspace_configuration(s,&issues) || issues!=0)
        goto done;
    if(lis_solver_set_mfamg_contact_partition(s,3,labels) ||
       lis_solver_get_mfamg_contact_configuration(s,&enabled,&n) ||
       enabled!=LIS_TRUE || n!=3)
        goto done;
    labels[2]=2; /* LIS owns a copy, not the caller's array. */
    if(lis_solver_compare_mfamg_contact_partition(s,3,different,&eq) || eq ||
       lis_solver_compare_mfamg_contact_partition(s,2,labels,&eq) || eq)
        goto done;
    labels[2]=1;
    if(lis_solver_compare_mfamg_contact_partition(s,3,labels,&eq) || !eq)
        goto done;
    eq=99;
    err=lis_solver_compare_mfamg_contact_partition(s,3,NULL,&eq);
    if(err!=LIS_ERR_ILL_ARG || eq!=99) goto done;
    if(lis_solver_set_mfamg_contact_partition(s,0,NULL) ||
       lis_solver_compare_mfamg_contact_partition(s,0,NULL,&eq) || !eq ||
       lis_solver_clear_near_nullspace(s) ||
       lis_solver_get_near_nullspace_dim(s,&dim) || dim!=0)
        goto done;
    failed=0;
done:
    if(v) lis_vector_destroy(v);
    if(s) lis_solver_destroy(s);
    if(failed) fprintf(stderr,"MFAMG configuration query regression failed\n");
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

    /* The native composite requires near-nullspace even though
     * the outer preconditioner type is not standalone MFAMG. */
    failed |= check_query(
        "FGMRES+MFAMG_ILUTP",
        "-i fgmres -p mfamg_ilutp",
        LIS_PRECON_TYPE_MFAMG_ILUTP,
        LIS_FALSE,
        LIS_SOLVER_CAP_FLEXIBLE_PRECON |
        LIS_SOLVER_CAP_NEAR_NULLSPACE);

    failed |= check_query(
        "FGCRODR+MFAMG_ILUTP",
        "-i fgcrodr -p mfamg_ilutp",
        LIS_PRECON_TYPE_MFAMG_ILUTP,
        LIS_FALSE,
        LIS_SOLVER_CAP_FLEXIBLE_PRECON |
        LIS_SOLVER_CAP_RECYCLING |
        LIS_SOLVER_CAP_NEAR_NULLSPACE);

    /* NONE must remain a real leaf, not a special parsing case. */
    failed |= check_query(
        "CG+NONE",
        "-i cg -p none",
        LIS_PRECON_TYPE_NONE,
        LIS_FALSE,
        0);

    failed |= check_recursive_hybrid_rejected();
    failed |= check_mfamg_configuration();

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
