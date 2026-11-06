#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"
#include <stdio.h>

/*
 * Exercise the real BSR ILUC setup path. In particular, this covers
 * lis_matrix_malloc_bsr() and the caller's success/ownership path.
 * Allocator-failure injection is covered separately by testblockallocoom.py.
 */
int main(int argc, char *argv[])
{
    LIS_MATRIX csr = NULL, bsr = NULL;
    LIS_SOLVER solver = NULL;
    LIS_PRECON precon = NULL;
    LIS_INT err, i, is, ie;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;

    err = lis_matrix_create(LIS_COMM_WORLD,&csr);
    if( err ) goto fail;
    err = lis_matrix_set_size(csr,0,8);
    if( err ) goto fail;

    /* The BSR ILUC regression is a single-process test. */
    if( csr->nprocs!=1 )
    {
        printf("ILUC_BSR_SMOKE SKIPPED (requires one rank)\n");
        goto cleanup;
    }

    err = lis_matrix_get_range(csr,&is,&ie);
    if( err ) goto fail;
    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)8.0,csr);
        if( err ) goto fail;
        if( i>0 )
        {
            err = lis_matrix_set_value(LIS_INS_VALUE,i,i-1,(LIS_SCALAR)-1.0,csr);
            if( err ) goto fail;
        }
        if( i<7 )
        {
            err = lis_matrix_set_value(LIS_INS_VALUE,i,i+1,(LIS_SCALAR)-1.0,csr);
            if( err ) goto fail;
        }
    }
    err = lis_matrix_set_type(csr,LIS_MATRIX_CSR);
    if( err ) goto fail;
    err = lis_matrix_assemble(csr);
    if( err ) goto fail;

    err = lis_matrix_duplicate(csr,&bsr);
    if( err ) goto fail;
    err = lis_matrix_set_blocksize(bsr,2,2,NULL,NULL);
    if( err ) goto fail;
    err = lis_matrix_set_type(bsr,LIS_MATRIX_BSR);
    if( err ) goto fail;
    err = lis_matrix_convert(csr,bsr);
    if( err ) goto fail;
    if( bsr->matrix_type!=LIS_MATRIX_BSR || bsr->bnr!=2 )
    {
        fprintf(stderr,"BSR conversion did not retain 2x2 blocks\n");
        err = LIS_FAILS;
        goto fail;
    }

    err = lis_solver_create(&solver);
    if( err ) goto fail;
    solver->A = bsr;
    err = lis_solver_set_option(
        "-p iluc -storage bsr -storage_block 2 "
        "-iluc_drop 0 -iluc_rate 5 -scale none",
        solver);
    if( err ) goto fail;

    err = lis_precon_create(solver,&precon);
    if( err )
    {
        /*
         * lis_precon_create() destroys its partial object on error.
         * Do not free the dangling output pointer again.
         */
        precon = NULL;
        goto fail;
    }
    if( precon==NULL || precon->L==NULL ||
        precon->U==NULL || precon->WD==NULL )
    {
        fprintf(stderr,"BSR ILUC did not create L, U, and WD\n");
        err = LIS_FAILS;
        goto fail;
    }

    printf("ILUC_BSR_SMOKE PASSED\n");
    goto cleanup;

fail:
    fprintf(stderr,"ILUC_BSR_SMOKE FAILED (err=%d)\n",(int)err);
    failed = 1;

cleanup:
    if( precon ) lis_precon_destroy(precon);
    if( solver ) lis_solver_destroy(solver);
    if( bsr ) lis_matrix_destroy(bsr);
    if( csr ) lis_matrix_destroy(csr);
    lis_finalize();
    return failed;
}
