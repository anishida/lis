#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include "lis.h"
#include <math.h>

static int check_generation(
	LIS_MATRIX A,
	LIS_UNSIGNED_INT expected_structure,
	LIS_UNSIGNED_INT expected_value,
	const char *label)
{
	LIS_UNSIGNED_INT sg = 0;
	LIS_UNSIGNED_INT vg = 0;
	LIS_INT err;

	err = lis_matrix_get_generation(A,&sg,&vg);
	if( err )
	{
		fprintf(stderr,"%s: get_generation failed: %d\n",
		        label,(int)err);
		return 1;
	}

	if( sg!=expected_structure || vg!=expected_value )
	{
		fprintf(
			stderr,
			"%s: expected=(%llu,%llu) actual=(%llu,%llu)\n",
			label,
			(unsigned long long)expected_structure,
			(unsigned long long)expected_value,
			(unsigned long long)sg,
			(unsigned long long)vg);
		return 1;
	}

	printf(
		"%s structure=%llu values=%llu\n",
		label,
		(unsigned long long)sg,
		(unsigned long long)vg);

	return 0;
}


/*
 * Incremental RCO assembly regression.
 *
 * Unique entries force row-capacity growth.  Repeated ADD updates
 * then exercise duplicate lookup without extending the structure.
 */
static int
check_incremental_add_assembly(void)
{
    enum
    {
        GLOBAL_N = 256,
        WIDTH = 64,
        REPEATS = 4
    };

    LIS_MATRIX A;

    LIS_VECTOR x;
    LIS_VECTOR y;

    LIS_INT err;
    LIS_INT i;
    LIS_INT j;
    LIS_INT pass;
    LIS_INT row;
    LIS_INT col;

    LIS_REAL ynorm;
    LIS_REAL expected_row;
    LIS_REAL expected_norm;
    LIS_REAL relerr;


    A = NULL;
    x = NULL;
    y = NULL;

    ynorm = 0.0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto fail;


    err =
        lis_matrix_set_size(
            A,
            0,
            GLOBAL_N);

    if( err )
        goto fail;


    for(i=0;i<A->n;i++)
    {
        row =
            A->is
            +
            i;

        for(j=0;j<WIDTH;j++)
        {
            col =
                (
                    row
                    +
                    j
                )
                %
                GLOBAL_N;

            err =
                lis_matrix_set_value(
                    LIS_ADD_VALUE,
                    row,
                    col,
                    (LIS_SCALAR)1.0,
                    A);

            if( err )
                goto fail;
        }
    }


    for(pass=0;pass<REPEATS;pass++)
    {
        for(i=0;i<A->n;i++)
        {
            row =
                A->is
                +
                i;

            for(j=0;j<WIDTH;j++)
            {
                col =
                    (
                        row
                        +
                        j
                    )
                    %
                    GLOBAL_N;

                err =
                    lis_matrix_set_value(
                        LIS_ADD_VALUE,
                        row,
                        col,
                        (LIS_SCALAR)0.25,
                        A);

                if( err )
                    goto fail;
            }
        }
    }


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if( err )
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &y);

    if( err )
        goto fail;


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        x);


    err =
        lis_matvec(
            A,
            x,
            y);

    if( err )
        goto fail;


    err =
        lis_vector_nrm2(
            y,
            &ynorm);

    if( err )
        goto fail;


    expected_row =
        (LIS_REAL)WIDTH
        *
        (
            1.0
            +
            0.25
            *
            (LIS_REAL)REPEATS
        );

    expected_norm =
        expected_row
        *
        sqrt(
            (double)GLOBAL_N);


    relerr =
        fabs(
            (double)(
                ynorm
                -
                expected_norm))
        /
        (
            1.0
            +
            fabs(
                (double)expected_norm)
        );


    if( relerr>1.0e-12 )
    {
        fprintf(
            stderr,
            "incremental ADD assembly mismatch: "
            "norm=% .17e expected=% .17e rel=% .3e\n",
            (double)ynorm,
            (double)expected_norm,
            (double)relerr);

        goto fail;
    }


    if(y)
        lis_vector_destroy(y);

    if(x)
        lis_vector_destroy(x);

    if(A)
        lis_matrix_destroy(A);


    printf(
        "LIS_INCREMENTAL_ADD_ASSEMBLY PASSED "
        "n=%d width=%d repeats=%d rel=% .3e\n",
        GLOBAL_N,
        WIDTH,
        REPEATS,
        (double)relerr);

    return 0;


fail:

    if(y)
        lis_vector_destroy(y);

    if(x)
        lis_vector_destroy(x);

    if(A)
        lis_matrix_destroy(A);

    return 1;
}


int main(int argc, char **argv)
{
	LIS_MATRIX A = NULL;
	LIS_MATRIX B = NULL;
	LIS_UNSIGNED_INT g = 0;
	LIS_INT err;
	int failed = 0;

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

	err = lis_matrix_create(LIS_COMM_WORLD,&A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	failed |= check_generation(A,0,0,"INITIAL");

	err = lis_matrix_mark_values_changed(A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}
	failed |= check_generation(A,0,1,"VALUES1");

	err = lis_matrix_mark_values_changed(A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}
	failed |= check_generation(A,0,2,"VALUES2");

	err = lis_matrix_mark_structure_changed(A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}
	failed |= check_generation(A,1,3,"STRUCTURE1");

	err = lis_matrix_set_size(A,0,2);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_matrix_set_value(
		LIS_INS_VALUE,0,0,(LIS_SCALAR)2.0,A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_matrix_set_value(
		LIS_INS_VALUE,1,1,(LIS_SCALAR)3.0,A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_matrix_assemble(A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	failed |= check_generation(A,1,3,"AFTER_ASSEMBLY");

	err = lis_matrix_duplicate(A,&B);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	failed |= check_generation(B,0,0,"DUPLICATE_INITIAL");

	err = lis_matrix_mark_values_changed(B);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	failed |= check_generation(B,0,1,"DUPLICATE_CHANGED");
	failed |= check_generation(A,1,3,"ORIGINAL_UNCHANGED");

	err = lis_matrix_get_generation(A,&g,NULL);
	if( err || g!=1 )
	{
		fprintf(stderr,"structure-only generation query failed\n");
		failed = 1;
	}

	g = 0;

	err = lis_matrix_get_generation(A,NULL,&g);
	if( err || g!=3 )
	{
		fprintf(stderr,"value-only generation query failed\n");
		failed = 1;
	}

	err = lis_matrix_get_generation(A,NULL,NULL);
	if( err!=LIS_ERR_ILL_ARG )
	{
		fprintf(
			stderr,
			"NULL output test: expected=%d actual=%d\n",
			(int)LIS_ERR_ILL_ARG,
			(int)err);
		failed = 1;
	}

cleanup:

	if( B ) lis_matrix_destroy(B);
	if( A ) lis_matrix_destroy(A);

	printf(
		"LIS_MATRIX_GENERATION_S1 %s\n",
		failed ? "FAILED" : "PASSED");

	failed |=
	    check_incremental_add_assembly();

	lis_finalize();

	return failed ? 1 : 0;
}
