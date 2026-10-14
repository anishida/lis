#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include "lis.h"
#include "lis_vector.h"

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
		fprintf(stderr,"%s: get_generation failed\n",label);
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

	return 0;
}

static int check_external_vector_storage(void)
{
	LIS_VECTOR v = NULL;
	LIS_VECTOR duplicate = NULL;
	LIS_SCALAR values0[3] = {
		(LIS_SCALAR)1.0,
		(LIS_SCALAR)2.0,
		(LIS_SCALAR)3.0
	};
	LIS_SCALAR values1[3] = {
		(LIS_SCALAR)4.0,
		(LIS_SCALAR)5.0,
		(LIS_SCALAR)6.0
	};
	LIS_SCALAR got = 0;
	LIS_INT err;
	int failed = 0;

	err = lis_vector_create(LIS_COMM_WORLD,&v);
	if( err ) return 1;

	err = lis_vector_set_size(v,0,3);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_set_external(v,values0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	if( v->value!=values0 || !v->value_external )
	{
		fprintf(stderr,"initial external vector attachment failed\n");
		failed = 1;
		goto cleanup;
	}

	if( !v->is_destroy )
	{
		fprintf(stderr,"external vector changed destroy policy\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_get_value(v,1,&got);
	if( err || got!=(LIS_SCALAR)2.0 )
	{
		fprintf(stderr,"external vector read failed\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_set_value(
		LIS_INS_VALUE,1,(LIS_SCALAR)8.0,v);
	if( err || values0[1]!=(LIS_SCALAR)8.0 )
	{
		fprintf(stderr,"external vector write-through failed\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_set_external(v,values1);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	if( v->value!=values1 || !v->value_external )
	{
		fprintf(stderr,"external vector rebind failed\n");
		failed = 1;
		goto cleanup;
	}

	if( values0[0]!=(LIS_SCALAR)1.0 ||
	    values0[1]!=(LIS_SCALAR)8.0 ||
	    values0[2]!=(LIS_SCALAR)3.0 )
	{
		fprintf(stderr,"previous external storage changed\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(v,&duplicate);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	if( duplicate->value==values1 ||
	    duplicate->value_external )
	{
		fprintf(stderr,"duplicate inherited external storage\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_unset(v);
	if( err )
	{
		fprintf(stderr,"external vector unset failed\n");
		failed = 1;
		goto cleanup;
	}

	if( v->value!=NULL || v->value_external ||
	    v->status!=LIS_VECTOR_NULL )
	{
		fprintf(stderr,"external vector unset state is invalid\n");
		failed = 1;
		goto cleanup;
	}

	/* Borrowed storage must remain application-owned after unset. */
	values1[0] = (LIS_SCALAR)9.0;

	err = lis_vector_reuse(&v);
	if( err )
	{
		fprintf(stderr,"vector reuse after external unset failed\n");
		failed = 1;
		goto cleanup;
	}

	if( v->value==values1 || v->value_external || !v->is_copy )
	{
		fprintf(stderr,"vector reuse ownership transition failed\n");
		failed = 1;
		goto cleanup;
	}

	printf("VECTOR_UNSET_REUSE PASSED\n");

cleanup:

	if( duplicate ) lis_vector_destroy(duplicate);
	if( v ) lis_vector_destroy(v);

        /*
         * Verify destruction while application-owned storage is still
         * attached.  In particular, lis_vector_destroy() must not free
         * the borrowed array.
         */
        if( !failed )
        {
                LIS_VECTOR borrowed = NULL;
                LIS_INT destroy_err;

                err = lis_vector_create(
                        LIS_COMM_WORLD,&borrowed);

                if( !err )
                {
                        err = lis_vector_set_size(
                                borrowed,0,3);
                }

                if( !err )
                {
                        err = lis_vector_set_external(
                                borrowed,values1);
                }

                if( err )
                {
                        fprintf(
                                stderr,
                                "external vector destroy setup failed rc=%d\n",
                                (int)err);

                        if( borrowed )
                        {
                                lis_vector_destroy(borrowed);
                        }

                        failed = 1;
                }
                else
                {
                        destroy_err =
                                lis_vector_destroy(borrowed);
                        borrowed = NULL;

                        if( destroy_err )
                        {
                                fprintf(
                                        stderr,
                                        "external vector destroy failed rc=%d\n",
                                        (int)destroy_err);
                                failed = 1;
                        }
                        else
                        {
                                /*
                                 * The array remains application-owned and
                                 * writable after destroying the LIS object.
                                 */
                                values1[0] = (LIS_SCALAR)12.0;

                                if( values1[0]!=(LIS_SCALAR)12.0 )
                                {
                                        fprintf(
                                                stderr,
                                                "external vector ownership lost after destroy\n");
                                        failed = 1;
                                }
                                else
                                {
                                        printf(
                                                "EXTERNAL_VECTOR_DESTROY_OWNERSHIP PASSED\n");
                                }
                        }
                }
        }

	values0[0] = (LIS_SCALAR)11.0;
	values1[0] = (LIS_SCALAR)12.0;

	printf(
		"LIS_EXTERNAL_VECTOR_S3 %s\n",
		failed ? "FAILED" : "PASSED");

	return failed;
}

int main(int argc, char **argv)
{
	LIS_MATRIX A = NULL;
	LIS_VECTOR x = NULL;
	LIS_VECTOR y = NULL;

	/*
	 * Deliberately use stack-owned CSR storage.
	 * lis_matrix_destroy() must never attempt to free these arrays.
	 */
	LIS_INT ptr[3] = {0,2,4};
	LIS_INT index[4] = {0,1,0,1};

	LIS_SCALAR values0[4] = {
		(LIS_SCALAR)2.0,
		(LIS_SCALAR)1.0,
		(LIS_SCALAR)1.0,
		(LIS_SCALAR)3.0
	};

	LIS_SCALAR values1[4] = {
		(LIS_SCALAR)4.0,
		(LIS_SCALAR)1.0,
		(LIS_SCALAR)1.0,
		(LIS_SCALAR)5.0
	};

	LIS_SCALAR y0 = 0;
	LIS_SCALAR y1 = 0;

	LIS_INT err;
	int failed = 0;
	int i;

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

	err = lis_matrix_create(LIS_COMM_WORLD,&A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_matrix_set_size(A,0,2);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_matrix_set_csr_external(
		4,ptr,index,values0,A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	if( !A->csr_external )
	{
		fprintf(stderr,"external CSR ownership flag not set\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_matrix_assemble(A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	if( !A->csr_external )
	{
		fprintf(stderr,"external CSR ownership lost during assembly\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_generation(A,0,0,"INITIAL");

	err = lis_vector_duplicate(A,&x);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A,&y);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_set_all((LIS_SCALAR)1.0,x);

	err = lis_matvec(A,x,y);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_get_value(y,0,&y0);
	lis_vector_get_value(y,1,&y1);

	if( y0!=(LIS_SCALAR)3.0 || y1!=(LIS_SCALAR)4.0 )
	{
		fprintf(stderr,"initial external CSR matvec failed\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_matrix_set_csr_values_external(values1,A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	failed |= check_generation(A,0,1,"VALUES_UPDATE");

	if( A->value!=values1 )
	{
		fprintf(stderr,"external CSR value pointer was not replaced\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_matvec(A,x,y);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_get_value(y,0,&y0);
	lis_vector_get_value(y,1,&y1);

	if( y0!=(LIS_SCALAR)5.0 || y1!=(LIS_SCALAR)6.0 )
	{
		fprintf(stderr,"updated external CSR matvec failed\n");
		failed = 1;
		goto cleanup;
	}

	/*
	 * Exercise a sequence of related matrices with identical CSR
	 * structure and repeatedly replaced external numerical values.
	 */
	for( i=0; i<100; i++ )
	{
		LIS_SCALAR *next_values;
		LIS_SCALAR expected0;
		LIS_SCALAR expected1;

		if( (i & 1)==0 )
		{
			next_values = values0;
			expected0 = (LIS_SCALAR)3.0;
			expected1 = (LIS_SCALAR)4.0;
		}
		else
		{
			next_values = values1;
			expected0 = (LIS_SCALAR)5.0;
			expected1 = (LIS_SCALAR)6.0;
		}

		err = lis_matrix_set_csr_values_external(next_values,A);
		if( err )
		{
			fprintf(stderr,"sequence update failed at %d\n",i);
			failed = 1;
			goto cleanup;
		}

		if( A->ptr!=ptr || A->index!=index || A->value!=next_values )
		{
			fprintf(stderr,"CSR pointer changed at update %d\n",i);
			failed = 1;
			goto cleanup;
		}

		if( !A->csr_external )
		{
			fprintf(stderr,"external ownership lost at update %d\n",i);
			failed = 1;
			goto cleanup;
		}

		err = lis_matvec(A,x,y);
		if( err )
		{
			failed = 1;
			goto cleanup;
		}

		lis_vector_get_value(y,0,&y0);
		lis_vector_get_value(y,1,&y1);

		if( y0!=expected0 || y1!=expected1 )
		{
			fprintf(
				stderr,
				"sequence matvec failed at update %d\n",
				i);
			failed = 1;
			goto cleanup;
		}
	}

	failed |= check_generation(A,0,101,"SEQUENCE_100");

	if( failed )
		goto cleanup;

	printf("SEQUENCE_100 PASSED\n");

	/*
	 * Change numerical values in-place without changing the external
	 * value pointer.  The application explicitly marks the new
	 * numerical matrix generation.
	 */
	if( A->value!=values1 )
	{
		fprintf(stderr,"unexpected value pointer before in-place update\n");
		failed = 1;
		goto cleanup;
	}

	values1[0] = (LIS_SCALAR)6.0;
	values1[1] = (LIS_SCALAR)2.0;
	values1[2] = (LIS_SCALAR)3.0;
	values1[3] = (LIS_SCALAR)7.0;

	err = lis_matrix_mark_values_changed(A);
	if( err )
	{
		fprintf(stderr,"in-place generation mark failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_generation(A,0,102,"INPLACE_UPDATE");
	if( failed )
		goto cleanup;

	if( A->ptr!=ptr || A->index!=index || A->value!=values1 )
	{
		fprintf(stderr,"CSR pointer changed during in-place update\n");
		failed = 1;
		goto cleanup;
	}

	if( !A->csr_external )
	{
		fprintf(stderr,"external ownership lost during in-place update\n");
		failed = 1;
		goto cleanup;
	}

	err = lis_matvec(A,x,y);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_get_value(y,0,&y0);
	lis_vector_get_value(y,1,&y1);

	if( y0!=(LIS_SCALAR)8.0 || y1!=(LIS_SCALAR)10.0 )
	{
		fprintf(stderr,"in-place external CSR matvec failed\n");
		failed = 1;
		goto cleanup;
	}

	printf("INPLACE_UPDATE PASSED\n");

	if( check_external_vector_storage() )
	{
		failed = 1;
		goto cleanup;
	}

cleanup:

	if( y ) lis_vector_destroy(y);
	if( x ) lis_vector_destroy(x);

	if( A )
	{
		lis_matrix_destroy(A);
		A = NULL;
	}

	/*
	 * These writes are intentional.  The arrays still belong to this
	 * application after LIS matrix destruction.
	 */
	ptr[0] = 0;
	index[0] = 0;
	values1[0] = (LIS_SCALAR)7.0;

	printf(
		"LIS_EXTERNAL_CSR_S2 %s\n",
		failed ? "FAILED" : "PASSED");

	lis_finalize();

	return failed ? 1 : 0;
}
