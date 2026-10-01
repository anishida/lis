#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lis.h"
#include <stdio.h>

static int scalar_close(LIS_SCALAR got, LIS_SCALAR expected)
{
	const LIS_REAL tol = (LIS_REAL)1.0e-10;
#ifdef _COMPLEX
	LIS_REAL dr = (LIS_REAL)creal(got-expected);
	LIS_REAL di = (LIS_REAL)cimag(got-expected);
	if( dr<0 ) dr = -dr;
	if( di<0 ) di = -di;
	return dr<=tol && di<=tol;
#else
	LIS_REAL d = (LIS_REAL)(got-expected);
	if( d<0 ) d = -d;
	return d<=tol;
#endif
}

static void print_mismatch(const char *name, LIS_INT global_i,
                           LIS_SCALAR got, LIS_SCALAR expected)
{
#ifdef _COMPLEX
	fprintf(stderr,
	        "%s mismatch at global index %d: got=(%.15e,%.15e) expected=(%.15e,%.15e)\n",
	        name,(int)global_i,
	        (double)creal(got),(double)cimag(got),
	        (double)creal(expected),(double)cimag(expected));
#else
	fprintf(stderr,
	        "%s mismatch at global index %d: got=%.15e expected=%.15e\n",
	        name,(int)global_i,(double)got,(double)expected);
#endif
}

static int check_result(const char *name, LIS_VECTOR y, LIS_INT global_n,
                        LIS_SCALAR alpha, LIS_SCALAR beta, int hermitian)
{
	LIS_INT i, gi;
	LIS_SCALAR a, b, expected;
	LIS_SCALAR aa = hermitian ? conj(alpha) : alpha;
	LIS_SCALAR bb = hermitian ? conj(beta)  : beta;

	for(i=0;i<y->n;i++)
	{
		gi = y->is + i;
		a = (LIS_SCALAR)(gi + 1);
		b = (LIS_SCALAR)(global_n - gi);
		expected = aa*a + bb*b;
		if( !scalar_close(y->value[i],expected) )
		{
			print_mismatch(name,gi,y->value[i],expected);
			return 1;
		}
	}
	return 0;
}

int main(int argc, char *argv[])
{
	LIS_MATRIX A = NULL, B = NULL, C = NULL;
	LIS_VECTOR x = NULL, y = NULL;
	LIS_INT i, is, ie, err;
	LIS_INT global_n = 4;
	LIS_INT local_fail = 0, global_fail = 0;
	LIS_SCALAR alpha, beta;

	lis_initialize(&argc,&argv);

	err = lis_matrix_create(LIS_COMM_WORLD,&A);
	if( err ) goto api_fail;

	/*
	 * Keep the canonical 4x4 test for serial/2-rank/4-rank runs.
	 * If somebody launches more than four MPI ranks, enlarge the diagonal
	 * problem so lis_matrix_set_size() remains valid.
	 */
	if( A->nprocs>global_n ) global_n = A->nprocs;

	err = lis_matrix_set_size(A,0,global_n);
	if( err ) goto api_fail;
	err = lis_matrix_get_range(A,&is,&ie);
	if( err ) goto api_fail;
	for(i=is;i<ie;i++)
	{
		err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)(i+1),A);
		if( err ) goto api_fail;
	}
	err = lis_matrix_assemble(A);
	if( err ) goto api_fail;

	err = lis_matrix_create(LIS_COMM_WORLD,&B);
	if( err ) goto api_fail;
	err = lis_matrix_set_size(B,0,global_n);
	if( err ) goto api_fail;
	err = lis_matrix_get_range(B,&is,&ie);
	if( err ) goto api_fail;
	for(i=is;i<ie;i++)
	{
		err = lis_matrix_set_value(LIS_INS_VALUE,i,i,
		                           (LIS_SCALAR)(global_n-i),B);
		if( err ) goto api_fail;
	}
	err = lis_matrix_assemble(B);
	if( err ) goto api_fail;

	/* Baseline Stage 2A case: C = 2*A - B. */
	alpha = (LIS_SCALAR)2.0;
	beta  = (LIS_SCALAR)-1.0;
	err = lis_matrix_create_operator(alpha,A,beta,B,&C);
	if( err ) goto api_fail;

	err = lis_vector_duplicate(C,&x);
	if( err ) goto api_fail;
	err = lis_vector_duplicate(C,&y);
	if( err ) goto api_fail;
	err = lis_vector_set_all((LIS_SCALAR)1.0,x);
	if( err ) goto api_fail;

	err = lis_matvec(C,x,y);
	if( err ) goto api_fail;
	local_fail |= check_result("matvec",y,global_n,alpha,beta,0);

	err = lis_matvech(C,x,y);
	if( err ) goto api_fail;
	local_fail |= check_result("matvech",y,global_n,alpha,beta,1);

	/*
	 * Destroying C must not destroy A or B. Verify that A still works.
	 */
	lis_matrix_destroy(C);
	C = NULL;
	err = lis_matvec(A,x,y);
	if( err ) goto api_fail;
	for(i=0;i<y->n;i++)
	{
		LIS_INT gi = y->is+i;
		LIS_SCALAR expected = (LIS_SCALAR)(gi+1);
		if( !scalar_close(y->value[i],expected) )
		{
			print_mismatch("A after operator destroy",gi,y->value[i],expected);
			local_fail = 1;
			break;
		}
	}

#ifdef _COMPLEX
	/*
	 * Complex-specific check. Real A/B make it easy to isolate whether
	 * lis_matvech() correctly conjugates alpha and beta.
	 */
	alpha = (LIS_SCALAR)(2.0 + 1.0*I);
	beta  = (LIS_SCALAR)(-1.0 + 2.0*I);

	err = lis_matrix_create_operator(alpha,A,beta,B,&C);
	if( err ) goto api_fail;

	err = lis_matvec(C,x,y);
	if( err ) goto api_fail;
	local_fail |= check_result("complex matvec",y,global_n,alpha,beta,0);

	err = lis_matvech(C,x,y);
	if( err ) goto api_fail;
	local_fail |= check_result("complex matvech",y,global_n,alpha,beta,1);

	lis_matrix_destroy(C);
	C = NULL;
#endif


#ifdef USE_MPI
	/*
	 * Exercise real MPI halo exchange, not just MPI execution.
	 *
	 * Use two rows per rank so A needs one remote column per rank while
	 * B needs two. This specifically checks that the operator vector layout
	 * is large enough for both operands:
	 *
	 *   A = I + 2*P1
	 *   B = 3*I + 4*P1 + 5*P2
	 *   C = 2*A - B
	 *
	 * For x = 1, both C*x and C^H*x are identically -6.
	 */
	if( A->nprocs>1 )
	{
		LIS_MATRIX Ah = NULL, Bh = NULL, Ch = NULL;
		LIS_VECTOR xh = NULL, yh = NULL;
		LIS_INT halo_n = 2*A->nprocs;
		LIS_INT his, hie, j1, j2;
		int halo_api_fail = 0;

		err = lis_matrix_create(LIS_COMM_WORLD,&Ah);
		if( err ) halo_api_fail = 1;
		if( !halo_api_fail )
		{
			err = lis_matrix_set_size(Ah,0,halo_n);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			err = lis_matrix_get_range(Ah,&his,&hie);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			for(i=his;i<hie;i++)
			{
				j1 = (i+1)%halo_n;
				err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,Ah);
				if( err ) { halo_api_fail = 1; break; }
				err = lis_matrix_set_value(LIS_INS_VALUE,i,j1,(LIS_SCALAR)2.0,Ah);
				if( err ) { halo_api_fail = 1; break; }
			}
		}
		if( !halo_api_fail )
		{
			err = lis_matrix_assemble(Ah);
			if( err ) halo_api_fail = 1;
		}

		if( !halo_api_fail )
		{
			err = lis_matrix_create(LIS_COMM_WORLD,&Bh);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			err = lis_matrix_set_size(Bh,0,halo_n);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			err = lis_matrix_get_range(Bh,&his,&hie);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			for(i=his;i<hie;i++)
			{
				j1 = (i+1)%halo_n;
				j2 = (i+2)%halo_n;
				err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)3.0,Bh);
				if( err ) { halo_api_fail = 1; break; }
				err = lis_matrix_set_value(LIS_INS_VALUE,i,j1,(LIS_SCALAR)4.0,Bh);
				if( err ) { halo_api_fail = 1; break; }
				err = lis_matrix_set_value(LIS_INS_VALUE,i,j2,(LIS_SCALAR)5.0,Bh);
				if( err ) { halo_api_fail = 1; break; }
			}
		}
		if( !halo_api_fail )
		{
			err = lis_matrix_assemble(Bh);
			if( err ) halo_api_fail = 1;
		}

		if( !halo_api_fail && Bh->np<=Ah->np )
		{
			fprintf(stderr,
			        "MPI halo setup failed on rank %d: A->np=%d B->np=%d\n",
			        (int)A->my_rank,(int)Ah->np,(int)Bh->np);
			local_fail = 1;
		}

		if( !halo_api_fail )
		{
			err = lis_matrix_create_operator((LIS_SCALAR)2.0,Ah,
			                                 (LIS_SCALAR)-1.0,Bh,&Ch);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			err = lis_vector_duplicate(Ch,&xh);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			err = lis_vector_duplicate(Ch,&yh);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			err = lis_vector_set_all((LIS_SCALAR)1.0,xh);
			if( err ) halo_api_fail = 1;
		}

		if( !halo_api_fail )
		{
			err = lis_matvec(Ch,xh,yh);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			for(i=0;i<yh->n;i++)
			{
				if( !scalar_close(yh->value[i],(LIS_SCALAR)-6.0) )
				{
					print_mismatch("MPI halo matvec",yh->is+i,
					               yh->value[i],(LIS_SCALAR)-6.0);
					local_fail = 1;
					break;
				}
			}
		}

		if( !halo_api_fail )
		{
			err = lis_matvech(Ch,xh,yh);
			if( err ) halo_api_fail = 1;
		}
		if( !halo_api_fail )
		{
			for(i=0;i<yh->n;i++)
			{
				if( !scalar_close(yh->value[i],(LIS_SCALAR)-6.0) )
				{
					print_mismatch("MPI halo matvech",yh->is+i,
					               yh->value[i],(LIS_SCALAR)-6.0);
					local_fail = 1;
					break;
				}
			}
		}

		if( halo_api_fail )
		{
			fprintf(stderr,"MPI halo API failure on rank %d: %d\n",
			        (int)A->my_rank,(int)err);
			local_fail = 1;
		}

		if( yh ) lis_vector_destroy(yh);
		if( xh ) lis_vector_destroy(xh);
		if( Ch ) lis_matrix_destroy(Ch);
		if( Bh ) lis_matrix_destroy(Bh);
		if( Ah ) lis_matrix_destroy(Ah);
	}
#endif

#ifdef USE_MPI
	{
		int lf = (int)local_fail;
		int gf = 0;
		MPI_Allreduce(&lf,&gf,1,MPI_INT,MPI_MAX,LIS_COMM_WORLD);
		global_fail = (LIS_INT)gf;
	}
#else
	global_fail = local_fail;
#endif

	if( A->my_rank==0 )
	{
		if( global_n==4 )
		{
			printf("C = 2*A - B, x = [1,1,1,1] -> [-2,1,4,7]\n");
		}
		printf("LIS_MATRIX_OPERATOR %s\n",global_fail ? "FAILED" : "PASSED");
	}

	lis_vector_destroy(y);
	lis_vector_destroy(x);
	lis_matrix_destroy(B);
	lis_matrix_destroy(A);
	lis_finalize();
	return global_fail ? 1 : 0;

api_fail:
	fprintf(stderr,"LIS_MATRIX_OPERATOR API failure: %d\n",(int)err);
	if( C ) lis_matrix_destroy(C);
	if( y ) lis_vector_destroy(y);
	if( x ) lis_vector_destroy(x);
	if( B ) lis_matrix_destroy(B);
	if( A ) lis_matrix_destroy(A);
	lis_finalize();
	return 1;
}
