#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lis.h"
#include <stdio.h>

typedef struct
{
	LIS_INT n;
	LIS_SCALAR diag;
} OPERATOR_USER_DIAG;

static LIS_INT operator_user_matvec(void *ctx,
                                    const LIS_SCALAR *x,
                                    LIS_SCALAR *y)
{
	OPERATOR_USER_DIAG *A = (OPERATOR_USER_DIAG *)ctx;
	LIS_INT i;

	for(i=0;i<A->n;i++) y[i] = A->diag*x[i];
	return LIS_SUCCESS;
}

static LIS_INT operator_user_matvech(void *ctx,
                                     const LIS_SCALAR *x,
                                     LIS_SCALAR *y)
{
#ifdef _COMPLEX
	OPERATOR_USER_DIAG *A = (OPERATOR_USER_DIAG *)ctx;
	LIS_INT i;

	for(i=0;i<A->n;i++) y[i] = conj(A->diag)*x[i];
	return LIS_SUCCESS;
#else
	return operator_user_matvec(ctx,x,y);
#endif
}

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

static int check_affine_result(const char *name, LIS_VECTOR y,
                               LIS_SCALAR slope, LIS_SCALAR offset)
{
	LIS_INT i, gi;
	LIS_SCALAR expected;

	for(i=0;i<y->n;i++)
	{
		gi = y->is + i;
		expected = slope*(LIS_SCALAR)gi + offset;
		if( !scalar_close(y->value[i],expected) )
		{
			print_mismatch(name,gi,y->value[i],expected);
			return 1;
		}
	}
	return 0;
}

static LIS_INT check_matrix_linear(const char *name, LIS_MATRIX A,
                                   LIS_VECTOR x, LIS_VECTOR y,
                                   LIS_INT global_n,
                                   LIS_SCALAR alpha, LIS_SCALAR beta,
                                   LIS_INT *fail)
{
	LIS_INT err;

	err = lis_matvec(A,x,y);
	if( err ) return err;
	*fail |= check_result(name,y,global_n,alpha,beta,0);

	err = lis_matvech(A,x,y);
	if( err ) return err;
	*fail |= check_result(name,y,global_n,alpha,beta,1);

	return LIS_SUCCESS;
}

static LIS_INT check_matrix_affine(const char *name, LIS_MATRIX A,
                                   LIS_VECTOR x, LIS_VECTOR y,
                                   LIS_SCALAR slope, LIS_SCALAR offset,
                                   LIS_INT *fail)
{
	LIS_INT err;

	err = lis_matvec(A,x,y);
	if( err ) return err;
	*fail |= check_affine_result(name,y,slope,offset);

	err = lis_matvech(A,x,y);
	if( err ) return err;
	*fail |= check_affine_result(name,y,slope,offset);

	return LIS_SUCCESS;
}

static LIS_INT check_matvec_linear(const char *name, LIS_MATRIX A,
                                   LIS_VECTOR x, LIS_VECTOR y,
                                   LIS_INT global_n,
                                   LIS_SCALAR alpha, LIS_SCALAR beta,
                                   LIS_INT *fail)
{
	LIS_INT err;

	err = lis_matvec(A,x,y);
	if( err ) return err;
	*fail |= check_result(name,y,global_n,alpha,beta,0);
	return LIS_SUCCESS;
}

static int expect_operator_error(const char *name,
                                 LIS_MATRIX A, LIS_MATRIX B,
                                 LIS_INT expected)
{
	LIS_MATRIX C = A ? A : B;
	LIS_INT err;

	err = lis_matrix_create_operator((LIS_SCALAR)1.0,A,
	                                 (LIS_SCALAR)1.0,B,&C);
	if( err!=expected )
	{
		fprintf(stderr,"%s: expected error %d, got %d\n",
		        name,(int)expected,(int)err);
		if( err==LIS_SUCCESS && C!=NULL && C!=A && C!=B )
			lis_matrix_destroy(C);
		return 1;
	}
	if( C!=NULL )
	{
		fprintf(stderr,"%s: output matrix was not cleared on failure\n",name);
		return 1;
	}
	return 0;
}

static LIS_INT create_user_diag(LIS_INT global_n, LIS_SCALAR diag,
                                OPERATOR_USER_DIAG *ctx, LIS_MATRIX *U)
{
	LIS_INT err, is, ie;

	*U = NULL;
	err = lis_matrix_create(LIS_COMM_WORLD,U);
	if( err ) return err;
	err = lis_matrix_set_size(*U,0,global_n);
	if( err ) return err;
	err = lis_matrix_get_range(*U,&is,&ie);
	if( err ) return err;

	ctx->n = ie-is;
	ctx->diag = diag;
	return lis_matrix_set_user(*U,ctx,
	                           operator_user_matvec,
	                           operator_user_matvech);
}

static void test_validation(LIS_MATRIX A, LIS_MATRIX B, LIS_INT *fail)
{
	LIS_INT err, saved;

	err = lis_matrix_create_operator((LIS_SCALAR)1.0,A,
	                                 (LIS_SCALAR)1.0,B,NULL);
	if( err!=LIS_ERR_ILL_ARG )
	{
		fprintf(stderr,"NULL output pointer: expected error %d, got %d\n",
		        (int)LIS_ERR_ILL_ARG,(int)err);
		*fail = 1;
	}

	*fail |= expect_operator_error("NULL A",NULL,B,LIS_ERR_ILL_ARG);
	*fail |= expect_operator_error("NULL B",A,NULL,LIS_ERR_ILL_ARG);

	saved = B->gn;
	B->gn = saved + 1;
	*fail |= expect_operator_error("global size mismatch",A,B,LIS_ERR_ILL_ARG);
	B->gn = saved;

	saved = B->n;
	B->n = saved + 1;
	*fail |= expect_operator_error("local size mismatch",A,B,LIS_ERR_ILL_ARG);
	B->n = saved;

	saved = B->is;
	B->is = saved + 1;
	*fail |= expect_operator_error("local start mismatch",A,B,LIS_ERR_ILL_ARG);
	B->is = saved;

	saved = B->ie;
	B->ie = saved + 1;
	*fail |= expect_operator_error("local end mismatch",A,B,LIS_ERR_ILL_ARG);
	B->ie = saved;

	saved = B->nprocs;
	B->nprocs = saved + 1;
	*fail |= expect_operator_error("nprocs mismatch",A,B,LIS_ERR_ILL_ARG);
	B->nprocs = saved;

	if( A->ranges!=NULL && B->ranges!=NULL )
	{
		saved = B->ranges[0];
		B->ranges[0] = saved + 1;
		*fail |= expect_operator_error("ranges mismatch",A,B,LIS_ERR_ILL_ARG);
		B->ranges[0] = saved;
	}
	else
	{
		fprintf(stderr,"ranges mismatch test skipped: ranges[] unavailable\n");
	}

}

static void test_nested(LIS_MATRIX A, LIS_MATRIX B,
                        LIS_VECTOR x, LIS_VECTOR y,
                        LIS_INT global_n, LIS_INT *fail)
{
	LIS_MATRIX C = NULL, D = NULL, E = NULL;
	LIS_INT err = LIS_SUCCESS;

	err = lis_matrix_create_operator((LIS_SCALAR)2.0,A,
	                                 (LIS_SCALAR)-1.0,B,&C);
	if( err ) goto fail;
	err = lis_matrix_create_operator((LIS_SCALAR)0.5,C,
	                                 (LIS_SCALAR)3.0,B,&D);
	if( err ) goto fail;
	err = lis_matrix_create_operator((LIS_SCALAR)1.0,C,
	                                 (LIS_SCALAR)1.0,D,&E);
	if( err ) goto fail;

	if( C->operator_work==NULL || D->operator_work==NULL ||
	    E->operator_work==NULL )
	{
		fprintf(stderr,"nested operator_work is NULL\n");
		*fail = 1;
	}
	if( C->operator_work==D->operator_work ||
	    C->operator_work==E->operator_work ||
	    D->operator_work==E->operator_work )
	{
		fprintf(stderr,"nested operators share operator_work\n");
		*fail = 1;
	}

	err = check_matrix_linear("nested C",C,x,y,global_n,
	                          (LIS_SCALAR)2.0,(LIS_SCALAR)-1.0,fail);
	if( err ) goto fail;
	err = check_matrix_linear("nested D",D,x,y,global_n,
	                          (LIS_SCALAR)1.0,(LIS_SCALAR)2.5,fail);
	if( err ) goto fail;
	err = check_matrix_linear("nested E",E,x,y,global_n,
	                          (LIS_SCALAR)3.0,(LIS_SCALAR)1.5,fail);
	if( err ) goto fail;

	lis_matrix_destroy(E);
	E = NULL;
	err = check_matvec_linear("D after E destroy",D,x,y,global_n,
	                          (LIS_SCALAR)1.0,(LIS_SCALAR)2.5,fail);
	if( err ) goto fail;

	lis_matrix_destroy(D);
	D = NULL;
	err = check_matvec_linear("C after D destroy",C,x,y,global_n,
	                          (LIS_SCALAR)2.0,(LIS_SCALAR)-1.0,fail);
	if( err ) goto fail;

	lis_matrix_destroy(C);
	C = NULL;
	err = check_matvec_linear("A after nested destroy",A,x,y,global_n,
	                          (LIS_SCALAR)1.0,(LIS_SCALAR)0.0,fail);
	if( err ) goto fail;
	err = check_matvec_linear("B after nested destroy",B,x,y,global_n,
	                          (LIS_SCALAR)0.0,(LIS_SCALAR)1.0,fail);
	if( err ) goto fail;
	goto cleanup;

fail:
	fprintf(stderr,"nested operator API failure: %d\n",(int)err);
	*fail = 1;

cleanup:
	if( E ) lis_matrix_destroy(E);
	if( D ) lis_matrix_destroy(D);
	if( C ) lis_matrix_destroy(C);
}

static void test_stress(LIS_MATRIX A, LIS_MATRIX B,
                        LIS_VECTOR x, LIS_VECTOR y,
                        LIS_INT global_n, LIS_INT *fail)
{
	LIS_INT k, err = LIS_SUCCESS;
	int stress_failed = 0;

	for(k=0;k<256;k++)
	{
		LIS_MATRIX C = NULL, D = NULL, E = NULL;

		err = lis_matrix_create_operator((LIS_SCALAR)2.0,A,
		                                 (LIS_SCALAR)-1.0,B,&C);
		if( err ) stress_failed = 1;
		if( !stress_failed )
		{
			err = lis_matrix_create_operator((LIS_SCALAR)0.5,C,
			                                 (LIS_SCALAR)3.0,B,&D);
			if( err ) stress_failed = 1;
		}
		if( !stress_failed )
		{
			err = lis_matrix_create_operator((LIS_SCALAR)1.0,C,
			                                 (LIS_SCALAR)1.0,D,&E);
			if( err ) stress_failed = 1;
		}
		if( !stress_failed )
		{
			err = lis_matvec(E,x,y);
			if( err ) stress_failed = 1;
		}
		if( !stress_failed &&
		    check_result("stress nested matvec",y,global_n,
		                 (LIS_SCALAR)3.0,(LIS_SCALAR)1.5,0) )
		{
			stress_failed = 1;
			*fail = 1;
		}

		if( E ) lis_matrix_destroy(E);
		if( D ) lis_matrix_destroy(D);
		if( C ) lis_matrix_destroy(C);

		if( stress_failed )
		{
			fprintf(stderr,"operator create/destroy stress failed: %d\n",(int)err);
			*fail = 1;
			break;
		}
	}

	if( !stress_failed )
	{
		err = check_matvec_linear("A after operator stress",A,x,y,global_n,
		                          (LIS_SCALAR)1.0,(LIS_SCALAR)0.0,fail);
		if( err )
		{
			fprintf(stderr,"A after operator stress API failure: %d\n",(int)err);
			*fail = 1;
		}
	}
}

static void test_user_combinations(LIS_MATRIX A, LIS_MATRIX B,
                                   LIS_VECTOR x, LIS_VECTOR y,
                                   LIS_INT global_n, LIS_INT *fail)
{
	LIS_MATRIX U = NULL, V = NULL, T = NULL;
	OPERATOR_USER_DIAG uctx, vctx;
	LIS_INT err = LIS_SUCCESS;

	err = create_user_diag(global_n,(LIS_SCALAR)5.0,&uctx,&U);
	if( err ) goto fail;
	err = create_user_diag(global_n,(LIS_SCALAR)-2.0,&vctx,&V);
	if( err ) goto fail;

	/* USER + CSR: 2*U - B = i + (10-global_n). */
	err = lis_matrix_create_operator((LIS_SCALAR)2.0,U,
	                                 (LIS_SCALAR)-1.0,B,&T);
	if( err ) goto fail;
	err = check_matrix_affine("USER + CSR",T,x,y,
	                          (LIS_SCALAR)1.0,
	                          (LIS_SCALAR)(10-global_n),fail);
	if( err ) goto fail;
	lis_matrix_destroy(T);
	T = NULL;

	/* CSR + USER: 2*A - U = 2*i - 3. */
	err = lis_matrix_create_operator((LIS_SCALAR)2.0,A,
	                                 (LIS_SCALAR)-1.0,U,&T);
	if( err ) goto fail;
	err = check_matrix_affine("CSR + USER",T,x,y,
	                          (LIS_SCALAR)2.0,(LIS_SCALAR)-3.0,fail);
	if( err ) goto fail;
	lis_matrix_destroy(T);
	T = NULL;

	/* USER + USER: 3*U + 0.5*V = 14*I. */
	err = lis_matrix_create_operator((LIS_SCALAR)3.0,U,
	                                 (LIS_SCALAR)0.5,V,&T);
	if( err ) goto fail;
	err = check_matrix_affine("USER + USER",T,x,y,
	                          (LIS_SCALAR)0.0,(LIS_SCALAR)14.0,fail);
	if( err ) goto fail;
	goto cleanup;

fail:
	fprintf(stderr,"USER/CSR operator test API failure: %d\n",(int)err);
	*fail = 1;

cleanup:
	if( T ) lis_matrix_destroy(T);
	if( V ) lis_matrix_destroy(V);
	if( U ) lis_matrix_destroy(U);
}

#ifdef USE_MPI
static LIS_INT create_identity_matrix(LIS_Comm comm,
                                      LIS_INT local_n, LIS_INT global_n,
                                      LIS_MATRIX *A)
{
	LIS_INT err, i, is, ie;

	*A = NULL;
	err = lis_matrix_create(comm,A);
	if( err ) return err;
	err = lis_matrix_set_size(*A,local_n,global_n);
	if( err ) return err;
	err = lis_matrix_get_range(*A,&is,&ie);
	if( err ) return err;

	for(i=is;i<ie;i++)
	{
		err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,*A);
		if( err ) return err;
	}
	return lis_matrix_assemble(*A);
}

static void test_mpi_validation(LIS_MATRIX A, LIS_MATRIX B, LIS_INT *fail)
{
	LIS_MATRIX Ap = NULL, Bp = NULL, T = NULL;
	LIS_INT part_n, local_b, err = LIS_SUCCESS;
	MPI_Comm dup_comm = MPI_COMM_NULL, rev_comm = MPI_COMM_NULL;
	LIS_Comm saved_comm;
	int mpi_err, comm_result = -1;
	int reverse_key;

	if( A->nprocs<=1 ) return;

	part_n = 4*A->nprocs;
	local_b = 4;
	if( A->my_rank==0 ) local_b = 3;
	else if( A->my_rank==1 ) local_b = 5;

	err = create_identity_matrix(LIS_COMM_WORLD,4,part_n,&Ap);
	if( err ) goto partition_fail;
	err = create_identity_matrix(LIS_COMM_WORLD,local_b,part_n,&Bp);
	if( err ) goto partition_fail;
	*fail |= expect_operator_error("MPI real partition mismatch",
	                               Ap,Bp,LIS_ERR_ILL_ARG);
	goto partition_cleanup;

partition_fail:
	fprintf(stderr,"MPI partition setup failure on rank %d: %d\n",
	        (int)A->my_rank,(int)err);
	*fail = 1;

partition_cleanup:
	if( Bp ) lis_matrix_destroy(Bp);
	if( Ap ) lis_matrix_destroy(Ap);

	mpi_err = MPI_Comm_dup(A->comm,&dup_comm);
	if( mpi_err!=MPI_SUCCESS )
	{
		fprintf(stderr,"MPI_Comm_dup failed on rank %d\n",(int)A->my_rank);
		*fail = 1;
	}
	else
	{
		mpi_err = MPI_Comm_compare(A->comm,dup_comm,&comm_result);
		if( mpi_err!=MPI_SUCCESS || comm_result!=MPI_CONGRUENT )
		{
			fprintf(stderr,
			        "MPI congruent communicator setup failed on rank %d: compare=%d\n",
			        (int)A->my_rank,comm_result);
			*fail = 1;
		}
		else
		{
			saved_comm = B->comm;
			B->comm = dup_comm;
			err = lis_matrix_create_operator((LIS_SCALAR)1.0,A,
			                                 (LIS_SCALAR)1.0,B,&T);
			B->comm = saved_comm;
			if( err!=LIS_SUCCESS || T==NULL )
			{
				fprintf(stderr,
				        "MPI_CONGRUENT operator rejected on rank %d: %d\n",
				        (int)A->my_rank,(int)err);
				*fail = 1;
			}
		}
		if( T ) lis_matrix_destroy(T);
		T = NULL;
		MPI_Comm_free(&dup_comm);
	}

	reverse_key = (int)(A->nprocs - 1 - A->my_rank);
	mpi_err = MPI_Comm_split(A->comm,0,reverse_key,&rev_comm);
	if( mpi_err!=MPI_SUCCESS )
	{
		fprintf(stderr,"MPI_Comm_split failed on rank %d\n",(int)A->my_rank);
		*fail = 1;
		return;
	}

	mpi_err = MPI_Comm_compare(A->comm,rev_comm,&comm_result);
	if( mpi_err!=MPI_SUCCESS || comm_result!=MPI_SIMILAR )
	{
		fprintf(stderr,
		        "MPI similar communicator setup failed on rank %d: compare=%d\n",
		        (int)A->my_rank,comm_result);
		*fail = 1;
	}
	else
	{
		saved_comm = B->comm;
		B->comm = rev_comm;
		*fail |= expect_operator_error("MPI similar communicator",
		                               A,B,LIS_ERR_ILL_ARG);
		B->comm = saved_comm;
	}
	MPI_Comm_free(&rev_comm);
}

static void test_mpi_halo(LIS_MATRIX A, LIS_INT *fail)
{
	LIS_MATRIX Ah = NULL, Bh = NULL, Ch = NULL;
	LIS_VECTOR xh = NULL, yh = NULL;
	LIS_INT halo_n, i, is, ie, j1, j2;
	LIS_INT err = LIS_SUCCESS;

	if( A->nprocs<=1 ) return;
	halo_n = 2*A->nprocs;

	err = lis_matrix_create(LIS_COMM_WORLD,&Ah);
	if( err ) goto fail;
	err = lis_matrix_set_size(Ah,0,halo_n);
	if( err ) goto fail;
	err = lis_matrix_get_range(Ah,&is,&ie);
	if( err ) goto fail;
	for(i=is;i<ie;i++)
	{
		j1 = (i+1)%halo_n;
		err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,Ah);
		if( err ) goto fail;
		err = lis_matrix_set_value(LIS_INS_VALUE,i,j1,(LIS_SCALAR)2.0,Ah);
		if( err ) goto fail;
	}
	err = lis_matrix_assemble(Ah);
	if( err ) goto fail;

	err = lis_matrix_create(LIS_COMM_WORLD,&Bh);
	if( err ) goto fail;
	err = lis_matrix_set_size(Bh,0,halo_n);
	if( err ) goto fail;
	err = lis_matrix_get_range(Bh,&is,&ie);
	if( err ) goto fail;
	for(i=is;i<ie;i++)
	{
		j1 = (i+1)%halo_n;
		j2 = (i+2)%halo_n;
		err = lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)3.0,Bh);
		if( err ) goto fail;
		err = lis_matrix_set_value(LIS_INS_VALUE,i,j1,(LIS_SCALAR)4.0,Bh);
		if( err ) goto fail;
		err = lis_matrix_set_value(LIS_INS_VALUE,i,j2,(LIS_SCALAR)5.0,Bh);
		if( err ) goto fail;
	}
	err = lis_matrix_assemble(Bh);
	if( err ) goto fail;

	if( Bh->np<=Ah->np )
	{
		fprintf(stderr,"MPI halo setup failed on rank %d: A->np=%d B->np=%d\n",
		        (int)A->my_rank,(int)Ah->np,(int)Bh->np);
		*fail = 1;
	}

	err = lis_matrix_create_operator((LIS_SCALAR)2.0,Ah,
	                                 (LIS_SCALAR)-1.0,Bh,&Ch);
	if( err ) goto fail;
	err = lis_vector_duplicate(Ch,&xh);
	if( err ) goto fail;
	err = lis_vector_duplicate(Ch,&yh);
	if( err ) goto fail;
	err = lis_vector_set_all((LIS_SCALAR)1.0,xh);
	if( err ) goto fail;

	err = lis_matvec(Ch,xh,yh);
	if( err ) goto fail;
	*fail |= check_affine_result("MPI halo matvec",yh,
	                             (LIS_SCALAR)0.0,(LIS_SCALAR)-6.0);

	err = lis_matvech(Ch,xh,yh);
	if( err ) goto fail;
	*fail |= check_affine_result("MPI halo matvech",yh,
	                             (LIS_SCALAR)0.0,(LIS_SCALAR)-6.0);
	goto cleanup;

fail:
	fprintf(stderr,"MPI halo API failure on rank %d: %d\n",
	        (int)A->my_rank,(int)err);
	*fail = 1;

cleanup:
	if( yh ) lis_vector_destroy(yh);
	if( xh ) lis_vector_destroy(xh);
	if( Ch ) lis_matrix_destroy(Ch);
	if( Bh ) lis_matrix_destroy(Bh);
	if( Ah ) lis_matrix_destroy(Ah);
}
#endif

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

	test_validation(A,B,&local_fail);

	/* Baseline: C = 2*A - B. */
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

	err = check_matrix_linear("baseline",C,x,y,global_n,
	                          alpha,beta,&local_fail);
	if( err ) goto api_fail;

	lis_matrix_destroy(C);
	C = NULL;
	err = check_matvec_linear("A after operator destroy",A,x,y,global_n,
	                          (LIS_SCALAR)1.0,(LIS_SCALAR)0.0,&local_fail);
	if( err ) goto api_fail;

#ifdef _COMPLEX
	alpha = (LIS_SCALAR)(2.0 + 1.0*I);
	beta  = (LIS_SCALAR)(-1.0 + 2.0*I);
	err = lis_matrix_create_operator(alpha,A,beta,B,&C);
	if( err ) goto api_fail;
	err = check_matrix_linear("complex",C,x,y,global_n,
	                          alpha,beta,&local_fail);
	if( err ) goto api_fail;
	lis_matrix_destroy(C);
	C = NULL;
#endif

	test_nested(A,B,x,y,global_n,&local_fail);
	test_stress(A,B,x,y,global_n,&local_fail);
	test_user_combinations(A,B,x,y,global_n,&local_fail);

#ifdef USE_MPI
	test_mpi_validation(A,B,&local_fail);
	test_mpi_halo(A,&local_fail);
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
			printf("C = 2*A - B, x = [1,1,1,1] -> [-2,1,4,7]\n");
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
