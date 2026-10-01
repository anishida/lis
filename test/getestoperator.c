#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lis.h"

typedef struct
{
	LIS_INT n;
	LIS_SCALAR *diag;
} DIAG_USER;

static LIS_INT diag_matvec(void *user_data, const LIS_SCALAR *x, LIS_SCALAR *y)
{
	DIAG_USER *D = (DIAG_USER *)user_data;
	LIS_INT i;

	for(i=0;i<D->n;i++) y[i] = D->diag[i] * x[i];
	return LIS_SUCCESS;
}

static LIS_INT diag_matvech(void *user_data, const LIS_SCALAR *x, LIS_SCALAR *y)
{
	DIAG_USER *D = (DIAG_USER *)user_data;
	LIS_INT i;

	for(i=0;i<D->n;i++) y[i] = conj(D->diag[i]) * x[i];
	return LIS_SUCCESS;
}

static int scalar_close(LIS_SCALAR a, LIS_SCALAR b, LIS_REAL tol)
{
#ifdef _COMPLEX
	return cabs(a-b) <= tol;
#else
	LIS_REAL d = (LIS_REAL)(a-b);
	if( d<0 ) d = -d;
	return d <= tol;
#endif
}

static int inputs_unchanged(LIS_MATRIX A, LIS_MATRIX B,
                            LIS_VECTOR x, LIS_VECTOR y)
{
	LIS_INT i, gi;
	LIS_SCALAR expected;

	lis_vector_set_all((LIS_SCALAR)1.0,x);

	if( lis_matvec(A,x,y) ) return 1;
	for(i=0;i<y->n;i++)
	{
		gi = y->is+i;
		expected = (LIS_SCALAR)((gi+1)*(gi+2));
		if( !scalar_close(y->value[i],expected,(LIS_REAL)1.0e-12) )
			return 1;
	}

	if( lis_matvec(B,x,y) ) return 1;
	for(i=0;i<y->n;i++)
	{
		gi = y->is+i;
		expected = (LIS_SCALAR)(gi+1);
		if( !scalar_close(y->value[i],expected,(LIS_REAL)1.0e-12) )
			return 1;
	}

	return A->matrix_type!=LIS_MATRIX_USER ||
	       B->matrix_type!=LIS_MATRIX_USER;
}

static int run_case(LIS_MATRIX A, LIS_MATRIX B, LIS_VECTOR x,
                    char *options, LIS_SCALAR expected, const char *name)
{
	LIS_ESOLVER esolver = NULL;
	LIS_SCALAR evalue = 0.0;
	LIS_INT err, status;

	err = lis_esolver_create(&esolver);
	if( err ) return 1;

	err = lis_esolver_set_option(options,esolver);
	if( err )
	{
		lis_esolver_destroy(esolver);
		return 1;
	}

	err = lis_gesolve(A,B,x,&evalue,esolver);
	if( err )
	{
		fprintf(stderr,"%s returned %d\n",name,(int)err);
		lis_esolver_destroy(esolver);
		return 1;
	}

	lis_esolver_get_status(esolver,&status);
	if( status!=LIS_SUCCESS )
	{
		fprintf(stderr,"%s status %d\n",name,(int)status);
		lis_esolver_destroy(esolver);
		return 1;
	}

	if( !scalar_close(evalue,expected,(LIS_REAL)1.0e-7) )
	{
#ifdef _COMPLEX
		fprintf(stderr,"%s unexpected eigenvalue (%e,%e)\n",
		        name,(double)creal(evalue),(double)cimag(evalue));
#else
		fprintf(stderr,"%s unexpected eigenvalue %e\n",
		        name,(double)evalue);
#endif
		lis_esolver_destroy(esolver);
		return 1;
	}

	lis_esolver_destroy(esolver);
	return 0;
}

int main(int argc, char *argv[])
{
	LIS_MATRIX A = NULL, B = NULL;
	LIS_VECTOR x = NULL, y = NULL;
	DIAG_USER Ad = {0,NULL}, Bd = {0,NULL};
	LIS_INT i, is, ie, err;
	LIS_INT global_n = 4;
	int fail = 0;

	lis_initialize(&argc,&argv);

	err = lis_matrix_create(LIS_COMM_WORLD,&A);
	if( err ) goto api_fail;
	if( A->nprocs>global_n ) global_n = A->nprocs;

	err = lis_matrix_set_size(A,0,global_n);
	if( err ) goto api_fail;
	err = lis_matrix_get_range(A,&is,&ie);
	if( err ) goto api_fail;

	Ad.n = A->n;
	Ad.diag = (LIS_SCALAR *)malloc((size_t)Ad.n*sizeof(LIS_SCALAR));
	if( !Ad.diag ) goto api_fail;
	for(i=0;i<Ad.n;i++)
	{
		LIS_INT gi = is+i;
		Ad.diag[i] = (LIS_SCALAR)((gi+1)*(gi+2));
	}
	err = lis_matrix_set_user(A,&Ad,diag_matvec,diag_matvech);
	if( err ) goto api_fail;

	err = lis_matrix_create(LIS_COMM_WORLD,&B);
	if( err ) goto api_fail;
	err = lis_matrix_set_size(B,0,global_n);
	if( err ) goto api_fail;
	err = lis_matrix_get_range(B,&is,&ie);
	if( err ) goto api_fail;

	Bd.n = B->n;
	Bd.diag = (LIS_SCALAR *)malloc((size_t)Bd.n*sizeof(LIS_SCALAR));
	if( !Bd.diag ) goto api_fail;
	for(i=0;i<Bd.n;i++)
	{
		LIS_INT gi = is+i;
		Bd.diag[i] = (LIS_SCALAR)(gi+1);
	}
	err = lis_matrix_set_user(B,&Bd,diag_matvec,diag_matvech);
	if( err ) goto api_fail;

	err = lis_vector_duplicate(A,&x);
	if( err ) goto api_fail;
	err = lis_vector_duplicate(A,&y);
	if( err ) goto api_fail;

	/* Generalized eigenvalues are 2,3,...,global_n+1. */
	fail |= run_case(A,B,x,
	                 "-e gpi -shift 1.0 -etol 1.0e-10 -emaxiter 1000",
	                 (LIS_SCALAR)(global_n+1),"shifted GPI");
	fail |= inputs_unchanged(A,B,x,y);

	/* 3.6 makes lambda=4 the unique nearest eigenvalue for inverse iteration. */
	fail |= run_case(A,B,x,
	                 "-e gii -shift 3.6 -etol 1.0e-10 -emaxiter 1000",
	                 (LIS_SCALAR)4.0,"shifted GII");
	fail |= inputs_unchanged(A,B,x,y);

#ifdef USE_MPI
	{
		int local_fail = fail, global_fail = 0;
		MPI_Allreduce(&local_fail,&global_fail,1,MPI_INT,MPI_MAX,LIS_COMM_WORLD);
		fail = global_fail;
	}
#endif

	if( A->my_rank==0 )
		printf("shifted USER GPI/GII %s\n",fail ? "FAILED" : "PASSED");

	lis_vector_destroy(y);
	lis_vector_destroy(x);
	lis_matrix_destroy(B);
	lis_matrix_destroy(A);
	free(Bd.diag);
	free(Ad.diag);
	lis_finalize();
	return fail ? 1 : 0;

api_fail:
	fprintf(stderr,"getestoperator API failure: %d\n",(int)err);
	if( y ) lis_vector_destroy(y);
	if( x ) lis_vector_destroy(x);
	if( B ) lis_matrix_destroy(B);
	if( A ) lis_matrix_destroy(A);
	free(Bd.diag);
	free(Ad.diag);
	lis_finalize();
	return 1;
}
