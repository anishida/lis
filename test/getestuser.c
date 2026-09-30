#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lis.h"

typedef struct
{
	LIS_INT n;
	LIS_INT is;
	LIS_INT which;
} USER_DIAG;

static LIS_REAL user_diag_value(USER_DIAG *ctx, LIS_INT global_i)
{
	static const LIS_REAL adiag[4] = {2.0, 6.0, 12.0, 20.0};
	static const LIS_REAL bdiag[4] = {1.0, 2.0, 3.0, 4.0};

	return ctx->which==0 ? adiag[global_i] : bdiag[global_i];
}

static void user_scale_real(LIS_SCALAR *dst,
                            const LIS_SCALAR *src,
                            LIS_REAL a)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
	(*dst)[0] = a * (*src)[0];
	(*dst)[1] = a * (*src)[1];
#else
	*dst = (LIS_SCALAR)a * (*src);
#endif
}

static LIS_REAL user_scalar_real(const LIS_SCALAR *value)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
	return (LIS_REAL)(*value)[0];
#elif defined(_COMPLEX)
	return (LIS_REAL)creal(*value);
#else
	return (LIS_REAL)(*value);
#endif
}

static LIS_REAL user_scalar_error(const LIS_SCALAR *value, LIS_REAL expected)
{
	LIS_REAL dr, di;

#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
	dr = (LIS_REAL)(*value)[0] - expected;
	di = (LIS_REAL)(*value)[1];
#elif defined(_COMPLEX)
	dr = (LIS_REAL)creal(*value) - expected;
	di = (LIS_REAL)cimag(*value);
#else
	dr = (LIS_REAL)(*value) - expected;
	di = 0.0;
#endif

	if( dr < 0.0 ) dr = -dr;
	if( di < 0.0 ) di = -di;
	return dr > di ? dr : di;
}

static LIS_INT user_diag_matvec(void *user_data,
                                const LIS_SCALAR *x,
                                LIS_SCALAR *y)
{
	USER_DIAG *ctx = (USER_DIAG *)user_data;
	LIS_INT i;

	for(i=0;i<ctx->n;i++)
	{
		LIS_REAL d = user_diag_value(ctx,ctx->is+i);
		user_scale_real(&y[i],&x[i],d);
	}
	return LIS_SUCCESS;
}

static LIS_INT user_diag_matvech(void *user_data,
                                 const LIS_SCALAR *x,
                                 LIS_SCALAR *y)
{
	return user_diag_matvec(user_data,x,y);
}

static LIS_INT create_user_diag(LIS_Comm comm,
                                LIS_INT global_n,
                                LIS_INT which,
                                USER_DIAG *ctx,
                                LIS_MATRIX *A)
{
	LIS_INT err, local_n, is, ie, gn;

	err = lis_matrix_create(comm,A);
	if( err ) return err;

	err = lis_matrix_set_size(*A,0,global_n);
	if( err ) goto fail;

	err = lis_matrix_get_size(*A,&local_n,&gn);
	if( err ) goto fail;

	err = lis_matrix_get_range(*A,&is,&ie);
	if( err ) goto fail;

	ctx->n = local_n;
	ctx->is = is;
	ctx->which = which;
	(void)ie;
	(void)gn;

	err = lis_matrix_set_user(*A,ctx,user_diag_matvec,user_diag_matvech);
	if( err ) goto fail;

	return LIS_SUCCESS;

fail:
	lis_matrix_destroy(*A);
	*A = NULL;
	return err;
}

static LIS_INT run_case(LIS_Comm comm,
                        LIS_MATRIX A,
                        LIS_MATRIX B,
                        LIS_VECTOR x,
                        const char *method,
                        LIS_REAL expected)
{
	LIS_ESOLVER esolver;
	LIS_SCALAR evalue;
	LIS_REAL actual, error;
	LIS_INT err;
	char option[128];

	err = lis_esolver_create(&esolver);
	if( err ) return err;

	snprintf(option,sizeof(option),
	         "-e %s -etol 1.0e-10 -emaxiter 1000 -eprint mem",
	         method);

	err = lis_esolver_set_option(option,esolver);
	if( err ) goto out;

	err = lis_gesolve(A,B,x,&evalue,esolver);
	if( err )
	{
		fprintf(stderr,"%s: lis_gesolve failed: %d\n",method,(int)err);
		goto out;
	}

	actual = user_scalar_real(&evalue);
	error = user_scalar_error(&evalue,expected);

	lis_printf(comm,
	           "LIS_MATRIX_USER generalized %s eigenvalue = %.15e, error = %.3e\n",
	           method,(double)actual,(double)error);

	if( error > (LIS_REAL)1.0e-7 )
	{
		fprintf(stderr,"%s: eigenvalue verification failed\n",method);
		err = LIS_FAILS;
	}

out:
	lis_esolver_destroy(esolver);
	return err;
}

LIS_INT main(int argc, char *argv[])
{
	LIS_Comm comm;
	LIS_MATRIX A = NULL, B = NULL;
	LIS_VECTOR x = NULL;
	USER_DIAG actx, bctx;
	LIS_INT err;
	const LIS_INT global_n = 4;

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

	comm = LIS_COMM_WORLD;

	err = create_user_diag(comm,global_n,0,&actx,&A);
	if( err ) goto fail;

	err = create_user_diag(comm,global_n,1,&bctx,&B);
	if( err ) goto fail;

	err = lis_vector_duplicate(A,&x);
	if( err ) goto fail;

	err = run_case(comm,A,B,x,"gpi",(LIS_REAL)5.0);
	if( err ) goto fail;

	err = run_case(comm,A,B,x,"gii",(LIS_REAL)2.0);
	if( err ) goto fail;

	lis_printf(comm,"LIS_MATRIX_USER generalized eigensolver test PASSED\n");

	lis_vector_destroy(x);
	lis_matrix_destroy(A);
	lis_matrix_destroy(B);
	lis_finalize();
	return 0;

fail:
	fprintf(stderr,"getestuser failed: %d\n",(int)err);
	if( x ) lis_vector_destroy(x);
	if( A ) lis_matrix_destroy(A);
	if( B ) lis_matrix_destroy(B);
	lis_finalize();
	return 1;
}
