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
} USER_DIAG;

static LIS_INT user_matvec(void *ctx,
                           const LIS_SCALAR *x,
                           LIS_SCALAR *y)
{
        USER_DIAG *A = (USER_DIAG *)ctx;
        LIS_INT i;

        for(i=0;i<A->n;i++)
                y[i] = A->diag*x[i];

        return LIS_SUCCESS;
}

static LIS_INT user_matvech(void *ctx,
                            const LIS_SCALAR *x,
                            LIS_SCALAR *y)
{
#ifdef _COMPLEX
        USER_DIAG *A = (USER_DIAG *)ctx;
        LIS_INT i;

        for(i=0;i<A->n;i++)
                y[i] = conj(A->diag)*x[i];

        return LIS_SUCCESS;
#else
        return user_matvec(ctx,x,y);
#endif
}

static LIS_INT user_get_diagonal(void *ctx, LIS_SCALAR *d)
{
        USER_DIAG *A = (USER_DIAG *)ctx;
        LIS_INT i;

        for(i=0;i<A->n;i++)
                d[i] = A->diag;

        return LIS_SUCCESS;
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

static void print_mismatch(const char *name,
                           LIS_INT gi,
                           LIS_SCALAR got,
                           LIS_SCALAR expected)
{
#ifdef _COMPLEX
        fprintf(stderr,
                "%s mismatch at global index %d: "
                "got=(%.15e,%.15e) expected=(%.15e,%.15e)\n",
                name,(int)gi,
                (double)creal(got),(double)cimag(got),
                (double)creal(expected),(double)cimag(expected));
#else
        fprintf(stderr,
                "%s mismatch at global index %d: "
                "got=%.15e expected=%.15e\n",
                name,(int)gi,(double)got,(double)expected);
#endif
}

static int check_affine(const char *name,
                        LIS_VECTOR d,
                        LIS_SCALAR slope,
                        LIS_SCALAR offset)
{
        LIS_INT i,gi;
        LIS_SCALAR expected;

        for(i=0;i<d->n;i++)
        {
                gi = d->is+i;
                expected = slope*(LIS_SCALAR)(gi+1) + offset;

                if( !scalar_close(d->value[i],expected) )
                {
                        print_mismatch(name,gi,d->value[i],expected);
                        return 1;
                }
        }

        return 0;
}

static int expect_error(const char *name,
                        LIS_INT got,
                        LIS_INT expected)
{
        if( got!=expected )
        {
                fprintf(stderr,
                        "%s: expected error %d, got %d\n",
                        name,(int)expected,(int)got);
                return 1;
        }

        return 0;
}

int main(int argc, char *argv[])
{
        LIS_MATRIX A = NULL;
        LIS_MATRIX U = NULL;
        LIS_MATRIX V = NULL;
        LIS_MATRIX C1 = NULL;
        LIS_MATRIX C2 = NULL;
        LIS_MATRIX C3 = NULL;
        LIS_MATRIX N = NULL;
#ifdef _COMPLEX
        LIS_MATRIX X = NULL;
#endif
        LIS_VECTOR d = NULL;
        USER_DIAG uctx,vctx;
        LIS_INT i,is,ie,err;
        LIS_INT global_n = 4;
        LIS_INT local_fail = 0;
        LIS_INT global_fail = 0;

        lis_initialize(&argc,&argv);

        /*
         * Baseline CSR matrix:
         *
         * diag(A)[i] = i + 1
         */
        err = lis_matrix_create(LIS_COMM_WORLD,&A);
        if( err ) goto api_fail;

        if( A->nprocs>global_n )
                global_n = A->nprocs;

        err = lis_matrix_set_size(A,0,global_n);
        if( err ) goto api_fail;

        err = lis_matrix_get_range(A,&is,&ie);
        if( err ) goto api_fail;

        for(i=is;i<ie;i++)
        {
                err = lis_matrix_set_value(
                        LIS_INS_VALUE,i,i,(LIS_SCALAR)(i+1),A);
                if( err ) goto api_fail;
        }

        err = lis_matrix_assemble(A);
        if( err ) goto api_fail;

        err = lis_vector_duplicate(A,&d);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(A,d);
        if( err ) goto api_fail;

        local_fail |= check_affine(
                "CSR diagonal",d,
                (LIS_SCALAR)1.0,
                (LIS_SCALAR)0.0);

        /*
         * Setting a USER diagonal callback on a non-USER matrix
         * must be rejected.
         */
        err = lis_matrix_set_user_diagonal(A,user_get_diagonal);
        local_fail |= expect_error(
                "diagonal callback on CSR",
                err,LIS_ERR_ILL_ARG);

        /*
         * U = 5*I as LIS_MATRIX_USER.
         */
        err = lis_matrix_create(LIS_COMM_WORLD,&U);
        if( err ) goto api_fail;

        err = lis_matrix_set_size(U,0,global_n);
        if( err ) goto api_fail;

        err = lis_matrix_get_range(U,&is,&ie);
        if( err ) goto api_fail;

        uctx.n = ie-is;
        uctx.diag = (LIS_SCALAR)5.0;

        err = lis_matrix_set_user(
                U,&uctx,user_matvec,user_matvech);
        if( err ) goto api_fail;

        err = lis_matrix_set_user_diagonal(
                U,user_get_diagonal);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(U,d);
        if( err ) goto api_fail;

        local_fail |= check_affine(
                "USER diagonal",d,
                (LIS_SCALAR)0.0,
                (LIS_SCALAR)5.0);

        /*
         * V = -2*I.
         *
         * First verify that USER without a diagonal callback
         * remains explicitly unsupported.
         */
        err = lis_matrix_create(LIS_COMM_WORLD,&V);
        if( err ) goto api_fail;

        err = lis_matrix_set_size(V,0,global_n);
        if( err ) goto api_fail;

        err = lis_matrix_get_range(V,&is,&ie);
        if( err ) goto api_fail;

        vctx.n = ie-is;
        vctx.diag = (LIS_SCALAR)-2.0;

        err = lis_matrix_set_user(
                V,&vctx,user_matvec,user_matvech);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(V,d);
        local_fail |= expect_error(
                "USER without diagonal callback",
                err,LIS_ERR_NOT_IMPLEMENTED);

        err = lis_matrix_set_user_diagonal(
                V,user_get_diagonal);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(V,d);
        if( err ) goto api_fail;

        local_fail |= check_affine(
                "second USER diagonal",d,
                (LIS_SCALAR)0.0,
                (LIS_SCALAR)-2.0);

        /*
         * USER + CSR:
         *
         * C1 = 2*U - A
         *
         * diag(C1) = 10 - (i+1)
         */
        err = lis_matrix_create_operator(
                (LIS_SCALAR)2.0,U,
                (LIS_SCALAR)-1.0,A,&C1);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(C1,d);
        if( err ) goto api_fail;

        local_fail |= check_affine(
                "USER + CSR diagonal",d,
                (LIS_SCALAR)-1.0,
                (LIS_SCALAR)10.0);

        /*
         * CSR + USER:
         *
         * C2 = 3*A + 0.5*U
         *
         * diag(C2) = 3*(i+1) + 2.5
         */
        err = lis_matrix_create_operator(
                (LIS_SCALAR)3.0,A,
                (LIS_SCALAR)0.5,U,&C2);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(C2,d);
        if( err ) goto api_fail;

        local_fail |= check_affine(
                "CSR + USER diagonal",d,
                (LIS_SCALAR)3.0,
                (LIS_SCALAR)2.5);

        /*
         * USER + USER:
         *
         * C3 = 3*U + 0.5*V = 14*I
         */
        err = lis_matrix_create_operator(
                (LIS_SCALAR)3.0,U,
                (LIS_SCALAR)0.5,V,&C3);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(C3,d);
        if( err ) goto api_fail;

        local_fail |= check_affine(
                "USER + USER diagonal",d,
                (LIS_SCALAR)0.0,
                (LIS_SCALAR)14.0);

        /*
         * Nested operator:
         *
         * N = 0.5*C1 + 2*C2
         *
         *   = 0.5*(10-A) + 2*(3*A+2.5*I)
         *   = 5.5*A + 10*I
         */
        err = lis_matrix_create_operator(
                (LIS_SCALAR)0.5,C1,
                (LIS_SCALAR)2.0,C2,&N);
        if( err ) goto api_fail;

        err = lis_matrix_get_diagonal(N,d);
        if( err ) goto api_fail;

        local_fail |= check_affine(
                "nested operator diagonal",d,
                (LIS_SCALAR)5.5,
                (LIS_SCALAR)10.0);

#ifdef _COMPLEX
        /*
         * get_diagonal() uses the coefficients of C itself.
         * Unlike matvech(), alpha and beta must NOT be conjugated.
         */
        {
                LIS_SCALAR alpha =
                        (LIS_SCALAR)(2.0 + 1.0*I);
                LIS_SCALAR beta =
                        (LIS_SCALAR)(-1.0 + 2.0*I);

                err = lis_matrix_create_operator(
                        alpha,U,beta,A,&X);
                if( err ) goto api_fail;

                err = lis_matrix_get_diagonal(X,d);
                if( err ) goto api_fail;

                local_fail |= check_affine(
                        "complex operator diagonal",d,
                        beta,
                        (LIS_SCALAR)5.0*alpha);
        }
#endif

#ifdef USE_MPI
        {
                int lf = (int)local_fail;
                int gf = 0;

                MPI_Allreduce(
                        &lf,&gf,1,MPI_INT,MPI_MAX,LIS_COMM_WORLD);

                global_fail = (LIS_INT)gf;
        }
#else
        global_fail = local_fail;
#endif

        if( A->my_rank==0 )
        {
                printf("LIS_MATRIX_DIAGONAL %s\n",
                       global_fail ? "FAILED" : "PASSED");
        }

#ifdef _COMPLEX
        if( X ) lis_matrix_destroy(X);
#endif
        if( N ) lis_matrix_destroy(N);
        if( C3 ) lis_matrix_destroy(C3);
        if( C2 ) lis_matrix_destroy(C2);
        if( C1 ) lis_matrix_destroy(C1);
        if( V ) lis_matrix_destroy(V);
        if( U ) lis_matrix_destroy(U);
        if( d ) lis_vector_destroy(d);
        if( A ) lis_matrix_destroy(A);

        lis_finalize();
        return global_fail ? 1 : 0;

api_fail:
        fprintf(stderr,
                "LIS_MATRIX_DIAGONAL API failure: %d\n",(int)err);

#ifdef _COMPLEX
        if( X ) lis_matrix_destroy(X);
#endif
        if( N ) lis_matrix_destroy(N);
        if( C3 ) lis_matrix_destroy(C3);
        if( C2 ) lis_matrix_destroy(C2);
        if( C1 ) lis_matrix_destroy(C1);
        if( V ) lis_matrix_destroy(V);
        if( U ) lis_matrix_destroy(U);
        if( d ) lis_vector_destroy(d);
        if( A ) lis_matrix_destroy(A);

        lis_finalize();
        return 1;
}
