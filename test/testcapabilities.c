#include "lis.h"
#include <stdio.h>


typedef struct
{
        LIS_INT n;
        LIS_INT gn;
        LIS_INT is;
} USER_CAP_CTX;


static void
scalar_set_real(
        LIS_SCALAR *x,
        LIS_REAL value)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        (*x)[0] = value;
        (*x)[1] = (LIS_REAL)0.0;
#else
        *x = (LIS_SCALAR)value;
#endif
}


static void
scalar_copy(
        LIS_SCALAR *dst,
        const LIS_SCALAR *src)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        (*dst)[0] = (*src)[0];
        (*dst)[1] = (*src)[1];
#else
        *dst = *src;
#endif
}


static LIS_INT
user_matvec(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        USER_CAP_CTX *ctx;
        LIS_INT i;

        ctx = (USER_CAP_CTX *)user_data;

        for(i=0;i<ctx->n;i++)
        {
                scalar_copy(&y[i],&x[i]);
        }

        return LIS_SUCCESS;
}


static LIS_INT
user_matvech(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        return user_matvec(
                user_data,
                x,
                y);
}


static LIS_INT
user_get_diagonal(
        void *user_data,
        LIS_SCALAR *d)
{
        USER_CAP_CTX *ctx;
        LIS_INT i;

        ctx = (USER_CAP_CTX *)user_data;

        for(i=0;i<ctx->n;i++)
        {
                scalar_set_real(
                        &d[i],
                        (LIS_REAL)1.0);
        }

        return LIS_SUCCESS;
}


static LIS_INT
user_get_row(
        void *user_data,
        LIS_INT local_row,
        LIS_INT *nnz,
        LIS_INT *cols,
        LIS_SCALAR *vals)
{
        USER_CAP_CTX *ctx;

        ctx = (USER_CAP_CTX *)user_data;

        if( nnz==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        if(
                local_row<0
                ||
                local_row>=ctx->n)
        {
                return LIS_ERR_ILL_ARG;
        }

        if( cols==NULL && vals==NULL )
        {
                *nnz = 1;
                return LIS_SUCCESS;
        }

        if( cols==NULL || vals==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        if( *nnz<1 )
        {
                *nnz = 1;
                return LIS_ERR_ILL_ARG;
        }

        cols[0] =
                ctx->is
                +
                local_row;

        scalar_set_real(
                &vals[0],
                (LIS_REAL)1.0);

        *nnz = 1;

        return LIS_SUCCESS;
}


static int
expect_capabilities(
        const char *name,
        LIS_MATRIX A,
        LIS_INT expected)
{
        LIS_INT err;
        LIS_INT got;

        got = -1;

        err = lis_matrix_get_capabilities(
                A,
                &got);

        if(
                err!=LIS_SUCCESS
                ||
                got!=expected)
        {
                fprintf(
                        stderr,
                        "%s: err=%d got=0x%x expected=0x%x\n",
                        name,
                        (int)err,
                        (unsigned int)got,
                        (unsigned int)expected);

                return 1;
        }

        return 0;
}


int
main(
        int argc,
        char *argv[])
{
        LIS_MATRIX A;
        LIS_MATRIX U;
        LIS_MATRIX V;
        LIS_MATRIX C1;
        LIS_MATRIX C2;
        LIS_MATRIX N;

        USER_CAP_CTX uctx;
        USER_CAP_CTX vctx;

        LIS_INT err;
        LIS_INT is;
        LIS_INT ie;
        LIS_INT global_n;
        LIS_INT caps;
        LIS_INT i;

        LIS_INT cap_matvec;
        LIS_INT cap_matvech;
        LIS_INT cap_diagonal;
        LIS_INT cap_row;
        LIS_INT cap_storage;

        LIS_SCALAR one;

        int fail;

#ifdef USE_MPI
        int nprocs;
#endif

        A = NULL;
        U = NULL;
        V = NULL;
        C1 = NULL;
        C2 = NULL;
        N = NULL;

        fail = 0;

        cap_matvec =
                LIS_MATRIX_CAP_MATVEC;

        cap_matvech =
                LIS_MATRIX_CAP_MATVECH;

        cap_diagonal =
                LIS_MATRIX_CAP_DIAGONAL;

        cap_row =
                LIS_MATRIX_CAP_ROW_ACCESS;

        cap_storage =
                cap_matvec
                |
                cap_matvech
                |
                cap_diagonal;

        scalar_set_real(
                &one,
                (LIS_REAL)1.0);

        global_n = 8;

        lis_initialize(
                &argc,
                &argv);

#ifdef USE_MPI
        MPI_Comm_size(
                LIS_COMM_WORLD,
                &nprocs);

        if(
                global_n
                <
                (LIS_INT)(2*nprocs))
        {
                global_n =
                        (LIS_INT)(2*nprocs);
        }
#endif


        /*
         * Unsized matrices do not yet have a valid operation layout.
         */
        err = lis_matrix_create(
                LIS_COMM_WORLD,
                &N);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        caps = -1;

        err = lis_matrix_get_capabilities(
                N,
                &caps);

        if( err!=LIS_ERR_ILL_ARG )
        {
                fprintf(
                        stderr,
                        "unsized matrix was not rejected: err=%d\n",
                        (int)err);

                fail = 1;
        }

        lis_matrix_destroy(N);
        N = NULL;


        /*
         * Baseline explicit CSR matrix.
         */
        err = lis_matrix_create(
                LIS_COMM_WORLD,
                &A);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        err = lis_matrix_set_size(
                A,
                0,
                global_n);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        err = lis_matrix_get_range(
                A,
                &is,
                &ie);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        for(i=is;i<ie;i++)
        {
                err = lis_matrix_set_value(
                        LIS_INS_VALUE,
                        i,
                        i,
                        one,
                        A);

                if( err )
                {
                        fail = 1;
                        goto cleanup;
                }
        }

        err = lis_matrix_assemble(A);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "CSR",
                A,
                cap_storage);


        /*
         * NULL output pointers must be rejected.
         */
        err = lis_matrix_get_capabilities(
                A,
                NULL);

        if( err!=LIS_ERR_ILL_ARG )
        {
                fprintf(
                        stderr,
                        "NULL capability output was not rejected: err=%d\n",
                        (int)err);

                fail = 1;
        }


        /*
         * USER matrix with matvec + matvech.
         */
        err = lis_matrix_create(
                LIS_COMM_WORLD,
                &U);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        err = lis_matrix_set_size(
                U,
                0,
                global_n);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        err = lis_matrix_get_range(
                U,
                &is,
                &ie);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        uctx.n = ie-is;
        uctx.gn = global_n;
        uctx.is = is;

        err = lis_matrix_set_user(
                U,
                &uctx,
                user_matvec,
                user_matvech);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "USER matvec+matvech",
                U,
                cap_matvec
                |
                cap_matvech);


        /*
         * USER capabilities must follow callbacks dynamically.
         */
        err = lis_matrix_set_user_diagonal(
                U,
                user_get_diagonal);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "USER + diagonal",
                U,
                cap_matvec
                |
                cap_matvech
                |
                cap_diagonal);

        err = lis_matrix_set_user_row(
                U,
                user_get_row);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "USER + diagonal + row",
                U,
                cap_matvec
                |
                cap_matvech
                |
                cap_diagonal
                |
                cap_row);

        err = lis_matrix_set_user_row(
                U,
                NULL);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "USER row cleared",
                U,
                cap_matvec
                |
                cap_matvech
                |
                cap_diagonal);

        err = lis_matrix_set_user_row(
                U,
                user_get_row);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        err = lis_matrix_set_user_diagonal(
                U,
                NULL);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "USER diagonal cleared",
                U,
                cap_matvec
                |
                cap_matvech
                |
                cap_row);

        err = lis_matrix_set_user_diagonal(
                U,
                user_get_diagonal);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }


        /*
         * A + U:
         *
         * both operands provide matvec, matvech and diagonal.
         * U also provides row access, but LIS_MATRIX_OPERATOR does not.
         */
        err = lis_matrix_create_operator(
                one,
                A,
                one,
                U,
                &C1);

        if( err )
        {
                fprintf(
                        stderr,
                        "operator creation failed: %d\n",
                        (int)err);

                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "CSR + full USER operator",
                C1,
                cap_storage);


        /*
         * V intentionally provides only matvec.
         */
        err = lis_matrix_create(
                LIS_COMM_WORLD,
                &V);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        err = lis_matrix_set_size(
                V,
                0,
                global_n);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        err = lis_matrix_get_range(
                V,
                &is,
                &ie);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        vctx.n = ie-is;
        vctx.gn = global_n;
        vctx.is = is;

        err = lis_matrix_set_user(
                V,
                &vctx,
                user_matvec,
                NULL);

        if( err )
        {
                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "USER matvec only",
                V,
                cap_matvec);


        /*
         * Intersection semantics:
         *
         * CSR provides three generic operations, but V only provides
         * matvec. Therefore A + V exposes only matvec.
         */
        err = lis_matrix_create_operator(
                one,
                A,
                one,
                V,
                &C2);

        if( err )
        {
                fprintf(
                        stderr,
                        "restricted operator creation failed: %d\n",
                        (int)err);

                fail = 1;
                goto cleanup;
        }

        fail |= expect_capabilities(
                "CSR + matvec-only USER operator",
                C2,
                cap_matvec);


cleanup:

        if( C2 )
                lis_matrix_destroy(C2);

        if( C1 )
                lis_matrix_destroy(C1);

        if( V )
                lis_matrix_destroy(V);

        if( U )
                lis_matrix_destroy(U);

        if( A )
                lis_matrix_destroy(A);

        if( N )
                lis_matrix_destroy(N);

        lis_finalize();

        if( fail )
        {
                fprintf(
                        stderr,
                        "LIS_MATRIX_CAPABILITIES FAILED\n");

                return 1;
        }

        printf(
                "LIS_MATRIX_CAPABILITIES PASSED\n");

        return 0;
}
