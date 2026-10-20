#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include <stdio.h>

#define BLOCK_NVECS 3
#define BLOCK_FACTOR ((LIS_REAL)2.0)

typedef struct
{
        LIS_INT n;

        LIS_INT scalar_matvec_calls;
        LIS_INT scalar_matvech_calls;

        LIS_INT block_matvec_calls;
        LIS_INT block_matvech_calls;

        LIS_INT fail_block_matvec;
        LIS_INT fail_block_matvech;
} USER_BLOCK_CTX;


static void scalar_set_real(
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


static void scalar_scale_real(
        LIS_SCALAR *dst,
        const LIS_SCALAR *src,
        LIS_REAL factor)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        (*dst)[0] = factor * (*src)[0];
        (*dst)[1] = factor * (*src)[1];
#else
        *dst = (LIS_SCALAR)(factor * (*src));
#endif
}


static int scalar_equal(
        const LIS_SCALAR *a,
        const LIS_SCALAR *b)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        return
                (*a)[0] == (*b)[0] &&
                (*a)[1] == (*b)[1];
#else
        return *a == *b;
#endif
}


static void reset_counts(
        USER_BLOCK_CTX *ctx)
{
        ctx->scalar_matvec_calls = 0;
        ctx->scalar_matvech_calls = 0;

        ctx->block_matvec_calls = 0;
        ctx->block_matvech_calls = 0;
}


static LIS_INT user_matvec(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        USER_BLOCK_CTX *ctx;
        LIS_INT i;

        ctx = (USER_BLOCK_CTX *)user_data;

        ctx->scalar_matvec_calls++;

        for(i=0;i<ctx->n;i++)
        {
                scalar_scale_real(
                        &y[i],
                        &x[i],
                        BLOCK_FACTOR);
        }

        return LIS_SUCCESS;
}


static LIS_INT user_matvech(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        USER_BLOCK_CTX *ctx;
        LIS_INT i;

        ctx = (USER_BLOCK_CTX *)user_data;

        ctx->scalar_matvech_calls++;

        for(i=0;i<ctx->n;i++)
        {
                scalar_scale_real(
                        &y[i],
                        &x[i],
                        BLOCK_FACTOR);
        }

        return LIS_SUCCESS;
}


static LIS_INT user_matvec_block(
        void *user_data,
        LIS_INT nvec,
        LIS_VECTOR X[],
        LIS_VECTOR Y[])
{
        USER_BLOCK_CTX *ctx;
        LIS_INT i,j;

        ctx = (USER_BLOCK_CTX *)user_data;

        ctx->block_matvec_calls++;

        if( ctx->fail_block_matvec )
        {
                return LIS_ERR_ILL_ARG;
        }

        for(j=0;j<nvec;j++)
        {
                for(i=0;i<ctx->n;i++)
                {
                        scalar_scale_real(
                                &Y[j]->value[i],
                                &X[j]->value[i],
                                BLOCK_FACTOR);
                }
        }

        return LIS_SUCCESS;
}


static LIS_INT user_matvech_block(
        void *user_data,
        LIS_INT nvec,
        LIS_VECTOR X[],
        LIS_VECTOR Y[])
{
        USER_BLOCK_CTX *ctx;
        LIS_INT i,j;

        ctx = (USER_BLOCK_CTX *)user_data;

        ctx->block_matvech_calls++;

        if( ctx->fail_block_matvech )
        {
                return LIS_ERR_ILL_ARG;
        }

        for(j=0;j<nvec;j++)
        {
                for(i=0;i<ctx->n;i++)
                {
                        scalar_scale_real(
                                &Y[j]->value[i],
                                &X[j]->value[i],
                                BLOCK_FACTOR);
                }
        }

        return LIS_SUCCESS;
}


static void fill_vectors(
        LIS_VECTOR X[],
        LIS_VECTOR Y[],
        LIS_INT nvec,
        LIS_INT n)
{
        LIS_INT i,j;

        for(j=0;j<nvec;j++)
        {
                for(i=0;i<n;i++)
                {
                        scalar_set_real(
                                &X[j]->value[i],
                                (LIS_REAL)(
                                        10*(j+1) +
                                        (i+1)));

                        scalar_set_real(
                                &Y[j]->value[i],
                                (LIS_REAL)-1.0);
                }
        }
}


static int check_values(
        LIS_VECTOR X[],
        LIS_VECTOR Y[],
        LIS_INT nvec,
        LIS_INT n)
{
        LIS_INT i,j;
        LIS_SCALAR expected;

        for(j=0;j<nvec;j++)
        {
                for(i=0;i<n;i++)
                {
                        scalar_scale_real(
                                &expected,
                                &X[j]->value[i],
                                BLOCK_FACTOR);

                        if( !scalar_equal(
                                &expected,
                                &Y[j]->value[i]) )
                        {
                                fprintf(
                                        stderr,
                                        "wrong block result: "
                                        "vec=%d row=%d\n",
                                        (int)j,
                                        (int)i);

                                return 1;
                        }
                }
        }

        return 0;
}


int main(int argc, char *argv[])
{
        LIS_MATRIX A;
        LIS_VECTOR X[BLOCK_NVECS];
        LIS_VECTOR Y[BLOCK_NVECS];

        USER_BLOCK_CTX ctx;

        LIS_INT err;
        LIS_INT i;
        LIS_INT is,ie;
        LIS_INT global_n;

        LIS_INT local_fail;
        LIS_INT global_fail;

        A = NULL;

        for(i=0;i<BLOCK_NVECS;i++)
        {
                X[i] = NULL;
                Y[i] = NULL;
        }

        ctx.n = 0;

        ctx.scalar_matvec_calls = 0;
        ctx.scalar_matvech_calls = 0;

        ctx.block_matvec_calls = 0;
        ctx.block_matvech_calls = 0;

        ctx.fail_block_matvec = 0;
        ctx.fail_block_matvech = 0;

        local_fail = 0;
        global_fail = 0;

        lis_initialize(&argc,&argv);


        /*
         * Create a simple USER operator:
         *
         *     A = 2 I
         *
         * Scalar and block callbacks implement the same operator.
         */
        err = lis_matrix_create(
                LIS_COMM_WORLD,
                &A);

        if( err )
        {
                fprintf(stderr,
                        "matrix create failed: %d\n",
                        (int)err);

                local_fail = 1;
                goto cleanup;
        }

        global_n = 4;

        if( A->nprocs > global_n )
        {
                global_n = A->nprocs;
        }

        err = lis_matrix_set_size(
                A,
                0,
                global_n);

        if( err )
        {
                fprintf(stderr,
                        "matrix size failed: %d\n",
                        (int)err);

                local_fail = 1;
                goto cleanup;
        }

        err = lis_matrix_get_range(
                A,
                &is,
                &ie);

        if( err )
        {
                fprintf(stderr,
                        "matrix range failed: %d\n",
                        (int)err);

                local_fail = 1;
                goto cleanup;
        }

        ctx.n = ie-is;

        err = lis_matrix_set_user(
                A,
                &ctx,
                user_matvec,
                user_matvech);

        if( err )
        {
                fprintf(stderr,
                        "set USER matrix failed: %d\n",
                        (int)err);

                local_fail = 1;
                goto cleanup;
        }


        for(i=0;i<BLOCK_NVECS;i++)
        {
                err = lis_vector_duplicate(
                        A,
                        &X[i]);

                if( err )
                {
                        fprintf(stderr,
                                "X duplicate failed: %d\n",
                                (int)err);

                        local_fail = 1;
                        goto cleanup;
                }

                err = lis_vector_duplicate(
                        A,
                        &Y[i]);

                if( err )
                {
                        fprintf(stderr,
                                "Y duplicate failed: %d\n",
                                (int)err);

                        local_fail = 1;
                        goto cleanup;
                }
        }


        /*
         * nvec == 0 must be a no-op.
         */
        err = lis_matvec_block(
                A,
                0,
                NULL,
                NULL);

        if( err!=LIS_SUCCESS )
        {
                fprintf(stderr,
                        "zero-length block failed: %d\n",
                        (int)err);

                local_fail = 1;
        }


        /*
         * No native block callback yet:
         * lis_matvec_block() must fall back to nvec scalar calls.
         */
        fill_vectors(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n);

        reset_counts(&ctx);

        err = lis_matvec_block(
                A,
                BLOCK_NVECS,
                X,
                Y);

        if( err!=LIS_SUCCESS )
        {
                fprintf(stderr,
                        "scalar matvec fallback failed: %d\n",
                        (int)err);

                local_fail = 1;
        }

        if( ctx.scalar_matvec_calls!=BLOCK_NVECS ||
            ctx.block_matvec_calls!=0 )
        {
                fprintf(
                        stderr,
                        "wrong matvec fallback dispatch: "
                        "scalar=%d block=%d\n",
                        (int)ctx.scalar_matvec_calls,
                        (int)ctx.block_matvec_calls);

                local_fail = 1;
        }

        if( check_values(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n) )
        {
                local_fail = 1;
        }


        /*
         * Same fallback contract for matvech.
         */
        fill_vectors(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n);

        reset_counts(&ctx);

        err = lis_matvech_block(
                A,
                BLOCK_NVECS,
                X,
                Y);

        if( err!=LIS_SUCCESS )
        {
                fprintf(stderr,
                        "scalar matvech fallback failed: %d\n",
                        (int)err);

                local_fail = 1;
        }

        if( ctx.scalar_matvech_calls!=BLOCK_NVECS ||
            ctx.block_matvech_calls!=0 )
        {
                fprintf(
                        stderr,
                        "wrong matvech fallback dispatch: "
                        "scalar=%d block=%d\n",
                        (int)ctx.scalar_matvech_calls,
                        (int)ctx.block_matvech_calls);

                local_fail = 1;
        }

        if( check_values(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n) )
        {
                local_fail = 1;
        }


        /*
         * Install native block callbacks.
         */
        err = lis_matrix_set_user_matvec_block(
                A,
                user_matvec_block,
                user_matvech_block);

        if( err )
        {
                fprintf(stderr,
                        "block callback setter failed: %d\n",
                        (int)err);

                local_fail = 1;
                goto cleanup;
        }


        /*
         * Native matvec must invoke one callback for the entire block.
         */
        fill_vectors(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n);

        reset_counts(&ctx);

        err = lis_matvec_block(
                A,
                BLOCK_NVECS,
                X,
                Y);

        if( err!=LIS_SUCCESS )
        {
                fprintf(stderr,
                        "native block matvec failed: %d\n",
                        (int)err);

                local_fail = 1;
        }

        if( ctx.block_matvec_calls!=1 ||
            ctx.scalar_matvec_calls!=0 )
        {
                fprintf(
                        stderr,
                        "wrong native matvec dispatch: "
                        "block=%d scalar=%d\n",
                        (int)ctx.block_matvec_calls,
                        (int)ctx.scalar_matvec_calls);

                local_fail = 1;
        }

        if( check_values(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n) )
        {
                local_fail = 1;
        }


        /*
         * Native matvech has the same one-call contract.
         */
        fill_vectors(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n);

        reset_counts(&ctx);

        err = lis_matvech_block(
                A,
                BLOCK_NVECS,
                X,
                Y);

        if( err!=LIS_SUCCESS )
        {
                fprintf(stderr,
                        "native block matvech failed: %d\n",
                        (int)err);

                local_fail = 1;
        }

        if( ctx.block_matvech_calls!=1 ||
            ctx.scalar_matvech_calls!=0 )
        {
                fprintf(
                        stderr,
                        "wrong native matvech dispatch: "
                        "block=%d scalar=%d\n",
                        (int)ctx.block_matvech_calls,
                        (int)ctx.scalar_matvech_calls);

                local_fail = 1;
        }

        if( check_values(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n) )
        {
                local_fail = 1;
        }


        /*
         * nvec == 1 still uses the native block callback.
         */
        fill_vectors(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n);

        reset_counts(&ctx);

        err = lis_matvec_block(
                A,
                1,
                X,
                Y);

        if( err!=LIS_SUCCESS )
        {
                fprintf(stderr,
                        "single-vector native block failed: %d\n",
                        (int)err);

                local_fail = 1;
        }

        if( ctx.block_matvec_calls!=1 ||
            ctx.scalar_matvec_calls!=0 )
        {
                fprintf(
                        stderr,
                        "nvec=1 did not use native block callback\n");

                local_fail = 1;
        }

        if( check_values(
                X,
                Y,
                1,
                ctx.n) )
        {
                local_fail = 1;
        }


        /*
         * Callback error must propagate unchanged.
         */
        reset_counts(&ctx);
        ctx.fail_block_matvec = 1;

        err = lis_matvec_block(
                A,
                BLOCK_NVECS,
                X,
                Y);

        ctx.fail_block_matvec = 0;

        if( err!=LIS_ERR_ILL_ARG ||
            ctx.block_matvec_calls!=1 )
        {
                fprintf(
                        stderr,
                        "matvec callback error was not propagated: "
                        "err=%d calls=%d\n",
                        (int)err,
                        (int)ctx.block_matvec_calls);

                local_fail = 1;
        }


        reset_counts(&ctx);
        ctx.fail_block_matvech = 1;

        err = lis_matvech_block(
                A,
                BLOCK_NVECS,
                X,
                Y);

        ctx.fail_block_matvech = 0;

        if( err!=LIS_ERR_ILL_ARG ||
            ctx.block_matvech_calls!=1 )
        {
                fprintf(
                        stderr,
                        "matvech callback error was not propagated: "
                        "err=%d calls=%d\n",
                        (int)err,
                        (int)ctx.block_matvech_calls);

                local_fail = 1;
        }


        /*
         * Clearing native callbacks restores scalar fallback.
         */
        err = lis_matrix_set_user_matvec_block(
                A,
                NULL,
                NULL);

        if( err )
        {
                fprintf(stderr,
                        "clearing block callbacks failed: %d\n",
                        (int)err);

                local_fail = 1;
                goto cleanup;
        }

        fill_vectors(
                X,
                Y,
                BLOCK_NVECS,
                ctx.n);

        reset_counts(&ctx);

        err = lis_matvec_block(
                A,
                BLOCK_NVECS,
                X,
                Y);

        if( err!=LIS_SUCCESS ||
            ctx.scalar_matvec_calls!=BLOCK_NVECS ||
            ctx.block_matvec_calls!=0 )
        {
                fprintf(
                        stderr,
                        "scalar fallback was not restored: "
                        "err=%d scalar=%d block=%d\n",
                        (int)err,
                        (int)ctx.scalar_matvec_calls,
                        (int)ctx.block_matvec_calls);

                local_fail = 1;
        }


        /*
         * Negative block length is invalid.
         */
        err = lis_matvec_block(
                A,
                -1,
                X,
                Y);

        if( err!=LIS_ERR_ILL_ARG )
        {
                fprintf(
                        stderr,
                        "negative nvec was not rejected: %d\n",
                        (int)err);

                local_fail = 1;
        }


cleanup:

#ifdef USE_MPI
        {
                int lf;
                int gf;

                lf = (int)local_fail;
                gf = 0;

                MPI_Allreduce(
                        &lf,
                        &gf,
                        1,
                        MPI_INT,
                        MPI_MAX,
                        LIS_COMM_WORLD);

                global_fail = (LIS_INT)gf;
        }
#else
        global_fail = local_fail;
#endif

        if( A==NULL || A->my_rank==0 )
        {
                printf(
                        "LIS_MATVEC_BLOCK %s\n",
                        global_fail ? "FAILED" : "PASSED");
        }

        for(i=0;i<BLOCK_NVECS;i++)
        {
                if( Y[i] ) lis_vector_destroy(Y[i]);
                if( X[i] ) lis_vector_destroy(X[i]);
        }

        if( A ) lis_matrix_destroy(A);

        lis_finalize();

        return global_fail ? 1 : 0;
}
