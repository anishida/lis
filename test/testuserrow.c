#include "lis.h"
#include <stdio.h>

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

static void scalar_copy(
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

static int scalar_equal(
        const LIS_SCALAR *a,
        const LIS_SCALAR *b)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        return
                (*a)[0]==(*b)[0] &&
                (*a)[1]==(*b)[1];
#else
        return *a==*b;
#endif
}

typedef struct
{
        LIS_INT n;
        LIS_INT gn;
        LIS_INT is;
        LIS_SCALAR diag;
        LIS_SCALAR offdiag;
} USER_ROW_CTX;

static LIS_INT user_matvec(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        USER_ROW_CTX *ctx;
        LIS_INT i;

        ctx = (USER_ROW_CTX *)user_data;

        for(i=0;i<ctx->n;i++)
        {
                scalar_copy(&y[i],&x[i]);
        }

        return LIS_SUCCESS;
}

static LIS_INT user_get_row(
        void *user_data,
        LIS_INT local_row,
        LIS_INT *nnz,
        LIS_INT *cols,
        LIS_SCALAR *vals)
{
        USER_ROW_CTX *ctx;
        LIS_INT global_row;
        LIS_INT required;
        LIS_INT capacity;
        LIS_INT k;

        ctx = (USER_ROW_CTX *)user_data;

        if( nnz==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        if( local_row<0 || local_row>=ctx->n )
        {
                return LIS_ERR_ILL_ARG;
        }

        global_row = ctx->is + local_row;

        required = 1;

        if( global_row>0 )
        {
                required++;
        }

        if( global_row+1<ctx->gn )
        {
                required++;
        }

        if( cols==NULL && vals==NULL )
        {
                *nnz = required;
                return LIS_SUCCESS;
        }

        if( cols==NULL || vals==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        capacity = *nnz;

        if( capacity<required )
        {
                *nnz = required;
                return LIS_ERR_ILL_ARG;
        }

        k = 0;

        if( global_row>0 )
        {
                cols[k] = global_row-1;
                scalar_copy(&vals[k],&ctx->offdiag);
                k++;
        }

        cols[k] = global_row;
        scalar_copy(&vals[k],&ctx->diag);
        k++;

        if( global_row+1<ctx->gn )
        {
                cols[k] = global_row+1;
                scalar_copy(&vals[k],&ctx->offdiag);
                k++;
        }

        *nnz = k;

        return LIS_SUCCESS;
}

static LIS_INT user_get_bad_row(
        void *user_data,
        LIS_INT local_row,
        LIS_INT *nnz,
        LIS_INT *cols,
        LIS_SCALAR *vals)
{
        USER_ROW_CTX *ctx;

        ctx = (USER_ROW_CTX *)user_data;

        if( local_row<0 || local_row>=ctx->n )
        {
                return LIS_ERR_ILL_ARG;
        }

        if( cols==NULL && vals==NULL )
        {
                *nnz = 1;
                return LIS_SUCCESS;
        }

        if( *nnz<1 )
        {
                return LIS_ERR_ILL_ARG;
        }

        cols[0] = ctx->gn;
        scalar_set_real(&vals[0],(LIS_REAL)1.0);
        *nnz = 1;

        return LIS_SUCCESS;
}

static int check_row(
        LIS_MATRIX A,
        USER_ROW_CTX *ctx,
        LIS_INT local_row)
{
        LIS_INT err;
        LIS_INT nnz;
        LIS_INT expected;
        LIS_INT cols[3];
        LIS_SCALAR vals[3];
        LIS_INT global_row;
        LIS_INT k;

        global_row = ctx->is + local_row;

        expected = 1;
        if( global_row>0 ) expected++;
        if( global_row+1<ctx->gn ) expected++;

        nnz = -1;

        err = lis_matrix_get_user_row(
                A,
                local_row,
                &nnz,
                NULL,
                NULL);

        if( err!=LIS_SUCCESS || nnz!=expected )
        {
                fprintf(stderr,
                        "row query failed: row=%d err=%d nnz=%d expected=%d\n",
                        (int)local_row,
                        (int)err,
                        (int)nnz,
                        (int)expected);
                return 1;
        }

        nnz = 3;

        err = lis_matrix_get_user_row(
                A,
                local_row,
                &nnz,
                cols,
                vals);

        if( err!=LIS_SUCCESS || nnz!=expected )
        {
                fprintf(stderr,"row fill failed\n");
                return 1;
        }

        k = 0;

        if( global_row>0 )
        {
                if( cols[k]!=global_row-1 ||
                    !scalar_equal(&vals[k],&ctx->offdiag) )
                {
                        fprintf(stderr,"wrong lower entry\n");
                        return 1;
                }
                k++;
        }

        if( cols[k]!=global_row ||
            !scalar_equal(&vals[k],&ctx->diag) )
        {
                fprintf(stderr,"wrong diagonal entry\n");
                return 1;
        }

        k++;

        if( global_row+1<ctx->gn )
        {
                if( cols[k]!=global_row+1 ||
                    !scalar_equal(&vals[k],&ctx->offdiag) )
                {
                        fprintf(stderr,"wrong upper entry\n");
                        return 1;
                }
        }

        return 0;
}

int main(int argc, char *argv[])
{
        LIS_MATRIX A;
        USER_ROW_CTX ctx;

        LIS_INT err;
        LIS_INT i;
        LIS_INT local_n;
        LIS_INT global_n;
        LIS_INT is;
        LIS_INT ie;
        LIS_INT nnz;
        LIS_INT cols[3];
        LIS_SCALAR vals[3];

        LIS_INT requested_global_n;
        int found_interior;

        requested_global_n = 8;

        lis_initialize(&argc,&argv);

#ifdef USE_MPI
        {
                int nprocs;

                MPI_Comm_size(LIS_COMM_WORLD,&nprocs);

                if( requested_global_n < (LIS_INT)(2*nprocs) )
                {
                        requested_global_n = (LIS_INT)(2*nprocs);
                }
        }
#endif

        A = NULL;

        err = lis_matrix_create(LIS_COMM_WORLD,&A);
        if( err ) return 1;

        err = lis_matrix_set_size(A,0,requested_global_n);
        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        lis_matrix_get_size(A,&local_n,&global_n);
        lis_matrix_get_range(A,&is,&ie);

        (void)ie;

        ctx.n       = local_n;
        ctx.gn      = global_n;
        ctx.is      = is;
        scalar_set_real(&ctx.diag,(LIS_REAL)2.0);
        scalar_set_real(&ctx.offdiag,(LIS_REAL)-1.0);

        err = lis_matrix_set_user(
                A,
                &ctx,
                user_matvec,
                NULL);

        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        nnz = 0;

        err = lis_matrix_get_user_row(
                A,0,&nnz,NULL,NULL);

        if( err!=LIS_ERR_NOT_IMPLEMENTED )
        {
                fprintf(stderr,
                        "missing row callback was not rejected\n");
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        err = lis_matrix_set_user_row(A,user_get_row);

        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        for(i=0;i<local_n;i++)
        {
                if( check_row(A,&ctx,i) )
                {
                        lis_matrix_destroy(A);
                        lis_finalize();
                        return 1;
                }
        }

        nnz = 3;

        err = lis_matrix_get_user_row(
                A,0,&nnz,cols,NULL);

        if( err!=LIS_ERR_ILL_ARG )
        {
                fprintf(stderr,
                        "mixed cols/vals pointers were not rejected\n");
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        found_interior = 0;

        for(i=0;i<local_n;i++)
        {
                LIS_INT global_row;

                global_row = is+i;

                if( global_row>0 &&
                    global_row+1<global_n )
                {
                        nnz = 2;

                        err = lis_matrix_get_user_row(
                                A,i,&nnz,cols,vals);

                        if( err!=LIS_ERR_ILL_ARG )
                        {
                                fprintf(stderr,
                                        "small row buffer was not rejected\n");
                                lis_matrix_destroy(A);
                                lis_finalize();
                                return 1;
                        }

                        found_interior = 1;
                        break;
                }
        }

        if( !found_interior )
        {
                fprintf(stderr,
                        "no interior row available for capacity test\n");
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        nnz = 0;

        err = lis_matrix_get_user_row(
                A,local_n,&nnz,NULL,NULL);

        if( err!=LIS_ERR_ILL_ARG )
        {
                fprintf(stderr,
                        "out-of-range local row was not rejected\n");
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        err = lis_matrix_set_user_row(
                A,user_get_bad_row);

        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        nnz = 1;

        err = lis_matrix_get_user_row(
                A,0,&nnz,cols,vals);

        if( err!=LIS_ERR_ILL_ARG )
        {
                fprintf(stderr,
                        "invalid global column was not rejected\n");
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        err = lis_matrix_set_user_row(A,NULL);

        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        nnz = 0;

        err = lis_matrix_get_user_row(
                A,0,&nnz,NULL,NULL);

        if( err!=LIS_ERR_NOT_IMPLEMENTED )
        {
                fprintf(stderr,
                        "cleared row callback was not rejected\n");
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        lis_matrix_destroy(A);
        lis_finalize();

        return 0;
}
