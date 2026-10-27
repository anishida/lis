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
        LIS_INT saw_parallel;
} USER_ROW_CTX;

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

#ifdef _OPENMP
        if( omp_in_parallel() && omp_get_num_threads()>1 )
        {
#pragma omp atomic write
                ctx->saw_parallel = 1;
        }
#endif

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

static LIS_INT user_get_duplicate_row(
        void *user_data,
        LIS_INT local_row,
        LIS_INT *nnz,
        LIS_INT *cols,
        LIS_SCALAR *vals)
{
        LIS_INT required;
        (void)user_data;

        if( nnz==NULL || local_row<0 || local_row>=2 ) return LIS_ERR_ILL_ARG;
        required = local_row==0 ? 3 : 1;
        if( cols==NULL && vals==NULL )
        {
                *nnz = required;
                return LIS_SUCCESS;
        }
        if( cols==NULL || vals==NULL || *nnz<required ) return LIS_ERR_ILL_ARG;

        if( local_row==0 )
        {
                cols[0]=0; scalar_set_real(&vals[0],(LIS_REAL)2.0);
                cols[1]=0; scalar_set_real(&vals[1],(LIS_REAL)-1.0);
                cols[2]=1; scalar_set_real(&vals[2],(LIS_REAL)3.0);
        }
        else
        {
                cols[0]=1; scalar_set_real(&vals[0],(LIS_REAL)1.0);
        }
        *nnz = required;
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
        ctx.saw_parallel = 0;

        /*
         * Row-only USER matrix. No explicit matvec, matvech or diagonal
         * callbacks are installed.
         */
        err = lis_matrix_set_user(
                A,
                &ctx,
                NULL,
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

        err = lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);

        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        /* Replacing the callback must clear the previous thread-safety opt-in. */
        err = lis_matrix_set_user_row(A,user_get_row);
        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }
        {
                LIS_INT capabilities = 0;
                err = lis_matrix_get_capabilities(A,&capabilities);
                if( err || (capabilities & LIS_MATRIX_CAP_ROW_THREAD_SAFE) )
                {
                        fprintf(stderr,"replacement USER row callback inherited thread-safe flag\n");
                        lis_matrix_destroy(A);
                        lis_finalize();
                        return 1;
                }
        }
        err = lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);
        if( err )
        {
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        {
                LIS_INT capabilities;
                LIS_INT expected;

                capabilities = 0;

                err = lis_matrix_get_capabilities(
                        A,
                        &capabilities);

                if( err )
                {
                        lis_matrix_destroy(A);
                        lis_finalize();
                        return 1;
                }

                expected = LIS_MATRIX_CAP_ROW_ACCESS |
                           LIS_MATRIX_CAP_ROW_THREAD_SAFE;

                if( local_n==global_n )
                {
                        expected |=
                                LIS_MATRIX_CAP_MATVEC
                                |
                                LIS_MATRIX_CAP_MATVECH
                                |
                                LIS_MATRIX_CAP_DIAGONAL;
                }

                if( capabilities!=expected )
                {
                        fprintf(
                                stderr,
                                "row-only capabilities: got=0x%x expected=0x%x\n",
                                (unsigned int)capabilities,
                                (unsigned int)expected);

                        lis_matrix_destroy(A);
                        lis_finalize();
                        return 1;
                }
        }

        if( local_n==global_n )
        {
                LIS_VECTOR x;
                LIS_VECTOR y;
                LIS_VECTOR yh;
                LIS_VECTOR d;
                LIS_SCALAR expected;
                LIS_INT global_row;

                x = NULL;
                y = NULL;
                yh = NULL;
                d = NULL;

                err = lis_vector_duplicate(A,&x);
                if( err ) goto serial_fallback_fail;

                err = lis_vector_duplicate(A,&y);
                if( err ) goto serial_fallback_fail;

                err = lis_vector_duplicate(A,&yh);
                if( err ) goto serial_fallback_fail;

                err = lis_vector_duplicate(A,&d);
                if( err ) goto serial_fallback_fail;

                for(i=0;i<local_n;i++)
                {
                        scalar_set_real(
                                &x->value[i],
                                (LIS_REAL)1.0);
                }

                err = lis_matrix_set_user_row_thread_safe(A,LIS_FALSE);
                if( err ) goto serial_fallback_fail;
                ctx.saw_parallel = 0;
                err = lis_matvec(A,x,y);
                if( err ) goto serial_fallback_fail;
#ifdef _OPENMP
                if( ctx.saw_parallel )
                {
                        fprintf(stderr,"non-thread-safe USER row callback entered OpenMP row path\n");
                        goto serial_fallback_fail;
                }
#endif

                err = lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);
                if( err ) goto serial_fallback_fail;
                ctx.saw_parallel = 0;
                err = lis_matvec(A,x,y);
                if( err ) goto serial_fallback_fail;
#ifdef _OPENMP
                if( omp_get_max_threads()>1 && !ctx.saw_parallel )
                {
                        fprintf(stderr,"thread-safe USER matvec did not use the OpenMP row path\n");
                        goto serial_fallback_fail;
                }
#endif

                err = lis_matvech(A,x,yh);
                if( err ) goto serial_fallback_fail;

                err = lis_matrix_get_diagonal(A,d);
                if( err ) goto serial_fallback_fail;

                {
                        LIS_REAL inf_norm;
                        LIS_REAL expected_norm;

                        inf_norm = (LIS_REAL)0.0;
                        err = lis_matrix_get_inf_norm(A,&inf_norm);
                        if( err ) goto serial_fallback_fail;

                        expected_norm =
                                global_n<=1 ? (LIS_REAL)2.0 :
                                global_n==2 ? (LIS_REAL)3.0 :
                                              (LIS_REAL)4.0;

                        if( fabs((double)(inf_norm-expected_norm))>1.0e-12 )
                        {
                                fprintf(stderr,
                                        "row-derived infinity norm mismatch: got=%.17g expected=%.17g\n",
                                        (double)inf_norm,(double)expected_norm);
                                goto serial_fallback_fail;
                        }
                }

                for(i=0;i<local_n;i++)
                {
                        global_row = is+i;

                        if( global_n==1 )
                        {
                                scalar_set_real(
                                        &expected,
                                        (LIS_REAL)2.0);
                        }
                        else if(
                                global_row==0
                                ||
                                global_row+1==global_n)
                        {
                                scalar_set_real(
                                        &expected,
                                        (LIS_REAL)1.0);
                        }
                        else
                        {
                                scalar_set_real(
                                        &expected,
                                        (LIS_REAL)0.0);
                        }

                        if(
                                !scalar_equal(&y->value[i],&expected)
                                ||
                                !scalar_equal(&yh->value[i],&expected))
                        {
                                fprintf(
                                        stderr,
                                        "row-derived matvec mismatch at row=%d\n",
                                        (int)i);

                                goto serial_fallback_fail;
                        }

                        scalar_set_real(
                                &expected,
                                (LIS_REAL)2.0);

                        if(
                                !scalar_equal(
                                        &d->value[i],
                                        &expected))
                        {
                                fprintf(
                                        stderr,
                                        "row-derived diagonal mismatch at row=%d\n",
                                        (int)i);

                                goto serial_fallback_fail;
                        }
                }

                lis_vector_destroy(d);
                lis_vector_destroy(yh);
                lis_vector_destroy(y);
                lis_vector_destroy(x);

                goto serial_fallback_done;

serial_fallback_fail:

                if( d ) lis_vector_destroy(d);
                if( yh ) lis_vector_destroy(yh);
                if( y ) lis_vector_destroy(y);
                if( x ) lis_vector_destroy(x);

                lis_matrix_destroy(A);
                lis_finalize();
                return 1;

serial_fallback_done:
                ;
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

        /* Infinity-norm regressions independent of the row-access API checks. */
        {
                LIS_MATRIX U=NULL,C=NULL,D=NULL,O=NULL;
                LIS_REAL norm;
                LIS_INT row_start,row_end,r,c;

                if( local_n==global_n )
                {
                        err = lis_matrix_create(LIS_COMM_WORLD,&U); if( err ) goto norm_fail;
                        err = lis_matrix_set_size(U,0,2); if( err ) goto norm_fail;
                        err = lis_matrix_set_user(U,NULL,NULL,NULL); if( err ) goto norm_fail;
                        err = lis_matrix_set_user_row(U,user_get_duplicate_row); if( err ) goto norm_fail;
                        norm = (LIS_REAL)0.0;
                        err = lis_matrix_get_inf_norm(U,&norm); if( err ) goto norm_fail;
                        if( fabs((double)(norm-(LIS_REAL)4.0))>1.0e-12 )
                        {
                                fprintf(stderr,"duplicate USER norm mismatch: %.17g\n",(double)norm);
                                goto norm_fail;
                        }
                        lis_matrix_destroy(U); U=NULL;
                }

                err = lis_matrix_create(LIS_COMM_WORLD,&C); if( err ) goto norm_fail;
                err = lis_matrix_set_size(C,0,4); if( err ) goto norm_fail;
                lis_matrix_get_range(C,&row_start,&row_end);
                for(r=row_start;r<row_end;r++)
                {
                        err = lis_matrix_set_value(LIS_INS_VALUE,r,r,(LIS_SCALAR)2.0,C); if( err ) goto norm_fail;
                        if(r>0){ err=lis_matrix_set_value(LIS_INS_VALUE,r,r-1,(LIS_SCALAR)-1.0,C); if(err) goto norm_fail; }
                        if(r+1<4){ err=lis_matrix_set_value(LIS_INS_VALUE,r,r+1,(LIS_SCALAR)-1.0,C); if(err) goto norm_fail; }
                }
                err = lis_matrix_set_type(C,LIS_MATRIX_CSR); if( err ) goto norm_fail;
                err = lis_matrix_assemble(C); if( err ) goto norm_fail;
                norm=(LIS_REAL)0.0;
                err=lis_matrix_get_inf_norm(C,&norm); if(err) goto norm_fail;
                if(fabs((double)(norm-(LIS_REAL)4.0))>1.0e-12) goto norm_fail;

                err=lis_matrix_duplicate(C,&D); if(err) goto norm_fail;
                err=lis_matrix_set_type(D,LIS_MATRIX_CSC); if(err) goto norm_fail;
                err=lis_matrix_convert(C,D); if(err) goto norm_fail;
                norm=(LIS_REAL)0.0;
                err=lis_matrix_get_inf_norm(D,&norm); if(err) goto norm_fail;
                if(fabs((double)(norm-(LIS_REAL)4.0))>1.0e-12) goto norm_fail;

                err=lis_matrix_create_operator((LIS_SCALAR)1.0,C,(LIS_SCALAR)0.0,C,&O); if(err) goto norm_fail;
                norm=(LIS_REAL)0.0;
                err=lis_matrix_get_inf_norm(O,&norm);
                if(err!=LIS_ERR_NOT_IMPLEMENTED) goto norm_fail;

                lis_matrix_destroy(O);
                lis_matrix_destroy(D);
                lis_matrix_destroy(C);
                goto norm_done;

norm_fail:
                if(O) lis_matrix_destroy(O);
                if(D) lis_matrix_destroy(D);
                if(C) lis_matrix_destroy(C);
                if(U) lis_matrix_destroy(U);
                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
norm_done:
                ;
        }

        lis_matrix_destroy(A);
        lis_finalize();

        return 0;
}
