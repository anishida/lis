#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lis.h"
#include "lis_precon.h"
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

static int scalar_close(LIS_SCALAR a, LIS_SCALAR b)
{
        LIS_SCALAR diff = a-b;
        LIS_REAL tol = (LIS_REAL)1.0e-10;

#ifdef _COMPLEX
        LIS_REAL dr = (LIS_REAL)creal(diff);
        LIS_REAL di = (LIS_REAL)cimag(diff);

        if( dr<0 ) dr = -dr;
        if( di<0 ) di = -di;

        return dr<=tol && di<=tol;
#else
        LIS_REAL d = (LIS_REAL)diff;

        if( d<0 ) d = -d;
        return d<=tol;
#endif
}

static int check_precon_diag(const char *name,
                             LIS_PRECON precon,
                             LIS_SCALAR expected)
{
        LIS_INT i;

        if( precon==NULL || precon->D==NULL )
        {
                fprintf(stderr,"%s: missing Jacobi diagonal\n",name);
                return 1;
        }

        for(i=0;i<precon->D->n;i++)
        {
                if( !scalar_close(precon->D->value[i],expected) )
                {
                        fprintf(stderr,
                                "%s: unexpected inverse diagonal at %d\n",
                                name,(int)i);
                        return 1;
                }
        }

        return 0;
}

static LIS_INT create_jacobi(LIS_SOLVER solver,
                             LIS_MATRIX A,
                             LIS_PRECON *precon)
{
        LIS_INT err;

        solver->A = A;

        err = lis_solver_set_option("-p jacobi",solver);
        if( err ) return err;

        *precon = NULL;
        err = lis_precon_create(solver,precon);
        if( err )
        {
                *precon = NULL;
                return err;
        }

        return LIS_SUCCESS;
}

int main(int argc, char *argv[])
{
        LIS_MATRIX U = NULL;
        LIS_MATRIX V = NULL;
        LIS_MATRIX C = NULL;
        LIS_MATRIX N = NULL;
        LIS_SOLVER solver = NULL;
        LIS_PRECON precon = NULL;
        USER_DIAG uctx,vctx,nctx;
        LIS_INT is,ie,err;
        LIS_INT global_n = 4;
        LIS_INT local_fail = 0;
        LIS_INT global_fail = 0;

        lis_initialize(&argc,&argv);

        err = lis_solver_create(&solver);
        if( err ) goto api_fail;

        /*
         * U = 2*I with a valid USER diagonal callback.
         */
        err = lis_matrix_create(LIS_COMM_WORLD,&U);
        if( err ) goto api_fail;

        if( U->nprocs>global_n )
                global_n = U->nprocs;

        err = lis_matrix_set_size(U,0,global_n);
        if( err ) goto api_fail;

        err = lis_matrix_get_range(U,&is,&ie);
        if( err ) goto api_fail;

        uctx.n = ie-is;
        uctx.diag = (LIS_SCALAR)2.0;

        err = lis_matrix_set_user(
                U,&uctx,user_matvec,user_matvech);
        if( err ) goto api_fail;

        err = lis_matrix_set_user_diagonal(
                U,user_get_diagonal);
        if( err ) goto api_fail;

        err = create_jacobi(solver,U,&precon);
        if( err ) goto api_fail;

        local_fail |= check_precon_diag(
                "USER Jacobi",
                precon,
                (LIS_SCALAR)0.5);

        lis_precon_destroy(precon);
        precon = NULL;

        /*
         * V = 4*I with a valid diagonal callback.
         */
        err = lis_matrix_create(LIS_COMM_WORLD,&V);
        if( err ) goto api_fail;

        err = lis_matrix_set_size(V,0,global_n);
        if( err ) goto api_fail;

        err = lis_matrix_get_range(V,&is,&ie);
        if( err ) goto api_fail;

        vctx.n = ie-is;
        vctx.diag = (LIS_SCALAR)4.0;

        err = lis_matrix_set_user(
                V,&vctx,user_matvec,user_matvech);
        if( err ) goto api_fail;

        err = lis_matrix_set_user_diagonal(
                V,user_get_diagonal);
        if( err ) goto api_fail;

        /*
         * C = 0.5*U + 0.5*V = 3*I.
         *
         * Jacobi must therefore contain 1/3 on its diagonal.
         */
        err = lis_matrix_create_operator(
                (LIS_SCALAR)0.5,U,
                (LIS_SCALAR)0.5,V,&C);
        if( err ) goto api_fail;

        err = create_jacobi(solver,C,&precon);
        if( err ) goto api_fail;

        local_fail |= check_precon_diag(
                "OPERATOR Jacobi",
                precon,
                (LIS_SCALAR)(1.0/3.0));

        lis_precon_destroy(precon);
        precon = NULL;

        /*
         * N = 7*I as USER without a diagonal callback.
         *
         * Jacobi creation must propagate LIS_ERR_NOT_IMPLEMENTED
         * from lis_matrix_get_diagonal().
         */
        err = lis_matrix_create(LIS_COMM_WORLD,&N);
        if( err ) goto api_fail;

        err = lis_matrix_set_size(N,0,global_n);
        if( err ) goto api_fail;

        err = lis_matrix_get_range(N,&is,&ie);
        if( err ) goto api_fail;

        nctx.n = ie-is;
        nctx.diag = (LIS_SCALAR)7.0;

        err = lis_matrix_set_user(
                N,&nctx,user_matvec,user_matvech);
        if( err ) goto api_fail;

        err = create_jacobi(solver,N,&precon);

        if( err!=LIS_ERR_NOT_IMPLEMENTED )
        {
                fprintf(stderr,
                        "USER without diagonal callback: "
                        "expected %d, got %d\n",
                        (int)LIS_ERR_NOT_IMPLEMENTED,
                        (int)err);
                local_fail = 1;

                if( precon )
                {
                        lis_precon_destroy(precon);
                        precon = NULL;
                }
        }
        else
        {
                /*
                 * lis_precon_create() destroys the partially created
                 * preconditioner when creation fails.
                 */
                precon = NULL;
        }

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

        if( U->my_rank==0 )
        {
                printf("LIS_JACOBI_USER_OPERATOR %s\n",
                       global_fail ? "FAILED" : "PASSED");
        }

        solver->A = NULL;

        if( N ) lis_matrix_destroy(N);
        if( C ) lis_matrix_destroy(C);
        if( V ) lis_matrix_destroy(V);
        if( U ) lis_matrix_destroy(U);
        if( solver ) lis_solver_destroy(solver);

        lis_finalize();
        return global_fail ? 1 : 0;

api_fail:
        fprintf(stderr,
                "LIS_JACOBI_USER_OPERATOR API failure: %d\n",
                (int)err);

        if( precon ) lis_precon_destroy(precon);

        if( solver )
                solver->A = NULL;

        if( N ) lis_matrix_destroy(N);
        if( C ) lis_matrix_destroy(C);
        if( V ) lis_matrix_destroy(V);
        if( U ) lis_matrix_destroy(U);
        if( solver ) lis_solver_destroy(solver);

        lis_finalize();
        return 1;
}
