/*
 * testgcrodrconfig.c
 *
 * GCRO-DR recycle-configuration lifecycle regression test.
 *
 * Covers:
 *   1. recycle dimension 0 disables recycling;
 *   2. recycle target >= restart is safely capped at restart-1;
 *   3. changing recycle dimension clears an existing recycle space;
 *   4. changing only quality thresholds preserves the recycle space.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lis.h"

#ifdef _COMPLEX
int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;
    printf("SKIP: GCRO-DR is not implemented for complex scalar builds\n");
    return 77;
}
#else

static LIS_INT build_matrix(LIS_Comm comm, LIS_INT m, LIS_INT n, LIS_MATRIX *Aout)
{
    LIS_MATRIX A;
    LIS_INT *ptr, *index;
    LIS_SCALAR *value;
    LIS_INT is, ie, ii, i, j, ctr, err;

    err = lis_matrix_create(comm,&A);
    if( err ) return err;

    err = lis_matrix_set_size(A,0,m*n);
    if( err )
    {
        lis_matrix_destroy(A);
        return err;
    }

    ptr   = (LIS_INT *)malloc((size_t)(A->n+1)*sizeof(LIS_INT));
    index = (LIS_INT *)malloc((size_t)(5*A->n)*sizeof(LIS_INT));
    value = (LIS_SCALAR *)malloc((size_t)(5*A->n)*sizeof(LIS_SCALAR));
    if( ptr==NULL || index==NULL || value==NULL )
    {
        free(ptr);
        free(index);
        free(value);
        lis_matrix_destroy(A);
        return LIS_ERR_OUT_OF_MEMORY;
    }

    lis_matrix_get_range(A,&is,&ie);

    ctr=0;
    ptr[0]=0;
    for(ii=is;ii<ie;ii++)
    {
        i=ii/m;
        j=ii-i*m;

        if( i>0 )
        {
            index[ctr]=ii-m;
            value[ctr++]=-1.0;
        }
        if( i<n-1 )
        {
            index[ctr]=ii+m;
            value[ctr++]=-1.0;
        }
        if( j>0 )
        {
            index[ctr]=ii-1;
            value[ctr++]=-1.0;
        }
        if( j<m-1 )
        {
            index[ctr]=ii+1;
            value[ctr++]=-1.0;
        }

        index[ctr]=ii;
        value[ctr++]=4.1;
        ptr[ii-is+1]=ctr;
    }

    err = lis_matrix_set_csr(ptr[ie-is],ptr,index,value,A);
    if( err )
    {
        free(ptr);
        free(index);
        free(value);
        lis_matrix_destroy(A);
        return err;
    }

    err = lis_matrix_assemble(A);
    if( err )
    {
        lis_matrix_destroy(A);
        return err;
    }

    *Aout=A;
    return LIS_SUCCESS;
}

static void set_exact_solution(LIS_VECTOR u)
{
    LIS_INT i;
    for(i=0;i<u->n;i++)
    {
        LIS_REAL g=(LIS_REAL)(u->is+i+1);
        u->value[i]=(LIS_SCALAR)(
            1.0 + 0.10*sin(0.031*g) + 0.04*cos(0.017*g)
        );
    }
}

static LIS_INT solve_once(
    LIS_Comm comm,
    LIS_SOLVER solver,
    LIS_INT m,
    LIS_INT n,
    LIS_REAL *true_resid)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,b=NULL,u=NULL,r=NULL;
    LIS_INT err,status;
    LIS_REAL rnrm,bnrm;

    err=build_matrix(comm,m,n,&A);
    if( err ) goto cleanup;

    if( (err=lis_vector_duplicate(A,&u)) ) goto cleanup;
    if( (err=lis_vector_duplicate(A,&b)) ) goto cleanup;
    if( (err=lis_vector_duplicate(A,&x)) ) goto cleanup;
    if( (err=lis_vector_duplicate(A,&r)) ) goto cleanup;

    set_exact_solution(u);
    lis_matvec(A,u,b);
    lis_vector_set_all(0.0,x);

    err=lis_solve(A,b,x,solver);
    if( err ) goto cleanup;

    lis_solver_get_status(solver,&status);
    if( status!=LIS_SUCCESS )
    {
        err=status;
        goto cleanup;
    }

    lis_matvec(A,x,r);
    lis_vector_xpay(b,-1.0,r);
    lis_vector_nrm2(r,&rnrm);
    lis_vector_nrm2(b,&bnrm);

    *true_resid=(bnrm>0.0) ? rnrm/bnrm : rnrm;
    if( *true_resid>1.0e-8 )
        err=LIS_BREAKDOWN;
    else
        err=LIS_SUCCESS;

cleanup:
    if( r ) lis_vector_destroy(r);
    if( x ) lis_vector_destroy(x);
    if( b ) lis_vector_destroy(b);
    if( u ) lis_vector_destroy(u);
    if( A ) lis_matrix_destroy(A);
    return err;
}

static int fail_case(LIS_Comm comm, const char *name, const char *message)
{
    lis_printf(comm,"CASE %-20s : FAIL (%s)\n",name,message);
    return 1;
}

int main(int argc, char *argv[])
{
    LIS_Comm comm;
    LIS_SOLVER solver=NULL;
    LIS_INT err,dim,action;
    LIS_INT dim_before;
    LIS_REAL quality,qmin,resid;
    LIS_REAL q_before,qmin_before;
    int failed=0;

    const LIS_INT m=24;
    const LIS_INT n=24;

    lis_initialize(&argc,&argv);
    comm=LIS_COMM_WORLD;

    err=lis_solver_create(&solver);
    if( err )
    {
        lis_finalize();
        return (int)err;
    }

    err=lis_solver_set_option(
        "-i gcrodr -p jacobi -restart 8 -maxiter 2000 -tol 1e-10 "
        "-gcrodr_k 0 -gcrodr_keep 0.90 -gcrodr_reject 0.35 "
        "-print none",
        solver);
    if( err )
    {
        lis_solver_destroy(solver);
        lis_finalize();
        return (int)err;
    }

    /*
     * CASE 1: k=0 disables recycling while leaving GCRO-DR usable.
     */
    err=solve_once(comm,solver,m,n,&resid);
    if( err )
    {
        failed |= fail_case(comm,"k=0","solve failed");
    }
    else
    {
        lis_solver_get_gcrodr_recycle(solver,&dim,&action,&quality,&qmin);
        if( dim!=0 || solver->recycle_valid )
            failed |= fail_case(comm,"k=0","recycle space was created");
        else
            lis_printf(comm,
                "CASE %-20s : PASS dim=%D true_resid=%.3e\n",
                "k=0",dim,(double)resid);
    }

    /*
     * CASE 2: a target larger than restart is accepted, but the effective
     * recycle space is capped at restart-1.
     */
    err=lis_solver_set_option("-gcrodr_k 20",solver);
    if( err )
    {
        failed |= fail_case(comm,"k>=restart","option rejected");
    }
    else
    {
        err=solve_once(comm,solver,m,n,&resid);
        if( err )
        {
            failed |= fail_case(comm,"k>=restart","solve failed");
        }
        else
        {
            lis_solver_get_gcrodr_recycle(solver,&dim,&action,&quality,&qmin);
            if( solver->recycle_target_dim!=20 || dim<=0 || dim>7 )
                failed |= fail_case(comm,"k>=restart","effective dimension not capped");
            else
                lis_printf(comm,
                    "CASE %-20s : PASS target=%D restart=%D dim=%D action=%D\n",
                    "k>=restart",
                    solver->recycle_target_dim,
                    solver->options[LIS_OPTIONS_RESTART],
                    dim,action);
        }
    }

    /*
     * CASE 3: changing k invalidates the old recycle space immediately.
     */
    dim_before=solver->recycle_dim;
    err=lis_solver_set_option("-gcrodr_k 4",solver);
    if( err )
    {
        failed |= fail_case(comm,"change-k-clears","option rejected");
    }
    else
    {
        lis_solver_get_gcrodr_recycle(solver,&dim,&action,&quality,&qmin);
        if( dim!=0 || solver->recycle_valid )
        {
            failed |= fail_case(comm,"change-k-clears","old recycle space survived");
        }
        else
        {
            err=solve_once(comm,solver,m,n,&resid);
            if( err )
            {
                failed |= fail_case(comm,"change-k-clears","rebuild solve failed");
            }
            else
            {
                lis_solver_get_gcrodr_recycle(
                    solver,&dim,&action,&quality,&qmin);
                if( dim<=0 || dim>4 ||
                    action!=LIS_GCRODR_RECYCLE_BUILD )
                {
                    failed |= fail_case(
                        comm,"change-k-clears","next solve did not rebuild");
                }
                else
                {
                    lis_printf(comm,
                        "CASE %-20s : PASS old_dim=%D cleared=0 new_dim=%D action=BUILD\n",
                        "change-k-clears",dim_before,dim);
                }
            }
        }
    }

    /*
     * CASE 4: changing only policy thresholds must not discard U/ZU/C.
     * Re-solving the identical operator should therefore reuse the space.
     */
    dim_before=solver->recycle_dim;
    q_before=solver->recycle_quality;
    qmin_before=solver->recycle_quality_min;

    err=lis_solver_set_option(
        "-gcrodr_keep 0.99 -gcrodr_reject 0.95",
        solver);
    if( err )
    {
        failed |= fail_case(comm,"threshold-preserve","option rejected");
    }
    else
    {
        lis_solver_get_gcrodr_recycle(solver,&dim,&action,&quality,&qmin);

        if( dim!=dim_before || !solver->recycle_valid )
        {
            failed |= fail_case(
                comm,"threshold-preserve","recycle space was cleared");
        }
        else if( solver->recycle_keep_threshold!=(LIS_REAL)0.99 ||
                 solver->recycle_reject_threshold!=(LIS_REAL)0.95 )
        {
            failed |= fail_case(
                comm,"threshold-preserve","thresholds were not updated");
        }
        else
        {
            err=solve_once(comm,solver,m,n,&resid);
            if( err )
            {
                failed |= fail_case(
                    comm,"threshold-preserve","reuse solve failed");
            }
            else
            {
                lis_solver_get_gcrodr_recycle(
                    solver,&dim,&action,&quality,&qmin);

                if( action!=LIS_GCRODR_RECYCLE_KEEP )
                {
                    failed |= fail_case(
                        comm,"threshold-preserve","identical operator was not kept");
                }
                else
                {
                    lis_printf(comm,
                        "CASE %-20s : PASS dim=%D action=KEEP q=%.4f qmin=%.4f\n",
                        "threshold-preserve",
                        dim,(double)quality,(double)qmin);
                }
            }
        }
    }

    /*
     * CASE 5: GCRO-DR currently supports only the default convergence
     * condition.  A non-default condition must be rejected cleanly.
     */
    err=lis_solver_set_option("-conv_cond nrm2_b",solver);
    if( err )
    {
        failed |= fail_case(comm,"conv-cond-guard","option parse failed");
    }
    else
    {
        err=solve_once(comm,solver,m,n,&resid);
        if( err!=LIS_ERR_ILL_ARG )
            failed |= fail_case(comm,"conv-cond-guard","unsupported condition was not rejected");
        else
            lis_printf(comm,
                "CASE %-20s : PASS non-default condition rejected\n",
                "conv-cond-guard");
    }

    /*
     * Keep otherwise-unused pre-change quality values visible to the
     * compiler in non-optimizing/debug builds; the important preserved
     * state in this case is the recycle space itself.
     */
    (void)q_before;
    (void)qmin_before;

    if( failed )
        lis_printf(comm,"\nGCRO-DR configuration lifecycle test: FAIL\n");
    else
        lis_printf(comm,"\nGCRO-DR configuration lifecycle test: PASS\n");

    lis_solver_destroy(solver);
    lis_finalize();
    return failed ? 1 : 0;
}
#endif /* !_COMPLEX */
