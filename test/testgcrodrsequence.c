/*
 * testgcrodrsequence.c
 *
 * Regression test for persistent GCRO-DR recycling across a short sequence
 * of closely related linear systems configured through textual LIS options.
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

static LIS_INT build_matrix(
    LIS_Comm comm, LIS_INT m, LIS_INT n, LIS_REAL shift, LIS_MATRIX *Aout)
{
    LIS_MATRIX A;
    LIS_INT *ptr,*index;
    LIS_SCALAR *value;
    LIS_INT is,ie,ii,i,j,ctr,err;

    err=lis_matrix_create(comm,&A);
    if( err ) return err;

    err=lis_matrix_set_size(A,0,m*n);
    if( err )
    {
        lis_matrix_destroy(A);
        return err;
    }

    ptr=(LIS_INT *)malloc((size_t)(A->n+1)*sizeof(LIS_INT));
    index=(LIS_INT *)malloc((size_t)(5*A->n)*sizeof(LIS_INT));
    value=(LIS_SCALAR *)malloc((size_t)(5*A->n)*sizeof(LIS_SCALAR));
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
        value[ctr++]=(LIS_SCALAR)(4.1+shift);
        ptr[ii-is+1]=ctr;
    }

    err=lis_matrix_set_csr(ptr[ie-is],ptr,index,value,A);
    if( err )
    {
        free(ptr);
        free(index);
        free(value);
        lis_matrix_destroy(A);
        return err;
    }

    err=lis_matrix_assemble(A);
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
            1.0 + 0.08*sin(0.021*g) + 0.03*cos(0.037*g)
        );
    }
}

static LIS_INT run_system(
    LIS_Comm comm,
    LIS_SOLVER solver,
    LIS_INT m,
    LIS_INT n,
    LIS_REAL shift,
    LIS_REAL *true_resid)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR u=NULL,b=NULL,x=NULL,r=NULL;
    LIS_REAL rnrm,bnrm;
    LIS_INT err,status;

    err=build_matrix(comm,m,n,shift,&A);
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
    err=(*true_resid<=1.0e-8) ? LIS_SUCCESS : LIS_BREAKDOWN;

cleanup:
    if( r ) lis_vector_destroy(r);
    if( x ) lis_vector_destroy(x);
    if( b ) lis_vector_destroy(b);
    if( u ) lis_vector_destroy(u);
    if( A ) lis_matrix_destroy(A);
    return err;
}

int main(int argc, char *argv[])
{
    LIS_Comm comm;
    LIS_SOLVER solver=NULL;
    LIS_INT err,dim,action;
    LIS_REAL q,qmin,resid;
    LIS_INT cross0,cross1;
    const LIS_INT m=24,n=24;

    lis_initialize(&argc,&argv);
    comm=LIS_COMM_WORLD;

    err=lis_solver_create(&solver);
    if( err ) goto fail_init;

    err=lis_solver_set_option(
        "-i gcrodr -p jacobi -restart 20 -maxiter 2000 -tol 1e-10 "
        "-gcrodr_k 4 -gcrodr_keep 0.80 -gcrodr_reject 0.20 "
        "-print none",
        solver);
    if( err ) goto fail;

    err=lis_solver_set_optionC(solver);
    if( err ) goto fail;

    if( solver->recycle_target_dim!=4 ||
        fabs((double)(solver->recycle_keep_threshold-(LIS_REAL)0.80))>1.0e-14 ||
        fabs((double)(solver->recycle_reject_threshold-(LIS_REAL)0.20))>1.0e-14 )
    {
        lis_printf(comm,"GCRO-DR sequence test: FAIL (text options)\n");
        goto fail;
    }

    err=run_system(comm,solver,m,n,(LIS_REAL)0.000,&resid);
    if( err ) goto fail;

    lis_solver_get_gcrodr_recycle(solver,&dim,&action,&q,&qmin);
    if( dim<=0 || dim>4 || action!=LIS_GCRODR_RECYCLE_BUILD )
    {
        lis_printf(comm,
            "GCRO-DR sequence test: FAIL (first solve dim=%D action=%D)\n",
            dim,action);
        goto fail;
    }
    cross0=solver->recycle_cross_uses;

    lis_printf(comm,
        "solve=0 action=BUILD dim=%D q=%.4f true_resid=%.3e cross=%D\n",
        dim,(double)q,(double)resid,cross0);

    err=run_system(comm,solver,m,n,(LIS_REAL)0.002,&resid);
    if( err ) goto fail;

    lis_solver_get_gcrodr_recycle(solver,&dim,&action,&q,&qmin);
    cross1=solver->recycle_cross_uses;

    if( dim<=0 || dim>4 ||
        action!=LIS_GCRODR_RECYCLE_KEEP ||
        cross1<=cross0 )
    {
        lis_printf(comm,
            "GCRO-DR sequence test: FAIL "
            "(second solve dim=%D action=%D cross=%D->%D q=%.4f)\n",
            dim,action,cross0,cross1,(double)q);
        goto fail;
    }

    lis_printf(comm,
        "solve=1 action=KEEP dim=%D q=%.4f qmin=%.4f "
        "true_resid=%.3e cross=%D\n",
        dim,(double)q,(double)qmin,(double)resid,cross1);

    cross0=cross1;
    err=run_system(comm,solver,m,n,(LIS_REAL)0.004,&resid);
    if( err ) goto fail;

    lis_solver_get_gcrodr_recycle(solver,&dim,&action,&q,&qmin);
    cross1=solver->recycle_cross_uses;

    if( dim<=0 || dim>4 ||
        action!=LIS_GCRODR_RECYCLE_KEEP ||
        cross1<=cross0 )
    {
        lis_printf(comm,
            "GCRO-DR sequence test: FAIL "
            "(third solve dim=%D action=%D cross=%D->%D q=%.4f)\n",
            dim,action,cross0,cross1,(double)q);
        goto fail;
    }

    lis_printf(comm,
        "solve=2 action=KEEP dim=%D q=%.4f qmin=%.4f "
        "true_resid=%.3e cross=%D\n",
        dim,(double)q,(double)qmin,(double)resid,cross1);

    lis_printf(comm,"GCRO-DR persistent sequence test: PASS\n");

    lis_solver_destroy(solver);
    lis_finalize();
    return 0;

fail:
    if( solver ) lis_solver_destroy(solver);
fail_init:
    lis_finalize();
    return 1;
}
#endif /* !_COMPLEX */
