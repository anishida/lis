#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include "lis_matrix_operator.h"
#include "lis_matrix_materialize.h"
#include <stdio.h>
#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif

/*
 * Run with mpirun -np 2. This exercises operator rebinds on distributed
 * CSR operands without requiring distributed USER/CSR materialization.
 */
static LIS_INT make_diagonal(LIS_Comm comm,LIS_INT n,LIS_SCALAR diagonal,
                             LIS_MATRIX *out)
{
    LIS_MATRIX A=NULL;
    LIS_INT i,is,ie,err;
    *out=NULL;
    err=lis_matrix_create(comm,&A);if(err)return err;
    err=lis_matrix_set_size(A,0,n);if(err)goto fail;
    err=lis_matrix_get_range(A,&is,&ie);if(err)goto fail;
    for(i=is;i<ie;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,diagonal,A);
        if(err)goto fail;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto fail;
    err=lis_matrix_assemble(A);if(err)goto fail;
    *out=A;
    return LIS_SUCCESS;
fail:
    lis_matrix_destroy(A);
    return err;
}
static int check_output(LIS_MATRIX op,LIS_VECTOR x,LIS_VECTOR y,
                        LIS_REAL factor)
{
    LIS_INT i,err;
    err=lis_matvec(op,x,y);
    if(err)return 1;
    for(i=0;i<op->n;i++)
    {
#ifdef _COMPLEX
        double error=cabs(y->value[i]-(LIS_SCALAR)factor*x->value[i]);
#else
        double error=fabs((double)(y->value[i]-(LIS_SCALAR)factor*x->value[i]));
#endif
        if(error>1e-11)return 1;
    }
    return 0;
}
#if defined(USE_MPI) && !defined(_COMPLEX) && !defined(USE_COMPLEX) && !defined(_LONG__DOUBLE)
static int check_distributed_pair(LIS_MATRIX A,LIS_MATRIX P,int rank)
{
    LIS_MATRIX_AUTO_SESSION session=NULL;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR b=NULL,x=NULL,r=NULL;
    LIS_PRECON kept=NULL;
    LIS_REAL nr,nb;
    LIS_INT err,status,pass;
    int failed=1,local,all;
    local=lis_matrix_auto_session_create(&session) ||
          lis_matrix_auto_session_get_solver(session,&solver);
    MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(all)goto done;
    local=lis_solver_set_option("-i fgmres -p jacobi -scale none -tol 1e-10 -restart 10 -maxiter 80 -print none",solver);
    MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(all)goto done;
    /* Matrix/vector creation follows native collective allocation rules. */
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&r) || lis_vector_set_all((LIS_SCALAR)1,x) ||
       lis_matvec(A,x,b) || lis_vector_nrm2(b,&nb) || !(nb>0))
        goto done;
    for(pass=0;pass<3;pass++)
    {
        if(pass==2)
        {
            /* Only one rank publishes a new generation. All ranks
             * must retire persistent state before the next solve. */
            if(rank==0)lis_matrix_mark_values_changed(P);
            kept=NULL;
        }
        if(lis_vector_set_all((LIS_SCALAR)0,x))goto done;
        err=lis_matrix_auto_session_solve_distributed_csr(session,A,P,b,x);
        local=err || lis_solver_get_status(solver,&status) ||
              status!=LIS_SUCCESS || !solver->sequence_precon ||
              solver->sequence_precon_matrix!=P ||
              (kept && solver->sequence_precon!=kept);
        MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        if(all)goto done;
        if(lis_matvec(A,x,r) || lis_vector_axpy((LIS_SCALAR)-1,b,r) ||
           lis_vector_nrm2(r,&nr) || !(nr/nb<1e-9))goto done;
        kept=solver->sequence_precon;
    }
    /* A rank-local invalid argument must not enter the numerical solver. */
    err=lis_matrix_auto_session_solve_distributed_csr(
                    session,A,P,rank==0?NULL:b,x);
    local=err!=LIS_ERR_ILL_ARG;
    MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(all)goto done;
    failed=0;
done:
    lis_matrix_auto_session_destroy(session);
    if(r)lis_vector_destroy(r);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(rank==0)printf("MPI_PAIRED_CSR_JACOBI %s\n",failed?"FAILED":"PASSED");
    return failed;
}
#endif

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,B=NULL,C=NULL,D=NULL,E=NULL,op=NULL,nested=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_UNSIGNED_INT old_s=0,old_v=0,new_s=0,new_v=0;
    LIS_INT i,is,ie,err,n=8;
    int failed=0,rank=0,size=1,global_failed=0;
#ifdef USE_MPI
    MPI_Comm reverse_comm=MPI_COMM_NULL;
    LIS_MATRIX reversed=NULL;
#endif
    err=lis_initialize(&argc,&argv);
    if(err)return 1;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
    err=make_diagonal(LIS_COMM_WORLD,n,(LIS_SCALAR)2.0,&A);
    if(err){failed=1;goto done;}
    err=make_diagonal(LIS_COMM_WORLD,n,(LIS_SCALAR)3.0,&B);
    if(err){failed=1;goto done;}
    err=make_diagonal(LIS_COMM_WORLD,n,(LIS_SCALAR)4.0,&C);
    if(err){failed=1;goto done;}
    err=make_diagonal(LIS_COMM_WORLD,n,(LIS_SCALAR)5.0,&D);
    if(err){failed=1;goto done;}
    err=make_diagonal(LIS_COMM_WORLD,n+2,(LIS_SCALAR)1.0,&E);
    if(err){failed=1;goto done;}
#if defined(USE_MPI) && !defined(_COMPLEX) && !defined(USE_COMPLEX) && !defined(_LONG__DOUBLE)
    if(size>1 && check_distributed_pair(A,B,rank)){failed=1;goto done;}
#endif
    err=lis_matrix_create_operator(
        (LIS_SCALAR)1.0,A,(LIS_SCALAR)1.0,B,&op);
    if(err){failed=1;goto done;}
    err=lis_matrix_create_operator(
        (LIS_SCALAR)0.5,op,(LIS_SCALAR)0.5,A,&nested);
    if(err){failed=1;goto done;}
    err=lis_vector_duplicate(op,&x);if(err){failed=1;goto done;}
    err=lis_vector_duplicate(op,&y);if(err){failed=1;goto done;}
    err=lis_matrix_get_range(op,&is,&ie);
    if(err){failed=1;goto done;}
    for(i=0;i<op->n;i++)
        x->value[i]=(LIS_SCALAR)(1.0+0.1*(double)(is+i));
    if(check_output(op,x,y,5.0))
    {
        fprintf(stderr,"rank %d initial distributed operator mismatch\n",rank);
        failed=1;goto done;
    }
    err=lis_matrix_operator_set_coefficients(
        op,(LIS_SCALAR)2.0,(LIS_SCALAR)-1.0);
    if(err||check_output(op,x,y,1.0)){failed=1;goto done;}
    err=lis_matrix_operator_rebind(op,C,D);
    if(err||check_output(op,x,y,3.0)){failed=1;goto done;}
    err=lis_matrix_get_generation(op,&old_s,&old_v);
    if(err){failed=1;goto done;}
    if(lis_matrix_operator_rebind(op,op,D)==LIS_SUCCESS ||
       lis_matrix_operator_rebind(op,nested,D)==LIS_SUCCESS ||
       lis_matrix_operator_rebind(op,E,D)==LIS_SUCCESS)
    {
        fprintf(stderr,"rank %d invalid distributed rebind accepted\n",rank);
        failed=1;goto done;
    }
#ifdef USE_MPI
    if(size>1)
    {
        int mpi_err=MPI_Comm_split(
            MPI_COMM_WORLD,0,size-1-rank,&reverse_comm);
        if(mpi_err!=MPI_SUCCESS){failed=1;goto done;}
        err=make_diagonal(reverse_comm,n,(LIS_SCALAR)7.0,&reversed);
        if(err){failed=1;goto done;}
        /* Same members but a different rank order (MPI_SIMILAR). */
        if(lis_matrix_operator_rebind(op,reversed,D)==LIS_SUCCESS)
        {
            fprintf(stderr,"rank %d incompatible communicator accepted\n",rank);
            failed=1;goto done;
        }
    }
#endif
    err=lis_matrix_get_generation(op,&new_s,&new_v);
    if(err||new_s!=old_s||new_v!=old_v||check_output(op,x,y,3.0))
    {
        fprintf(stderr,"rank %d distributed rollback mismatch\n",rank);
        failed=1;goto done;
    }
done:
#ifdef USE_MPI
    if(reversed)lis_matrix_destroy(reversed);
    if(reverse_comm!=MPI_COMM_NULL)MPI_Comm_free(&reverse_comm);
#endif
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(nested)lis_matrix_destroy(nested);
    if(op)lis_matrix_destroy(op);
    if(E)lis_matrix_destroy(E);
    if(D)lis_matrix_destroy(D);
    if(C)lis_matrix_destroy(C);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
#ifdef USE_MPI
    MPI_Allreduce(&failed,&global_failed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#else
    global_failed=failed;
#endif
    lis_finalize();
    if(rank==0)
        printf("LIS_MATRIX_OPERATOR_MPI %s (ranks=%d)\n",
               global_failed?"FAILED":"PASSED",size);
    return global_failed?1:0;
}
