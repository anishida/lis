/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Mutable lifecycle helpers for non-owning matrix operators.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lislib.h"
#include "lis_matrix_operator.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

static LIS_INT
lis_matrix_operator_same_distribution(LIS_MATRIX A,LIS_MATRIX B)
{
    LIS_INT i;
#ifdef USE_MPI
    int mpi_err,comm_result;
#endif

    if(A==NULL || B==NULL)
        return LIS_FALSE;
    if(A->gn!=B->gn || A->n!=B->n || A->is!=B->is || A->ie!=B->ie ||
       A->nprocs!=B->nprocs)
        return LIS_FALSE;
    if((A->ranges==NULL)!=(B->ranges==NULL))
        return LIS_FALSE;
    if(A->ranges!=NULL)
        for(i=0;i<=A->nprocs;i++)
            if(A->ranges[i]!=B->ranges[i])
                return LIS_FALSE;
#ifdef USE_MPI
    mpi_err=MPI_Comm_compare(A->comm,B->comm,&comm_result);
    if(mpi_err!=MPI_SUCCESS ||
       (comm_result!=MPI_IDENT && comm_result!=MPI_CONGRUENT))
        return LIS_FALSE;
#endif
    return LIS_TRUE;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_operator_set_coefficients"
LIS_INT
lis_matrix_operator_set_coefficients(
    LIS_MATRIX C,
    LIS_SCALAR alpha,
    LIS_SCALAR beta)
{
    LIS_INT err;

    err=lis_matrix_check(C,LIS_MATRIX_CHECK_ALL);
    if(err) return err;
    if(C->matrix_type!=LIS_MATRIX_OPERATOR || C->operator_scale!=NULL ||
       C->operator_A==NULL || C->operator_B==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "matrix is not a linear-combination operator\n");
        return LIS_ERR_ILL_ARG;
    }

    C->operator_alpha=alpha;
    C->operator_beta=beta;
    return lis_matrix_mark_values_changed(C);
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_operator_rebind"
LIS_INT
lis_matrix_operator_rebind(
    LIS_MATRIX C,
    LIS_MATRIX A,
    LIS_MATRIX B)
{
    LIS_MATRIX candidate=NULL;
    LIS_VECTOR old_work;
    LIS_INT err;

    err=lis_matrix_check(C,LIS_MATRIX_CHECK_ALL);
    if(err) return err;
    if(C->matrix_type!=LIS_MATRIX_OPERATOR || C->operator_scale!=NULL ||
       C->operator_A==NULL || C->operator_B==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "matrix is not a linear-combination operator\n");
        return LIS_ERR_ILL_ARG;
    }

    /* Reuse the constructor as the validation/allocation transaction. */
    err=lis_matrix_create_operator(
        C->operator_alpha,A,C->operator_beta,B,&candidate);
    if(err) return err;

    if(!lis_matrix_operator_same_distribution(C,candidate))
    {
        lis_matrix_destroy(candidate);
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "replacement operator changes the matrix distribution\n");
        return LIS_ERR_ILL_ARG;
    }

    old_work=C->operator_work;
    C->operator_A=A;
    C->operator_B=B;
    C->np=candidate->np;
    C->pad=candidate->pad;
    C->origin=candidate->origin;
    C->operator_work=candidate->operator_work;
    candidate->operator_work=NULL;
    lis_matrix_destroy(candidate);
    if(old_work) lis_vector_destroy(old_work);

    return lis_matrix_mark_structure_changed(C);
}
