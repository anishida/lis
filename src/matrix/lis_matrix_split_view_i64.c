/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Fixed-width external-index import for split-triangular views.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <limits.h>
#include <stdint.h>
#include <stddef.h>

#include "lislib.h"
#include "lis_matrix_split_view.h"

static LIS_INT
lis_matrix_split_view_i64_fits_lis_int(int64_t value)
{
    if(sizeof(LIS_INT)>=sizeof(int64_t))
        return LIS_TRUE;

    return value>=(int64_t)INT_MIN && value<=(int64_t)INT_MAX;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_split_triangular_view_i64"
LIS_INT
lis_matrix_create_split_triangular_view_i64(
    LIS_Comm comm,
    int64_t n,
    int64_t nnz,
    const int64_t ptr[],
    const int64_t index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    int64_t origin,
    LIS_MATRIX *A,
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW *view)
{
    LIS_INT *ptr_lis=NULL;
    LIS_INT *index_lis=NULL;
    LIS_INT n_lis,nnz_lis,origin_lis;
    LIS_INT err;
    int64_t i;

    if(A==NULL || view==NULL)
        return LIS_ERR_ILL_ARG;

    *A=NULL;
    *view=NULL;

    if(n<=0 || nnz<0 || ptr==NULL || diag==NULL ||
       (nnz>0 && (index==NULL || lower_values==NULL)) ||
       (origin!=0 && origin!=1))
        return LIS_ERR_ILL_ARG;

    if(!lis_matrix_split_view_i64_fits_lis_int(n) ||
       !lis_matrix_split_view_i64_fits_lis_int(nnz) ||
       !lis_matrix_split_view_i64_fits_lis_int(origin))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "split-triangular 64-bit dimensions exceed LIS_INT range\n");
        return LIS_ERR_ILL_ARG;
    }

    if((uint64_t)n+UINT64_C(1)>(uint64_t)(SIZE_MAX/sizeof(LIS_INT)) ||
       (uint64_t)nnz>(uint64_t)(SIZE_MAX/sizeof(LIS_INT)))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "split-triangular 64-bit index arrays exceed address range\n");
        return LIS_ERR_ILL_ARG;
    }

    n_lis=(LIS_INT)n;
    nnz_lis=(LIS_INT)nnz;
    origin_lis=(LIS_INT)origin;

    ptr_lis=(LIS_INT *)lis_malloc(
        ((size_t)n_lis+1u)*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_view_i64::ptr");
    index_lis=nnz_lis>0 ? (LIS_INT *)lis_malloc(
        (size_t)nnz_lis*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_view_i64::index") : NULL;

    if(ptr_lis==NULL || (nnz_lis>0 && index_lis==NULL))
    {
        lis_free(index_lis);
        lis_free(ptr_lis);
        return LIS_ERR_OUT_OF_MEMORY;
    }

    for(i=0;i<=n;i++)
    {
        if(!lis_matrix_split_view_i64_fits_lis_int(ptr[i]))
        {
            err=LIS_ERR_ILL_ARG;
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "split-triangular 64-bit pointer exceeds LIS_INT range\n");
            goto cleanup;
        }
        ptr_lis[i]=(LIS_INT)ptr[i];
    }

    for(i=0;i<nnz;i++)
    {
        if(!lis_matrix_split_view_i64_fits_lis_int(index[i]))
        {
            err=LIS_ERR_ILL_ARG;
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "split-triangular 64-bit index exceeds LIS_INT range\n");
            goto cleanup;
        }
        index_lis[i]=(LIS_INT)index[i];
    }

    err=lis_matrix_create_split_triangular_view(
        comm,n_lis,nnz_lis,ptr_lis,index_lis,
        diag,lower_values,upper_values,origin_lis,A,view);

cleanup:
    lis_free(index_lis);
    lis_free(ptr_lis);
    return err;
}
