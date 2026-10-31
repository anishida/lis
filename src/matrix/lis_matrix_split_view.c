/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Non-owning split-triangular external matrix view.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "lislib.h"
#include "lis_matrix_split_view.h"

struct LIS_MATRIX_SPLIT_TRIANGULAR_VIEW_STRUCT
{
    LIS_MATRIX A;
    LIS_INT n;
    LIS_INT nnz;
    LIS_INT symmetric;
    const LIS_SCALAR *diag;
    const LIS_SCALAR *lower;
    const LIS_SCALAR *upper;
    LIS_INT *row_ptr;
    LIS_INT *row_col;
    LIS_INT *row_src;
};

static LIS_INT
lis_matrix_split_view_value(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    LIS_INT src,
    LIS_SCALAR *value)
{
    uintmax_t offset;

    if(view==NULL || value==NULL)
        return LIS_ERR_ILL_ARG;

    if(src<0)
    {
        LIS_INT i = -src-1;
        if(i<0 || i>=view->n || view->diag==NULL)
            return LIS_ERR_ILL_ARG;
        *value = view->diag[i];
        return LIS_SUCCESS;
    }

    if(src<=0)
        return LIS_ERR_ILL_ARG;

    offset = (uintmax_t)(src-1);
    if(offset<(uintmax_t)view->nnz)
    {
        if(view->lower==NULL)
            return LIS_ERR_ILL_ARG;
        *value = view->lower[(size_t)offset];
        return LIS_SUCCESS;
    }

    offset -= (uintmax_t)view->nnz;
    if(offset>=(uintmax_t)view->nnz)
        return LIS_ERR_ILL_ARG;

    if(view->symmetric)
    {
        if(view->lower==NULL)
            return LIS_ERR_ILL_ARG;
        *value = view->lower[(size_t)offset];
    }
    else
    {
        if(view->upper==NULL)
            return LIS_ERR_ILL_ARG;
        *value = view->upper[(size_t)offset];
    }

    return LIS_SUCCESS;
}

/*
 * Native row kernels avoid the generic USER-row scratch path.  The expanded
 * row metadata is immutable between structure rebinds, so each output row is
 * independent and can be evaluated safely in parallel.
 */
static LIS_INT
lis_matrix_split_view_matvec(
    void *user_data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view;
    LIS_INT i,p;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_VIEW)user_data;
    if(view==NULL || x==NULL || y==NULL || view->diag==NULL ||
       (view->nnz>0 && view->lower==NULL))
        return LIS_ERR_ILL_ARG;

#ifdef _OPENMP
#pragma omp parallel for private(p) schedule(static)
#endif
    for(i=0;i<view->n;i++)
    {
        LIS_SCALAR sum=(LIS_SCALAR)0.0;

        for(p=view->row_ptr[i];p<view->row_ptr[i+1];p++)
        {
            LIS_INT src=view->row_src[p];
            LIS_SCALAR value;

            if(src<0)
            {
                value=view->diag[-src-1];
            }
            else if(src<=view->nnz)
            {
                value=view->lower[src-1];
            }
            else
            {
                LIS_INT k=src-view->nnz-1;
                value=view->symmetric ? view->lower[k] : view->upper[k];
            }

            sum+=value*x[view->row_col[p]];
        }

        y[i]=sum;
    }

    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_split_view_matvech(
    void *user_data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view;
    LIS_INT i,p;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_VIEW)user_data;
    if(view==NULL || x==NULL || y==NULL || view->diag==NULL ||
       (view->nnz>0 && view->lower==NULL))
        return LIS_ERR_ILL_ARG;

#ifdef _OPENMP
#pragma omp parallel for private(p) schedule(static)
#endif
    for(i=0;i<view->n;i++)
    {
        LIS_SCALAR sum=(LIS_SCALAR)0.0;

        for(p=view->row_ptr[i];p<view->row_ptr[i+1];p++)
        {
            LIS_INT src=view->row_src[p];
            LIS_SCALAR value;

            if(src<0)
            {
                value=conj(view->diag[-src-1]);
            }
            else if(src<=view->nnz)
            {
                LIS_INT k=src-1;
                value=conj(view->symmetric ? view->lower[k] : view->upper[k]);
            }
            else
            {
                LIS_INT k=src-view->nnz-1;
                value=conj(view->lower[k]);
            }

            sum+=value*x[view->row_col[p]];
        }

        y[i]=sum;
    }

    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_split_view_get_row(
    void *user_data,
    LIS_INT local_row,
    LIS_INT *row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view;
    LIS_INT p0,p1,count,k,capacity,err;

    view = (LIS_MATRIX_SPLIT_TRIANGULAR_VIEW)user_data;
    if(view==NULL || row_nnz==NULL ||
       local_row<0 || local_row>=view->n)
        return LIS_ERR_ILL_ARG;

    p0 = view->row_ptr[local_row];
    p1 = view->row_ptr[local_row+1];
    if(p0<0 || p1<p0)
        return LIS_ERR_ILL_ARG;

    count = p1-p0;
    if(cols==NULL && vals==NULL)
    {
        *row_nnz = count;
        return LIS_SUCCESS;
    }

    if(cols==NULL || vals==NULL)
        return LIS_ERR_ILL_ARG;

    capacity = *row_nnz;
    if(capacity<count)
    {
        *row_nnz = count;
        return LIS_ERR_ILL_ARG;
    }

    for(k=0;k<count;k++)
    {
        LIS_INT p = p0+k;
        cols[k] = view->row_col[p];
        err = lis_matrix_split_view_value(view,view->row_src[p],&vals[k]);
        if(err)
            return err;
    }

    *row_nnz = count;
    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_split_view_get_diagonal(
    void *user_data,
    LIS_SCALAR *diag)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view;
    LIS_INT i;

    view = (LIS_MATRIX_SPLIT_TRIANGULAR_VIEW)user_data;
    if(view==NULL || diag==NULL || view->diag==NULL)
        return LIS_ERR_ILL_ARG;

    for(i=0;i<view->n;i++)
        diag[i] = view->diag[i];

    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_view_build_rows"
static LIS_INT
lis_matrix_split_view_build_rows(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    const LIS_INT ptr[],
    const LIS_INT index[],
    LIS_INT origin)
{
    LIS_INT *next=NULL;
    LIS_INT i,j,k,p,start,end;
    LIS_INT entries;
    uintmax_t total;

    if(view==NULL || ptr==NULL ||
       (view->nnz>0 && index==NULL) ||
       (origin!=0 && origin!=1))
        return LIS_ERR_ILL_ARG;

    total = (uintmax_t)view->n + 2u*(uintmax_t)view->nnz;
    if(total>(sizeof(LIS_INT)==4 ? (uintmax_t)INT_MAX : (uintmax_t)LLONG_MAX))
        return LIS_ERR_ILL_ARG;
    entries = (LIS_INT)total;

    view->row_ptr = (LIS_INT *)lis_malloc(
        (view->n+1)*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_view::row_ptr");
    view->row_col = entries ? (LIS_INT *)lis_malloc(
        entries*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_view::row_col") : NULL;
    view->row_src = entries ? (LIS_INT *)lis_malloc(
        entries*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_view::row_src") : NULL;
    next = view->n ? (LIS_INT *)lis_malloc(
        view->n*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_view::next") : NULL;

    if(view->row_ptr==NULL ||
       (entries>0 && (view->row_col==NULL || view->row_src==NULL)) ||
       (view->n>0 && next==NULL))
    {
        lis_free(next);
        return LIS_ERR_OUT_OF_MEMORY;
    }

    if(ptr[0]!=origin || ptr[view->n]-origin!=view->nnz)
    {
        lis_free(next);
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "split-triangular pointer range is inconsistent with nnz\n");
        return LIS_ERR_ILL_ARG;
    }

    view->row_ptr[0]=0;
    for(j=0;j<view->n;j++)
    {
        start=ptr[j]-origin;
        end=ptr[j+1]-origin;
        if(start<0 || end<start || end>view->nnz)
        {
            lis_free(next);
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "split-triangular column pointer is invalid\n");
            return LIS_ERR_ILL_ARG;
        }

        /*
         * The public view contract is one strictly lower triangle stored by
         * columns. Enforcing row>column and strictly increasing row indices
         * makes the expanded rows deterministic and sorted around the
         * diagonal, which is also the useful order for row-based ILU.
         */
        for(k=start;k<end;k++)
        {
            i=index[k]-origin;
            if(i<=j || i>=view->n ||
               (k>start && index[k]<=index[k-1]))
            {
                lis_free(next);
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "split-triangular indices must be sorted strict-lower columns\n");
                return LIS_ERR_ILL_ARG;
            }
        }

        view->row_ptr[j+1]=1+(end-start);
    }

    for(k=0;k<view->nnz;k++)
    {
        i=index[k]-origin;
        if(i<0 || i>=view->n)
        {
            lis_free(next);
            return LIS_ERR_ILL_ARG;
        }
        view->row_ptr[i+1]++;
    }

    for(j=0;j<view->n;j++)
    {
        view->row_ptr[j+1]+=view->row_ptr[j];
        next[j]=view->row_ptr[j];
    }

    for(j=0;j<view->n;j++)
    {
        start=ptr[j]-origin;
        end=ptr[j+1]-origin;
        for(k=start;k<end;k++)
        {
            i=index[k]-origin;
            p=next[i]++;
            view->row_col[p]=j;
            view->row_src[p]=k+1;
        }
    }

    for(j=0;j<view->n;j++)
    {
        p=next[j]++;
        view->row_col[p]=j;
        view->row_src[p]=-(j+1);

        start=ptr[j]-origin;
        end=ptr[j+1]-origin;
        for(k=start;k<end;k++)
        {
            i=index[k]-origin;
            p=next[j]++;
            view->row_col[p]=i;
            view->row_src[p]=view->nnz+k+1;
        }
    }

    for(j=0;j<view->n;j++)
    {
        if(next[j]!=view->row_ptr[j+1])
        {
            lis_free(next);
            return LIS_ERR_ILL_ARG;
        }
    }

    lis_free(next);
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_split_triangular_view"
LIS_INT
lis_matrix_create_split_triangular_view(
    LIS_Comm comm,
    LIS_INT n,
    LIS_INT nnz,
    const LIS_INT ptr[],
    const LIS_INT index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    LIS_INT origin,
    LIS_MATRIX *A,
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW *view_out)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view=NULL;
    LIS_MATRIX M=NULL;
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    if(A==NULL || view_out==NULL || n<=0 || nnz<0 || ptr==NULL ||
       diag==NULL || (nnz>0 && (index==NULL || lower_values==NULL)))
        return LIS_ERR_ILL_ARG;

    *A=NULL;
    *view_out=NULL;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_VIEW)lis_malloc(
        sizeof(*view),"lis_matrix_create_split_triangular_view::view");
    if(view==NULL)
        return LIS_ERR_OUT_OF_MEMORY;
    memset(view,0,sizeof(*view));

    view->n=n;
    view->nnz=nnz;
    view->symmetric=upper_values==NULL;
    view->diag=diag;
    view->lower=lower_values;
    view->upper=upper_values;

    err=lis_matrix_split_view_build_rows(view,ptr,index,origin);
    if(err) goto fail;

    err=lis_matrix_create(comm,&M);
    if(err) goto fail;
    err=lis_matrix_set_size(M,0,n);
    if(err) goto fail;
    if(M->nprocs!=1 || M->n!=M->gn)
    {
        err=LIS_ERR_NOT_IMPLEMENTED;
        goto fail;
    }

    err=lis_matrix_set_user(
        M,view,lis_matrix_split_view_matvec,lis_matrix_split_view_matvech);
    if(err) goto fail;
    err=lis_matrix_set_user_row(M,lis_matrix_split_view_get_row);
    if(err) goto fail;
    err=lis_matrix_set_user_row_thread_safe(M,LIS_TRUE);
    if(err) goto fail;
    err=lis_matrix_set_user_diagonal(M,lis_matrix_split_view_get_diagonal);
    if(err) goto fail;

    view->A=M;
    *A=M;
    *view_out=view;
    return LIS_SUCCESS;

fail:
    if(M) lis_matrix_destroy(M);
    lis_matrix_split_triangular_view_destroy(view);
    return err;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_triangular_view_set_values"
LIS_INT
lis_matrix_split_triangular_view_set_values(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[])
{
    if(view==NULL || view->A==NULL || diag==NULL ||
       (view->nnz>0 && lower_values==NULL))
        return LIS_ERR_ILL_ARG;

    if((upper_values==NULL)!=(view->symmetric!=0))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "split-triangular symmetry mode cannot change in a value update\n");
        return LIS_ERR_ILL_ARG;
    }

    view->diag=diag;
    view->lower=lower_values;
    view->upper=upper_values;

    return lis_matrix_mark_values_changed(view->A);
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_triangular_view_rebind"
LIS_INT
lis_matrix_split_triangular_view_rebind(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    LIS_INT nnz,
    const LIS_INT ptr[],
    const LIS_INT index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    LIS_INT origin)
{
    struct LIS_MATRIX_SPLIT_TRIANGULAR_VIEW_STRUCT candidate;
    LIS_INT err;

    if(view==NULL || view->A==NULL || nnz<0 || ptr==NULL || diag==NULL ||
       (nnz>0 && (index==NULL || lower_values==NULL)))
        return LIS_ERR_ILL_ARG;

    memset(&candidate,0,sizeof(candidate));
    candidate.n=view->n;
    candidate.nnz=nnz;
    candidate.symmetric=upper_values==NULL;
    candidate.diag=diag;
    candidate.lower=lower_values;
    candidate.upper=upper_values;

    err=lis_matrix_split_view_build_rows(&candidate,ptr,index,origin);
    if(err) goto fail;

    /*
     * Publish the structural generation before swapping storage.  After the
     * candidate metadata has been built successfully the remaining swap is
     * allocation-free and cannot fail.
     */
    err=lis_matrix_mark_structure_changed(view->A);
    if(err) goto fail;

    lis_free(view->row_src);
    lis_free(view->row_col);
    lis_free(view->row_ptr);

    view->nnz=candidate.nnz;
    view->symmetric=candidate.symmetric;
    view->diag=candidate.diag;
    view->lower=candidate.lower;
    view->upper=candidate.upper;
    view->row_ptr=candidate.row_ptr;
    view->row_col=candidate.row_col;
    view->row_src=candidate.row_src;

    return LIS_SUCCESS;

fail:
    lis_free(candidate.row_src);
    lis_free(candidate.row_col);
    lis_free(candidate.row_ptr);
    return err;
}

static LIS_INT
lis_matrix_split_view_i64_fits_lis_int(int64_t value)
{
    if(sizeof(LIS_INT)>=sizeof(int64_t))
        return LIS_TRUE;

    return value>=(int64_t)INT_MIN && value<=(int64_t)INT_MAX;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_triangular_view_rebind_i64"
LIS_INT
lis_matrix_split_triangular_view_rebind_i64(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view,
    int64_t nnz,
    const int64_t ptr[],
    const int64_t index[],
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    int64_t origin)
{
    LIS_INT *ptr_lis=NULL;
    LIS_INT *index_lis=NULL;
    LIS_INT nnz_lis,origin_lis,err;
    LIS_INT i;

    if(view==NULL || view->A==NULL || nnz<0 || ptr==NULL || diag==NULL ||
       (nnz>0 && (index==NULL || lower_values==NULL)) ||
       (origin!=0 && origin!=1))
        return LIS_ERR_ILL_ARG;

    if(!lis_matrix_split_view_i64_fits_lis_int(nnz) ||
       (uint64_t)view->n+UINT64_C(1)>(uint64_t)(SIZE_MAX/sizeof(LIS_INT)) ||
       (uint64_t)nnz>(uint64_t)(SIZE_MAX/sizeof(LIS_INT)))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "split-triangular 64-bit rebind exceeds LIS_INT or address range\n");
        return LIS_ERR_ILL_ARG;
    }

    nnz_lis=(LIS_INT)nnz;
    origin_lis=(LIS_INT)origin;
    ptr_lis=(LIS_INT *)lis_malloc(
        ((size_t)view->n+1u)*sizeof(LIS_INT),
        "lis_matrix_split_triangular_view_rebind_i64::ptr");
    index_lis=nnz_lis>0 ? (LIS_INT *)lis_malloc(
        (size_t)nnz_lis*sizeof(LIS_INT),
        "lis_matrix_split_triangular_view_rebind_i64::index") : NULL;
    if(ptr_lis==NULL || (nnz_lis>0 && index_lis==NULL))
    {
        lis_free(index_lis);
        lis_free(ptr_lis);
        return LIS_ERR_OUT_OF_MEMORY;
    }

    for(i=0;i<=view->n;i++)
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

    for(i=0;i<nnz_lis;i++)
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

    err=lis_matrix_split_triangular_view_rebind(
        view,nnz_lis,ptr_lis,index_lis,
        diag,lower_values,upper_values,origin_lis);

cleanup:
    lis_free(index_lis);
    lis_free(ptr_lis);
    return err;
}

void
lis_matrix_split_triangular_view_destroy(
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view)
{
    if(view==NULL)
        return;

    lis_free(view->row_src);
    lis_free(view->row_col);
    lis_free(view->row_ptr);
    memset(view,0,sizeof(*view));
    lis_free(view);
}
