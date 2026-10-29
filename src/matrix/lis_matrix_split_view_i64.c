/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Fixed-width external-index import and shared topology support for
 * split-triangular views.
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
#include <string.h>

#include "lislib.h"
#include "lis_matrix_split_view.h"

struct LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY_STRUCT
{
    LIS_INT refs;
    LIS_INT n;
    LIS_INT nnz;
    LIS_INT *row_ptr;
    LIS_INT *row_col;
    LIS_INT *row_src;
};

struct LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW_STRUCT
{
    LIS_MATRIX A;
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology;
    LIS_INT symmetric;
    const LIS_SCALAR *diag;
    const LIS_SCALAR *lower;
    const LIS_SCALAR *upper;
};

static LIS_INT
lis_matrix_split_view_i64_fits_lis_int(int64_t value)
{
    if(sizeof(LIS_INT)>=sizeof(int64_t))
        return LIS_TRUE;

    return value>=(int64_t)INT_MIN && value<=(int64_t)INT_MAX;
}

static void
lis_matrix_split_topology_free(
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology)
{
    if(topology==NULL)
        return;

    lis_free(topology->row_src);
    lis_free(topology->row_col);
    lis_free(topology->row_ptr);
    memset(topology,0,sizeof(*topology));
    lis_free(topology);
}

static LIS_INT
lis_matrix_split_topology_retain(
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology)
{
    if(topology==NULL || topology->refs<=0)
        return LIS_ERR_ILL_ARG;

    topology->refs++;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_split_triangular_topology"
LIS_INT
lis_matrix_create_split_triangular_topology(
    LIS_INT n,
    LIS_INT nnz,
    const LIS_INT ptr[],
    const LIS_INT index[],
    LIS_INT origin,
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY *topology_out)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology=NULL;
    LIS_INT *next=NULL;
    LIS_INT i,j,k,p,start,end,entries;
    uintmax_t total;

    LIS_DEBUG_FUNC_IN;

    if(topology_out==NULL || n<=0 || nnz<0 || ptr==NULL ||
       (nnz>0 && index==NULL) || (origin!=0 && origin!=1))
        return LIS_ERR_ILL_ARG;

    *topology_out=NULL;

    total=(uintmax_t)n+2u*(uintmax_t)nnz;
    if(total>(sizeof(LIS_INT)==4 ? (uintmax_t)INT_MAX :
                                   (uintmax_t)LLONG_MAX))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "split-triangular topology exceeds LIS_INT range\n");
        return LIS_ERR_ILL_ARG;
    }
    entries=(LIS_INT)total;

    topology=(LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY)lis_malloc(
        sizeof(*topology),
        "lis_matrix_create_split_triangular_topology::topology");
    if(topology==NULL)
        return LIS_ERR_OUT_OF_MEMORY;
    memset(topology,0,sizeof(*topology));

    topology->n=n;
    topology->nnz=nnz;
    topology->row_ptr=(LIS_INT *)lis_malloc(
        ((size_t)n+1u)*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_topology::row_ptr");
    topology->row_col=entries>0 ? (LIS_INT *)lis_malloc(
        (size_t)entries*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_topology::row_col") : NULL;
    topology->row_src=entries>0 ? (LIS_INT *)lis_malloc(
        (size_t)entries*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_topology::row_src") : NULL;
    next=n>0 ? (LIS_INT *)lis_malloc(
        (size_t)n*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_topology::next") : NULL;

    if(topology->row_ptr==NULL ||
       (entries>0 && (topology->row_col==NULL || topology->row_src==NULL)) ||
       (n>0 && next==NULL))
    {
        lis_free(next);
        lis_matrix_split_topology_free(topology);
        return LIS_ERR_OUT_OF_MEMORY;
    }

    if(ptr[0]!=origin || ptr[n]-origin!=nnz)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "split-triangular pointer range is inconsistent with nnz\n");
        goto invalid;
    }

    topology->row_ptr[0]=0;
    for(j=0;j<n;j++)
    {
        start=ptr[j]-origin;
        end=ptr[j+1]-origin;
        if(start<0 || end<start || end>nnz)
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "split-triangular column pointer is invalid\n");
            goto invalid;
        }

        for(k=start;k<end;k++)
        {
            i=index[k]-origin;
            if(i<=j || i>=n ||
               (k>start && index[k]<=index[k-1]))
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "split-triangular indices must be sorted strict-lower columns\n");
                goto invalid;
            }
        }

        topology->row_ptr[j+1]=1+(end-start);
    }

    for(k=0;k<nnz;k++)
    {
        i=index[k]-origin;
        if(i<0 || i>=n)
            goto invalid;
        topology->row_ptr[i+1]++;
    }

    for(j=0;j<n;j++)
    {
        topology->row_ptr[j+1]+=topology->row_ptr[j];
        next[j]=topology->row_ptr[j];
    }

    for(j=0;j<n;j++)
    {
        start=ptr[j]-origin;
        end=ptr[j+1]-origin;
        for(k=start;k<end;k++)
        {
            i=index[k]-origin;
            p=next[i]++;
            topology->row_col[p]=j;
            topology->row_src[p]=k+1;
        }
    }

    for(j=0;j<n;j++)
    {
        p=next[j]++;
        topology->row_col[p]=j;
        topology->row_src[p]=-(j+1);

        start=ptr[j]-origin;
        end=ptr[j+1]-origin;
        for(k=start;k<end;k++)
        {
            i=index[k]-origin;
            p=next[j]++;
            topology->row_col[p]=i;
            topology->row_src[p]=nnz+k+1;
        }
    }

    for(j=0;j<n;j++)
    {
        if(next[j]!=topology->row_ptr[j+1])
            goto invalid;
    }

    lis_free(next);
    topology->refs=1;
    *topology_out=topology;
    return LIS_SUCCESS;

invalid:
    lis_free(next);
    lis_matrix_split_topology_free(topology);
    return LIS_ERR_ILL_ARG;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_split_triangular_topology_i64"
LIS_INT
lis_matrix_create_split_triangular_topology_i64(
    int64_t n,
    int64_t nnz,
    const int64_t ptr[],
    const int64_t index[],
    int64_t origin,
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY *topology)
{
    LIS_INT *ptr_lis=NULL;
    LIS_INT *index_lis=NULL;
    LIS_INT n_lis,nnz_lis,origin_lis;
    LIS_INT err;
    int64_t i;

    if(topology==NULL)
        return LIS_ERR_ILL_ARG;

    *topology=NULL;

    if(n<=0 || nnz<0 || ptr==NULL ||
       (nnz>0 && index==NULL) || (origin!=0 && origin!=1))
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
        "lis_matrix_create_split_triangular_topology_i64::ptr");
    index_lis=nnz_lis>0 ? (LIS_INT *)lis_malloc(
        (size_t)nnz_lis*sizeof(LIS_INT),
        "lis_matrix_create_split_triangular_topology_i64::index") : NULL;

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

    err=lis_matrix_create_split_triangular_topology(
        n_lis,nnz_lis,ptr_lis,index_lis,origin_lis,topology);

cleanup:
    lis_free(index_lis);
    lis_free(ptr_lis);
    return err;
}

void
lis_matrix_split_triangular_topology_destroy(
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology)
{
    if(topology==NULL)
        return;

    if(topology->refs<=0)
        return;

    topology->refs--;
    if(topology->refs==0)
        lis_matrix_split_topology_free(topology);
}

static LIS_INT
lis_matrix_split_shared_value(
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view,
    LIS_INT src,
    LIS_SCALAR *value)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology;
    uintmax_t offset;

    if(view==NULL || value==NULL || view->topology==NULL)
        return LIS_ERR_ILL_ARG;

    topology=view->topology;

    if(src<0)
    {
        LIS_INT i=-src-1;
        if(i<0 || i>=topology->n || view->diag==NULL)
            return LIS_ERR_ILL_ARG;
        *value=view->diag[i];
        return LIS_SUCCESS;
    }

    if(src<=0)
        return LIS_ERR_ILL_ARG;

    offset=(uintmax_t)(src-1);
    if(offset<(uintmax_t)topology->nnz)
    {
        if(view->lower==NULL)
            return LIS_ERR_ILL_ARG;
        *value=view->lower[(size_t)offset];
        return LIS_SUCCESS;
    }

    offset-=(uintmax_t)topology->nnz;
    if(offset>=(uintmax_t)topology->nnz)
        return LIS_ERR_ILL_ARG;

    if(view->symmetric)
    {
        if(view->lower==NULL)
            return LIS_ERR_ILL_ARG;
        *value=view->lower[(size_t)offset];
    }
    else
    {
        if(view->upper==NULL)
            return LIS_ERR_ILL_ARG;
        *value=view->upper[(size_t)offset];
    }

    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_split_shared_matvec(
    void *user_data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view;
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology;
    LIS_INT i,p;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW)user_data;
    if(view==NULL || x==NULL || y==NULL || view->topology==NULL ||
       view->diag==NULL ||
       (view->topology->nnz>0 && view->lower==NULL))
        return LIS_ERR_ILL_ARG;

    topology=view->topology;

#ifdef _OPENMP
#pragma omp parallel for private(p) schedule(static)
#endif
    for(i=0;i<topology->n;i++)
    {
        LIS_SCALAR sum=(LIS_SCALAR)0.0;

        for(p=topology->row_ptr[i];p<topology->row_ptr[i+1];p++)
        {
            LIS_INT src=topology->row_src[p];
            LIS_SCALAR value;

            if(src<0)
            {
                value=view->diag[-src-1];
            }
            else if(src<=topology->nnz)
            {
                value=view->lower[src-1];
            }
            else
            {
                LIS_INT k=src-topology->nnz-1;
                value=view->symmetric ? view->lower[k] : view->upper[k];
            }

            sum+=value*x[topology->row_col[p]];
        }

        y[i]=sum;
    }

    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_split_shared_matvech(
    void *user_data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view;
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology;
    LIS_INT i,p;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW)user_data;
    if(view==NULL || x==NULL || y==NULL || view->topology==NULL ||
       view->diag==NULL ||
       (view->topology->nnz>0 && view->lower==NULL))
        return LIS_ERR_ILL_ARG;

    topology=view->topology;

#ifdef _OPENMP
#pragma omp parallel for private(p) schedule(static)
#endif
    for(i=0;i<topology->n;i++)
    {
        LIS_SCALAR sum=(LIS_SCALAR)0.0;

        for(p=topology->row_ptr[i];p<topology->row_ptr[i+1];p++)
        {
            LIS_INT src=topology->row_src[p];
            LIS_SCALAR value;

            if(src<0)
            {
                value=conj(view->diag[-src-1]);
            }
            else if(src<=topology->nnz)
            {
                LIS_INT k=src-1;
                value=conj(view->symmetric ?
                    view->lower[k] : view->upper[k]);
            }
            else
            {
                LIS_INT k=src-topology->nnz-1;
                value=conj(view->lower[k]);
            }

            sum+=value*x[topology->row_col[p]];
        }

        y[i]=sum;
    }

    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_split_shared_get_row(
    void *user_data,
    LIS_INT local_row,
    LIS_INT *row_nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view;
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology;
    LIS_INT p0,p1,count,k,capacity,err;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW)user_data;
    if(view==NULL || view->topology==NULL || row_nnz==NULL)
        return LIS_ERR_ILL_ARG;

    topology=view->topology;
    if(local_row<0 || local_row>=topology->n)
        return LIS_ERR_ILL_ARG;

    p0=topology->row_ptr[local_row];
    p1=topology->row_ptr[local_row+1];
    if(p0<0 || p1<p0)
        return LIS_ERR_ILL_ARG;

    count=p1-p0;
    if(cols==NULL && vals==NULL)
    {
        *row_nnz=count;
        return LIS_SUCCESS;
    }

    if(cols==NULL || vals==NULL)
        return LIS_ERR_ILL_ARG;

    capacity=*row_nnz;
    if(capacity<count)
    {
        *row_nnz=count;
        return LIS_ERR_ILL_ARG;
    }

    for(k=0;k<count;k++)
    {
        LIS_INT pos=p0+k;
        cols[k]=topology->row_col[pos];
        err=lis_matrix_split_shared_value(
            view,topology->row_src[pos],&vals[k]);
        if(err)
            return err;
    }

    *row_nnz=count;
    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_split_shared_get_diagonal(
    void *user_data,
    LIS_SCALAR *diag)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view;
    LIS_INT i;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW)user_data;
    if(view==NULL || view->topology==NULL ||
       diag==NULL || view->diag==NULL)
        return LIS_ERR_ILL_ARG;

    for(i=0;i<view->topology->n;i++)
        diag[i]=view->diag[i];

    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_split_triangular_shared_view"
LIS_INT
lis_matrix_create_split_triangular_shared_view(
    LIS_Comm comm,
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology,
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[],
    LIS_MATRIX *A,
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW *view_out)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view=NULL;
    LIS_MATRIX M=NULL;
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    if(A==NULL || view_out==NULL || topology==NULL ||
       topology->refs<=0 || topology->n<=0 || topology->nnz<0 ||
       topology->row_ptr==NULL || diag==NULL ||
       (topology->nnz>0 && lower_values==NULL))
        return LIS_ERR_ILL_ARG;

    *A=NULL;
    *view_out=NULL;

    view=(LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW)lis_malloc(
        sizeof(*view),
        "lis_matrix_create_split_triangular_shared_view::view");
    if(view==NULL)
        return LIS_ERR_OUT_OF_MEMORY;
    memset(view,0,sizeof(*view));

    view->topology=topology;
    view->symmetric=upper_values==NULL;
    view->diag=diag;
    view->lower=lower_values;
    view->upper=upper_values;

    err=lis_matrix_create(comm,&M);
    if(err) goto fail;
    err=lis_matrix_set_size(M,0,topology->n);
    if(err) goto fail;
    if(M->nprocs!=1 || M->n!=M->gn)
    {
        err=LIS_ERR_NOT_IMPLEMENTED;
        goto fail;
    }

    err=lis_matrix_set_user(
        M,view,lis_matrix_split_shared_matvec,
        lis_matrix_split_shared_matvech);
    if(err) goto fail;
    err=lis_matrix_set_user_row(M,lis_matrix_split_shared_get_row);
    if(err) goto fail;
    err=lis_matrix_set_user_row_thread_safe(M,LIS_TRUE);
    if(err) goto fail;
    err=lis_matrix_set_user_diagonal(
        M,lis_matrix_split_shared_get_diagonal);
    if(err) goto fail;

    err=lis_matrix_split_topology_retain(topology);
    if(err) goto fail;

    view->A=M;
    *A=M;
    *view_out=view;
    return LIS_SUCCESS;

fail:
    if(M) lis_matrix_destroy(M);
    memset(view,0,sizeof(*view));
    lis_free(view);
    return err;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_triangular_shared_view_set_values"
LIS_INT
lis_matrix_split_triangular_shared_view_set_values(
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view,
    const LIS_SCALAR diag[],
    const LIS_SCALAR lower_values[],
    const LIS_SCALAR upper_values[])
{
    if(view==NULL || view->A==NULL || view->topology==NULL ||
       diag==NULL ||
       (view->topology->nnz>0 && lower_values==NULL))
        return LIS_ERR_ILL_ARG;

    if((upper_values==NULL)!=(view->symmetric!=0))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "shared split-triangular symmetry mode cannot change in a value update\n");
        return LIS_ERR_ILL_ARG;
    }

    view->diag=diag;
    view->lower=lower_values;
    view->upper=upper_values;

    return lis_matrix_mark_values_changed(view->A);
}

void
lis_matrix_split_triangular_shared_view_destroy(
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view)
{
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology;

    if(view==NULL)
        return;

    topology=view->topology;
    memset(view,0,sizeof(*view));
    lis_free(view);
    lis_matrix_split_triangular_topology_destroy(topology);
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
