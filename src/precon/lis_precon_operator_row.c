/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Row access for non-owning linear-combination matrix operators.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include "lislib.h"
#include "lis_precon.h"

typedef struct
{
    LIS_INT col;
    LIS_SCALAR value;
}
LIS_PRECON_OPERATOR_ENTRY;

static int
lis_precon_operator_entry_compare(const void *a,const void *b)
{
    const LIS_PRECON_OPERATOR_ENTRY *ea;
    const LIS_PRECON_OPERATOR_ENTRY *eb;

    ea=(const LIS_PRECON_OPERATOR_ENTRY *)a;
    eb=(const LIS_PRECON_OPERATOR_ENTRY *)b;
    if(ea->col<eb->col) return -1;
    if(ea->col>eb->col) return 1;
    return 0;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_operator_reserve"
static LIS_INT
lis_precon_operator_reserve(
    LIS_PRECON_ROW_SOURCE *source,
    LIS_INT needed)
{
    LIS_INT *new_index=NULL;
    LIS_SCALAR *new_value=NULL;

    if(source==NULL || needed<0)
        return LIS_ERR_ILL_ARG;

    if(needed<=source->capacity)
        return LIS_SUCCESS;

    if((uintmax_t)needed>
       (uintmax_t)SIZE_MAX/sizeof(LIS_INT) ||
       (uintmax_t)needed>
       (uintmax_t)SIZE_MAX/sizeof(LIS_SCALAR))
        return LIS_ERR_ILL_ARG;

    if(needed>0)
    {
        new_index=(LIS_INT *)lis_precon_row_source_alloc(source,
            needed*sizeof(LIS_INT),
            "lis_precon_operator_reserve::index");
        if(new_index==NULL)
        {
            LIS_SETERR_MEM(needed*sizeof(LIS_INT));
            return LIS_OUT_OF_MEMORY;
        }

        new_value=(LIS_SCALAR *)lis_precon_row_source_alloc(source,
            needed*sizeof(LIS_SCALAR),
            "lis_precon_operator_reserve::value");
        if(new_value==NULL)
        {
            lis_precon_row_source_release(source,new_index);
            LIS_SETERR_MEM(needed*sizeof(LIS_SCALAR));
            return LIS_OUT_OF_MEMORY;
        }
    }

    lis_precon_row_source_release(source,source->index_buffer);
    lis_precon_row_source_release(source,source->value_buffer);
    source->index_buffer=new_index;
    source->value_buffer=new_value;
    source->capacity=needed;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_operator_entries_reserve"
static LIS_INT
lis_precon_operator_entries_reserve(
    LIS_PRECON_ROW_SOURCE *source,
    LIS_INT needed)
{
    LIS_PRECON_OPERATOR_ENTRY *new_entries=NULL;

    if(source==NULL || needed<0)
        return LIS_ERR_ILL_ARG;

    if(needed<=source->operator_entries_capacity)
        return LIS_SUCCESS;

    if((uintmax_t)needed>
       (uintmax_t)SIZE_MAX/sizeof(LIS_PRECON_OPERATOR_ENTRY))
        return LIS_ERR_ILL_ARG;

    if(needed>0)
    {
        new_entries=(LIS_PRECON_OPERATOR_ENTRY *)lis_precon_row_source_alloc(source,
            needed*sizeof(*new_entries),
            "lis_precon_operator_entries_reserve::entries");
        if(new_entries==NULL)
        {
            LIS_SETERR_MEM(needed*sizeof(*new_entries));
            return LIS_OUT_OF_MEMORY;
        }
    }

    lis_precon_row_source_release(source,source->operator_entries);
    source->operator_entries=new_entries;
    source->operator_entries_capacity=needed;
    return LIS_SUCCESS;
}

LIS_INT
lis_precon_operator_row_source_init(
    LIS_PRECON_ROW_SOURCE *source)
{
    LIS_MATRIX A;
    LIS_INT err;

    if(source==NULL || source->A==NULL)
        return LIS_ERR_ILL_ARG;

    A=source->A;
    if(A->matrix_type!=LIS_MATRIX_OPERATOR || A->operator_A==NULL)
        return LIS_ERR_NOT_IMPLEMENTED;
    if(A->nprocs!=1 || A->n!=A->gn ||
       A->operator_A->n!=A->n || A->operator_A->gn!=A->gn)
        return LIS_ERR_NOT_IMPLEMENTED;

    source->operator_A=(LIS_PRECON_ROW_SOURCE *)lis_precon_row_source_alloc(source,
        sizeof(LIS_PRECON_ROW_SOURCE),
        "lis_precon_operator_row_source_init::A");
    if(source->operator_A==NULL)
        return LIS_OUT_OF_MEMORY;
    memset(source->operator_A,0,sizeof(*source->operator_A));
    err=lis_precon_row_source_init_mode(A->operator_A,source->operator_A,source->raw_scratch);
    if(err) goto fail;

    if(A->operator_scale!=NULL)
    {
        if(A->operator_B!=NULL || A->operator_scale->n!=A->n ||
           A->operator_scale->gn!=A->gn || A->operator_scale->value==NULL)
        {
            err=LIS_ERR_ILL_ARG;
            goto fail;
        }
        return LIS_SUCCESS;
    }

    if(A->operator_B==NULL || A->operator_B->n!=A->n ||
       A->operator_B->gn!=A->gn)
    {
        err=LIS_ERR_NOT_IMPLEMENTED;
        goto fail;
    }

    source->operator_B=(LIS_PRECON_ROW_SOURCE *)lis_precon_row_source_alloc(source,
        sizeof(LIS_PRECON_ROW_SOURCE),
        "lis_precon_operator_row_source_init::B");
    if(source->operator_B==NULL)
    {
        err=LIS_OUT_OF_MEMORY;
        goto fail;
    }
    memset(source->operator_B,0,sizeof(*source->operator_B));
    err=lis_precon_row_source_init_mode(A->operator_B,source->operator_B,source->raw_scratch);
    if(err) goto fail;

    return LIS_SUCCESS;

fail:
    lis_precon_operator_row_source_destroy(source);
    return err;
}

void
lis_precon_operator_row_source_destroy(
    LIS_PRECON_ROW_SOURCE *source)
{
    if(source==NULL)
        return;

    if(source->operator_A)
    {
        lis_precon_row_source_destroy(source->operator_A);
        lis_precon_row_source_release(source,source->operator_A);
        source->operator_A=NULL;
    }
    if(source->operator_B)
    {
        lis_precon_row_source_destroy(source->operator_B);
        lis_precon_row_source_release(source,source->operator_B);
        source->operator_B=NULL;
    }

    lis_precon_row_source_release(source,source->operator_entries);
    source->operator_entries=NULL;
    source->operator_entries_capacity=0;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_operator_row_source_get"
LIS_INT
lis_precon_operator_row_source_get(
    LIS_PRECON_ROW_SOURCE *source,
    LIS_INT local_row)
{
    LIS_PRECON_OPERATOR_ENTRY *entries;
    LIS_MATRIX A;
    LIS_INT na,nb,total,k,j,out,err;

    if(source==NULL || source->A==NULL || source->operator_A==NULL)
        return LIS_ERR_ILL_ARG;

    A=source->A;
    err=lis_precon_row_source_get(source->operator_A,local_row);
    if(err) return err;

    if(A->operator_scale!=NULL)
    {
        LIS_SCALAR left;

        na=source->operator_A->nnz;
        if(na<0 || local_row<0 || local_row>=A->n ||
           A->operator_scale->value==NULL)
            return LIS_ERR_ILL_ARG;

        err=lis_precon_operator_reserve(source,na);
        if(err) return err;

        left=A->operator_scale->value[local_row];
        for(k=0;k<na;k++)
        {
            LIS_INT col=source->operator_A->index[k];
            if(col<0 || col>=A->n)
                return LIS_ERR_ILL_ARG;
            source->index_buffer[k]=col;
            source->value_buffer[k]=left*source->operator_A->value[k]*
                                    A->operator_scale->value[col];
        }

        source->index=source->index_buffer;
        source->value=source->value_buffer;
        source->nnz=na;
        return LIS_SUCCESS;
    }

    if(source->operator_B==NULL)
        return LIS_ERR_ILL_ARG;
    err=lis_precon_row_source_get(source->operator_B,local_row);
    if(err) return err;

    na=source->operator_A->nnz;
    nb=source->operator_B->nnz;

    if(na<0 || nb<0)
        return LIS_ERR_ILL_ARG;

    {
        uintmax_t total_u=(uintmax_t)na+(uintmax_t)nb;
        uintmax_t lis_int_max=
            sizeof(LIS_INT)==4
            ? (uintmax_t)INT_MAX
            : (uintmax_t)LLONG_MAX;

        if(total_u>lis_int_max ||
           total_u>(uintmax_t)SIZE_MAX/sizeof(LIS_PRECON_OPERATOR_ENTRY))
            return LIS_ERR_ILL_ARG;

        total=(LIS_INT)total_u;
    }

    err=lis_precon_operator_entries_reserve(source,total);
    if(err)
        return err;
    entries=(LIS_PRECON_OPERATOR_ENTRY *)source->operator_entries;

    for(k=0;k<na;k++)
    {
        entries[k].col=source->operator_A->index[k];
        entries[k].value=A->operator_alpha*source->operator_A->value[k];
    }
    for(k=0;k<nb;k++)
    {
        entries[na+k].col=source->operator_B->index[k];
        entries[na+k].value=A->operator_beta*source->operator_B->value[k];
    }

    if(total>1)
        qsort(entries,(size_t)total,sizeof(*entries),
              lis_precon_operator_entry_compare);

    /* First pass determines the algebraically combined column count. */
    out=0;
    k=0;
    while(k<total)
    {
        LIS_INT col=entries[k].col;
        j=k+1;
        while(j<total && entries[j].col==col) j++;
        out++;
        k=j;
    }

    err=lis_precon_operator_reserve(source,out);
    if(err)
        return err;

    out=0;
    k=0;
    while(k<total)
    {
        LIS_INT col=entries[k].col;
        LIS_SCALAR value=entries[k].value;
        j=k+1;
        while(j<total && entries[j].col==col)
        {
            value+=entries[j].value;
            j++;
        }
        if(col<0 || col>=A->n)
            return LIS_ERR_ILL_ARG;
        /* Preserve the structural union even when values cancel to zero. */
        source->index_buffer[out]=col;
        source->value_buffer[out]=value;
        out++;
        k=j;
    }

    source->index=source->index_buffer;
    source->value=source->value_buffer;
    source->nnz=out;
    return LIS_SUCCESS;
}
