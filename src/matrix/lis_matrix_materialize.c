/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Generic cached CSR materialization for row-backed matrices.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stdint.h>
#include <string.h>

#include "lislib.h"
#include "lis_matrix_materialize.h"

#undef __FUNC__
#define __FUNC__ "lis_matrix_materialize_csr_cached"

struct LIS_MATRIX_CSR_CACHE_STRUCT
{
    LIS_MATRIX source;
    LIS_MATRIX csr;
    uint64_t structure_stamp;
    uint64_t value_stamp;
    LIS_INT valid;
};

#define LIS_MATERIALIZE_HASH_OFFSET UINT64_C(1469598103934665603)
#define LIS_MATERIALIZE_HASH_PRIME  UINT64_C(1099511628211)

static uint64_t
lis_matrix_materialize_hash_word(uint64_t hash,uint64_t value)
{
    hash^=value;
    hash*=LIS_MATERIALIZE_HASH_PRIME;
    return hash;
}

static uint64_t
lis_matrix_materialize_hash_bytes(
    uint64_t hash,
    const void *data,
    size_t size)
{
    const unsigned char *p=(const unsigned char *)data;
    size_t i;

    for(i=0;i<size;i++)
    {
        hash^=(uint64_t)p[i];
        hash*=LIS_MATERIALIZE_HASH_PRIME;
    }
    return hash;
}

static LIS_INT
lis_matrix_materialize_source_stamp(
    LIS_MATRIX A,
    uint64_t *structure_stamp,
    uint64_t *value_stamp,
    LIS_INT depth)
{
    LIS_UNSIGNED_INT sg=0,vg=0;
    uint64_t hs=LIS_MATERIALIZE_HASH_OFFSET;
    uint64_t hv=LIS_MATERIALIZE_HASH_OFFSET;
    LIS_INT err;

    if(A==NULL || structure_stamp==NULL || value_stamp==NULL || depth>64)
        return LIS_ERR_ILL_ARG;

    err=lis_matrix_get_generation(A,&sg,&vg);
    if(err) return err;

    hs=lis_matrix_materialize_hash_word(hs,(uint64_t)(uintptr_t)A);
    hs=lis_matrix_materialize_hash_word(hs,(uint64_t)A->matrix_type);
    hs=lis_matrix_materialize_hash_word(hs,(uint64_t)sg);

    hv=lis_matrix_materialize_hash_word(hv,(uint64_t)(uintptr_t)A);
    hv=lis_matrix_materialize_hash_word(hv,(uint64_t)A->matrix_type);
    hv=lis_matrix_materialize_hash_word(hv,(uint64_t)vg);

    if(A->matrix_type==LIS_MATRIX_OPERATOR)
    {
        uint64_t as,avs,bs,bvs;

        if(A->operator_A==NULL)
            return LIS_ERR_ILL_ARG;

        err=lis_matrix_materialize_source_stamp(
            A->operator_A,&as,&avs,depth+1);
        if(err) return err;
        hs=lis_matrix_materialize_hash_word(hs,as);
        hv=lis_matrix_materialize_hash_word(hv,avs);

        if(A->operator_scale!=NULL)
        {
            if(A->operator_B!=NULL || A->operator_scale->n!=A->n ||
               A->operator_scale->gn!=A->gn)
                return LIS_ERR_ILL_ARG;

            /*
             * A scale vector can be changed in place without replacing
             * its address. Fingerprint its values to invalidate cached CSR
             * even when the application has not marked the operator.
             */
            if(A->n<0 || (size_t)A->n>SIZE_MAX/sizeof(LIS_SCALAR) ||
               (A->n>0 && A->operator_scale->value==NULL))
                return LIS_ERR_ILL_ARG;
            hv=lis_matrix_materialize_hash_word(
                hv,(uint64_t)(uintptr_t)A->operator_scale);
            hv=lis_matrix_materialize_hash_bytes(
                hv,A->operator_scale->value,
                (size_t)A->n*sizeof(LIS_SCALAR));
        }
        else
        {
            if(A->operator_B==NULL)
                return LIS_ERR_ILL_ARG;

            err=lis_matrix_materialize_source_stamp(
                A->operator_B,&bs,&bvs,depth+1);
            if(err) return err;

            hs=lis_matrix_materialize_hash_word(hs,bs);
            hv=lis_matrix_materialize_hash_word(hv,bvs);
            hv=lis_matrix_materialize_hash_bytes(
                hv,&A->operator_alpha,sizeof(A->operator_alpha));
            hv=lis_matrix_materialize_hash_bytes(
                hv,&A->operator_beta,sizeof(A->operator_beta));
        }
    }

    *structure_stamp=hs;
    *value_stamp=hv;
    return LIS_SUCCESS;
}

static LIS_INT
lis_matrix_materialize_build_csr(
    LIS_MATRIX source,
    LIS_MATRIX *csr_out)
{
    LIS_PRECON_ROW_SOURCE row;
    LIS_MATRIX B=NULL;
    LIS_INT *ptr=NULL,*index=NULL;
    LIS_SCALAR *value=NULL;
    LIS_INT nnz=0,cursor=0,i,k,err;
    LIS_INT row_ready=LIS_FALSE;
    LIS_INT storage_allocated=LIS_FALSE;
    LIS_INT storage_attached=LIS_FALSE;

    if(source==NULL || csr_out==NULL)
        return LIS_ERR_ILL_ARG;
    *csr_out=NULL;

    if(source->nprocs!=1 || source->n!=source->gn)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "cached CSR materialization is serial-only\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }

    err=lis_precon_row_source_count_nnz(source,&nnz);
    if(err) return err;

    err=lis_matrix_duplicate(source,&B);
    if(err) goto fail;

    err=lis_matrix_malloc_csr(source->n,nnz,&ptr,&index,&value);
    if(err) goto fail;
    storage_allocated=LIS_TRUE;

    err=lis_precon_row_source_init(source,&row);
    if(err) goto fail;
    row_ready=LIS_TRUE;

    ptr[0]=0;
    for(i=0;i<source->n;i++)
    {
        err=lis_precon_row_source_get(&row,i);
        if(err) goto fail;
        if(row.nnz<0 || row.nnz>nnz-cursor)
        {
            err=LIS_ERR_ILL_ARG;
            goto fail;
        }

        for(k=0;k<row.nnz;k++)
        {
            index[cursor]=row.index[k];
            value[cursor]=row.value[k];
            cursor++;
        }
        ptr[i+1]=cursor;
    }

    if(cursor!=nnz)
    {
        err=LIS_ERR_ILL_ARG;
        goto fail;
    }

    lis_precon_row_source_destroy(&row);
    row_ready=LIS_FALSE;

    err=lis_matrix_set_csr(nnz,ptr,index,value,B);
    if(err) goto fail;
    storage_attached=LIS_TRUE;

    err=lis_matrix_assemble(B);
    if(err) goto fail;

    *csr_out=B;
    return LIS_SUCCESS;

fail:
    if(row_ready)
        lis_precon_row_source_destroy(&row);
    if(B)
        lis_matrix_destroy(B);
    if(storage_allocated && !storage_attached)
    {
        lis_free(ptr);
        lis_free(index);
        lis_free(value);
    }
    return err;
}

static LIS_INT
lis_matrix_materialize_pattern_matches(
    LIS_MATRIX source,
    LIS_MATRIX csr,
    LIS_INT *matches)
{
    LIS_PRECON_ROW_SOURCE row;
    LIS_INT i,k,p0,p1,err;

    if(source==NULL || csr==NULL || matches==NULL)
        return LIS_ERR_ILL_ARG;

    *matches=LIS_FALSE;
    if(csr->matrix_type!=LIS_MATRIX_CSR ||
       source->n!=csr->n || source->gn!=csr->gn)
        return LIS_SUCCESS;

    err=lis_precon_row_source_init(source,&row);
    if(err) return err;

    for(i=0;i<source->n;i++)
    {
        err=lis_precon_row_source_get(&row,i);
        if(err) goto done;

        p0=csr->ptr[i];
        p1=csr->ptr[i+1];
        if(row.nnz!=p1-p0)
            goto mismatch;

        for(k=0;k<row.nnz;k++)
            if(row.index[k]!=csr->index[p0+k])
                goto mismatch;
    }

    *matches=LIS_TRUE;
    err=LIS_SUCCESS;
    goto done;

mismatch:
    *matches=LIS_FALSE;
    err=LIS_SUCCESS;

done:
    lis_precon_row_source_destroy(&row);
    return err;
}

static LIS_INT
lis_matrix_materialize_refresh_values(
    LIS_MATRIX source,
    LIS_MATRIX csr)
{
    LIS_PRECON_ROW_SOURCE row;
    LIS_INT i,k,p0,err;

    err=lis_precon_row_source_init(source,&row);
    if(err) return err;

    for(i=0;i<source->n;i++)
    {
        err=lis_precon_row_source_get(&row,i);
        if(err) goto done;

        p0=csr->ptr[i];
        for(k=0;k<row.nnz;k++)
            csr->value[p0+k]=row.value[k];
    }

    /*
     * A solver/preconditioner may have split this cache-owned CSR matrix
     * (for example, SSOR uses the D/L/U representation). Updating the
     * canonical CSR values does not update those derived arrays. Discard
     * them so subsequent matvecs and solves use the refreshed coefficients.
     * Keep ptr/index/value and the CSR matrix object itself allocated.
     */
    if(csr->is_splited)
    {
        err=lis_matrix_DLU_destroy(csr);
        if(err) goto done;
        csr->is_splited=LIS_FALSE;
    }
    if(csr->WD!=NULL)
    {
        lis_matrix_diag_destroy(csr->WD);
        csr->WD=NULL;
    }
    csr->use_wd=0;
    csr->is_scaled=LIS_FALSE;

    err=lis_matrix_mark_values_changed(csr);

done:
    lis_precon_row_source_destroy(&row);
    return err;
}

LIS_INT
lis_matrix_materialize_csr_cached(
    LIS_MATRIX source,
    LIS_MATRIX_CSR_CACHE *cache_io,
    LIS_MATRIX *csr,
    LIS_INT *action)
{
    LIS_MATRIX_CSR_CACHE cache;
    LIS_MATRIX fresh=NULL;
    uint64_t structure_stamp,value_stamp;
    LIS_INT matches,err;

    LIS_DEBUG_FUNC_IN;

    if(source==NULL || cache_io==NULL || csr==NULL)
        return LIS_ERR_ILL_ARG;

    *csr=NULL;
    if(action) *action=LIS_MATRIX_MATERIALIZE_REBUILD;

    err=lis_matrix_materialize_source_stamp(
        source,&structure_stamp,&value_stamp,0);
    if(err) return err;

    cache=*cache_io;
    if(cache==NULL)
    {
        err=lis_matrix_materialize_build_csr(source,&fresh);
        if(err) return err;

        cache=(LIS_MATRIX_CSR_CACHE)lis_malloc(
            sizeof(*cache),"lis_matrix_materialize_csr_cached::cache");
        if(cache==NULL)
        {
            lis_matrix_destroy(fresh);
            return LIS_OUT_OF_MEMORY;
        }
        memset(cache,0,sizeof(*cache));
        cache->source=source;
        cache->csr=fresh;
        cache->structure_stamp=structure_stamp;
        cache->value_stamp=value_stamp;
        cache->valid=LIS_TRUE;
        *cache_io=cache;
        *csr=fresh;
        return LIS_SUCCESS;
    }

    if(cache->valid && cache->source==source &&
       cache->structure_stamp==structure_stamp &&
       cache->value_stamp==value_stamp)
    {
        *csr=cache->csr;
        if(action) *action=LIS_MATRIX_MATERIALIZE_REUSE;
        return LIS_SUCCESS;
    }

    if(cache->valid && cache->source==source &&
       cache->structure_stamp==structure_stamp)
    {
        err=lis_matrix_materialize_pattern_matches(
            source,cache->csr,&matches);
        if(err) return err;

        if(matches)
        {
            /*
             * The pattern was validated in a separate pass. If a row callback
             * fails during the value pass, invalidate the cache so the next
             * call rebuilds instead of exposing a partially refreshed CSR.
             */
            cache->valid=LIS_FALSE;
            err=lis_matrix_materialize_refresh_values(source,cache->csr);
            if(err) return err;

            cache->value_stamp=value_stamp;
            cache->valid=LIS_TRUE;
            *csr=cache->csr;
            if(action) *action=LIS_MATRIX_MATERIALIZE_VALUES;
            return LIS_SUCCESS;
        }
    }

    err=lis_matrix_materialize_build_csr(source,&fresh);
    if(err) return err;

    if(cache->csr)
        lis_matrix_destroy(cache->csr);
    cache->source=source;
    cache->csr=fresh;
    cache->structure_stamp=structure_stamp;
    cache->value_stamp=value_stamp;
    cache->valid=LIS_TRUE;

    *csr=fresh;
    if(action) *action=LIS_MATRIX_MATERIALIZE_REBUILD;
    return LIS_SUCCESS;
}

void
lis_matrix_csr_cache_destroy(
    LIS_MATRIX_CSR_CACHE cache)
{
    if(cache==NULL)
        return;

    if(cache->csr)
        lis_matrix_destroy(cache->csr);
    memset(cache,0,sizeof(*cache));
    lis_free(cache);
}
