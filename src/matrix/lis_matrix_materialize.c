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
#include <limits.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "lislib.h"
#include "lis_matrix_materialize.h"
#include "lis_solver_capabilities.h"

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


/*
 * Concurrent row reads are only legal when every USER leaf opts in.
 * Each OpenMP worker owns a separate LIS_PRECON_ROW_SOURCE and its scratch.
 * The fallback remains serial for unsafe callbacks and small problems.
 */
/*
 * USER matrices in MPI builds need not have a communication table.
 * Allocate an independent cache-owned CSR target with the same distribution.
 */
#ifdef USE_MPI
static LIS_INT
lis_matrix_materialize_create_target(LIS_MATRIX source,LIS_MATRIX *out)
{
    LIS_INT err;
    *out=NULL;
    err=lis_matrix_create(source->comm,out);
    if(err) return err;
    err=lis_matrix_set_size(*out,source->n,source->gn);
    if(err)
    {
        lis_matrix_destroy(*out);
        *out=NULL;
    }
    return err;
}
#endif

#ifdef _OPENMP
static LIS_INT
lis_matrix_materialize_rows_parallel_safe(LIS_MATRIX A,LIS_INT depth)
{
    if(A==NULL || depth>64 || A->nprocs!=1 || A->n!=A->gn)
        return LIS_FALSE;

    switch(A->matrix_type)
    {
    case LIS_MATRIX_CSR:
        return LIS_TRUE;
    case LIS_MATRIX_USER:
        return A->user_get_row!=NULL &&
               A->user_get_row_thread_safe==LIS_TRUE;
    case LIS_MATRIX_OPERATOR:
        if(A->operator_A==NULL ||
           !lis_matrix_materialize_rows_parallel_safe(A->operator_A,depth+1))
            return LIS_FALSE;
        if(A->operator_scale!=NULL)
            return A->operator_B==NULL &&
                   A->operator_scale->value!=NULL;
        return lis_matrix_materialize_rows_parallel_safe(
            A->operator_B,depth+1);
    default:
        return LIS_FALSE;
    }
}

static LIS_INT
lis_matrix_materialize_use_parallel_rows(LIS_MATRIX A)
{
    return A!=NULL && A->n>=1024 && omp_get_max_threads()>1 &&
           !omp_in_parallel() &&
           lis_matrix_materialize_rows_parallel_safe(A,0);
}

/*
 * mode=0: count entries per row; mode=1: fill CSR.
 * On fill, row bounds must match the counted structure.
 */
static LIS_INT
lis_matrix_materialize_parallel_rows(
    LIS_MATRIX A,
    LIS_INT mode,
    LIS_INT *counts,
    const LIS_INT *ptr,
    LIS_INT *index,
    LIS_SCALAR *value)
{
    LIS_INT error=LIS_SUCCESS;

    #pragma omp parallel
    {
        LIS_PRECON_ROW_SOURCE row;
        LIS_INT i,local_error,ready;
        local_error=lis_precon_row_source_init_mode(A,&row,LIS_TRUE);
        ready=(local_error==LIS_SUCCESS);

        #pragma omp for schedule(static)
        for(i=0;i<A->n;i++)
        {
            LIS_INT k;
            if(local_error!=LIS_SUCCESS)
                continue;

            local_error=lis_precon_row_source_get(&row,i);
            if(local_error!=LIS_SUCCESS)
                continue;
            if(row.nnz<0)
            {
                local_error=LIS_ERR_ILL_ARG;
                continue;
            }
            if(mode==0)
            {
                counts[i]=row.nnz;
                continue;
            }
            if(ptr==NULL || row.nnz!=ptr[i+1]-ptr[i])
            {
                local_error=LIS_ERR_ILL_ARG;
                continue;
            }
            if(row.nnz>0 &&
               (row.index==NULL || row.value==NULL || value==NULL || index==NULL))
            {
                local_error=LIS_ERR_ILL_ARG;
                continue;
            }
            for(k=0;k<row.nnz;k++)
            {
                index[ptr[i]+k]=row.index[k];
                value[ptr[i]+k]=row.value[k];
            }
        }

        if(ready)
            lis_precon_row_source_destroy(&row);

        if(local_error!=LIS_SUCCESS)
        {
            #pragma omp critical(lis_materialize_parallel_error)
            {
                if(error==LIS_SUCCESS)
                    error=local_error;
            }
        }
    }
    return error;
}

static LIS_INT
lis_matrix_materialize_build_csr_parallel(
    LIS_MATRIX source,
    LIS_MATRIX *csr_out)
{
    LIS_MATRIX B=NULL;
    LIS_INT *counts=NULL,*ptr=NULL,*index=NULL;
    LIS_SCALAR *value=NULL;
    LIS_INT i,nnz=0,err;
    LIS_INT attached=LIS_FALSE;
    uintmax_t max_nnz=sizeof(LIS_INT)==4 ?
        (uintmax_t)INT_MAX : (uintmax_t)LLONG_MAX;

    if((uintmax_t)source->n>(uintmax_t)SIZE_MAX/sizeof(LIS_INT))
        return LIS_ERR_ILL_ARG;

    counts=(LIS_INT *)lis_malloc(
        (size_t)source->n*sizeof(LIS_INT),
        "lis_matrix_materialize_build_csr_parallel::counts");
    if(counts==NULL)
        return LIS_OUT_OF_MEMORY;

    err=lis_matrix_materialize_parallel_rows(
        source,0,counts,NULL,NULL,NULL);
    if(err) goto fail;

    for(i=0;i<source->n;i++)
    {
        if(counts[i]<0 || (uintmax_t)counts[i]>max_nnz-(uintmax_t)nnz)
        {
            err=LIS_ERR_ILL_ARG;
            goto fail;
        }
        nnz+=counts[i];
    }

#ifdef USE_MPI
    err=lis_matrix_materialize_create_target(source,&B);
#else
    /* Retain existing duplicate semantics in non-MPI builds. */
    err=lis_matrix_duplicate(source,&B);
#endif
    if(err) goto fail;

    /* The cache owns these CSR buffers even if the source does not own its data. */
    B->is_destroy=LIS_TRUE;
    err=lis_matrix_malloc_csr(source->n,nnz,&ptr,&index,&value);
    if(err) goto fail;

    ptr[0]=0;
    for(i=0;i<source->n;i++)
        ptr[i+1]=ptr[i]+counts[i];

    err=lis_matrix_materialize_parallel_rows(
        source,1,NULL,ptr,index,value);
    if(err) goto fail;

    err=lis_matrix_set_csr(nnz,ptr,index,value,B);
    if(err) goto fail;
    attached=LIS_TRUE;
    err=lis_matrix_assemble(B);
    if(err) goto fail;

    *csr_out=B;
    lis_free(counts);
    return LIS_SUCCESS;

fail:
    if(B) lis_matrix_destroy(B);
    if(ptr!=NULL && !attached)
    {
        lis_free(ptr);
        lis_free(index);
        lis_free(value);
    }
    lis_free(counts);
    return err;
}
#endif

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

#ifdef _OPENMP
    if(lis_matrix_materialize_use_parallel_rows(source))
        return lis_matrix_materialize_build_csr_parallel(source,csr_out);
#endif

    err=lis_precon_row_source_count_nnz(source,&nnz);
    if(err) return err;

#ifdef USE_MPI
    err=lis_matrix_materialize_create_target(source,&B);
#else
    /* Retain existing duplicate semantics in non-MPI builds. */
    err=lis_matrix_duplicate(source,&B);
#endif
    if(err) goto fail;

    /* The cache owns these CSR buffers even if the source does not own its data. */
    B->is_destroy=LIS_TRUE;

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

    /*
     * Refresh stays serial: the preceding full structural verification is
     * serial too, and benchmarks did not show a consistent OpenMP benefit.
     */
    err=lis_precon_row_source_init(source,&row);
    if(err) return err;
    for(i=0;i<source->n;i++)
    {
        err=lis_precon_row_source_get(&row,i);
        if(err) break;
        p0=csr->ptr[i];
        for(k=0;k<row.nnz;k++)
            csr->value[p0+k]=row.value[k];
    }
    lis_precon_row_source_destroy(&row);
    if(err) return err;

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


/*
 * Check whether every leaf of a USER/operator expression exposes row data
 * usable by the existing row-backed preconditioner path. This is a
 * capability query only; it does not fetch or copy matrix values.
 */
static LIS_INT
lis_matrix_policy_has_row_access(
    LIS_MATRIX A,
    LIS_INT depth)
{
    if(A==NULL || depth>64)
        return LIS_FALSE;

    if(A->nprocs!=1 || A->n!=A->gn)
        return LIS_FALSE;

    switch(A->matrix_type)
    {
    case LIS_MATRIX_CSR:
        return A->ptr!=NULL ? LIS_TRUE : LIS_FALSE;

    case LIS_MATRIX_USER:
        return A->user_get_row!=NULL ? LIS_TRUE : LIS_FALSE;

    case LIS_MATRIX_OPERATOR:
        if(A->operator_A==NULL)
            return LIS_FALSE;
        if(A->operator_scale!=NULL)
            return lis_matrix_policy_has_row_access(
                A->operator_A,depth+1);
        if(A->operator_B==NULL)
            return LIS_FALSE;
        return lis_matrix_policy_has_row_access(
                   A->operator_A,depth+1)
            && lis_matrix_policy_has_row_access(
                   A->operator_B,depth+1);

    default:
        return LIS_FALSE;
    }
}

static LIS_INT
lis_matrix_policy_has_matvec(
    LIS_MATRIX A,
    LIS_INT depth)
{
    if(A==NULL || depth>64)
        return LIS_FALSE;

    switch(A->matrix_type)
    {
    case LIS_MATRIX_CSR:
        return LIS_TRUE;

    case LIS_MATRIX_USER:
        return A->user_matvec!=NULL ? LIS_TRUE : LIS_FALSE;

    case LIS_MATRIX_OPERATOR:
        if(A->operator_A==NULL)
            return LIS_FALSE;
        if(A->operator_scale!=NULL)
            return lis_matrix_policy_has_matvec(
                A->operator_A,depth+1);
        if(A->operator_B==NULL)
            return LIS_FALSE;
        return lis_matrix_policy_has_matvec(
                   A->operator_A,depth+1)
            && lis_matrix_policy_has_matvec(
                   A->operator_B,depth+1);

    default:
        return LIS_TRUE;
    }
}

static LIS_INT
lis_matrix_policy_direct_supported(
    LIS_MATRIX A,
    LIS_INT precon_type)
{
    LIS_INT has_matvec;

    if(A==NULL)
        return LIS_FALSE;

    if(A->matrix_type==LIS_MATRIX_USER ||
       A->matrix_type==LIS_MATRIX_OPERATOR)
        has_matvec=lis_matrix_policy_has_matvec(A,0);
    else
        return LIS_TRUE;

    if(!has_matvec)
        return LIS_FALSE;

    if(precon_type==LIS_PRECON_TYPE_NONE)
        return LIS_TRUE;

    if(precon_type==LIS_PRECON_TYPE_JACOBI)
    {
        return A->matrix_type==LIS_MATRIX_USER
            && A->user_get_diagonal!=NULL
            ? LIS_TRUE : LIS_FALSE;
    }

    if(precon_type==LIS_PRECON_TYPE_ILU ||
       precon_type==LIS_PRECON_TYPE_ILUC ||
       precon_type==LIS_PRECON_TYPE_ILUT)
    {
        return lis_matrix_policy_has_row_access(A,0);
    }

    return LIS_FALSE;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_prepare_solver_matrix"
LIS_INT
lis_matrix_prepare_solver_matrix(
    LIS_MATRIX source,
    LIS_SOLVER solver,
    LIS_MATRIX_CSR_CACHE *cache,
    LIS_MATRIX *selected,
    LIS_INT *policy,
    LIS_INT *action)
{
    LIS_INT precon_type,nested,err;

    LIS_DEBUG_FUNC_IN;

    if(source==NULL || solver==NULL || selected==NULL || policy==NULL)
        return LIS_ERR_ILL_ARG;

    *selected=NULL;
    *policy=LIS_MATRIX_POLICY_DIRECT;
    if(action) *action=LIS_MATRIX_POLICY_ACTION_DIRECT;

    /*
     * This policy is deliberately limited to external USER matrices and
     * linear-combination operators. All other formats retain LIS's existing
     * solver/preconditioner conversion behavior.
     */
    if(source->matrix_type!=LIS_MATRIX_USER &&
       source->matrix_type!=LIS_MATRIX_OPERATOR)
    {
        *selected=source;
        return LIS_SUCCESS;
    }

    err=lis_solver_get_effective_preconditioner(
        solver,&precon_type,&nested);
    if(err) return err;

    if(lis_matrix_policy_direct_supported(source,precon_type))
    {
        *selected=source;
        return LIS_SUCCESS;
    }

    if(cache==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "automatic matrix policy requires a CSR cache for fallback\n");
        return LIS_ERR_ILL_ARG;
    }

    err=lis_matrix_materialize_csr_cached(
        source,cache,selected,action);
    if(err) return err;

    *policy=LIS_MATRIX_POLICY_CACHED_CSR;
    return LIS_SUCCESS;
}


/*
 * Own the solver and the optional materialized CSR as one lifetime unit.
 * selected is borrowed: it aliases source on DIRECT and cache->csr on CSR.
 */
struct LIS_MATRIX_AUTO_SESSION_STRUCT
{
    LIS_SOLVER solver;
    LIS_MATRIX_CSR_CACHE cache;
    LIS_MATRIX selected;
};

#undef __FUNC__
#define __FUNC__ "lis_matrix_auto_session_create"
LIS_INT
lis_matrix_auto_session_create(LIS_MATRIX_AUTO_SESSION *session_out)
{
    LIS_MATRIX_AUTO_SESSION session;
    LIS_INT err;

    if(session_out==NULL)
        return LIS_ERR_ILL_ARG;
    *session_out=NULL;

    session=(LIS_MATRIX_AUTO_SESSION)lis_malloc(
        sizeof(*session),"lis_matrix_auto_session_create::session");
    if(session==NULL)
        return LIS_OUT_OF_MEMORY;
    memset(session,0,sizeof(*session));

    err=lis_solver_create(&session->solver);
    if(err) goto fail;
    err=lis_solver_set_sequence_policy(
        session->solver,LIS_SEQUENCE_POLICY_AUTO);
    if(err) goto fail;

    *session_out=session;
    return LIS_SUCCESS;

fail:
    if(session->solver) lis_solver_destroy(session->solver);
    lis_free(session);
    return err;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_auto_session_get_solver"
LIS_INT
lis_matrix_auto_session_get_solver(
    LIS_MATRIX_AUTO_SESSION session,
    LIS_SOLVER *solver)
{
    if(session==NULL || solver==NULL || session->solver==NULL)
        return LIS_ERR_ILL_ARG;
    *solver=session->solver;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_auto_session_solve"
LIS_INT
lis_matrix_auto_session_solve(
    LIS_MATRIX_AUTO_SESSION session,
    LIS_MATRIX source,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_INT *policy,
    LIS_INT *action)
{
    LIS_MATRIX selected=NULL;
    LIS_INT chosen_policy=LIS_MATRIX_POLICY_DIRECT;
    LIS_INT chosen_action=LIS_MATRIX_POLICY_ACTION_DIRECT;
    LIS_INT err;

    if(session==NULL || session->solver==NULL ||
       source==NULL || b==NULL || x==NULL)
        return LIS_ERR_ILL_ARG;

    if(policy) *policy=chosen_policy;
    if(action) *action=chosen_action;

    /*
     * A cached CSR may be updated in place or DESTROYED by the following
     * policy call. An existing sequence preconditioner may still reference
     * that CSR (or its split D/L/U buffers), so release the solver state
     * before any cache mutation, including a failed prior refresh.
     *
     * Inspect the same effective nested-source stamps as the CSR cache;
     * a value-only update must reset the solver too, since refresh may
     * invalidate arrays referenced by the old preconditioner.
     */
    if(session->cache!=NULL &&
       session->selected==session->cache->csr)
    {
        uint64_t structure_stamp,value_stamp;
        err=lis_matrix_materialize_source_stamp(
            source,&structure_stamp,&value_stamp,0);
        if(err) return err;

        if(!session->cache->valid ||
           session->cache->source!=source ||
           session->cache->structure_stamp!=structure_stamp ||
           session->cache->value_stamp!=value_stamp)
        {
            err=lis_solver_reset_sequence(session->solver);
            if(err) return err;
            session->selected=NULL;
        }
    }

    err=lis_matrix_prepare_solver_matrix(
        source,session->solver,&session->cache,&selected,
        &chosen_policy,&chosen_action);
    if(err)
    {
        /* A failure may follow partial value refresh: force a fresh CSR. */
        if(session->cache) session->cache->valid=LIS_FALSE;
        return err;
    }

    /*
     * DIRECT <-> CSR selection or source identity changed. Retain
     * persistent state only when the actual solver matrix is unchanged.
     */
    if(session->selected!=NULL && session->selected!=selected)
    {
        err=lis_solver_reset_sequence(session->solver);
        if(err) return err;
    }

    session->selected=selected;
    if(policy) *policy=chosen_policy;
    if(action) *action=chosen_action;
    return lis_solve_sequence(selected,b,x,session->solver);
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_auto_session_reset"
LIS_INT
lis_matrix_auto_session_reset(LIS_MATRIX_AUTO_SESSION session)
{
    LIS_INT err;

    if(session==NULL || session->solver==NULL)
        return LIS_ERR_ILL_ARG;

    /* The preconditioner may reference the cached CSR: destroy it first. */
    err=lis_solver_reset_sequence(session->solver);
    if(err) return err;

    lis_matrix_csr_cache_destroy(session->cache);
    session->cache=NULL;
    session->selected=NULL;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_auto_session_destroy"
void
lis_matrix_auto_session_destroy(LIS_MATRIX_AUTO_SESSION session)
{
    if(session==NULL)
        return;

    /* Always destroy solver/preconditioner state before cache-owned CSR. */
    if(session->solver)
        lis_solver_destroy(session->solver);
    lis_matrix_csr_cache_destroy(session->cache);
    memset(session,0,sizeof(*session));
    lis_free(session);
}
