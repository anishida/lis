/* Conservative, generic AUTO matrix-pair eigensolver session.
 *
 * Both CSR caches belong to one session.  An old eigensolver (including
 * nested linear solvers and borrowed matrix references) is ALWAYS destroyed
 * before either cache is refreshed. No algorithm-specific work reuse is
 * attempted. The public materializer retains its own structural validation.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <string.h>
#include "lislib.h"
#include "lis_matrix_materialize.h"

struct LIS_GESOLVE_MATRIX_SESSION_STRUCT
{
    LIS_ESOLVER esolver;
    LIS_MATRIX_CSR_CACHE cache_a,cache_b;
    LIS_VECTOR single_vector;
    LIS_SCALAR single_evalue;
    LIS_INT computed_nev;
    LIS_INT committed;
};

/*
 * Transfer only user-supplied eigensolver CONFIGURATION to a new instance.
 * Never copy A/B/work/evalue/evector or nested solver ownership.  In
 * particular, option strings are deep-copied using the existing setters.
 */
static LIS_INT
lis_gesolve_session_clone_config(LIS_ESOLVER old,LIS_ESOLVER *out)
{
    LIS_ESOLVER next=NULL;
    LIS_INT err;

    if(out==NULL) return LIS_ERR_ILL_ARG;
    *out=NULL;
    err=lis_esolver_create(&next);
    if(err) return err;
    if(old!=NULL)
    {
        memcpy(next->options,old->options,sizeof(next->options));
        memcpy(next->params,old->params,sizeof(next->params));
        next->nev=old->requested_nev_active ?
            old->requested_nev : old->nev;
        next->requested_nev=old->requested_nev;
        next->requested_nev_active=old->requested_nev_active;
        next->ishift=old->ishift;
        next->transform=old->transform;
        next->feast_emin=old->feast_emin;
        next->feast_emax=old->feast_emax;
        next->feast_ratio=old->feast_ratio;
        next->feast_points=old->feast_points;
        next->feast_solver=old->feast_solver;
        next->feast_interval_mask=old->feast_interval_mask;
        next->block_size=old->block_size;
        next->refine_gap_fraction=old->refine_gap_fraction;
        next->subspace_retry_max=old->subspace_retry_max;
        next->subspace_retry_growth=old->subspace_retry_growth;
        next->subspace_retry_min=old->subspace_retry_min;
        if(old->linear_options!=NULL)
        {
            err=lis_esolver_set_linear_options(old->linear_options,next);
            if(err) goto fail;
        }
        if(old->refine_linear_options!=NULL)
        {
            err=lis_esolver_set_refine_linear_options(
                old->refine_linear_options,next);
            if(err) goto fail;
        }
    }
    *out=next;
    return LIS_SUCCESS;
fail:
    lis_esolver_destroy(next);
    return err;
}

/* Retire all references to A/B BEFORE either cache is changed or freed. */
static LIS_INT
lis_gesolve_session_retire(LIS_GESOLVE_MATRIX_SESSION s,LIS_INT drop_caches)
{
    LIS_ESOLVER next=NULL;
    LIS_INT err;

    err=lis_gesolve_session_clone_config(s->esolver,&next);
    if(err) return err;
    lis_esolver_destroy(s->esolver);
    s->esolver=next;
    if(s->single_vector)
    {
        lis_vector_destroy(s->single_vector);
        s->single_vector=NULL;
    }
    s->committed=LIS_FALSE;
    s->computed_nev=0;
    if(drop_caches)
    {
        lis_matrix_csr_cache_destroy(s->cache_a);
        lis_matrix_csr_cache_destroy(s->cache_b);
        s->cache_a=NULL;
        s->cache_b=NULL;
    }
    return LIS_SUCCESS;
}

static LIS_INT
lis_gesolve_session_prepare(LIS_MATRIX source,
    LIS_MATRIX_CSR_CACHE *cache,LIS_MATRIX *selected,LIS_INT *policy)
{
    if(source->matrix_type==LIS_MATRIX_USER ||
       source->matrix_type==LIS_MATRIX_OPERATOR)
    {
        LIS_INT action;
        LIS_INT err=lis_matrix_materialize_csr_cached(
            source,cache,selected,&action);
        if(err) return err;
        *policy=LIS_MATRIX_POLICY_CACHED_CSR;
    }
    else
    {
        *selected=source;
        *policy=LIS_MATRIX_POLICY_DIRECT;
    }
    return LIS_SUCCESS;
}

LIS_INT
lis_gesolve_matrix_session_create(LIS_GESOLVE_MATRIX_SESSION *out)
{
    LIS_GESOLVE_MATRIX_SESSION s;
    LIS_INT err;
    if(out==NULL) return LIS_ERR_ILL_ARG;
    *out=NULL;
    s=(LIS_GESOLVE_MATRIX_SESSION)lis_malloc(
        sizeof(*s),"lis_gesolve_matrix_session_create");
    if(s==NULL) return LIS_OUT_OF_MEMORY;
    memset(s,0,sizeof(*s));
    err=lis_esolver_create(&s->esolver);
    if(err)
    {
        lis_free(s);
        return err;
    }
    *out=s;
    return LIS_SUCCESS;
}

LIS_INT
lis_gesolve_matrix_session_get_esolver(
    LIS_GESOLVE_MATRIX_SESSION s,LIS_ESOLVER *esolver)
{
    if(s==NULL || esolver==NULL || s->esolver==NULL)
        return LIS_ERR_ILL_ARG;
    *esolver=s->esolver;
    return LIS_SUCCESS;
}

LIS_INT
lis_gesolve_matrix_session_solve(
    LIS_GESOLVE_MATRIX_SESSION s,
    LIS_MATRIX A,LIS_MATRIX B,LIS_VECTOR x,
    LIS_INT *policy_A,LIS_INT *policy_B)
{
    LIS_MATRIX a=NULL,b=NULL;
    LIS_VECTOR single=NULL;
    LIS_SCALAR lambda;
    LIS_INT pa=LIS_MATRIX_POLICY_DIRECT,pb=LIS_MATRIX_POLICY_DIRECT;
    LIS_INT err,count,i;
    LIS_ESOLVER e;

    if(s==NULL || s->esolver==NULL || A==NULL || B==NULL ||
       x==NULL || A->n<=0 || A->gn!=B->gn || A->n!=B->n ||
       A->n!=x->n || A->gn!=x->gn ||
       A->nprocs!=1 || B->nprocs!=1 || x->nprocs!=1)
        return LIS_ERR_ILL_ARG;

    /* Retire old eigensolver and all its nested references even if A/B
     * appear unchanged: algorithm-specific retention has not been audited.
     */
    err=lis_gesolve_session_retire(s,LIS_FALSE);
    if(err) return err;
    err=lis_gesolve_session_prepare(A,&s->cache_a,&a,&pa);
    if(err) goto fail;
    err=lis_gesolve_session_prepare(B,&s->cache_b,&b,&pb);
    if(err) goto fail;

    memset(&lambda,0,sizeof(lambda));
    err=lis_gesolve(a,b,x,&lambda,s->esolver);
    if(err) goto fail;
    e=s->esolver;
    if(e->retcode!=LIS_SUCCESS)
    {
        err=e->retcode;
        goto fail;
    }

    /*
     * Subspace methods own one stored vector for each available result.
     * Other generalized methods typically return one eigenpair in x.
     */
    count=0;
    if(e->evector!=NULL && e->evalue!=NULL)
    {
        count=e->evector_owned_count;
        if(count>e->options[LIS_EOPTIONS_SUBSPACE])
            count=e->options[LIS_EOPTIONS_SUBSPACE];
        for(i=0;i<count;i++)
            if(e->evector[i]==NULL)
            {
                count=i;
                break;
            }
    }
    if(e->options[LIS_EOPTIONS_ESOLVER]==LIS_ESOLVER_FEAST &&
       e->nev>=0 && count>e->nev)
        count=e->nev;
    if(count==0 && e->options[LIS_EOPTIONS_ESOLVER]!=LIS_ESOLVER_FEAST)
    {
        err=lis_vector_duplicate(x,&single);
        if(err) goto fail;
        err=lis_vector_copy(x,single);
        if(err) goto fail;
        s->single_vector=single;
        single=NULL;
        memcpy(&s->single_evalue,&lambda,sizeof(lambda));
        count=1;
    }
    s->computed_nev=count;
    s->committed=LIS_TRUE;
    if(policy_A) *policy_A=pa;
    if(policy_B) *policy_B=pb;
    return LIS_SUCCESS;

fail:
    if(single) lis_vector_destroy(single);
    /* Invalidate BOTH caches after any partial preparation or solve
     * failure. The old esolver must be retired before freeing them. */
    {
        LIS_INT reset_err=lis_gesolve_session_retire(s,LIS_TRUE);
        if(reset_err)
        {
            lis_esolver_destroy(s->esolver);
            s->esolver=NULL;
            if(s->single_vector)
            {
                lis_vector_destroy(s->single_vector);
                s->single_vector=NULL;
            }
            lis_matrix_csr_cache_destroy(s->cache_a);
            lis_matrix_csr_cache_destroy(s->cache_b);
            s->cache_a=NULL;
            s->cache_b=NULL;
            s->committed=LIS_FALSE;
            s->computed_nev=0;
        }
    }
    return err;
}

LIS_INT
lis_gesolve_matrix_session_get_nev(
    LIS_GESOLVE_MATRIX_SESSION s,LIS_INT *computed_nev)
{
    if(s==NULL || computed_nev==NULL || !s->committed)
        return LIS_ERR_ILL_ARG;
    *computed_nev=s->computed_nev;
    return LIS_SUCCESS;
}

LIS_INT
lis_gesolve_matrix_session_get_eigenpair(
    LIS_GESOLVE_MATRIX_SESSION s,LIS_INT mode,
    LIS_SCALAR *eigenvalue,LIS_VECTOR eigenvector_out)
{
    LIS_VECTOR src;
    const LIS_SCALAR *value;
    LIS_INT err;
    if(s==NULL || !s->committed || mode<0 ||
       mode>=s->computed_nev || eigenvalue==NULL ||
       eigenvector_out==NULL)
        return LIS_ERR_ILL_ARG;
    if(s->single_vector!=NULL)
    {
        src=s->single_vector;
        value=&s->single_evalue;
    }
    else
    {
        if(s->esolver==NULL || s->esolver->evector==NULL ||
           s->esolver->evalue==NULL)
            return LIS_ERR_ILL_ARG;
        src=s->esolver->evector[mode];
        value=&s->esolver->evalue[mode];
    }
    if(src==NULL || src->value==NULL ||
       eigenvector_out->value==NULL ||
       src->precision!=LIS_PRECISION_DEFAULT ||
       eigenvector_out->precision!=LIS_PRECISION_DEFAULT ||
       src->n!=eigenvector_out->n ||
       src->gn!=eigenvector_out->gn ||
       src->is!=eigenvector_out->is ||
       src->ie!=eigenvector_out->ie)
        return LIS_ERR_ILL_ARG;

    err=lis_vector_copy(src,eigenvector_out);
    if(err) return err;
    memcpy(eigenvalue,value,sizeof(*eigenvalue));
    return LIS_SUCCESS;
}

LIS_INT
lis_gesolve_matrix_session_reset(LIS_GESOLVE_MATRIX_SESSION s)
{
    if(s==NULL || s->esolver==NULL) return LIS_ERR_ILL_ARG;
    return lis_gesolve_session_retire(s,LIS_TRUE);
}

void
lis_gesolve_matrix_session_destroy(LIS_GESOLVE_MATRIX_SESSION s)
{
    if(s==NULL) return;
    lis_esolver_destroy(s->esolver);
    if(s->single_vector) lis_vector_destroy(s->single_vector);
    lis_matrix_csr_cache_destroy(s->cache_a);
    lis_matrix_csr_cache_destroy(s->cache_b);
    lis_free(s);
}
