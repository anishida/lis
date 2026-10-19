/*
 * Flexible GCRO-DR.
 *
 * Stage 1B provides solver registration, persistent private recycle
 * state, configuration APIs and lifecycle handling.
 *
 * The actual flexible recycling iteration is added separately.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <string.h>

#include "lislib.h"


#define FGCRODR_DEFAULT_RECYCLE_DIM      ((LIS_INT)8)
#define FGCRODR_DEFAULT_KEEP_QUALITY     ((LIS_REAL)0.90)
#define FGCRODR_DEFAULT_REJECT_QUALITY   ((LIS_REAL)0.10)
#define FGCRODR_BUILD_TINY               ((LIS_REAL)1.0e-14)


struct LIS_FGCRODR_STATE_STRUCT
{
    /*
     * Flexible recycle relation used by the numerical core:
     *
     *     A Z = C
     *     C^T C = I
     *
     * W:
     *   recycle/search-space representation used by the update.
     *
     * Z:
     *   vectors after application of the (possibly changing)
     *   preconditioner.
     *
     * C:
     *   orthonormalized images A*Z.
     *
     * A later preconditioner must never be applied to W in order to
     * reconstruct Z.  Z is persistent state in its own right.
     */
    LIS_VECTOR *w;
    LIS_VECTOR *z;
    LIS_VECTOR *c;

    LIS_INT dim;
    LIS_INT alloc;
    LIS_INT valid;

    /* Owner of the exact recycle relation A Z = C. */
    LIS_MATRIX matrix;
    LIS_UNSIGNED_INT structure_generation;
    LIS_UNSIGNED_INT value_generation;

    LIS_INT target_dim;

    LIS_REAL keep_threshold;
    LIS_REAL reject_threshold;

    LIS_INT update;

    LIS_REAL quality;
    LIS_REAL quality_min;

    LIS_INT last_action;

    /*
     * Cross-system action selected before entering the flexible core.
     * This is transient for one solve, while last_action remains
     * available as the public diagnostic result.
     */
    LIS_INT cross_action;

    LIS_UNSIGNED_INT builds;
    LIS_UNSIGNED_INT updates;
    LIS_UNSIGNED_INT uses;
    LIS_UNSIGNED_INT rejects;

    LIS_UNSIGNED_INT cross_uses;

    LIS_UNSIGNED_INT policy_keeps;
    LIS_UNSIGNED_INT policy_transports;
    LIS_UNSIGNED_INT policy_rejects;
};


static void
lis_fgcrodr_state_defaults(
    LIS_FGCRODR_STATE state)
{
    memset(state,0,sizeof(*state));

    state->target_dim =
        FGCRODR_DEFAULT_RECYCLE_DIM;

    state->keep_threshold =
        FGCRODR_DEFAULT_KEEP_QUALITY;

    state->reject_threshold =
        FGCRODR_DEFAULT_REJECT_QUALITY;

    state->update = LIS_TRUE;

    state->last_action =
        LIS_FGCRODR_RECYCLE_NONE;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_state_ensure"
static LIS_INT
lis_fgcrodr_state_ensure(
    LIS_SOLVER solver)
{
    LIS_FGCRODR_STATE state;

    if( solver==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "solver is NULL\n");

        return LIS_ERR_ILL_ARG;
    }

    if( solver->fgcrodr_state!=NULL )
        return LIS_SUCCESS;

    state = (LIS_FGCRODR_STATE)lis_malloc(
        sizeof(struct LIS_FGCRODR_STATE_STRUCT),
        "lis_fgcrodr_state_ensure::state");

    if( state==NULL )
    {
        LIS_SETERR_MEM(
            sizeof(struct LIS_FGCRODR_STATE_STRUCT));

        return LIS_OUT_OF_MEMORY;
    }

    lis_fgcrodr_state_defaults(state);

    solver->fgcrodr_state = state;

    return LIS_SUCCESS;
}


static void
lis_fgcrodr_vectors_destroy(
    LIS_FGCRODR_STATE state)
{
    LIS_INT i;

    if( state==NULL )
        return;

    if( state->w!=NULL )
    {
        for(i=0;i<state->alloc;i++)
        {
            if( state->w[i] )
                lis_vector_destroy(state->w[i]);
        }

        lis_free(state->w);
    }

    if( state->z!=NULL )
    {
        for(i=0;i<state->alloc;i++)
        {
            if( state->z[i] )
                lis_vector_destroy(state->z[i]);
        }

        lis_free(state->z);
    }

    if( state->c!=NULL )
    {
        for(i=0;i<state->alloc;i++)
        {
            if( state->c[i] )
                lis_vector_destroy(state->c[i]);
        }

        lis_free(state->c);
    }

    state->w = NULL;
    state->z = NULL;
    state->c = NULL;

    state->dim = 0;
    state->alloc = 0;
    state->valid = LIS_FALSE;

    /*
     * No A-dependent metadata may survive recycle destruction.
     * Configuration and lifetime counters intentionally remain.
     */
    state->matrix = NULL;
    state->structure_generation = 0;
    state->value_generation = 0;

    state->quality = 0.0;
    state->quality_min = 0.0;

    state->last_action =
        LIS_FGCRODR_RECYCLE_NONE;

    state->cross_action =
        LIS_FGCRODR_RECYCLE_NONE;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_clear_recycle"
LIS_INT
lis_fgcrodr_clear_recycle(
    LIS_SOLVER solver)
{
    LIS_DEBUG_FUNC_IN;

    /*
     * A sequence reset is allowed before FGCRO-DR has ever been
     * selected, so an absent private state is not an error.
     */
    if( solver!=NULL &&
        solver->fgcrodr_state!=NULL )
    {
        lis_fgcrodr_vectors_destroy(
            solver->fgcrodr_state);
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_state_destroy"
LIS_INT
lis_fgcrodr_state_destroy(
    LIS_SOLVER solver)
{
    LIS_DEBUG_FUNC_IN;

    if( solver!=NULL &&
        solver->fgcrodr_state!=NULL )
    {
        lis_fgcrodr_vectors_destroy(
            solver->fgcrodr_state);

        lis_free(
            solver->fgcrodr_state);

        solver->fgcrodr_state = NULL;
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_set_fgcrodr_recycle"
LIS_INT
lis_solver_set_fgcrodr_recycle(
    LIS_SOLVER solver,
    LIS_INT dim,
    LIS_REAL keep_quality,
    LIS_REAL reject_quality)
{
    LIS_INT err;
    LIS_FGCRODR_STATE state;

    LIS_DEBUG_FUNC_IN;

    if( solver==NULL || dim<0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid FGCRO-DR recycle configuration\n");

        return LIS_ERR_ILL_ARG;
    }

    if( keep_quality!=keep_quality ||
        reject_quality!=reject_quality ||
        keep_quality<0.0 ||
        keep_quality>1.0 ||
        reject_quality<0.0 ||
        reject_quality>1.0 ||
        reject_quality>keep_quality )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "FGCRO-DR quality thresholds must satisfy "
            "0 <= reject <= keep <= 1\n");

        return LIS_ERR_ILL_ARG;
    }

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    state = solver->fgcrodr_state;

    if( state->target_dim!=dim )
    {
        lis_fgcrodr_vectors_destroy(state);
    }

    state->target_dim = dim;
    state->keep_threshold = keep_quality;
    state->reject_threshold = reject_quality;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_get_fgcrodr_recycle_config"
LIS_INT
lis_solver_get_fgcrodr_recycle_config(
    LIS_SOLVER solver,
    LIS_INT *dim,
    LIS_REAL *keep_quality,
    LIS_REAL *reject_quality,
    LIS_INT *update)
{
    LIS_INT err;
    LIS_FGCRODR_STATE state;

    LIS_DEBUG_FUNC_IN;

    if( solver==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "solver is NULL\n");

        return LIS_ERR_ILL_ARG;
    }

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    state = solver->fgcrodr_state;

    if( dim )
        *dim = state->target_dim;

    if( keep_quality )
        *keep_quality = state->keep_threshold;

    if( reject_quality )
        *reject_quality = state->reject_threshold;

    if( update )
        *update = state->update;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_get_fgcrodr_recycle"
LIS_INT
lis_solver_get_fgcrodr_recycle(
    LIS_SOLVER solver,
    LIS_INT *dim,
    LIS_INT *action,
    LIS_REAL *quality,
    LIS_REAL *quality_min)
{
    LIS_INT err;
    LIS_FGCRODR_STATE state;

    LIS_DEBUG_FUNC_IN;

    if( solver==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "solver is NULL\n");

        return LIS_ERR_ILL_ARG;
    }

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    state = solver->fgcrodr_state;

    if( dim )
        *dim = state->dim;

    if( action )
        *action = state->last_action;

    if( quality )
        *quality = state->quality;

    if( quality_min )
        *quality_min = state->quality_min;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_set_fgcrodr_recycle_update"
LIS_INT
lis_solver_set_fgcrodr_recycle_update(
    LIS_SOLVER solver,
    LIS_INT update)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    if( solver==NULL ||
        (update!=LIS_FALSE &&
         update!=LIS_TRUE) )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "FGCRO-DR recycle update must be true or false\n");

        return LIS_ERR_ILL_ARG;
    }

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    solver->fgcrodr_state->update = update;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_get_fgcrodr_recycle_update"
LIS_INT
lis_solver_get_fgcrodr_recycle_update(
    LIS_SOLVER solver,
    LIS_INT *update)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    if( solver==NULL || update==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid FGCRO-DR update query\n");

        return LIS_ERR_ILL_ARG;
    }

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    *update =
        solver->fgcrodr_state->update;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_set_option"
LIS_INT
lis_fgcrodr_set_option(
    LIS_SOLVER solver,
    const char *name,
    const char *value)
{
    LIS_INT err;
    LIS_INT ivalue;
    LIS_FGCRODR_STATE state;

    double rvalue;
    char extra;

    LIS_DEBUG_FUNC_IN;

    if( solver==NULL ||
        name==NULL ||
        value==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid FGCRO-DR option\n");

        return LIS_ERR_ILL_ARG;
    }

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    state = solver->fgcrodr_state;


    /* -------------------------------------------------------- */
    /* -fgcrodr_k                                               */
    /* -------------------------------------------------------- */

    if( strcmp(name,"-fgcrodr_k")==0 )
    {
#ifdef _LONG__LONG
        if( sscanf(
                value,
                "%lld %c",
                &ivalue,
                &extra)!=1 ||
            ivalue<0 )
#else
        if( sscanf(
                value,
                "%d %c",
                &ivalue,
                &extra)!=1 ||
            ivalue<0 )
#endif
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "Parameter -fgcrodr_k must be a "
                "nonnegative integer\n");

            return LIS_ERR_ILL_ARG;
        }

        if( state->target_dim!=ivalue )
        {
            lis_fgcrodr_vectors_destroy(state);
        }

        state->target_dim = ivalue;

        LIS_DEBUG_FUNC_OUT;

        return LIS_SUCCESS;
    }


    /* -------------------------------------------------------- */
    /* -fgcrodr_keep                                            */
    /* -------------------------------------------------------- */

    if( strcmp(name,"-fgcrodr_keep")==0 )
    {
        if( sscanf(
                value,
                "%lg %c",
                &rvalue,
                &extra)!=1 ||
            rvalue!=rvalue ||
            rvalue<0.0 ||
            rvalue>1.0 )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "Parameter -fgcrodr_keep must satisfy "
                "0 <= value <= 1\n");

            return LIS_ERR_ILL_ARG;
        }

        state->keep_threshold =
            (LIS_REAL)rvalue;

        LIS_DEBUG_FUNC_OUT;

        return LIS_SUCCESS;
    }


    /* -------------------------------------------------------- */
    /* -fgcrodr_reject                                          */
    /* -------------------------------------------------------- */

    if( strcmp(name,"-fgcrodr_reject")==0 )
    {
        if( sscanf(
                value,
                "%lg %c",
                &rvalue,
                &extra)!=1 ||
            rvalue!=rvalue ||
            rvalue<0.0 ||
            rvalue>1.0 )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "Parameter -fgcrodr_reject must satisfy "
                "0 <= value <= 1\n");

            return LIS_ERR_ILL_ARG;
        }

        state->reject_threshold =
            (LIS_REAL)rvalue;

        LIS_DEBUG_FUNC_OUT;

        return LIS_SUCCESS;
    }


    /* -------------------------------------------------------- */
    /* -fgcrodr_update                                          */
    /* -------------------------------------------------------- */

    if( strcmp(name,"-fgcrodr_update")==0 )
    {
        if( strcmp(value,"true")==0 ||
            strcmp(value,"TRUE")==0 ||
            strcmp(value,"1")==0 )
        {
            state->update = LIS_TRUE;
        }
        else if(
            strcmp(value,"false")==0 ||
            strcmp(value,"FALSE")==0 ||
            strcmp(value,"0")==0 )
        {
            state->update = LIS_FALSE;
        }
        else
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "Parameter -fgcrodr_update must be "
                "true or false\n");

            return LIS_ERR_ILL_ARG;
        }

        LIS_DEBUG_FUNC_OUT;

        return LIS_SUCCESS;
    }


    LIS_SETERR1(
        LIS_ERR_ILL_ARG,
        "unknown FGCRO-DR option %s\n",
        name);

    return LIS_ERR_ILL_ARG;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_check_params"
LIS_INT
lis_fgcrodr_check_params(
    LIS_SOLVER solver)
{
    LIS_INT err;
    LIS_FGCRODR_STATE state;

    LIS_DEBUG_FUNC_IN;

    err = lis_fgmres_check_params(solver);
    if( err )
        return err;

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    state = solver->fgcrodr_state;


#ifdef _COMPLEX
    /*
     * The current recycling implementation uses a real small-dense
     * harmonic Ritz eigensystem.  In particular, its projected
     * algebra is formed and interpreted as the real
     *
     *     G^T G y = theta G^T F y
     *
     * problem.
     *
     * A complex recycling core requires Hermitian projected products
     * and a complex small eigensystem.  Do not silently run the real
     * harmonic path for complex scalars.
     *
     * k=0 remains the existing FGMRES reference path.
     */
    if( state->target_dim>0 )
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "complex FGCRO-DR recycling core is not implemented yet; "
            "use -fgcrodr_k 0 for the FGMRES reference path\n");

        return LIS_ERR_NOT_IMPLEMENTED;
    }
#endif


    if(
        state->target_dim>0 &&
        solver->options[LIS_OPTIONS_RESTART]<2 )
    {
        LIS_SETERR1(
            LIS_ERR_ILL_ARG,
            "Parameter LIS_OPTIONS_RESTART(=%D) must be "
            "at least 2 for FGCRO-DR recycling\n",
            solver->options[LIS_OPTIONS_RESTART]);

        return LIS_ERR_ILL_ARG;
    }

    if( state->target_dim<0 )
    {
        LIS_SETERR1(
            LIS_ERR_ILL_ARG,
            "FGCRO-DR recycle dimension %D is negative\n",
            state->target_dim);

        return LIS_ERR_ILL_ARG;
    }

    if( state->keep_threshold<0.0 ||
        state->keep_threshold>1.0 ||
        state->reject_threshold<0.0 ||
        state->reject_threshold>1.0 ||
        state->reject_threshold >
            state->keep_threshold )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "FGCRO-DR quality thresholds must satisfy "
            "0 <= reject <= keep <= 1\n");

        return LIS_ERR_ILL_ARG;
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_malloc_work"
LIS_INT
lis_fgcrodr_malloc_work(
    LIS_SOLVER solver)
{
    /*
     * Until the augmented Arnoldi core is installed, the working
     * vectors have exactly the FGMRES layout.
     */
    return lis_fgmres_malloc_work(solver);
}


/*
 * Flexible Arnoldi core used before recycle augmentation is applied.
 *
 * This intentionally follows the existing FGMRES numerical path.
 * V stores the orthonormal Arnoldi basis and Z stores the vectors
 * produced by the preconditioner:
 *
 *     z_j = M_j^{-1} v_j
 *     w   = A z_j
 *
 * Later FGCRO-DR stages insert the C-space projection between these
 * two operations and the ordinary V-space Arnoldi orthogonalization.
 */

#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_ensure_recycle_vectors"
static LIS_INT
lis_fgcrodr_ensure_recycle_vectors(
    LIS_SOLVER solver,
    LIS_INT nvec)
{
    LIS_FGCRODR_STATE state;

    LIS_INT i;
    LIS_INT err;


    if( solver==NULL || nvec<0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid FGCRO-DR recycle allocation request\n");

        return LIS_ERR_ILL_ARG;
    }


    state =
        solver->fgcrodr_state;


    if( state==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "FGCRO-DR state is not initialized\n");

        return LIS_ERR_ILL_ARG;
    }


    if( nvec==0 )
    {
        lis_fgcrodr_vectors_destroy(
            state);

        return LIS_SUCCESS;
    }


    /*
     * Reuse a completely allocated empty recycle container.
     *
     * Stage 2B does not yet transport recycle vectors between changed
     * operators.  Compatibility/generation handling is added with the
     * cross-system policy.
     */
    if(
        state->alloc==nvec &&
        state->w!=NULL &&
        state->z!=NULL &&
        state->c!=NULL )
    {
        return LIS_SUCCESS;
    }


    lis_fgcrodr_vectors_destroy(
        state);


    state->w =
        (LIS_VECTOR *)lis_malloc(
            (size_t)nvec*sizeof(LIS_VECTOR),
            "lis_fgcrodr_ensure_recycle_vectors::w");

    state->z =
        (LIS_VECTOR *)lis_malloc(
            (size_t)nvec*sizeof(LIS_VECTOR),
            "lis_fgcrodr_ensure_recycle_vectors::z");

    state->c =
        (LIS_VECTOR *)lis_malloc(
            (size_t)nvec*sizeof(LIS_VECTOR),
            "lis_fgcrodr_ensure_recycle_vectors::c");


    if(
        state->w==NULL ||
        state->z==NULL ||
        state->c==NULL )
    {
        state->alloc = nvec;

        if( state->w )
            for(i=0;i<nvec;i++)
                state->w[i] = NULL;

        if( state->z )
            for(i=0;i<nvec;i++)
                state->z[i] = NULL;

        if( state->c )
            for(i=0;i<nvec;i++)
                state->c[i] = NULL;

        lis_fgcrodr_vectors_destroy(
            state);

        LIS_SETERR_MEM(
            (size_t)3*nvec*sizeof(LIS_VECTOR));

        return LIS_ERR_OUT_OF_MEMORY;
    }


    state->alloc = nvec;


    for(i=0;i<nvec;i++)
    {
        state->w[i] = NULL;
        state->z[i] = NULL;
        state->c[i] = NULL;
    }


    for(i=0;i<nvec;i++)
    {
        err = lis_vector_duplicate(
            solver->A,
            &state->w[i]);

        if( err )
            goto fail;


        err = lis_vector_duplicate(
            solver->A,
            &state->z[i]);

        if( err )
            goto fail;


        err = lis_vector_duplicate(
            solver->A,
            &state->c[i]);

        if( err )
            goto fail;
    }


    state->dim = 0;
    state->valid = LIS_FALSE;

    return LIS_SUCCESS;


fail:

    lis_fgcrodr_vectors_destroy(
        state);

    return err;
}



/* ---------------------------------------------------------------------- */

/*
 * Return false if NaN or infinity entered a small harmonic
 * projection.  Harmonic selection is optional, so callers can
 * retain their existing verified fallback/recycle space.
 */
static LIS_INT
lis_fgcrodr_dense_finite(
    LIS_INT n,
    const LIS_SCALAR *a)
{
    LIS_INT i;

    for(i=0;i<n*n;i++)
    {
        LIS_SCALAR d =
            a[i]-a[i];

        if( d!=d )
        {
            return LIS_FALSE;
        }
    }

    return LIS_TRUE;
}


/* Harmonic Ritz support for the flexible recycle basis                    */
/* ---------------------------------------------------------------------- */

typedef struct
{
    LIS_REAL r;
    LIS_REAL i;
} LIS_FGCRODR_CPX;


typedef struct
{
    LIS_REAL lr;
    LIS_REAL li;
    LIS_REAL score;
    LIS_INT pair;
} LIS_FGCRODR_EIG;


static LIS_FGCRODR_CPX
lis_fgcrodr_cpx(
    LIS_REAL r,
    LIS_REAL i)
{
    LIS_FGCRODR_CPX z;

    z.r = r;
    z.i = i;

    return z;
}


static LIS_REAL
lis_fgcrodr_cabs2(
    LIS_FGCRODR_CPX z)
{
    return
        z.r*z.r
        + z.i*z.i;
}


static LIS_FGCRODR_CPX
lis_fgcrodr_csub(
    LIS_FGCRODR_CPX a,
    LIS_FGCRODR_CPX b)
{
    return
        lis_fgcrodr_cpx(
            a.r-b.r,
            a.i-b.i);
}


static LIS_FGCRODR_CPX
lis_fgcrodr_cmul(
    LIS_FGCRODR_CPX a,
    LIS_FGCRODR_CPX b)
{
    return
        lis_fgcrodr_cpx(
            a.r*b.r-a.i*b.i,
            a.r*b.i+a.i*b.r);
}


static LIS_FGCRODR_CPX
lis_fgcrodr_cdiv(
    LIS_FGCRODR_CPX a,
    LIS_FGCRODR_CPX b)
{
    LIS_REAL d;

    d =
        b.r*b.r
        + b.i*b.i;

    if(
        d<=
        FGCRODR_BUILD_TINY
        *FGCRODR_BUILD_TINY )
    {
        return
            lis_fgcrodr_cpx(
                0.0,
                0.0);
    }

    return
        lis_fgcrodr_cpx(
            (a.r*b.r+a.i*b.i)/d,
            (a.i*b.r-a.r*b.i)/d);
}


static LIS_INT
lis_fgcrodr_csolve(
    LIS_INT n,
    LIS_FGCRODR_CPX *a,
    LIS_FGCRODR_CPX *b,
    LIS_FGCRODR_CPX *x)
{
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT p;

    LIS_REAL vmax;
    LIS_REAL v;

    LIS_FGCRODR_CPX tmp;
    LIS_FGCRODR_CPX factor;
    LIS_FGCRODR_CPX sum;


    for(k=0;k<n;k++)
    {
        p = k;

        vmax =
            lis_fgcrodr_cabs2(
                a[k*n+k]);


        for(i=k+1;i<n;i++)
        {
            v =
                lis_fgcrodr_cabs2(
                    a[i*n+k]);

            if( v>vmax )
            {
                vmax = v;
                p = i;
            }
        }


        if(
            vmax<=
            FGCRODR_BUILD_TINY
            *FGCRODR_BUILD_TINY )
        {
            return LIS_BREAKDOWN;
        }


        if( p!=k )
        {
            for(j=k;j<n;j++)
            {
                tmp =
                    a[k*n+j];

                a[k*n+j] =
                    a[p*n+j];

                a[p*n+j] =
                    tmp;
            }

            tmp = b[k];
            b[k] = b[p];
            b[p] = tmp;
        }


        for(i=k+1;i<n;i++)
        {
            factor =
                lis_fgcrodr_cdiv(
                    a[i*n+k],
                    a[k*n+k]);

            a[i*n+k] =
                lis_fgcrodr_cpx(
                    0.0,
                    0.0);


            for(j=k+1;j<n;j++)
            {
                a[i*n+j] =
                    lis_fgcrodr_csub(
                        a[i*n+j],
                        lis_fgcrodr_cmul(
                            factor,
                            a[k*n+j]));
            }


            b[i] =
                lis_fgcrodr_csub(
                    b[i],
                    lis_fgcrodr_cmul(
                        factor,
                        b[k]));
        }
    }


    for(i=n-1;i>=0;i--)
    {
        sum = b[i];


        for(j=i+1;j<n;j++)
        {
            sum =
                lis_fgcrodr_csub(
                    sum,
                    lis_fgcrodr_cmul(
                        a[i*n+j],
                        x[j]));
        }


        if(
            lis_fgcrodr_cabs2(
                a[i*n+i])
            <=
            FGCRODR_BUILD_TINY
            *FGCRODR_BUILD_TINY )
        {
            return LIS_BREAKDOWN;
        }


        x[i] =
            lis_fgcrodr_cdiv(
                sum,
                a[i*n+i]);
    }


    return LIS_SUCCESS;
}


static void
lis_fgcrodr_sort_eigs(
    LIS_FGCRODR_EIG *eig,
    LIS_INT n)
{
    LIS_INT i;
    LIS_INT j;

    LIS_FGCRODR_EIG x;


    for(i=1;i<n;i++)
    {
        x =
            eig[i];

        j =
            i-1;


        while(
            j>=0 &&
            eig[j].score<x.score )
        {
            eig[j+1] =
                eig[j];

            j--;
        }


        eig[j+1] =
            x;
    }
}


static LIS_INT
lis_fgcrodr_extract_eigs(
    LIS_INT n,
    LIS_SCALAR *tqr,
    LIS_FGCRODR_EIG *eig,
    LIS_INT *neig)
{
    LIS_INT i;
    LIS_INT m;

    LIS_REAL a;
    LIS_REAL b;
    LIS_REAL c;
    LIS_REAL d;

    LIS_REAL tr;
    LIS_REAL det;
    LIS_REAL disc;
    LIS_REAL s;

    LIS_REAL sub;
    LIS_REAL scale;


    m = 0;
    i = 0;


    while( i<n )
    {
        scale =
            (LIS_REAL)1.0
            + fabs(
                (LIS_REAL)
                tqr[i+i*n]);

        sub =
            (i+1<n)
            ? fabs(
                (LIS_REAL)
                tqr[(i+1)+i*n])
            : 0.0;


        if(
            i+1<n &&
            sub>
            (LIS_REAL)1.0e-10
            *scale )
        {
            a =
                (LIS_REAL)
                tqr[i+i*n];

            b =
                (LIS_REAL)
                tqr[i+(i+1)*n];

            c =
                (LIS_REAL)
                tqr[(i+1)+i*n];

            d =
                (LIS_REAL)
                tqr[
                    (i+1)
                    +(i+1)*n];


            tr =
                a+d;

            det =
                a*d-b*c;

            disc =
                tr*tr
                -(LIS_REAL)4.0*det;


            if( disc>=0.0 )
            {
                s =
                    sqrt(disc);


                eig[m].lr =
                    ((LIS_REAL)0.5)
                    *(tr+s);

                eig[m].li =
                    0.0;

                eig[m].score =
                    fabs(
                        eig[m].lr);

                eig[m].pair =
                    0;

                m++;


                eig[m].lr =
                    ((LIS_REAL)0.5)
                    *(tr-s);

                eig[m].li =
                    0.0;

                eig[m].score =
                    fabs(
                        eig[m].lr);

                eig[m].pair =
                    0;

                m++;
            }
            else
            {
                eig[m].lr =
                    ((LIS_REAL)0.5)
                    *tr;

                eig[m].li =
                    ((LIS_REAL)0.5)
                    *sqrt(-disc);

                eig[m].score =
                    sqrt(
                        eig[m].lr
                        *eig[m].lr
                        +
                        eig[m].li
                        *eig[m].li);

                eig[m].pair =
                    1;

                m++;
            }


            i += 2;
        }
        else
        {
            eig[m].lr =
                (LIS_REAL)
                tqr[i+i*n];

            eig[m].li =
                0.0;

            eig[m].score =
                fabs(
                    eig[m].lr);

            eig[m].pair =
                0;

            m++;
            i++;
        }
    }


    *neig =
        m;


    lis_fgcrodr_sort_eigs(
        eig,
        m);


    return LIS_SUCCESS;
}


static LIS_INT
lis_fgcrodr_inverse_eigenvector(
    LIS_INT n,
    const LIS_SCALAR *T,
    LIS_REAL lr,
    LIS_REAL li,
    LIS_REAL *vr,
    LIS_REAL *vi)
{
    LIS_FGCRODR_CPX *a = NULL;
    LIS_FGCRODR_CPX *rhs = NULL;
    LIS_FGCRODR_CPX *x = NULL;

    LIS_FGCRODR_CPX sigma;

    LIS_REAL delta;
    LIS_REAL nrm;
    LIS_REAL mag;

    LIS_INT i;
    LIS_INT j;
    LIS_INT it;
    LIS_INT err;


    a =
        (LIS_FGCRODR_CPX *)lis_malloc(
            (size_t)n
            *(size_t)n
            *sizeof(*a),
            "lis_fgcrodr_inverse_eigenvector::a");

    rhs =
        (LIS_FGCRODR_CPX *)lis_malloc(
            (size_t)n
            *sizeof(*rhs),
            "lis_fgcrodr_inverse_eigenvector::rhs");

    x =
        (LIS_FGCRODR_CPX *)lis_malloc(
            (size_t)n
            *sizeof(*x),
            "lis_fgcrodr_inverse_eigenvector::x");


    if(
        a==NULL ||
        rhs==NULL ||
        x==NULL )
    {
        lis_free(a);
        lis_free(rhs);
        lis_free(x);

        return LIS_ERR_OUT_OF_MEMORY;
    }


    mag =
        sqrt(
            lr*lr
            + li*li);

    delta =
        (LIS_REAL)1.0e-8
        *((LIS_REAL)1.0+mag);


    sigma.r =
        lr+delta;

    sigma.i =
        (li==0.0)
        ? 0.0
        : li
          + (
              (li>0.0)
              ? (LIS_REAL)0.125*delta
              : -(LIS_REAL)0.125*delta
            );


    for(i=0;i<n;i++)
    {
        rhs[i].r =
            (LIS_REAL)1.0
            /(LIS_REAL)(i+1);

        rhs[i].i =
            (li!=0.0)
            ? (LIS_REAL)0.125
              /(LIS_REAL)(n-i)
            : 0.0;
    }


    err =
        LIS_SUCCESS;


    for(it=0;it<5;it++)
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                a[i*n+j].r =
                    (LIS_REAL)
                    T[i+j*n];

                a[i*n+j].i =
                    0.0;
            }


            a[i*n+i] =
                lis_fgcrodr_csub(
                    a[i*n+i],
                    sigma);

            x[i] =
                lis_fgcrodr_cpx(
                    0.0,
                    0.0);
        }


        err =
            lis_fgcrodr_csolve(
                n,
                a,
                rhs,
                x);

        if( err )
            break;


        nrm =
            0.0;


        for(i=0;i<n;i++)
        {
            nrm +=
                lis_fgcrodr_cabs2(
                    x[i]);
        }


        nrm =
            sqrt(nrm);


        if(
            nrm<=FGCRODR_BUILD_TINY ||
            nrm!=nrm )
        {
            err =
                LIS_BREAKDOWN;

            break;
        }


        for(i=0;i<n;i++)
        {
            rhs[i].r =
                x[i].r/nrm;

            rhs[i].i =
                x[i].i/nrm;
        }
    }


    if( err==LIS_SUCCESS )
    {
        for(i=0;i<n;i++)
        {
            vr[i] =
                rhs[i].r;

            vi[i] =
                rhs[i].i;
        }
    }


    lis_free(a);
    lis_free(rhs);
    lis_free(x);


    return err;
}


/*
 * Build a candidate in the real solution-space basis Z.
 *
 * W receives the same Ritz coefficients applied to V only as an
 * auxiliary search-space representation.  Persistent Z is always the
 * authoritative flexible direction.
 */
static LIS_INT
lis_fgcrodr_build_harmonic_vector(
    LIS_VECTOR *v,
    LIS_VECTOR *z,
    LIS_INT jdim,
    const LIS_REAL *coef,
    LIS_VECTOR wout,
    LIS_VECTOR zout)
{
    LIS_INT i;
    LIS_INT err;

    LIS_REAL nrm;


    err =
        lis_vector_set_all(
            0.0,
            wout);

    if( err )
        return err;


    err =
        lis_vector_set_all(
            0.0,
            zout);

    if( err )
        return err;


    for(i=0;i<jdim;i++)
    {
        err =
            lis_vector_axpy(
                (LIS_SCALAR)coef[i],
                v[i],
                wout);

        if( err )
            return err;


        err =
            lis_vector_axpy(
                (LIS_SCALAR)coef[i],
                z[i],
                zout);

        if( err )
            return err;
    }


    err =
        lis_vector_nrm2(
            zout,
            &nrm);

    if( err )
        return err;


    if(
        nrm<=FGCRODR_BUILD_TINY ||
        nrm!=nrm )
    {
        return LIS_BREAKDOWN;
    }


    return LIS_SUCCESS;
}


static LIS_INT
lis_fgcrodr_install_recycle(
    LIS_SOLVER solver,
    LIS_VECTOR *v,
    LIS_VECTOR *z,
    LIS_INT jdim);


/*
 * First-generation harmonic Ritz BUILD.
 *
 * Stage 2D-A intentionally handles oldk=0 only.  Existing persistent
 * recycle vectors are consumed by the Stage-2C path and are not replaced
 * here.  The augmented oldk>0 harmonic update is added next.
 *
 * Algebra:
 *
 *     P = Z_m
 *     Q = V_{m+1}
 *
 *     A P = Q G
 *
 *     F = Q^T P = V_{m+1}^T Z_m
 *
 *     G^T G y = theta G^T F y
 *
 *     T = (G^T G)^-1 (G^T F)
 *
 * eig(T) = 1/theta.
 */
static LIS_INT
lis_fgcrodr_try_harmonic_build(
    LIS_SOLVER solver,
    LIS_VECTOR *v,
    LIS_VECTOR *z,
    LIS_INT jdim,
    const LIS_SCALAR *hraw,
    LIS_INT hld)
{
    LIS_FGCRODR_STATE state;

    LIS_SCALAR *G = NULL;
    LIS_SCALAR *F = NULL;

    LIS_SCALAR *AP = NULL;
    LIS_SCALAR *BP = NULL;

    LIS_SCALAR *T = NULL;
    LIS_SCALAR *TQR = NULL;

    LIS_SCALAR *QW = NULL;
    LIS_SCALAR *RW = NULL;
    LIS_SCALAR *AINV = NULL;

    LIS_FGCRODR_EIG *eig = NULL;

    LIS_REAL *vr = NULL;
    LIS_REAL *vi = NULL;

    LIS_VECTOR *neww = NULL;
    LIS_VECTOR *newz = NULL;

    LIS_SCALAR ag;
    LIS_SCALAR bg;
    LIS_SCALAR s;

    LIS_INT target;
    LIS_INT psmall;
    LIS_INT qsmall;

    LIS_INT i;
    LIS_INT r;
    LIS_INT c;
    LIS_INT l;

    LIS_INT neig;
    LIS_INT nsel;
    LIS_INT eidx;

    LIS_INT qriter;
    LIS_REAL qrerr;

    LIS_INT err;


    if(
        solver==NULL ||
        v==NULL ||
        z==NULL ||
        hraw==NULL ||
        jdim<=0 ||
        hld<jdim+1 )
    {
        return LIS_SUCCESS;
    }


    state =
        solver->fgcrodr_state;


    if(
        state==NULL ||
        !state->update ||
        state->valid ||
        state->target_dim<=0 )
    {
        return LIS_SUCCESS;
    }


    target =
        solver->options[
            LIS_OPTIONS_RESTART]
        - 1;

    if( target>state->target_dim )
        target = state->target_dim;

    if( target>jdim )
        target = jdim;


    if( target<=0 )
        return LIS_SUCCESS;


    psmall =
        jdim;

    qsmall =
        jdim+1;


    G =
        (LIS_SCALAR *)lis_malloc(
            (size_t)qsmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::G");

    F =
        (LIS_SCALAR *)lis_malloc(
            (size_t)qsmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::F");

    AP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::AP");

    BP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::BP");

    T =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::T");

    TQR =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::TQR");

    QW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::QW");

    RW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::RW");

    AINV =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::AINV");

    eig =
        (LIS_FGCRODR_EIG *)lis_malloc(
            (size_t)psmall
            *sizeof(*eig),
            "lis_fgcrodr_try_harmonic_build::eig");

    vr =
        (LIS_REAL *)lis_malloc(
            (size_t)psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_build::vr");

    vi =
        (LIS_REAL *)lis_malloc(
            (size_t)psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_build::vi");

    neww =
        (LIS_VECTOR *)lis_malloc(
            (size_t)target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_build::neww");

    newz =
        (LIS_VECTOR *)lis_malloc(
            (size_t)target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_build::newz");


    if( neww )
    {
        for(i=0;i<target;i++)
            neww[i] = NULL;
    }

    if( newz )
    {
        for(i=0;i<target;i++)
            newz[i] = NULL;
    }


    if(
        G==NULL ||
        F==NULL ||
        AP==NULL ||
        BP==NULL ||
        T==NULL ||
        TQR==NULL ||
        QW==NULL ||
        RW==NULL ||
        AINV==NULL ||
        eig==NULL ||
        vr==NULL ||
        vi==NULL ||
        neww==NULL ||
        newz==NULL )
    {
        err =
            LIS_ERR_OUT_OF_MEMORY;

        goto cleanup;
    }


    memset(
        G,
        0,
        (size_t)qsmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));

    memset(
        F,
        0,
        (size_t)qsmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));


    /*
     * Unrotated flexible Arnoldi relation:
     *
     *     A Z_m = V_{m+1} Hbar_m
     */
    for(c=0;c<psmall;c++)
    {
        for(r=0;r<qsmall;r++)
        {
            G[
                r+c*qsmall] =
                hraw[
                    r+c*hld];
        }
    }


    /*
     * Flexible-domain overlap:
     *
     *     F = V_{m+1}^T Z_m
     *
     * Unlike ordinary GCRO-DR, this is not an identity block because
     * flexible Z_j are not the Arnoldi vectors V_j.
     */
    for(c=0;c<psmall;c++)
    {
        for(r=0;r<qsmall;r++)
        {
            err =
                lis_vector_dot(
                    v[r],
                    z[c],
                    &F[
                        r
                        +c*qsmall]);

            if( err )
                goto cleanup;
        }
    }


    /*
     * AP = G^T G
     * BP = G^T F
     */
    for(c=0;c<psmall;c++)
    {
        for(r=0;r<psmall;r++)
        {
            ag =
                0.0;

            bg =
                0.0;


            for(l=0;l<qsmall;l++)
            {
                ag +=
                    G[
                        l+r*qsmall]
                    *G[
                        l+c*qsmall];

                bg +=
                    G[
                        l+r*qsmall]
                    *F[
                        l+c*qsmall];
            }


            AP[
                r+c*psmall] =
                ag;

            BP[
                r+c*psmall] =
                bg;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            psmall,
            AP)
        ||
        !lis_fgcrodr_dense_finite(
            psmall,
            BP) )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        AINV,
        AP,
        (size_t)psmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));


    err =
        lis_array_ge(
            psmall,
            AINV);

    if(
        err
        ||
        !lis_fgcrodr_dense_finite(
            psmall,
            AINV) )
    {
        /*
         * Harmonic selection is optional.  The caller will fall back
         * to the already verified Stage-2B candidate selection.
         */
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    /*
     * T = AP^-1 BP
     */
    for(c=0;c<psmall;c++)
    {
        for(r=0;r<psmall;r++)
        {
            s =
                0.0;


            for(l=0;l<psmall;l++)
            {
                s +=
                    AINV[
                        r+l*psmall]
                    *BP[
                        l+c*psmall];
            }


            T[
                r+c*psmall] =
                s;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            psmall,
            T) )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        TQR,
        T,
        (size_t)psmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));


    /*
     * Harmonic recycle extraction is optional.  Bound the
     * unshifted QR iteration so a difficult small eigensystem
     * cannot dominate the linear solve.
     */
    qriter =
        1024;

    err =
        lis_array_qr_limit(
            psmall,
            TQR,
            QW,
            RW,
            qriter,
            &qriter,
            &qrerr);

    if( err )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    (void)qriter;
    (void)qrerr;


    err =
        lis_fgcrodr_extract_eigs(
            psmall,
            TQR,
            eig,
            &neig);

    if( err )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    for(i=0;i<target;i++)
    {
        err =
            lis_vector_duplicate(
                solver->A,
                &neww[i]);

        if( err )
            goto cleanup;


        err =
            lis_vector_duplicate(
                solver->A,
                &newz[i]);

        if( err )
            goto cleanup;
    }


    nsel =
        0;


    for(
        eidx=0;
        eidx<neig && nsel<target;
        eidx++ )
    {
        if(
            eig[eidx].pair &&
            nsel+2>target )
        {
            continue;
        }


        err =
            lis_fgcrodr_inverse_eigenvector(
                psmall,
                T,
                eig[eidx].lr,
                eig[eidx].li,
                vr,
                vi);


        if( err==LIS_ERR_OUT_OF_MEMORY )
            goto cleanup;

        if( err )
        {
            err =
                LIS_SUCCESS;

            continue;
        }


        if( eig[eidx].pair )
        {
            err =
                lis_fgcrodr_build_harmonic_vector(
                    v,
                    z,
                    jdim,
                    vr,
                    neww[nsel],
                    newz[nsel]);

            if( err )
            {
                if(
                    err==
                    LIS_ERR_OUT_OF_MEMORY )
                {
                    goto cleanup;
                }

                err =
                    LIS_SUCCESS;

                continue;
            }


            err =
                lis_fgcrodr_build_harmonic_vector(
                    v,
                    z,
                    jdim,
                    vi,
                    neww[nsel+1],
                    newz[nsel+1]);

            if( err )
            {
                if(
                    err==
                    LIS_ERR_OUT_OF_MEMORY )
                {
                    goto cleanup;
                }

                err =
                    LIS_SUCCESS;

                continue;
            }


            nsel +=
                2;
        }
        else
        {
            err =
                lis_fgcrodr_build_harmonic_vector(
                    v,
                    z,
                    jdim,
                    vr,
                    neww[nsel],
                    newz[nsel]);

            if( err )
            {
                if(
                    err==
                    LIS_ERR_OUT_OF_MEMORY )
                {
                    goto cleanup;
                }

                err =
                    LIS_SUCCESS;

                continue;
            }


            nsel++;
        }
    }


    /*
     * Require a complete requested space at this stage.
     *
     * If the small eigensystem is rank-deficient or an odd target
     * cannot be filled because only complex pairs remain, retain the
     * proven Stage-2B fallback.
     */
    if( nsel!=target )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    err =
        lis_fgcrodr_install_recycle(
            solver,
            neww,
            newz,
            nsel);

    if( err )
        goto cleanup;


    state =
        solver->fgcrodr_state;


    if(
        state!=NULL &&
        state->valid &&
        state->dim>0 )
    {
    }


cleanup:

    if( neww )
    {
        for(i=0;i<target;i++)
        {
            if( neww[i] )
                lis_vector_destroy(
                    neww[i]);
        }
    }


    if( newz )
    {
        for(i=0;i<target;i++)
        {
            if( newz[i] )
                lis_vector_destroy(
                    newz[i]);
        }
    }


    lis_free(neww);
    lis_free(newz);

    lis_free(G);
    lis_free(F);

    lis_free(AP);
    lis_free(BP);

    lis_free(T);
    lis_free(TQR);

    lis_free(QW);
    lis_free(RW);
    lis_free(AINV);

    lis_free(eig);

    lis_free(vr);
    lis_free(vi);


    return err;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_install_recycle"
static LIS_INT
lis_fgcrodr_install_recycle(
    LIS_SOLVER solver,
    LIS_VECTOR *v,
    LIS_VECTOR *z,
    LIS_INT jdim)
{
    LIS_FGCRODR_STATE state;

    LIS_INT err;

    LIS_INT target;
    LIS_INT candidate;
    LIS_INT keep;
    LIS_INT j;
    LIS_INT pass;

    LIS_SCALAR alpha;

    LIS_REAL nrm;


    if(
        solver==NULL ||
        v==NULL ||
        z==NULL ||
        jdim<0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid FGCRO-DR recycle install arguments\n");

        return LIS_ERR_ILL_ARG;
    }


    state =
        solver->fgcrodr_state;


    if( state==NULL )
        return LIS_SUCCESS;




    /*
     * update=false means that the current solve is not allowed to
     * create or replace the persistent recycle space.
     */
    if( !state->update )
        return LIS_SUCCESS;


    /*
     * Installation is only allowed when no persistent recycle space
     * is currently active.  Existing state is consumed by the reuse
     * path and is updated by a later augmented harmonic stage.
     */
    if(
        state->valid &&
        state->dim>0 )
    {
        return LIS_SUCCESS;
    }


    if(
        state->target_dim<=0 ||
        jdim<=0 )
    {
        return LIS_SUCCESS;
    }


    /*
     * Reserve at least one Arnoldi direction outside the requested
     * recycle space, matching the existing GCRO-DR dimensional rule.
     */
    target =
        solver->options[LIS_OPTIONS_RESTART]-1;

    if( target>state->target_dim )
        target = state->target_dim;

    if( target>jdim )
        target = jdim;


    if( target<=0 )
        return LIS_SUCCESS;


    err =
        lis_fgcrodr_ensure_recycle_vectors(
            solver,
            target);

    if( err )
        return err;


    state =
        solver->fgcrodr_state;


    keep = 0;


    /*
     * Install the supplied candidate pairs.  They may come from the
     * harmonic-Ritz selector or from the provisional fallback.
     */
    for(
        candidate=0;
        candidate<jdim && keep<target;
        candidate++ )
    {
        /*
         * W is retained as the corresponding Arnoldi-space
         * representation.
         *
         * Z is the actual preconditioned direction and is the
         * authoritative solution-space recycle vector.
         */
        err =
            lis_vector_copy(
                v[candidate],
                state->w[keep]);

        if( err )
            goto fail;


        err =
            lis_vector_copy(
                z[candidate],
                state->z[keep]);

        if( err )
            goto fail;


        /*
         * Start with the exact operator image.
         *
         *     C_i = A Z_i
         */
        lis_matvec(
            solver->A,
            state->z[keep],
            state->c[keep]);


        /*
         * Two-pass modified Gram-Schmidt.
         *
         * Every operation applied to C is also applied to Z.
         * Therefore the relation
         *
         *     A Z_i = C_i
         *
         * is preserved through orthogonalization.
         *
         * W receives the same small-space transformation but is
         * never used to reconstruct Z.
         */
        for(pass=0;pass<2;pass++)
        {
            for(j=0;j<keep;j++)
            {
                lis_vector_dot(
                    state->c[keep],
                    state->c[j],
                    &alpha);


                lis_vector_axpy(
                    -alpha,
                    state->c[j],
                    state->c[keep]);


                lis_vector_axpy(
                    -alpha,
                    state->z[j],
                    state->z[keep]);


                lis_vector_axpy(
                    -alpha,
                    state->w[j],
                    state->w[keep]);
            }
        }


        lis_vector_nrm2(
            state->c[keep],
            &nrm);


        /*
         * Drop rank-deficient candidates rather than poisoning the
         * persistent space.
         */
        if(
            nrm<=FGCRODR_BUILD_TINY ||
            nrm!=nrm )
        {

            continue;
        }


        lis_vector_scale(
            1.0/nrm,
            state->c[keep]);


        lis_vector_scale(
            1.0/nrm,
            state->z[keep]);


        lis_vector_scale(
            1.0/nrm,
            state->w[keep]);


        keep++;

    }




    state->dim =
        keep;

    state->valid =
        (keep>0)
        ? LIS_TRUE
        : LIS_FALSE;


    if( keep>0 )
    {
        state->matrix =
            solver->A;

        state->structure_generation =
            solver->A->structure_generation;

        state->value_generation =
            solver->A->value_generation;
    }


    if( keep>0 )
    {
        state->last_action =
            LIS_FGCRODR_RECYCLE_BUILD;

        state->builds++;

        /*
         * Cross-system quality is introduced later.
         *
         * A freshly constructed exact A*Z=C relation starts at
         * nominal quality 1.
         */
        state->quality =
            1.0;

        state->quality_min =
            1.0;
    }
    else
    {
        state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->rejects++;

        state->quality =
            0.0;

        state->quality_min =
            0.0;
    }


    return LIS_SUCCESS;


fail:

    lis_fgcrodr_vectors_destroy(
        state);

    state->last_action =
        LIS_FGCRODR_RECYCLE_REJECT;

    state->rejects++;

    return err;
}





/*
 * Build one harmonic candidate from the augmented flexible domain basis
 *
 *     [ Z_old , Z_krylov ].
 *
 * W follows the same small coefficients only as an auxiliary
 * search-space representation.  Z remains authoritative.
 */
static LIS_INT
lis_fgcrodr_build_augmented_vector(
    LIS_FGCRODR_STATE state,
    LIS_VECTOR *v,
    LIS_VECTOR *z,
    LIS_INT jdim,
    LIS_INT oldk,
    const LIS_REAL *coef,
    LIS_VECTOR wout,
    LIS_VECTOR zout)
{
    LIS_INT i;
    LIS_INT err;

    LIS_REAL nrm;


    err =
        lis_vector_set_all(
            0.0,
            wout);

    if( err )
        return err;


    err =
        lis_vector_set_all(
            0.0,
            zout);

    if( err )
        return err;


    /*
     * Persistent part.
     */
    for(i=0;i<oldk;i++)
    {
        err =
            lis_vector_axpy(
                (LIS_SCALAR)coef[i],
                state->w[i],
                wout);

        if( err )
            return err;


        err =
            lis_vector_axpy(
                (LIS_SCALAR)coef[i],
                state->z[i],
                zout);

        if( err )
            return err;
    }


    /*
     * Current flexible Krylov part.
     */
    for(i=0;i<jdim;i++)
    {
        err =
            lis_vector_axpy(
                (LIS_SCALAR)coef[oldk+i],
                v[i],
                wout);

        if( err )
            return err;


        err =
            lis_vector_axpy(
                (LIS_SCALAR)coef[oldk+i],
                z[i],
                zout);

        if( err )
            return err;
    }


    err =
        lis_vector_nrm2(
            zout,
            &nrm);

    if( err )
        return err;


    if(
        nrm<=FGCRODR_BUILD_TINY ||
        nrm!=nrm )
    {
        return LIS_BREAKDOWN;
    }


    return LIS_SUCCESS;
}


/*
 * Atomically replace an existing recycle space.
 *
 * Candidate installation is first performed in a private trial state.
 * The currently usable recycle relation is therefore retained if the
 * candidate set turns out to be rank deficient.
 */
static LIS_INT
lis_fgcrodr_replace_recycle_basis(
    LIS_SOLVER solver,
    LIS_VECTOR *neww,
    LIS_VECTOR *newz,
    LIS_INT nsel,
    LIS_INT count_update,
    LIS_INT action,
    LIS_REAL quality,
    LIS_REAL quality_min)
{
    LIS_FGCRODR_STATE oldstate;
    LIS_FGCRODR_STATE trial;

    LIS_INT saved_cross_action;

    LIS_REAL saved_quality;
    LIS_REAL saved_quality_min;

    LIS_INT err;


    oldstate =
        solver->fgcrodr_state;


    saved_cross_action =
        oldstate
        ? oldstate->cross_action
        : LIS_FGCRODR_RECYCLE_NONE;

    saved_quality =
        oldstate
        ? oldstate->quality
        : 0.0;

    saved_quality_min =
        oldstate
        ? oldstate->quality_min
        : 0.0;


    if(
        oldstate==NULL ||
        nsel<=0 )
    {
        return LIS_SUCCESS;
    }


    trial =
        (LIS_FGCRODR_STATE)lis_malloc(
            sizeof(
                struct LIS_FGCRODR_STATE_STRUCT),
            "lis_fgcrodr_replace_recycle_basis::trial");

    if( trial==NULL )
    {
        LIS_SETERR_MEM(
            sizeof(
                struct LIS_FGCRODR_STATE_STRUCT));

        return LIS_ERR_OUT_OF_MEMORY;
    }


    lis_fgcrodr_state_defaults(
        trial);


    /*
     * Preserve user configuration in the trial installer.
     */
    trial->target_dim =
        oldstate->target_dim;

    trial->keep_threshold =
        oldstate->keep_threshold;

    trial->reject_threshold =
        oldstate->reject_threshold;

    trial->update =
        oldstate->update;


    /*
     * Install candidates into an isolated state.
     */
    solver->fgcrodr_state =
        trial;


    err =
        lis_fgcrodr_install_recycle(
            solver,
            neww,
            newz,
            nsel);


    solver->fgcrodr_state =
        oldstate;


    if( err )
    {
        lis_fgcrodr_vectors_destroy(
            trial);

        lis_free(
            trial);

        return err;
    }


    /*
     * A failed/rank-deficient harmonic candidate is nonfatal.
     * Keep the previous valid recycle space.
     */
    if(
        !trial->valid ||
        trial->dim<=0 )
    {
        lis_fgcrodr_vectors_destroy(
            trial);

        lis_free(
            trial);

        return LIS_SUCCESS;
    }


    /*
     * Commit atomically.
     *
     * lis_fgcrodr_vectors_destroy() deliberately leaves configuration
     * and lifetime counters intact.
     */
    lis_fgcrodr_vectors_destroy(
        oldstate);


    oldstate->w =
        trial->w;

    oldstate->z =
        trial->z;

    oldstate->c =
        trial->c;

    oldstate->dim =
        trial->dim;

    oldstate->alloc =
        trial->alloc;

    oldstate->valid =
        trial->valid;


    oldstate->matrix =
        trial->matrix;

    oldstate->structure_generation =
        trial->structure_generation;

    oldstate->value_generation =
        trial->value_generation;


    /*
     * Preserve cross-system diagnostics across same-solve harmonic
     * replacements.  For ordinary same-generation updates use the
     * diagnostics supplied by the caller.
     */
    if(
        saved_cross_action
        !=LIS_FGCRODR_RECYCLE_NONE )
    {
        oldstate->quality =
            saved_quality;

        oldstate->quality_min =
            saved_quality_min;

        oldstate->last_action =
            saved_cross_action;

        oldstate->cross_action =
            saved_cross_action;
    }
    else
    {
        oldstate->quality =
            quality;

        oldstate->quality_min =
            quality_min;

        oldstate->last_action =
            action;

        oldstate->cross_action =
            LIS_FGCRODR_RECYCLE_NONE;
    }


    if( count_update )
        oldstate->updates++;


    /*
     * Detach vectors before freeing the temporary state.
     */
    trial->w = NULL;
    trial->z = NULL;
    trial->c = NULL;

    trial->dim = 0;
    trial->alloc = 0;
    trial->valid = LIS_FALSE;


    lis_free(
        trial);


    return LIS_SUCCESS;
}



static LIS_INT
lis_fgcrodr_commit_recycle_update(
    LIS_SOLVER solver,
    LIS_VECTOR *neww,
    LIS_VECTOR *newz,
    LIS_INT nsel)
{
    return
        lis_fgcrodr_replace_recycle_basis(
            solver,
            neww,
            newz,
            nsel,
            LIS_TRUE,
            LIS_FGCRODR_RECYCLE_KEEP,
            1.0,
            1.0);
}


/*
 * Augmented harmonic Ritz UPDATE.
 *
 * Existing flexible recycle basis:
 *
 *     Z_o = [Z_1,...,Z_k]
 *
 * Current flexible Krylov basis:
 *
 *     Z_m = [z_0,...,z_{m-1}]
 *
 * Orthogonal image basis:
 *
 *     Q = [C_o,V_{m+1}]
 *
 * with
 *
 *     A [Z_o,Z_m] = Q G
 *
 * and
 *
 *             [ I  B ]
 *     G =     [      ].
 *             [ 0  H ]
 *
 * The flexible harmonic overlap is computed explicitly:
 *
 *     F = Q^T [Z_o,Z_m].
 *
 * No current or future preconditioner is applied to reconstruct any
 * historical Z vector.
 */
static LIS_INT
lis_fgcrodr_try_harmonic_update(
    LIS_SOLVER solver,
    LIS_VECTOR *v,
    LIS_VECTOR *z,
    LIS_INT jdim,
    const LIS_SCALAR *hraw,
    LIS_INT hld,
    const LIS_SCALAR *bcoef,
    LIS_INT bld)
{
    LIS_FGCRODR_STATE state;

    LIS_SCALAR *G = NULL;
    LIS_SCALAR *F = NULL;

    LIS_SCALAR *AP = NULL;
    LIS_SCALAR *BP = NULL;

    LIS_SCALAR *T = NULL;
    LIS_SCALAR *TQR = NULL;

    LIS_SCALAR *QW = NULL;
    LIS_SCALAR *RW = NULL;
    LIS_SCALAR *AINV = NULL;

    LIS_FGCRODR_EIG *eig = NULL;

    LIS_REAL *vr = NULL;
    LIS_REAL *vi = NULL;

    LIS_VECTOR *neww = NULL;
    LIS_VECTOR *newz = NULL;

    LIS_SCALAR ag;
    LIS_SCALAR bg;
    LIS_SCALAR s;

    LIS_INT oldk;
    LIS_INT target;

    LIS_INT psmall;
    LIS_INT qsmall;

    LIS_INT i;
    LIS_INT r;
    LIS_INT c;
    LIS_INT l;

    LIS_INT neig;
    LIS_INT nsel;
    LIS_INT eidx;

    LIS_INT qriter;
    LIS_REAL qrerr;

    LIS_INT err;


    if(
        solver==NULL ||
        v==NULL ||
        z==NULL ||
        hraw==NULL ||
        bcoef==NULL ||
        jdim<=0 )
    {
        return LIS_SUCCESS;
    }


    state =
        solver->fgcrodr_state;


    if(
        state==NULL ||
        !state->update ||
        !state->valid ||
        state->dim<=0 )
    {
        return LIS_SUCCESS;
    }


    oldk =
        state->dim;


    if(
        bld<oldk ||
        hld<jdim+1 )
    {
        return LIS_SUCCESS;
    }


    target =
        solver->options[
            LIS_OPTIONS_RESTART]
        - 1;

    if( target>state->target_dim )
        target = state->target_dim;


    /*
     * The update should preserve a useful nonzero recycle rank.
     */
    if( target<=0 )
        return LIS_SUCCESS;


    psmall =
        oldk+jdim;

    qsmall =
        oldk+jdim+1;


    if( target>psmall )
        target = psmall;


    G =
        (LIS_SCALAR *)lis_malloc(
            (size_t)qsmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::G");

    F =
        (LIS_SCALAR *)lis_malloc(
            (size_t)qsmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::F");

    AP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::AP");

    BP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::BP");

    T =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::T");

    TQR =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::TQR");

    QW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::QW");

    RW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::RW");

    AINV =
        (LIS_SCALAR *)lis_malloc(
            (size_t)psmall
            *(size_t)psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::AINV");

    eig =
        (LIS_FGCRODR_EIG *)lis_malloc(
            (size_t)psmall
            *sizeof(*eig),
            "lis_fgcrodr_try_harmonic_update::eig");

    vr =
        (LIS_REAL *)lis_malloc(
            (size_t)psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_update::vr");

    vi =
        (LIS_REAL *)lis_malloc(
            (size_t)psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_update::vi");

    neww =
        (LIS_VECTOR *)lis_malloc(
            (size_t)target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_update::neww");

    newz =
        (LIS_VECTOR *)lis_malloc(
            (size_t)target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_update::newz");


    if( neww )
    {
        for(i=0;i<target;i++)
            neww[i] = NULL;
    }

    if( newz )
    {
        for(i=0;i<target;i++)
            newz[i] = NULL;
    }


    if(
        G==NULL ||
        F==NULL ||
        AP==NULL ||
        BP==NULL ||
        T==NULL ||
        TQR==NULL ||
        QW==NULL ||
        RW==NULL ||
        AINV==NULL ||
        eig==NULL ||
        vr==NULL ||
        vi==NULL ||
        neww==NULL ||
        newz==NULL )
    {
        err =
            LIS_ERR_OUT_OF_MEMORY;

        goto cleanup;
    }


    memset(
        G,
        0,
        (size_t)qsmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));

    memset(
        F,
        0,
        (size_t)qsmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));


    /*
     * Old recycle columns:
     *
     *     A Z_old = C_old.
     */
    for(i=0;i<oldk;i++)
    {
        G[
            i+i*qsmall] =
            1.0;
    }


    /*
     * Current flexible Krylov columns:
     *
     *     A z_j =
     *         C B(:,j)
     *       + Vnext H(:,j).
     */
    for(c=0;c<jdim;c++)
    {
        for(r=0;r<oldk;r++)
        {
            G[
                r
                +(oldk+c)*qsmall] =
                bcoef[
                    r+c*bld];
        }


        for(r=0;r<=jdim;r++)
        {
            G[
                oldk+r
                +(oldk+c)*qsmall] =
                hraw[
                    r+c*hld];
        }
    }


    /*
     * Full flexible overlap:
     *
     *     F = [C,Vnext]^T [Z_old,Z_krylov].
     *
     * We intentionally compute every block explicitly.
     */

    /* old domain columns */
    for(c=0;c<oldk;c++)
    {
        for(r=0;r<oldk;r++)
        {
            err =
                lis_vector_dot(
                    state->c[r],
                    state->z[c],
                    &F[
                        r+c*qsmall]);

            if( err )
                goto cleanup;
        }


        for(r=0;r<=jdim;r++)
        {
            err =
                lis_vector_dot(
                    v[r],
                    state->z[c],
                    &F[
                        oldk+r
                        +c*qsmall]);

            if( err )
                goto cleanup;
        }
    }


    /* new flexible Krylov domain columns */
    for(c=0;c<jdim;c++)
    {
        for(r=0;r<oldk;r++)
        {
            err =
                lis_vector_dot(
                    state->c[r],
                    z[c],
                    &F[
                        r
                        +(oldk+c)*qsmall]);

            if( err )
                goto cleanup;
        }


        for(r=0;r<=jdim;r++)
        {
            err =
                lis_vector_dot(
                    v[r],
                    z[c],
                    &F[
                        oldk+r
                        +(oldk+c)*qsmall]);

            if( err )
                goto cleanup;
        }
    }


    /*
     * AP = G^T G
     * BP = G^T F
     */
    for(c=0;c<psmall;c++)
    {
        for(r=0;r<psmall;r++)
        {
            ag = 0.0;
            bg = 0.0;


            for(l=0;l<qsmall;l++)
            {
                ag +=
                    G[
                        l+r*qsmall]
                    *G[
                        l+c*qsmall];

                bg +=
                    G[
                        l+r*qsmall]
                    *F[
                        l+c*qsmall];
            }


            AP[
                r+c*psmall] =
                ag;

            BP[
                r+c*psmall] =
                bg;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            psmall,
            AP)
        ||
        !lis_fgcrodr_dense_finite(
            psmall,
            BP) )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        AINV,
        AP,
        (size_t)psmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));


    err =
        lis_array_ge(
            psmall,
            AINV);

    if(
        err
        ||
        !lis_fgcrodr_dense_finite(
            psmall,
            AINV) )
    {
        /*
         * Failure of recycle enhancement is nonfatal.
         */
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    /*
     * T = (G^T G)^-1 (G^T F)
     */
    for(c=0;c<psmall;c++)
    {
        for(r=0;r<psmall;r++)
        {
            s = 0.0;


            for(l=0;l<psmall;l++)
            {
                s +=
                    AINV[
                        r+l*psmall]
                    *BP[
                        l+c*psmall];
            }


            T[
                r+c*psmall] =
                s;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            psmall,
            T) )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        TQR,
        T,
        (size_t)psmall
        *(size_t)psmall
        *sizeof(LIS_SCALAR));


    /*
     * Harmonic recycle extraction is optional.  Bound the
     * unshifted QR iteration so a difficult small eigensystem
     * cannot dominate the linear solve.
     */
    qriter =
        512;

    err =
        lis_array_qr_limit(
            psmall,
            TQR,
            QW,
            RW,
            qriter,
            &qriter,
            &qrerr);

    if( err )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    (void)qriter;
    (void)qrerr;


    err =
        lis_fgcrodr_extract_eigs(
            psmall,
            TQR,
            eig,
            &neig);

    if( err )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    for(i=0;i<target;i++)
    {
        err =
            lis_vector_duplicate(
                solver->A,
                &neww[i]);

        if( err )
            goto cleanup;


        err =
            lis_vector_duplicate(
                solver->A,
                &newz[i]);

        if( err )
            goto cleanup;
    }


    nsel = 0;


    for(
        eidx=0;
        eidx<neig && nsel<target;
        eidx++ )
    {
        if(
            eig[eidx].pair &&
            nsel+2>target )
        {
            continue;
        }


        err =
            lis_fgcrodr_inverse_eigenvector(
                psmall,
                T,
                eig[eidx].lr,
                eig[eidx].li,
                vr,
                vi);


        if( err==LIS_ERR_OUT_OF_MEMORY )
            goto cleanup;

        if( err )
        {
            err =
                LIS_SUCCESS;

            continue;
        }


        if( eig[eidx].pair )
        {
            err =
                lis_fgcrodr_build_augmented_vector(
                    state,
                    v,
                    z,
                    jdim,
                    oldk,
                    vr,
                    neww[nsel],
                    newz[nsel]);

            if( err )
            {
                if(
                    err==
                    LIS_ERR_OUT_OF_MEMORY )
                {
                    goto cleanup;
                }

                err =
                    LIS_SUCCESS;

                continue;
            }


            err =
                lis_fgcrodr_build_augmented_vector(
                    state,
                    v,
                    z,
                    jdim,
                    oldk,
                    vi,
                    neww[nsel+1],
                    newz[nsel+1]);

            if( err )
            {
                if(
                    err==
                    LIS_ERR_OUT_OF_MEMORY )
                {
                    goto cleanup;
                }

                err =
                    LIS_SUCCESS;

                continue;
            }


            nsel += 2;
        }
        else
        {
            err =
                lis_fgcrodr_build_augmented_vector(
                    state,
                    v,
                    z,
                    jdim,
                    oldk,
                    vr,
                    neww[nsel],
                    newz[nsel]);

            if( err )
            {
                if(
                    err==
                    LIS_ERR_OUT_OF_MEMORY )
                {
                    goto cleanup;
                }

                err =
                    LIS_SUCCESS;

                continue;
            }


            nsel++;
        }
    }


    /*
     * Do not replace a healthy persistent recycle space by a
     * partially filled harmonic selection.
     */
    if( nsel!=target )
    {
        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    err =
        lis_fgcrodr_commit_recycle_update(
            solver,
            neww,
            newz,
            nsel);

    if( err )
        goto cleanup;


    state =
        solver->fgcrodr_state;


    if(
        state!=NULL &&
        state->valid &&
        state->dim>0 )
    {
    }


cleanup:

    if( neww )
    {
        for(i=0;i<target;i++)
        {
            if( neww[i] )
                lis_vector_destroy(
                    neww[i]);
        }
    }


    if( newz )
    {
        for(i=0;i<target;i++)
        {
            if( newz[i] )
                lis_vector_destroy(
                    newz[i]);
        }
    }


    lis_free(neww);
    lis_free(newz);

    lis_free(G);
    lis_free(F);

    lis_free(AP);
    lis_free(BP);

    lis_free(T);
    lis_free(TQR);

    lis_free(QW);
    lis_free(RW);
    lis_free(AINV);

    lis_free(eig);

    lis_free(vr);
    lis_free(vi);


    return err;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_build_recycle"
static LIS_INT
lis_fgcrodr_build_recycle(
    LIS_SOLVER solver,
    LIS_VECTOR *v,
    LIS_VECTOR *z,
    LIS_INT jdim,
    const LIS_SCALAR *hraw,
    LIS_INT hld,
    const LIS_SCALAR *bcoef,
    LIS_INT bld)
{
    LIS_FGCRODR_STATE state;

    LIS_INT err;


    if( solver==NULL )
        return LIS_ERR_ILL_ARG;


    state =
        solver->fgcrodr_state;


    if(
        state==NULL ||
        !state->update )
    {
        return LIS_SUCCESS;
    }


    /*
     * Existing same-operator recycle space:
     *
     * perform the augmented flexible harmonic update.
     */
    if(
        state->valid &&
        state->dim>0 )
    {
        return
            lis_fgcrodr_try_harmonic_update(
                solver,
                v,
                z,
                jdim,
                hraw,
                hld,
                bcoef,
                bld);
    }


    /*
     * No recycle space yet:
     *
     * preferred Stage-2D harmonic BUILD.
     */
    err =
        lis_fgcrodr_try_harmonic_build(
            solver,
            v,
            z,
            jdim,
            hraw,
            hld);

    if( err )
        return err;


    state =
        solver->fgcrodr_state;


    if(
        state!=NULL &&
        state->valid &&
        state->dim>0 )
    {
        return LIS_SUCCESS;
    }


    /*
     * Nonfatal fallback to the verified Stage-2B initial selection.
     */
    return
        lis_fgcrodr_install_recycle(
            solver,
            v,
            z,
            jdim);
}




LIS_INT
lis_fgcrodr_get_lifecycle_state(
    LIS_SOLVER solver,
    LIS_INT *valid,
    LIS_MATRIX *matrix,
    LIS_UNSIGNED_INT *structure_generation,
    LIS_UNSIGNED_INT *value_generation)
{
    LIS_FGCRODR_STATE state;


    if( solver==NULL )
        return LIS_ERR_ILL_ARG;


    state =
        solver->fgcrodr_state;


    if( valid )
        *valid =
            state
            ? state->valid
            : LIS_FALSE;

    if( matrix )
        *matrix =
            state
            ? state->matrix
            : NULL;

    if( structure_generation )
        *structure_generation =
            state
            ? state->structure_generation
            : 0;

    if( value_generation )
        *value_generation =
            state
            ? state->value_generation
            : 0;


    return LIS_SUCCESS;
}


LIS_INT
lis_fgcrodr_get_stats(
    LIS_SOLVER solver,
    LIS_UNSIGNED_INT *builds,
    LIS_UNSIGNED_INT *updates,
    LIS_UNSIGNED_INT *uses,
    LIS_UNSIGNED_INT *rejects)
{
    LIS_FGCRODR_STATE state;


    if( solver==NULL )
        return LIS_ERR_ILL_ARG;


    state =
        solver->fgcrodr_state;


    if( builds )
        *builds =
            state
            ? state->builds
            : 0;

    if( updates )
        *updates =
            state
            ? state->updates
            : 0;

    if( uses )
        *uses =
            state
            ? state->uses
            : 0;

    if( rejects )
        *rejects =
            state
            ? state->rejects
            : 0;


    return LIS_SUCCESS;
}



LIS_INT
lis_fgcrodr_get_policy_stats(
    LIS_SOLVER solver,
    LIS_UNSIGNED_INT *cross_uses,
    LIS_UNSIGNED_INT *keeps,
    LIS_UNSIGNED_INT *transports,
    LIS_UNSIGNED_INT *rejects)
{
    LIS_FGCRODR_STATE state;


    if( solver==NULL )
        return LIS_ERR_ILL_ARG;


    state =
        solver->fgcrodr_state;


    if( cross_uses )
        *cross_uses =
            state
            ? state->cross_uses
            : 0;

    if( keeps )
        *keeps =
            state
            ? state->policy_keeps
            : 0;

    if( transports )
        *transports =
            state
            ? state->policy_transports
            : 0;

    if( rejects )
        *rejects =
            state
            ? state->policy_rejects
            : 0;


    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_validate_recycle"
LIS_INT
lis_fgcrodr_validate_recycle(
    LIS_SOLVER solver,
    LIS_REAL *orth_error,
    LIS_REAL *relation_error)
{
    LIS_FGCRODR_STATE state;

    LIS_VECTOR az = NULL;

    LIS_SCALAR dot;

    LIS_REAL err_orth;
    LIS_REAL err_relation;
    LIS_REAL value;
    LIS_REAL nrm;

    LIS_INT i;
    LIS_INT j;
    LIS_INT err;


    LIS_DEBUG_FUNC_IN;


    if(
        solver==NULL ||
        orth_error==NULL ||
        relation_error==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid FGCRO-DR recycle validation arguments\n");

        return LIS_ERR_ILL_ARG;
    }


    *orth_error = 0.0;
    *relation_error = 0.0;


    state =
        solver->fgcrodr_state;


    if(
        state==NULL ||
        !state->valid ||
        state->dim<=0 )
    {
        LIS_DEBUG_FUNC_OUT;

        return LIS_SUCCESS;
    }


    /*
     * Check
     *
     *     C^T C = I
     */
    err_orth = 0.0;


    for(i=0;i<state->dim;i++)
    {
        for(j=0;j<state->dim;j++)
        {
            err =
                lis_vector_dot(
                    state->c[i],
                    state->c[j],
                    &dot);

            if( err )
                return err;


            value =
                fabs(
                    (double)(
                        dot
                        - (i==j
                           ? (LIS_SCALAR)1.0
                           : (LIS_SCALAR)0.0)));


            if( value>err_orth )
                err_orth = value;
        }
    }


    /*
     * Check
     *
     *     A Z_i = C_i
     *
     * Since every C_i is normalized, the absolute residual is also
     * directly meaningful relative to the unit image scale.
     */
    err =
        lis_vector_duplicate(
            solver->A,
            &az);

    if( err )
        return err;


    err_relation = 0.0;


    for(i=0;i<state->dim;i++)
    {
        lis_matvec(
            solver->A,
            state->z[i],
            az);


        err =
            lis_vector_axpy(
                (LIS_SCALAR)-1.0,
                state->c[i],
                az);

        if( err )
        {
            lis_vector_destroy(az);
            return err;
        }


        err =
            lis_vector_nrm2(
                az,
                &nrm);

        if( err )
        {
            lis_vector_destroy(az);
            return err;
        }


        if( nrm>err_relation )
            err_relation = nrm;
    }


    lis_vector_destroy(az);


    *orth_error =
        err_orth;

    *relation_error =
        err_relation;


    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}



/*
 * Prepare a persistent flexible recycle basis for a changed matrix
 * value generation.
 *
 * Historical Z is authoritative:
 *
 *     C_new = A_new Z_old
 *
 * Historical directions are never reconstructed by applying the
 * current preconditioner.
 *
 * Quality follows the GCRO-DR cross-solve policy:
 *
 *     q_i = || C_old^T (A_new Z_i) ||_2
 *           / || A_new Z_i ||_2.
 *
 * KEEP      : qmin >= keep threshold
 * REJECT    : qavg < reject threshold, or fewer than half of the
 *             directions reach the reject threshold
 * TRANSPORT : otherwise
 *
 * KEEP and TRANSPORT both rebuild A_new Z = C and apply a
 * consistent QR transformation to W/Z/C.
 */
static LIS_INT
lis_fgcrodr_prepare_recycle(
    LIS_SOLVER solver)
{
    LIS_FGCRODR_STATE state;

    LIS_VECTOR *scratch_c = NULL;

    LIS_SCALAR alpha;

    LIS_REAL nrm;
    LIS_REAL proj2;
    LIS_REAL qi;
    LIS_REAL qsum;
    LIS_REAL qmin;
    LIS_REAL qavg;

    LIS_INT oldk;
    LIS_INT i;
    LIS_INT j;
    LIS_INT good;
    LIS_INT need;
    LIS_INT action;
    LIS_INT err;


    state =
        solver->fgcrodr_state;


    if(
        state==NULL ||
        !state->valid ||
        state->dim<=0 )
    {
        return LIS_SUCCESS;
    }


    state->cross_action =
        LIS_FGCRODR_RECYCLE_NONE;


    /*
     * Different matrix objects and structure generations remain
     * non-transportable in Stage 2E-A.
     */
    if(
        state->matrix!=solver->A ||
        state->structure_generation!=
            solver->A->structure_generation )
    {
        err =
            lis_fgcrodr_clear_recycle(
                solver);

        if( err )
            return err;


        state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->cross_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->rejects++;
        state->policy_rejects++;

        state->quality =
            0.0;

        state->quality_min =
            0.0;


        return LIS_SUCCESS;
    }


    /*
     * Exact generation: normal same-system reuse.
     */
    if(
        state->value_generation==
        solver->A->value_generation )
    {
        state->last_action =
            LIS_FGCRODR_RECYCLE_NONE;

        return LIS_SUCCESS;
    }


    oldk =
        state->dim;


    scratch_c =
        (LIS_VECTOR *)lis_malloc(
            (size_t)oldk
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_prepare_recycle::scratch_c");

    if( scratch_c==NULL )
    {
        LIS_SETERR_MEM(
            (size_t)oldk
            *sizeof(LIS_VECTOR));

        return LIS_ERR_OUT_OF_MEMORY;
    }


    for(i=0;i<oldk;i++)
        scratch_c[i] = NULL;


    for(i=0;i<oldk;i++)
    {
        err =
            lis_vector_duplicate(
                solver->A,
                &scratch_c[i]);

        if( err )
            goto cleanup;
    }


    qsum =
        0.0;

    qmin =
        1.0;

    good =
        0;


    /*
     * Read-only assessment against the OLD orthonormal C basis.
     */
    for(i=0;i<oldk;i++)
    {
        lis_matvec(
            solver->A,
            state->z[i],
            scratch_c[i]);


        lis_vector_nrm2(
            scratch_c[i],
            &nrm);


        if(
            nrm<=FGCRODR_BUILD_TINY ||
            nrm!=nrm )
        {
            qi =
                0.0;
        }
        else
        {
            proj2 =
                0.0;


            for(j=0;j<oldk;j++)
            {
                lis_vector_dot(
                    state->c[j],
                    scratch_c[i],
                    &alpha);


                proj2 +=
                    (LIS_REAL)(
                        fabs(alpha)
                        *fabs(alpha));
            }


            if( proj2<0.0 )
                proj2 = 0.0;


            qi =
                sqrt(proj2)
                /nrm;


            if( qi<0.0 )
                qi = 0.0;

            if( qi>1.0 )
                qi = 1.0;
        }


        qsum +=
            qi;


        if( qi<qmin )
            qmin = qi;


        if(
            qi>=
            state->reject_threshold )
        {
            good++;
        }
    }


    qavg =
        qsum
        /(LIS_REAL)oldk;


    need =
        _max(
            (LIS_INT)1,
            (oldk+1)/2);


    if(
        qmin>=
        state->keep_threshold )
    {
        action =
            LIS_FGCRODR_RECYCLE_KEEP;
    }
    else if(
        qavg<
            state->reject_threshold ||
        good<need )
    {
        action =
            LIS_FGCRODR_RECYCLE_REJECT;
    }
    else
    {
        action =
            LIS_FGCRODR_RECYCLE_TRANSPORT;
    }


    if(
        action==
        LIS_FGCRODR_RECYCLE_REJECT )
    {
        err =
            lis_fgcrodr_clear_recycle(
                solver);

        if( err )
            goto cleanup;


        state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->cross_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->rejects++;
        state->policy_rejects++;

        state->quality =
            qavg;

        state->quality_min =
            qmin;


        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    /*
     * Preserve measured quality before the atomic basis replacement.
     */
    state->quality =
        qavg;

    state->quality_min =
        qmin;

    state->cross_action =
        action;


    err =
        lis_fgcrodr_replace_recycle_basis(
            solver,
            state->w,
            state->z,
            oldk,
            LIS_FALSE,
            action,
            qavg,
            qmin);

    if( err )
        goto cleanup;


    state =
        solver->fgcrodr_state;


    if(
        !state->valid ||
        state->dim<=0 ||
        state->matrix!=solver->A ||
        state->structure_generation!=
            solver->A->structure_generation ||
        state->value_generation!=
            solver->A->value_generation )
    {
        err =
            lis_fgcrodr_clear_recycle(
                solver);

        if( err )
            goto cleanup;


        state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->cross_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->rejects++;
        state->policy_rejects++;

        state->quality =
            qavg;

        state->quality_min =
            qmin;

        state->cross_action =
            LIS_FGCRODR_RECYCLE_NONE;


        err =
            LIS_SUCCESS;

        goto cleanup;
    }


    state->cross_uses++;

    state->last_action =
        action;

    state->quality =
        qavg;

    state->quality_min =
        qmin;


    if(
        action==
        LIS_FGCRODR_RECYCLE_KEEP )
    {
        state->policy_keeps++;
    }
    else
    {
        state->policy_transports++;
    }


cleanup:

    if( scratch_c )
    {
        for(i=0;i<oldk;i++)
        {
            if( scratch_c[i] )
                lis_vector_destroy(
                    scratch_c[i]);
        }
    }


    lis_free(
        scratch_c);


    return err;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_flexible_core"
static LIS_INT
lis_fgcrodr_flexible_core(
    LIS_SOLVER solver)
{
    LIS_Comm comm;

    LIS_MATRIX A;

    LIS_VECTOR b;
    LIS_VECTOR x;

    LIS_VECTOR s;
    LIS_VECTOR *z;
    LIS_VECTOR *v;

    LIS_SCALAR *h;
    LIS_SCALAR *hraw;
    LIS_SCALAR *bcoef;

    LIS_SCALAR aa;
    LIS_SCALAR bb;
    LIS_SCALAR rr;

#ifdef _COMPLEX
    LIS_SCALAR a2;
    LIS_SCALAR b2;
#else
    LIS_REAL aa_abs;
    LIS_REAL bb_abs;
    LIS_REAL givens_scale;
    LIS_REAL givens_a;
    LIS_REAL givens_b;
#endif

    LIS_SCALAR t;
    LIS_SCALAR alpha;
    LIS_SCALAR gamma;

    LIS_REAL bnrm2;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL min_nrm2;
    LIS_REAL rnorm;

    LIS_INT err;
    LIS_INT iter;
    LIS_INT maxiter;
    LIS_INT output;
    LIS_INT maxiter_noimp;
    LIS_INT noimp_count;
    LIS_INT noimp_stop;

    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT m;

    LIS_INT ii;
    LIS_INT i1;
    LIS_INT iiv;
    LIS_INT i1v;
    LIS_INT iih;
    LIS_INT jj;

    LIS_INT h_dim;
    LIS_INT cs;
    LIS_INT sn;

    LIS_FGCRODR_STATE state;

    LIS_INT use_recycle;
    LIS_INT kcur;
    LIS_INT ri;

    double time;
    double ptime;


    comm = LIS_COMM_WORLD;

    A = solver->A;
    b = solver->b;
    x = solver->x;

    state =
        solver->fgcrodr_state;

    use_recycle =
        (
            state!=NULL &&
            state->valid &&
            state->dim>0 &&
            state->matrix==A &&
            state->structure_generation==
                A->structure_generation &&
            state->value_generation==
                A->value_generation
        )
        ? LIS_TRUE
        : LIS_FALSE;

    kcur =
        use_recycle
        ? state->dim
        : 0;

    hraw = NULL;
    bcoef = NULL;

    maxiter =
        solver->options[LIS_OPTIONS_MAXITER];

    output =
        solver->options[LIS_OPTIONS_OUTPUT];

    m =
        solver->options[LIS_OPTIONS_RESTART];

    maxiter_noimp =
        solver->options[LIS_OPTIONS_MAXITER_NO_IMP];

    noimp_count = 0;
    min_nrm2 = LIS_SCALAR_MAX;
    noimp_stop = LIS_FALSE;

    h_dim = m+1;

    ptime = 0.0;


    /*
     * Keep the Stage-2A workspace exactly compatible with FGMRES.
     *
     * work[0]       : small residual/Givens vector
     * work[2..]     : Z basis
     * work[m+2..]   : V basis
     */
    s = solver->work[0];
    z = &solver->work[2];
    v = &solver->work[m+2];


    h = (LIS_SCALAR *)lis_malloc(
        sizeof(LIS_SCALAR)
        *(h_dim+1)
        *(h_dim+2),
        "lis_fgcrodr_flexible_core::h");

    if( h==NULL )
    {
        LIS_SETERR_MEM(
            sizeof(LIS_SCALAR)
            *(h_dim+1)
            *(h_dim+2));

        return LIS_ERR_OUT_OF_MEMORY;
    }


    cs = (m+1)*h_dim;
    sn = (m+2)*h_dim;


    hraw =
        (LIS_SCALAR *)lis_malloc(
            (size_t)h_dim
            *(size_t)m
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_flexible_core::hraw");

    if( hraw==NULL )
    {
        if( hraw )
            lis_free(hraw);

        lis_free(h);

        LIS_SETERR_MEM(
            (size_t)h_dim
            *(size_t)m
            *sizeof(LIS_SCALAR));

        return LIS_ERR_OUT_OF_MEMORY;
    }


    memset(
        hraw,
        0,
        (size_t)h_dim
        *(size_t)m
        *sizeof(LIS_SCALAR));


    if( use_recycle )
    {
        bcoef =
            (LIS_SCALAR *)lis_malloc(
                (size_t)kcur
                *(size_t)m
                *sizeof(LIS_SCALAR),
                "lis_fgcrodr_flexible_core::bcoef");

        if( bcoef==NULL )
        {
            if( bcoef )
            lis_free(bcoef);

        if( hraw )
            lis_free(hraw);

        lis_free(h);

            LIS_SETERR_MEM(
                (size_t)kcur
                *(size_t)m
                *sizeof(LIS_SCALAR));

            return LIS_ERR_OUT_OF_MEMORY;
        }

        memset(
            bcoef,
            0,
            (size_t)kcur
            *(size_t)m
            *sizeof(LIS_SCALAR));
    }


    /* -------------------------------------------------------- */
    /* Initial residual                                         */
    /* -------------------------------------------------------- */

    if( lis_solver_get_initial_residual(
            solver,
            NULL,
            NULL,
            v[0],
            &bnrm2) )
    {
        if( bcoef )
            lis_free(bcoef);

        if( hraw )
            lis_free(hraw);

        lis_free(h);

        return LIS_SUCCESS;
    }


    tol = solver->tol;


    if( use_recycle )
    {
        for(ri=0;ri<kcur;ri++)
        {
            lis_vector_dot(
                state->c[ri],
                v[0],
                &gamma);

            lis_vector_axpy(
                gamma,
                state->z[ri],
                x);

            lis_vector_axpy(
                -gamma,
                state->c[ri],
                v[0]);
        }


        lis_vector_nrm2(
            v[0],
            &rnorm);


        if(
            rnorm==0.0 ||
            rnorm!=rnorm )
        {
            solver->retcode =
                (rnorm==0.0)
                ? LIS_SUCCESS
                : LIS_BREAKDOWN;

            solver->iter = 0;
            solver->resid = rnorm;
            solver->ptime = ptime;

            if( bcoef )
                lis_free(bcoef);

            if( hraw )
                lis_free(hraw);

            lis_free(h);

            return solver->retcode;
        }


        /*
         * Current FGCRO-DR core follows the FGMRES/default residual
         * normalization path: v0 is scaled by 1/||r||.
         */
        bnrm2 =
            1.0/rnorm;

        if(
            state->cross_action==
            LIS_FGCRODR_RECYCLE_NONE )
        {
            state->last_action =
                LIS_FGCRODR_RECYCLE_KEEP;
        }
        else
        {
            state->last_action =
                state->cross_action;
        }

        state->uses++;
    }
    else
    {
        rnorm =
            1.0/bnrm2;
    }


    iter = 0;


    while( iter<maxiter )
    {
        /*
         * First Arnoldi vector:
         *
         *     v_0 = r / ||r||
         */
        lis_vector_scale(
            bnrm2,
            v[0]);


        /*
         * Least-squares RHS:
         *
         *     s = ||r|| e_1
         */
        lis_vector_set_all(
            0,
            s);

        s->value[0] = rnorm;


        i = 0;


        do
        {
            iter++;
            i++;

            ii  = i-1;
            i1  = i;

            iiv = i-1;
            i1v = i;

            iih = (i-1)*h_dim;


            /* ------------------------------------------------ */
            /* FLEXIBLE STEP                                    */
            /*                                                  */
            /*     z_j = M_j^{-1} v_j                           */
            /* ------------------------------------------------ */

            time = lis_wtime();

            err = lis_psolve(
                solver,
                v[iiv],
                z[iiv]);

            ptime +=
                lis_wtime()-time;


            if( err )
            {
                solver->retcode = err;
                solver->iter = iter;
                solver->ptime = ptime;

                if( bcoef )
            lis_free(bcoef);

        if( hraw )
            lis_free(hraw);

        lis_free(h);

                return err;
            }


            /* ------------------------------------------------ */
            /* Operator image                                   */
            /*                                                  */
            /*     w = A z_j                                    */
            /*                                                  */
            /* Stage 2B will insert here:                       */
            /*                                                  */
            /*     B(:,j) = C^T w                               */
            /*     w      = w - C B(:,j)                        */
            /* ------------------------------------------------ */

            lis_matvec(
                A,
                z[iiv],
                v[i1v]);


            /*
             * Active FGCRO projection:
             *
             *     B(:,j) = C^T A z_j
             *     w      = A z_j - C B(:,j)
             */
            if( use_recycle )
            {
                for(ri=0;ri<kcur;ri++)
                {
                    lis_vector_dot(
                        state->c[ri],
                        v[i1v],
                        &alpha);

                    bcoef[
                        ri
                        + iiv*kcur] =
                        alpha;

                    lis_vector_axpy(
                        -alpha,
                        state->c[ri],
                        v[i1v]);
                }
            }


            /* ------------------------------------------------ */
            /* Arnoldi orthogonalization against V              */
            /* ------------------------------------------------ */

            for(k=0;k<i;k++)
            {
                lis_vector_dot(
                    v[i1v],
                    v[k],
                    &t);

                h[k+iih] = t;

                lis_vector_axpy(
                    -t,
                    v[k],
                    v[i1v]);
            }


            /*
             * h(i+1,i) = ||w||
             * v(i+1)   = w / h(i+1,i)
             */
            lis_vector_nrm2(
                v[i1v],
                (LIS_REAL *)&t);

            h[i1+iih] = t;


            /*
             * Preserve the unrotated flexible Arnoldi column before
             * Givens transformations modify h.
             */
            for(k=0;k<=i;k++)
            {
                hraw[
                    k
                    +iiv*h_dim] =
                    h[
                        k
                        +iih];
            }



            /*
             * Preserve a zero norm so the existing lucky-breakdown
             * behavior remains identical to FGMRES.
             */
            if( t!=0.0 && t==t )
            {
                lis_vector_scale(
                    1.0/t,
                    v[i1v]);
            }


            /* ------------------------------------------------ */
            /* Apply previous Givens rotations                  */
            /* ------------------------------------------------ */

            for(k=1;k<=ii;k++)
            {
                jj = k-1;

                t = h[jj+iih];

                aa =
                    h[jj+cs]*t;

                aa +=
                    h[jj+sn]
                    *h[k+iih];

                bb =
                    -h[jj+sn]*t;

                bb +=
                    h[jj+cs]
                    *h[k+iih];

                h[jj+iih] = aa;
                h[k+iih] = bb;
            }


            /* ------------------------------------------------ */
            /* New Givens rotation                              */
            /* ------------------------------------------------ */

            aa = h[ii+iih];
            bb = h[i1+iih];

#ifdef _COMPLEX
            /*
             * Complex recycling is rejected before this core is
             * entered.  Keep the existing expression here so the
             * complex build remains source-compatible.
             */
            a2 = aa*aa;
            b2 = bb*bb;

            rr = sqrt(a2+b2);
#else
            /*
             * Scale the Euclidean norm before squaring.
             *
             * Direct evaluation of
             *
             *     sqrt(aa*aa + bb*bb)
             *
             * can overflow when aa or bb is large and underflow when
             * both are small, even though the corresponding Givens
             * rotation is perfectly representable.
             *
             * After division by max(|aa|,|bb|), both squared terms
             * are bounded by one.
             */
            aa_abs =
                fabs(aa);

            bb_abs =
                fabs(bb);

            givens_scale =
                (aa_abs>bb_abs)
                ? aa_abs
                : bb_abs;

            if( givens_scale==0.0 )
            {
                rr =
                    0.0;
            }
            else
            {
                givens_a =
                    aa_abs
                    /givens_scale;

                givens_b =
                    bb_abs
                    /givens_scale;

                rr =
                    givens_scale
                    *sqrt(
                        givens_a*givens_a
                        +
                        givens_b*givens_b);
            }
#endif


            if( rr==0.0 || rr!=rr )
            {
                nrm2 =
                    fabs(s->value[ii]);

                solver->retcode =
                    LIS_BREAKDOWN;

                solver->iter =
                    iter;

                solver->resid =
                    nrm2;

                solver->ptime =
                    ptime;

                if( bcoef )
            lis_free(bcoef);

        if( hraw )
            lis_free(hraw);

        lis_free(h);

                return LIS_BREAKDOWN;
            }


            h[ii+cs] =
                aa/rr;

            h[ii+sn] =
                bb/rr;


            s->value[i1] =
                -h[ii+sn]
                *s->value[ii];

            s->value[ii] =
                h[ii+cs]
                *s->value[ii];


            aa =
                h[ii+cs]
                *h[ii+iih];

            aa +=
                h[ii+sn]
                *h[i1+iih];

            h[ii+iih] = aa;


            /* ------------------------------------------------ */
            /* Convergence                                      */
            /* ------------------------------------------------ */

            nrm2 =
                fabs(s->value[i1]);


            if( output )
            {
                if( output & LIS_PRINT_MEM )
                {
                    solver->rhistory[iter] =
                        nrm2;
                }

                if( output & LIS_PRINT_OUT )
                {
                    lis_print_rhistory(
                        comm,
                        iter,
                        nrm2);
                }
            }


            if( tol>=nrm2 )
                break;


            if( maxiter_noimp )
            {
                if( min_nrm2>nrm2 )
                {
                    min_nrm2 =
                        nrm2;

                    noimp_count =
                        0;
                }
                else if(
                    ++noimp_count >
                    maxiter_noimp )
                {
                    /*
                     * Preserve the approximation generated by the
                     * current Arnoldi cycle before terminating for
                     * lack of further residual improvement.
                     */
                    noimp_stop =
                        LIS_TRUE;

                    break;
                }
            }


        } while(
            i<m &&
            iter<maxiter );


        /* ---------------------------------------------------- */
        /* Solve upper-triangular H y = s                       */
        /* ---------------------------------------------------- */

        if(
            h[ii+iih]==0.0 ||
            h[ii+iih]!=h[ii+iih] )
        {
            solver->retcode =
                LIS_BREAKDOWN;

            solver->iter =
                iter;

            solver->resid =
                nrm2;

            solver->ptime =
                ptime;

            if( bcoef )
            lis_free(bcoef);

        if( hraw )
            lis_free(hraw);

        lis_free(h);

            return LIS_BREAKDOWN;
        }


        s->value[ii] =
            s->value[ii]
            /h[ii+iih];


        /*
         * A nonzero pivot can still produce an unrepresentable
         * reduced coefficient.  Reject it before applying Z*y to
         * the physical solution vector.
         */
        if(
            s->value[ii]!=s->value[ii] ||
            fabs(s->value[ii])>LIS_SCALAR_MAX )
        {
            solver->retcode =
                LIS_BREAKDOWN;

            solver->iter =
                iter;

            solver->resid =
                nrm2;

            solver->ptime =
                ptime;

            if( bcoef )
                lis_free(
                    bcoef);

            if( hraw )
                lis_free(
                    hraw);

            lis_free(
                h);

            return LIS_BREAKDOWN;
        }


        for(k=1;k<=ii;k++)
        {
            jj = ii-k;

            t =
                s->value[jj];


            for(j=jj+1;j<=ii;j++)
            {
                t -=
                    h[jj+j*h_dim]
                    *s->value[j];
            }


            if(
                h[jj+jj*h_dim]==0.0 ||
                h[jj+jj*h_dim]
                    !=h[jj+jj*h_dim] )
            {
                solver->retcode =
                    LIS_BREAKDOWN;

                solver->iter =
                    iter;

                solver->resid =
                    nrm2;

                solver->ptime =
                    ptime;

                if( bcoef )
            lis_free(bcoef);

        if( hraw )
            lis_free(hraw);

        lis_free(h);

                return LIS_BREAKDOWN;
            }


            s->value[jj] =
                t
                /h[jj+jj*h_dim];


            /*
             * Preserve x until the complete reduced solution is
             * representable.  A finite nonzero pivot alone does
             * not guarantee a finite quotient.
             */
            if(
                s->value[jj]!=s->value[jj] ||
                fabs(s->value[jj])>LIS_SCALAR_MAX )
            {
                solver->retcode =
                    LIS_BREAKDOWN;

                solver->iter =
                    iter;

                solver->resid =
                    nrm2;

                solver->ptime =
                    ptime;

                if( bcoef )
                    lis_free(
                        bcoef);

                if( hraw )
                    lis_free(
                        hraw);

                lis_free(
                    h);

                return LIS_BREAKDOWN;
            }
        }


        /* ---------------------------------------------------- */
        /* Flexible correction                                  */
        /*                                                      */
        /*     x = x_0 + Z y                                    */
        /*                                                      */
        /* Never substitute V for Z here.                       */
        /* ---------------------------------------------------- */

        for(j=0;j<=ii;j++)
        {
            lis_vector_axpy(
                s->value[j],
                z[j],
                x);
        }


        if( use_recycle )
        {
            for(ri=0;ri<kcur;ri++)
            {
                gamma = 0.0;

                for(j=0;j<=ii;j++)
                {
                    gamma +=
                        bcoef[
                            ri
                            + j*kcur]
                        *s->value[j];
                }

                lis_vector_axpy(
                    -gamma,
                    state->z[ri],
                    x);
            }
        }


        /*
         * Preserve the first usable flexible Arnoldi cycle before
         * restart overwrites V/Z.
         *
         * The first cycle can BUILD a harmonic recycle space.
         * Subsequent cycles with active recycling perform an augmented
         * harmonic UPDATE using [Z_old,Z_krylov].
         */
        err =
            lis_fgcrodr_build_recycle(
                solver,
                v,
                z,
                ii+1,
                hraw,
                h_dim,
                bcoef,
                kcur);

        if( err )
        {
            solver->retcode =
                err;

            solver->iter =
                iter;

            solver->resid =
                nrm2;

            solver->ptime =
                ptime;

            if( bcoef )
            lis_free(bcoef);

        if( hraw )
            lis_free(hraw);

        lis_free(h);

            return err;
        }


        /*
         * The no-improvement stop was requested during Arnoldi.
         * The current flexible correction and recycle-space
         * maintenance are now complete.
         */
        if( noimp_stop )
        {
            solver->retcode =
                LIS_MAXITER;

            solver->iter =
                iter;

            solver->resid =
                nrm2;

            solver->ptime =
                ptime;

            if( bcoef )
                lis_free(
                    bcoef);

            if( hraw )
                lis_free(
                    hraw);

            lis_free(
                h);

            return LIS_MAXITER;
        }


        /*
         * The projected FGCRO-DR residual is only a convergence
         * candidate.  Form the explicit residual of the updated
         * iterate before accepting success.
         *
         * Keep this residual in v[0] so a rejected candidate can
         * continue directly through the existing restart path.
         */
        lis_matvec(
            A,
            x,
            v[0]);

        lis_vector_xpay(
            b,
            -1.0,
            v[0]);


        if( tol>=nrm2 )
        {
            lis_solver_get_residual[0](
                v[0],
                solver,
                &nrm2);


            if(
                nrm2!=nrm2 ||
                fabs(nrm2)>LIS_SCALAR_MAX )
            {
                solver->retcode =
                    LIS_BREAKDOWN;

                solver->iter =
                    iter;

                solver->resid =
                    nrm2;

                solver->ptime =
                    ptime;

                if( bcoef )
                    lis_free(
                        bcoef);

                if( hraw )
                    lis_free(
                        hraw);

                lis_free(
                    h);

                return LIS_BREAKDOWN;
            }


            if( tol>=nrm2 )
            {
                solver->retcode =
                    LIS_SUCCESS;

                solver->iter =
                    iter;

                solver->resid =
                    nrm2;

                solver->ptime =
                    ptime;

                if( bcoef )
                    lis_free(
                        bcoef);

                if( hraw )
                    lis_free(
                        hraw);

                lis_free(
                    h);

                return LIS_SUCCESS;
            }
        }


        /* ---------------------------------------------------- */
        /* Restart residual                                     */
        /* ---------------------------------------------------- */


        if( use_recycle )
        {
            for(ri=0;ri<kcur;ri++)
            {
                lis_vector_dot(
                    state->c[ri],
                    v[0],
                    &gamma);

                lis_vector_axpy(
                    gamma,
                    state->z[ri],
                    x);

                lis_vector_axpy(
                    -gamma,
                    state->c[ri],
                    v[0]);
            }
        }


        lis_vector_nrm2(
            v[0],
            &rnorm);


        /*
         * The physical restart residual may vanish even when the
         * least-squares estimate which led to this restart was not
         * exactly zero.
         *
         * Do not form 1/||r|| for an exact restart solution.  Keep
         * the same zero/NaN contract as the initial recycled-residual
         * projection above.
         */
        if(
            rnorm==0.0 ||
            rnorm!=rnorm )
        {
            solver->retcode =
                (rnorm==0.0)
                ? LIS_SUCCESS
                : LIS_BREAKDOWN;

            solver->iter =
                iter;

            solver->resid =
                rnorm;

            solver->ptime =
                ptime;

            if( bcoef )
                lis_free(
                    bcoef);

            if( hraw )
                lis_free(
                    hraw);

            lis_free(
                h);


            return
                solver->retcode;
        }


        bnrm2 =
            1.0/rnorm;
    }


    solver->retcode =
        LIS_MAXITER;

    solver->iter =
        iter+1;

    solver->resid =
        nrm2;

    if( hraw )
        lis_free(hraw);

    lis_free(h);

    return LIS_MAXITER;
}


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr"
LIS_INT
lis_fgcrodr(
    LIS_SOLVER solver)
{
    LIS_INT err;
    LIS_INT solve_cross_action;

    LIS_FGCRODR_STATE state;

    LIS_DEBUG_FUNC_IN;


    err =
        lis_fgcrodr_state_ensure(
            solver);

    if( err )
        return err;


    state =
        solver->fgcrodr_state;


    solve_cross_action =
        LIS_FGCRODR_RECYCLE_NONE;


    /*
     * Keep k=0 as an exact call into the existing FGMRES
     * implementation.  It is our numerical reference path.
     */
    if( state->target_dim==0 )
    {
        LIS_DEBUG_FUNC_OUT;

        return lis_fgmres(
            solver);
    }


    /*
     * Make the persistent recycle relation compatible with the
     * current matrix generation before entering the flexible core.
     */
    err =
        lis_fgcrodr_prepare_recycle(
            solver);

    if( err )
        return err;


    solve_cross_action =
        state->cross_action;


    err =
        lis_fgcrodr_flexible_core(
            solver);


    /*
     * cross_action is transient for one solve.
     *
     * Use the local snapshot because a REJECT may be followed by a
     * fresh BUILD whose recycle initialization legitimately resets
     * state->cross_action.
     */
    if(
        solve_cross_action
        !=LIS_FGCRODR_RECYCLE_NONE )
    {
        state->last_action =
            solve_cross_action;
    }


    state->cross_action =
        LIS_FGCRODR_RECYCLE_NONE;


    LIS_DEBUG_FUNC_OUT;

    return err;
}


#ifdef USE_QUAD_PRECISION

#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_quad"
LIS_INT
lis_fgcrodr_quad(
    LIS_SOLVER solver)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_fgcrodr_state_ensure(solver);
    if( err )
        return err;

    if( solver->fgcrodr_state->target_dim==0 )
    {
        LIS_DEBUG_FUNC_OUT;

        return lis_fgmres_quad(solver);
    }

    LIS_SETERR(
        LIS_ERR_NOT_IMPLEMENTED,
        "quad FGCRO-DR recycling core is not implemented yet; "
        "use -fgcrodr_k 0 for the FGMRES reference path\n");

    return LIS_ERR_NOT_IMPLEMENTED;
}

#endif


#undef FGCRODR_DEFAULT_RECYCLE_DIM
#undef FGCRODR_DEFAULT_KEEP_QUALITY
#undef FGCRODR_DEFAULT_REJECT_QUALITY
#undef FGCRODR_BUILD_TINY
