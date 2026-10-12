/*
 * Flexible GCRO-DR.
 *
 * Stage 1B provides solver registration, persistent private recycle
 * state, configuration APIs and lifecycle handling.
 *
 * The actual flexible recycling iteration is added separately.
 */

#include "lislib.h"
#include <string.h>
#include <math.h>


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

    state->quality = 0.0;
    state->quality_min = 0.0;

    state->last_action =
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


#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_build_recycle"
static LIS_INT
lis_fgcrodr_build_recycle(
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
            "invalid FGCRO-DR recycle BUILD arguments\n");

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
     * Stage 2B performs only the first BUILD.
     *
     * Reuse/projection of an already valid C space is introduced in
     * Stage 2C.
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
     * Stage 2B uses the first independent flexible directions from
     * the final Arnoldi cycle.
     *
     * Harmonic-Ritz selection replaces this provisional choice in a
     * later stage.
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
    LIS_SCALAR *bcoef;

    LIS_SCALAR aa;
    LIS_SCALAR bb;
    LIS_SCALAR rr;
    LIS_SCALAR a2;
    LIS_SCALAR b2;
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

            lis_free(h);

            return solver->retcode;
        }


        /*
         * Current FGCRO-DR core follows the FGMRES/default residual
         * normalization path: v0 is scaled by 1/||r||.
         */
        bnrm2 =
            1.0/rnorm;

        state->last_action =
            LIS_FGCRODR_RECYCLE_KEEP;

        state->uses++;
        state->policy_keeps++;
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

            a2 = aa*aa;
            b2 = bb*bb;

            rr = sqrt(a2+b2);


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
                    solver->retcode =
                        LIS_MAXITER;

                    solver->iter =
                        iter;

                    solver->resid =
                        nrm2;

                    solver->ptime =
                        ptime;

                    if( bcoef )
            lis_free(bcoef);

        lis_free(h);

                    return LIS_MAXITER;
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

        lis_free(h);

            return LIS_BREAKDOWN;
        }


        s->value[ii] =
            s->value[ii]
            /h[ii+iih];


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

        lis_free(h);

                return LIS_BREAKDOWN;
            }


            s->value[jj] =
                t
                /h[jj+jj*h_dim];
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
         * lis_fgcrodr_build_recycle() is intentionally a no-op once
         * a valid initial recycle space exists, so later cycles do
         * not replace this Stage-2B BUILD.
         */
        err =
            lis_fgcrodr_build_recycle(
                solver,
                v,
                z,
                ii+1);

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

        lis_free(h);

            return err;
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
            lis_free(bcoef);

        lis_free(h);

            return LIS_SUCCESS;
        }


        /* ---------------------------------------------------- */
        /* Restart residual                                     */
        /* ---------------------------------------------------- */

        lis_matvec(
            A,
            x,
            v[0]);

        lis_vector_xpay(
            b,
            -1.0,
            v[0]);


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

        bnrm2 =
            1.0/rnorm;
    }


    solver->retcode =
        LIS_MAXITER;

    solver->iter =
        iter+1;

    solver->resid =
        nrm2;

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
    LIS_FGCRODR_STATE state;

    LIS_DEBUG_FUNC_IN;


    err =
        lis_fgcrodr_state_ensure(
            solver);

    if( err )
        return err;


    state =
        solver->fgcrodr_state;


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
     * Stage 2C consumes recycle vectors only for the exact operator
     * generation that created A*Z=C.
     *
     * Changed-operator transport is added in a later stage.
     */
    if(
        state->valid &&
        state->dim>0 &&
        (
            state->matrix!=solver->A ||
            state->structure_generation!=
                solver->A->structure_generation ||
            state->value_generation!=
                solver->A->value_generation
        )
      )
    {
        err =
            lis_fgcrodr_clear_recycle(
                solver);

        if( err )
            return err;

        state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        state->rejects++;
    }
    else
    {
        state->last_action =
            LIS_FGCRODR_RECYCLE_NONE;
    }


    err =
        lis_fgcrodr_flexible_core(
            solver);


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
