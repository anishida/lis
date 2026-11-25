/*
 * Experimental MFAMG -> residual -> native ILUTP preconditioner.
 * Global ILUTP is serial; explicit local RAS ILUTP uses OpenMP across domains.
 *
 * Both children are constructed exactly once and are owned by this
 * preconditioner.  No recursive Hybrid and no silent numerical fallback.
 * The outer Krylov solver always operates on the original matrix.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lislib.h"

struct LIS_MFAMG_ILUTP_STATE
{
    LIS_SOLVER mfamg;
    LIS_SOLVER ilutp;
    LIS_VECTOR rhs;
    LIS_VECTOR ax;
    LIS_VECTOR correction;
    LIS_MATRIX matrix; /* borrowed, same operator used for both factors */
    LIS_UNSIGNED_INT structure_generation;
    LIS_UNSIGNED_INT value_generation;
};

static void
lis_mfamg_ilutp_release_child(LIS_SOLVER child)
{
    if(child==NULL) return;
    if(child->precon!=NULL)
    {
        lis_precon_destroy(child->precon);
        child->precon=NULL;
    }
    lis_solver_destroy(child);
}

LIS_INT
lis_precon_destroy_mfamg_ilutp(LIS_PRECON precon)
{
    struct LIS_MFAMG_ILUTP_STATE *state;
    if(precon==NULL || precon->user_data==NULL) return LIS_SUCCESS;
    state=(struct LIS_MFAMG_ILUTP_STATE *)precon->user_data;
    precon->user_data=NULL;
    if(state->correction!=NULL) lis_vector_destroy(state->correction);
    if(state->ax!=NULL) lis_vector_destroy(state->ax);
    if(state->rhs!=NULL) lis_vector_destroy(state->rhs);
    lis_mfamg_ilutp_release_child(state->ilutp);
    lis_mfamg_ilutp_release_child(state->mfamg);
    lis_free(state);
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_mfamg_ilutp_create_child"
static LIS_INT
lis_mfamg_ilutp_create_child(LIS_SOLVER outer, LIS_INT kind,
                              LIS_SOLVER *created)
{
    LIS_SOLVER child=NULL;
    LIS_INT err;
    *created=NULL;
    err=lis_solver_create(&child);
    if(err) return err;

    /*
     * Copy only numerical configuration, never solver-owned state.
     * Near-nullspace vectors are deep-copied by their public setter.
     */
    memcpy(child->options,outer->options,sizeof(child->options));
    memcpy(child->params,outer->params,sizeof(child->params));
    child->options[LIS_OPTIONS_PRECON]=kind;
    child->options[LIS_OPTIONS_ADDS]=LIS_FALSE;
    child->options[LIS_OPTIONS_SCALE]=LIS_SCALE_NONE;
    child->A=outer->A;
    child->precision=outer->precision;
    child->mfamg_block_size=outer->mfamg_block_size;
    child->mfamg_fine_agg_factor=outer->mfamg_fine_agg_factor;
    child->mfamg_vcycle_smooth_steps=outer->mfamg_vcycle_smooth_steps;
    child->mfamg_deep_sibling_merge=outer->mfamg_deep_sibling_merge;
    child->precon_half_mode=(kind==LIS_PRECON_TYPE_MFAMG) ?
                            outer->precon_half_mode : LIS_PRECON_HALF_OFF;

    if(kind==LIS_PRECON_TYPE_MFAMG)
    {
        err=lis_solver_set_near_nullspace(child,
                                         lis_solver_amg_basis_dim(outer),
                                         lis_solver_amg_basis_vectors(outer));
        if(err) goto fail;

        /* CA-MFAMG is a property of the MFAMG child, not the ILUTP child.
         * Keep partition metadata intact when composing the two native
         * preconditioners. The setter owns a private copy of the labels.
         * Never guess labels from node numbers, fixed DOF counts or matrices. */
        if(outer->mfamg_contact_aware)
        {
            if(outer->mfamg_contact_partition==NULL ||
               outer->mfamg_contact_partition_n!=outer->A->n)
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                           "CA-MFAMG_ILUTP requires partition labels for every local DOF\n");
                err=LIS_ERR_ILL_ARG;
                goto fail;
            }
            child->mfamg_contact_aware=1;
            err=lis_solver_set_mfamg_contact_partition(
                child,outer->mfamg_contact_partition_n,
                outer->mfamg_contact_partition);
            if(err) goto fail;
        }
    }
    else
    {
        /* The composite always uses pivoted native ILUTP, not legacy ILUT. */
        child->ilutp_enabled=1;
        child->ilutp_threshold=outer->ilutp_threshold;
        child->ilutp_full=outer->ilutp_full;
        child->ilutp_full_limit=outer->ilutp_full_limit;
        child->ilut_rcm=outer->ilut_rcm;
        child->ilut_row_scale=outer->ilut_row_scale;
        child->ilut_equil=outer->ilut_equil;
        /* Explicit RAS mode gives each domain its own native pivoted LU.
         * Global ILUTP (overlap=-1) retains the original serial policy. */
        child->ilut_ras_overlap=outer->ilut_ras_overlap;
        child->ilut_ras_domains=outer->ilut_ras_domains;
        child->ilut_schur=outer->ilut_schur;
    }

    /* Metadata setters invalidate retained sequence state, including A.
     * Restore the borrowed operator only after all child setters finish. */
    child->A=outer->A;

    err=lis_precon_create(child,&child->precon);
    if(err)
    {
        /* lis_precon_create() frees its output on failure, but leaves it
         * non-NULL; do not try to free that dangling pointer again. */
        child->precon=NULL;
        goto fail;
    }
    *created=child;
    return LIS_SUCCESS;

fail:
    lis_mfamg_ilutp_release_child(child);
    return err;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_mfamg_ilutp"
LIS_INT
lis_precon_create_mfamg_ilutp(LIS_SOLVER solver,LIS_PRECON precon)
{
    LIS_MATRIX A;
    struct LIS_MFAMG_ILUTP_STATE *state;
    LIS_INT err;
    if(solver==NULL || precon==NULL || solver->A==NULL)
        return LIS_ERR_ILL_ARG;
    A=solver->A;

    if(A->matrix_type!=LIS_MATRIX_CSR || A->nprocs!=1 ||
       A->n!=A->gn || A->is_splited)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                   "MFAMG_ILUTP requires an unsplit, serial CSR matrix\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
#ifdef _COMPLEX
    LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
               "MFAMG_ILUTP requires a real-double build\n");
    return LIS_ERR_NOT_IMPLEMENTED;
#endif
    if(solver->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE ||
       solver->options[LIS_OPTIONS_SCALE]!=LIS_SCALE_NONE ||
       solver->options[LIS_OPTIONS_ADDS])
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                   "MFAMG_ILUTP requires double, scale none, no Additive Schwarz\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
#ifdef _OPENMP
    if(solver->ilut_ras_overlap<0 && omp_get_max_threads()!=1)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                   "MFAMG_ILUTP global ILUTP requires OMP_NUM_THREADS=1; "
                   "enable -ilut_ras_overlap for domain-parallel ILUTP\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
#endif
    /*
     * Threshold zero selects the legacy ILUT factorization in limited-fill
     * mode, contradicting the native-pivoting contract of this composite.
     * Never silently downgrade a requested MFAMG+ILUTP configuration.
     */
    if(!isfinite((double)solver->ilutp_threshold) ||
       solver->ilutp_threshold<=0.0 || solver->ilutp_threshold>1.0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "MFAMG_ILUTP requires 0 < ilutp_threshold <= 1\n");
        return LIS_ERR_ILL_ARG;
    }
    /* The ILUTP child owns its optional approximate interface Schur
     * correction, applied to the residual remaining after MFAMG. */
    if(solver->ilut_schur && solver->ilut_ras_overlap<0)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                   "MFAMG_ILUTP Schur correction requires explicit RAS\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
    if(solver->ilut_ras_overlap>=0)
    {
        if(solver->ilut_ras_overlap>2 ||
           solver->ilut_ras_domains<2 || solver->ilut_ras_domains>16 ||
           solver->ilut_ras_domains>A->n)
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                       "MFAMG_ILUTP RAS requires overlap 0..2 and 2..16 "
                       "domains, no more than the number of rows\n");
            return LIS_ERR_ILL_ARG;
        }
    }
    else if(solver->ilut_ras_overlap!=-1 || solver->ilut_ras_domains!=2)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "MFAMG_ILUTP RAS domains require an enabled RAS overlap\n");
        return LIS_ERR_ILL_ARG;
    }
    if(lis_solver_amg_basis_dim(solver)<1 || lis_solver_amg_basis_vectors(solver)==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "MFAMG_ILUTP requires supplied interpolation basis vectors\n");
        return LIS_ERR_ILL_ARG;
    }

    state=(struct LIS_MFAMG_ILUTP_STATE *)
        lis_malloc(sizeof(*state),"lis_precon_create_mfamg_ilutp::state");
    if(state==NULL)
    {
        LIS_SETERR_MEM(sizeof(*state));
        return LIS_OUT_OF_MEMORY;
    }
    memset(state,0,sizeof(*state));
    /* Built-in type 13 owns user_data; the registered USERDEF types
     * retain their existing user_data/user_destroy semantics. */
    precon->user_data=state;

    err=lis_mfamg_ilutp_create_child(solver,LIS_PRECON_TYPE_MFAMG,
                                     &state->mfamg);
    if(err) return err;
    err=lis_mfamg_ilutp_create_child(solver,LIS_PRECON_TYPE_ILUT,
                                     &state->ilutp);
    if(err) return err;

    err=lis_vector_duplicate(A,&state->rhs);
    if(err) return err;
    err=lis_vector_duplicate(A,&state->ax);
    if(err) return err;
    err=lis_vector_duplicate(A,&state->correction);
    if(err) return err;
    state->matrix=A; /* borrowed */
    err=lis_matrix_get_generation(A,&state->structure_generation,
                                  &state->value_generation);
    if(err) return err;
    precon->A=A; /* borrowed */
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_mfamg_ilutp"
LIS_INT
lis_psolve_mfamg_ilutp(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    struct LIS_MFAMG_ILUTP_STATE *state;
    LIS_UNSIGNED_INT structure_generation,value_generation;
    LIS_INT err;
    if(solver==NULL || solver->precon==NULL || b==NULL || x==NULL)
        return LIS_ERR_ILL_ARG;
    state=(struct LIS_MFAMG_ILUTP_STATE *)solver->precon->user_data;
    if(state==NULL || state->mfamg==NULL || state->ilutp==NULL ||
       state->rhs==NULL || state->ax==NULL || state->correction==NULL)
        return LIS_ERR_ILL_ARG;

    /* The setup matrix is P. The physical A may differ only when this
     * object is the sequence-owned preconditioner explicitly built on P.
     * For standalone/caller-owned applications retain the old identity
     * check. In either case verify the generation of P, not A. */
    if(solver->A!=state->matrix &&
       !(solver->sequence_precon==solver->precon &&
         solver->sequence_precon_matrix==state->matrix))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "MFAMG_ILUTP operator no longer matches setup matrix\n");
        return LIS_ERR_ILL_ARG;
    }
    err=lis_matrix_get_generation(state->matrix,
                                  &structure_generation,&value_generation);
    if(err) return err;
    if(structure_generation!=state->structure_generation ||
       value_generation!=state->value_generation)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "MFAMG_ILUTP preconditioner is stale; rebuild after matrix change\n");
        return LIS_ERR_ILL_ARG;
    }

    /* Save b first, including when caller aliases input and output. */
    err=lis_vector_copy(b,state->rhs);
    if(err) return err;
    /* Use the saved input so b==x is also safe. */
    err=lis_psolve_mfamg(state->mfamg,state->rhs,x);
    if(err) return err;

    err=lis_matvec(state->matrix,x,state->ax);
    if(err) return err;
    err=lis_vector_axpy((LIS_SCALAR)-1.0,state->ax,state->rhs);
    if(err) return err;

    /* Native ILUTP dispatches to domain-parallel RAS when explicitly
     * requested, or to the single global factor in legacy serial mode. */
    err=lis_psolve_ilut_csr(state->ilutp,state->rhs,state->correction);
    if(err) return err;
    return lis_vector_axpy((LIS_SCALAR)1.0,state->correction,x);
}
