/* Generated from pinned native lis_krylov.c; Copyright (C) 2026, contributors to LIS.
 * Redistribution and use are governed by the LIS project BSD license. */
struct frame_apply_metric { int pc,parent;
    struct krylov_workspace * s;
    const double * x;
    double * y;
    const lis_krylov_problem64 * p;
    lis_krylov_matvec64 fn;
};
struct frame_verify_inverse { int pc,parent;
    struct krylov_workspace * s;
    const double * rhs;
    const double * solution;
    const lis_krylov_problem64 * p;
    double numerator;
    double denominator;
    int64_t i;
};
struct frame_apply_operator { int pc,parent;
    struct krylov_workspace * s;
    const double * x;
    const double * metric_x;
    double * y;
    const lis_krylov_problem64 * p;
    int err;
};
struct frame_orthogonalize { int pc,parent;
    struct krylov_workspace * s;
    int count;
    double * x;
    double * metric_x;
    double initial;
    double square;
    double norm;
    int pass;
    int j;
    int err;
    double previous_square;
    const double * v;
    const double * bv;
    double coefficient;
    double inverse_norm;
};
struct frame_append_range_seed { int pc,parent;
    struct krylov_workspace * s;
    int count;
    double * x;
    int err;
};
struct frame_append_vector { int pc,parent;
    struct krylov_workspace * s;
    int count;
    double * x;
    int err;
};
struct frame_append_extra_seeds { int pc,parent;
    struct krylov_workspace * s;
    int * count;
    int capacity;
    int additional;
    unsigned start_seed;
    int accepted;
    int extra;
    int trial;
    int found;
    int rc;
};
struct frame_eigenpair_residual { int pc,parent;
    struct krylov_workspace * s;
    int m;
    int index;
    double mu;
    double * lambda;
    double * relative;
    const lis_krylov_problem64 * p;
    double norm_a;
    double norm_b;
    double norm_r;
    double norm_x;
    double denom;
};
struct frame_lis_krylov_symmetric_core { int pc,parent;
    lis_krylov_workspace * reuse;
    const lis_krylov_problem64 * p;
    const lis_krylov_options * o;
    const lis_krylov_ext_options * ext;
    double * values;
    double * vectors;
    double * residuals;
    lis_krylov_stats * stats;
    lis_krylov_ext_stats * ext_stats;
    struct krylov_workspace s;
    int m;
    int previous;
    int iter;
    int ncv;
    int capacity;
    int nev;
    int maxiter;
    int rc;
    int restart_directions;
    int last_nconv;
    int stalled_restarts;
    int i;
    int j;
    double tol;
    unsigned seed;
    const char * which;
    int output_sort;
    int both_ends;
    int be_high_locked;
    int be_low_locked;
    size_t which_len;
    double norm_x;
    double norm_a;
    int attempt;
    int old_m;
    int extensions;
    int reuse_verified;
    double norm_a_before;
    int captured;
    int idx;
    double lambda;
    double rel;
    int err;
    int slot;
    int trial;
    int ok;
    int added;
    int from_cache;
    int cache_slot;
    int increment;
    int remaining;
    int keep;
    int result;
    int extra;
    double mu;
    double norm_w;
    double norm_v;
    double norm_r;
    double threshold;
    double left;
    double right;
    double mu_left;
    double mu_right;
    double tmp;
    int64_t row;
};
struct lis_rci_krylov_restore {
 int active,retval,status,native_code; uint64_t token; lis_rci_krylov_restore_request request; lis_krylov_budget budget; lis_krylov_budget_stats counts;
 struct frame_apply_metric f_apply_metric;
 struct frame_verify_inverse f_verify_inverse;
 struct frame_apply_operator f_apply_operator;
 struct frame_orthogonalize f_orthogonalize;
 struct frame_append_range_seed f_append_range_seed;
 struct frame_append_vector f_append_vector;
 struct frame_append_extra_seeds f_append_extra_seeds;
 struct frame_eigenpair_residual f_eigenpair_residual;
 struct frame_lis_krylov_symmetric_core f_lis_krylov_symmetric_core;
};
#define K_STEP_AGAIN 1
#define K_STEP_YIELD 2
#define K_STEP_STOP 3
static int k_return(lis_rci_krylov_restore *,int,int);
static int k_suspend(lis_rci_krylov_restore *,uint32_t,int64_t,const double *,double *);
static int step_apply_metric(lis_rci_krylov_restore *K)
{
 struct frame_apply_metric *F=&K->f_apply_metric;
 switch(F->pc){case 0:break;case 1:goto resume_1;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    F->p = F->s->problem;
    ;
    if (F->p->transform == LIS_KRYLOV_STANDARD) {
        memcpy(F->y, F->x, (size_t)F->s->n * sizeof(double));
    } else {
        F->fn = F->p->transform == LIS_KRYLOV_SHIFT_INVERT_A ?
             F->p->apply_a : F->p->apply_b;
        KRYLOV_INC(F->s->stats->metric_count);
        F->pc = 1; return k_suspend(K, (F->p->transform == LIS_KRYLOV_SHIFT_INVERT_A ? LIS_RCI_OP_A : LIS_RCI_OP_B), F->s->n, F->x, F->y);
        resume_1: F->pc = 0;
        if (K->retval) return k_return(K,F->parent,-1);
    }
    return k_return(K,F->parent,all_finite(F->s->n, F->y) ? 0 : -2);

}
static int step_verify_inverse(lis_rci_krylov_restore *K)
{
 struct frame_verify_inverse *F=&K->f_verify_inverse;
 switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    F->p = F->s->problem;
    ;
    ;

    if (F->s->options->inner_tol == 0.0) return k_return(K,F->parent,0);
    F->pc = 1; return k_suspend(K, LIS_RCI_OP_A, F->s->n, F->solution, F->s->check);
    resume_1: F->pc = 0;
    if (K->retval) return k_return(K,F->parent,-1);
    F->pc = 2; return k_suspend(K, LIS_RCI_OP_B, F->s->n, F->solution, F->s->rhs);
    resume_2: F->pc = 0;
    if (K->retval) return k_return(K,F->parent,-1);
    for (F->i = 0; F->i < F->s->n; ++F->i)
        F->s->check[F->i] -= F->p->sigma * F->s->rhs[F->i] + F->rhs[F->i];
    F->numerator = stable_norm(F->s->n, F->s->check);
    F->denominator = stable_norm(F->s->n, F->rhs);
    if (!isfinite(F->numerator) || !isfinite(F->denominator)) return k_return(K,F->parent,-2);
    if (F->numerator > F->s->options->inner_tol * F->denominator)
        return k_return(K,F->parent,-3);
    return k_return(K,F->parent,0);

}
static int step_apply_operator(lis_rci_krylov_restore *K)
{
 struct frame_apply_operator *F=&K->f_apply_operator;
 switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    F->p = F->s->problem;
    ;
    KRYLOV_INC(F->s->stats->op_count);
    if (F->p->transform == LIS_KRYLOV_STANDARD) {
        F->pc = 1; return k_suspend(K, LIS_RCI_OP_A, F->s->n, F->x, F->y);
        resume_1: F->pc = 0;
        F->err = K->retval;
    } else {
        KRYLOV_INC(F->s->stats->solve_count);
        F->pc = 2; return k_suspend(K, LIS_RCI_OP_SHIFTED_SOLVE, F->s->n, F->metric_x, F->y);
        resume_2: F->pc = 0;
        F->err = K->retval;
        if (!F->err && all_finite(F->s->n, F->y)) {
            K->f_verify_inverse.s = F->s;
            K->f_verify_inverse.rhs = F->metric_x;
            K->f_verify_inverse.solution = F->y;
            K->f_verify_inverse.pc = 0; K->f_verify_inverse.parent = 3;
            F->pc = 3; K->active = 2; return K_STEP_AGAIN;
            resume_3: F->pc = 0;
            F->err = K->retval;
        }
    }
    if (F->err) return k_return(K,F->parent,F->err == -3 ? -3 : -1);
    return k_return(K,F->parent,all_finite(F->s->n, F->y) ? 0 : -2);

}
static int step_orthogonalize(lis_rci_krylov_restore *K)
{
 struct frame_orthogonalize *F=&K->f_orthogonalize;
 switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    ;
    ;

    F->norm = stable_norm(F->s->n, F->x);
    if (!isfinite(F->norm)) return k_return(K,F->parent,-2);
    if (F->norm == 0.0) return k_return(K,F->parent,1);
    vector_scale(F->s->n, 1.0 / F->norm, F->x);

    K->f_apply_metric.s = F->s;
    K->f_apply_metric.x = F->x;
    K->f_apply_metric.y = F->metric_x;
    K->f_apply_metric.pc = 0; K->f_apply_metric.parent = 4;
    F->pc = 1; K->active = 1; return K_STEP_AGAIN;
    resume_1: F->pc = 0;
    F->err = K->retval;
    if (F->err) return k_return(K,F->parent,F->err);
    F->initial = compensated_dot(F->s->n, F->x, F->metric_x);
    if (!isfinite(F->initial)) return k_return(K,F->parent,-2);
    if (F->initial <= 0.0) {
        if (F->s->options->allow_semidefinite_metric && F->initial == 0.0)
            return k_return(K,F->parent,1);
        return k_return(K,F->parent,-2);
    }





    {
        F->previous_square = F->initial;
        for (F->pass = 0; F->pass < 3; ++F->pass) {
            for (F->j = 0; F->j < F->s->nlocked; ++F->j) {
                F->v = F->s->locked_basis + (size_t)F->j * F->s->n;
                F->bv = F->s->locked_metric + (size_t)F->j * F->s->n;
                F->coefficient = compensated_dot(F->s->n, F->bv, F->x);
                if (!isfinite(F->coefficient)) return k_return(K,F->parent,-2);
                vector_axpy(F->s->n, -F->coefficient, F->v, F->x);
            }
            for (F->j = 0; F->j < F->count; ++F->j) {
                F->v = F->s->basis + (size_t)F->j * F->s->n;
                F->bv = F->s->metric_basis + (size_t)F->j * F->s->n;
                F->coefficient = compensated_dot(F->s->n, F->bv, F->x);
                if (!isfinite(F->coefficient)) return k_return(K,F->parent,-2);
                vector_axpy(F->s->n, -F->coefficient, F->v, F->x);
            }
            K->f_apply_metric.s = F->s;
            K->f_apply_metric.x = F->x;
            K->f_apply_metric.y = F->metric_x;
            K->f_apply_metric.pc = 0; K->f_apply_metric.parent = 4;
            F->pc = 2; K->active = 1; return K_STEP_AGAIN;
            resume_2: F->pc = 0;
            F->err = K->retval;
            if (F->err) return k_return(K,F->parent,F->err);
            F->square = compensated_dot(F->s->n, F->x, F->metric_x);
            if (!isfinite(F->square)) return k_return(K,F->parent,-2);

            if (F->square <= F->initial * KRYLOV_ORTH_REJECT *
                          KRYLOV_ORTH_REJECT) {
                if (F->square >= 0.0 ||
                    (F->s->options->allow_semidefinite_metric &&
                     F->square >= -64.0 * DBL_EPSILON * F->initial))
                    return k_return(K,F->parent,1);
                return k_return(K,F->parent,-2);
            }
            if (F->square >= KRYLOV_DGKS_ETA * KRYLOV_DGKS_ETA *
                          F->previous_square)
                break;
            F->previous_square = F->square;
        }
    }
    {
        F->inverse_norm = 1.0 / sqrt(F->square);
        vector_scale_pair(F->s->n, F->inverse_norm, F->x, F->metric_x);
    }
    for (F->j = 0; F->j < F->s->nlocked; ++F->j) {
        F->bv = F->s->locked_metric + (size_t)F->j * F->s->n;
        if (fabs(compensated_dot(F->s->n, F->bv, F->x)) > 1e-8)
            return k_return(K,F->parent,1);
    }
    for (F->j = 0; F->j < F->count; ++F->j) {
        F->bv = F->s->metric_basis + (size_t)F->j * F->s->n;
        if (fabs(compensated_dot(F->s->n, F->bv, F->x)) > 1.0e-8)
            return k_return(K,F->parent,1);
    }
    return k_return(K,F->parent,0);

}
static int step_append_range_seed(lis_rci_krylov_restore *K)
{
 struct frame_append_range_seed *F=&K->f_append_range_seed;
 switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    ;
    if (F->s->options->allow_semidefinite_metric) {
        K->f_apply_metric.s = F->s;
        K->f_apply_metric.x = F->x;
        K->f_apply_metric.y = F->s->metric_x;
        K->f_apply_metric.pc = 0; K->f_apply_metric.parent = 5;
        F->pc = 1; K->active = 1; return K_STEP_AGAIN;
        resume_1: F->pc = 0;
        F->err = K->retval;
        if (F->err) return k_return(K,F->parent,F->err);
        if (stable_norm(F->s->n, F->s->metric_x) == 0.0) return k_return(K,F->parent,1);


        K->f_apply_operator.s = F->s;
        K->f_apply_operator.x = F->x;
        K->f_apply_operator.metric_x = F->s->metric_x;
        K->f_apply_operator.y = F->s->scratch;
        K->f_apply_operator.pc = 0; K->f_apply_operator.parent = 5;
        F->pc = 2; K->active = 3; return K_STEP_AGAIN;
        resume_2: F->pc = 0;
        F->err = K->retval;
        if (F->err) return k_return(K,F->parent,F->err);
        memcpy(F->x, F->s->scratch, (size_t)F->s->n * sizeof(double));
    }
    K->f_append_vector.s = F->s;
    K->f_append_vector.count = F->count;
    K->f_append_vector.x = F->x;
    K->f_append_vector.pc = 0; K->f_append_vector.parent = 5;
    F->pc = 3; K->active = 6; return K_STEP_AGAIN;
    resume_3: F->pc = 0;
    return k_return(K,F->parent,K->retval);

}
static int step_append_vector(lis_rci_krylov_restore *K)
{
 struct frame_append_vector *F=&K->f_append_vector;
 switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    K->f_orthogonalize.s = F->s;
    K->f_orthogonalize.count = F->count;
    K->f_orthogonalize.x = F->x;
    K->f_orthogonalize.metric_x = F->s->metric_x;
    K->f_orthogonalize.pc = 0; K->f_orthogonalize.parent = 6;
    F->pc = 1; K->active = 4; return K_STEP_AGAIN;
    resume_1: F->pc = 0;
    F->err = K->retval;
    if (F->err) return k_return(K,F->parent,F->err);
    memcpy(F->s->basis + (size_t)F->count * F->s->n, F->x,
           (size_t)F->s->n * sizeof(double));
    if (F->s->metric_basis != F->s->basis)
        memcpy(F->s->metric_basis + (size_t)F->count * F->s->n, F->s->metric_x,
               (size_t)F->s->n * sizeof(double));
    K->f_apply_operator.s = F->s;
    K->f_apply_operator.x = F->x;
    K->f_apply_operator.metric_x = F->s->metric_x;
    K->f_apply_operator.y = F->s->images + (size_t)F->count * F->s->n;
    K->f_apply_operator.pc = 0; K->f_apply_operator.parent = 6;
    F->pc = 2; K->active = 3; return K_STEP_AGAIN;
    resume_2: F->pc = 0;
    F->err = K->retval;
    if (!F->err && F->s->ext_stats) KRYLOV_INC(F->s->ext_stats->basis_vectors);
    return k_return(K,F->parent,F->err);

}
static int step_append_extra_seeds(lis_rci_krylov_restore *K)
{
 struct frame_append_extra_seeds *F=&K->f_append_extra_seeds;
 switch(F->pc){case 0:break;case 1:goto resume_1;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    F->accepted = 0;
    for (F->extra = 0; F->extra < F->additional && *F->count < F->capacity; ++F->extra) {
        F->found = 0;
        for (F->trial = 0; F->trial < 12; ++F->trial) {
            ;
            seed_vector(F->s->n, F->s->x, F->start_seed +
                        (unsigned)F->extra * 41u + (unsigned)F->trial);
            if (F->s->ext_stats) KRYLOV_INC(F->s->ext_stats->extra_seed_attempts);
            K->f_append_range_seed.s = F->s;
            K->f_append_range_seed.count = *F->count;
            K->f_append_range_seed.x = F->s->x;
            K->f_append_range_seed.pc = 0; K->f_append_range_seed.parent = 7;
            F->pc = 1; K->active = 5; return K_STEP_AGAIN;
            resume_1: F->pc = 0;
            F->rc = K->retval;
            if (F->rc < 0) return k_return(K,F->parent,F->rc);
            if (F->rc == 0) {
                ++*F->count;
                ++F->accepted;
                if (F->s->ext_stats) KRYLOV_INC(F->s->ext_stats->extra_seed_accepted);
                F->found = 1;
                break;
            }
        }
        if (!F->found) break;
    }
    return k_return(K,F->parent,F->accepted);

}
static int step_eigenpair_residual(lis_rci_krylov_restore *K)
{
 struct frame_eigenpair_residual *F=&K->f_eigenpair_residual;
 switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    F->p = F->s->problem;
    ;

    rotate_column(F->s, F->s->basis, F->m, F->index, F->s->x);
    if (F->p->transform == LIS_KRYLOV_STANDARD) {
        *F->lambda = F->mu;
    } else if (F->p->transform == LIS_KRYLOV_SHIFT_INVERT_B) {
        if (fabs(F->mu) < DBL_MIN) return k_return(K,F->parent,1);
        *F->lambda = F->p->sigma + 1.0 / F->mu;
    } else {


        if (fabs(F->mu - 1.0) <= 64.0 * DBL_EPSILON *
                                fmax(1.0, fabs(F->mu))) return k_return(K,F->parent,1);
        *F->lambda = F->p->sigma * F->mu / (F->mu - 1.0);
    }
    if (!isfinite(*F->lambda)) return k_return(K,F->parent,1);

    F->pc = 1; return k_suspend(K, LIS_RCI_OP_A, F->s->n, F->s->x, F->s->metric_x);
    resume_1: F->pc = 0;
    if (K->retval) return k_return(K,F->parent,-1);
    F->norm_a = stable_norm(F->s->n, F->s->metric_x);
    if (F->p->transform == LIS_KRYLOV_STANDARD) {
        memcpy(F->s->rhs, F->s->x, (size_t)F->s->n * sizeof(double));
    } else {
        F->pc = 2; return k_suspend(K, LIS_RCI_OP_B, F->s->n, F->s->x, F->s->rhs);
        resume_2: F->pc = 0;
        if (K->retval) return k_return(K,F->parent,-1);
    }
    F->norm_b = stable_norm(F->s->n, F->s->rhs);
    F->norm_x = stable_norm(F->s->n, F->s->x);
    if (!isfinite(F->norm_a) || !isfinite(F->norm_b) ||
        !isfinite(F->norm_x)) return k_return(K,F->parent,-2);

    if (F->s->options->residual_scale == 0.0 && F->norm_x > 0.0 &&
        F->norm_a / F->norm_x > F->s->norm_a)
        F->s->norm_a = F->norm_a / F->norm_x;

    vector_axpy(F->s->n, -*F->lambda, F->s->rhs, F->s->metric_x);
    F->norm_r = stable_norm(F->s->n, F->s->metric_x);
    F->denom = F->norm_a + fabs(*F->lambda) * F->norm_b;
    if (F->s->norm_a * F->norm_x > F->denom)
        F->denom = F->s->norm_a * F->norm_x;
    if (!isfinite(F->denom) || !isfinite(F->norm_r)) return k_return(K,F->parent,-2);
    *F->relative = F->denom > 0.0 ? F->norm_r / F->denom : F->norm_r;
    return k_return(K,F->parent,isfinite(*F->relative) ? 0 : -2);

}
static int step_lis_krylov_symmetric_core(lis_rci_krylov_restore *K)
{
 struct frame_lis_krylov_symmetric_core *F=&K->f_lis_krylov_symmetric_core;
 switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;case 4:goto resume_4;case 5:goto resume_5;case 6:goto resume_6;case 7:goto resume_7;case 8:goto resume_8;case 9:goto resume_9;case 10:goto resume_10;case 11:goto resume_11;case 12:goto resume_12;default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}

    F->s = (struct krylov_workspace){0};
    F->m = 0;
    F->previous = 0;
    F->restart_directions = 1;
    F->last_nconv = 0;
    F->stalled_restarts = 0;
    ;
    ;
    ;
    ;
    F->output_sort = 0;
    F->both_ends = 0;
    F->be_high_locked = 0;
    F->be_low_locked = 0;

    if (F->stats) memset(F->stats, 0, sizeof(*F->stats));
    if (F->ext_stats) {
        if (F->ext_stats->struct_size != (unsigned int)sizeof(*F->ext_stats) ||
            F->ext_stats->version != LIS_KRYLOV_EXT_VERSION)
            return k_return(K,F->parent,LIS_KRYLOV_BAD_ARGUMENT);
        lis_krylov_ext_stats_default(F->ext_stats);
    }
    if (!F->p || !F->o || !F->values || !F->vectors || !F->residuals || !F->stats ||
        F->p->n < 3 || F->o->nev < 1 || F->o->nev >= F->p->n || !F->p->apply_a ||
        F->p->transform < LIS_KRYLOV_STANDARD ||
        F->p->transform > LIS_KRYLOV_SHIFT_INVERT_A ||
        (F->p->transform != LIS_KRYLOV_STANDARD &&
         (!F->p->apply_b || !F->p->solve)) ||
        (F->p->transform == LIS_KRYLOV_SHIFT_INVERT_A && F->p->sigma == 0.0) ||
        !isfinite(F->p->sigma) || !isfinite(F->o->tol) || F->o->tol < 0.0 ||
        (F->o->hard_locking != 0 && F->o->hard_locking != 1) ||
        (F->o->allow_semidefinite_metric != 0 &&
         F->o->allow_semidefinite_metric != 1) ||
        (F->o->allow_semidefinite_metric &&
         F->p->transform != LIS_KRYLOV_SHIFT_INVERT_B) ||
        !isfinite(F->o->inner_tol) || F->o->inner_tol < 0.0 ||
        !isfinite(F->o->residual_scale) || F->o->residual_scale < 0.0 ||
        F->o->maxiter < 0)
        return k_return(K,F->parent,LIS_KRYLOV_BAD_ARGUMENT);
    F->which = F->o->which ? F->o->which : "LM";
    {
        F->which_len = strlen(F->which);
        if ((F->which_len != 2 && F->which_len != 6 && F->which_len != 7) ||
            (strncmp(F->which, "LM", 2) && strncmp(F->which, "SM", 2) &&
             strncmp(F->which, "LA", 2) && strncmp(F->which, "SA", 2) &&
             strncmp(F->which, "BE", 2)))
            return k_return(K,F->parent,LIS_KRYLOV_BAD_ARGUMENT);
        if (F->which_len > 2) {
            if (!strcmp(F->which + 2, ":ASC"))
                F->output_sort = 1;
            else if (!strcmp(F->which + 2, ":DESC"))
                F->output_sort = -1;
            else
                return k_return(K,F->parent,LIS_KRYLOV_BAD_ARGUMENT);
        }
        F->both_ends = F->which[0] == 'B' && F->which[1] == 'E';
    }

    F->nev = F->o->nev;
    F->maxiter = F->o->maxiter ? F->o->maxiter : 300;
    F->tol = F->o->tol ? F->o->tol : 1.0e-10;
    if(F->nev >= INT_MAX/8) return k_return(K,F->parent,LIS_KRYLOV_BAD_ARGUMENT);
    F->ncv = F->o->ncv ? F->o->ncv :
          (F->nev < (F->p->n - 8) / 3 ? 3 * F->nev + 8 : (int)F->p->n);
    if (F->ncv <= F->nev || F->ncv > F->p->n || F->ncv > INT_MAX / 8) return k_return(K,F->parent,LIS_KRYLOV_BAD_ARGUMENT);
    F->rc = validate_extension(F->ext, F->p->n, F->ncv, &F->capacity, &F->restart_directions);
    if (F->rc != LIS_KRYLOV_SUCCESS) return k_return(K,F->parent,F->rc);
    if (F->reuse) {
        if (F->reuse->n != F->p->n || F->reuse->ncv != F->ncv ||
            F->reuse->capacity_ncv != F->capacity ||
            F->reuse->nev != F->nev || F->reuse->transform != F->p->transform ||
            F->reuse->hard_locking != F->o->hard_locking)
            return k_return(K,F->parent,LIS_KRYLOV_BAD_ARGUMENT);
        F->s = F->reuse->buffers;
    } else {
        F->rc = allocate_workspace_buffers(&F->s, F->p->n, F->capacity, F->nev,
                                        F->p->transform, F->o->hard_locking);
        if (F->rc != LIS_KRYLOV_SUCCESS) return k_return(K,F->parent,F->rc);
    }
    if (F->ext_stats)
        F->ext_stats->ncv_initial = F->ext_stats->ncv_final = F->ncv;
    F->s.problem = F->p;
    F->s.options = F->o;
    F->s.stats = F->stats;
    F->s.ext_stats = F->ext_stats;
    F->s.norm_a = F->o->residual_scale;
    F->s.nlocked = 0;
    F->seed = F->o->random_seed ? F->o->random_seed : 1u;
    F->rc = LIS_KRYLOV_NO_CONVERGENCE;

    if (F->o->initial_vector) {
        memcpy(F->s.x, F->o->initial_vector, (size_t)F->s.n * sizeof(double));
    } else {
        seed_vector(F->s.n, F->s.x, F->seed);
    }
    if (F->o->residual_scale == 0.0) {
        F->norm_x = stable_norm(F->s.n, F->s.x);
        ;
        if (!isfinite(F->norm_x) || F->norm_x == 0.0) {
            F->rc = LIS_KRYLOV_BAD_ARGUMENT;
            goto done;
        }
        F->pc = 1; return k_suspend(K, LIS_RCI_OP_A, F->s.n, F->s.x, F->s.rhs);
        resume_1: F->pc = 0;
        if (K->retval) {
            F->rc = LIS_KRYLOV_CALLBACK_ERROR;
            goto done;
        }
        F->norm_a = stable_norm(F->s.n, F->s.rhs);
        if (!isfinite(F->norm_a)) {
            F->rc = LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
        F->s.norm_a = F->norm_a / F->norm_x;
    }




    {
        ;
        for (F->attempt = 0; F->attempt <
             (F->o->allow_semidefinite_metric ? 16 : 1); ++F->attempt) {
            if (F->attempt > 0)
                seed_vector(F->s.n, F->s.x, F->seed + (unsigned)F->attempt);
            K->f_append_range_seed.s = &F->s;
            K->f_append_range_seed.count = F->m;
            K->f_append_range_seed.x = F->s.x;
            K->f_append_range_seed.pc = 0; K->f_append_range_seed.parent = 9;
            F->pc = 2; K->active = 5; return K_STEP_AGAIN;
            resume_2: F->pc = 0;
            F->i = K->retval;
            if (F->i == 0) break;
            if (F->i < 0) {
                F->rc = F->i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     F->i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
        }
        if (F->i != 0) {
            F->rc = F->o->allow_semidefinite_metric ?
                 LIS_KRYLOV_NO_CONVERGENCE : LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
    }
    F->m = 1;

    for (F->iter = 0; F->iter < F->maxiter; ++F->iter) {
        F->old_m = F->m;
        F->extensions = 0;
        F->reuse_verified = 0;
        F->norm_a_before = F->s.norm_a;
        F->stats->iterations = F->iter + 1;
        update_projected(&F->s, F->previous, F->m);


        for (F->i = 0; F->i < F->m; ++F->i)
            memcpy(F->s.dense_work + (size_t)F->i * F->s.ncv,
                   F->s.projected + (size_t)F->i * F->s.ncv,
                   (size_t)F->m * sizeof(double));
        if (F->ext_stats) KRYLOV_INC(F->ext_stats->projected_solves);
        if (diagonalize(F->m, F->s.dense_work, F->s.ritz, F->s.eval, F->s.ncv)) {
            F->rc = LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
        sort_ritz(F->m, F->s.eval, F->s.order, F->which,
                  F->be_high_locked, F->be_low_locked, F->nev, F->s.locked);
        if (F->o->hard_locking && F->s.nlocked < F->nev) {
            F->captured = 0;
            memset(F->s.verified, 0, (size_t)F->nev * sizeof(*F->s.verified));


            for (F->j = 0; F->j < F->m && F->j < F->nev - F->s.nlocked; ++F->j) {
                F->idx = F->s.order[F->j];
                F->lambda = NAN;
                F->rel = INFINITY;
                K->f_eigenpair_residual.s = &F->s;
                K->f_eigenpair_residual.m = F->m;
                K->f_eigenpair_residual.index = F->idx;
                K->f_eigenpair_residual.mu = F->s.eval[F->idx];
                K->f_eigenpair_residual.lambda = &F->lambda;
                K->f_eigenpair_residual.relative = &F->rel;
                K->f_eigenpair_residual.pc = 0; K->f_eigenpair_residual.parent = 9;
                F->pc = 3; K->active = 8; return K_STEP_AGAIN;
                resume_3: F->pc = 0;
                F->err = K->retval;
                if (F->err < 0) {
                    F->rc = F->err == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }





                F->s.verified[F->j] = F->err == 0 ? 1 : -1;
                if (F->err == 0) {
                    F->slot = F->s.nlocked + F->j;
                    F->values[F->slot] = F->lambda;
                    F->residuals[F->slot] = F->rel;
                    memcpy(F->vectors + (size_t)F->slot * F->s.n, F->s.x,
                           (size_t)F->s.n * sizeof(double));
                }
                if (F->err != 0 || F->rel > F->tol) continue;
                K->f_apply_metric.s = &F->s;
                K->f_apply_metric.x = F->s.x;
                K->f_apply_metric.y = F->s.metric_x;
                K->f_apply_metric.pc = 0; K->f_apply_metric.parent = 9;
                F->pc = 4; K->active = 1; return K_STEP_AGAIN;
                resume_4: F->pc = 0;
                F->err = K->retval;
                if (F->err) {
                    F->rc = F->err == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                memcpy(F->s.locked_basis + (size_t)F->s.nlocked * F->s.n, F->s.x,
                       (size_t)F->s.n * sizeof(double));
                if (F->s.locked_metric != F->s.locked_basis)
                    memcpy(F->s.locked_metric + (size_t)F->s.nlocked * F->s.n,
                           F->s.metric_x, (size_t)F->s.n * sizeof(double));
                F->s.locked_values[F->s.nlocked] = F->lambda;
                F->s.locked_residuals[F->s.nlocked] = F->rel;
                ++F->s.nlocked;
                if (F->both_ends) {
                    if (F->s.locked[F->j]) ++F->be_high_locked;
                    else ++F->be_low_locked;
                }
                F->captured = 1;
                break;
            }
            if (F->captured) {
                F->ok = 0;
                F->stats->nreturned = F->stats->nconv = F->s.nlocked;
                F->stats->max_residual = 0.0;
                for (F->j = 0; F->j < F->s.nlocked; ++F->j) {
                    F->values[F->j] = F->s.locked_values[F->j];
                    F->residuals[F->j] = F->s.locked_residuals[F->j];
                    if (F->residuals[F->j] > F->stats->max_residual)
                        F->stats->max_residual = F->residuals[F->j];
                    memcpy(F->vectors + (size_t)F->j * F->s.n,
                           F->s.locked_basis + (size_t)F->j * F->s.n,
                           (size_t)F->s.n * sizeof(double));
                }
                if (F->s.nlocked == F->nev) {
                    F->rc = LIS_KRYLOV_SUCCESS;
                    goto done;
                }


                F->m = F->previous = 0;
                F->last_nconv = F->stats->nconv;
                F->stalled_restarts = 0;
                for (F->trial = 0; F->trial < 32; ++F->trial) {
                    ;
                    seed_vector(F->s.n, F->s.x,
                                F->seed + (unsigned)F->iter * 577u +
                                (unsigned)F->trial + 1001u);
                    K->f_append_range_seed.s = &F->s;
                    K->f_append_range_seed.count = 0;
                    K->f_append_range_seed.x = F->s.x;
                    K->f_append_range_seed.pc = 0; K->f_append_range_seed.parent = 9;
                    F->pc = 5; K->active = 5; return K_STEP_AGAIN;
                    resume_5: F->pc = 0;
                    F->attempt = K->retval;
                    if (F->attempt < 0) {
                        F->rc = F->attempt == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                             F->attempt == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                             LIS_KRYLOV_BREAKDOWN;
                        goto done;
                    }
                    if (F->attempt == 0) { F->m = 1; F->ok = 1; break; }
                }
                if (!F->ok) {
                    F->rc = LIS_KRYLOV_NO_CONVERGENCE;
                    goto done;
                }
                if (F->restart_directions > 1) {
                    K->f_append_extra_seeds.s = &F->s;
                    K->f_append_extra_seeds.count = &F->m;
                    K->f_append_extra_seeds.capacity = F->ncv;
                    K->f_append_extra_seeds.additional = F->restart_directions - 1;
                    K->f_append_extra_seeds.start_seed = F->seed + (unsigned)F->iter * 577u + 12001u;
                    K->f_append_extra_seeds.pc = 0; K->f_append_extra_seeds.parent = 9;
                    F->pc = 6; K->active = 7; return K_STEP_AGAIN;
                    resume_6: F->pc = 0;
                    F->added = K->retval;
                    if (F->added < 0) {
                        F->rc = F->added == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                             F->added == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                             LIS_KRYLOV_BREAKDOWN;
                        goto done;
                    }
                }
                KRYLOV_INC(F->stats->restarts);
                continue;
            }



            F->reuse_verified = (F->s.norm_a == F->norm_a_before);
        }
        F->stats->nreturned = F->s.nlocked;
        F->stats->nconv = F->s.nlocked;
        F->stats->max_residual = 0.0;
        for (F->j = 0; F->j < F->s.nlocked; ++F->j) {
            F->values[F->j] = F->s.locked_values[F->j];
            F->residuals[F->j] = F->s.locked_residuals[F->j];
            memcpy(F->vectors + (size_t)F->j * F->s.n,
                   F->s.locked_basis + (size_t)F->j * F->s.n,
                   (size_t)F->s.n * sizeof(double));
            if (F->residuals[F->j] > F->stats->max_residual)
                F->stats->max_residual = F->residuals[F->j];
        }
        memset(F->s.locked, 0, (size_t)F->m * sizeof(int));
        for (F->j = 0; F->j < F->m && F->stats->nreturned < F->nev; ++F->j) {
            F->idx = F->s.order[F->j];
            F->lambda = NAN;
            F->rel = INFINITY;
            F->from_cache = F->reuse_verified && F->j < F->nev &&
                             F->s.verified[F->j] != 0;
            F->cache_slot = F->s.nlocked + F->j;
            ;
            if (F->from_cache) {
                F->err = F->s.verified[F->j] < 0 ? 1 : 0;
                if (F->err == 0) {
                    F->lambda = F->values[F->cache_slot];
                    F->rel = F->residuals[F->cache_slot];
                }
            } else {
                K->f_eigenpair_residual.s = &F->s;
                K->f_eigenpair_residual.m = F->m;
                K->f_eigenpair_residual.index = F->idx;
                K->f_eigenpair_residual.mu = F->s.eval[F->idx];
                K->f_eigenpair_residual.lambda = &F->lambda;
                K->f_eigenpair_residual.relative = &F->rel;
                K->f_eigenpair_residual.pc = 0; K->f_eigenpair_residual.parent = 9;
                F->pc = 7; K->active = 8; return K_STEP_AGAIN;
                resume_7: F->pc = 0;
                F->err = K->retval;
            }
            if (F->err < 0) {
                F->rc = F->err == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }


            if (F->err == 1) continue;
            F->s.locked[F->j] = (F->rel <= F->tol);
            if (F->s.locked[F->j]) ++F->stats->nconv;
            if (F->rel > F->stats->max_residual) F->stats->max_residual = F->rel;
            F->values[F->stats->nreturned] = F->lambda;
            F->residuals[F->stats->nreturned] = F->rel;
            if (F->from_cache) {


                memmove(F->vectors + (size_t)F->stats->nreturned * F->s.n,
                        F->vectors + (size_t)F->cache_slot * F->s.n,
                        (size_t)F->s.n * sizeof(double));
            } else {
                memcpy(F->vectors + (size_t)F->stats->nreturned * F->s.n, F->s.x,
                       (size_t)F->s.n * sizeof(double));
            }
            ++F->stats->nreturned;
        }
        if (F->stats->nconv == F->nev) {
            F->rc = LIS_KRYLOV_SUCCESS;
            goto done;
        }
        if (F->m >= F->ncv) {


            if (F->stats->nconv > F->last_nconv)
                F->stalled_restarts = 0;
            else if (F->stalled_restarts < INT_MAX)
                ++F->stalled_restarts;
            F->last_nconv = F->stats->nconv;




            if (F->ext && F->ext->adaptive_ncv && F->ncv < F->capacity) {
                if (F->stalled_restarts >= 2) {
                    F->increment = F->ncv / 2;
                    F->remaining = F->capacity - F->ncv;
                    if (F->increment < 2) F->increment = 2;
                    if (F->increment > F->remaining) F->increment = F->remaining;
                    F->ncv += F->increment;
                    if (F->ext_stats) {
                        F->ext_stats->ncv_final = F->ncv;
                        KRYLOV_INC(F->ext_stats->ncv_growths);
                    }
                    F->stalled_restarts = 0;
                }
            }
            if (F->m >= F->ncv) {
                F->keep = F->m / 2;
                if (F->keep < F->nev) F->keep = F->nev;
                if (F->keep >= F->m) F->keep = F->m - 1;
                restart_basis(&F->s, F->m, F->keep);
                KRYLOV_INC(F->stats->restarts);
                F->m = F->keep;




                F->previous = F->restart_directions > 1 ? 0 : F->keep;

                for (F->i = 0; F->i < 4 && F->m < F->ncv; ++F->i) {
                    ;
                    seed_vector(F->s.n, F->s.x,
                                F->seed + (unsigned)F->iter * 997u + (unsigned)F->i + 17u);
                    K->f_append_range_seed.s = &F->s;
                    K->f_append_range_seed.count = F->m;
                    K->f_append_range_seed.x = F->s.x;
                    K->f_append_range_seed.pc = 0; K->f_append_range_seed.parent = 9;
                    F->pc = 8; K->active = 5; return K_STEP_AGAIN;
                    resume_8: F->pc = 0;
                    F->result = K->retval;
                    if (F->result < 0) {
                        F->rc = F->result == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                             F->result == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                             LIS_KRYLOV_BREAKDOWN;
                        goto done;
                    }
                    if (F->result == 0) { ++F->m; break; }
                }
                {
                    F->extra = F->restart_directions - 1;






                    if (!F->extra && F->o->allow_semidefinite_metric &&
                        F->stalled_restarts >= 1 && F->o->residual_scale == 0.0)
                        F->extra = 1;
                    if (F->extra) {
                        K->f_append_extra_seeds.s = &F->s;
                        K->f_append_extra_seeds.count = &F->m;
                        K->f_append_extra_seeds.capacity = F->ncv;
                        K->f_append_extra_seeds.additional = F->extra;
                        K->f_append_extra_seeds.start_seed = F->seed + (unsigned)F->iter * 997u + 13001u;
                        K->f_append_extra_seeds.pc = 0; K->f_append_extra_seeds.parent = 9;
                        F->pc = 9; K->active = 7; return K_STEP_AGAIN;
                        resume_9: F->pc = 0;
                        F->added = K->retval;
                        if (F->added < 0) {
                            F->rc = F->added == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                                 F->added == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                                 LIS_KRYLOV_BREAKDOWN;
                            goto done;
                        }
                    }
                }
                continue;
            }
        }

        for (F->j = 0; F->j < F->old_m && F->j < F->nev - F->s.nlocked && F->m < F->ncv; ++F->j) {
            F->idx = F->s.order[F->j];
            F->mu = F->s.eval[F->idx];
            ;
            if (F->s.locked[F->j]) continue;
            rotate_column(&F->s, F->s.images, F->old_m, F->idx, F->s.x);
            rotate_column(&F->s, F->s.basis, F->old_m, F->idx, F->s.metric_x);
            F->norm_w = stable_norm(F->s.n, F->s.x);
            F->norm_v = stable_norm(F->s.n, F->s.metric_x);
            for (F->row = 0; F->row < F->s.n; ++F->row)
                F->s.x[F->row] -= F->mu * F->s.metric_x[F->row];
            F->norm_r = stable_norm(F->s.n, F->s.x);
            F->threshold = KRYLOV_DIRECTION_REJECT *
                        (F->norm_w + fabs(F->mu) * F->norm_v);
            if (F->norm_r <= F->threshold) continue;
            K->f_append_vector.s = &F->s;
            K->f_append_vector.count = F->m;
            K->f_append_vector.x = F->s.x;
            K->f_append_vector.pc = 0; K->f_append_vector.parent = 9;
            F->pc = 10; K->active = 6; return K_STEP_AGAIN;
            resume_10: F->pc = 0;
            F->i = K->retval;
            if (F->i < 0) {
                F->rc = F->i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     F->i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
            if (F->i == 0) { ++F->m; ++F->extensions; }
        }
        F->previous = F->old_m;
        if (!F->extensions && F->m < F->ncv) {
            ;
            for (F->trial = 0; F->trial < 6 && !F->extensions; ++F->trial) {
                seed_vector(F->s.n, F->s.x,
                            F->seed + (unsigned)F->iter * 7u + (unsigned)F->trial + 1u);
                K->f_append_range_seed.s = &F->s;
                K->f_append_range_seed.count = F->m;
                K->f_append_range_seed.x = F->s.x;
                K->f_append_range_seed.pc = 0; K->f_append_range_seed.parent = 9;
                F->pc = 11; K->active = 5; return K_STEP_AGAIN;
                resume_11: F->pc = 0;
                F->i = K->retval;
                if (F->i < 0) {
                    F->rc = F->i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         F->i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                if (F->i == 0) { ++F->m; ++F->extensions; }
            }
            if (F->extensions && F->restart_directions > 1) {
                K->f_append_extra_seeds.s = &F->s;
                K->f_append_extra_seeds.count = &F->m;
                K->f_append_extra_seeds.capacity = F->ncv;
                K->f_append_extra_seeds.additional = F->restart_directions - 1;
                K->f_append_extra_seeds.start_seed = F->seed + (unsigned)F->iter * 7u + 14001u;
                K->f_append_extra_seeds.pc = 0; K->f_append_extra_seeds.parent = 9;
                F->pc = 12; K->active = 7; return K_STEP_AGAIN;
                resume_12: F->pc = 0;
                F->added = K->retval;
                if (F->added < 0) {
                    F->rc = F->added == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         F->added == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                F->extensions += F->added;
            }
            if (!F->extensions) {



                F->rc = F->o->allow_semidefinite_metric ?
                     LIS_KRYLOV_NO_CONVERGENCE : LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
        }
    }
    F->rc = LIS_KRYLOV_NO_CONVERGENCE;

done:



    if (F->s.scratch &&
        (F->rc == LIS_KRYLOV_SUCCESS || F->rc == LIS_KRYLOV_NO_CONVERGENCE)) {
        if (F->output_sort || F->both_ends) {
            reorder_results(F->s.n, F->stats->nreturned, F->p->transform, F->p->sigma,
                            F->both_ends, F->output_sort, F->values, F->residuals,
                            F->vectors, F->s.scratch, F->s.order, F->s.locked);
        } else if (F->o && F->o->hard_locking) {
            for (F->i = 1; F->i < F->stats->nreturned; ++F->i) {
                for (F->j = F->i; F->j > 0; --F->j) {
                    F->left = F->values[F->j];
                    F->right = F->values[F->j-1];
                    ;
                    if (F->p->transform == LIS_KRYLOV_STANDARD) {
                        F->mu_left = F->left;
                        F->mu_right = F->right;
                    } else if (F->p->transform == LIS_KRYLOV_SHIFT_INVERT_B) {
                        F->mu_left = 1.0 / (F->left - F->p->sigma);
                        F->mu_right = 1.0 / (F->right - F->p->sigma);
                    } else {
                        F->mu_left = F->left / (F->left - F->p->sigma);
                        F->mu_right = F->right / (F->right - F->p->sigma);
                    }
                    if (!earlier(F->mu_left, F->mu_right, F->which)) break;
                    F->tmp = F->values[F->j-1];
                    F->values[F->j-1] = F->values[F->j];
                    F->values[F->j] = F->tmp;
                    F->tmp = F->residuals[F->j-1];
                    F->residuals[F->j-1] = F->residuals[F->j];
                    F->residuals[F->j] = F->tmp;
                    memcpy(F->s.scratch, F->vectors + (size_t)(F->j-1) * F->s.n,
                           (size_t)F->s.n * sizeof(double));
                    memcpy(F->vectors + (size_t)(F->j-1) * F->s.n,
                           F->vectors + (size_t)F->j * F->s.n,
                           (size_t)F->s.n * sizeof(double));
                    memcpy(F->vectors + (size_t)F->j * F->s.n, F->s.scratch,
                           (size_t)F->s.n * sizeof(double));
                }
            }
        }
    }
    if (F->rc != LIS_KRYLOV_SUCCESS && F->rc != LIS_KRYLOV_NO_CONVERGENCE) {
        F->stats->nreturned = 0;
        F->stats->nconv = 0;
    }
    if (!F->reuse) release_workspace(&F->s);
    return k_return(K,F->parent,F->rc);

}
static int k_step(lis_rci_krylov_restore *K)
{ switch(K->active){
case 1:return step_apply_metric(K);
case 2:return step_verify_inverse(K);
case 3:return step_apply_operator(K);
case 4:return step_orthogonalize(K);
case 5:return step_append_range_seed(K);
case 6:return step_append_vector(K);
case 7:return step_append_extra_seeds(K);
case 8:return step_eigenpair_residual(K);
case 9:return step_lis_krylov_symmetric_core(K);
default:return K_STEP_STOP;
}}
