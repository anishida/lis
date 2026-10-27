/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Generic exact/near-nullspace candidate classification.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <float.h>
#include <math.h>
#include <string.h>

#include "lislib.h"
#include "lis_nullspace_candidates.h"


static LIS_INT
lis_solver_candidate_vector_compatible(
    LIS_MATRIX A,
    LIS_VECTOR v)
{
    if(A==NULL || v==NULL)
        return LIS_FALSE;

    if(v->label!=LIS_LABEL_VECTOR)
        return LIS_FALSE;

    if(v->gn!=A->gn ||
       v->n!=A->n ||
       v->np!=A->np ||
       v->pad!=A->pad ||
       v->origin!=A->origin ||
       v->my_rank!=A->my_rank ||
       v->nprocs!=A->nprocs ||
       v->is!=A->is ||
       v->ie!=A->ie)
        return LIS_FALSE;

    return LIS_TRUE;
}


static void
lis_solver_candidate_destroy_array(
    LIS_VECTOR *vectors,
    LIS_INT count)
{
    LIS_INT i;

    if(vectors==NULL)
        return;

    for(i=0;i<count;i++)
    {
        if(vectors[i])
            lis_vector_destroy(vectors[i]);
    }

    lis_free(vectors);
}


static LIS_INT
lis_solver_candidate_append_copy(
    LIS_MATRIX A,
    LIS_VECTOR source,
    LIS_VECTOR *vectors,
    LIS_INT *count)
{
    LIS_VECTOR v;
    LIS_INT err;

    if(A==NULL || source==NULL || vectors==NULL || count==NULL)
        return LIS_ERR_ILL_ARG;

    v = NULL;

    err = lis_vector_duplicate(A,&v);
    if(err)
        return err;

    err = lis_vector_copy(source,v);
    if(err)
    {
        lis_vector_destroy(v);
        return err;
    }

    vectors[*count] = v;
    (*count)++;

    return LIS_SUCCESS;
}


static LIS_INT
lis_solver_candidate_orthogonalize(
    LIS_VECTOR vector,
    LIS_VECTOR basis[],
    LIS_INT rank,
    LIS_REAL *norm0,
    LIS_REAL *norm)
{
    LIS_INT i,pass,err;
    LIS_SCALAR alpha;

    err = lis_vector_nrm2(vector,norm0);
    if(err)
        return err;

    if(!isfinite((double)(*norm0)) || *norm0<=0.0)
    {
        *norm = 0.0;
        return LIS_SUCCESS;
    }

    for(pass=0;pass<2;pass++)
    {
        for(i=0;i<rank;i++)
        {
            err = lis_vector_dot(basis[i],vector,&alpha);
            if(err)
                return err;

            err = lis_vector_axpy(-alpha,basis[i],vector);
            if(err)
                return err;
        }
    }

    err = lis_vector_nrm2(vector,norm);
    if(err)
        return err;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_candidate_jacobi_pair"
static LIS_INT
lis_solver_candidate_jacobi_pair(
    LIS_VECTOR left_basis,
    LIS_VECTOR right_basis,
    LIS_VECTOR left_image,
    LIS_VECTOR right_image,
    LIS_VECTOR tmp_basis,
    LIS_VECTOR tmp_image,
    LIS_INT *changed)
{
    LIS_REAL alpha,beta,gabs,threshold;
    LIS_REAL zeta,t,c,s;
    LIS_SCALAR gamma,phase;
    LIS_INT err;

    err = lis_vector_nrm2(left_image,&alpha);
    if(err)
        return err;

    err = lis_vector_nrm2(right_image,&beta);
    if(err)
        return err;

    alpha *= alpha;
    beta *= beta;

    err = lis_vector_dot(left_image,right_image,&gamma);
    if(err)
        return err;

    gabs = (LIS_REAL)fabs(gamma);

    if(!isfinite((double)alpha) ||
       !isfinite((double)beta) ||
       !isfinite((double)gabs))
    {
        LIS_SETERR(
            LIS_BREAKDOWN,
            "nullspace candidate Jacobi state became non-finite\n");
        return LIS_BREAKDOWN;
    }

    if(alpha<=DBL_MIN || beta<=DBL_MIN || gabs<=DBL_MIN)
        return LIS_SUCCESS;

    threshold =
        (LIS_REAL)DBL_EPSILON
        *
        sqrt(alpha*beta);

    if(gabs<=threshold)
        return LIS_SUCCESS;

    /*
     * Rotate the right column by a unit phase so that the cross term
     * becomes real and nonnegative.  The following Jacobi rotation is
     * then the same stable real two-column update used by one-sided SVD.
     */
    phase = conj(gamma)/(LIS_SCALAR)gabs;

    err = lis_vector_scale(phase,right_basis);
    if(err)
        return err;

    err = lis_vector_scale(phase,right_image);
    if(err)
        return err;

    zeta =
        (beta-alpha)
        /
        ((LIS_REAL)2.0*gabs);

    t =
        copysign((LIS_REAL)1.0,zeta)
        /
        (fabs(zeta)+hypot((LIS_REAL)1.0,zeta));

    c =
        (LIS_REAL)1.0
        /
        sqrt((LIS_REAL)1.0+t*t);

    s = c*t;

    err = lis_vector_copy(left_basis,tmp_basis);
    if(err)
        return err;

    err = lis_vector_copy(left_image,tmp_image);
    if(err)
        return err;

    err = lis_vector_scale((LIS_SCALAR)c,left_basis);
    if(err)
        return err;

    err = lis_vector_axpy((LIS_SCALAR)(-s),right_basis,left_basis);
    if(err)
        return err;

    err = lis_vector_scale((LIS_SCALAR)c,right_basis);
    if(err)
        return err;

    err = lis_vector_axpy((LIS_SCALAR)s,tmp_basis,right_basis);
    if(err)
        return err;

    err = lis_vector_scale((LIS_SCALAR)c,left_image);
    if(err)
        return err;

    err = lis_vector_axpy((LIS_SCALAR)(-s),right_image,left_image);
    if(err)
        return err;

    err = lis_vector_scale((LIS_SCALAR)c,right_image);
    if(err)
        return err;

    err = lis_vector_axpy((LIS_SCALAR)s,tmp_image,right_image);
    if(err)
        return err;

    *changed = LIS_TRUE;

    return LIS_SUCCESS;
}


static LIS_INT
lis_solver_candidate_group_seen(
    const LIS_INT groups[],
    LIS_INT index)
{
    LIS_INT i;
    LIS_INT group;

    if(groups==NULL)
        return index>0 ? LIS_TRUE : LIS_FALSE;

    group = groups[index];

    for(i=0;i<index;i++)
    {
        if(groups[i]==group)
            return LIS_TRUE;
    }

    return LIS_FALSE;
}


#undef __FUNC__
#define __FUNC__ "lis_solver_classify_nullspace_candidates"
LIS_INT
lis_solver_classify_nullspace_candidates(
    LIS_MATRIX A,
    LIS_SOLVER solver,
    LIS_INT nvec,
    LIS_VECTOR candidates[],
    const LIS_INT groups[],
    LIS_REAL operator_scale,
    LIS_REAL exact_tolerance,
    LIS_REAL near_tolerance,
    LIS_INT flags,
    LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS *diagnostics)
{
    const LIS_REAL dependence_tol = (LIS_REAL)1.0e-12;
    const LIS_INT max_sweeps = 64;

    LIS_VECTOR *basis = NULL;
    LIS_VECTOR *image = NULL;
    LIS_VECTOR *exact = NULL;
    LIS_VECTOR *near = NULL;
    LIS_VECTOR *left = NULL;
    LIS_VECTOR tmp_basis = NULL;
    LIS_VECTOR tmp_image = NULL;
    LIS_VECTOR left_image = NULL;

    LIS_INT exact_count = 0;
    LIS_INT near_count = 0;
    LIS_INT left_count = 0;
    LIS_INT independent_count = 0;
    LIS_INT i,j,k,rank,sweep,changed,err;
    LIS_INT group;

    LIS_REAL norm0,norm,defect;
    LIS_REAL max_exact_defect = 0.0;
    LIS_REAL max_near_defect = 0.0;

    if(diagnostics)
        memset(diagnostics,0,sizeof(*diagnostics));

    if(A==NULL || solver==NULL || nvec<0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid nullspace candidate arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if(nvec==0)
    {
        err = lis_solver_clear_nullspace(solver);
        if(err)
            return err;

        err = lis_solver_clear_near_nullspace(solver);
        if(err)
            return err;

        err = lis_solver_clear_left_nullspace(solver);
        if(err)
            return err;

        return LIS_SUCCESS;
    }

    if(candidates==NULL ||
       !isfinite((double)operator_scale) ||
       operator_scale<=0.0 ||
       !isfinite((double)exact_tolerance) ||
       exact_tolerance<0.0 ||
       !isfinite((double)near_tolerance) ||
       near_tolerance<exact_tolerance)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid nullspace candidate configuration\n");
        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<nvec;i++)
    {
        if(!lis_solver_candidate_vector_compatible(A,candidates[i]))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,"nullspace candidate layout is incompatible with matrix\n");
            return LIS_ERR_ILL_ARG;
        }
    }

    basis = (LIS_VECTOR *)lis_malloc(
        nvec*sizeof(LIS_VECTOR),
        "lis_solver_classify_nullspace_candidates::basis");
    image = (LIS_VECTOR *)lis_malloc(
        nvec*sizeof(LIS_VECTOR),
        "lis_solver_classify_nullspace_candidates::image");
    exact = (LIS_VECTOR *)lis_malloc(
        nvec*sizeof(LIS_VECTOR),
        "lis_solver_classify_nullspace_candidates::exact");
    near = (LIS_VECTOR *)lis_malloc(
        nvec*sizeof(LIS_VECTOR),
        "lis_solver_classify_nullspace_candidates::near");
    left = (LIS_VECTOR *)lis_malloc(
        nvec*sizeof(LIS_VECTOR),
        "lis_solver_classify_nullspace_candidates::left");

    if(basis==NULL || image==NULL || exact==NULL || near==NULL || left==NULL)
    {
        err = LIS_OUT_OF_MEMORY;
        goto cleanup;
    }

    for(i=0;i<nvec;i++)
    {
        basis[i] = NULL;
        image[i] = NULL;
        exact[i] = NULL;
        near[i] = NULL;
        left[i] = NULL;
    }

    err = lis_vector_duplicate(A,&tmp_basis);
    if(err)
        goto cleanup;

    err = lis_vector_duplicate(A,&tmp_image);
    if(err)
        goto cleanup;

    for(i=0;i<nvec;i++)
    {
        if(lis_solver_candidate_group_seen(groups,i))
            continue;

        group = groups ? groups[i] : 0;
        rank = 0;

        for(j=0;j<nvec;j++)
        {
            if(groups && groups[j]!=group)
                continue;

            if(!groups && i!=0)
                continue;

            err = lis_vector_duplicate(A,&basis[rank]);
            if(err)
                goto cleanup;

            err = lis_vector_copy(candidates[j],basis[rank]);
            if(err)
                goto cleanup;

            err = lis_solver_candidate_orthogonalize(
                basis[rank],basis,rank,&norm0,&norm);
            if(err)
                goto cleanup;

            if(!isfinite((double)norm0) ||
               norm0<=0.0 ||
               !isfinite((double)norm) ||
               norm<=dependence_tol*norm0)
            {
                lis_vector_destroy(basis[rank]);
                basis[rank] = NULL;
                continue;
            }

            err = lis_vector_scale(
                (LIS_SCALAR)((LIS_REAL)1.0/norm),
                basis[rank]);
            if(err)
                goto cleanup;

            err = lis_vector_duplicate(A,&image[rank]);
            if(err)
                goto cleanup;

            rank++;
        }

        independent_count += rank;

        if(rank<=0)
            continue;

        err = lis_matvec_block(A,rank,basis,image);
        if(err)
            goto cleanup;

        for(j=0;j<rank;j++)
        {
            err = lis_vector_scale(
                (LIS_SCALAR)((LIS_REAL)1.0/operator_scale),
                image[j]);
            if(err)
                goto cleanup;
        }

        for(sweep=0;sweep<max_sweeps;sweep++)
        {
            changed = LIS_FALSE;

            for(j=0;j<rank;j++)
            {
                for(k=j+1;k<rank;k++)
                {
                    err = lis_solver_candidate_jacobi_pair(
                        basis[j],basis[k],image[j],image[k],
                        tmp_basis,tmp_image,&changed);
                    if(err)
                        goto cleanup;
                }
            }

            if(!changed)
                break;
        }

        for(j=0;j<rank;j++)
        {
            LIS_INT is_exact;

            err = lis_vector_nrm2(image[j],&defect);
            if(err)
                goto cleanup;

            if(!isfinite((double)defect))
            {
                err = LIS_BREAKDOWN;
                LIS_SETERR(LIS_BREAKDOWN,"nullspace candidate defect is non-finite\n");
                goto cleanup;
            }

            is_exact = defect<=exact_tolerance;

            if(is_exact)
            {
                err = lis_solver_candidate_append_copy(
                    A,basis[j],exact,&exact_count);
                if(err)
                    goto cleanup;

                if(defect>max_exact_defect)
                    max_exact_defect = defect;
            }

            if((flags & LIS_NULLSPACE_CANDIDATE_BIND_NEAR) &&
               defect<=near_tolerance &&
               (is_exact==LIS_FALSE ||
                (flags & LIS_NULLSPACE_CANDIDATE_NEAR_INCLUDE_EXACT)))
            {
                err = lis_solver_candidate_append_copy(
                    A,basis[j],near,&near_count);
                if(err)
                    goto cleanup;

                if(defect>max_near_defect)
                    max_near_defect = defect;
            }
        }

        for(j=0;j<rank;j++)
        {
            if(basis[j])
            {
                lis_vector_destroy(basis[j]);
                basis[j] = NULL;
            }

            if(image[j])
            {
                lis_vector_destroy(image[j]);
                image[j] = NULL;
            }
        }
    }

    /*
     * Left and right nullspaces are independent for nonsymmetric
     * operators. Keep every right-exact vector selected above. When a
     * left basis is requested, classify the same application candidates
     * against A^H into separate storage instead of shrinking the right
     * basis. Applications with distinct physical left candidates can call
     * this routine separately for the adjoint-side candidate set.
     */
    if((flags & LIS_NULLSPACE_CANDIDATE_BIND_LEFT) && exact_count>0)
    {
        err = lis_vector_duplicate(A,&left_image);
        if(err)
            goto cleanup;

        for(i=0;i<exact_count;i++)
        {
            err = lis_matvech(A,exact[i],left_image);
            if(err)
                goto cleanup;

            err = lis_vector_nrm2(left_image,&defect);
            if(err)
                goto cleanup;

            defect /= operator_scale;

            if(isfinite((double)defect) && defect<=exact_tolerance)
            {
                err = lis_solver_candidate_append_copy(
                    A,exact[i],left,&left_count);
                if(err)
                    goto cleanup;
            }
        }
    }

    if(exact_count>0)
        err = lis_solver_set_nullspace(solver,exact_count,exact);
    else
        err = lis_solver_clear_nullspace(solver);

    if(err)
        goto cleanup;

    if(flags & LIS_NULLSPACE_CANDIDATE_BIND_NEAR)
    {
        if(near_count>0)
            err = lis_solver_set_near_nullspace(solver,near_count,near);
        else
            err = lis_solver_clear_near_nullspace(solver);
    }
    else
    {
        err = lis_solver_clear_near_nullspace(solver);
    }

    if(err)
        goto cleanup;

    if(flags & LIS_NULLSPACE_CANDIDATE_BIND_LEFT)
    {
        if(left_count>0)
            err = lis_solver_set_left_nullspace(solver,left_count,left);
        else
            err = lis_solver_clear_left_nullspace(solver);
    }
    else
    {
        err = lis_solver_clear_left_nullspace(solver);
    }

    if(err)
        goto cleanup;

    if(diagnostics)
    {
        diagnostics->input_count = nvec;
        diagnostics->independent_count = independent_count;
        diagnostics->exact_count = exact_count;
        diagnostics->near_count = near_count;
        diagnostics->left_exact_count = left_count;
        diagnostics->operator_scale = operator_scale;
        diagnostics->max_exact_defect = max_exact_defect;
        diagnostics->max_near_defect = max_near_defect;
    }

    err = LIS_SUCCESS;

cleanup:

    if(tmp_basis)
        lis_vector_destroy(tmp_basis);
    if(tmp_image)
        lis_vector_destroy(tmp_image);
    if(left_image)
        lis_vector_destroy(left_image);

    if(basis)
    {
        for(i=0;i<nvec;i++)
        {
            if(basis[i])
                lis_vector_destroy(basis[i]);
        }
        lis_free(basis);
    }

    if(image)
    {
        for(i=0;i<nvec;i++)
        {
            if(image[i])
                lis_vector_destroy(image[i]);
        }
        lis_free(image);
    }

    lis_solver_candidate_destroy_array(exact,nvec);
    lis_solver_candidate_destroy_array(near,nvec);
    lis_solver_candidate_destroy_array(left,nvec);

    return err;
}
