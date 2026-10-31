/* Copyright (C) 2005 The Scalable Software Infrastructure Project.
 *
 * Generic component-localized nullspace candidate classification.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "lislib.h"
#include "lis_matrix_graph.h"
#include "lis_nullspace_candidates.h"


static void
lis_solver_component_candidate_destroy_array(
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


#undef __FUNC__
#define __FUNC__ "lis_solver_classify_component_nullspace_candidates"
LIS_INT
lis_solver_classify_component_nullspace_candidates(
    LIS_MATRIX A,
    LIS_SOLVER solver,
    LIS_INT nvec,
    LIS_VECTOR candidates[],
    LIS_REAL graph_drop_tolerance,
    LIS_REAL operator_scale,
    LIS_REAL exact_tolerance,
    LIS_REAL near_tolerance,
    LIS_INT flags,
    LIS_COMPONENT_NULLSPACE_DIAGNOSTICS *diagnostics)
{
    LIS_VECTOR *localized = NULL;
    LIS_INT *groups = NULL;
    LIS_INT *labels = NULL;
    LIS_INT component_count = 0;
    LIS_INT localized_count = 0;
    LIS_INT capacity = 0;
    LIS_INT c,j,i,err;
    LIS_REAL norm;
    LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS classifier;

    LIS_DEBUG_FUNC_IN;

    memset(&classifier,0,sizeof(classifier));
    if(diagnostics)
        memset(diagnostics,0,sizeof(*diagnostics));

    if(A==NULL || solver==NULL || nvec<0 ||
       !(graph_drop_tolerance>=(LIS_REAL)0.0))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "invalid component nullspace candidate arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if(nvec==0)
    {
        err = lis_solver_classify_nullspace_candidates(
            A,solver,0,NULL,NULL,
            operator_scale,exact_tolerance,near_tolerance,
            flags,&classifier);
        if(!err && diagnostics)
            diagnostics->classifier = classifier;
        return err;
    }

    if(candidates==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "component nullspace candidates are missing\n");
        return LIS_ERR_ILL_ARG;
    }

    if(A->nprocs!=1 || A->n!=A->gn)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "distributed component nullspace localization is not implemented\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }

    for(j=0;j<nvec;j++)
    {
        if(candidates[j]==NULL ||
           candidates[j]->label!=LIS_LABEL_VECTOR ||
           candidates[j]->n!=A->n ||
           candidates[j]->gn!=A->gn)
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "component nullspace candidate layout is incompatible with matrix\n");
            return LIS_ERR_ILL_ARG;
        }
    }

    labels = (LIS_INT *)lis_malloc(
        A->n*sizeof(LIS_INT),
        "lis_solver_classify_component_nullspace_candidates::labels");
    if(labels==NULL && A->n>0)
    {
        LIS_SETERR_MEM(A->n*sizeof(LIS_INT));
        return LIS_ERR_OUT_OF_MEMORY;
    }

    err = lis_matrix_connected_components(
        A,graph_drop_tolerance,labels,&component_count);
    if(err)
        goto cleanup;

    if(component_count<0 ||
       (component_count>0 &&
        (uintmax_t)nvec>(uintmax_t)SIZE_MAX/(uintmax_t)component_count))
    {
        err = LIS_ERR_ILL_ARG;
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "component nullspace candidate count overflows address range\n");
        goto cleanup;
    }

    capacity = nvec*component_count;

    if(capacity>0)
    {
        localized = (LIS_VECTOR *)lis_malloc(
            capacity*sizeof(LIS_VECTOR),
            "lis_solver_classify_component_nullspace_candidates::localized");
        groups = (LIS_INT *)lis_malloc(
            capacity*sizeof(LIS_INT),
            "lis_solver_classify_component_nullspace_candidates::groups");

        if(localized==NULL || groups==NULL)
        {
            err = LIS_ERR_OUT_OF_MEMORY;
            LIS_SETERR_MEM(
                capacity*(sizeof(LIS_VECTOR)+sizeof(LIS_INT)));
            goto cleanup;
        }

        for(i=0;i<capacity;i++)
            localized[i] = NULL;
    }

    for(c=0;c<component_count;c++)
    {
        for(j=0;j<nvec;j++)
        {
            LIS_VECTOR v = NULL;

            err = lis_vector_duplicate(A,&v);
            if(err)
                goto cleanup;

            err = lis_vector_set_all((LIS_SCALAR)0.0,v);
            if(err)
            {
                lis_vector_destroy(v);
                goto cleanup;
            }

            for(i=0;i<A->n;i++)
            {
                if(labels[i]==c)
                    v->value[i] = candidates[j]->value[i];
            }

            err = lis_vector_nrm2(v,&norm);
            if(err)
            {
                lis_vector_destroy(v);
                goto cleanup;
            }

            if(!isfinite((double)norm))
            {
                lis_vector_destroy(v);
                err = LIS_BREAKDOWN;
                LIS_SETERR(LIS_BREAKDOWN,
                    "component nullspace candidate norm is non-finite\n");
                goto cleanup;
            }

            if(norm<=0.0)
            {
                lis_vector_destroy(v);
                continue;
            }

            localized[localized_count] = v;
            groups[localized_count] = c;
            localized_count++;
        }
    }

    err = lis_solver_classify_nullspace_candidates(
        A,solver,localized_count,localized,groups,
        operator_scale,exact_tolerance,near_tolerance,
        flags,&classifier);
    if(err)
        goto cleanup;

    if(diagnostics)
    {
        diagnostics->component_count = component_count;
        diagnostics->localized_candidate_count = localized_count;
        diagnostics->classifier = classifier;
    }

cleanup:
    lis_solver_component_candidate_destroy_array(localized,localized_count);
    lis_free(groups);
    lis_free(labels);

    LIS_DEBUG_FUNC_OUT;
    return err;
}
