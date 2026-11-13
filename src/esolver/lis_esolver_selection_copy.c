/* Generic copying and streaming of selected converged eigensolver pairs. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lislib.h"
#include "lis_esolver_selection.h"

#undef __FUNC__
#define __FUNC__ "lis_esolver_copy_selected_eigenpairs"
LIS_INT
lis_esolver_copy_selected_eigenpairs(
    LIS_ESOLVER esolver,
    LIS_INT scan_count,
    LIS_INT requested_count,
    LIS_REAL residual_tolerance,
    LIS_INT flags,
    LIS_SCALAR values[],
    LIS_VECTOR vectors[],
    LIS_INT *available_count)
{
    LIS_ESOLVER_SELECTED_PAIR *pairs = NULL;
    LIS_INT i,available=0,err;

    LIS_DEBUG_FUNC_IN;
    if(esolver==NULL || scan_count<0 || requested_count<0 ||
       (requested_count>0 && values==NULL && vectors==NULL))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid selected eigenpair copy arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if(requested_count>0)
    {
        pairs=(LIS_ESOLVER_SELECTED_PAIR *)lis_malloc(
            requested_count*sizeof(*pairs),
            "lis_esolver_copy_selected_eigenpairs::pairs");
        if(pairs==NULL)
        {
            LIS_SETERR_MEM(requested_count*sizeof(*pairs));
            return LIS_OUT_OF_MEMORY;
        }
    }

    err=lis_esolver_select_converged_pairs(
        esolver,scan_count,requested_count,residual_tolerance,
        flags,pairs,&available);
    if(err) goto cleanup;
    if(available_count) *available_count=available;

    for(i=0;i<requested_count;i++)
    {
        if(values) values[i]=pairs[i].value;
        if(vectors)
        {
            if(vectors[i]==NULL)
            {
                err=LIS_ERR_ILL_ARG;
                LIS_SETERR(LIS_ERR_ILL_ARG,"selected eigenpair output vector is NULL\n");
                goto cleanup;
            }
            err=lis_esolver_get_specific_evector(esolver,pairs[i].mode,vectors[i]);
            if(err) goto cleanup;
        }
    }

cleanup:
    lis_free(pairs);
    LIS_DEBUG_FUNC_OUT;
    return err;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_visit_converged_pairs"
LIS_INT
lis_esolver_visit_converged_pairs(
    LIS_ESOLVER esolver,
    LIS_INT scan_count,
    LIS_INT requested_count,
    LIS_REAL residual_tolerance,
    LIS_INT flags,
    LIS_VECTOR work_vector,
    LIS_ESOLVER_SELECTED_PAIR_VISITOR visitor,
    void *context,
    LIS_INT *available_count)
{
    LIS_ESOLVER_SELECTED_PAIR *pairs=NULL;
    LIS_INT i,available=0,err=LIS_SUCCESS;

    LIS_DEBUG_FUNC_IN;
    if(esolver==NULL || scan_count<0 || requested_count<0 ||
       (requested_count>0 && (work_vector==NULL || visitor==NULL)))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid selected eigenpair visitor arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if(requested_count>0)
    {
        pairs=(LIS_ESOLVER_SELECTED_PAIR *)lis_malloc(
            requested_count*sizeof(*pairs),
            "lis_esolver_visit_converged_pairs::pairs");
        if(pairs==NULL)
        {
            LIS_SETERR_MEM(requested_count*sizeof(*pairs));
            return LIS_OUT_OF_MEMORY;
        }
    }

    err=lis_esolver_select_converged_pairs(
        esolver,scan_count,requested_count,residual_tolerance,
        flags,pairs,&available);
    if(err) goto cleanup;
    if(available_count) *available_count=available;

    for(i=0;i<requested_count;i++)
    {
        err=lis_esolver_get_specific_evector(
            esolver,pairs[i].mode,work_vector);
        if(err) goto cleanup;

        err=visitor(
            context,i,pairs[i].mode,pairs[i].value,
            pairs[i].residual,work_vector);
        if(err) goto cleanup;
    }

cleanup:
    lis_free(pairs);
    LIS_DEBUG_FUNC_OUT;
    return err;
}
