/* Generic converged eigensolver-pair selection. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include "lislib.h"
#include "lis_esolver_selection.h"

static LIS_REAL
lis_esolver_selected_real(LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (LIS_REAL)creal(value);
#else
    return (LIS_REAL)value;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_select_converged_pairs"
LIS_INT
lis_esolver_select_converged_pairs(
    LIS_ESOLVER esolver,
    LIS_INT scan_count,
    LIS_INT requested_count,
    LIS_REAL residual_tolerance,
    LIS_INT flags,
    LIS_ESOLVER_SELECTED_PAIR pairs[],
    LIS_INT *available_count)
{
    LIS_ESOLVER_SELECTED_PAIR *work = NULL;
    LIS_INT i,j,best,count,err;

    if(available_count)
        *available_count = 0;

    if(esolver==NULL ||
       scan_count<0 ||
       requested_count<0 ||
       requested_count>scan_count ||
       !isfinite((double)residual_tolerance) ||
       residual_tolerance<0.0 ||
       (requested_count>0 && pairs==NULL))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid eigensolver selection arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if(scan_count==0)
        return requested_count==0 ? LIS_SUCCESS : LIS_BREAKDOWN;

    work = (LIS_ESOLVER_SELECTED_PAIR *)lis_malloc(
        scan_count*sizeof(LIS_ESOLVER_SELECTED_PAIR),
        "lis_esolver_select_converged_pairs::work");
    if(work==NULL)
    {
        LIS_SETERR_MEM(scan_count*sizeof(LIS_ESOLVER_SELECTED_PAIR));
        return LIS_OUT_OF_MEMORY;
    }

    count = 0;

    for(i=0;i<scan_count;i++)
    {
        LIS_SCALAR value;
        LIS_REAL residual;
        LIS_REAL real_value;

        err = lis_esolver_get_specific_evalue(esolver,i,&value);
        if(err!=LIS_SUCCESS)
        {
            lis_free(work);
            return err;
        }

        err = lis_esolver_get_specific_residualnorm(esolver,i,&residual);
        if(err!=LIS_SUCCESS)
        {
            lis_free(work);
            return err;
        }

        real_value = lis_esolver_selected_real(value);

        if(!isfinite((double)real_value) ||
           !isfinite((double)residual) ||
           residual<0.0 ||
           residual>residual_tolerance)
            continue;

        if((flags & LIS_ESOLVER_SELECT_POSITIVE_REAL) &&
           real_value<=0.0)
            continue;

        work[count].mode = i;
        work[count].value = value;
        work[count].residual = residual;
        count++;
    }

    if(flags & LIS_ESOLVER_SELECT_ASCENDING_REAL)
    {
        for(i=0;i<count-1;i++)
        {
            best = i;

            for(j=i+1;j<count;j++)
            {
                if(lis_esolver_selected_real(work[j].value) <
                   lis_esolver_selected_real(work[best].value))
                    best = j;
            }

            if(best!=i)
            {
                LIS_ESOLVER_SELECTED_PAIR tmp = work[i];
                work[i] = work[best];
                work[best] = tmp;
            }
        }
    }

    if(available_count)
        *available_count = count;

    if(count<requested_count)
    {
        lis_free(work);
        return LIS_BREAKDOWN;
    }

    for(i=0;i<requested_count;i++)
        pairs[i] = work[i];

    lis_free(work);
    return LIS_SUCCESS;
}
