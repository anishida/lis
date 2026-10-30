/* Public application real64 vector transfer helpers. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif
#include <math.h>
#include <stdint.h>
#if defined(_COMPLEX) && defined(HAVE_COMPLEX_H)
#include <complex.h>
#endif
#include "lislib.h"

static LIS_INT
lis_vector_real64_validate(LIS_VECTOR v,LIS_INT length,LIS_INT stride)
{
    if(v==NULL || v->precision!=LIS_PRECISION_DEFAULT ||
       length<0 || length!=v->n || stride<1 ||
       (length>0 && v->value==NULL))
        return LIS_ERR_ILL_ARG;
    if(length>1 && (size_t)(length-1)>SIZE_MAX/(size_t)stride)
        return LIS_ERR_ILL_ARG;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_import_real64"
LIS_INT
lis_vector_import_real64(
    LIS_VECTOR v,const double *src,LIS_INT length,LIS_INT stride)
{
    LIS_INT i,err=lis_vector_real64_validate(v,length,stride);
    if(err)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid real64 import vector or stride\n");
        return err;
    }
    if(length==0) return LIS_SUCCESS;
    if(src==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"real64 import source is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    /* Check the entire source first so failed conversions are atomic. */
    for(i=0;i<length;i++)
    {
        double x=src[(size_t)i*(size_t)stride];
        LIS_REAL y=(LIS_REAL)x;
        if(!isfinite(x) || !isfinite((double)y))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,"real64 import requires finite representable values\n");
            return LIS_ERR_ILL_ARG;
        }
    }
    for(i=0;i<length;i++)
    {
        LIS_REAL x=(LIS_REAL)src[(size_t)i*(size_t)stride];
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        v->value[i][0]=x;
        v->value[i][1]=(LIS_REAL)0;
#else
        v->value[i]=(LIS_SCALAR)x;
#endif
    }
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_export_real64"
LIS_INT
lis_vector_export_real64(
    LIS_VECTOR v,double *dst,LIS_INT length,LIS_INT stride,
    LIS_INT loss_policy)
{
    LIS_INT i,err=lis_vector_real64_validate(v,length,stride);
    if(err ||
       (loss_policy!=LIS_VECTOR_REAL64_REJECT_LOSS &&
        loss_policy!=LIS_VECTOR_REAL64_ALLOW_LOSS))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid real64 export layout or policy\n");
        return LIS_ERR_ILL_ARG;
    }
    if(length==0) return LIS_SUCCESS;
    if(dst==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"real64 export destination is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    /* Validate first: a failed strict export must not partially write dst. */
    for(i=0;i<length;i++)
    {
        LIS_REAL real_value,imag_value=(LIS_REAL)0;
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        real_value=v->value[i][0];
        imag_value=v->value[i][1];
#elif defined(_COMPLEX)
        real_value=(LIS_REAL)creal(v->value[i]);
        imag_value=(LIS_REAL)cimag(v->value[i]);
#else
        real_value=(LIS_REAL)v->value[i];
#endif
        /* Reject nonfinite imaginary parts even when ALLOW_LOSS ignores them. */
        if(!isfinite((double)real_value) ||
           !isfinite((double)imag_value))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,"real64 export requires finite values\n");
            return LIS_ERR_ILL_ARG;
        }
        if(loss_policy==LIS_VECTOR_REAL64_REJECT_LOSS &&
           (imag_value!=(LIS_REAL)0 ||
            (LIS_REAL)(double)real_value!=real_value))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,"strict real64 export would lose information\n");
            return LIS_ERR_ILL_ARG;
        }
    }
    for(i=0;i<length;i++)
    {
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        dst[(size_t)i*(size_t)stride]=(double)v->value[i][0];
#elif defined(_COMPLEX)
        dst[(size_t)i*(size_t)stride]=(double)creal(v->value[i]);
#else
        dst[(size_t)i*(size_t)stride]=(double)v->value[i];
#endif
    }
    return LIS_SUCCESS;
}
