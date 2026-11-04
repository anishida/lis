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

/*
 * Atomic export implementation shared by the legacy and tolerance-aware
 * APIs.  The first pass validates the complete vector; destination bytes
 * are touched only after every check has passed.
 */
static LIS_INT
lis_vector_export_real64_impl(
    LIS_VECTOR v,double *dst,LIS_INT length,LIS_INT stride,
    LIS_INT loss_policy,LIS_INT check_imag_tolerance,
    LIS_REAL absolute_imag_tolerance,LIS_REAL relative_imag_tolerance)
{
    LIS_INT i,err=lis_vector_real64_validate(v,length,stride);
    LIS_REAL max_real=(LIS_REAL)0,max_imag=(LIS_REAL)0;

    if(err ||
       (loss_policy!=LIS_VECTOR_REAL64_REJECT_LOSS &&
        loss_policy!=LIS_VECTOR_REAL64_ALLOW_LOSS) ||
       (check_imag_tolerance &&
        (!isfinite((double)relative_imag_tolerance) ||
         !isfinite((double)absolute_imag_tolerance) ||
         relative_imag_tolerance<(LIS_REAL)0 ||
         absolute_imag_tolerance<(LIS_REAL)0)))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid real64 export layout, policy or tolerances\n");
        return LIS_ERR_ILL_ARG;
    }
    if(length==0) return LIS_SUCCESS;
    if(dst==NULL)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"real64 export destination is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<length;i++)
    {
        LIS_REAL real_value,imag_value=(LIS_REAL)0,abs_real,abs_imag;
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
        real_value=v->value[i][0];
        imag_value=v->value[i][1];
#elif defined(_COMPLEX)
        real_value=(LIS_REAL)creal(v->value[i]);
        imag_value=(LIS_REAL)cimag(v->value[i]);
#else
        real_value=(LIS_REAL)v->value[i];
#endif
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
        if(check_imag_tolerance)
        {
            abs_real=(LIS_REAL)fabs(real_value);
            abs_imag=(LIS_REAL)fabs(imag_value);
            if(abs_real>max_real) max_real=abs_real;
            if(abs_imag>max_imag) max_imag=abs_imag;
        }
    }

    /*
     * Compare maxima over the LOCAL vector (v->n); no MPI collective is\n     * performed. Compare global-to-this-buffer maxima, not per-element\n     * ratios. Use
     * subtraction/division to avoid overflow in atol+rtol*max_real.
     */
    if(check_imag_tolerance &&
       max_imag>absolute_imag_tolerance &&
       (max_real==(LIS_REAL)0 ||
        (max_imag-absolute_imag_tolerance)/max_real>
            relative_imag_tolerance))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "real64 export maximum imaginary magnitude exceeds tolerance\n");
        return LIS_ERR_ILL_ARG;
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

#undef __FUNC__
#define __FUNC__ "lis_vector_export_real64"
LIS_INT
lis_vector_export_real64(
    LIS_VECTOR v,double *dst,LIS_INT length,LIS_INT stride,
    LIS_INT loss_policy)
{
    return lis_vector_export_real64_impl(
        v,dst,length,stride,loss_policy,LIS_FALSE,
        (LIS_REAL)0,(LIS_REAL)0);
}

#undef __FUNC__
#define __FUNC__ "lis_vector_export_real64_tolerant"
LIS_INT
lis_vector_export_real64_tolerant(
    LIS_VECTOR v,double *dst,LIS_INT length,LIS_INT stride,
    LIS_INT loss_policy,LIS_REAL absolute_imag_tolerance,
    LIS_REAL relative_imag_tolerance)
{
    return lis_vector_export_real64_impl(
        v,dst,length,stride,loss_policy,LIS_TRUE,
        absolute_imag_tolerance,relative_imag_tolerance);
}
