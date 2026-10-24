/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "lislib.h"

/*
 * Native FEAST eigensolver foundation.
 *
 * FEAST 4.0 is used as an algorithmic reference.
 * This implementation has no external FEAST library dependency.
 *
 * Stage F0 provides registration, interval configuration and
 * parameter validation. The contour projector and shifted systems
 * are implemented in later stages.
 */

#define LIS_FEAST_INTERVAL_EMIN 0x01
#define LIS_FEAST_INTERVAL_EMAX 0x02

#define LIS_FEAST_PI \
        ((LIS_REAL)3.141592653589793238462643383279502884)


/*
 * Compute the n-point Gauss-Legendre rule on [-1,1].
 *
 * The implementation uses Newton refinement of the Legendre roots
 * and does not depend on a table of precomputed quadrature nodes.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_gauss_legendre"

LIS_INT lis_efeast_gauss_legendre(
        LIS_INT n,
        LIS_REAL *x,
        LIS_REAL *w)
{
        LIS_INT i;
        LIS_INT j;
        LIS_INT m;
        LIS_INT iter;

        LIS_REAL z;
        LIS_REAL zold;
        LIS_REAL p1;
        LIS_REAL p2;
        LIS_REAL p3;
        LIS_REAL pp;
        LIS_REAL wi;
        LIS_REAL tol;

        LIS_DEBUG_FUNC_IN;

        if( n<1 || x==NULL || w==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Gauss-Legendre quadrature requires n >= 1 and valid arrays\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        m = (n+1)/2;

        tol = (LIS_REAL)(64.0*DBL_EPSILON);

        for(i=0;i<m;i++)
        {
                /*
                 * Standard asymptotic initial estimate for the
                 * positive roots.
                 */
                z = cos(
                        LIS_FEAST_PI *
                        ((LIS_REAL)i + (LIS_REAL)0.75) /
                        ((LIS_REAL)n + (LIS_REAL)0.5));

                for(iter=0;iter<100;iter++)
                {
                        p1 = 1.0;
                        p2 = 0.0;

                        for(j=1;j<=n;j++)
                        {
                                p3 = p2;
                                p2 = p1;

                                p1 =
                                        (((LIS_REAL)(2*j-1))*z*p2 -
                                         ((LIS_REAL)(j-1))*p3) /
                                        (LIS_REAL)j;
                        }

                        pp =
                                ((LIS_REAL)n) *
                                (z*p1-p2) /
                                (z*z-(LIS_REAL)1.0);

                        zold = z;
                        z = zold-p1/pp;

                        if( fabs(z-zold) <=
                            tol*((LIS_REAL)1.0+fabs(z)) )
                        {
                                break;
                        }
                }

                if( iter==100 )
                {
                        LIS_SETERR1(
                                LIS_BREAKDOWN,
                                "Gauss-Legendre root %D did not converge\n",
                                i);

                        LIS_DEBUG_FUNC_OUT;
                        return LIS_BREAKDOWN;
                }

                /*
                 * Re-evaluate P_n and its derivative at the
                 * converged root before forming the weight.
                 */
                p1 = 1.0;
                p2 = 0.0;

                for(j=1;j<=n;j++)
                {
                        p3 = p2;
                        p2 = p1;

                        p1 =
                                (((LIS_REAL)(2*j-1))*z*p2 -
                                 ((LIS_REAL)(j-1))*p3) /
                                (LIS_REAL)j;
                }

                pp =
                        ((LIS_REAL)n) *
                        (z*p1-p2) /
                        (z*z-(LIS_REAL)1.0);

                wi =
                        ((LIS_REAL)2.0) /
                        (((LIS_REAL)1.0-z*z)*pp*pp);

                x[i] = -z;
                x[n-1-i] = z;

                w[i] = wi;
                w[n-1-i] = wi;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


/*
 * Generate the upper-half elliptic contour used by the native
 * real-symmetric FEAST path.
 *
 * The complex values are represented as separate real and imaginary
 * arrays so this geometry kernel is available in both real and
 * complex LIS builds.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_contour_gauss"

LIS_INT lis_efeast_contour_gauss(
        LIS_REAL emin,
        LIS_REAL emax,
        LIS_INT n,
        LIS_REAL ratio,
        LIS_REAL *z_real,
        LIS_REAL *z_imag,
        LIS_REAL *w_real,
        LIS_REAL *w_imag)
{
        LIS_REAL x[1024];
        LIS_REAL w[1024];

        LIS_REAL radius;
        LIS_REAL center;
        LIS_REAL theta;
        LIS_REAL ctheta;
        LIS_REAL stheta;

        LIS_INT i;
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        if( !(emin<emax) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour requires emin < emax\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( n<1 || n>1024 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour point count must be in [1,1024]\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !(ratio>0.0) || !isfinite((double)ratio) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour ratio must be finite and positive\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( z_real==NULL ||
            z_imag==NULL ||
            w_real==NULL ||
            w_imag==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour output arrays are undefined\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        err = lis_efeast_gauss_legendre(
                n,
                x,
                w);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        radius =
                (emax-emin) /
                (LIS_REAL)2.0;

        center = emin+radius;

        for(i=0;i<n;i++)
        {
                /*
                 * Map [-1,1] to the upper half ellipse:
                 *
                 * theta in [pi,0].
                 */
                theta =
                        -(LIS_FEAST_PI/(LIS_REAL)2.0)*x[i] +
                        LIS_FEAST_PI/(LIS_REAL)2.0;

                ctheta = cos(theta);
                stheta = sin(theta);

                z_real[i] =
                        center +
                        radius*ctheta;

                z_imag[i] =
                        radius*ratio*stheta;

                /*
                 * FEAST projector weight:
                 *
                 *   (1/4) * w_i * dz/dtheta
                 *
                 * after the half-contour parameter mapping.
                 */
                w_real[i] =
                        (LIS_REAL)0.25 *
                        w[i] *
                        radius *
                        ratio *
                        ctheta;

                w_imag[i] =
                        (LIS_REAL)0.25 *
                        w[i] *
                        radius *
                        stheta;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


/*
 * Apply the complex shifted operator
 *
 *     y = (z B - A) x
 *
 * for real A and B using split-complex vectors:
 *
 *     x = x_real + i x_imag
 *     z = z_real + i z_imag
 *
 * If B is NULL, the identity matrix is used.
 *
 * This representation is intentionally available in both real
 * and complex LIS builds.  It provides the foundation for the
 * real 2N block formulation required by applications that keep
 * the original matrices real.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_shifted_matvec_split"

LIS_INT lis_efeast_shifted_matvec_split(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_REAL z_real,
        LIS_REAL z_imag,
        LIS_VECTOR x_real,
        LIS_VECTOR x_imag,
        LIS_VECTOR y_real,
        LIS_VECTOR y_imag,
        LIS_VECTOR work_real,
        LIS_VECTOR work_imag)
{
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        if( A==NULL ||
            x_real==NULL ||
            x_imag==NULL ||
            y_real==NULL ||
            y_imag==NULL ||
            work_real==NULL ||
            work_imag==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted operator received an undefined argument\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !isfinite((double)z_real) ||
            !isfinite((double)z_imag) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted operator requires a finite complex shift\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        /*
         * y_real = -A*x_real
         */
        err = lis_matvec(
                A,
                x_real,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_scale(
                (LIS_SCALAR)-1.0,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        /*
         * y_imag = -A*x_imag
         */
        err = lis_matvec(
                A,
                x_imag,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_scale(
                (LIS_SCALAR)-1.0,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        /*
         * work_real = B*x_real
         * work_imag = B*x_imag
         *
         * B == NULL means B = I.
         */
        if( B!=NULL )
        {
                err = lis_matvec(
                        B,
                        x_real,
                        work_real);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }

                err = lis_matvec(
                        B,
                        x_imag,
                        work_imag);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }
        }
        else
        {
                err = lis_vector_copy(
                        x_real,
                        work_real);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }

                err = lis_vector_copy(
                        x_imag,
                        work_imag);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }
        }

        /*
         * Real part:
         *
         * y_real =
         *     -A*x_real
         *     + z_real*B*x_real
         *     - z_imag*B*x_imag
         */
        err = lis_vector_axpy(
                (LIS_SCALAR)z_real,
                work_real,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_axpy(
                (LIS_SCALAR)(-z_imag),
                work_imag,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        /*
         * Imaginary part:
         *
         * y_imag =
         *     -A*x_imag
         *     + z_imag*B*x_real
         *     + z_real*B*x_imag
         */
        err = lis_vector_axpy(
                (LIS_SCALAR)z_imag,
                work_real,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_axpy(
                (LIS_SCALAR)z_real,
                work_imag,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


static LIS_INT lis_efeast_parse_real(
        const char *text,
        LIS_REAL *value)
{
        char *endptr;
        double parsed;

        if( text==NULL || value==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        endptr = NULL;
        parsed = strtod(text,&endptr);

        if( endptr==text ||
            endptr==NULL ||
            *endptr!='\0' ||
            !isfinite(parsed) )
        {
                return LIS_ERR_ILL_ARG;
        }

        *value = (LIS_REAL)parsed;

        return LIS_SUCCESS;
}


static LIS_INT lis_efeast_parse_int(
        const char *text,
        LIS_INT *value)
{
        char *endptr;
        long parsed;

        if( text==NULL || value==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        endptr = NULL;
        parsed = strtol(text,&endptr,10);

        if( endptr==text ||
            endptr==NULL ||
            *endptr!='\0' ||
            parsed<1 ||
            parsed>1024 )
        {
                return LIS_ERR_ILL_ARG;
        }

        *value = (LIS_INT)parsed;

        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast_set_option"

LIS_INT lis_efeast_set_option(
        char *arg1,
        char *arg2,
        LIS_ESOLVER esolver)
{
        LIS_REAL value;
        LIS_INT ivalue;
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL ||
            arg1==NULL ||
            arg2==NULL )
        {
                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( strcmp(arg1,"-feast_emin")==0 )
        {
                err = lis_efeast_parse_real(arg2,&value);

                if( err )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_emin must be a finite real value\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_emin = value;
                esolver->feast_interval_mask |=
                        LIS_FEAST_INTERVAL_EMIN;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( strcmp(arg1,"-feast_emax")==0 )
        {
                err = lis_efeast_parse_real(arg2,&value);

                if( err )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_emax must be a finite real value\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_emax = value;
                esolver->feast_interval_mask |=
                        LIS_FEAST_INTERVAL_EMAX;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( strcmp(arg1,"-feast_ratio")==0 )
        {
                err = lis_efeast_parse_real(arg2,&value);

                if( err ||
                    !(value>0.0) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_ratio must be finite and positive\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_ratio = value;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( strcmp(arg1,"-feast_points")==0 )
        {
                err = lis_efeast_parse_int(arg2,&ivalue);

                if( err )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_points must be an integer in [1,1024]\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_points = ivalue;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        LIS_SETERR1(
                LIS_ERR_ILL_ARG,
                "Unknown native FEAST option %s\n",
                arg1);

        LIS_DEBUG_FUNC_OUT;
        return LIS_ERR_ILL_ARG;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast_check_params"

LIS_INT lis_efeast_check_params(
        LIS_ESOLVER esolver)
{
        LIS_INT ss;

        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL )
        {
                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

        if( ss<=0 )
        {
                LIS_SETERR1(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST requires positive subspace size -ss (got %D)\n",
                        ss);

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( esolver->feast_interval_mask !=
            (LIS_FEAST_INTERVAL_EMIN |
             LIS_FEAST_INTERVAL_EMAX) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST requires both -feast_emin and -feast_emax\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !(esolver->feast_emin <
              esolver->feast_emax) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST requires -feast_emin < -feast_emax\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !(esolver->feast_ratio>0.0) ||
            !isfinite((double)esolver->feast_ratio) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST contour ratio must be finite and positive\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( esolver->feast_points<1 ||
            esolver->feast_points>1024 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST contour point count is invalid\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast_malloc_work"

LIS_INT lis_efeast_malloc_work(
        LIS_ESOLVER esolver)
{
        LIS_DEBUG_FUNC_IN;

        esolver->work = NULL;
        esolver->worklen = 0;

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast"

LIS_INT lis_efeast(
        LIS_ESOLVER esolver)
{
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        err = lis_efeast_check_params(esolver);

        if( err )
        {
                esolver->retcode = err;
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

#ifndef _COMPLEX

        LIS_SETERR(
                LIS_ERR_NOT_IMPLEMENTED,
                "Native FEAST requires complex shifted systems; "
                "the real-build bridge is not implemented yet\n");

#else

        LIS_SETERR(
                LIS_ERR_NOT_IMPLEMENTED,
                "Native FEAST contour projector is not implemented in stage F0\n");

#endif

        esolver->retcode = LIS_ERR_NOT_IMPLEMENTED;

        LIS_DEBUG_FUNC_OUT;
        return LIS_ERR_NOT_IMPLEMENTED;
}
