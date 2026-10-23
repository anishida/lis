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
