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
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_esolver_create
 * lis_esolver_destroy
 * lis_esolver_set_option
 * lis_esolver_get_option
 * lis_esolve
 * lis_gesolve
 ************************************************/

#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_esolver_create_f"
void lis_esolver_create_f(LIS_ESOLVER_F *esolver, LIS_INT *ierr)
{
	LIS_ESOLVER	s;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_create(&s);
	if( *ierr )	return;

	*esolver = LIS_P2V(s);
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_destroy_f"
void lis_esolver_destroy_f(LIS_ESOLVER_F *esolver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_destroy((LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolve_f"
void lis_esolve_f(LIS_MATRIX_F *A, LIS_VECTOR_F x, LIS_SCALAR *evalue0, LIS_ESOLVER_F *esolver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolve((LIS_MATRIX)LIS_V2P(A), (LIS_VECTOR)LIS_V2P(x), evalue0, (LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_gesolve_f"
void lis_gesolve_f(LIS_MATRIX_F *A, LIS_MATRIX_F *B, LIS_VECTOR_F x, LIS_SCALAR *evalue0, LIS_ESOLVER_F *esolver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_gesolve((LIS_MATRIX)LIS_V2P(A), (LIS_MATRIX)LIS_V2P(B), (LIS_VECTOR)LIS_V2P(x), evalue0, (LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_f"
void lis_esolver_set_option_f(char *text, LIS_ESOLVER_F *esolver, LIS_INT *ierr, LIS_INT len)
{
	char	buf[1024];
	LIS_DEBUG_FUNC_IN;

	strncpy(buf,text,len);
	buf[len] = '\0';
	*ierr = lis_esolver_set_option(buf,(LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_optionC_f"
void lis_esolver_set_optionC_f(LIS_ESOLVER esolver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_set_optionC((LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;
	
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iter_f"
void lis_esolver_get_iter_f(LIS_ESOLVER_F *esolver, LIS_INT *iter, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_iter((LIS_ESOLVER)LIS_V2P(esolver),iter);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iterex_f"
void lis_esolver_get_iterex_f(LIS_ESOLVER_F *esolver, LIS_INT *iter, LIS_INT *iter_double, LIS_INT *iter_quad, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_iterex((LIS_ESOLVER)LIS_V2P(esolver),iter,iter_double,iter_quad);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_time_f"
void lis_esolver_get_time_f(LIS_ESOLVER_F *esolver, double *time, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_time((LIS_ESOLVER)LIS_V2P(esolver),time);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_timeex_f"
void lis_esolver_get_timeex_f(LIS_ESOLVER_F *esolver, double *time, double *itime, double *ptime, double *p_c_time, double *p_i_time, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_timeex((LIS_ESOLVER)LIS_V2P(esolver),time,itime,ptime,p_c_time,p_i_time);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_residualnorm_f"
void lis_esolver_get_residualnorm_f(LIS_ESOLVER_F *esolver, LIS_REAL *residual, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_residualnorm((LIS_ESOLVER)LIS_V2P(esolver),residual);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_rhistory_f"
void lis_esolver_get_rhistory_f(LIS_SOLVER_F *esolver, LIS_VECTOR_F v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_rhistory((LIS_ESOLVER)LIS_V2P(esolver),(LIS_VECTOR)LIS_V2P(v));
	if( *ierr )	return;
	
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evalues_f"
void lis_esolver_get_evalues_f(LIS_ESOLVER_F *esolver, LIS_VECTOR v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_evalues((LIS_ESOLVER)LIS_V2P(esolver),(LIS_VECTOR)LIS_V2P(v));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_evalue_f"
void lis_esolver_get_specific_evalue_f(LIS_ESOLVER_F *esolver, LIS_INT *mode, LIS_SCALAR *evalue, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_specific_evalue((LIS_ESOLVER)LIS_V2P(esolver),*mode,evalue);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evectors_f"
void lis_esolver_get_evectors_f(LIS_ESOLVER_F *esolver, LIS_MATRIX_F A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_evectors((LIS_ESOLVER)LIS_V2P(esolver),(LIS_MATRIX)LIS_V2P(A));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_evector_f"
void lis_esolver_get_specific_evector_f(LIS_ESOLVER_F *esolver, LIS_INT *mode, LIS_VECTOR x, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_specific_evector((LIS_ESOLVER)LIS_V2P(esolver),*mode,(LIS_VECTOR)LIS_V2P(x));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_residualnorms_f"
void lis_esolver_get_residualnorms_f(LIS_ESOLVER_F *esolver, LIS_VECTOR v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_residualnorms((LIS_ESOLVER)LIS_V2P(esolver),(LIS_VECTOR)LIS_V2P(v));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_residualnorm_f"
void lis_esolver_get_specific_residualnorm_f(LIS_ESOLVER_F *esolver, LIS_INT *mode, LIS_REAL *residual, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_specific_residualnorm((LIS_ESOLVER)LIS_V2P(esolver),*mode,residual);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iters_f"
void lis_esolver_get_iters_f(LIS_ESOLVER_F *esolver, LIS_VECTOR v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_iters((LIS_ESOLVER)LIS_V2P(esolver),(LIS_VECTOR)LIS_V2P(v));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_specific_iter_f"
void lis_esolver_get_specific_iter_f(LIS_ESOLVER_F *esolver, LIS_INT *mode, LIS_INT *iter, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_specific_iter((LIS_ESOLVER)LIS_V2P(esolver),*mode,iter);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_esolver_f"
void lis_esolver_get_esolver_f(LIS_ESOLVER_F *esolver, LIS_INT *nsol, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_esolver((LIS_ESOLVER)LIS_V2P(esolver),nsol);
	if( *ierr )	return;
	
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_esolvername_f"
void lis_esolver_get_esolvername_f(LIS_INT *esolver, char *name, LIS_INT *ierr, LIS_INT len)
{
	char	buf[1024];
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_esolvername(*esolver, buf);
	if( *ierr )	return;
	strncpy(name,buf,len);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif
