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
 * lis_precon_create
 * lis_precon_psd_create
 * lis_precon_psd_update
 * lis_precon_destroy
 * lis_precon_set_option
 * lis_precon_get_option
 * lis_solve
 ************************************************/

#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_precon_create_f"
void lis_precon_create_f(LIS_SOLVER_F solver, LIS_PRECON_F *precon, LIS_INT *ierr)
{
	LIS_PRECON p;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_precon_create((LIS_SOLVER)LIS_V2P(solver), &p);
	if( *ierr )	return;

	*precon = LIS_P2V(p);
	LIS_DEBUG_FUNC_OUT;
	return;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_create_f"
void lis_precon_psd_create_f(LIS_SOLVER_F solver, LIS_PRECON_F *precon, LIS_INT *ierr)
{
	LIS_PRECON p;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_precon_psd_create((LIS_SOLVER)LIS_V2P(solver), &p);
	if( *ierr )	return;

	*precon = LIS_P2V(p);
	LIS_DEBUG_FUNC_OUT;
	return;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_update_f"
void lis_precon_psd_update_f(LIS_SOLVER_F solver, LIS_PRECON_F precon, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_precon_psd_update((LIS_SOLVER)LIS_V2P(solver),(LIS_PRECON)LIS_V2P(precon));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_destroy_f"
void lis_precon_destroy_f(LIS_PRECON_F *precon, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_precon_destroy((LIS_PRECON)LIS_V2P(precon));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif
