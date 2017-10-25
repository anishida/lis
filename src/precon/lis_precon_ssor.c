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
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_precon_create
 * lis_psolve
 * lis_psolveh
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_precon_create_ssor"
LIS_INT lis_precon_create_ssor(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	err;
	LIS_SCALAR w;
	LIS_MATRIX A;

	LIS_DEBUG_FUNC_IN;

	A   = solver->A;
	w   = solver->params[LIS_PARAMS_SSOR_OMEGA-LIS_OPTIONS_LEN];

	err = lis_matrix_convert_self(solver);
	if( err ) return err;

	err = lis_matrix_split(A);
	if( err )
	{
		return err;
	}
	if( A->use_wd!=LIS_SOLVER_SOR )
	{
		if( !A->WD )
		{
			err = lis_matrix_diag_duplicate(A->D,&A->WD);
			if( err ) return err;
		}
		lis_matrix_diag_copy(A->D,A->WD);
		lis_matrix_diag_scale(w,A->WD);
		lis_matrix_diag_inverse(A->WD);
		A->use_wd = LIS_SOLVER_SOR;
	}

	precon->A       = A;
	precon->is_copy = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_ssor"
LIS_INT lis_psolve_ssor(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_MATRIX A;

	/*
	 *  Mx = b
	 *  M  = (D/w + L) * (I + w*D^-1 * U)
	 */

	LIS_DEBUG_FUNC_IN;

	A = solver->precon->A;

	lis_matrix_solve(A,B,X,LIS_MATRIX_SSOR);


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_ssor"
LIS_INT lis_psolveh_ssor(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_MATRIX A;

	/*
	 *  M'x = b
	 *  M'  = (I + U' * w*D^-1) * (D/w + L')
	 */

	LIS_DEBUG_FUNC_IN;

	A = solver->precon->A;

	lis_matrix_solveh(A,B,X,LIS_MATRIX_SSOR);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
