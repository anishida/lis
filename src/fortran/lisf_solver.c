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
 * lis_solver_create
 * lis_solver_destroy
 * lis_solver_set_option
 * lis_solver_get_option
 * lis_solve
 ************************************************/

#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_solver_create_f"
void lis_solver_create_f(LIS_SOLVER_F *solver, LIS_INT *ierr)
{
	LIS_SOLVER s;

	LIS_DEBUG_FUNC_IN;

    *ierr = lis_solver_create(&s);
	if( *ierr )	return;

	*solver = LIS_P2V(s);
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_destroy_f"
void lis_solver_destroy_f(LIS_SOLVER_F *solver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_destroy((LIS_SOLVER)LIS_V2P(solver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_solver_set_matrix_f"
void lis_solver_set_matrix_f(LIS_MATRIX_F *A, LIS_SOLVER_F *solver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_set_matrix((LIS_MATRIX)LIS_V2P(A), (LIS_SOLVER)LIS_V2P(solver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solve_f"
void lis_solve_f(LIS_MATRIX_F *A, LIS_VECTOR_F b, LIS_VECTOR_F x, LIS_SOLVER_F *solver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solve((LIS_MATRIX)LIS_V2P(A), (LIS_VECTOR)LIS_V2P(b), (LIS_VECTOR)LIS_V2P(x), (LIS_SOLVER)LIS_V2P(solver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_f"
void lis_solver_set_option_f(char *text, LIS_SOLVER_F *solver, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	LIS_DEBUG_FUNC_IN;

	strncpy(buf,text,len);
	buf[len] = '\0';
	*ierr = lis_solver_set_option(buf,(LIS_SOLVER)LIS_V2P(solver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_optionC_f"
void lis_solver_set_optionC_f(LIS_SOLVER solver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_set_optionC((LIS_SOLVER)LIS_V2P(solver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_iter_f"
void lis_solver_get_iter_f(LIS_SOLVER_F *solver, LIS_INT *iter, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_iter((LIS_SOLVER)LIS_V2P(solver),iter);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_iterex_f"
void lis_solver_get_iterex_f(LIS_SOLVER_F *solver, LIS_INT *iter, LIS_INT *iter_double, LIS_INT *iter_quad, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_iterex((LIS_SOLVER)LIS_V2P(solver),iter,iter_double,iter_quad);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_time_f"
void lis_solver_get_time_f(LIS_SOLVER_F *solver, double *time, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_time((LIS_SOLVER)LIS_V2P(solver),time);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_timeex_f"
void lis_solver_get_timeex_f(LIS_SOLVER_F *solver, double *time, double *itime, double *ptime, double *p_c_time, double *p_i_time, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_timeex((LIS_SOLVER)LIS_V2P(solver),time,itime,ptime,p_c_time,p_i_time);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_residualnorm_f"
void lis_solver_get_residualnorm_f(LIS_SOLVER_F *solver, LIS_REAL *residual, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_residualnorm((LIS_SOLVER)LIS_V2P(solver),residual);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_solver_f"
void lis_solver_get_solver_f(LIS_SOLVER_F *solver, LIS_INT *nsol, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_solver((LIS_SOLVER)LIS_V2P(solver),nsol);
	if( *ierr )	return;
	
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_rhistory_f"
void lis_solver_get_rhistory_f(LIS_SOLVER_F *solver, LIS_VECTOR_F v, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_rhistory((LIS_SOLVER)LIS_V2P(solver),(LIS_VECTOR)LIS_V2P(v));
	if( *ierr )	return;
	
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_solvername_f"
void lis_solver_get_solvername_f(LIS_INT *solver, char *name, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_solvername(*solver, buf);
	if( *ierr )	return;	
	strncpy(name,buf,len);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solve_kernel_f"
void lis_solve_kernel_f(LIS_MATRIX_F *A, LIS_VECTOR_F b, LIS_VECTOR_F x, LIS_SOLVER_F *solver, LIS_PRECON *precon, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solve_kernel((LIS_MATRIX)LIS_V2P(A), (LIS_VECTOR)LIS_V2P(b), (LIS_VECTOR)LIS_V2P(x), (LIS_SOLVER)LIS_V2P(solver), (LIS_PRECON)LIS_V2P(precon));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solve_setup_f"
void lis_solve_setup_f(LIS_MATRIX_F *A, LIS_SOLVER_F *solver, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solve_setup((LIS_MATRIX)LIS_V2P(A), (LIS_SOLVER)LIS_V2P(solver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_precon_f"
void lis_solver_get_precon_f(LIS_SOLVER_F *solver, LIS_INT *nsol, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_precon((LIS_SOLVER)LIS_V2P(solver),nsol);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_preconname_f"
void lis_solver_get_preconname_f(LIS_INT *solver, char *name, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_preconname(*solver, buf);
	if( *ierr )	return;	
	strncpy(name,buf,len);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_status_f"
void lis_solver_get_status_f(LIS_SOLVER_F *solver, LIS_INT *status, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_solver_get_status((LIS_SOLVER)LIS_V2P(solver),status);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif
