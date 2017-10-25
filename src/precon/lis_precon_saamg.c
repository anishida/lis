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

#if defined(USE_SAAMG)

extern void F77_FUNC_(finit_data_creation,FINIT_DATA_CREATION)(void *c_data_creation_ptr_bar);
extern void F77_FUNC_(finit_data_creation_unsym,FINIT_DATA_CREATION_UNSYM)(void *c_data_creation_unsym_ptr_bar);
extern void F77_FUNC_(finit_v_cycle,FINIT_V_CYCLE)(void *c_v_cycle_ptr_bar);
extern void F77_FUNC_(finit_clear_matrix,FINIT_CLEAR_MATRIX)(void *c_clear_matrix_ptr_bar);

char *f_data_creation_ptr;
char *f_data_creation_unsym_ptr;
char *f_v_cycle_ptr;
char *f_clear_matrix_ptr;
void c_data_creation_ptr_bar(char *bar)
{
	f_data_creation_ptr = bar;
}
void c_data_creation_unsym_ptr_bar(char *bar)
{
	f_data_creation_unsym_ptr = bar;
}
void c_v_cycle_ptr_bar(char *bar)
{
	f_v_cycle_ptr = bar;
}
void c_clear_matrix_ptr_bar(char *bar)
{
	f_clear_matrix_ptr = bar;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_precon_create_saamg"
LIS_INT lis_precon_create_saamg(LIS_SOLVER solver, LIS_PRECON precon)
{
#if defined(USE_SAAMG)
	LIS_MATRIX A,B;
	LIS_COMMTABLE table;
	LIS_INT	unsym,sol;
	LIS_INT	err;
	LIS_REAL theta; 
	#ifdef USE_MPI
		LIS_MPI_Fint comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	if( solver->A->matrix_type!=LIS_MATRIX_CSR )
	{
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		lis_matrix_destroy(B);
		solver->A = A;
	}
	precon->A       = solver->A;
	precon->is_copy = LIS_FALSE;
	A               = precon->A;
	sol             = solver->options[LIS_OPTIONS_SOLVER];
	unsym           = solver->options[LIS_OPTIONS_SAAMG_UNSYM];
	theta           = solver->params[LIS_PARAMS_SAAMG_THETA - LIS_OPTIONS_LEN];


#if 0
	if( sol!=LIS_SOLVER_CG && !unsym )
	{
		unsym = LIS_TRUE;
	}
#endif

	err = lis_vector_duplicate(A,&precon->temp);
	if( err )
	{
		return err;
	}
	F77_FUNC_(finit_data_creation,FINIT_DATA_CREATION)(c_data_creation_ptr_bar);
	F77_FUNC_(finit_data_creation_unsym,FINIT_DATA_CREATION_UNSYM)(c_data_creation_unsym_ptr_bar);
	F77_FUNC_(finit_v_cycle,FINIT_V_CYCLE)(c_v_cycle_ptr_bar);
	F77_FUNC_(finit_clear_matrix,FINIT_CLEAR_MATRIX)(c_clear_matrix_ptr_bar);
	
	lis_matrix_split(A);

	#ifdef USE_MPI
		comm = MPI_Comm_c2f(A->comm);
		lis_send_recv(A->commtable,A->D->value);
		table = A->commtable;
		if( !unsym )
		{
			(*(void (*)())f_data_creation_ptr)(&A->n,&A->np,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index,
				&table->neibpetot, table->neibpe, table->import_ptr,
				table->import_index, table->export_ptr, table->export_index,
				&table->imnnz,&table->exnnz,
				&comm, &precon->level_num,&precon->wsize, &theta);
		}
		else
		{
			(*(void (*)())f_data_creation_unsym_ptr)(&A->n,&A->np,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index,
				&table->neibpetot, table->neibpe, table->import_ptr,
				table->import_index, table->export_ptr, table->export_index,
				&table->imnnz,&table->exnnz,
				&comm, &precon->level_num,&precon->wsize, &theta);
		}
	#else
		if( !unsym )
		{
			(*(void (*)())f_data_creation_ptr)(&A->n,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index, &precon->level_num, &theta);
		}
		else
		{
			(*(void (*)())f_data_creation_unsym_ptr)(&A->n,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index, &precon->level_num, &theta);
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
#else
	LIS_DEBUG_FUNC_IN;

    precon->A       = solver->A;
    precon->is_copy = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
#endif
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_create_saamg"
LIS_INT lis_precon_psd_create_saamg(LIS_SOLVER solver, LIS_PRECON precon)
{
#if defined(USE_SAAMG)
	LIS_MATRIX A,B;
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	if( solver->A->matrix_type!=LIS_MATRIX_CSR )
	{
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		lis_matrix_destroy(B);
		solver->A = A;
	}
	precon->A       = solver->A;
	precon->is_copy = LIS_FALSE;
	A               = precon->A;


	err = lis_vector_duplicate(A,&precon->temp);
	if( err )
	{
		return err;
	}
	F77_FUNC_(finit_data_creation,FINIT_DATA_CREATION)(c_data_creation_ptr_bar);
	F77_FUNC_(finit_data_creation_unsym,FINIT_DATA_CREATION_UNSYM)(c_data_creation_unsym_ptr_bar);
	F77_FUNC_(finit_v_cycle,FINIT_V_CYCLE)(c_v_cycle_ptr_bar);
	F77_FUNC_(finit_clear_matrix,FINIT_CLEAR_MATRIX)(c_clear_matrix_ptr_bar);
	
    lis_matrix_split_create(A);

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
#else
	LIS_DEBUG_FUNC_IN;

    precon->A       = solver->A;
    precon->is_copy = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
#endif
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_update_saamg"
LIS_INT lis_precon_psd_update_saamg(LIS_SOLVER solver, LIS_PRECON precon)
{
#if defined(USE_SAAMG)
	LIS_MATRIX A;
	LIS_COMMTABLE table;
    LIS_INT	unsym,sol;
	LIS_INT	err;
	LIS_REAL theta; 
	#ifdef USE_MPI
		LIS_MPI_Fint comm;
	#endif

	LIS_DEBUG_FUNC_IN;

    A               = precon->A;
    sol             = solver->options[LIS_OPTIONS_SOLVER];
    unsym           = solver->options[LIS_OPTIONS_SAAMG_UNSYM];
    theta           = solver->params[LIS_PARAMS_SAAMG_THETA - LIS_OPTIONS_LEN];


#if 0
    if( sol!=LIS_SOLVER_CG && !unsym )
    {
        unsym = LIS_TRUE;
    }
#endif

    lis_matrix_split_update(A);

	#ifdef USE_MPI
		comm = MPI_Comm_c2f(A->comm);
		lis_send_recv(A->commtable,A->D->value); // this line should probably only be in the "update" version . . .
		table = A->commtable;
		if( !unsym )
		{
			(*(void (*)())f_data_creation_ptr)(&A->n,&A->np,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index,
				&table->neibpetot, table->neibpe, table->import_ptr,
				table->import_index, table->export_ptr, table->export_index,
				&table->imnnz,&table->exnnz,
				&comm, &precon->level_num,&precon->wsize, &theta);
		}
		else
		{
			(*(void (*)())f_data_creation_unsym_ptr)(&A->n,&A->np,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index,
				&table->neibpetot, table->neibpe, table->import_ptr,
				table->import_index, table->export_ptr, table->export_index,
				&table->imnnz,&table->exnnz,
				&comm, &precon->level_num,&precon->wsize, &theta);
		}
	#else
		if( !unsym )
		{
			(*(void (*)())f_data_creation_ptr)(&A->n,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index, &precon->level_num, &theta);
		}
		else
		{
			(*(void (*)())f_data_creation_unsym_ptr)(&A->n,&A->L->nnz,&A->U->nnz,
				A->D->value,A->L->value,A->L->ptr,A->L->index,
				A->U->value, A->U->ptr, A->U->index, &precon->level_num, &theta);
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_saamg"
LIS_INT lis_psolve_saamg(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
#if defined(USE_SAAMG)
	LIS_INT	n;
	LIS_PRECON precon;
	LIS_MATRIX A;
	#ifdef USE_MPI
		LIS_MPI_Fint comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	A      = solver->A;

	#ifdef USE_MPI
		comm = MPI_Comm_c2f(A->comm);
		n = b->np;
		 (*(void (*)())f_v_cycle_ptr)(b->value,x->value,precon->temp->value,&precon->level_num,
			 &comm,A->commtable->ws,A->commtable->wr,&n,&precon->wsize);
	#else
		n = b->n;
		 (*(void (*)())f_v_cycle_ptr)(&n,b->value,x->value,&precon->level_num,precon->temp->value);
	#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_DEBUG_FUNC_IN;

	lis_vector_copy(b,x);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_saamg"
LIS_INT lis_psolveh_saamg(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
#if defined(USE_SAAMG)
	LIS_INT	n;
	LIS_PRECON precon;
	LIS_MATRIX A;
	#ifdef USE_MPI
		LIS_MPI_Fint comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	A      = solver->A;
	precon = solver->precon;

	#ifdef USE_MPI
		comm = MPI_Comm_c2f(A->comm);
		n = b->np;
		 (*(void (*)())f_v_cycle_ptr)(b->value,x->value,precon->temp->value,&precon->level_num,
			 &comm,A->commtable->ws,A->commtable->wr,&n,&precon->wsize);
	#else
		n = b->n;
		 (*(void (*)())f_v_cycle_ptr)(&n,b->value,x->value,&precon->level_num,precon->temp->value);
	#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	LIS_DEBUG_FUNC_IN;

	lis_vector_copy(b,x);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}
