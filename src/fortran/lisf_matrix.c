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
 * lis_matrix_create_f
 * lis_matrix_duplicate_f
 * lis_matrix_destroy_f
 ************************************************/

#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_f"
void lis_matrix_create_f(LIS_Comm_f *comm, LIS_MATRIX_F *Amat, LIS_INT *ierr)
{
	LIS_MATRIX A;
	LIS_Comm c_comm;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		if( *comm==lis_comm_world_f )
		{
			c_comm = MPI_COMM_WORLD;
		}
		else
		{
			c_comm = MPI_Comm_f2c(*comm);
		}
	#else
		c_comm = *comm;
	#endif

	*ierr = lis_matrix_create(c_comm,&A);
	if( *ierr )	return;

	A->origin = LIS_ORIGIN_1;
	*Amat = LIS_P2V(A);
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_size_f"
void lis_matrix_set_size_f(LIS_MATRIX_F *A, LIS_INT *local_n, LIS_INT *global_n, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_size((LIS_MATRIX)LIS_V2P(A),*local_n,*global_n);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_duplicate_f"
void lis_matrix_duplicate_f(LIS_MATRIX_F *Ain, LIS_MATRIX_F *Aout, LIS_INT *ierr)
{
	LIS_MATRIX A;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_duplicate((LIS_MATRIX)LIS_V2P(Ain),&A);
	if( *ierr )	return;	
	*Aout = LIS_P2V(A);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_destroy_f"
void lis_matrix_destroy_f(LIS_MATRIX_F *Amat, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_destroy((LIS_MATRIX)LIS_V2P(Amat));
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_range_f"
void lis_matrix_get_range_f(LIS_MATRIX_F *A, LIS_INT *is, LIS_INT *ie, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_get_range((LIS_MATRIX)LIS_V2P(A),is,ie);
	if( *ierr )	return;	
	(*is)++;
	(*ie)++;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_size_f"
void lis_matrix_get_size_f(LIS_MATRIX_F *A, LIS_INT *local_n, LIS_INT *global_n, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_get_size((LIS_MATRIX)LIS_V2P(A),local_n,global_n);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_nnz_f"
void lis_matrix_get_nnz_f(LIS_MATRIX_F *A, LIS_INT *nnz, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_get_nnz((LIS_MATRIX)LIS_V2P(A),nnz);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_value_f"
void lis_matrix_set_value_f(LIS_INT *flag, LIS_INT *i, LIS_INT *j, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_value(*flag,*i,*j,*value,(LIS_MATRIX)LIS_V2P(A));
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_psd_set_value_f"
void lis_matrix_psd_set_value_f(LIS_INT *flag, LIS_INT *i, LIS_INT *j, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_psd_set_value(*flag,*i,*j,*value,(LIS_MATRIX)LIS_V2P(A));
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_type_f"
void lis_matrix_set_type_f(LIS_MATRIX_F *A, LIS_INT *matrix_type, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_type((LIS_MATRIX)LIS_V2P(A),*matrix_type);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_type_f"
void lis_matrix_get_type_f(LIS_MATRIX_F *A, LIS_INT *matrix_type, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_get_type((LIS_MATRIX)LIS_V2P(A),matrix_type);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_f"
void lis_matrix_malloc_f(LIS_MATRIX_F *A, LIS_INT *nnz_row, LIS_INT *nnz, LIS_INT *ierr)
{
	LIS_MATRIX AA;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	*ierr = lis_matrix_malloc(AA,*nnz_row,nnz);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_csr_f"
void lis_matrix_malloc_csr_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_csr(*n, *nnz, &ptr, &index, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_csc_f"
void lis_matrix_malloc_csc_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_csc(*n, *nnz, &ptr, &index, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_msr_f"
void lis_matrix_malloc_msr_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *ndz, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_msr(*n, *nnz, *ndz, &index, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_dia_f"
void lis_matrix_malloc_dia_f(LIS_INT *n, LIS_INT *nnd, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_dia(*n, *nnd, &index, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_ell_f"
void lis_matrix_malloc_ell_f(LIS_INT *n, LIS_INT *maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_ell(*n, *maxnzr, &index, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_jad_f"
void lis_matrix_malloc_jad_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_jad(*n, *nnz, *maxnzr, &perm, &ptr, &index, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_bsr_f"
void lis_matrix_malloc_bsr_f(LIS_INT *n, LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_bsr(*n, *bnr, *bnc, *bnnz, &bptr, &bindex, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_malloc_malloc_bsc_f"
void lis_matrix_malloc_bsc_f(LIS_INT *n, LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_bsc(*n, *bnr, *bnc, *bnnz, &bptr, &bindex, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_coo_f"
void lis_matrix_malloc_coo_f(LIS_INT *nnz, LIS_INT *row, LIS_INT *col, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_coo(*nnz, &row, &col, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_dns_f"
void lis_matrix_malloc_dns_f(LIS_INT *n, LIS_INT *np, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_dns(*n, *np, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_vbr_f"
void lis_matrix_malloc_vbr_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *nr, LIS_INT *nc, LIS_INT *bnnz, LIS_INT *row, LIS_INT *col, LIS_INT *ptr, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_vbr(*n, *nnz, *nr, *nc, *bnnz, &row, &col, &ptr, &bptr, &bindex, &value);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_csr_f"
void lis_matrix_set_csr_f(LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *Amat, LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_csr(*nnz,ptr,index,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_csc_f"
void lis_matrix_set_csc_f(LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_csc(*nnz,ptr,index,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_msr_f"
void lis_matrix_set_msr_f(LIS_INT *nnz, LIS_INT *ndz, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;	

	*ierr = lis_matrix_set_msr(*nnz,*ndz,index,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_dia_f"
void lis_matrix_set_dia_f(LIS_INT *nnd, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_dia(*nnd,index,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_ell_f"
void lis_matrix_set_ell_f(LIS_INT *maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_ell(*maxnzr,index,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_jad_f"
void lis_matrix_set_jad_f(LIS_INT *nnz, LIS_INT *maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_jad(*nnz,*maxnzr,perm,ptr,index,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_bsr_f"
void lis_matrix_set_bsr_f(LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_bsr(*bnr,*bnc,*bnnz,bptr,bindex,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_bsc_f"
void lis_matrix_set_bsc_f(LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_bsc(*bnr,*bnc,*bnnz,bptr,bindex,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_coo_f"
void lis_matrix_set_coo_f(LIS_INT *nnz, LIS_INT *row, LIS_INT *col, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_coo(*nnz,row,col,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_dns_f"
void lis_matrix_set_dns_f(LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_dns(value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_vbr_f"
void lis_matrix_set_vbr_f(LIS_INT *nnz, LIS_INT *nr, LIS_INT *nc, LIS_INT *bnnz, LIS_INT *row, LIS_INT *col, LIS_INT *ptr, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX_F *Amat,LIS_INT *ierr)
{
	LIS_MATRIX A;
	
	LIS_DEBUG_FUNC_IN;

	A = (LIS_MATRIX)LIS_V2P(Amat);
	A->is_fallocated = 1;
	
	*ierr = lis_matrix_set_vbr(*nnz,*nr,*nc,*bnnz,row,col,ptr,bptr,bindex,value,A);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_assemble_f"
void lis_matrix_assemble_f(LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_MATRIX AA;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	*ierr = lis_matrix_assemble(AA);
	if( *ierr )	return;	
	*A    = LIS_P2V(AA);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_f"
void lis_matrix_convert_f(LIS_MATRIX_F *Ain, LIS_MATRIX_F *Aout, LIS_INT *ierr)
{
	LIS_MATRIX AAin,AAout;

	LIS_DEBUG_FUNC_IN;

	AAin    = (LIS_MATRIX)LIS_V2P(Ain);
	AAout   = (LIS_MATRIX)LIS_V2P(Aout);
	*ierr   = lis_matrix_convert(AAin,AAout);
	if( *ierr )	return;	
	*Aout   = LIS_P2V(AAout);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_f"
void lis_matrix_copy_f(LIS_MATRIX_F *Ain, LIS_MATRIX_F *Aout, LIS_INT *ierr)
{
	LIS_MATRIX AAin,AAout;

	LIS_DEBUG_FUNC_IN;

	AAin    = (LIS_MATRIX)LIS_V2P(Ain);
	AAout   = (LIS_MATRIX)LIS_V2P(Aout);
	*ierr   = lis_matrix_copy(AAin,AAout);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_axpy_f"
void lis_matrix_axpy_f(LIS_SCALAR_F *alpha, LIS_MATRIX_F *A, LIS_MATRIX_F *B, LIS_INT *ierr)
{
	LIS_MATRIX AA,BB;
	LIS_SCALAR ss;

	LIS_DEBUG_FUNC_IN;

	ss      = (LIS_SCALAR)LIS_V2P(alpha);
	AA      = (LIS_MATRIX)LIS_V2P(A);
	BB      = (LIS_MATRIX)LIS_V2P(B);
	*ierr   = lis_matrix_axpy(ss,AA,BB);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_xpay_f"
void lis_matrix_xpay_f(LIS_SCALAR_F *alpha, LIS_MATRIX_F *A, LIS_MATRIX_F *B, LIS_INT *ierr)
{
	LIS_MATRIX AA,BB;
	LIS_SCALAR ss;

	LIS_DEBUG_FUNC_IN;

	ss      = (LIS_SCALAR)LIS_V2P(alpha);
	AA      = (LIS_MATRIX)LIS_V2P(A);
	BB      = (LIS_MATRIX)LIS_V2P(B);
	*ierr   = lis_matrix_xpay(ss,AA,BB);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_axpyz_f"
void lis_matrix_axpyz_f(LIS_SCALAR_F *alpha, LIS_MATRIX_F *A, LIS_MATRIX_F *B, LIS_MATRIX_F *C, LIS_INT *ierr)
{
	LIS_MATRIX AA,BB,CC;
	LIS_SCALAR ss;

	LIS_DEBUG_FUNC_IN;

	ss      = (LIS_SCALAR)LIS_V2P(alpha);
	AA      = (LIS_MATRIX)LIS_V2P(A);
	BB      = (LIS_MATRIX)LIS_V2P(B);
	CC      = (LIS_MATRIX)LIS_V2P(C);	
	*ierr   = lis_matrix_axpyz(ss,AA,BB,CC);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scale_f"
void lis_matrix_scale_f(LIS_MATRIX_F *A, LIS_VECTOR_F *b, LIS_VECTOR_F *d, LIS_INT *action, LIS_INT *ierr)
{
	LIS_MATRIX AA;
	LIS_VECTOR dd,bb;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	bb    = (LIS_VECTOR)LIS_V2P(b);
	dd    = (LIS_VECTOR)LIS_V2P(d);
	*ierr = lis_matrix_scale(AA,bb,dd,*action);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_matrix_psd_reset_scale_f"
void lis_matrix_psd_reset_scale_f(LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_MATRIX AA;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	*ierr = lis_matrix_psd_reset_scale(AA);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_f"
void lis_matrix_get_diagonal_f(LIS_MATRIX_F *A, LIS_VECTOR_F *d, LIS_INT *ierr)
{
	LIS_MATRIX AA;
	LIS_VECTOR dd;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	dd    = (LIS_VECTOR)LIS_V2P(d);
	*ierr = lis_matrix_get_diagonal(AA,dd);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_f"
void lis_matrix_shift_diagonal_f(LIS_MATRIX_F *A, LIS_SCALAR_F *sigma, LIS_INT *ierr)
{
	LIS_MATRIX AA;
	LIS_SCALAR ss;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	ss    = (LIS_SCALAR)LIS_V2P(sigma);
	*ierr = lis_matrix_shift_diagonal(AA,ss);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_matrix_f"
void lis_matrix_shift_matrix_f(LIS_MATRIX_F *A, LIS_MATRIX_F *B, LIS_SCALAR_F *sigma, LIS_INT *ierr)
{
	LIS_MATRIX AA;
	LIS_MATRIX BB;	
	LIS_SCALAR ss;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	BB    = (LIS_MATRIX)LIS_V2P(B);	
	ss    = (LIS_SCALAR)LIS_V2P(sigma);
	*ierr = lis_matrix_shift_matrix(AA,BB,ss);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_blocksize_f"
void lis_matrix_set_blocksize_f(LIS_MATRIX_F *A, LIS_INT *bnr, LIS_INT *bnc, LIS_INT *row, LIS_INT *col, LIS_INT *ierr)
{
	LIS_MATRIX AA;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	if( *row==0 )
	{
		*ierr = lis_matrix_set_blocksize(AA,*bnr,*bnc,NULL,NULL);
	}
	else
	{
		*ierr = lis_matrix_set_blocksize(AA,*bnr,*bnc,row,col);
	}
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_unset_f"
void lis_matrix_unset_f(LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_unset(((LIS_MATRIX)LIS_V2P(A)));
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}


#endif


