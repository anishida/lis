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

LIS_MATVEC_FUNC LIS_MATVEC  = lis_matvec;
LIS_MATVEC_FUNC LIS_MATVECH = lis_matvech;

#undef __FUNC__
#define __FUNC__ "lis_matvec"
LIS_INT lis_matvec(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	if( X->precision==LIS_PRECISION_DEFAULT )
	{
		x = X->value;
		y = Y->value;

		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_csr(A, x, y);
			break;
		case LIS_MATRIX_BSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			if( A->bnr<=4 && A->bnc<=4 )
			{
				lis_matvec_bsr_xxx[A->bnr-1][A->bnc-1](A, x, y);
			}
			else
			{
				lis_matvec_bsr(A, x, y);
			}
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_csc(A, x, y);
			break;
		case LIS_MATRIX_BSC:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_bsc(A, x, y);
			break;
		case LIS_MATRIX_MSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_msr(A, x, y);
			break;
		case LIS_MATRIX_ELL:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_ell(A, x, y);
			break;
		case LIS_MATRIX_DIA:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_dia(A, x, y);
			break;
		case LIS_MATRIX_JAD:
			#ifdef USE_MPI
				#ifndef USE_OVERLAP
					LIS_MATVEC_SENDRECV;
				#else
					LIS_MATVEC_REALLOC;
				#endif
			#endif
			lis_matvec_jad(A, x, y);
			break;
		case LIS_MATRIX_VBR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_vbr(A, x, y);
			break;
		case LIS_MATRIX_DNS:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_dns(A, x, y);
			break;
		case LIS_MATRIX_COO:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_coo(A, x, y);
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#ifdef USE_QUAD_PRECISION
	else
	{

		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				lis_send_recv_mp(A->commtable,X);
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvec_csr_mp(A, X, Y);
			#else
				lis_matvec_csr_mp2(A, X, Y);
			#endif
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				lis_send_recv_mp(A->commtable,X);
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvec_csc_mp(A, X, Y);
			#else
				lis_matvec_csc_mp2(A, X, Y);
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matvech"
LIS_INT lis_matvech(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	x = X->value;
	y = Y->value;

	if( X->precision==LIS_PRECISION_DEFAULT )
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_csr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_BSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_bsr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_csc(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_BSC:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_bsc(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_MSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_msr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_ELL:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_ell(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_JAD:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_jad(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_DIA:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_dia(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_VBR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_vbr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_DNS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_dns(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_COO:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_coo(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#ifdef USE_QUAD_PRECISION
	else
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvech_csr_mp(A, X, Y);
			#else
				lis_matvech_csr_mp2(A, X, Y);
			#endif
			#ifdef USE_MPI
				lis_reduce_mp(A->commtable,Y);
			#endif
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvech_csc_mp(A, X, Y);
			#else
				lis_matvech_csc_mp2(A, X, Y);
			#endif
			#ifdef USE_MPI
				lis_reduce_mp(A->commtable,Y);
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matvec_optimize"
LIS_INT lis_matvec_optimize(LIS_MATRIX A, LIS_INT *matrix_type_maxperf)
{
  int    	   nprocs,my_rank;
  LIS_INT	   err,iter,i,matrix_type,ss,se;
  double	   time,time2,convtime;
  LIS_REAL         val;
  double	   commtime,comptime,flops,flops_maxperf;
  LIS_MATRIX       A1;
  LIS_VECTOR       X, Y;
  char             *lis_storagename2[]   = {"CSR", "CSC", "MSR", "DIA", "ELL", "JAD", "BSR", "BSC", "VBR", "COO", "DNS"};

	LIS_DEBUG_FUNC_IN;

#ifdef USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#else
	nprocs  = 1;
	my_rank = 0;
#endif

	ss = 1;
	se = 11;

	lis_vector_duplicate(A,&X);
	lis_vector_duplicate(A,&Y);
	err = lis_vector_set_all(1.0,X);

	printf("\nmeasuring matvec performance...\n");
	iter = (int)(10000000 / A->nnz) + 1; 
	flops_maxperf = 0.0;
#ifdef _LONG__LONG
	printf("number of iterations = 1e7 / %lld + 1 = %lld\n", A->nnz, iter);
#else
	printf("number of iterations = 1e7 / %d + 1 = %d\n", A->nnz, iter);
#endif
	for (matrix_type=ss;matrix_type<se;matrix_type++)
	  {
	    if ( nprocs>1 && matrix_type==9 ) continue;
	    lis_matrix_duplicate(A,&A1);
	    lis_matrix_set_type(A1,matrix_type);
	    err = lis_matrix_convert(A,A1);
	    if( err ) CHKERR(err);
		    
	    comptime = 0.0;
	    commtime = 0.0;

	    for(i=0;i<iter;i++)
	      {
#ifdef USE_MPI
		MPI_Barrier(A1->comm);
		time = lis_wtime();
		lis_send_recv(A1->commtable,X->value);
		commtime += lis_wtime() - time;
#endif
		time2 = lis_wtime();
		lis_matvec(A1,X,Y);
		comptime += lis_wtime() - time2;
	      }
	    lis_vector_nrm2(Y,&val);

	    if( my_rank==0 )
	      {
		flops = 2.0*A->nnz*iter*1.0e-6 / comptime;
#ifdef USE_MPI
#ifdef _LONG__DOUBLE
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#endif
#else
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#endif
#endif
#else
#ifdef _LONG__DOUBLE
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#endif
#else
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#endif
#endif
#endif
	      }
	    lis_matrix_destroy(A1);

	    if ( flops > flops_maxperf )
	      {
		*matrix_type_maxperf=matrix_type;
		flops_maxperf = flops;
	      }
	  }

	printf("matrix format is set to %s\n\n", lis_storagename2[*matrix_type_maxperf-1]);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
