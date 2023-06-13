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

#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/**************************************************************
 * lis_vector_dot		v   <- x^H * y (Hermitian)
 * lis_vector_nhdot		v   <- x^T * y (non-Hermitian)
 * lis_vector_nrm2		v   <- ||x||_2
 * lis_vector_nrm1		v   <- ||x||_1
 * lis_vector_nrmi		v   <- ||x||_infinity
 * lis_vector_sum		v   <- sum x_i
 **************************************************************/

/**********************/
/* v <- x^H * y       */
/**********************/
#undef __FUNC__
#define __FUNC__ "lis_vector_dot"
LIS_INT lis_vector_dot(LIS_VECTOR vx, LIS_VECTOR vy, LIS_SCALAR *value)
{
	LIS_INT i,n;
	LIS_SCALAR dot;
	LIS_SCALAR *x,*y;
	LIS_SCALAR tmp;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	n = vx->n;
	#ifndef NO_ERROR_CHECK
		if( n!=vy->n )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"length of vector x and y is not equal\n");
			return LIS_ERR_ILL_ARG;
		}
	#endif

	x      = vx->value;
	y      = vy->value;
	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,tmp,my_rank)
		{
			my_rank = omp_get_thread_num();
			tmp     = 0.0;
			#ifdef USE_VEC_COMP
		    #pragma cdir nodep
			#endif
			#pragma omp for
			for(i=0; i<n; i++)
			{
				tmp += conj(x[i])*y[i];
			}
			lis_vec_tmp[my_rank*LIS_VEC_TMP_PADD] = tmp;
		}
		dot = 0.0;
		for(i=0;i<nprocs;i++)
		{
			dot += lis_vec_tmp[i*LIS_VEC_TMP_PADD];
		}
	#else
		dot  = 0.0;
		#ifdef USE_VEC_COMP
	    #pragma cdir nodep
		#endif
		for(i=0; i<n; i++)
		{
			dot += conj(x[i])*y[i];
		}
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(&dot,&tmp,1,LIS_MPI_SCALAR,MPI_SUM,comm);
		*value = tmp;
	#else
		*value = dot;
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/**********************/
/* v <- x^T * y       */
/**********************/
#undef __FUNC__
#define __FUNC__ "lis_vector_nhdot"
LIS_INT lis_vector_nhdot(LIS_VECTOR vx, LIS_VECTOR vy, LIS_SCALAR *value)
{
	LIS_INT i,n;
	LIS_SCALAR dot;
	LIS_SCALAR *x,*y;
	LIS_SCALAR tmp;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	n = vx->n;
	#ifndef NO_ERROR_CHECK
		if( n!=vy->n )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"length of vector x and y is not equal\n");
			return LIS_ERR_ILL_ARG;
		}
	#endif

	x      = vx->value;
	y      = vy->value;
	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,tmp,my_rank)
		{
			my_rank = omp_get_thread_num();
			tmp     = 0.0;
			#ifdef USE_VEC_COMP
		    #pragma cdir nodep
			#endif
			#pragma omp for
			for(i=0; i<n; i++)
			{
				tmp += x[i]*y[i];
			}
			lis_vec_tmp[my_rank*LIS_VEC_TMP_PADD] = tmp;
		}
		dot = 0.0;
		for(i=0;i<nprocs;i++)
		{
			dot += lis_vec_tmp[i*LIS_VEC_TMP_PADD];
		}
	#else
		dot  = 0.0;
		#ifdef USE_VEC_COMP
	    #pragma cdir nodep
		#endif
		for(i=0; i<n; i++)
		{
			dot += x[i]*y[i];
		}
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(&dot,&tmp,1,LIS_MPI_SCALAR,MPI_SUM,comm);
		*value = tmp;
	#else
		*value = dot;
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/**********************/
/* v <- ||x||_2       */
/**********************/
#undef __FUNC__
#define __FUNC__ "lis_vector_nrm2"
LIS_INT lis_vector_nrm2(LIS_VECTOR vx, LIS_REAL *value)
{
	LIS_INT i,n;
	LIS_SCALAR dot;
	LIS_SCALAR *x;
	LIS_SCALAR tmp;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	n      = vx->n;

	x      = vx->value;
	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,tmp,my_rank)
		{
			my_rank = omp_get_thread_num();
			tmp     = 0.0;
			#ifdef USE_VEC_COMP
		    #pragma cdir nodep
			#endif
			#pragma omp for
			for(i=0; i<n; i++)
			{
				tmp += conj(x[i])*x[i];
			}
			lis_vec_tmp[my_rank*LIS_VEC_TMP_PADD] = tmp;
		}
		dot = 0.0;
		for(i=0;i<nprocs;i++)
		{
			dot += lis_vec_tmp[i*LIS_VEC_TMP_PADD];
		}
	#else
		dot  = 0.0;
		#ifdef USE_VEC_COMP
	    #pragma cdir nodep
		#endif
		for(i=0; i<n; i++)
		{
			dot += conj(x[i])*x[i];
		}
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(&dot,&tmp,1,LIS_MPI_REAL,MPI_SUM,comm);
		*value = sqrt(tmp);
	#else
		*value = sqrt(dot);
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/**********************/
/* v <- ||x||_1       */
/**********************/
#undef __FUNC__
#define __FUNC__ "lis_vector_nrm1"
LIS_INT lis_vector_nrm1(LIS_VECTOR vx, LIS_REAL *value)
{
	LIS_INT i,n;
	LIS_SCALAR sum;
	LIS_SCALAR *x;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
		LIS_SCALAR tmp;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	n   = vx->n;

	x      = vx->value;
	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
	sum    = 0.0;
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,tmp,my_rank)
		{
			my_rank = omp_get_thread_num();
			tmp     = 0.0;
			#ifdef USE_VEC_COMP
		    #pragma cdir nodep
			#endif
			#pragma omp for
			for(i=0; i<n; i++)
			{
				tmp += fabs(x[i]);
			}
			lis_vec_tmp[my_rank*LIS_VEC_TMP_PADD] = tmp;
		}
		for(i=0;i<nprocs;i++)
		{
			sum += lis_vec_tmp[i*LIS_VEC_TMP_PADD];
		}
	#else
		#ifdef USE_VEC_COMP
	    #pragma cdir nodep
		#endif
		for(i=0; i<n; i++)
		{
			sum += fabs(x[i]);
		}
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(&sum,value,1,LIS_MPI_REAL,MPI_SUM,comm);
	#else
		*value = sum;
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/***********************/
/* v <- ||x||_infinity */
/***********************/
#undef __FUNC__
#define __FUNC__ "lis_vector_nrmi"
LIS_INT lis_vector_nrmi(LIS_VECTOR vx, LIS_REAL *value)
{
	LIS_INT i,n;
	LIS_SCALAR *x;
	LIS_REAL tmp;

	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	n   = vx->n;

	x      = vx->value;
	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
		tmp     = 0.0;		
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,tmp,my_rank)
		{
			my_rank = omp_get_thread_num();
			#ifdef USE_VEC_COMP
		    #pragma cdir nodep
			#endif
			#pragma omp for
			for(i=0; i<n; i++)
			{
			  if (fabs(x[i]) > tmp) 
			    {
				  tmp = fabs(x[i]);
			    }
			}
			lis_vec_tmp[my_rank*LIS_VEC_TMP_PADD] = tmp;
		}
		for(i=0;i<nprocs;i++)
		{
		  if ((LIS_REAL)lis_vec_tmp[i*LIS_VEC_TMP_PADD] > tmp) 
		    {
		      tmp = lis_vec_tmp[i*LIS_VEC_TMP_PADD];
		    }
		}
	#else
		#ifdef USE_VEC_COMP
	    #pragma cdir nodep
		#endif
		for(i=0; i<n; i++)
		{
		  if (fabs(x[i]) > tmp)
		    {
		      tmp = fabs(x[i]);
		    }
		}
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(&tmp,value,1,LIS_MPI_REAL,MPI_MAX,comm);
	#else
		*value = tmp;
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/**********************/
/* v <- sum x_i       */
/**********************/
#undef __FUNC__
#define __FUNC__ "lis_vector_sum"
LIS_INT lis_vector_sum(LIS_VECTOR vx, LIS_SCALAR *value)
{
	LIS_INT i,n;
	LIS_SCALAR sum;
	LIS_SCALAR *x;
	#ifdef _OPENMP
		LIS_INT nprocs,my_rank;
		LIS_SCALAR tmp;
	#endif
	#ifdef USE_MPI
		MPI_Comm comm;
	#endif

	LIS_DEBUG_FUNC_IN;

	n   = vx->n;

	x      = vx->value;
	#ifdef USE_MPI
		comm   = vx->comm;
	#endif
	sum    = 0.0;
	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
		#pragma omp parallel private(i,tmp,my_rank)
		{
			my_rank = omp_get_thread_num();
			tmp     = 0.0;
			#ifdef USE_VEC_COMP
		    #pragma cdir nodep
			#endif
			#pragma omp for
			for(i=0; i<n; i++)
			{
				tmp += x[i];
			}
			lis_vec_tmp[my_rank*LIS_VEC_TMP_PADD] = tmp;
		}
		for(i=0;i<nprocs;i++)
		{
			sum += lis_vec_tmp[i*LIS_VEC_TMP_PADD];
		}
	#else
		#ifdef USE_VEC_COMP
	    #pragma cdir nodep
		#endif
		for(i=0; i<n; i++)
		{
			sum += x[i];
		}
	#endif
	#ifdef USE_MPI
		MPI_Allreduce(&sum,value,1,LIS_MPI_SCALAR,MPI_SUM,comm);
	#else
		*value = sum;
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

