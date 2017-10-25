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

LIS_MATVEC_XXX lis_matvec_bsc_xxx[4][4] = { 
	{lis_matvec_bsc_1x1, lis_matvec_bsc_1x2, lis_matvec_bsc_1x3, lis_matvec_bsc_1x4},
	{lis_matvec_bsc_2x1, lis_matvec_bsc_2x2, lis_matvec_bsc_2x3, lis_matvec_bsc_2x4},
	{lis_matvec_bsc_3x1, lis_matvec_bsc_3x2, lis_matvec_bsc_3x3, lis_matvec_bsc_3x4},
	{lis_matvec_bsc_4x1, lis_matvec_bsc_4x2, lis_matvec_bsc_4x3, lis_matvec_bsc_4x4}
};

void lis_matvec_bsc(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k;
	LIS_INT bi,bj,bc,bs;
	LIS_INT nr,nc,bnr,bnc;
	LIS_INT n;
	LIS_SCALAR t;

	n   = A->n;
	nr  = A->nr;
	nc  = A->nc;
	bnr = A->bnr;
	bnc = A->bnc;
	bs  = bnr*bnc;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,i,j,t,k)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k = 0;
			for(i=0;i<bnr;i++)
			{
				t = 0.0;
				for(j=0;j<bnc;j++)
				{
					t += A->D->value[bi*bs+ j*bnr + i] * x[bi*bnr+j];
					k++;
				}
				y[bi*bnr+i] = t;
			}
		}
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,k,bi,bj,bc)
		#endif
		for(bi=0;bi<nc;bi++)
		{
			for(bc=A->L->bptr[bi];bc<A->L->bptr[bi+1];bc++)
			{
				bj   = A->L->bindex[bc] * bnr;
				k    = bc*bs;
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						y[bj+i] += A->L->value[k] * x[bi*bnc+j];
						k++;
					}
				}
			}
			for(bc=A->U->bptr[bi];bc<A->U->bptr[bi+1];bc++)
			{
				bj   = A->U->bindex[bc] * bnr;
				k    = bc*bs;
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						y[bj+i] += A->U->value[k] * x[bi*bnc+j];
						k++;
					}
				}
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0; i<n; i++)
		{
			y[i] = 0.0;
		}
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,k,bi,bj,bc)
		#endif
		for(bi=0;bi<nc;bi++)
		{
			for(bc=A->bptr[bi];bc<A->bptr[bi+1];bc++)
			{
				bj   = A->bindex[bc] * bnr;
				k    = bc*bs;
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						y[bj+i] += A->value[k] * x[bi*bnc+j];
						k++;
					}
				}
			}
		}
	}
}

void lis_matvec_bsc_1x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0;

	nr   = A->nr;
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 += A->value[j] * x[jj];
		}
		y[i] = t0;
	}
}

void lis_matvec_bsc_1x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 += A->value[j*2+0] * x[jj*2+0];
			t0 += A->value[j*2+1] * x[jj*2+1];
		}
		y[i] = t0;
	}
}

void lis_matvec_bsc_1x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 += A->value[j*3+0] * x[jj*3+0];
			t0 += A->value[j*3+1] * x[jj*3+1];
			t0 += A->value[j*3+2] * x[jj*3+2];
		}
		y[i] = t0;
	}
}

void lis_matvec_bsc_1x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 += A->value[j*4+0] * x[jj*4+0];
			t0 += A->value[j*4+1] * x[jj*4+1];
			t0 += A->value[j*4+2] * x[jj*4+2];
			t0 += A->value[j*4+3] * x[jj*4+3];
		}
		y[i] = t0;
	}
}

void lis_matvec_bsc_2x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1;

	nr   = A->nr;

	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,jj,js,je,t0,t1)
		#endif
		for(i=0; i<nr; i++)
		{
			t0 = A->D->value[4*i+0] * x[2*i+0] + A->D->value[4*i+2] * x[2*i+1];
			t1 = A->D->value[4*i+1] * x[2*i+0] + A->D->value[4*i+3] * x[2*i+1];
			js = A->L->bptr[i];
			je = A->L->bptr[i+1];
			for(j=js;j<je;j++)
			{
				jj  = A->L->bindex[j];
				t0 += A->L->value[j*4+0] * x[jj*2+0];
				t1 += A->L->value[j*4+1] * x[jj*2+0];
				t0 += A->L->value[j*4+2] * x[jj*2+1];
				t1 += A->L->value[j*4+3] * x[jj*2+1];
			}
			js = A->U->bptr[i];
			je = A->U->bptr[i+1];
			for(j=js;j<je;j++)
			{
				jj  = A->U->bindex[j];
				t0 += A->U->value[j*4+0] * x[jj*2+0];
				t1 += A->U->value[j*4+1] * x[jj*2+0];
				t0 += A->U->value[j*4+2] * x[jj*2+1];
				t1 += A->U->value[j*4+3] * x[jj*2+1];
			}
			y[2*i+0] = t0;
			y[2*i+1] = t1;
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,jj,js,je,t0,t1)
		#endif
		for(i=0; i<nr; i++)
		{
			js = A->bptr[i];
			je = A->bptr[i+1];
			t0 = 0.0;
			t1 = 0.0;
			for(j=js;j<je;j++)
			{
				jj  = A->bindex[j];
				t0 += A->value[j*4+0] * x[jj*2+0];
				t1 += A->value[j*4+1] * x[jj*2+0];
				t0 += A->value[j*4+2] * x[jj*2+1];
				t1 += A->value[j*4+3] * x[jj*2+1];
			}
			y[2*i+0] = t0;
			y[2*i+1] = t1;
		}
	}
}

void lis_matvec_bsc_2x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 += A->value[j*2+0] * x[jj];
			t1 += A->value[j*2+1] * x[jj];
		}
		y[2*i+0] = t0;
		y[2*i+1] = t1;
	}
}

void lis_matvec_bsc_2x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 += A->value[j*6+0] * x[jj*3+0];
			t1 += A->value[j*6+1] * x[jj*3+0];
			t0 += A->value[j*6+2] * x[jj*3+1];
			t1 += A->value[j*6+3] * x[jj*3+1];
			t0 += A->value[j*6+4] * x[jj*3+2];
			t1 += A->value[j*6+5] * x[jj*3+2];
		}
		y[2*i+0] = t0;
		y[2*i+1] = t1;
	}
}

void lis_matvec_bsc_2x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*8+0] * x[jj*4+0];
			t1 +=  A->value[j*8+1] * x[jj*4+0];
			t0 +=  A->value[j*8+2] * x[jj*4+1];
			t1 +=  A->value[j*8+3] * x[jj*4+1];
			t0 +=  A->value[j*8+4] * x[jj*4+2];
			t1 +=  A->value[j*8+5] * x[jj*4+2];
			t0 +=  A->value[j*8+6] * x[jj*4+3];
			t1 +=  A->value[j*8+7] * x[jj*4+3];
		}
		y[2*i+0] = t0;
		y[2*i+1] = t1;
	}
}

void lis_matvec_bsc_3x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*9+0] * x[jj*3+0];
			t1 +=  A->value[j*9+1] * x[jj*3+0];
			t2 +=  A->value[j*9+2] * x[jj*3+0];
			t0 +=  A->value[j*9+3] * x[jj*3+1];
			t1 +=  A->value[j*9+4] * x[jj*3+1];
			t2 +=  A->value[j*9+5] * x[jj*3+1];
			t0 +=  A->value[j*9+6] * x[jj*3+2];
			t1 +=  A->value[j*9+7] * x[jj*3+2];
			t2 +=  A->value[j*9+8] * x[jj*3+2];
		}
		y[3*i+0] = t0;
		y[3*i+1] = t1;
		y[3*i+2] = t2;
	}
}

void lis_matvec_bsc_3x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj,ii;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2,ii)
	#endif
	for(i=0; i<nr; i++)
	{
		ii = 3*i;
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*3+0] * x[jj];
			t1 +=  A->value[j*3+1] * x[jj];
			t2 +=  A->value[j*3+2] * x[jj];
		}
		y[ii+0] = t0;
		y[ii+1] = t1;
		y[ii+2] = t2;
	}
}

void lis_matvec_bsc_3x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*6+0] * x[jj*2+0];
			t1 +=  A->value[j*6+1] * x[jj*2+0];
			t2 +=  A->value[j*6+2] * x[jj*2+0];
			t0 +=  A->value[j*6+3] * x[jj*2+1];
			t1 +=  A->value[j*6+4] * x[jj*2+1];
			t2 +=  A->value[j*6+5] * x[jj*2+1];
		}
		y[3*i+0] = t0;
		y[3*i+1] = t1;
		y[3*i+2] = t2;
	}
}

void lis_matvec_bsc_3x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*12+ 0] * x[jj*4+0];
			t1 +=  A->value[j*12+ 1] * x[jj*4+0];
			t2 +=  A->value[j*12+ 2] * x[jj*4+0];
			t0 +=  A->value[j*12+ 3] * x[jj*4+1];
			t1 +=  A->value[j*12+ 4] * x[jj*4+1];
			t2 +=  A->value[j*12+ 5] * x[jj*4+1];
			t0 +=  A->value[j*12+ 6] * x[jj*4+2];
			t1 +=  A->value[j*12+ 7] * x[jj*4+2];
			t2 +=  A->value[j*12+ 8] * x[jj*4+2];
			t0 +=  A->value[j*12+ 9] * x[jj*4+3];
			t1 +=  A->value[j*12+10] * x[jj*4+3];
			t2 +=  A->value[j*12+11] * x[jj*4+3];
		}
		y[3*i+0] = t0;
		y[3*i+1] = t1;
		y[3*i+2] = t2;
	}
}

void lis_matvec_bsc_4x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2,t3;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2,t3)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		t3 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*16+ 0] * x[jj*4+0];
			t1 +=  A->value[j*16+ 1] * x[jj*4+0];
			t2 +=  A->value[j*16+ 2] * x[jj*4+0];
			t3 +=  A->value[j*16+ 3] * x[jj*4+0];
			t0 +=  A->value[j*16+ 4] * x[jj*4+1];
			t1 +=  A->value[j*16+ 5] * x[jj*4+1];
			t2 +=  A->value[j*16+ 6] * x[jj*4+1];
			t3 +=  A->value[j*16+ 7] * x[jj*4+1];
			t0 +=  A->value[j*16+ 8] * x[jj*4+2];
			t1 +=  A->value[j*16+ 9] * x[jj*4+2];
			t2 +=  A->value[j*16+10] * x[jj*4+2];
			t3 +=  A->value[j*16+11] * x[jj*4+2];
			t0 +=  A->value[j*16+12] * x[jj*4+3];
			t1 +=  A->value[j*16+13] * x[jj*4+3];
			t2 +=  A->value[j*16+14] * x[jj*4+3];
			t3 +=  A->value[j*16+15] * x[jj*4+3];
		}
		y[4*i+0] = t0;
		y[4*i+1] = t1;
		y[4*i+2] = t2;
		y[4*i+3] = t3;
	}
}

void lis_matvec_bsc_4x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2,t3;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2,t3)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		t3 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*4+ 0] * x[jj];
			t1 +=  A->value[j*4+ 1] * x[jj];
			t2 +=  A->value[j*4+ 2] * x[jj];
			t3 +=  A->value[j*4+ 3] * x[jj];
		}
		y[4*i+0] = t0;
		y[4*i+1] = t1;
		y[4*i+2] = t2;
		y[4*i+3] = t3;
	}
}

void lis_matvec_bsc_4x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2,t3;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2,t3)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		t3 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*8+ 0] * x[jj*2+0];
			t1 +=  A->value[j*8+ 1] * x[jj*2+0];
			t2 +=  A->value[j*8+ 2] * x[jj*2+0];
			t3 +=  A->value[j*8+ 3] * x[jj*2+0];
			t0 +=  A->value[j*8+ 4] * x[jj*2+1];
			t1 +=  A->value[j*8+ 5] * x[jj*2+1];
			t2 +=  A->value[j*8+ 6] * x[jj*2+1];
			t3 +=  A->value[j*8+ 7] * x[jj*2+1];
		}
		y[4*i+0] = t0;
		y[4*i+1] = t1;
		y[4*i+2] = t2;
		y[4*i+3] = t3;
	}
}

void lis_matvec_bsc_4x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,js,je,jj;
	LIS_INT nr;
	LIS_SCALAR t0,t1,t2,t3;

	nr   = A->nr;

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,jj,js,je,t0,t1,t2,t3)
	#endif
	for(i=0; i<nr; i++)
	{
		js = A->bptr[i];
		je = A->bptr[i+1];
		t0 = 0.0;
		t1 = 0.0;
		t2 = 0.0;
		t3 = 0.0;
		for(j=js;j<je;j++)
		{
			jj  = A->bindex[j];
			t0 +=  A->value[j*12+ 0] * x[jj*3+0];
			t1 +=  A->value[j*12+ 1] * x[jj*3+0];
			t2 +=  A->value[j*12+ 2] * x[jj*3+0];
			t3 +=  A->value[j*12+ 3] * x[jj*3+0];
			t0 +=  A->value[j*12+ 4] * x[jj*3+1];
			t1 +=  A->value[j*12+ 5] * x[jj*3+1];
			t2 +=  A->value[j*12+ 6] * x[jj*3+1];
			t3 +=  A->value[j*12+ 7] * x[jj*3+1];
			t0 +=  A->value[j*12+ 8] * x[jj*3+2];
			t1 +=  A->value[j*12+ 9] * x[jj*3+2];
			t2 +=  A->value[j*12+10] * x[jj*3+2];
			t3 +=  A->value[j*12+11] * x[jj*3+2];
		}
		y[4*i+0] = t0;
		y[4*i+1] = t1;
		y[4*i+2] = t2;
		y[4*i+3] = t3;
	}
}

void lis_matvech_bsc(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[])
{
	LIS_INT i,j,k;
	LIS_INT bi,bj,bc,bs;
	LIS_INT nr,nc,bnr,bnc;
	LIS_INT n,np;
	#ifdef _OPENMP
		LIS_INT nprocs;
		LIS_SCALAR t;
		LIS_SCALAR *w;
	#endif

	n   = A->n;
	np  = A->np;
	nr  = A->nr;
	nc  = A->nc;
	bnr = A->bnr;
	bnc = A->bnc;
	bs  = bnr*bnc;
	if( A->is_splited )
	{
		for(i=0;i<n;i++)
		{
			y[i] = 0.0;
		}
		for(bi=0;bi<nr;bi++)
		{
			k    = bi*bs;
			for(j=0;j<bnc;j++)
			{
				for(i=0;i<bnr;i++)
				{
					y[bi*bnr+j] += conj(A->D->value[k++]) * x[bi*bnr+i];
				}
			}
		}
		for(bi=0;bi<nc;bi++)
		{
			for(bc=A->L->bptr[bi];bc<A->L->bptr[bi+1];bc++)
			{
				bj   = A->L->bindex[bc] * bnr;
				k    = bc*bs;
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						y[bj+j] += conj(A->L->value[k]) * x[bi*bnr+i];
						k++;
					}
				}
			}
			for(bc=A->U->bptr[bi];bc<A->U->bptr[bi+1];bc++)
			{
				bj   = A->U->bindex[bc] * bnr;
				k    = bc*bs;
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						y[bj+j] += conj(A->U->value[k]) * x[bi*bnr+i];
						k++;
					}
				}
			}
		}
	}
	else
	{
		#ifdef _OPENMP
			nprocs = omp_get_max_threads();
			w = (LIS_SCALAR *)lis_malloc( nprocs*np*sizeof(LIS_SCALAR),"lis_matvech_bsc::w" );

			#pragma omp parallel private(bi,bc,bj,i,j,k)
			{
				#pragma omp for
				for(j=0;j<nprocs;j++)
				{
					memset( &w[j*np], 0, np*sizeof(LIS_SCALAR) );
				}
				#pragma omp for 
				for(bi=0;bi<nc;bi++)
				{
					for(bc=A->bptr[bi];bc<A->bptr[bi+1];bc++)
					{
						bj   = A->bindex[bc] * bnr;
						k    = bc*bs;
						for(j=0;j<bnc;j++)
						{
							for(i=0;i<bnr;i++)
							{
								w[bi*bnc+j] += conj(A->value[k]) * x[bj+i];
								k++;
							}
						}
					}
				}
				#pragma omp barrier
				#pragma omp for 
				for(i=0;i<np;i++)
				{
					t = 0.0;
					for(j=0;j<nprocs;j++)
					{
						t += w[j*np+i];
					}
					y[i] = t;
				}
			}
			lis_free(w);
		#else
			for(i=0; i<n; i++)
			{
				y[i] = 0.0;
			}
			for(bi=0;bi<nc;bi++)
			{
				for(bc=A->bptr[bi];bc<A->bptr[bi+1];bc++)
				{
					bj   = A->bindex[bc] * bnr;
					k    = bc*bs;
					for(j=0;j<bnc;j++)
					{
						for(i=0;i<bnr;i++)
						{
							y[bi*bnc+j] += conj(A->value[k]) * x[bj+i];
							k++;
						}
					}
				}
			}
		#endif
	}
}
