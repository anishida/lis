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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/*************************************************************************
 * lis_array_swap		x <-> y
 * lis_array_copy		y <- x
 * lis_array_axpy		y <- y + alpha * x
 * lis_array_xpay		y <- x + alpha * y
 * lis_array_axpyz		z <- y + alpha * x
 * lis_array_scale		x <- alpha * x
 * lis_array_pmul		z_i <- x_i * y_i
 * lis_array_pdiv		z_i <- x_i / y_i
 * lis_array_set_all		x_i <- alpha
 * lis_array_abs		x_i <- |x_i|
 * lis_array_reciprocal		x_i <- 1 / x_i
 * lis_array_conjugate		x_i <- conj(x_i)
 * lis_array_shift		x_i <- x_i - sigma
 * lis_array_dot		v <- x^H * y
 * lis_array_nhdot		v <- x^T * y
 * lis_array_nrm1		v <- ||x||_1
 * lis_array_nrm2		v <- ||x||_2
 * lis_array_nrmi		v <- ||x||_infinity
 * lis_array_sum		v <- sum x_i
 * lis_array_matvec		c <- A * b
 * lis_array_matvech		c <- A^H * b
 * lis_array_matvec_ns		c <- A * b where A is not square
 * lis_array_matmat		C <- A * B
 * lis_array_matmat_ns		C <- A * B where A and B are not square
 * lis_array_ge			A <- A^-1 with Gaussian elimination
 * lis_array_solve		x <- A^-1 b
 * lis_array_cgs		Q * R <- A with classical Gram-Schmidt
 * lis_array_mgs		Q * R <- A with modified Gram-Schmidt
 * lis_array_qr			QR algorithm
 * lis_array_hqr_limit          bounded shifted Hessenberg QR
 *************************************************************************/

#undef __FUNC__
#define __FUNC__ "lis_array_swap"
LIS_INT lis_array_swap(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y)
{
  LIS_INT i;
  LIS_SCALAR t;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      t = y[i];
      y[i] = x[i];
      x[i] = t;
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_copy"
LIS_INT lis_array_copy(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      y[i] = x[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_axpy"
LIS_INT lis_array_axpy(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x, LIS_SCALAR *y)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      y[i] = alpha * x[i] + y[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_xpay"
LIS_INT lis_array_xpay(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR alpha, LIS_SCALAR *y)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      y[i] = x[i] + alpha * y[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_axpyz"
LIS_INT lis_array_axpyz(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      z[i] = alpha * x[i] + y[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_scale"
LIS_INT lis_array_scale(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      x[i] = alpha * x[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_pmul"
LIS_INT lis_array_pmul(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      z[i] = x[i]*y[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_pdiv"
LIS_INT lis_array_pdiv(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      z[i] = x[i]/y[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_set_all"
LIS_INT lis_array_set_all(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      x[i] = alpha; 
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_abs"
LIS_INT lis_array_abs(LIS_INT n, LIS_SCALAR *x)
{
  LIS_INT i;
  
  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      x[i] = fabs(x[i]);
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_reciprocal"
LIS_INT lis_array_reciprocal(LIS_INT n, LIS_SCALAR *x)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      x[i] = 1/x[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_conjugate"
LIS_INT lis_array_conjugate(LIS_INT n, LIS_SCALAR *x)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
#ifdef _COMPLEX
      x[i] = conj(x[i]);
#endif      
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_shift"
LIS_INT lis_array_shift(LIS_INT n, LIS_SCALAR sigma, LIS_SCALAR *x)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  for(i=0;i<n;i++)
    {
      x[i] = x[i] - sigma;
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_dot"
LIS_INT lis_array_dot(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *value)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  *value = 0;
  for(i=0;i<n;i++)
    {
      *value = *value + conj(x[i]) * y[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nhdot"
LIS_INT lis_array_nhdot(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *value)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  *value = 0;
  for(i=0;i<n;i++)
    {
      *value = *value + x[i] * y[i];
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nrm1"
LIS_INT lis_array_nrm1(LIS_INT n, LIS_SCALAR *x, LIS_REAL *value)
{
	LIS_INT i;
	LIS_SCALAR t;

	LIS_DEBUG_FUNC_IN;

	t = 0.0;
	for(i=0;i<n;i++)
	{
		t += fabs(x[i]);
	}
	*value = t;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nrm2"
LIS_INT lis_array_nrm2(LIS_INT n, LIS_SCALAR *x, LIS_REAL *value)
{
	LIS_INT i;
	LIS_SCALAR t;

	LIS_DEBUG_FUNC_IN;

	t = 0.0;
	for(i=0;i<n;i++)
	{
		t += conj(x[i]) * x[i];
	}
	*value = sqrt(t);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_nrmi"
LIS_INT lis_array_nrmi(LIS_INT n, LIS_SCALAR *x, LIS_REAL *value)
{
	LIS_INT i;
	LIS_REAL t;

	LIS_DEBUG_FUNC_IN;

	t = 0.0;
	for(i=0;i<n;i++)
	{
	  if (t < fabs(x[i]))
	    {
	      t = fabs(x[i]);
	    }
	}
	*value = t;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_sum"
LIS_INT lis_array_sum(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *value)
{
	LIS_INT i;
	LIS_SCALAR t;

	LIS_DEBUG_FUNC_IN;

	t = 0.0;
	for(i=0;i<n;i++)
	{
		t += x[i];
	}
	*value = t;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matvec"
LIS_INT lis_array_matvec(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT op)
{
	LIS_INT			i,j;
	LIS_SCALAR	t;

	LIS_DEBUG_FUNC_IN;

	/* y = A*x */

	if( op==LIS_INS_VALUE )
	{
		switch(n)
		{
		case 1:
			y[0] = a[0]*x[0];
			break;
		case 2:
			y[0] = a[0]*x[0] + a[2]*x[1];
			y[1] = a[1]*x[0] + a[3]*x[1];
			break;
		case 3:
			y[0] = a[0]*x[0] + a[3]*x[1] + a[6]*x[2];
			y[1] = a[1]*x[0] + a[4]*x[1] + a[7]*x[2];
			y[2] = a[2]*x[0] + a[5]*x[1] + a[8]*x[2];
			break;
		default:
			for(i=0;i<n;i++)
			{
				t = 0.0;
				for(j=0;j<n;j++)
				{
					t += a[i+j*n] * x[j];
				}
				y[i] = t;
			}
			break;
		}
	}
	else if( op==LIS_SUB_VALUE )
	{
		switch(n)
		{
		case 1:
			y[0] -= a[0]*x[0];
			break;
		case 2:
			y[0] -= a[0]*x[0] + a[2]*x[1];
			y[1] -= a[1]*x[0] + a[3]*x[1];
			break;
		case 3:
			y[0] -= a[0]*x[0] + a[3]*x[1] + a[6]*x[2];
			y[1] -= a[1]*x[0] + a[4]*x[1] + a[7]*x[2];
			y[2] -= a[2]*x[0] + a[5]*x[1] + a[8]*x[2];
			break;
		default:
			for(i=0;i<n;i++)
			{
				t = 0.0;
				for(j=0;j<n;j++)
				{
					t += a[i+j*n] * x[j];
				}
				y[i] -= t;
			}
			break;
		}
	}
	else
	{
		switch(n)
		{
		case 1:
			y[0] += a[0]*x[0];
			break;
		case 2:
			y[0] += a[0]*x[0] + a[2]*x[1];
			y[1] += a[1]*x[0] + a[3]*x[1];
			break;
		case 3:
			y[0] += a[0]*x[0] + a[3]*x[1] + a[6]*x[2];
			y[1] += a[1]*x[0] + a[4]*x[1] + a[7]*x[2];
			y[2] += a[2]*x[0] + a[5]*x[1] + a[8]*x[2];
			break;
		default:
			for(i=0;i<n;i++)
			{
				t = 0.0;
				for(j=0;j<n;j++)
				{
					t += a[i+j*n] * x[j];
				}
				y[i] += t;
			}
			break;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matvech"
LIS_INT lis_array_matvech(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT op)
{
	LIS_INT			i,j;
	LIS_SCALAR	t;

	LIS_DEBUG_FUNC_IN;

	/* y = A*x */

	if( op==LIS_INS_VALUE )
	{
		switch(n)
		{
		case 1:
			y[0] = conj(a[0])*x[0];
			break;
		case 2:
			y[0] = conj(a[0])*x[0] + conj(a[1])*x[1];
			y[1] = conj(a[2])*x[0] + conj(a[3])*x[1];
			break;
		case 3:
			y[0] = conj(a[0])*x[0] + conj(a[1])*x[1] + conj(a[2])*x[2];
			y[1] = conj(a[3])*x[0] + conj(a[4])*x[1] + conj(a[5])*x[2];
			y[2] = conj(a[6])*x[0] + conj(a[7])*x[1] + conj(a[8])*x[2];
			break;
		default:
			for(i=0;i<n;i++)
			{
				t = 0.0;
				for(j=0;j<n;j++)
				{
					t += conj(a[i*n+j]) * x[j];
				}
				y[i] = t;
			}
			break;
		}
	}
	else if( op==LIS_SUB_VALUE )
	{
		switch(n)
		{
		case 1:
			y[0] -= conj(a[0])*x[0];
			break;
		case 2:
			y[0] -= conj(a[0])*x[0] + conj(a[1])*x[1];
			y[1] -= conj(a[2])*x[0] + conj(a[3])*x[1];
			break;
		case 3:
			y[0] -= conj(a[0])*x[0] + conj(a[1])*x[1] + conj(a[2])*x[2];
			y[1] -= conj(a[3])*x[0] + conj(a[4])*x[1] + conj(a[5])*x[2];
			y[2] -= conj(a[6])*x[0] + conj(a[7])*x[1] + conj(a[8])*x[2];
			break;
		default:
			for(i=0;i<n;i++)
			{
				t = 0.0;
				for(j=0;j<n;j++)
				{
					t += conj(a[i*n+j]) * x[j];
				}
				y[i] -= t;
			}
			break;
		}
	}
	else
	{
		switch(n)
		{
		case 1:
			y[0] += conj(a[0])*x[0];
			break;
		case 2:
			y[0] += conj(a[0])*x[0] + conj(a[1])*x[1];
			y[1] += conj(a[2])*x[0] + conj(a[3])*x[1];
			break;
		case 3:
			y[0] += conj(a[0])*x[0] + conj(a[1])*x[1] + conj(a[2])*x[2];
			y[1] += conj(a[3])*x[0] + conj(a[4])*x[1] + conj(a[5])*x[2];
			y[2] += conj(a[6])*x[0] + conj(a[7])*x[1] + conj(a[8])*x[2];
			break;
		default:
			for(i=0;i<n;i++)
			{
				t = 0.0;
				for(j=0;j<n;j++)
				{
					t += conj(a[i*n+j]) * x[j];
				}
				y[i] += t;
			}
			break;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matvec_ns"
LIS_INT lis_array_matvec_ns(LIS_INT m, LIS_INT n, LIS_SCALAR *a, LIS_INT lda, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT op)
{
	LIS_INT			i,j;
	LIS_SCALAR	t;

	LIS_DEBUG_FUNC_IN;

	/* y = A*x */

	if( op==LIS_INS_VALUE )
	{
		for(i=0;i<m;i++)
		{
			t = 0.0;
			for(j=0;j<n;j++)
			{
				t += a[i+j*lda] * x[j];
			}
			y[i] = t;
		}
	}
	else if( op==LIS_SUB_VALUE )
	{
		for(i=0;i<m;i++)
		{
			t = 0.0;
			for(j=0;j<n;j++)
			{
				t += a[i+j*lda] * x[j];
			}
			y[i] -= t;
		}
	}
	else if( op==LIS_ADD_VALUE )
	{
		for(i=0;i<m;i++)
		{
			t = 0.0;
			for(j=0;j<n;j++)
			{
				t += a[i+j*lda] * x[j];
			}
			y[i] += t;
		}
	}
	else
	{
		switch(n)
		{
		case 1:
			y[0] += a[0]*x[0];
			break;
		case 2:
			y[0] += a[0]*x[0] + a[2]*x[1];
			y[1] += a[1]*x[0] + a[3]*x[1];
			break;
		case 3:
			y[0] += a[0]*x[0] + a[3]*x[1] + a[6]*x[2];
			y[1] += a[1]*x[0] + a[4]*x[1] + a[7]*x[2];
			y[2] += a[2]*x[0] + a[5]*x[1] + a[8]*x[2];
			break;
		default:
			for(i=0;i<n;i++)
			{
				t = 0.0;
				for(j=0;j<n;j++)
				{
					t += a[i+j*n] * x[j];
				}
				y[i] += t;
			}
			break;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matmat"
LIS_INT lis_array_matmat(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, LIS_INT op)
{
	LIS_INT			i,j,l;

	LIS_DEBUG_FUNC_IN;

	/* C = A*B */

	if( op==LIS_INS_VALUE )
	{
		switch(n)
		{
		case 1:
			c[0] = a[0]*b[0];
			break;
		case 2:
			c[0] = a[0]*b[0] + a[2]*b[1];
			c[1] = a[1]*b[0] + a[3]*b[1];
			c[2] = a[0]*b[2] + a[2]*b[3];
			c[3] = a[1]*b[2] + a[3]*b[3];
			break;
		case 3:
			c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
			c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
			c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
			c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
			c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
			c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
			c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
			c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
			c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
			break;
		default:
			for(j=0;j<n;j++)
			  {
			    for(i=0;i<n;i++)
			      {
				c[i+j*n] = 0.0;
			      }
			    for(l=0;l<n;l++)
			      {
				for(i=0;i<n;i++)
				  {
				    c[i+j*n] += a[i+l*n] * b[l+j*n];
				  }
			      }
			  }
			break;
		}
	}
	else if( op==LIS_SUB_VALUE )
	{
		switch(n)
		{
		case 1:
			c[0] -= a[0]*b[0];
			break;
		case 2:
			c[0] -= a[0]*b[0] + a[2]*b[1];
			c[1] -= a[1]*b[0] + a[3]*b[1];
			c[2] -= a[0]*b[2] + a[2]*b[3];
			c[3] -= a[1]*b[2] + a[3]*b[3];
			break;
		case 3:
			c[0] -= a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
			c[1] -= a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
			c[2] -= a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
			c[3] -= a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
			c[4] -= a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
			c[5] -= a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
			c[6] -= a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
			c[7] -= a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
			c[8] -= a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
			break;
		default:
			for(j=0;j<n;j++)
			{
				for(l=0;l<n;l++)
				{
					for(i=0;i<n;i++)
					{
						c[i+j*n] -= a[i+l*n] * b[l+j*n];
					}
				}
			}
			break;
		}
	}
	else
	{
		switch(n)
		{
		case 1:
			c[0] += a[0]*b[0];
			break;
		case 2:
			c[0] += a[0]*b[0] + a[2]*b[1];
			c[1] += a[1]*b[0] + a[3]*b[1];
			c[2] += a[0]*b[2] + a[2]*b[3];
			c[3] += a[1]*b[2] + a[3]*b[3];
			break;
		case 3:
			c[0] += a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
			c[1] += a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
			c[2] += a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
			c[3] += a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
			c[4] += a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
			c[5] += a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
			c[6] += a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
			c[7] += a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
			c[8] += a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
			break;
		default:
			for(j=0;j<n;j++)
			{
				for(l=0;l<n;l++)
				{
					for(i=0;i<n;i++)
			        	{
						c[i+j*n] += a[i+l*n] * b[l+j*n];
					}
				}
			}
			break;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_matmat_ns"
LIS_INT lis_array_matmat_ns(LIS_INT l, LIS_INT m, LIS_INT n, LIS_SCALAR *a, LIS_INT lda, LIS_SCALAR *b, LIS_INT ldb, LIS_SCALAR *c, LIS_INT ldc, LIS_INT op)
{
	LIS_INT			i,j,k;

	LIS_DEBUG_FUNC_IN;

	/* C = A*B */

	if( op==LIS_INS_VALUE )
	{
		for(j=0;j<m;j++)
		{
			for(i=0;i<l;i++)
			{
				c[i+j*ldc] = 0.0;
			}
			for(k=0;k<n;k++)
			{
				for(i=0;i<l;i++)
				{
					c[i+j*ldc] += a[i+k*lda] * b[k+j*ldb];
				}
			}
		}
	}
	else if( op==LIS_SUB_VALUE )
	{
		for(j=0;j<m;j++)
		{
			for(k=0;k<n;k++)
			{
				for(i=0;i<l;i++)
				{
					c[i+j*ldc] -= a[i+k*lda] * b[k+j*ldb];
				}
			}
		}
	}
	else
	{
		for(j=0;j<m;j++)
		{
			for(k=0;k<n;k++)
			{
				for(i=0;i<l;i++)
				{
					c[i+j*ldc] += a[i+k*lda] * b[k+j*ldb];
				}
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_ge"
LIS_INT lis_array_ge(LIS_INT n, LIS_SCALAR *a)
{
	LIS_INT i,j,k;
	LIS_SCALAR t,*lu;

	LIS_DEBUG_FUNC_IN;

	/* compute inverse matrix with Gaussian elimination */	

	lu = (LIS_SCALAR *)lis_malloc(n*n*sizeof(LIS_SCALAR), "lis_array_ge::lu");
	memcpy(lu,a,n*n*sizeof(LIS_SCALAR));
	for(k=0;k<n;k++)
	{
		lu[k+k*n] = 1.0 / lu[k+k*n];
		for(i=k+1;i<n;i++)
		{
			t = lu[i+k*n] * lu[k+k*n];
			for(j=k+1;j<n;j++)
			{
				lu[i+j*n] -= t * lu[k+j*n];
			}
			lu[i+k*n] = t;
		}
	}
	for(k=0;k<n;k++)
	{
		for(i=0;i<n;i++)
		{
			 t = (i==k);
			 for(j=0;j<i;j++)
			 {
				 t -= lu[i+j*n] * a[j+k*n];
			 }
			 a[i+k*n] = t;
		}
		for(i=n-1;i>=0;i--)
		{
			t = a[i+k*n];
			for(j=i+1;j<n;j++)
			{
				t -= lu[i+j*n] * a[j+k*n];
			}
			a[k*n+i] = t * lu[i+i*n];
		}
	}
	lis_free(lu);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_solve"
LIS_INT lis_array_solve(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *x, LIS_SCALAR *w)
{
	LIS_INT i,j,k;
	LIS_SCALAR t;

	LIS_DEBUG_FUNC_IN;

	for(i=0;i<n*n;i++) w[i] = a[i];

	switch( n )
	{
	case 1:
		x[0] = b[0] / w[0];
		break;
	case 2:
		w[0]  = 1.0 / w[0];
		w[1] *= w[0];
		w[3] -= w[1] * w[2];
		w[3]  = 1.0 / w[3];
		/* forward sub */
		x[0] = b[0];
		x[1] = b[1] - w[1] * x[0];
		/* backward sub */
		x[1] *= w[3];
		x[0] -= w[2] * x[1];
		x[0] *= w[0];
		break;
	default:
		for(k=0;k<n;k++)
		{
			w[k+k*n] = 1.0 / w[k+k*n];
			for(i=k+1;i<n;i++)
			{
				t = w[i+k*n] * w[k+k*n];
				for(j=k+1;j<n;j++)
				{
					w[i+j*n] -= t * w[k+j*n];
				}
				w[i+k*n] = t;
			}
		}

		/* forward sub */
		for(i=0;i<n;i++)
		{
			x[i] = b[i];
			for(j=0;j<i;j++)
			{
				x[i] -= w[i+j*n] * x[j];
			}
		}
		/* backward sub */
		for(i=n-1;i>=0;i--)
		{
			for(j=i+1;j<n;j++)
			{
				x[i] -= w[i+j*n] * x[j];
			}
			x[i] *= w[i+i*n];
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_array_cgs"
LIS_INT lis_array_cgs(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r)
{
  LIS_INT i, j, k; 
  LIS_SCALAR *a_k;
  LIS_REAL tol, nrm2;

  LIS_DEBUG_FUNC_IN;

  tol = 1e-12;

  a_k = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR), "lis_array_cgs::a_k");

  for (i=0;i<n*n;i++)
    {
    q[i] = 0.0;
    r[i] = 0.0;
    }

  for (k=0;k<n;k++)
    {
      for (i=0;i<n;i++)
	{
	  a_k[i] = a[i+k*n];
	}
      for (j=0;j<k;j++)
	{
	  r[j+k*n] = 0; 
	  for (i=0;i<n;i++)
	    {
	      r[j+k*n] += q[i+j*n] * a[i+k*n];
	    }
	  for (i=0;i<n;i++)
	    {
	      a_k[i] -= r[j+k*n] * q[i+j*n];
	    }
	}
      lis_array_nrm2(n, &a_k[0], &nrm2);
      r[k+k*n] = nrm2;
      if (nrm2<tol) break; 
      for (i=0;i<n;i++)
	{
	  q[i+k*n] = a_k[i] / nrm2;
	}
      
    }

  lis_free(a_k);


  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
} 

#undef __FUNC__
#define __FUNC__ "lis_array_mgs"
LIS_INT lis_array_mgs(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r)
{
  LIS_INT i, j, k; 
  LIS_SCALAR *a_j;
  LIS_REAL tol, nrm2;

  LIS_DEBUG_FUNC_IN;

  tol = 1e-12;

  a_j = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR), "lis_array_mgs::a_j");

  for (i=0;i<n*n;i++)
    {
    q[i] = 0.0;
    r[i] = 0.0;
    }

  for (j=0;j<n;j++)
    {
      for (i=0;i<n;i++)
	{
	  a_j[i] = a[i+j*n];
	}
      lis_array_nrm2(n, &a_j[0], &nrm2);
      r[j+j*n] = nrm2;
      for (i=0;i<n;i++)
	{
	  if (nrm2<tol) break; 
	  q[i+j*n] = a_j[i] / nrm2;
	}
      for (k=j+1;k<n;k++)
	{
	  r[j+k*n] = 0; 
	  for (i=0;i<n;i++)
	    {
	      r[j+k*n] += q[i+j*n] * a[i+k*n];
	    }
	  for (i=0;i<n;i++)
	    {
	      a[i+k*n] -= r[j+k*n] * q[i+j*n];
	    }
	}
    }
  lis_free(a_j);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__

#ifndef _COMPLEX
/*
 * Return LIS_TRUE when A has converged to a real quasi-triangular
 * form suitable for the recycling eigensystem extractors.
 *
 * Isolated nonzero first-subdiagonal entries represent real 2x2
 * Schur blocks.  Entries below the first subdiagonal must have
 * converged to zero at the requested relative tolerance.
 */
static LIS_INT
lis_array_qr_real_quasi_schur(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_REAL tol,
    LIS_REAL *residual)
{
    LIS_INT i,j;
    LIS_INT active;
    LIS_INT previous_active;

    LIS_REAL scale;
    LIS_REAL value;
    LIS_REAL max_below;
    LIS_REAL finite_test;

    max_below =
        0.0;


    /*
     * Reject NaN and infinity before classifying the block
     * structure.  For every finite real scalar x, x-x is zero;
     * NaN or infinity produces a non-finite result.
     */
    for(i=0;i<n*n;i++)
    {
        finite_test =
            (LIS_REAL)a[i];

        finite_test =
            finite_test
            -finite_test;

        if(
            finite_test
            !=finite_test )
        {
            *residual =
                finite_test;

            return LIS_FALSE;
        }
    }


    /*
     * A real Schur matrix is upper quasi-triangular.
     */
    for(j=0;j<n;j++)
    {
        for(i=j+2;i<n;i++)
        {
            scale =
                1.0
                +fabs(
                    (LIS_REAL)a[
                        j+j*n])
                +fabs(
                    (LIS_REAL)a[
                        i+i*n]);

            value =
                fabs(
                    (LIS_REAL)a[
                        i+j*n])
                /scale;

            if( value!=value )
            {
                *residual =
                    value;

                return LIS_FALSE;
            }

            if( value>max_below )
            {
                max_below =
                    value;
            }
        }
    }


    if( max_below>=tol )
    {
        *residual =
            max_below;

        return LIS_FALSE;
    }


    /*
     * Isolated active first-subdiagonal entries are valid 2x2
     * real Schur blocks.  Adjacent active entries would create
     * an unreduced block of order greater than two.
     */
    previous_active =
        LIS_FALSE;

    for(i=0;i<n-1;i++)
    {
        scale =
            1.0
            +fabs(
                (LIS_REAL)a[
                    i+i*n])
            +fabs(
                (LIS_REAL)a[
                    (i+1)+(i+1)*n]);

        value =
            fabs(
                (LIS_REAL)a[
                    (i+1)+i*n])
            /scale;

        if( value!=value )
        {
            *residual =
                value;

            return LIS_FALSE;
        }

        active =
            value>=tol;

        if(
            active
            &&
            previous_active )
        {
            *residual =
                value;

            return LIS_FALSE;
        }

        previous_active =
            active;
    }


    *residual =
        max_below;

    return LIS_TRUE;
}
#endif


#define __FUNC__ "lis_array_qr"
static LIS_INT lis_array_qr_core(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_SCALAR *q,
    LIS_SCALAR *r,
    LIS_INT maxiter,
    LIS_INT quasi_schur,
    LIS_INT *qriter,
    LIS_REAL *qrerr)
{
    LIS_INT i,j,k,iter;
    LIS_INT converged;
    LIS_SCALAR *a0;
    LIS_REAL tol,err;

    LIS_DEBUG_FUNC_IN;

    tol = 1e-12;

    a0 = (LIS_SCALAR *)lis_malloc(
        n*n*sizeof(LIS_SCALAR),
        "lis_array_qr::a0");

    iter = 0;
    converged = LIS_FALSE;

    while( iter<maxiter )
    {
        iter++;

        lis_array_cgs(
            n,
            a,
            q,
            r);

        for(j=0;j<n;j++)
        {
            for(i=0;i<n;i++)
            {
                a[i+j*n] = 0;

                for(k=0;k<n;k++)
                {
                    a[i+j*n] +=
                        r[i+k*n]
                        *q[k+j*n];
                }
            }
        }

        err = fabs(a[1]);

        if( quasi_schur )
        {
#ifndef _COMPLEX
            if(
                lis_array_qr_real_quasi_schur(
                    n,
                    a,
                    tol,
                    &err) )
            {
                converged =
                    LIS_TRUE;

                break;
            }
#else
            if( err<tol )
            {
                converged =
                    LIS_TRUE;

                break;
            }
#endif
        }
        else if( err<tol )
        {
            converged =
                LIS_TRUE;

            break;
        }
    }

    /*
     * For the bounded quasi-Schur path, exhaustion of maxiter
     * must remain distinguishable from successful convergence
     * even when a[1] happened to deflate earlier.
     */
    if(
        quasi_schur
        &&
        !converged
        &&
        err==err
        &&
        err<tol )
    {
        err =
            tol;
    }

    *qriter = iter;
    *qrerr = err;

    lis_free(a0);

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}



#ifndef _COMPLEX

static LIS_INT
lis_array_hqr_reflector(
    LIS_REAL *v,
    LIS_INT n,
    LIS_REAL *beta,
    LIS_REAL *alpha)
{
    LIS_INT i;
    LIS_REAL norm = 0.0;
    LIS_REAL vv = 0.0;
    LIS_REAL x;

    for(i=0;i<n;i++)
    {
        x = v[i];
        if( (x-x)!=(x-x) ) return LIS_BREAKDOWN;
        norm += x*x;
    }

    norm = sqrt(norm);
    if( (norm-norm)!=(norm-norm) ) return LIS_BREAKDOWN;

    if( norm==0.0 )
    {
        *alpha = 0.0;
        *beta = 0.0;
        return LIS_SUCCESS;
    }

    *alpha = (v[0]>=0.0) ? -norm : norm;
    v[0] -= *alpha;

    for(i=0;i<n;i++) vv += v[i]*v[i];

    if( (vv-vv)!=(vv-vv) ) return LIS_BREAKDOWN;

    *beta = (vv==0.0) ? 0.0 : 2.0/vv;

    return LIS_SUCCESS;
}


static LIS_INT
lis_array_hqr_hessenberg(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_REAL *v)
{
    LIS_INT i,j,k,m,err;
    LIS_REAL alpha,beta,dot;

    for(k=0;k<n-2;k++)
    {
        m = n-k-1;

        for(i=0;i<m;i++)
            v[i] = (LIS_REAL)a[(k+1+i)+k*n];

        err = lis_array_hqr_reflector(v,m,&beta,&alpha);
        if( err ) return err;
        if( beta==0.0 ) continue;

        for(j=k;j<n;j++)
        {
            dot = 0.0;
            for(i=0;i<m;i++)
                dot += v[i]*(LIS_REAL)a[(k+1+i)+j*n];

            dot *= beta;

            for(i=0;i<m;i++)
                a[(k+1+i)+j*n] -= (LIS_SCALAR)(v[i]*dot);
        }

        for(i=0;i<n;i++)
        {
            dot = 0.0;
            for(j=0;j<m;j++)
                dot += (LIS_REAL)a[i+(k+1+j)*n]*v[j];

            dot *= beta;

            for(j=0;j<m;j++)
                a[i+(k+1+j)*n] -= (LIS_SCALAR)(dot*v[j]);
        }

        for(i=k+2;i<n;i++)
            a[i+k*n] = (LIS_SCALAR)0.0;
    }

    return LIS_SUCCESS;
}


static LIS_INT
lis_array_hqr_householder_qr(
    LIS_INT n,
    LIS_SCALAR *work,
    LIS_SCALAR *q,
    LIS_SCALAR *r,
    LIS_REAL *v)
{
    LIS_INT i,j,k,m,err;
    LIS_REAL alpha,beta,dot;

    for(j=0;j<n;j++)
    {
        for(i=0;i<n;i++)
        {
            r[i+j*n] = work[i+j*n];
            q[i+j*n] =
                (i==j) ? (LIS_SCALAR)1.0 : (LIS_SCALAR)0.0;
        }
    }

    for(k=0;k<n;k++)
    {
        m = n-k;

        for(i=0;i<m;i++)
            v[i] = (LIS_REAL)r[(k+i)+k*n];

        err = lis_array_hqr_reflector(v,m,&beta,&alpha);
        if( err ) return err;
        if( beta==0.0 ) continue;

        for(j=k;j<n;j++)
        {
            dot = 0.0;
            for(i=0;i<m;i++)
                dot += v[i]*(LIS_REAL)r[(k+i)+j*n];

            dot *= beta;

            for(i=0;i<m;i++)
                r[(k+i)+j*n] -= (LIS_SCALAR)(v[i]*dot);
        }

        for(i=0;i<n;i++)
        {
            dot = 0.0;
            for(j=0;j<m;j++)
                dot += (LIS_REAL)q[i+(k+j)*n]*v[j];

            dot *= beta;

            for(j=0;j<m;j++)
                q[i+(k+j)*n] -= (LIS_SCALAR)(dot*v[j]);
        }

        r[k+k*n] = (LIS_SCALAR)alpha;

        for(i=k+1;i<n;i++)
            r[i+k*n] = (LIS_SCALAR)0.0;
    }

    return LIS_SUCCESS;
}


static void
lis_array_hqr_deflate(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_REAL tol)
{
    LIS_INT i;
    LIS_REAL scale;

    for(i=1;i<n;i++)
    {
        scale =
            1.0
            + fabs((LIS_REAL)a[(i-1)+(i-1)*n])
            + fabs((LIS_REAL)a[i+i*n]);

        if( fabs((LIS_REAL)a[i+(i-1)*n])<=tol*scale )
            a[i+(i-1)*n] = (LIS_SCALAR)0.0;
    }
}


static LIS_INT
lis_array_hqr_active_block(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_INT *lo,
    LIS_INT *hi)
{
    LIS_INT i,j;

    i = n-1;

    while(i>0)
    {
        if( a[i+(i-1)*n]==(LIS_SCALAR)0.0 )
        {
            i--;
            continue;
        }

        if(
            i==1
            ||
            a[(i-1)+(i-2)*n]==(LIS_SCALAR)0.0 )
        {
            i -= 2;
            continue;
        }

        *hi = i;
        j = i-1;

        while(
            j>0
            &&
            a[j+(j-1)*n]!=(LIS_SCALAR)0.0 )
            j--;

        *lo = j;

        return LIS_TRUE;
    }

    return LIS_FALSE;
}


/*
 * One explicit shifted QR similarity step on an unreduced
 * Hessenberg block.  Complex-conjugate shifts use the real
 * quadratic H^2-sH+tI.
 */
static LIS_INT
lis_array_hqr_shift_step(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_INT lo,
    LIS_INT hi,
    LIS_INT double_shift,
    LIS_REAL shift_s,
    LIS_REAL shift_t,
    LIS_SCALAR *q,
    LIS_SCALAR *r,
    LIS_SCALAR *work,
    LIS_REAL *v)
{
    LIS_INT i,j,k,m,err;
    LIS_REAL sum;

    m = hi-lo+1;

    for(j=0;j<m;j++)
    {
        for(i=0;i<m;i++)
        {
            if( double_shift )
            {
                sum = 0.0;

                for(k=0;k<m;k++)
                    sum +=
                        (LIS_REAL)a[(lo+i)+(lo+k)*n]
                        *(LIS_REAL)a[(lo+k)+(lo+j)*n];

                sum -=
                    shift_s
                    *(LIS_REAL)a[(lo+i)+(lo+j)*n];

                if( i==j ) sum += shift_t;
            }
            else
            {
                sum = (LIS_REAL)a[(lo+i)+(lo+j)*n];
                if( i==j ) sum -= shift_s;
            }

            work[i+j*m] = (LIS_SCALAR)sum;
        }
    }

    err =
        lis_array_hqr_householder_qr(
            m,work,q,r,v);

    if( err ) return err;

    for(j=0;j<m;j++)
    {
        for(i=0;i<n;i++)
        {
            sum = 0.0;

            for(k=0;k<m;k++)
                sum +=
                    (LIS_REAL)a[i+(lo+k)*n]
                    *(LIS_REAL)q[k+j*m];

            work[i+(lo+j)*n] = (LIS_SCALAR)sum;
        }
    }

    for(j=0;j<m;j++)
        for(i=0;i<n;i++)
            a[i+(lo+j)*n] = work[i+(lo+j)*n];

    for(j=0;j<n;j++)
    {
        for(i=0;i<m;i++)
        {
            sum = 0.0;

            for(k=0;k<m;k++)
                sum +=
                    (LIS_REAL)q[k+i*m]
                    *(LIS_REAL)a[(lo+k)+j*n];

            work[(lo+i)+j*n] = (LIS_SCALAR)sum;
        }
    }

    for(j=0;j<n;j++)
        for(i=0;i<m;i++)
            a[(lo+i)+j*n] = work[(lo+i)+j*n];

    for(j=0;j<n;j++)
        for(i=j+2;i<n;i++)
            a[i+j*n] = (LIS_SCALAR)0.0;

    return LIS_SUCCESS;
}

#endif /* !_COMPLEX */


LIS_INT lis_array_qr(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_SCALAR *q,
    LIS_SCALAR *r,
    LIS_INT *qriter,
    LIS_REAL *qrerr)
{
    return lis_array_qr_core(
        n,
        a,
        q,
        r,
        100000,
        LIS_FALSE,
        qriter,
        qrerr);
}


LIS_INT lis_array_qr_limit(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_SCALAR *q,
    LIS_SCALAR *r,
    LIS_INT maxiter,
    LIS_INT *qriter,
    LIS_REAL *qrerr)
{
    LIS_INT err;

    if( n<=0 || maxiter<=0 )
    {
        return LIS_ERR_ILL_ARG;
    }

    /*
     * A 1x1 matrix is already in real Schur form.
     * Keep the bounded path from accessing a[1].
     */
    if( n==1 )
    {
        q[0] = 1.0;
        r[0] = a[0];

        *qriter = 1;
        *qrerr = 0.0;

        return LIS_SUCCESS;
    }

    err = lis_array_qr_core(
        n,
        a,
        q,
        r,
        maxiter,
        LIS_TRUE,
        qriter,
        qrerr);

    if( err )
    {
        return err;
    }

    if(
        *qrerr!=*qrerr ||
        *qrerr>=1e-12)
    {
        return LIS_BREAKDOWN;
    }

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_array_hqr_limit"
LIS_INT lis_array_hqr_limit(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_SCALAR *q,
    LIS_SCALAR *r,
    LIS_INT maxiter,
    LIS_INT *qriter,
    LIS_REAL *qrerr)
{
#ifdef _COMPLEX

    return lis_array_qr_limit(
        n,a,q,r,maxiter,qriter,qrerr);

#else

    LIS_INT i,j,iter,lo,hi,err;
    LIS_INT converged,double_shift,exceptional;

    LIS_REAL tol,residual;
    LIS_REAL aa,bb,cc,dd;
    LIS_REAL trace,det,disc,root;
    LIS_REAL lambda1,lambda2;
    LIS_REAL shift_s,shift_t,shift_scale;
    LIS_REAL finite_test;

    LIS_SCALAR *work;
    LIS_REAL *v;

    LIS_DEBUG_FUNC_IN;

    if(
        n<=0
        ||
        maxiter<=0
        ||
        a==NULL
        ||
        q==NULL
        ||
        r==NULL
        ||
        qriter==NULL
        ||
        qrerr==NULL )
    {
        LIS_DEBUG_FUNC_OUT;
        return LIS_ERR_ILL_ARG;
    }

    tol = 1.0e-12;
    iter = 0;
    residual = 0.0;
    converged = LIS_FALSE;

    for(j=0;j<n;j++)
    {
        for(i=0;i<n;i++)
        {
            q[i+j*n] =
                (i==j) ? (LIS_SCALAR)1.0 : (LIS_SCALAR)0.0;
            r[i+j*n] = a[i+j*n];
        }
    }

    if( lis_array_qr_real_quasi_schur(n,a,tol,&residual) )
    {
        *qriter = 0;
        *qrerr = residual;
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
    }

    if( residual!=residual )
    {
        *qriter = 0;
        *qrerr = residual;
        LIS_DEBUG_FUNC_OUT;
        return LIS_BREAKDOWN;
    }

    work =
        (LIS_SCALAR *)lis_malloc(
            (size_t)n*(size_t)n*sizeof(LIS_SCALAR),
            "lis_array_hqr_limit::work");

    v =
        (LIS_REAL *)lis_malloc(
            (size_t)n*sizeof(LIS_REAL),
            "lis_array_hqr_limit::v");

    if( work==NULL || v==NULL )
    {
        lis_free(work);
        lis_free(v);
        LIS_DEBUG_FUNC_OUT;
        return LIS_ERR_OUT_OF_MEMORY;
    }

    err = lis_array_hqr_hessenberg(n,a,v);

    if( err )
    {
        lis_free(work);
        lis_free(v);
        *qriter = 0;
        *qrerr = residual;
        LIS_DEBUG_FUNC_OUT;
        return err;
    }

    lis_array_hqr_deflate(n,a,tol);

    converged =
        lis_array_qr_real_quasi_schur(
            n,a,tol,&residual);

    while( !converged && iter<maxiter )
    {
        if( !lis_array_hqr_active_block(n,a,&lo,&hi) )
        {
            converged =
                lis_array_qr_real_quasi_schur(
                    n,a,tol,&residual);
            break;
        }

        aa = (LIS_REAL)a[(hi-1)+(hi-1)*n];
        bb = (LIS_REAL)a[(hi-1)+hi*n];
        cc = (LIS_REAL)a[hi+(hi-1)*n];
        dd = (LIS_REAL)a[hi+hi*n];

        finite_test =
            aa-aa + bb-bb + cc-cc + dd-dd;

        if( finite_test!=finite_test )
        {
            residual = finite_test;
            break;
        }

        trace = aa+dd;
        det = aa*dd-bb*cc;
        disc = trace*trace-4.0*det;

        shift_scale =
            1.0
            +fabs(aa)
            +fabs(bb)
            +fabs(cc)
            +fabs(dd);

        double_shift = LIS_FALSE;
        exceptional = (((iter+1)%32)==0);
        shift_s = 0.0;
        shift_t = 0.0;

        if( disc>=0.0 )
        {
            root = sqrt(disc);
            lambda1 = 0.5*(trace+root);
            lambda2 = 0.5*(trace-root);

            shift_s =
                (
                    fabs(lambda1-dd)
                    <=
                    fabs(lambda2-dd)
                )
                ? lambda1
                : lambda2;

            if(
                fabs(shift_s-dd)<=1.0e-14*shift_scale
                &&
                fabs(cc)>1.0e-14*shift_scale )
                exceptional = LIS_TRUE;

            if( exceptional )
                shift_s =
                    dd
                    +0.75*fabs(cc)
                    +0.25*fabs(aa-dd);
        }
        else if( exceptional )
        {
            shift_s =
                dd
                +0.75*fabs(cc)
                +0.25*fabs(aa-dd);
        }
        else
        {
            double_shift = LIS_TRUE;
            shift_s = trace;
            shift_t = det;
        }

        finite_test =
            shift_s-shift_s
            +shift_t-shift_t;

        if( finite_test!=finite_test )
        {
            residual = finite_test;
            break;
        }

        err =
            lis_array_hqr_shift_step(
                n,a,lo,hi,double_shift,
                shift_s,shift_t,q,r,work,v);

        if( err )
        {
            residual = (LIS_REAL)err;
            break;
        }

        iter++;

        lis_array_hqr_deflate(n,a,tol);

        converged =
            lis_array_qr_real_quasi_schur(
                n,a,tol,&residual);

        if( residual!=residual )
            break;
    }

    *qriter = iter;
    *qrerr = residual;

    lis_free(work);
    lis_free(v);

    if(
        !converged
        ||
        residual!=residual
        ||
        residual>=tol )
    {
        LIS_DEBUG_FUNC_OUT;
        return LIS_BREAKDOWN;
    }

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;

#endif
}
