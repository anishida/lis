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
	free(lu);

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
#define __FUNC__ "lis_array_qr"
LIS_INT lis_array_qr(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r, LIS_INT *qriter, LIS_REAL *qrerr)
{
  LIS_INT i, j, k, iter, maxiter; 
  LIS_SCALAR *a0;
  LIS_REAL tol, err;

  LIS_DEBUG_FUNC_IN;

  maxiter = 100000;
  tol = 1e-12;

  a0 = (LIS_SCALAR *)lis_malloc(n*n*sizeof(LIS_SCALAR), "lis_array_qr::a0");
  iter = 0;
  while (iter < maxiter)
    {
      iter = iter + 1;
      lis_array_cgs(n, a, q, r);
      for (j=0;j<n;j++)
	{
	  for (i=0;i<n;i++)
	    {
	      a[i+j*n] = 0;
	      for (k=0;k<n;k++)
		{
		  a[i+j*n] += r[i+k*n] * q[k+j*n];
		}
	    }
	}
      err = fabs(a[1]);
      if (err<tol) break;
    }

  *qriter = iter;
  *qrerr = err;

  lis_free(a0);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}
