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

#include "lislib.h"

/*
 * Shared internal dense eigensystem support for block eigensolvers.
 *
 * The first implementation is the existing real symmetric Jacobi
 * projected eigensolver extracted into shared internal infrastructure.
 */

#undef __FUNC__
#define __FUNC__ "lis_eigen_projected_solve"
LIS_INT lis_eigen_projected_solve(
	LIS_INT n,
	LIS_SCALAR *a,
	LIS_SCALAR *eval,
	LIS_SCALAR *evec,
	LIS_INT *niter,
	LIS_REAL *error)
{
#ifdef _COMPLEX
	LIS_SETERR(
		LIS_ERR_NOT_IMPLEMENTED,
		"projected eigensystem solve is not yet implemented "
		"for complex arithmetic\n");
	return LIS_ERR_NOT_IMPLEMENTED;
#else
	LIS_INT i,j,p,q,sweep;
	LIS_INT max_sweeps;
	LIS_REAL app,aqq,apq;
	LIS_REAL akp,akq;
	LIS_REAL vip,viq;
	LIS_REAL tau,t,c,s;
	LIS_REAL off,maxoff,scale,threshold;

	LIS_DEBUG_FUNC_IN;

	if( n<=0 || a==NULL || eval==NULL || evec==NULL )
	{
		LIS_SETERR(
			LIS_ERR_ILL_ARG,
			"invalid projected eigensystem arguments\n");
		return LIS_ERR_ILL_ARG;
	}

	for(i=0;i<n*n;i++)
	{
		evec[i] = 0.0;
	}

	for(i=0;i<n;i++)
	{
		evec[i*n+i] = 1.0;
	}

	scale = 0.0;
	for(i=0;i<n;i++)
	{
		if( fabs(a[i*n+i]) > scale )
		{
			scale = fabs(a[i*n+i]);
		}
	}

	if( scale < 1.0 )
	{
		scale = 1.0;
	}

	threshold = 1.0e-14*scale;
	max_sweeps = 50;

	maxoff = 0.0;

	for(sweep=0;sweep<max_sweeps;sweep++)
	{
		maxoff = 0.0;

		for(p=0;p<n-1;p++)
		{
			for(q=p+1;q<n;q++)
			{
				apq = a[p*n+q];
				off = fabs(apq);

				if( off > maxoff )
				{
					maxoff = off;
				}

				if( off <= threshold )
				{
					continue;
				}

				app = a[p*n+p];
				aqq = a[q*n+q];

				tau = (aqq-app)/(2.0*apq);

				if( tau>=0.0 )
				{
					t = 1.0/
						(tau+sqrt(1.0+tau*tau));
				}
				else
				{
					t = -1.0/
						(-tau+sqrt(1.0+tau*tau));
				}

				c = 1.0/sqrt(1.0+t*t);
				s = t*c;

				for(j=0;j<n;j++)
				{
					if( j==p || j==q )
					{
						continue;
					}

					akp = a[j*n+p];
					akq = a[j*n+q];

					a[j*n+p] = c*akp-s*akq;
					a[p*n+j] = a[j*n+p];

					a[j*n+q] = s*akp+c*akq;
					a[q*n+j] = a[j*n+q];
				}

				a[p*n+p] = app-t*apq;
				a[q*n+q] = aqq+t*apq;
				a[p*n+q] = 0.0;
				a[q*n+p] = 0.0;

				for(j=0;j<n;j++)
				{
					vip = evec[j*n+p];
					viq = evec[j*n+q];

					evec[j*n+p] = c*vip-s*viq;
					evec[j*n+q] = s*vip+c*viq;
				}
			}
		}

		if( maxoff <= threshold )
		{
			break;
		}
	}

	for(i=0;i<n;i++)
	{
		eval[i] = a[i*n+i];
	}

	if( niter )
	{
		*niter = sweep+1;
	}

	if( error )
	{
		*error = maxoff;
	}

	if( maxoff > threshold )
	{
		LIS_SETERR(
			LIS_MAXITER,
			"projected Jacobi eigensolver did not converge\n");
		return LIS_MAXITER;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}


#undef __FUNC__
#define __FUNC__ "lis_eigen_block_gram"
LIS_INT lis_eigen_block_gram(
    LIS_INT nx,
    LIS_VECTOR X[],
    LIS_INT ny,
    LIS_VECTOR Y[],
    LIS_SCALAR *gram)
{
    LIS_INT i,j,err;

    LIS_DEBUG_FUNC_IN;

    if( nx<=0 || ny<=0 ||
        X==NULL || Y==NULL ||
        gram==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block Gram arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<nx;i++)
    {
        if( X[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null left block vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(j=0;j<ny;j++)
    {
        if( Y[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null right block vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            err = lis_vector_dot(
                X[i],
                Y[j],
                &gram[i*ny+j]);

            if( err )
            {
                return err;
            }
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_eigen_block_bgram"
LIS_INT lis_eigen_block_bgram(
    LIS_MATRIX B,
    LIS_INT nx,
    LIS_VECTOR X[],
    LIS_INT ny,
    LIS_VECTOR Y[],
    LIS_VECTOR BY[],
    LIS_SCALAR *gram)
{
    LIS_INT i,j,err;

    LIS_DEBUG_FUNC_IN;

    if( nx<=0 || ny<=0 ||
        X==NULL || Y==NULL ||
        gram==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block B-Gram arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    if( B==NULL )
    {
        return lis_eigen_block_gram(
            nx,
            X,
            ny,
            Y,
            gram);
    }

    if( BY==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "missing B-Gram work vectors\n");

        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<nx;i++)
    {
        if( X[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null left block vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(j=0;j<ny;j++)
    {
        if( Y[j]==NULL || BY[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null right or B-work block vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    err = lis_matvec_block(
        B,
        ny,
        Y,
        BY);

    if( err )
    {
        return err;
    }

    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            err = lis_vector_dot(
                X[i],
                BY[j],
                &gram[i*ny+j]);

            if( err )
            {
                return err;
            }
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}
