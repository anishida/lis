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
#include <string.h>

#ifdef _COMPLEX
#include <complex.h>
#endif

#include "lislib.h"

/*
 * Shared internal dense eigensystem support for block eigensolvers.
 *
 * The projected dense eigensystem helper supports real symmetric
 * and complex Hermitian matrices.
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
        LIS_INT i,j,p,q,sweep;
        LIS_INT max_sweeps;
        LIS_REAL app,aqq;
        LIS_REAL mag,tau,t,c,sabs;
        LIS_REAL off,maxoff,scale,threshold;
        LIS_SCALAR apq,phase,s;
        LIS_SCALAR akp,akq;
        LIS_SCALAR vip,viq;

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
                app = creal(a[i*n+i]);

                if( fabs(app)>scale )
                {
                        scale = fabs(app);
                }
        }

        if( scale<1.0 )
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
                                off = cabs(apq);

                                if( off>maxoff )
                                {
                                        maxoff = off;
                                }

                                if( off<=threshold )
                                {
                                        continue;
                                }

                                app = creal(a[p*n+p]);
                                aqq = creal(a[q*n+q]);

                                mag = off;
                                phase = apq/mag;

                                tau =
                                        (aqq-app)/
                                        (2.0*mag);

                                if( tau>=0.0 )
                                {
                                        t =
                                                1.0/
                                                (
                                                        tau+
                                                        sqrt(
                                                                1.0+
                                                                tau*tau)
                                                );
                                }
                                else
                                {
                                        t =
                                                -1.0/
                                                (
                                                        -tau+
                                                        sqrt(
                                                                1.0+
                                                                tau*tau)
                                                );
                                }

                                c =
                                        1.0/
                                        sqrt(
                                                1.0+
                                                t*t);

                                sabs = t*c;
                                s = sabs*phase;

                                for(j=0;j<n;j++)
                                {
                                        if( j==p || j==q )
                                        {
                                                continue;
                                        }

                                        akp =
                                                a[j*n+p];

                                        akq =
                                                a[j*n+q];

                                        a[j*n+p] =
                                                c*akp-
                                                conj(s)*akq;

                                        a[p*n+j] =
                                                conj(
                                                        a[j*n+p]);

                                        a[j*n+q] =
                                                s*akp+
                                                c*akq;

                                        a[q*n+j] =
                                                conj(
                                                        a[j*n+q]);
                                }

                                a[p*n+p] =
                                        app-t*mag;

                                a[q*n+q] =
                                        aqq+t*mag;

                                a[p*n+q] = 0.0;
                                a[q*n+p] = 0.0;

                                for(j=0;j<n;j++)
                                {
                                        vip =
                                                evec[j*n+p];

                                        viq =
                                                evec[j*n+q];

                                        evec[j*n+p] =
                                                c*vip-
                                                conj(s)*viq;

                                        evec[j*n+q] =
                                                s*vip+
                                                c*viq;
                                }
                        }
                }

                if( maxoff<=threshold )
                {
                        break;
                }
        }

        for(i=0;i<n;i++)
        {
                eval[i] =
                        creal(
                                a[i*n+i]);
        }

        if( niter )
        {
                *niter = sweep+1;
        }

        if( error )
        {
                *error = maxoff;
        }

        if( maxoff>threshold )
        {
                LIS_SETERR(
                        LIS_MAXITER,
                        "projected Hermitian Jacobi eigensolver "
                        "did not converge\n");

                return LIS_MAXITER;
        }

        LIS_DEBUG_FUNC_OUT;

        return LIS_SUCCESS;
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

#undef __FUNC__
#define __FUNC__ "lis_eigen_block_orthogonalize"
LIS_INT lis_eigen_block_orthogonalize(
    LIS_INT nq,
    LIS_VECTOR Q[],
    LIS_INT nw,
    LIS_VECTOR W[],
    LIS_INT passes,
    LIS_SCALAR *coeff,
    LIS_SCALAR *work)
{
    LIS_INT i,j,k,p,err;
    LIS_SCALAR alpha;

    LIS_DEBUG_FUNC_IN;

    if( nq<=0 || nw<=0 ||
        Q==NULL || W==NULL ||
        passes<=0 ||
        coeff==NULL || work==NULL ||
        coeff==work )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block orthogonalization arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<nq;i++)
    {
        if( Q[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null orthogonalization basis vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(j=0;j<nw;j++)
    {
        if( W[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null orthogonalization candidate vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(k=0;k<nq*nw;k++)
    {
        coeff[k] = 0.0;
    }

    for(p=0;p<passes;p++)
    {
        err = lis_eigen_block_gram(
            nq,
            Q,
            nw,
            W,
            work);

        if( err )
        {
            return err;
        }

        for(i=0;i<nq;i++)
        {
            for(j=0;j<nw;j++)
            {
                k = i*nw+j;

                coeff[k] += work[k];

                alpha = -work[k];

                err = lis_vector_axpy(
                    alpha,
                    Q[i],
                    W[j]);

                if( err )
                {
                    return err;
                }
            }
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_eigen_block_borthogonalize"
LIS_INT lis_eigen_block_borthogonalize(
    LIS_MATRIX B,
    LIS_INT nq,
    LIS_VECTOR Q[],
    LIS_INT nw,
    LIS_VECTOR W[],
    LIS_VECTOR BW[],
    LIS_INT passes,
    LIS_SCALAR *coeff,
    LIS_SCALAR *work)
{
    LIS_INT i,j,k,p,err;
    LIS_SCALAR alpha;

    LIS_DEBUG_FUNC_IN;

    if( B==NULL )
    {
        return lis_eigen_block_orthogonalize(
            nq,
            Q,
            nw,
            W,
            passes,
            coeff,
            work);
    }

    if( nq<=0 || nw<=0 ||
        Q==NULL || W==NULL ||
        BW==NULL ||
        passes<=0 ||
        coeff==NULL || work==NULL ||
        coeff==work )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block B-orthogonalization arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<nq;i++)
    {
        if( Q[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null B-orthogonalization basis vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(j=0;j<nw;j++)
    {
        if( W[j]==NULL || BW[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null B-orthogonalization candidate or work vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(k=0;k<nq*nw;k++)
    {
        coeff[k] = 0.0;
    }

    for(p=0;p<passes;p++)
    {
        err = lis_eigen_block_bgram(
            B,
            nq,
            Q,
            nw,
            W,
            BW,
            work);

        if( err )
        {
            return err;
        }

        for(i=0;i<nq;i++)
        {
            for(j=0;j<nw;j++)
            {
                k = i*nw+j;

                coeff[k] += work[k];

                alpha = -work[k];

                err = lis_vector_axpy(
                    alpha,
                    Q[i],
                    W[j]);

                if( err )
                {
                    return err;
                }
            }
        }
    }

    /*
     * Keep BW synchronized with the final orthogonalized W.
     */
    err = lis_matvec_block(
        B,
        nw,
        W,
        BW);

    if( err )
    {
        return err;
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}

static LIS_REAL
lis_eigen_block_scalar_real(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (LIS_REAL)creal(value);
#else
    return (LIS_REAL)value;
#endif
}


#undef __FUNC__
#define __FUNC__ "lis_eigen_block_linear_combination"
LIS_INT lis_eigen_block_linear_combination(
    LIS_INT nx,
    LIS_VECTOR X[],
    LIS_INT ny,
    LIS_VECTOR Y[],
    LIS_SCALAR *coeff,
    LIS_INT ldcoeff)
{
    LIS_INT i,j,k,err;

    LIS_DEBUG_FUNC_IN;

    if( nx<=0 || ny<=0 ||
        X==NULL || Y==NULL ||
        coeff==NULL ||
        ldcoeff<ny )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block linear-combination arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<nx;i++)
    {
        if( X[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null block linear-combination input vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(j=0;j<ny;j++)
    {
        if( Y[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null block linear-combination output vector\n");

            return LIS_ERR_ILL_ARG;
        }

        for(k=0;k<nx;k++)
        {
            if( Y[j]==X[k] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "block linear combination requires "
                    "out-of-place output vectors\n");

                return LIS_ERR_ILL_ARG;
            }
        }

        for(k=0;k<j;k++)
        {
            if( Y[j]==Y[k] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "block linear combination requires "
                    "distinct output vectors\n");

                return LIS_ERR_ILL_ARG;
            }
        }
    }

    for(j=0;j<ny;j++)
    {
        err = lis_vector_copy(
            X[0],
            Y[j]);

        if( err )
        {
            return err;
        }

        err = lis_vector_scale(
            coeff[j],
            Y[j]);

        if( err )
        {
            return err;
        }

        for(i=1;i<nx;i++)
        {
            err = lis_vector_axpy(
                coeff[i*ldcoeff+j],
                X[i],
                Y[j]);

            if( err )
            {
                return err;
            }
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


static LIS_INT
lis_eigen_block_rank_transform(
    LIS_INT n,
    LIS_REAL rank_tol,
    LIS_SCALAR *gram,
    LIS_SCALAR *eval,
    LIS_SCALAR *evec,
    LIS_SCALAR *transform,
    LIS_INT *active_rank)
{
    LIS_INT i,j,k,best;
    LIS_INT niter;
    LIS_INT err;
    LIS_INT rank;

    LIS_REAL error;
    LIS_REAL best_value;
    LIS_REAL value;
    LIS_REAL lambda_max;
    LIS_REAL threshold;
    LIS_REAL inverse_root;

    LIS_SCALAR swap;

    if( n<=0 ||
        !isfinite((double)rank_tol) ||
        rank_tol<0.0 ||
        rank_tol>=1.0 ||
        gram==NULL ||
        eval==NULL ||
        evec==NULL ||
        transform==NULL ||
        active_rank==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block rank-transform arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    *active_rank = 0;

    niter = 0;
    error = 0.0;

    err = lis_eigen_projected_solve(
        n,
        gram,
        eval,
        evec,
        &niter,
        &error);

    if( err )
    {
        return err;
    }

    /*
     * Sort the Hermitian Gram eigenpairs in descending order.
     *
     * Keeping the retained eigenpairs at the front gives the
     * following stages a deterministic active subspace layout.
     */
    for(i=0;i<n-1;i++)
    {
        best = i;

        best_value =
            lis_eigen_block_scalar_real(
                eval[i]);

        if( !isfinite((double)best_value) )
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "non-finite block Gram eigenvalue\n");

            return LIS_BREAKDOWN;
        }

        for(j=i+1;j<n;j++)
        {
            value =
                lis_eigen_block_scalar_real(
                    eval[j]);

            if( !isfinite((double)value) )
            {
                LIS_SETERR(
                    LIS_BREAKDOWN,
                    "non-finite block Gram eigenvalue\n");

                return LIS_BREAKDOWN;
            }

            if( value>best_value )
            {
                best = j;
                best_value = value;
            }
        }

        if( best!=i )
        {
            swap = eval[i];
            eval[i] = eval[best];
            eval[best] = swap;

            for(k=0;k<n;k++)
            {
                swap =
                    evec[k*n+i];

                evec[k*n+i] =
                    evec[k*n+best];

                evec[k*n+best] =
                    swap;
            }
        }
    }

    value =
        lis_eigen_block_scalar_real(
            eval[n-1]);

    if( !isfinite((double)value) )
    {
        LIS_SETERR(
            LIS_BREAKDOWN,
            "non-finite block Gram eigenvalue\n");

        return LIS_BREAKDOWN;
    }

    for(i=0;i<n*n;i++)
    {
        transform[i] = 0.0;
    }

    lambda_max =
        lis_eigen_block_scalar_real(
            eval[0]);

    if( !isfinite((double)lambda_max) )
    {
        LIS_SETERR(
            LIS_BREAKDOWN,
            "non-finite maximum block Gram eigenvalue\n");

        return LIS_BREAKDOWN;
    }

    if( lambda_max<=0.0 )
    {
        *active_rank = 0;

        return LIS_SUCCESS;
    }

    threshold =
        rank_tol*lambda_max;

    if( !isfinite((double)threshold) )
    {
        LIS_SETERR(
            LIS_BREAKDOWN,
            "block rank threshold overflowed\n");

        return LIS_BREAKDOWN;
    }

    rank = 0;

    for(j=0;j<n;j++)
    {
        value =
            lis_eigen_block_scalar_real(
                eval[j]);

        if( !isfinite((double)value) )
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "non-finite block Gram eigenvalue\n");

            return LIS_BREAKDOWN;
        }

        if( value>threshold &&
            value>0.0 )
        {
            rank++;
        }
        else
        {
            break;
        }
    }

    for(j=0;j<rank;j++)
    {
        value =
            lis_eigen_block_scalar_real(
                eval[j]);

        inverse_root =
            1.0/sqrt(value);

        if( !isfinite((double)inverse_root) )
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "block rank normalization overflowed\n");

            return LIS_BREAKDOWN;
        }

        for(i=0;i<n;i++)
        {
            transform[i*n+j] =
                evec[i*n+j] *
                inverse_root;
        }
    }

    *active_rank = rank;

    return LIS_SUCCESS;
}


static LIS_INT
lis_eigen_block_zero_outputs(
    LIS_INT capacity,
    LIS_INT start,
    LIS_VECTOR source,
    LIS_VECTOR Y[])
{
    LIS_INT j,err;

    if( capacity<0 ||
        start<0 ||
        start>capacity ||
        source==NULL ||
        Y==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block zero-output arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    for(j=start;j<capacity;j++)
    {
        if( Y[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null inactive block output vector\n");

            return LIS_ERR_ILL_ARG;
        }

        err = lis_vector_copy(
            source,
            Y[j]);

        if( err )
        {
            return err;
        }

        err = lis_vector_scale(
            0.0,
            Y[j]);

        if( err )
        {
            return err;
        }
    }

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_eigen_block_rank_deflate"
LIS_INT lis_eigen_block_rank_deflate(
    LIS_INT capacity,
    LIS_INT nominal_rank,
    LIS_VECTOR W[],
    LIS_VECTOR Q[],
    LIS_REAL rank_tol,
    LIS_INT *active_rank,
    LIS_SCALAR *gram,
    LIS_SCALAR *eval,
    LIS_SCALAR *evec,
    LIS_SCALAR *transform)
{
    LIS_INT i,j,err;
    LIS_INT rank;

    LIS_DEBUG_FUNC_IN;

    if( capacity<=0 ||
        nominal_rank<=0 ||
        nominal_rank>capacity ||
        W==NULL ||
        Q==NULL ||
        active_rank==NULL ||
        gram==NULL ||
        eval==NULL ||
        evec==NULL ||
        transform==NULL ||
        !isfinite((double)rank_tol) ||
        rank_tol<0.0 ||
        rank_tol>=1.0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block rank-deflation arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    *active_rank = 0;

    for(i=0;i<nominal_rank;i++)
    {
        if( W[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null rank-deflation candidate vector\n");

            return LIS_ERR_ILL_ARG;
        }
    }

    for(j=0;j<capacity;j++)
    {
        if( Q[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null rank-deflation output vector\n");

            return LIS_ERR_ILL_ARG;
        }

        for(i=0;i<nominal_rank;i++)
        {
            if( Q[j]==W[i] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "block rank deflation requires "
                    "out-of-place output vectors\n");

                return LIS_ERR_ILL_ARG;
            }
        }
    }

    err = lis_eigen_block_gram(
        nominal_rank,
        W,
        nominal_rank,
        W,
        gram);

    if( err )
    {
        return err;
    }

    err = lis_eigen_block_rank_transform(
        nominal_rank,
        rank_tol,
        gram,
        eval,
        evec,
        transform,
        &rank);

    if( err )
    {
        return err;
    }

    if( rank>0 )
    {
        err =
            lis_eigen_block_linear_combination(
                nominal_rank,
                W,
                rank,
                Q,
                transform,
                nominal_rank);

        if( err )
        {
            return err;
        }
    }

    err = lis_eigen_block_zero_outputs(
        capacity,
        rank,
        W[0],
        Q);

    if( err )
    {
        return err;
    }

    *active_rank = rank;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_eigen_block_brank_deflate"
LIS_INT lis_eigen_block_brank_deflate(
    LIS_MATRIX B,
    LIS_INT capacity,
    LIS_INT nominal_rank,
    LIS_VECTOR W[],
    LIS_VECTOR BW[],
    LIS_VECTOR Q[],
    LIS_VECTOR BQ[],
    LIS_REAL rank_tol,
    LIS_INT *active_rank,
    LIS_SCALAR *gram,
    LIS_SCALAR *eval,
    LIS_SCALAR *evec,
    LIS_SCALAR *transform)
{
    LIS_INT i,j,k,err;
    LIS_INT rank;

    LIS_DEBUG_FUNC_IN;

    if( B==NULL ||
        capacity<=0 ||
        nominal_rank<=0 ||
        nominal_rank>capacity ||
        W==NULL ||
        BW==NULL ||
        Q==NULL ||
        BQ==NULL ||
        active_rank==NULL ||
        gram==NULL ||
        eval==NULL ||
        evec==NULL ||
        transform==NULL ||
        !isfinite((double)rank_tol) ||
        rank_tol<0.0 ||
        rank_tol>=1.0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block B-rank-deflation arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    *active_rank = 0;

    for(i=0;i<nominal_rank;i++)
    {
        if( W[i]==NULL ||
            BW[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null B-rank-deflation candidate or work vector\n");

            return LIS_ERR_ILL_ARG;
        }

        for(k=0;k<nominal_rank;k++)
        {
            if( BW[i]==W[k] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "B-rank deflation requires separate BW vectors\n");

                return LIS_ERR_ILL_ARG;
            }
        }
    }

    for(j=0;j<capacity;j++)
    {
        if( Q[j]==NULL ||
            BQ[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null B-rank-deflation output vector\n");

            return LIS_ERR_ILL_ARG;
        }

        if( Q[j]==BQ[j] )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "B-rank deflation requires separate Q and BQ vectors\n");

            return LIS_ERR_ILL_ARG;
        }

        for(i=0;i<nominal_rank;i++)
        {
            if( Q[j]==W[i] ||
                Q[j]==BW[i] ||
                BQ[j]==W[i] ||
                BQ[j]==BW[i] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "B-rank deflation requires "
                    "out-of-place output vectors\n");

                return LIS_ERR_ILL_ARG;
            }
        }

        for(k=0;k<j;k++)
        {
            if( Q[j]==Q[k] ||
                BQ[j]==BQ[k] ||
                Q[j]==BQ[k] ||
                BQ[j]==Q[k] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "B-rank deflation requires distinct outputs\n");

                return LIS_ERR_ILL_ARG;
            }
        }
    }

    /*
     * This computes BW = B*W exactly once and forms W^H*BW.
     */
    err = lis_eigen_block_bgram(
        B,
        nominal_rank,
        W,
        nominal_rank,
        W,
        BW,
        gram);

    if( err )
    {
        return err;
    }

    err = lis_eigen_block_rank_transform(
        nominal_rank,
        rank_tol,
        gram,
        eval,
        evec,
        transform,
        &rank);

    if( err )
    {
        return err;
    }

    if( rank>0 )
    {
        err =
            lis_eigen_block_linear_combination(
                nominal_rank,
                W,
                rank,
                Q,
                transform,
                nominal_rank);

        if( err )
        {
            return err;
        }

        /*
         * BQ = BW * U_r * Lambda_r^(-1/2)
         *
         * Reusing BW avoids an additional B*Q block matvec.
         */
        err =
            lis_eigen_block_linear_combination(
                nominal_rank,
                BW,
                rank,
                BQ,
                transform,
                nominal_rank);

        if( err )
        {
            return err;
        }
    }

    err = lis_eigen_block_zero_outputs(
        capacity,
        rank,
        W[0],
        Q);

    if( err )
    {
        return err;
    }

    err = lis_eigen_block_zero_outputs(
        capacity,
        rank,
        BW[0],
        BQ);

    if( err )
    {
        return err;
    }

    *active_rank = rank;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}

static LIS_SCALAR
lis_eigen_block_scalar_conjugate(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return conj(value);
#else
    return value;
#endif
}


#undef __FUNC__
#define __FUNC__ "lis_eigen_block_rectangular_coupling"
LIS_INT lis_eigen_block_rectangular_coupling(
    LIS_INT nominal_rank,
    LIS_INT active_rank,
    LIS_SCALAR *eval,
    LIS_SCALAR *evec,
    LIS_SCALAR *coupling)
{
    LIS_INT i,j;

    LIS_REAL value;
    LIS_REAL root;

    LIS_DEBUG_FUNC_IN;

    if( nominal_rank<=0 ||
        active_rank<0 ||
        active_rank>nominal_rank ||
        eval==NULL ||
        evec==NULL ||
        coupling==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid rectangular coupling arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    /*
     * The output storage is nominal_rank by nominal_rank.
     *
     * Only the first active_rank rows are active.  Keeping the
     * remaining rows zero gives callers deterministic storage while
     * the logical matrix remains active_rank by nominal_rank.
     *
     * For retained Gram eigenpairs:
     *
     *     Q = W U_r Lambda_r^(-1/2)
     *
     * this helper builds:
     *
     *     C = Lambda_r^(1/2) U_r^H
     *
     * so that:
     *
     *     W = Q C
     *
     * on the retained candidate subspace.
     */
    for(i=0;i<nominal_rank*nominal_rank;i++)
    {
        coupling[i] = 0.0;
    }

    for(j=0;j<active_rank;j++)
    {
        value =
            lis_eigen_block_scalar_real(
                eval[j]);

        if( !isfinite((double)value) ||
            value<=0.0 )
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "invalid retained block Gram eigenvalue\n");

            return LIS_BREAKDOWN;
        }

        root = sqrt(value);

        if( !isfinite((double)root) )
        {
            LIS_SETERR(
                LIS_BREAKDOWN,
                "rectangular coupling normalization overflowed\n");

            return LIS_BREAKDOWN;
        }

        for(i=0;i<nominal_rank;i++)
        {
            coupling[j*nominal_rank+i] =
                root *
                lis_eigen_block_scalar_conjugate(
                    evec[i*nominal_rank+j]);
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_eigen_ritz_reconstruct"
LIS_INT lis_eigen_ritz_reconstruct(
    LIS_INT basis_rank,
    LIS_VECTOR basis[],
    LIS_INT nritz,
    LIS_SCALAR *coeff,
    LIS_INT ldcoeff,
    LIS_VECTOR ritz[])
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    if( basis_rank<=0 ||
        nritz<=0 ||
        basis==NULL ||
        coeff==NULL ||
        ritz==NULL ||
        ldcoeff<nritz )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid Ritz reconstruction arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    /*
     * Reconstruct physical Ritz vectors from projected
     * eigenvector coefficients:
     *
     *     ritz[j] =
     *         sum_i basis[i] * coeff[i*ldcoeff+j]
     *
     * Coefficients are applied directly.  Complex
     * reconstruction does not conjugate the projected
     * eigenvector coefficients.
     */
    err = lis_eigen_block_linear_combination(
        basis_rank,
        basis,
        nritz,
        ritz,
        coeff,
        ldcoeff);

    if( err )
    {
        return err;
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


/*
 * General complex dense eigensystem support.
 *
 * This path is intentionally separate from
 * lis_eigen_projected_solve(), whose contract is symmetric /
 * Hermitian.
 *
 * The general-complex path is used for small non-Hermitian
 * projected problems such as harmonic Ritz extraction.
 *
 * No LAPACK or external dense eigensolver is required.
 *
 * Storage follows the native general-complex kernel convention:
 *
 *     a[i + j*n]
 *
 * is matrix element (i,j), and right eigenvectors are returned
 * column-wise with the same leading dimension.
 */
#ifdef _COMPLEX

static LIS_REAL
zabs1(
    LIS_SCALAR z)
{
    return (LIS_REAL)cabs(z);
}


static void
identity(
    LIS_INT n,
    LIS_SCALAR *q)
{
    LIS_INT i,j;

    for(j=0;j<n;j++)
    {
        for(i=0;i<n;i++)
        {
            q[i+j*n] =
                (i==j)
                ? 1.0
                : 0.0;
        }
    }
}


static void
matmul(
    LIS_INT n,
    const LIS_SCALAR *a,
    const LIS_SCALAR *b,
    LIS_SCALAR *c)
{
    LIS_INT i,j,k;

    for(j=0;j<n;j++)
    {
        for(i=0;i<n;i++)
        {
            LIS_SCALAR s = 0.0;

            for(k=0;k<n;k++)
            {
                s +=
                    a[i+k*n]
                    *b[k+j*n];
            }

            c[i+j*n] = s;
        }
    }
}


static LIS_REAL
subdiag_norm(
    LIS_INT n,
    const LIS_SCALAR *a)
{
    LIS_INT i,j;
    LIS_REAL s = 0.0;

    for(j=0;j<n;j++)
    {
        for(i=j+1;i<n;i++)
        {
            LIS_REAL v =
                zabs1(
                    a[i+j*n]);

            s += v*v;
        }
    }

    return sqrt(s);
}


/*
 * Build a unitary Householder vector for x.
 *
 * The reflector is
 *
 *     H = I - tau v v^H
 *
 * with v[0] = 1.
 */
static LIS_INT
householder(
    LIS_INT n,
    const LIS_SCALAR *x,
    LIS_SCALAR *v,
    LIS_REAL *tau)
{
    LIS_INT i;

    LIS_REAL sigma;
    LIS_REAL x0abs;
    LIS_REAL norm;

    LIS_SCALAR alpha;
    LIS_SCALAR beta;
    LIS_SCALAR denom;

    if( n<=0 )
        return LIS_ERR_ILL_ARG;

    sigma = 0.0;

    for(i=1;i<n;i++)
    {
        LIS_REAL a =
            zabs1(x[i]);

        sigma += a*a;
    }

    x0abs =
        zabs1(x[0]);

    if( sigma==0.0 )
    {
        v[0] = 1.0;

        for(i=1;i<n;i++)
            v[i] = 0.0;

        *tau = 0.0;

        return LIS_SUCCESS;
    }

    norm =
        sqrt(
            x0abs*x0abs
            +sigma);

    if( x0abs==0.0 )
    {
        alpha = 1.0;
    }
    else
    {
        alpha =
            x[0]/x0abs;
    }

    beta =
        -alpha*norm;

    denom =
        x[0]-beta;

    if( zabs1(denom)==0.0 )
        return LIS_BREAKDOWN;

    v[0] = 1.0;

    for(i=1;i<n;i++)
    {
        v[i] =
            x[i]/denom;
    }

    {
        LIS_REAL vv = 1.0;

        for(i=1;i<n;i++)
        {
            LIS_REAL a =
                zabs1(v[i]);

            vv += a*a;
        }

        *tau =
            2.0/vv;
    }

    return LIS_SUCCESS;
}


/*
 * QR factorization using unitary complex Householder reflectors.
 *
 * A = Q R
 */
static LIS_INT
complex_householder_qr(
    LIS_INT n,
    const LIS_SCALAR *a,
    LIS_SCALAR *q,
    LIS_SCALAR *r)
{
    LIS_INT i,j,k,m;
    LIS_INT err;

    LIS_SCALAR *x;
    LIS_SCALAR *v;

    LIS_REAL tau;

    x =
        (LIS_SCALAR *)malloc(
            (size_t)n*sizeof(*x));

    v =
        (LIS_SCALAR *)malloc(
            (size_t)n*sizeof(*v));

    if( x==NULL || v==NULL )
    {
        free(x);
        free(v);
        return LIS_OUT_OF_MEMORY;
    }

    memcpy(
        r,
        a,
        (size_t)n*n*sizeof(*r));

    identity(n,q);

    for(k=0;k<n;k++)
    {
        LIS_INT row;

        m = n-k;

        for(i=0;i<m;i++)
        {
            x[i] =
                r[(k+i)+k*n];
        }

        err =
            householder(
                m,
                x,
                v,
                &tau);

        if( err )
        {
            free(x);
            free(v);
            return err;
        }

        if( tau==0.0 )
            continue;

        /*
         * R <- H R
         */
        for(j=k;j<n;j++)
        {
            LIS_SCALAR dot = 0.0;

            for(i=0;i<m;i++)
            {
                dot +=
                    conj(v[i])
                    *r[(k+i)+j*n];
            }

            dot *= tau;

            for(i=0;i<m;i++)
            {
                r[(k+i)+j*n] -=
                    v[i]*dot;
            }
        }

        /*
         * Q <- Q H
         *
         * H is Hermitian and unitary.
         */
        for(row=0;row<n;row++)
        {
            LIS_SCALAR dot = 0.0;

            for(i=0;i<m;i++)
            {
                dot +=
                    q[row+(k+i)*n]
                    *v[i];
            }

            dot *= tau;

            for(i=0;i<m;i++)
            {
                q[row+(k+i)*n] -=
                    dot*conj(v[i]);
            }
        }
    }

    free(x);
    free(v);

    return LIS_SUCCESS;
}


static LIS_SCALAR
wilkinson_shift_complex(
    LIS_INT n,
    const LIS_SCALAR *a)
{
    LIS_SCALAR aa;
    LIS_SCALAR bb;
    LIS_SCALAR cc;
    LIS_SCALAR dd;

    LIS_SCALAR tr;
    LIS_SCALAR det;
    LIS_SCALAR disc;

    LIS_SCALAR l1;
    LIS_SCALAR l2;

    if( n<2 )
        return a[0];

    aa =
        a[(n-2)+(n-2)*n];

    bb =
        a[(n-2)+(n-1)*n];

    cc =
        a[(n-1)+(n-2)*n];

    dd =
        a[(n-1)+(n-1)*n];

    tr =
        aa+dd;

    det =
        aa*dd-bb*cc;

    disc =
        csqrt(
            tr*tr
            -4.0*det);

    l1 =
        0.5*(tr+disc);

    l2 =
        0.5*(tr-disc);

    if(
        zabs1(l1-dd)
        <
        zabs1(l2-dd) )
    {
        return l1;
    }

    return l2;
}


static LIS_INT
complex_shifted_qr_schur(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_SCALAR *z,
    LIS_INT maxiter,
    LIS_INT *niter,
    LIS_REAL *final_resid)
{
    LIS_INT i,j,iter;
    LIS_INT err;

    LIS_REAL scale;
    LIS_REAL resid;
    LIS_REAL tol;

    LIS_SCALAR mu;

    LIS_SCALAR *work;
    LIS_SCALAR *q;
    LIS_SCALAR *r;
    LIS_SCALAR *rq;
    LIS_SCALAR *znew;

    if(
        n<=0 ||
        a==NULL ||
        z==NULL )
    {
        return LIS_ERR_ILL_ARG;
    }

    work =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(*work));

    q =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(*q));

    r =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(*r));

    rq =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(*rq));

    znew =
        (LIS_SCALAR *)malloc(
            (size_t)n*n*sizeof(*znew));

    if(
        work==NULL ||
        q==NULL ||
        r==NULL ||
        rq==NULL ||
        znew==NULL )
    {
        free(work);
        free(q);
        free(r);
        free(rq);
        free(znew);

        return LIS_OUT_OF_MEMORY;
    }

    identity(n,z);

    scale = 1.0;

    for(i=0;i<n*n;i++)
    {
        LIS_REAL v =
            zabs1(a[i]);

        if( v>scale )
            scale = v;
    }

    tol =
        1.0e-12*scale;

    resid =
        subdiag_norm(
            n,
            a);

    for(iter=0;iter<maxiter;iter++)
    {
        /*
         * Deflate all sufficiently small subdiagonal entries.
         */
        for(i=1;i<n;i++)
        {
            LIS_REAL local_scale =
                1.0
                +zabs1(a[(i-1)+(i-1)*n])
                +zabs1(a[i+i*n]);

            if(
                zabs1(a[i+(i-1)*n])
                <=
                1.0e-12*local_scale )
            {
                a[i+(i-1)*n] = 0.0;
            }
        }

        resid =
            subdiag_norm(
                n,
                a);

        if( resid<=tol )
            break;

        mu =
            wilkinson_shift_complex(
                n,
                a);

        memcpy(
            work,
            a,
            (size_t)n*n*sizeof(*work));

        for(i=0;i<n;i++)
        {
            work[i+i*n] -=
                mu;
        }

        err =
            complex_householder_qr(
                n,
                work,
                q,
                r);

        if( err )
        {
            free(work);
            free(q);
            free(r);
            free(rq);
            free(znew);

            return err;
        }

        matmul(
            n,
            r,
            q,
            rq);

        for(j=0;j<n;j++)
        {
            for(i=0;i<n;i++)
            {
                a[i+j*n] =
                    rq[i+j*n];
            }
        }

        for(i=0;i<n;i++)
        {
            a[i+i*n] +=
                mu;
        }

        matmul(
            n,
            z,
            q,
            znew);

        memcpy(
            z,
            znew,
            (size_t)n*n*sizeof(*z));
    }

    if( niter )
        *niter = iter;

    if( final_resid )
        *final_resid = resid;

    free(work);
    free(q);
    free(r);
    free(rq);
    free(znew);

    if( resid>tol )
        return LIS_MAXITER;

    return LIS_SUCCESS;
}



static LIS_REAL
matrix_frobenius_norm(
    LIS_INT n,
    const LIS_SCALAR *a)
{
    LIS_INT i;
    LIS_REAL sum = 0.0;

    for(i=0;i<n*n;i++)
    {
        LIS_REAL v =
            zabs1(a[i]);

        sum += v*v;
    }

    return sqrt(sum);
}


static LIS_REAL
vector_norm2_local(
    LIS_INT n,
    const LIS_SCALAR *x)
{
    LIS_INT i;
    LIS_REAL sum = 0.0;

    for(i=0;i<n;i++)
    {
        LIS_REAL v =
            zabs1(x[i]);

        sum += v*v;
    }

    return sqrt(sum);
}


/*
 * Compute one right eigenvector of an upper triangular
 * complex Schur matrix.
 *
 * For lambda = T[k,k], choose y[k] = 1 and solve
 *
 *     (T - lambda I)y = 0
 *
 * upward through rows 0..k-1.
 */
static LIS_INT
triangular_right_eigenvector(
    LIS_INT n,
    const LIS_SCALAR *t,
    LIS_INT k,
    LIS_SCALAR *y)
{
    LIS_INT i,j;

    LIS_SCALAR lambda;

    if(
        n<=0 ||
        t==NULL ||
        y==NULL ||
        k<0 ||
        k>=n )
    {
        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<n;i++)
        y[i] = 0.0;

    lambda =
        t[k+k*n];

    y[k] = 1.0;

    for(i=k-1;i>=0;i--)
    {
        LIS_SCALAR sum = 0.0;
        LIS_SCALAR denom;

        LIS_REAL amag;
        LIS_REAL floor;

        for(j=i+1;j<=k;j++)
        {
            sum +=
                t[i+j*n]
                *y[j];
        }

        denom =
            t[i+i*n]
            -lambda;

        amag =
            zabs1(denom);

        floor =
            1.0e-14*
            (
                1.0
                +zabs1(t[i+i*n])
                +zabs1(lambda)
            );

        /*
         * Avoid an uncontrolled division by a nearly zero
         * diagonal difference.  The final eigenpair residual
         * decides whether the regularized vector is acceptable.
         */
        if( amag<floor )
        {
            if( amag>0.0 )
            {
                denom *=
                    floor/amag;
            }
            else
            {
                denom = floor;
            }
        }

        y[i] =
            -sum/denom;

        /*
         * Local overflow protection for highly non-normal
         * triangular systems.
         */
        if( zabs1(y[i])>1.0e100 )
        {
            for(j=i;j<=k;j++)
                y[j] *= 1.0e-100;
        }
    }

    {
        LIS_REAL norm =
            vector_norm2_local(
                n,
                y);

        if( norm<=0.0 )
            return LIS_BREAKDOWN;

        for(i=0;i<n;i++)
            y[i] /= norm;
    }

    return LIS_SUCCESS;
}


static LIS_REAL
right_eigenpair_residual(
    LIS_INT n,
    const LIS_SCALAR *a,
    const LIS_SCALAR *v,
    LIS_SCALAR lambda)
{
    LIS_INT i,j;

    LIS_REAL nr = 0.0;
    LIS_REAL nv;
    LIS_REAL na;
    LIS_REAL denom;

    for(i=0;i<n;i++)
    {
        LIS_SCALAR r =
            -lambda*v[i];

        for(j=0;j<n;j++)
        {
            r +=
                a[i+j*n]
                *v[j];
        }

        {
            LIS_REAL x =
                zabs1(r);

            nr += x*x;
        }
    }

    nr = sqrt(nr);

    nv =
        vector_norm2_local(
            n,
            v);

    na =
        matrix_frobenius_norm(
            n,
            a);

    denom =
        (
            na
            +zabs1(lambda)
        )*nv;

    if( denom<1.0e-30 )
        denom = 1.0;

    return nr/denom;
}


static LIS_INT
schur_right_eigensystem(
    LIS_INT n,
    const LIS_SCALAR *a0,
    const LIS_SCALAR *t,
    const LIS_SCALAR *q,
    LIS_SCALAR *eval,
    LIS_SCALAR *evec,
    LIS_REAL *max_resid)
{
    LIS_INT i,j,k;
    LIS_INT err;

    LIS_REAL worst = 0.0;

    LIS_SCALAR *y;
    LIS_SCALAR *v;

    y =
        (LIS_SCALAR *)malloc(
            (size_t)n*sizeof(*y));

    v =
        (LIS_SCALAR *)malloc(
            (size_t)n*sizeof(*v));

    if( y==NULL || v==NULL )
    {
        free(y);
        free(v);
        return LIS_OUT_OF_MEMORY;
    }

    for(k=0;k<n;k++)
    {
        LIS_REAL resid;
        LIS_REAL norm;

        eval[k] =
            t[k+k*n];

        err =
            triangular_right_eigenvector(
                n,
                t,
                k,
                y);

        if( err )
        {
            free(y);
            free(v);
            return err;
        }

        /*
         * v = Q y
         */
        for(i=0;i<n;i++)
        {
            LIS_SCALAR sum = 0.0;

            for(j=0;j<n;j++)
            {
                sum +=
                    q[i+j*n]
                    *y[j];
            }

            v[i] = sum;
        }

        norm =
            vector_norm2_local(
                n,
                v);

        if( norm<=0.0 )
        {
            free(y);
            free(v);
            return LIS_BREAKDOWN;
        }

        for(i=0;i<n;i++)
        {
            v[i] /= norm;
            evec[i+k*n] = v[i];
        }

        resid =
            right_eigenpair_residual(
                n,
                a0,
                v,
                eval[k]);

        if( resid>worst )
            worst = resid;

    }

    free(y);
    free(v);

    if( max_resid )
        *max_resid = worst;

    return LIS_SUCCESS;
}


#endif


#undef __FUNC__
#define __FUNC__ "lis_eigen_general_complex_solve"
LIS_INT
lis_eigen_general_complex_solve(
    LIS_INT n,
    LIS_SCALAR *a,
    LIS_SCALAR *eval,
    LIS_SCALAR *evec,
    LIS_INT *niter,
    LIS_REAL *error)
{
#ifndef _COMPLEX

    LIS_DEBUG_FUNC_IN;

    if( n<=0 ||
        a==NULL ||
        eval==NULL ||
        evec==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid general complex eigensystem arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    if( niter )
        *niter = 0;

    if( error )
        *error = 0.0;

    LIS_SETERR(
        LIS_ERR_NOT_IMPLEMENTED,
        "general complex dense eigensystem requires complex arithmetic\n");

    LIS_DEBUG_FUNC_OUT;

    return LIS_ERR_NOT_IMPLEMENTED;

#else

    LIS_INT i;
    LIS_INT err;
    LIS_INT iter;

    LIS_REAL lower;
    LIS_REAL residual;

    LIS_SCALAR *original;
    LIS_SCALAR *q;

    LIS_DEBUG_FUNC_IN;

    if( n<=0 ||
        a==NULL ||
        eval==NULL ||
        evec==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid general complex eigensystem arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    iter = 0;
    lower = 0.0;
    residual = 0.0;

    if( niter )
        *niter = 0;

    if( error )
        *error = 0.0;

    original =
        (LIS_SCALAR *)lis_malloc(
            n*n*sizeof(LIS_SCALAR),
            "lis_eigen_general_complex_solve::original");

    q =
        (LIS_SCALAR *)lis_malloc(
            n*n*sizeof(LIS_SCALAR),
            "lis_eigen_general_complex_solve::q");

    if( original==NULL ||
        q==NULL )
    {
        if( q )
            lis_free(q);

        if( original )
            lis_free(original);

        LIS_SETERR_MEM(
            n*n*sizeof(LIS_SCALAR));

        return LIS_OUT_OF_MEMORY;
    }

    for(i=0;i<n*n;i++)
        original[i] = a[i];

    err =
        complex_shifted_qr_schur(
            n,
            a,
            q,
            20000,
            &iter,
            &lower);

    if( err )
    {
        lis_free(q);
        lis_free(original);

        if( niter )
            *niter = iter;

        if( error )
            *error = lower;

        return err;
    }

    err =
        schur_right_eigensystem(
            n,
            original,
            a,
            q,
            eval,
            evec,
            &residual);

    lis_free(q);
    lis_free(original);

    if( niter )
        *niter = iter;

    if( error )
        *error =
            residual>lower
                ? residual
                : lower;

    if( err )
        return err;

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;

#endif
}
