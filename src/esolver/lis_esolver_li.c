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
#include <math.h>
#include <string.h>
#include <stdarg.h>
#ifdef USE_SSE2
	#include <emmintrin.h>
#endif
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/***************************************
 * Lanczos Iteration                   *
 ***************************************
 v(0)    = (0,...,0)^T
 ***************************************
 r = v(0)
 beta(0) = ||r||_2 
 for j=1,2,...
   v(j)      = r / beta(j-1)
   r         = A * v(j)
   r         = r - beta(j-1) * v(j-1)
   alpha(j)  = <v(j), r>
   r         = r - alpha(j) * v(j)
   reorthogonalization 
   beta(j)   = ||r||_2
 end for
 compute eigenvalues of a symmetric tridiagonal matrix 
 T(j) = ST'(j)S^*, where
     (alpha(1) beta(1)                      )
     (beta(1)  alpha(2)                     )
 T = (               ...                    )
     (                  alpha(j-1) beta(j-1))                       
     (                  beta(j-1)  alpha(j) )
 compute refined eigenpairs
 ***************************************/

#define NWORK 2
#undef __FUNC__
#define __FUNC__ "lis_eli_check_params"
LIS_INT lis_eli_check_params(LIS_ESOLVER esolver)
{
        LIS_INT ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
	if( ss<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%D) is less than 0\n",ss);
		return LIS_ERR_ILL_ARG;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_eli_malloc_work"
LIS_INT lis_eli_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_eli_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( esolver->eprecision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(esolver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
			if( err ) break;
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	esolver->worklen = worklen;
	esolver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*
 * Orthogonalize one Lanczos residual against an active basis span.
 * Keeping the active range explicit separates basis length from the
 * requested Ritz count and allows the operation to be reused by
 * multi-vector recurrences.
 */
static LIS_INT
lis_lanczos_reorthogonalize_residual(
    LIS_VECTOR *basis,
    LIS_INT first,
    LIS_INT last,
    LIS_VECTOR residual)
{
    LIS_INT err,k;
    LIS_SCALAR dot;

    for(k=first;k<=last;k++)
    {
        err = lis_vector_dot(
            basis[k],
            residual,
            &dot);
        if( err )
            return err;

        err = lis_vector_axpy(
            -dot,
            basis[k],
            residual);
        if( err )
            return err;
    }

    return LIS_SUCCESS;
}


/*
 * Orthogonalize the generalized Lanczos dual residual.
 *
 * With dual_basis[k] = B*basis[k], subtracting
 * <basis[k],residual>*dual_basis[k] makes B^-1*residual
 * orthogonal to the active basis in the B inner product.
 */
static LIS_INT
lis_lanczos_reorthogonalize_dual_residual(
    LIS_VECTOR *basis,
    LIS_VECTOR *dual_basis,
    LIS_INT first,
    LIS_INT last,
    LIS_VECTOR residual)
{
    LIS_INT err,k;
    LIS_SCALAR dot;

    for(k=first;k<=last;k++)
    {
        err = lis_vector_dot(
            basis[k],
            residual,
            &dot);
        if( err )
            return err;

        err = lis_vector_axpy(
            -dot,
            dual_basis[k],
            residual);
        if( err )
            return err;
    }

    return LIS_SUCCESS;
}


/*
 * Build the logical block layout over the existing linear Lanczos
 * basis storage.
 *
 * The subspace size remains the maximum total basis dimension.
 * A final block may therefore have fewer columns than block_size.
 *
 * Rank reduction caused by numerical deflation is intentionally not
 * handled here.  This foundation records only the nominal layout;
 * rank-revealing orthogonalization is supplied by the shared block
 * eigensolver layer.
 */
static LIS_INT
lis_lanczos_block_layout(
    LIS_INT ss,
    LIS_INT block_size,
    LIS_INT *block_count,
    LIS_INT **block_offset,
    LIS_INT **block_rank)
{
    LIS_INT i;
    LIS_INT count;
    LIS_INT offset;
    LIS_INT rank;
    LIS_INT *offsets;
    LIS_INT *ranks;

    if( ss<1 || block_size<1 ||
        block_count==NULL ||
        block_offset==NULL ||
        block_rank==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid Block Lanczos layout arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if( block_size>ss )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "Block Lanczos block size exceeds subspace size\n");
        return LIS_ERR_ILL_ARG;
    }

    count = (ss + block_size - 1) / block_size;

    offsets = (LIS_INT *)lis_malloc(
        ss*sizeof(LIS_INT),
        "lis_lanczos_block_layout::offsets");

    if( offsets==NULL )
    {
        LIS_SETERR_MEM(
            ss*sizeof(LIS_INT));
        return LIS_ERR_OUT_OF_MEMORY;
    }

    ranks = (LIS_INT *)lis_malloc(
        ss*sizeof(LIS_INT),
        "lis_lanczos_block_layout::ranks");

    if( ranks==NULL )
    {
        lis_free(offsets);

        LIS_SETERR_MEM(
            ss*sizeof(LIS_INT));
        return LIS_ERR_OUT_OF_MEMORY;
    }

    /*
     * Numerical rank deflation can create more blocks than the
     * nominal ceil(ss/block_size) partition.  Reserve metadata for
     * the maximum possible case: ss rank-one blocks.
     */
    for(i=0;i<ss;i++)
    {
        offsets[i] = -1;
        ranks[i] = 0;
    }

    offset = 0;

    for(i=0;i<count;i++)
    {
        rank = block_size;

        if( offset + rank > ss )
        {
            rank = ss - offset;
        }

        offsets[i] = offset;
        ranks[i] = rank;

        offset += rank;
    }

    if( offset!=ss )
    {
        lis_free(ranks);
        lis_free(offsets);

        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid Block Lanczos basis layout\n");
        return LIS_ERR_ILL_ARG;
    }

    *block_count = count;
    *block_offset = offsets;
    *block_rank = ranks;

    return LIS_SUCCESS;
}


/*
 * Build the deterministic raw initial block used by Block Lanczos.
 *
 * Global row g belongs to column j when
 *
 *     g mod block_size == j.
 *
 * The raw block is intentionally not normalized here.  Numerical
 * rank detection, deflation, compaction, and normalization belong
 * to the shared block eigensolver layer.
 */
#undef __FUNC__
#define __FUNC__ "lis_lanczos_block_initial_raw"
LIS_INT
lis_lanczos_block_initial_raw(
    LIS_INT block_size,
    LIS_INT block_rank,
    LIS_VECTOR W[])
{
    LIS_INT err;
    LIS_INT i,j;
    LIS_INT is,ie;
    LIS_INT index;

    LIS_DEBUG_FUNC_IN;

    if( block_size<1 ||
        block_rank<1 ||
        block_rank>block_size ||
        W==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid Block Lanczos initial block arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    for(j=0;j<block_rank;j++)
    {
        if( W[j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null Block Lanczos initial block vector\n");
            return LIS_ERR_ILL_ARG;
        }

        err = lis_vector_set_all(
            0.0,
            W[j]);

        if( err )
        {
            return err;
        }

        err = lis_vector_get_range(
            W[j],
            &is,
            &ie);

        if( err )
        {
            return err;
        }

        for(i=is;i<ie;i++)
        {
            if( (i % block_size)==j )
            {
                index = i + W[j]->origin;

                err = lis_vector_set_value(
                    LIS_INS_VALUE,
                    index,
                    1.0,
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



/*
 * Insert one diagonal Alpha block into the dense projected matrix.
 *
 * The projected matrix uses the full Lanczos subspace capacity as
 * its leading dimension.  Active block offsets and ranks determine
 * the actual occupied entries.
 */
#undef __FUNC__
#define __FUNC__ "lis_lanczos_block_projected_set_alpha"
LIS_INT
lis_lanczos_block_projected_set_alpha(
    LIS_INT ld,
    LIS_SCALAR *t,
    LIS_INT offset,
    LIS_INT rank,
    const LIS_SCALAR *alpha)
{
    LIS_INT i,j;

    LIS_DEBUG_FUNC_IN;

    if( ld<1 ||
        t==NULL ||
        alpha==NULL ||
        offset<0 ||
        rank<1 ||
        offset+rank>ld )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid Block Lanczos Alpha placement arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<rank;i++)
    {
        for(j=0;j<rank;j++)
        {
            t[
                (offset+i)*ld+
                (offset+j)
            ] =
                alpha[i*rank+j];
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


/*
 * Insert one rectangular Beta block and its Hermitian transpose
 * into the dense projected matrix.
 *
 * beta has row-major dimensions:
 *
 *     next_rank x previous_rank
 *
 * and represents the lower off-diagonal block
 *
 *     T(next, previous) = Beta.
 *
 * The matching upper block is Beta^H.
 */
#undef __FUNC__
#define __FUNC__ "lis_lanczos_block_projected_set_beta"
LIS_INT
lis_lanczos_block_projected_set_beta(
    LIS_INT ld,
    LIS_SCALAR *t,
    LIS_INT previous_offset,
    LIS_INT previous_rank,
    LIS_INT next_offset,
    LIS_INT next_rank,
    const LIS_SCALAR *beta)
{
    LIS_INT i,j;
    LIS_SCALAR value;

    LIS_DEBUG_FUNC_IN;

    if( ld<1 ||
        t==NULL ||
        beta==NULL ||
        previous_offset<0 ||
        next_offset<0 ||
        previous_rank<1 ||
        next_rank<1 ||
        previous_offset+previous_rank>ld ||
        next_offset+next_rank>ld )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid Block Lanczos Beta placement arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    for(i=0;i<next_rank;i++)
    {
        for(j=0;j<previous_rank;j++)
        {
            value =
                beta[
                    i*previous_rank+j
                ];

            t[
                (next_offset+i)*ld+
                (previous_offset+j)
            ] = value;

#ifdef _COMPLEX
            t[
                (previous_offset+j)*ld+
                (next_offset+i)
            ] = conj(value);
#else
            t[
                (previous_offset+j)*ld+
                (next_offset+i)
            ] = value;
#endif
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return LIS_SUCCESS;
}


/*
 * Initial Block Lanczos driver.
 *
 * BL0-BL2 establish:
 *   - block-size semantics,
 *   - scalar/block dispatch,
 *   - linear basis capacity,
 *   - nominal block offsets and ranks.
 *
 * The block recurrence is deliberately not approximated here.
 * BL3/BL4 will connect Alpha/Beta block recurrence after the shared
 * block orthogonalization and rank-deflation contract is available.
 */
/*
 * Apply the scalar Lanczos happy-breakdown tolerance to a block
 * residual.
 *
 * Scalar Lanczos stops when beta = ||r||_2 < tol.
 * For a block residual, require every active residual column to
 * satisfy the same tolerance:
 *
 *     max_j ||W_j||_2 < tol.
 *
 * When this occurs, snap the round-off residual block to exact
 * zero. The existing common rank-deflation path then reports
 * active_rank == 0 and terminates without inventing directions
 * from relative rank tests on a near-zero Gram matrix.
 */
static LIS_INT
lis_lanczos_block_apply_breakdown_tol(
    LIS_INT width,
    LIS_VECTOR W[],
    LIS_REAL tol,
    LIS_REAL *max_norm)
{
    LIS_INT err;
    LIS_INT i;
    LIS_REAL nrm;
    LIS_REAL maximum;

    if( width<1 ||
        W==NULL ||
        max_norm==NULL ||
        !isfinite((double)tol) ||
        tol<0.0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid Block Lanczos breakdown tolerance arguments\n");

        return LIS_ERR_ILL_ARG;
    }

    maximum = 0.0;

    for(i=0;i<width;i++)
    {
        if( W[i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null Block Lanczos residual vector\n");

            return LIS_ERR_ILL_ARG;
        }

        err = lis_vector_nrm2(
            W[i],
            &nrm);

        if( err )
        {
            return err;
        }

        if( nrm>maximum )
        {
            maximum = nrm;
        }
    }

    *max_norm = maximum;

    if( maximum<tol )
    {
        for(i=0;i<width;i++)
        {
            err = lis_vector_set_all(
                0.0,
                W[i]);

            if( err )
            {
                return err;
            }
        }
    }

    return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_eli_block"
static LIS_INT
lis_eli_block(LIS_ESOLVER esolver)
{
    LIS_INT err;
    LIS_INT final_status;
    LIS_REAL cert_resid,cert_tol;
    LIS_INT i,j;
    LIS_INT ss;
    LIS_INT block_size;
    LIS_INT block_count;
    LIS_INT active_block_count;
    LIS_INT basis_count;
    LIS_INT nominal_rank;
    LIS_INT active_rank;
    LIS_INT current_block;
    LIS_INT previous_offset;
    LIS_INT previous_rank;
    LIS_INT current_offset;
    LIS_INT current_rank;
    LIS_INT full_next_rank;
    LIS_INT next_rank;
    LIS_INT next_offset;
    LIS_INT remaining;
    LIS_REAL rank_tol;
    LIS_REAL tol;
    LIS_REAL residual_max;
    LIS_INT *block_offset;
    LIS_INT *block_rank;
    LIS_SCALAR *t;
    LIS_SCALAR *coeff;
    LIS_SCALAR *gram_work;
    LIS_SCALAR *small;
    LIS_SCALAR *eval;
    LIS_VECTOR *candidate;
    LIS_VECTOR *partial_q;
    LIS_VECTOR *next_output;

    ss = esolver->options[
        LIS_EOPTIONS_SUBSPACE];

    block_size = esolver->block_size;

    final_status = LIS_SUCCESS;
    cert_tol = esolver->params[
        LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN];

    block_count = 0;
    active_block_count = 0;
    basis_count = 0;
    nominal_rank = 0;
    active_rank = 0;
    current_block = 0;
    previous_offset = 0;
    previous_rank = 0;
    current_offset = 0;
    current_rank = 0;
    full_next_rank = 0;
    next_rank = 0;
    next_offset = 0;
    remaining = 0;
    rank_tol = (LIS_REAL)1.0e-12;

    tol = esolver->params[
        LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN];

    residual_max = 0.0;
    block_offset = NULL;
    block_rank = NULL;
    t = NULL;
    coeff = NULL;
    gram_work = NULL;
    small = NULL;
    eval = NULL;
    candidate = NULL;
    partial_q = NULL;
    next_output = NULL;

    err = lis_lanczos_block_layout(
        ss,
        block_size,
        &block_count,
        &block_offset,
        &block_rank);

    if( err )
        return err;

    /*
     * BL2 intentionally reuses the existing linear basis capacity.
     * No additional vector workspace is needed merely to describe
     * the block partition.
     */
    if( esolver->work==NULL ||
        esolver->worklen < NWORK + ss )
    {
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "insufficient Block Lanczos workspace\n");
        return LIS_ERR_ILL_ARG;
    }

    /*
     * Validate the nominal storage layout before numerical rank
     * deflation.  The actual block ranks and offsets become dynamic
     * once C4 constructs the first active basis block.
     */
    for(i=0;i<block_count;i++)
    {
        if( block_offset[i]!=basis_count ||
            block_rank[i]<1 ||
            block_rank[i]>block_size ||
            block_offset[i]+block_rank[i]>ss )
        {
            lis_free(block_rank);
            lis_free(block_offset);

            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "invalid Block Lanczos block layout\n");
            return LIS_ERR_ILL_ARG;
        }

        basis_count += block_rank[i];
    }

    if( basis_count!=ss )
    {
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "Block Lanczos basis layout does not cover subspace\n");
        return LIS_ERR_ILL_ARG;
    }

    nominal_rank = block_rank[0];

    /*
     * BL4-B1 gives the block path ownership of one dense projected
     * matrix over the full Lanczos subspace capacity.  Alpha and
     * rectangular Beta blocks will be inserted into this matrix by
     * the BL4-B0 placement helpers once the numerical recurrence is
     * connected.
     */
    t = (LIS_SCALAR *)lis_malloc(
        ss*ss*sizeof(LIS_SCALAR),
        "lis_eli_block::t");

    if( t==NULL )
    {
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR_MEM(
            ss*ss*sizeof(LIS_SCALAR));

        return LIS_ERR_OUT_OF_MEMORY;
    }

    for(i=0;i<ss*ss;i++)
    {
        t[i] = 0.0;
    }

    /*
     * BL4-B3 reserves capacities for the future variable-rank
     * recurrence.  Runtime kernels will use actual active ranks;
     * these allocations only define upper bounds.
     *
     * C3 requires coeff and work to be distinct nq-by-nw arrays.
     * Their worst-case Block Lanczos capacity is
     *
     *     ss x block_size.
     *
     * The small square buffer is reserved for current-block
     * quantities such as Alpha and the future C4 normalization
     * factor.
     */
    coeff = (LIS_SCALAR *)lis_malloc(
        ss*block_size*sizeof(LIS_SCALAR),
        "lis_eli_block::coeff");

    if( coeff==NULL )
    {
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR_MEM(
            ss*block_size*sizeof(LIS_SCALAR));

        return LIS_ERR_OUT_OF_MEMORY;
    }

    gram_work = (LIS_SCALAR *)lis_malloc(
        ss*block_size*sizeof(LIS_SCALAR),
        "lis_eli_block::gram_work");

    if( gram_work==NULL )
    {
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR_MEM(
            ss*block_size*sizeof(LIS_SCALAR));

        return LIS_ERR_OUT_OF_MEMORY;
    }

    small = (LIS_SCALAR *)lis_malloc(
        block_size*block_size*sizeof(LIS_SCALAR),
        "lis_eli_block::small");

    if( small==NULL )
    {
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR_MEM(
            block_size*block_size*sizeof(LIS_SCALAR));

        return LIS_ERR_OUT_OF_MEMORY;
    }

    candidate = (LIS_VECTOR *)lis_malloc(
        block_size*sizeof(LIS_VECTOR),
        "lis_eli_block::candidate");

    if( candidate==NULL )
    {
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR_MEM(
            block_size*sizeof(LIS_VECTOR));

        return LIS_ERR_OUT_OF_MEMORY;
    }

    for(i=0;i<block_size;i++)
    {
        candidate[i] = NULL;
    }

    for(i=0;i<block_size;i++)
    {
        if( esolver->eprecision==LIS_PRECISION_DEFAULT )
        {
            err = lis_vector_duplicate(
                esolver->A,
                &candidate[i]);
        }
        else
        {
            err = lis_vector_duplicateex(
                LIS_PRECISION_QUAD,
                esolver->A,
                &candidate[i]);
        }

        if( err )
        {
            LIS_INT j;

            for(j=0;j<i;j++)
            {
                lis_vector_destroy(
                    candidate[j]);
            }

            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }
    }

    /*
     * C4 rank-revealing normalization requires one eigenvalue array
     * in addition to the three dense block-size square workspaces.
     */
    eval = (LIS_SCALAR *)lis_malloc(
        block_size*sizeof(LIS_SCALAR),
        "lis_eli_block::eval");

    if( eval==NULL )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR_MEM(
            block_size*sizeof(LIS_SCALAR));

        return LIS_ERR_OUT_OF_MEMORY;
    }

    /*
     * Build the deterministic raw candidate block out of place.
     * C4 owns numerical rank detection, compaction, and
     * normalization into the persistent Lanczos basis.
     */
    err = lis_lanczos_block_initial_raw(
        block_size,
        nominal_rank,
        candidate);

    if( err )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        return err;
    }

    err = lis_eigen_block_rank_deflate(
        block_size,
        nominal_rank,
        candidate,
        &esolver->work[NWORK],
        rank_tol,
        &active_rank,
        coeff,
        eval,
        gram_work,
        small);

    if( err )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        return err;
    }

    if( active_rank==0 )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        LIS_SETERR(
            LIS_BREAKDOWN,
            "Block Lanczos initial block has zero numerical rank\n");

        esolver->retcode = LIS_BREAKDOWN;

        return LIS_BREAKDOWN;
    }

    /*
     * From this point block_rank/block_offset describe the numerical
     * basis rather than the nominal storage partition.
     */
    block_rank[0] = active_rank;
    block_offset[0] = 0;
    basis_count = active_rank;
    active_block_count = 1;

    /*
     * The metadata arrays have ss-entry numerical capacity.
     * All blocks after Q0 remain unbuilt until recurrence creates
     * them, regardless of the original nominal block count.
     */
    for(i=1;i<ss;i++)
    {
        block_rank[i] = 0;
        block_offset[i] = -1;
    }

    /*
     * First numerical Block Lanczos projection:
     *
     *     Y0     = A Q0
     *     Alpha0 = Q0^H Y0
     *
     * candidate[] is no longer needed for the raw initial block,
     * so it becomes the transient A*Q0 block.  The C4 transform in
     * small[] is also complete, allowing the same dense workspace
     * to hold Alpha0.
     */
    err = lis_matvec_block(
        esolver->A,
        active_rank,
        &esolver->work[NWORK],
        candidate);

    if( err )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        return err;
    }

    err = lis_eigen_block_gram(
        active_rank,
        &esolver->work[NWORK],
        active_rank,
        candidate,
        small);

    if( err )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        return err;
    }

    err =
        lis_lanczos_block_projected_set_alpha(
            ss,
            t,
            block_offset[0],
            active_rank,
            small);

    if( err )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        return err;
    }

    /*
     * First Block Lanczos residual:
     *
     *     W1 = A*Q0 - Q0*Alpha0
     *
     * candidate[] currently owns A*Q0 and small[] owns Alpha0.
     * Alpha0 uses row-major indexing:
     *
     *     Alpha0(i,j) = small[i*active_rank+j]
     *
     * so candidate[j] receives the j-th residual column.
     */
    for(j=0;j<active_rank;j++)
    {
        for(i=0;i<active_rank;i++)
        {
            err = lis_vector_axpy(
                -small[i*active_rank+j],
                esolver->work[NWORK+i],
                candidate[j]);

            if( err )
            {
                for(i=0;i<block_size;i++)
                {
                    lis_vector_destroy(
                        candidate[i]);
                }

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                return err;
            }
        }
    }

    /*
     * Remove the round-off projection left after the Alpha0
     * recurrence update.
     *
     * The explicit Alpha0 subtraction above is the first projection.
     * One C3 pass therefore provides the full reorthogonalization
     * pass used by the scalar Lanczos path.
     *
     * coeff[] records the removed round-off projection for inspection,
     * but Alpha0/T remains the projection computed in BL4-D.
     */
    err = lis_eigen_block_orthogonalize(
        active_rank,
        &esolver->work[NWORK],
        active_rank,
        candidate,
        1,
        coeff,
        gram_work);

    if( err )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        return err;
    }

    /*
     * First complete two-block Lanczos transition.
     *
     * candidate[] contains the orthogonalized W1 from BL4-E.
     *
     * C4 requires:
     *
     *     nominal_rank <= capacity
     *
     * so rank deflation always receives capacity=previous_rank.
     *
     * If the remaining linear basis capacity is smaller than
     * previous_rank, use a temporary full-width output block and
     * retain only the directions that fit in the final partial block.
     */
    /*
     * Scalar-compatible happy breakdown for the first residual
     * block. C4 remains responsible for relative rank reduction
     * whenever the residual block is not absolutely negligible.
     */
    err =
        lis_lanczos_block_apply_breakdown_tol(
            active_rank,
            candidate,
            tol,
            &residual_max);

    if( err )
    {
        for(i=0;i<block_size;i++)
        {
            lis_vector_destroy(
                candidate[i]);
        }

        lis_free(eval);
        lis_free(candidate);
        lis_free(small);
        lis_free(gram_work);
        lis_free(coeff);
        lis_free(t);
        lis_free(block_rank);
        lis_free(block_offset);

        return err;
    }

    previous_rank = block_rank[0];
    next_offset = basis_count;
    remaining = ss - basis_count;

    full_next_rank = 0;
    next_rank = 0;
    partial_q = NULL;
    next_output = NULL;

    if( remaining>0 )
    {
        if( remaining>=previous_rank )
        {
            next_output =
                &esolver->work[NWORK+next_offset];
        }
        else
        {
            partial_q =
                (LIS_VECTOR *)lis_malloc(
                    previous_rank*sizeof(LIS_VECTOR),
                    "lis_eli_block::partial_q");

            if( partial_q==NULL )
            {
                for(i=0;i<block_size;i++)
                {
                    lis_vector_destroy(
                        candidate[i]);
                }

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                LIS_SETERR_MEM(
                    previous_rank*sizeof(LIS_VECTOR));

                return LIS_ERR_OUT_OF_MEMORY;
            }

            for(i=0;i<previous_rank;i++)
            {
                partial_q[i] = NULL;
            }

            for(i=0;i<previous_rank;i++)
            {
                if( esolver->eprecision==
                    LIS_PRECISION_DEFAULT )
                {
                    err = lis_vector_duplicate(
                        esolver->A,
                        &partial_q[i]);
                }
                else
                {
                    err = lis_vector_duplicateex(
                        LIS_PRECISION_QUAD,
                        esolver->A,
                        &partial_q[i]);
                }

                if( err )
                {
                    for(j=0;j<previous_rank;j++)
                    {
                        if( partial_q[j]!=NULL )
                        {
                            lis_vector_destroy(
                                partial_q[j]);
                        }
                    }

                    lis_free(partial_q);

                    for(j=0;j<block_size;j++)
                    {
                        lis_vector_destroy(
                            candidate[j]);
                    }

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }
            }

            next_output = partial_q;
        }

        err = lis_eigen_block_rank_deflate(
            previous_rank,
            previous_rank,
            candidate,
            next_output,
            rank_tol,
            &full_next_rank,
            coeff,
            eval,
            gram_work,
            small);

        if( err )
        {
            if( partial_q!=NULL )
            {
                for(i=0;i<previous_rank;i++)
                {
                    lis_vector_destroy(
                        partial_q[i]);
                }

                lis_free(partial_q);
            }

            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        next_rank = full_next_rank;

        if( next_rank>remaining )
        {
            next_rank = remaining;
        }

        /*
         * For the final partial block the C4 output was intentionally
         * built in temporary full-width storage.  Copy only the
         * retained directions that still fit in the linear basis.
         */
        if( partial_q!=NULL )
        {
            for(i=0;i<next_rank;i++)
            {
                err = lis_vector_copy(
                    partial_q[i],
                    esolver->work[
                        NWORK+next_offset+i]);

                if( err )
                {
                    for(j=0;j<previous_rank;j++)
                    {
                        lis_vector_destroy(
                            partial_q[j]);
                    }

                    lis_free(partial_q);

                    for(j=0;j<block_size;j++)
                    {
                        lis_vector_destroy(
                            candidate[j]);
                    }

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }
            }

            for(i=0;i<previous_rank;i++)
            {
                lis_vector_destroy(
                    partial_q[i]);
            }

            lis_free(partial_q);
            partial_q = NULL;
            next_output = NULL;
        }

        /*
         * active_rank==0 after W1 rank deflation is a valid Lanczos
         * happy breakdown.  In that case there is no Q1/Beta1 block.
         */
        if( next_rank>0 )
        {
            /*
             * C4 stores the coupling as:
             *
             *     Beta1 = Lambda^(1/2) U^H
             *
             * with logical dimensions:
             *
             *     next_rank x previous_rank.
             *
             * small[] no longer needs the rank-transform data, so it
             * is reused for this rectangular coupling.
             */
            err =
                lis_eigen_block_rectangular_coupling(
                    previous_rank,
                    next_rank,
                    eval,
                    gram_work,
                    small);

            if( err )
            {
                for(i=0;i<block_size;i++)
                {
                    lis_vector_destroy(
                        candidate[i]);
                }

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                return err;
            }

            err =
                lis_lanczos_block_projected_set_beta(
                    ss,
                    t,
                    block_offset[0],
                    previous_rank,
                    next_offset,
                    next_rank,
                    small);

            if( err )
            {
                for(i=0;i<block_size;i++)
                {
                    lis_vector_destroy(
                        candidate[i]);
                }

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                return err;
            }

            block_offset[active_block_count] =
                next_offset;

            block_rank[active_block_count] =
                next_rank;

            basis_count += next_rank;
            active_block_count++;
        }
    }

    /*
     * General variable-rank Block Lanczos recurrence.
     *
     * BL4-F remains the validated bootstrap:
     *
     *     Q0 -> Alpha0 -> W1 -> Q1 -> Beta1
     *
     * Continue from numerical block Q1.  For j >= 1:
     *
     *     W = A*Qj - Q(j-1)*Beta_j^H
     *     Alpha_j = Qj^H*W
     *     W -= Qj*Alpha_j
     *     full reorthogonalization
     *     W -> Q(j+1), Beta_(j+1)
     *
     * block_rank[] and block_offset[] describe the numerical,
     * rank-deflated basis.
     */
    current_block = 1;

    while( current_block<active_block_count )
    {
        previous_offset =
            block_offset[current_block-1];

        previous_rank =
            block_rank[current_block-1];

        current_offset =
            block_offset[current_block];

        current_rank =
            block_rank[current_block];

        if( previous_offset<0 ||
            previous_rank<1 ||
            current_offset<0 ||
            current_rank<1 ||
            previous_offset+previous_rank>basis_count ||
            current_offset+current_rank>basis_count )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "invalid numerical Block Lanczos recurrence layout\n");

            return LIS_ERR_ILL_ARG;
        }

        /*
         * Y_j = A * Q_j
         */
        err = lis_matvec_block(
            esolver->A,
            current_rank,
            &esolver->work[
                NWORK+current_offset],
            candidate);

        if( err )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        /*
         * Remove the previous block coupling:
         *
         *     Y_j -= Q_(j-1) * Beta_j^H
         *
         * set_beta() already stored Beta_j^H in the upper
         * projected block T[previous,current].  Reading that
         * block directly therefore preserves the correct
         * conjugation in complex arithmetic.
         */
        for(j=0;j<current_rank;j++)
        {
            for(i=0;i<previous_rank;i++)
            {
                err = lis_vector_axpy(
                    -t[
                        (previous_offset+i)*ss+
                        (current_offset+j)
                    ],
                    esolver->work[
                        NWORK+previous_offset+i],
                    candidate[j]);

                if( err )
                {
                    for(i=0;i<block_size;i++)
                    {
                        lis_vector_destroy(
                            candidate[i]);
                    }

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }
            }
        }

        /*
         * Alpha_j = Q_j^H * Y_j
         */
        err = lis_eigen_block_gram(
            current_rank,
            &esolver->work[
                NWORK+current_offset],
            current_rank,
            candidate,
            small);

        if( err )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        err =
            lis_lanczos_block_projected_set_alpha(
                ss,
                t,
                current_offset,
                current_rank,
                small);

        if( err )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        /*
         * A basis-capacity stop still requires the diagonal block
         * Alpha_j of the final active block. Once Alpha_j is stored,
         * the projected matrix is complete and no residual or
         * Q_(j+1) construction is required.
         */
        if( basis_count>=ss )
        {
            break;
        }

        /*
         * W_(j+1) = Y_j - Q_j*Alpha_j
         */
        for(j=0;j<current_rank;j++)
        {
            for(i=0;i<current_rank;i++)
            {
                err = lis_vector_axpy(
                    -small[
                        i*current_rank+j],
                    esolver->work[
                        NWORK+current_offset+i],
                    candidate[j]);

                if( err )
                {
                    for(i=0;i<block_size;i++)
                    {
                        lis_vector_destroy(
                            candidate[i]);
                    }

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }
            }
        }

        /*
         * The recurrence has explicitly removed the previous
         * coupling and Alpha_j projections.  One C3 pass removes
         * the remaining round-off projection against the entire
         * active basis Q0..Qj.
         *
         * As in the scalar Lanczos path, this correction is not
         * added back to the already stored projected coefficients.
         */
        err = lis_eigen_block_orthogonalize(
            basis_count,
            &esolver->work[NWORK],
            current_rank,
            candidate,
            1,
            coeff,
            gram_work);

        if( err )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        /*
         * Apply the same absolute happy-breakdown criterion after
         * full-basis reorthogonalization. Relative C4 rank
         * deflation is used only when a non-negligible residual
         * block remains.
         */
        err =
            lis_lanczos_block_apply_breakdown_tol(
                current_rank,
                candidate,
                tol,
                &residual_max);

        if( err )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        next_offset = basis_count;
        remaining = ss - basis_count;

        full_next_rank = 0;
        next_rank = 0;

        partial_q = NULL;
        next_output = NULL;

        /*
         * basis_count<ss is part of the loop condition, so
         * remaining is positive here.
         *
         * C4 requires capacity >= nominal_rank.  The nominal
         * residual width is current_rank.  For a final partial
         * block, first build the full rank-deflated block in
         * temporary storage and copy only the directions that fit.
         */
        if( remaining>=current_rank )
        {
            next_output =
                &esolver->work[
                    NWORK+next_offset];
        }
        else
        {
            partial_q =
                (LIS_VECTOR *)lis_malloc(
                    current_rank*sizeof(LIS_VECTOR),
                    "lis_eli_block::general_partial_q");

            if( partial_q==NULL )
            {
                for(i=0;i<block_size;i++)
                {
                    lis_vector_destroy(
                        candidate[i]);
                }

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                LIS_SETERR_MEM(
                    current_rank*sizeof(LIS_VECTOR));

                return LIS_ERR_OUT_OF_MEMORY;
            }

            for(i=0;i<current_rank;i++)
            {
                partial_q[i] = NULL;
            }

            for(i=0;i<current_rank;i++)
            {
                if( esolver->eprecision==
                    LIS_PRECISION_DEFAULT )
                {
                    err = lis_vector_duplicate(
                        esolver->A,
                        &partial_q[i]);
                }
                else
                {
                    err = lis_vector_duplicateex(
                        LIS_PRECISION_QUAD,
                        esolver->A,
                        &partial_q[i]);
                }

                if( err )
                {
                    for(j=0;j<current_rank;j++)
                    {
                        if( partial_q[j]!=NULL )
                        {
                            lis_vector_destroy(
                                partial_q[j]);
                        }
                    }

                    lis_free(partial_q);

                    for(j=0;j<block_size;j++)
                    {
                        lis_vector_destroy(
                            candidate[j]);
                    }

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }
            }

            next_output = partial_q;
        }

        err = lis_eigen_block_rank_deflate(
            current_rank,
            current_rank,
            candidate,
            next_output,
            rank_tol,
            &full_next_rank,
            coeff,
            eval,
            gram_work,
            small);

        if( err )
        {
            if( partial_q!=NULL )
            {
                for(i=0;i<current_rank;i++)
                {
                    lis_vector_destroy(
                        partial_q[i]);
                }

                lis_free(partial_q);
            }

            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        next_rank = full_next_rank;

        if( next_rank>remaining )
        {
            next_rank = remaining;
        }

        if( partial_q!=NULL )
        {
            for(i=0;i<next_rank;i++)
            {
                err = lis_vector_copy(
                    partial_q[i],
                    esolver->work[
                        NWORK+next_offset+i]);

                if( err )
                {
                    for(j=0;j<current_rank;j++)
                    {
                        lis_vector_destroy(
                            partial_q[j]);
                    }

                    lis_free(partial_q);

                    for(j=0;j<block_size;j++)
                    {
                        lis_vector_destroy(
                            candidate[j]);
                    }

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }
            }

            for(i=0;i<current_rank;i++)
            {
                lis_vector_destroy(
                    partial_q[i]);
            }

            lis_free(partial_q);

            partial_q = NULL;
            next_output = NULL;
        }

        /*
         * Zero numerical rank is a valid invariant-subspace
         * termination.  Keep the actual projected dimension equal
         * to basis_count and leave the loop without creating a
         * zero-rank metadata block.
         */
        if( next_rank==0 )
        {
            break;
        }

        /*
         * Beta_(j+1) has logical dimensions:
         *
         *     next_rank x current_rank
         */
        err =
            lis_eigen_block_rectangular_coupling(
                current_rank,
                next_rank,
                eval,
                gram_work,
                small);

        if( err )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        err =
            lis_lanczos_block_projected_set_beta(
                ss,
                t,
                current_offset,
                current_rank,
                next_offset,
                next_rank,
                small);

        if( err )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        if( active_block_count>=ss )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "Block Lanczos numerical block metadata overflow\n");

            return LIS_ERR_ILL_ARG;
        }

        block_offset[active_block_count] =
            next_offset;

        block_rank[active_block_count] =
            next_rank;

        basis_count += next_rank;
        active_block_count++;

        current_block++;
    }

    if( esolver->options[LIS_EOPTIONS_OUTPUT] )
    {
        lis_printf(
            LIS_COMM_WORLD,
            "Block Lanczos block size     : %D\n",
            block_size);

        lis_printf(
            LIS_COMM_WORLD,
            "Block Lanczos block count    : %D\n",
            block_count);

        lis_printf(
            LIS_COMM_WORLD,
            "Block Lanczos active blocks  : %D\n",
            active_block_count);

        lis_printf(
            LIS_COMM_WORLD,
            "Block Lanczos basis capacity : %D\n",
            basis_count);
    }

    /*
     * Solve the active projected eigensystem.
     *
     * T keeps its physical ss-by-ss stride, while basis_count is
     * the numerical projected dimension after rank deflation or
     * happy breakdown.  The common projected solver modifies its
     * matrix argument in place, so copy only the active leading
     * block into compact storage.
     *
     * No Ritz ordering or selection policy is applied here.
     */
    {
        LIS_INT projected_dim;
        LIS_INT projected_niter;
        LIS_INT pi,pj;
        LIS_REAL projected_error;
        LIS_SCALAR *projected_matrix;
        LIS_SCALAR *projected_eval;
        LIS_SCALAR *projected_evec;

        projected_dim = basis_count;
        projected_niter = 0;
        projected_error = 0.0;

        projected_matrix = NULL;
        projected_eval = NULL;
        projected_evec = NULL;

        if( projected_dim<1 ||
            projected_dim>ss )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "invalid Block Lanczos projected dimension\n");

            return LIS_ERR_ILL_ARG;
        }

        /*
         * The active projected dimension is the number of
         * physically available Ritz results. Keep it separate
         * from the nominal -ss capacity so public result getters
         * do not expose inactive slots after rank deflation.
         */
        esolver->nev = projected_dim;

        projected_matrix =
            (LIS_SCALAR *)lis_malloc(
                projected_dim*projected_dim*
                sizeof(LIS_SCALAR),
                "lis_eli_block::projected_matrix");

        if( projected_matrix==NULL )
        {
            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            LIS_SETERR_MEM(
                projected_dim*projected_dim*
                sizeof(LIS_SCALAR));

            return LIS_ERR_OUT_OF_MEMORY;
        }

        projected_eval =
            (LIS_SCALAR *)lis_malloc(
                projected_dim*sizeof(LIS_SCALAR),
                "lis_eli_block::projected_eval");

        if( projected_eval==NULL )
        {
            lis_free(projected_matrix);

            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            LIS_SETERR_MEM(
                projected_dim*sizeof(LIS_SCALAR));

            return LIS_ERR_OUT_OF_MEMORY;
        }

        projected_evec =
            (LIS_SCALAR *)lis_malloc(
                projected_dim*projected_dim*
                sizeof(LIS_SCALAR),
                "lis_eli_block::projected_evec");

        if( projected_evec==NULL )
        {
            lis_free(projected_eval);
            lis_free(projected_matrix);

            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            LIS_SETERR_MEM(
                projected_dim*projected_dim*
                sizeof(LIS_SCALAR));

            return LIS_ERR_OUT_OF_MEMORY;
        }

        for(pi=0;pi<projected_dim;pi++)
        {
            for(pj=0;pj<projected_dim;pj++)
            {
                projected_matrix[
                    pi*projected_dim+pj
                ] =
                    t[
                        pi*ss+pj
                    ];
            }
        }

        err = lis_eigen_projected_solve(
            projected_dim,
            projected_matrix,
            projected_eval,
            projected_evec,
            &projected_niter,
            &projected_error);

        if( err )
        {
            lis_free(projected_evec);
            lis_free(projected_eval);
            lis_free(projected_matrix);

            for(i=0;i<block_size;i++)
            {
                lis_vector_destroy(
                    candidate[i]);
            }

            lis_free(eval);
            lis_free(candidate);
            lis_free(small);
            lis_free(gram_work);
            lis_free(coeff);
            lis_free(t);
            lis_free(block_rank);
            lis_free(block_offset);

            return err;
        }

        /*
         * Preserve the unsorted common projected eigensystem.
         * BL5 does not impose a Ritz ordering policy.
         */
        for(pi=0;pi<projected_dim;pi++)
        {
            esolver->evalue[pi] =
                projected_eval[pi];
        }

        if( esolver->options[LIS_EOPTIONS_OUTPUT] )
        {
            lis_printf(
                LIS_COMM_WORLD,
                "Block Lanczos projected dimension : %D\n",
                projected_dim);

            lis_printf(
                LIS_COMM_WORLD,
                "projected eigensolver iterations  : %D\n",
                projected_niter);

            lis_printf(
                LIS_COMM_WORLD,
                "projected eigensolver error       : %e\n",
                (double)projected_error);
        }

        /*
         * BL7 reconstruct physical Ritz vectors only when the
         * caller requested eigenvectors. The common helper uses
         * the active Lanczos basis directly.
         */
        if( !esolver->options[LIS_EOPTIONS_RVAL] )
        {
            if( esolver->evector==NULL )
            {
                lis_free(projected_evec);
                lis_free(projected_eval);
                lis_free(projected_matrix);

                for(i=0;i<block_size;i++)
                    lis_vector_destroy(candidate[i]);

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "Block Lanczos eigenvector storage is unavailable\n");

                return LIS_ERR_ILL_ARG;
            }

            for(pi=0;pi<projected_dim;pi++)
            {
                if( esolver->evector[pi]!=NULL )
                {
                    for(pj=0;pj<pi;pj++)
                    {
                        lis_vector_destroy(esolver->evector[pj]);
                        esolver->evector[pj] = NULL;
                    }

                    lis_free(projected_evec);
                    lis_free(projected_eval);
                    lis_free(projected_matrix);

                    for(i=0;i<block_size;i++)
                        lis_vector_destroy(candidate[i]);

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Block Lanczos eigenvector slot is already in use\n");

                    return LIS_ERR_ILL_ARG;
                }

                if( esolver->eprecision==LIS_PRECISION_DEFAULT )
                {
                    err = lis_vector_duplicate(
                        esolver->A,
                        &esolver->evector[pi]);
                }
                else
                {
                    err = lis_vector_duplicateex(
                        LIS_PRECISION_QUAD,
                        esolver->A,
                        &esolver->evector[pi]);
                }

                if( err )
                {
                    for(pj=0;pj<pi;pj++)
                    {
                        lis_vector_destroy(esolver->evector[pj]);
                        esolver->evector[pj] = NULL;
                    }

                    lis_free(projected_evec);
                    lis_free(projected_eval);
                    lis_free(projected_matrix);

                    for(i=0;i<block_size;i++)
                        lis_vector_destroy(candidate[i]);

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }
            }

            err = lis_eigen_ritz_reconstruct(
                projected_dim,
                &esolver->work[NWORK],
                projected_dim,
                projected_evec,
                projected_dim,
                esolver->evector);

            if( err )
            {
                for(pi=0;pi<projected_dim;pi++)
                {
                    lis_vector_destroy(esolver->evector[pi]);
                    esolver->evector[pi] = NULL;
                }

                lis_free(projected_evec);
                lis_free(projected_eval);
                lis_free(projected_matrix);

                for(i=0;i<block_size;i++)
                    lis_vector_destroy(candidate[i]);

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                return err;
            }

            err = lis_vector_copy(
                esolver->evector[0],
                esolver->x);

            if( err )
            {
                for(pi=0;pi<projected_dim;pi++)
                {
                    lis_vector_destroy(esolver->evector[pi]);
                    esolver->evector[pi] = NULL;
                }

                lis_free(projected_evec);
                lis_free(projected_eval);
                lis_free(projected_matrix);

                for(i=0;i<block_size;i++)
                    lis_vector_destroy(candidate[i]);

                lis_free(eval);
                lis_free(candidate);
                lis_free(small);
                lis_free(gram_work);
                lis_free(coeff);
                lis_free(t);
                lis_free(block_rank);
                lis_free(block_offset);

                return err;
            }

            /*
             * BL8 certify reconstructed Ritz pairs with the same
             * normalized standard residual used by scalar Lanczos.
             * candidate[0] and candidate[1] are no longer needed by
             * the completed recurrence and provide two live scratch
             * vectors because Block Lanczos requires block_size > 1.
             */
            for(pi=0;pi<projected_dim;pi++)
            {
                err = lis_esolver_standard_residual(
                    esolver->A,
                    esolver->evector[pi],
                    esolver->evalue[pi],
                    candidate[0],
                    candidate[1],
                    &cert_resid);

                if( err )
                {
                    for(pj=0;pj<projected_dim;pj++)
                    {
                        lis_vector_destroy(esolver->evector[pj]);
                        esolver->evector[pj] = NULL;
                    }

                    lis_free(projected_evec);
                    lis_free(projected_eval);
                    lis_free(projected_matrix);

                    for(i=0;i<block_size;i++)
                        lis_vector_destroy(candidate[i]);

                    lis_free(eval);
                    lis_free(candidate);
                    lis_free(small);
                    lis_free(gram_work);
                    lis_free(coeff);
                    lis_free(t);
                    lis_free(block_rank);
                    lis_free(block_offset);

                    return err;
                }

                esolver->resid[pi] = cert_resid;

                if( cert_resid>cert_tol &&
                    final_status==LIS_SUCCESS )
                {
                    final_status = LIS_MAXITER;
                }
            }

        }

        lis_free(projected_evec);
        lis_free(projected_eval);
        lis_free(projected_matrix);
        }

for(i=0;i<block_size;i++)
    {
        lis_vector_destroy(
            candidate[i]);
    }

    lis_free(eval);
    lis_free(candidate);
    lis_free(small);
    lis_free(gram_work);
    lis_free(coeff);
    lis_free(t);
    lis_free(block_rank);
    lis_free(block_offset);

    /*
     * Match scalar Lanczos raw-Ritz behavior.
     * Preserve the common projected eigensolver ordering exactly:
     * no explicit sort, reorder, or -m selection is applied here.
     * With -rval true, raw Ritz values are the completed result.
     */
    if( esolver->options[LIS_EOPTIONS_RVAL] )
    {
        return LIS_SUCCESS;
    }

/* Physical Ritz vectors are available for -rval false. */
esolver->retcode = final_status;
    return final_status;
}


#undef __FUNC__
#define __FUNC__ "lis_eli"
LIS_INT lis_eli(LIS_ESOLVER esolver)
{
  LIS_Comm comm; 
  LIS_MATRIX A;
  LIS_INT err;
  LIS_INT ss,ic;
  LIS_INT emaxiter,iter0,qriter;
  LIS_REAL tol,qrerr;
  double time,time0;    
  LIS_INT i,j,k;
  LIS_INT output, niesolver;
  LIS_REAL nrm2,resid0,cert_resid;
  LIS_SCALAR dot;
  LIS_VECTOR *v,r;
  LIS_SCALAR *t,*tq,*tr,evalue,evalue0;
  LIS_SOLVER solver;
  LIS_ESOLVER esolver2;
  char esolvername[128],solvername[128],preconname[128];
  LIS_INT nsol,precon_type;
  LIS_INT rval,final_status;

  LIS_DEBUG_FUNC_IN;

  if( esolver->block_size>1 )
    {
      err = lis_eli_block(esolver);
      LIS_DEBUG_FUNC_OUT;
      return err;
    }

  comm = LIS_COMM_WORLD;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];
  rval = esolver->options[LIS_EOPTIONS_RVAL];  

  t = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::t");
  tq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tq");
  tr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tr");
  
  A = esolver->A;
  r = esolver->work[0];
  v = &esolver->work[1];
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);

  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  err = lis_solver_set_optionC(solver);
  CHKERR(err);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  lis_esolver_get_esolvername(niesolver, esolvername);
  if( output )
    {
      lis_printf(comm,"inner eigensolver     : %s\n", esolvername);
      lis_printf(comm,"linear solver         : %s\n", solvername);
      lis_printf(comm,"preconditioner        : %s\n", preconname);
    }

  for (i=0;i<ss*ss;i++) t[i] = 0.0;

  j=0;
  while (j<ss-1)
    {
      j = j+1;
      lis_vector_copy(r, v[j]);

      /* v(j) = r / beta(j-1) */
      /* r = A * v(j) */
      /* r = r - beta(j-1) * v(j-1) */
      if (j==1) {
	lis_vector_scale(1.0/nrm2, v[j]);
	lis_matvec(A, v[j], r);
      }
      else {
	lis_vector_scale(1.0/t[(j-2)*ss+j-1], v[j]);
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-t[(j-2)*ss+j-1], v[j-1], r); 
      }

      /* alpha(j)  = <v(j), r> */
      lis_vector_dot(v[j], r, &t[(j-1)*ss+j-1]); 

      /* r = r - alpha(j) * v(j) */
      lis_vector_axpy(-t[(j-1)*ss+j-1], v[j], r); 

      /*
       * Reorthogonalize the residual that generates v(j+1).
       * Include the current vector as well to remove round-off
       * left after the alpha(j) recurrence update.
       */
      err = lis_lanczos_reorthogonalize_residual(
          v,
          1,
          j,
          r);
      if( err )
          return err;

      /* beta(j) = ||r||_2 */
      lis_vector_nrm2(r, (LIS_REAL *)&t[(j-1)*ss+j]);

      /* convergence check */
      if (fabs(t[(j-1)*ss+j])<tol) break;  
      t[j*ss+j-1] = t[(j-1)*ss+j];
    }

  /* compute eigenvalues of a symmetric tridiagonal matrix 
     T(j) = ST'(j)S^* */
  time0 = lis_wtime();
  lis_array_qr(ss,t,tq,tr,&qriter,&qrerr);
  time = lis_wtime() - time0;

  for (i=0;i<ss;i++)
    {
      esolver->evalue[i] = t[i*ss+i];
    }

  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"Ritz values:\n\n");
      for (i=0;i<ss;i++)
	{
	  lis_printf(comm,"Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Lanczos: Ritz value           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else
	  lis_printf(comm,"Lanczos: Ritz value           = %e\n", (double)esolver->evalue[i]);
#endif	  
	}
      lis_printf(comm,"Lanczos: elapsed time         = %e sec.\n\n", time);
    }
  if( rval )
    {
      lis_free(t); 
      lis_free(tq);
      lis_free(tr);
      lis_solver_destroy(solver);
      return LIS_SUCCESS;
    }
  
  if( output ) 
    {
      lis_printf(comm,"computing refined eigenpairs using inner eigensolver:\n\n");      
    }

  lis_esolver_create(&esolver2);
  esolver2->options[LIS_EOPTIONS_ESOLVER] = niesolver;
  esolver2->options[LIS_EOPTIONS_SUBSPACE] = 1;
  esolver2->options[LIS_EOPTIONS_MAXITER] = emaxiter;
  esolver2->options[LIS_EOPTIONS_OUTPUT] = esolver->options[LIS_EOPTIONS_OUTPUT];
  esolver2->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = tol; 

  /* compute refined eigenpairs */
  final_status = LIS_SUCCESS;
  for (i=0;i<ss;i++)
    {
      lis_vector_duplicate(A, &esolver->evector[i]); 
      esolver2->ishift = esolver->evalue[i];
      lis_esolve(A, esolver->evector[i], &evalue, esolver2);
      lis_esolver_work_destroy(esolver2); 
      esolver->evalue[i] = evalue;
      esolver->iter[i] = esolver2->iter[0];

      err = lis_esolver_standard_residual(
              A, esolver->evector[i], esolver->evalue[i],
              r, v[0], &cert_resid);
      if( err )
        {
          esolver->retcode = err;
          lis_esolver_destroy(esolver2);
          lis_free(t);
          lis_free(tq);
          lis_free(tr);
          lis_solver_destroy(solver);
          LIS_DEBUG_FUNC_OUT;
          return err;
        }

      esolver->resid[i] = cert_resid;
      if( cert_resid>tol && final_status==LIS_SUCCESS )
        final_status = LIS_MAXITER;

      if (i==0) 
	{
	  evalue0 = esolver->evalue[0];
	  iter0 = esolver2->iter[0];
	  resid0 = esolver->resid[0];
	  if( output & LIS_EPRINT_MEM ) 
	    {
	      for (ic=0;ic<iter0+1;ic++)
		{
		  esolver->rhistory[ic] = esolver2->rhistory[ic]; 
		}
	    }
	  esolver->time = esolver2->time;	  
	  esolver->ptime += esolver2->ptime;
	  esolver->itime += esolver2->itime;
	  esolver->p_c_time += esolver2->p_c_time;
	  esolver->p_i_time += esolver2->p_i_time;
	}

      if (output) 
	{
	  lis_printf(comm,"Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX
	  lis_printf(comm,"Lanczos: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else	  
	  lis_printf(comm,"Lanczos: eigenvalue           = %e\n", (double)esolver->evalue[i]);
#endif
	  lis_printf(comm,"Lanczos: elapsed time         = %e sec.\n", esolver2->time);	  	  
	  lis_printf(comm,"Lanczos: number of iterations = %D\n",esolver2->iter[0]);
	  lis_printf(comm,"Lanczos: relative residual    = %e\n\n",(double)esolver->resid[i]);
	}
    }
  esolver->evalue[0] = evalue0; 
  esolver->iter[0] = iter0;
  esolver->resid[0] = resid0;

  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_esolver_destroy(esolver2); 

  lis_free(t); 
  lis_free(tq);
  lis_free(tr);

  lis_solver_destroy(solver);

  LIS_DEBUG_FUNC_OUT;
  esolver->retcode = final_status;
  return final_status;
}

/***************************************
 * Generalized Lanczos Iteration       *
 ***************************************
 q       = (1,...,1)^T
 ***************************************
 r = B * q
 beta(0) = |<q,r>|^1/2
 for j=1,2,...
   w(j)      = r / beta(j-1)
   v(j)      = q / beta(j-1)
   r         = A * v(j)
   r         = r - beta(j-1) * w(j-1)
   alpha(j)  = <v(j), r>
   r         = r - alpha(j) * w(j)
   reorthogonalization 
   solve B * q = r for q
   beta(j)   = |<q,r>|^1/2
 end for
 compute eigenvalues of a symmetric tridiagonal matrix 
 T(j) = ST'(j)S^*, where
     (alpha(1) beta(1)                      )
     (beta(1)  alpha(2)                     )
 T = (               ...                    )
     (                  alpha(j-1) beta(j-1))                       
     (                  beta(j-1)  alpha(j) )
 compute refined eigenpairs
 ***************************************/

#undef NWORK
#define NWORK 4
#undef __FUNC__
#define __FUNC__ "lis_egli_check_params"
LIS_INT lis_egli_check_params(LIS_ESOLVER esolver)
{
        LIS_INT ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
	if( ss<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%D) is less than 0\n",ss);
		return LIS_ERR_ILL_ARG;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_egli_malloc_work"
LIS_INT lis_egli_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR *work;
	LIS_INT	i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	/*
	 * The untransformed generalized Lanczos recurrence keeps both
	 * w_j = B*v_j and v_j for full reorthogonalization.  They must
	 * not alias.  The transformed path uses one Lanczos basis.
	 */
	if( esolver->transform==LIS_EIGEN_TRANSFORM_AUTO )
	{
		worklen = NWORK + 2*ss;
	}
	else
	{
		worklen = NWORK + ss;
	}
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_egli_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( esolver->eprecision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(esolver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
			if( err ) break;
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	esolver->worklen = worklen;
	esolver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

static LIS_INT
lis_egli_solve_with_retry(
        LIS_MATRIX S,
        LIS_VECTOR b,
        LIS_VECTOR x,
        LIS_SOLVER solver,
        LIS_PRECON *precon,
        LIS_REAL acceptance_tol)
{
        LIS_INT err,status;

        status =
                lis_esolver_inner_solve_attempt(
                        S,b,x,solver,*precon);
        if( status==LIS_SUCCESS )
        {
                return LIS_SUCCESS;
        }

        if( !lis_esolver_inner_solve_retryable(status) )
        {
                return status;
        }

        err = lis_esolver_rebuild_precon(
                S,solver,precon);
        if( err )
        {
                return err;
        }

        lis_vector_set_all(0.0,x);

        status =
                lis_esolver_inner_solve_attempt(
                        S,b,x,solver,*precon);

        return lis_esolver_inner_solve_certify(
                S,b,x,solver,
                status,
                acceptance_tol);
}


static LIS_INT
lis_egli_replace_shifted_operator(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_SCALAR shift,
        LIS_SOLVER solver,
        LIS_PRECON *precon,
        LIS_MATRIX *S,
        LIS_MATRIX *Sop)
{
        LIS_INT err;
        LIS_MATRIX trial_operator;

        trial_operator = NULL;

        err = lis_matrix_create_operator(
                1.0,A,-shift,B,
                &trial_operator);
        if( err )
        {
                return err;
        }

        err = lis_esolver_rebuild_precon(
                trial_operator,
                solver,
                precon);
        if( err )
        {
                /*
                 * The rebuild helper has detached work/preconditioner
                 * state tied to the previous operator.  Keep solver->A
                 * valid for cleanup before discarding the failed trial.
                 */
                solver->A = *S;
                lis_matrix_destroy(trial_operator);
                return err;
        }

        if( *Sop )
        {
                lis_matrix_destroy(*Sop);
        }

        *Sop = trial_operator;
        *S = trial_operator;

        return LIS_SUCCESS;
}


static LIS_INT
lis_egli_first_solve_with_recovery(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_VECTOR b,
        LIS_VECTOR x,
        LIS_SOLVER solver,
        LIS_PRECON *precon,
        LIS_MATRIX *S,
        LIS_MATRIX *Sop,
        LIS_SCALAR requested_shift,
        LIS_REAL acceptance_tol,
        LIS_SCALAR *active_shift)
{
        LIS_INT err,status,attempt;
        LIS_SCALAR perturbation;
        LIS_SCALAR trial_shift;

        *active_shift = requested_shift;

        status =
                lis_esolver_inner_solve_attempt(
                        *S,b,x,solver,*precon);
        if( status==LIS_SUCCESS )
        {
                return LIS_SUCCESS;
        }

        if( !lis_esolver_inner_solve_retryable(status) )
        {
                return status;
        }

        err = lis_esolver_rebuild_precon(
                *S,solver,precon);
        if( err )
        {
                return err;
        }

        lis_vector_set_all(0.0,x);

        status =
                lis_esolver_inner_solve_attempt(
                        *S,b,x,solver,*precon);

                status =
                        lis_esolver_inner_solve_certify(
                                *S,b,x,solver,
                                status,
                                acceptance_tol);
        if( status==LIS_SUCCESS )
        {
                return LIS_SUCCESS;
        }

        if( !lis_esolver_inner_solve_retryable(status) )
        {
                return status;
        }

        perturbation =
                lis_esolver_shift_perturbation(
                        requested_shift);
        if( perturbation==0.0 )
        {
                return status;
        }

        for(attempt=0;attempt<2;attempt++)
        {
                trial_shift =
                        attempt==0 ?
                        requested_shift+perturbation :
                        requested_shift-perturbation;

                err = lis_egli_replace_shifted_operator(
                        A,B,trial_shift,
                        solver,precon,
                        S,Sop);
                if( err )
                {
                        return err;
                }

                lis_vector_set_all(0.0,x);

                status =
                        lis_esolver_inner_solve_attempt(
                                *S,b,x,solver,*precon);

                status =
                        lis_esolver_inner_solve_certify(
                                *S,b,x,solver,
                                status,
                                acceptance_tol);

                if( status==LIS_SUCCESS )
                {
                        *active_shift = trial_shift;
                        return LIS_SUCCESS;
                }

                if( !lis_esolver_inner_solve_retryable(status) )
                {
                        return status;
                }
        }

        return status;
}


#undef __FUNC__
#define __FUNC__ "lis_egli_transformed"
static LIS_INT lis_egli_transformed(LIS_ESOLVER esolver)
{
	LIS_Comm comm;
	LIS_MATRIX A,B,M,S,Sop;
	LIS_VECTOR q,r,z,*v;
	LIS_SOLVER solver;
	LIS_PRECON precon;
	LIS_ESOLVER esolver2;
	LIS_SCALAR *t,*tq,*tr;
	LIS_SCALAR oshift,ishift,active_shift;
	LIS_SCALAR alpha,dot,tmp,nu;
	LIS_SCALAR evalue,evalue0;
	LIS_REAL beta,beta_prev,mag,maxmag;
	LIS_REAL tol,projerr;
	LIS_REAL resid,resid0,nrm2;
	LIS_INT transform,ss,nev,output,rval;
	LIS_INT nigesolver,emaxiter;
	LIS_INT i,j,k,imax,ic;
	LIS_INT projiter,iter0;
	LIS_INT err,inner_err;
	LIS_INT nsol,precon_type;
	LIS_INT final_status;
	char solvername[128],preconname[128],esolvername[128];
	double time,time0;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A = esolver->A;
	B = esolver->B;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
        nev = esolver->nev>0 ? esolver->nev : ss;
	emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
	tol = esolver->params[
		LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN];

	output = esolver->options[LIS_EOPTIONS_OUTPUT];
	rval = esolver->options[LIS_EOPTIONS_RVAL];

	nigesolver =
		esolver->options[
			LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER];

	transform = esolver->transform;

#ifdef _COMPLEX
	oshift =
		esolver->params[
			LIS_EPARAMS_SHIFT-LIS_EOPTIONS_LEN] +
		esolver->params[
			LIS_EPARAMS_SHIFT_IM-LIS_EOPTIONS_LEN]
			*_Complex_I;
#else
	oshift =
		esolver->params[
			LIS_EPARAMS_SHIFT-LIS_EOPTIONS_LEN];
#endif

	ishift = esolver->ishift;

	if( ishift!=0.0 )
	{
		oshift = ishift;
	}

	active_shift = oshift;

	if( ss<1 )
	{
		LIS_SETERR1(
			LIS_ERR_ILL_ARG,
			"Parameter LIS_EOPTIONS_SUBSPACE is %D "
			"(Set greater than 0 for transformed Lanczos)\n",
			ss);
		return LIS_ERR_ILL_ARG;
	}

	if( nev<1 || nev>ss )
	{
	        LIS_SETERR2(
	                LIS_ERR_ILL_ARG,
	                "Parameter -nev is %D "
	                "(Set between 1 and subspace size %D)\n",
	                nev,
	                ss);
	        return LIS_ERR_ILL_ARG;
	}

	if( transform==LIS_EIGEN_TRANSFORM_BUCKLING )
	{
		M = A;
	}
	else
	{
		M = B;
	}

	S = A;
	Sop = NULL;
	solver = NULL;
	precon = NULL;
	esolver2 = NULL;

	t = NULL;
	tq = NULL;
	tr = NULL;

	t = (LIS_SCALAR *)lis_malloc(
		ss*ss*sizeof(LIS_SCALAR),
		"lis_egli_transformed::t");

	tq = (LIS_SCALAR *)lis_malloc(
		ss*ss*sizeof(LIS_SCALAR),
		"lis_egli_transformed::tq");

	tr = (LIS_SCALAR *)lis_malloc(
		ss*ss*sizeof(LIS_SCALAR),
		"lis_egli_transformed::tr");

	if( t==NULL || tq==NULL || tr==NULL )
	{
		LIS_SETERR_MEM(ss*ss*sizeof(LIS_SCALAR));
		err = LIS_ERR_OUT_OF_MEMORY;
		goto cleanup;
	}

	for(i=0;i<ss*ss;i++)
	{
		t[i] = 0.0;
		tq[i] = 0.0;
		tr[i] = 0.0;
	}

	if( oshift!=0.0 )
	{
		err = lis_matrix_create_operator(
			1.0,A,-oshift,B,&Sop);

		if( err ) goto cleanup;

		S = Sop;
	}

	err = lis_solver_create(&solver);
	if( err ) goto cleanup;

	lis_solver_set_option("-i bicg -p none",solver);

	err = lis_solver_set_optionC(solver);
	if( err ) goto cleanup;

	lis_solver_get_solver(solver,&nsol);
	lis_solver_get_precon(solver,&precon_type);

	lis_solver_get_solvername(
		nsol,solvername);

	lis_solver_get_preconname(
		precon_type,preconname);

	lis_esolver_get_esolvername(
		nigesolver,esolvername);

	if( output )
	{
#ifdef _COMPLEX
		lis_printf(
			comm,
			"shift                 : (%e, %e)\n",
			(double)creal(oshift),
			(double)cimag(oshift));
#else
		lis_printf(
			comm,
			"shift                 : %e\n",
			(double)oshift);
#endif

		if( transform==
		    LIS_EIGEN_TRANSFORM_SHIFT_INVERT )
		{
			lis_printf(
				comm,
				"spectral transform    : shift-invert\n");
		}
		else
		{
			lis_printf(
				comm,
				"spectral transform    : buckling\n");
		}

		lis_printf(
			comm,
			"linear solver         : %s\n",
			solvername);

		lis_printf(
			comm,
			"preconditioner        : %s\n",
			preconname);

		lis_printf(
			comm,
			"size of subspace      : %D\n",
			ss);
	}

	solver->A = S;

	err = lis_precon_create(
		solver,&precon);

	if( err ) goto cleanup;

	q = esolver->work[0];
	r = esolver->work[1];
	z = esolver->work[2];
	v = &esolver->work[3];

	lis_vector_set_all(0.0,v[0]);
	lis_vector_set_all(1.0,q);

	/*
	 * Normalize the starting vector in the M inner product.
	 */
	err = lis_matvec(M,q,r);
	if( err ) goto cleanup;

	err = lis_vector_dot(q,r,&dot);
	if( err ) goto cleanup;

#ifdef _COMPLEX
	beta = sqrt(cabs(dot));
#else
	beta = sqrt(fabs(dot));
#endif

	if( !isfinite((double)beta) || beta<=0.0 )
	{
		LIS_SETERR(
			LIS_BREAKDOWN,
			"invalid transformed Lanczos starting norm\n");
		err = LIS_BREAKDOWN;
		goto cleanup;
	}

	lis_vector_copy(q,v[1]);
	lis_vector_scale(1.0/beta,v[1]);

	beta_prev = 0.0;

	time0 = lis_wtime();

	for(j=1;j<=ss;j++)
	{
		/*
		 * r = M v_j
		 */
		err = lis_matvec(M,v[j],r);
		if( err ) goto cleanup;

		/*
		 * q = (A-sigma B)^-1 M v_j
		 *
		 * Only the first transformed solve may adopt a bounded
		 * perturbed shift.  Once the Lanczos basis is started, all
		 * remaining solves use the same transformed operator.
		 */
		if( j==1 )
		{
			err = lis_egli_first_solve_with_recovery(
				A,B,r,q,solver,&precon,
				&S,&Sop,oshift,tol,&active_shift);
		}
		else
		{
			err = lis_egli_solve_with_retry(
				S,r,q,solver,&precon,tol);
		}

		if( err ) goto cleanup;

		esolver->ptime += solver->ptime;
		esolver->itime += solver->itime;
		esolver->p_c_time += solver->p_c_time;
		esolver->p_i_time += solver->p_i_time;

		/*
		 * alpha_j = <v_j,Tv_j>_M
		 */
		err = lis_vector_dot(r,q,&alpha);
		if( err ) goto cleanup;

		t[(j-1)*ss+(j-1)] = alpha;

		if( j>1 )
		{
			lis_vector_axpy(
				-beta_prev,
				v[j-1],
				q);
		}

		lis_vector_axpy(
			-alpha,
			v[j],
			q);

		/*
		 * One full M-orthogonalization pass.
		 */
		err = lis_matvec(M,q,r);
		if( err ) goto cleanup;

		for(k=1;k<=j;k++)
		{
			err = lis_vector_dot(
				v[k],r,&dot);

			if( err ) goto cleanup;

			lis_vector_axpy(
				-dot,
				v[k],
				q);
		}

		if( j==ss )
		{
			break;
		}

		err = lis_matvec(M,q,r);
		if( err ) goto cleanup;

		err = lis_vector_dot(q,r,&dot);
		if( err ) goto cleanup;

#ifdef _COMPLEX
		beta = sqrt(cabs(dot));
#else
		beta = sqrt(fabs(dot));
#endif

		if( !isfinite((double)beta) ||
		    beta<=tol )
		{
			LIS_SETERR(
				LIS_BREAKDOWN,
				"transformed Lanczos invariant subspace "
				"closed before requested subspace size\n");
			err = LIS_BREAKDOWN;
			goto cleanup;
		}

		t[(j-1)*ss+j] = (LIS_SCALAR)beta;
		t[j*ss+j-1] = (LIS_SCALAR)beta;

		lis_vector_copy(
			q,
			v[j+1]);

		lis_vector_scale(
			1.0/beta,
			v[j+1]);

		beta_prev = beta;
	}

	time = lis_wtime()-time0;

	/*
	 * Compute the complete eigensystem of the projected
	 * symmetric tridiagonal matrix.
	 */
	for(i=0;i<ss*ss;i++)
	{
		tq[i] = t[i];
	}

	err = lis_eigen_projected_solve(
		ss,
		tq,
		esolver->evalue,
		tr,
		&projiter,
		&projerr);

	if( err ) goto cleanup;

	/*
	 * Sort transformed eigenpairs by decreasing |nu|.
	 */
	for(i=0;i<ss-1;i++)
	{
		imax = i;

#ifdef _COMPLEX
		maxmag = cabs(esolver->evalue[i]);
#else
		maxmag = fabs(esolver->evalue[i]);
#endif

		for(j=i+1;j<ss;j++)
		{
#ifdef _COMPLEX
			mag = cabs(esolver->evalue[j]);
#else
			mag = fabs(esolver->evalue[j]);
#endif

			if( mag>maxmag )
			{
				maxmag = mag;
				imax = j;
			}
		}

		if( imax!=i )
		{
			tmp = esolver->evalue[i];
			esolver->evalue[i] =
				esolver->evalue[imax];
			esolver->evalue[imax] = tmp;

			for(k=0;k<ss;k++)
			{
				tmp = tr[k*ss+i];
				tr[k*ss+i] =
					tr[k*ss+imax];
				tr[k*ss+imax] = tmp;
			}
		}
	}

	/*
	 * Map nu to the original generalized eigenvalue lambda.
	 */
	for(i=0;i<ss;i++)
	{
		nu = esolver->evalue[i];

		if( transform==
		    LIS_EIGEN_TRANSFORM_SHIFT_INVERT )
		{
#ifdef _COMPLEX
			mag = cabs(nu);
#else
			mag = fabs(nu);
#endif

			if( mag<=1.0e-30 )
			{
				LIS_SETERR(
					LIS_BREAKDOWN,
					"zero shift-invert Ritz value\n");
				err = LIS_BREAKDOWN;
				goto cleanup;
			}

			esolver->evalue[i] =
				active_shift+1.0/nu;
		}
		else
		{
#ifdef _COMPLEX
			mag = cabs(nu-1.0);
#else
			mag = fabs(nu-1.0);
#endif

			if( mag<=1.0e-30 )
			{
				LIS_SETERR(
					LIS_BREAKDOWN,
					"buckling Ritz value is too close to one\n");
				err = LIS_BREAKDOWN;
				goto cleanup;
			}

			esolver->evalue[i] =
				active_shift*nu/(nu-1.0);
		}

		esolver->iter[i] = ss;
		esolver->resid[i] = 0.0;
	}

	if( output )
	{
		lis_printf(
			comm,
			"projected eigensolver iterations : %D\n",
			projiter);

		lis_printf(
			comm,
			"projected eigensolver error      : %e\n\n",
			(double)projerr);

		lis_printf(
			comm,
			"transformed Lanczos Ritz values:\n\n");

		for(i=0;i<ss;i++)
		{
#ifdef _COMPLEX
			lis_printf(
				comm,
				"Generalized Lanczos: mode number = %D, "
				"eigenvalue = (%e, %e)\n",
				i,
				(double)creal(esolver->evalue[i]),
				(double)cimag(esolver->evalue[i]));
#else
			lis_printf(
				comm,
				"Generalized Lanczos: mode number = %D, "
				"eigenvalue = %e\n",
				i,
				(double)esolver->evalue[i]);
#endif
		}

		lis_printf(
			comm,
			"Generalized Lanczos: elapsed time = %e sec.\n\n",
			time);
	}

	/*
	 * -rval keeps the established GLI behaviour: return the
	 * Ritz values without constructing eigenvectors.
	 */
	if( rval )
	{
		esolver->retcode = LIS_SUCCESS;
		err = LIS_SUCCESS;
		goto cleanup;
	}

	/*
	 * Build physical Ritz vectors:
	 *
	 *     x_i = V y_i
	 *
	 * where columns of tr are projected eigenvectors.
	 */
	final_status = LIS_SUCCESS;
	evalue0 = 0.0;
	iter0 = ss;
	resid0 = 0.0;

	err = lis_esolver_create(&esolver2);
	if( err ) goto cleanup;

	esolver2->options[LIS_EOPTIONS_ESOLVER] =
		nigesolver;

	esolver2->options[LIS_EOPTIONS_SUBSPACE] =
		1;

	esolver2->options[LIS_EOPTIONS_MAXITER] =
		emaxiter;

	esolver2->options[LIS_EOPTIONS_OUTPUT] =
		esolver->options[LIS_EOPTIONS_OUTPUT];

	esolver2->options[LIS_EOPTIONS_INITGUESS_ONES] =
		0;

	esolver2->params[
		LIS_EPARAMS_RESID-LIS_EOPTIONS_LEN] =
		tol;

	esolver2->transform =
		LIS_EIGEN_TRANSFORM_AUTO;

	if( output )
	{
		lis_printf(
			comm,
			"inner eigensolver     : %s\n\n",
			esolvername);
	}

	for(i=0;i<nev;i++)
	{
		err = lis_vector_duplicate(
			A,
			&esolver->evector[i]);

		if( err ) goto cleanup;

		err = lis_eigen_ritz_reconstruct(
			ss,
			&v[1],
			1,
			&tr[i],
			ss,
			&esolver->evector[i]);

		if( err ) goto cleanup;

		err = lis_vector_nrm2(
			esolver->evector[i],
			&nrm2);

		if( err ) goto cleanup;

		if( !isfinite((double)nrm2) ||
		    nrm2<=0.0 )
		{
			LIS_SETERR(
				LIS_BREAKDOWN,
				"invalid transformed generalized Lanczos Ritz vector norm\n");
			err = LIS_BREAKDOWN;
			goto cleanup;
		}

		lis_vector_scale(
			1.0/nrm2,
			esolver->evector[i]);

		err = lis_esolver_generalized_residual(
			A,
			B,
			esolver->evector[i],
			esolver->evalue[i],
			q,
			r,
			z,
			&resid);

		if( err ) goto cleanup;


		/*
		 * Refine only when the physical generalized residual
		 * is not already converged.
		 */
		if( resid>tol )
		{
			esolver2->ishift =
				esolver->evalue[i];

			/*
			 * Keep the inner inverse-iteration system away from the
			 * nearly singular Ritz pole.  Use one percent of the
			 * nearest distinct Ritz-value spacing and move toward
			 * the original spectral target.
			 */
#ifndef _COMPLEX
			{
			        LIS_REAL gap;
			        LIS_REAL d;
			        LIS_REAL delta;
			        LIS_INT jj;

			        gap = 0.0;

			        for(jj=0;jj<ss;jj++)
			        {
			                if( jj==i )
			                        continue;

			                d = fabs(
			                        (double)(
			                                esolver->evalue[jj] -
			                                esolver->evalue[i]));

			                if( d>0.0 &&
			                    (gap==0.0 || d<gap) )
			                {
			                        gap = d;
			                }
			        }

			        if( gap>0.0 )
			        {
			                delta = 0.01*gap;

			                if( oshift<esolver->evalue[i] )
			                        delta = -delta;

			                esolver2->ishift +=
			                        (LIS_SCALAR)delta;
			        }
			}
#endif


			inner_err = lis_gesolve(
				A,
				B,
				esolver->evector[i],
				&evalue,
				esolver2);

			if( inner_err!=LIS_SUCCESS &&
			    inner_err!=LIS_MAXITER )
			{
				err = inner_err;
				goto cleanup;
			}

			esolver->evalue[i] =
				evalue;

			esolver->iter[i] =
				esolver2->iter[0];

			err = lis_esolver_generalized_residual(
				A,
				B,
				esolver->evector[i],
				esolver->evalue[i],
				q,
				r,
				z,
				&resid);

			if( err ) goto cleanup;

			if( inner_err==LIS_MAXITER &&
			    resid>tol )
			{
				final_status = LIS_MAXITER;
			}

			if( i==0 )
			{
				esolver->ptime +=
					esolver2->ptime;

				esolver->itime +=
					esolver2->itime;

				esolver->p_c_time +=
					esolver2->p_c_time;

				esolver->p_i_time +=
					esolver2->p_i_time;

				if( output & LIS_EPRINT_MEM )
				{
					for(ic=0;
					    ic<esolver2->iter[0]+1;
					    ic++)
					{
						esolver->rhistory[ic] =
							esolver2->rhistory[ic];
					}
				}
			}

			lis_esolver_work_destroy(
				esolver2);
		}

		esolver->resid[i] = resid;

		if( resid>tol )
		{
			final_status = LIS_MAXITER;
		}

		if( output )
		{
			lis_printf(
				comm,
				"Generalized Lanczos: mode number          = %D\n",
				i);

#ifdef _COMPLEX
			lis_printf(
				comm,
				"Generalized Lanczos: eigenvalue           = (%e, %e)\n",
				(double)creal(esolver->evalue[i]),
				(double)cimag(esolver->evalue[i]));
#else
			lis_printf(
				comm,
				"Generalized Lanczos: eigenvalue           = %e\n",
				(double)esolver->evalue[i]);
#endif

			lis_printf(
				comm,
				"Generalized Lanczos: true residual        = %e\n\n",
				(double)esolver->resid[i]);
		}

		if( i==0 )
		{
			evalue0 = esolver->evalue[0];
			iter0 = esolver->iter[0];
			resid0 = esolver->resid[0];
		}
	}

	esolver->evalue[0] = evalue0;
	esolver->iter[0] = iter0;
	esolver->resid[0] = resid0;

	lis_vector_copy(
		esolver->evector[0],
		esolver->x);

	esolver->time = time;
	esolver->retcode = final_status;

	err = final_status;

cleanup:

	if( err!=LIS_SUCCESS )
	{
		esolver->retcode = err;
	}

	if( esolver2 )
	{
		lis_esolver_destroy(esolver2);
	}

	if( precon )
	{
		lis_precon_destroy(precon);
	}

	if( solver )
	{
		lis_solver_destroy(solver);
	}

	if( Sop )
	{
		lis_matrix_destroy(Sop);
	}

	if( t ) lis_free(t);
	if( tq ) lis_free(tq);
	if( tr ) lis_free(tr);

	LIS_DEBUG_FUNC_OUT;
	return err;
}


/*
 * Generalized Block Lanczos dispatch boundary.
 *
 * F1 introduces the internal block entry point only.  The numerical
 * generalized block recurrence is added in later stages so that scalar
 * generalized Lanczos remains unchanged while block requests have a
 * dedicated, auditable path.
 *
 * Transformed generalized block modes are intentionally rejected here.
 * In particular, a future buckling block path must preserve the existing
 * transformed-operator semantics rather than treating the geometric
 * stiffness matrix as a positive-definite block metric.
 */
#undef __FUNC__
#define __FUNC__ "lis_egli_block"
LIS_INT lis_egli_block(LIS_ESOLVER esolver)
{
    LIS_DEBUG_FUNC_IN;

    if( esolver==NULL ||
        esolver->block_size<=1 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid generalized Block Lanczos dispatch arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    if( esolver->transform!=LIS_EIGEN_TRANSFORM_AUTO )
    {
        LIS_SETERR(
            LIS_ERR_NOT_IMPLEMENTED,
            "transformed generalized Block Lanczos is not implemented yet\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }

    LIS_SETERR(
        LIS_ERR_NOT_IMPLEMENTED,
        "generalized Block Lanczos recurrence is not implemented yet\n");

    return LIS_ERR_NOT_IMPLEMENTED;
}


#undef __FUNC__
#define __FUNC__ "lis_egli"
LIS_INT lis_egli(LIS_ESOLVER esolver)
{
  LIS_Comm comm;

  if( esolver->block_size>1 )
    {
      return lis_egli_block(esolver);
    }

  if( esolver->transform != LIS_EIGEN_TRANSFORM_AUTO )
    {
      return lis_egli_transformed(esolver);
    }
  LIS_INT err;
  LIS_MATRIX A,B;
  LIS_INT ss,ic;
  LIS_INT emaxiter,iter0,qriter;
  LIS_REAL tol,qrerr;
  double time,time0;    
  LIS_INT i,j,k;
  LIS_INT output, nigesolver;
  LIS_REAL resid0;
  LIS_SCALAR beta,dot;
  LIS_VECTOR q,r,*w,*v;
  LIS_SCALAR *t,*tq,*tr,evalue,evalue0;
  LIS_SOLVER solver;
  LIS_PRECON precon;
  LIS_ESOLVER esolver2;
  char esolvername[128],solvername[128],preconname[128];
  LIS_INT nsol,precon_type;
  LIS_INT rval;  

  LIS_DEBUG_FUNC_IN;

  comm = LIS_COMM_WORLD;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  nigesolver = esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER];
  rval = esolver->options[LIS_EOPTIONS_RVAL];    

  t = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egli::t");
  tq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egli::tq");
  tr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_egli::tr");
  
  A = esolver->A;
  B = esolver->B;
  q = esolver->work[0];  
  r = esolver->work[1];
  w = &esolver->work[2];
  v = &esolver->work[2+ss];
  lis_vector_set_all(0.0,w[0]);  
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,q);

  /* create linear solver */
  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p none",solver);
  err = lis_solver_set_optionC(solver);
  CHKERR(err);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
  lis_esolver_get_esolvername(nigesolver, esolvername);
  if( output )
    {
      lis_printf(comm,"inner eigensolver     : %s\n", esolvername);
      lis_printf(comm,"linear solver         : %s\n", solvername);
      lis_printf(comm,"preconditioner        : %s\n", preconname);
    }

  /* create preconditioner */
  solver->A = B;
  err = lis_precon_create(solver, &precon);
  if( err )
    {
      lis_solver_work_destroy(solver);
      solver->retcode = err;
      return err;
    }

  for (i=0;i<ss*ss;i++) t[i] = 0.0;

  j=0;
  while (j<ss-1)
    {
      j = j+1;

      /* w(j) = r / beta(j-1) */      
      /* v(j) = r / beta(j-1) */
      /* r = A * v(j) */
      /* r = r - beta(j-1) * w(j-1) */
      
      if (j==1) {
 	lis_matvec(B, q, r);
	lis_vector_dot(q, r, &beta);
	beta = sqrt(fabs(beta));
	lis_vector_copy(r, w[j]);
	lis_vector_copy(q, v[j]);	
	lis_vector_scale(1.0/beta, w[j]);
	lis_vector_scale(1.0/beta, v[j]);	
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-beta, w[j-1], r); 	
      }
      else {
	lis_vector_copy(r, w[j]);
	lis_vector_copy(q, v[j]);	
	lis_vector_scale(1.0/t[(j-2)*ss+j-1], w[j]);	
	lis_vector_scale(1.0/t[(j-2)*ss+j-1], v[j]);
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-t[(j-2)*ss+j-1], w[j-1], r); 
      }

      /* alpha(j)  = <v(j), r> */
      lis_vector_dot(v[j], r, &t[(j-1)*ss+j-1]); 

      /* r = r - alpha(j) * w(j) */
      lis_vector_axpy(-t[(j-1)*ss+j-1], w[j], r); 

      /*
       * Reorthogonalize the dual residual.  With w(k)=B*v(k),
       * subtracting <v(k),r>*w(k) makes q=B^-1*r orthogonal
       * to the existing basis in the B inner product.
       */
      err = lis_lanczos_reorthogonalize_dual_residual(
          v,
          w,
          1,
          j,
          r);
      if( err )
          return err;

      /* solve B * q = r */
      err = lis_solve_kernel(B, r, q, solver, precon);
      if( err )
	{
	  lis_solver_work_destroy(solver);	  
	  solver->retcode = err;
	  lis_precon_destroy(precon);
	  return err;
	}

      /* beta(j) = |<q,r>|^1/2 */
      lis_vector_dot(q, r, &beta);
      beta = sqrt(fabs(beta));
      t[(j-1)*ss+j] = (LIS_REAL)beta;
      
      /* convergence check */
      if (fabs(t[(j-1)*ss+j])<tol) break;  
      t[j*ss+j-1] = t[(j-1)*ss+j];
    }

  /* compute eigenvalues of a symmetric tridiagonal matrix 
     T(j) = ST'(j)S^* */
  time0 = lis_wtime();
  lis_array_qr(ss,t,tq,tr,&qriter,&qrerr);
  time = lis_wtime() - time0;

  for (i=0;i<ss;i++)
    {
      esolver->evalue[i] = t[i*ss+i];
    }

  if( output ) 
    {
      lis_printf(comm,"size of subspace      : %D\n\n", ss);
      lis_printf(comm,"Ritz values :\n\n");
      for (i=0;i<ss;i++)
	{
	  lis_printf(comm,"Generalized Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX	  
	  lis_printf(comm,"Generalized Lanczos: Ritz value           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else
	  lis_printf(comm,"Generalized Lanczos: Ritz value           = %e\n", (double)esolver->evalue[i]);
#endif	  
	}
      lis_printf(comm,"Lanczos: elapsed time         = %e sec.\n\n", time);  
    }
  if( rval ) return LIS_SUCCESS;
      
  if( output ) 
    {
      lis_printf(comm,"computing refined eigenpairs using inner eigensolver:\n\n");      
    }

  lis_esolver_create(&esolver2);
  esolver2->options[LIS_EOPTIONS_ESOLVER] = nigesolver;
  esolver2->options[LIS_EOPTIONS_SUBSPACE] = 1;
  esolver2->options[LIS_EOPTIONS_MAXITER] = emaxiter;
  esolver2->options[LIS_EOPTIONS_OUTPUT] = esolver->options[LIS_EOPTIONS_OUTPUT];
  esolver2->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = tol; 

  /* compute refined eigenpairs */
  for (i=0;i<ss;i++)
    {
      lis_vector_duplicate(A, &esolver->evector[i]); 
      esolver2->ishift = esolver->evalue[i];
      lis_gesolve(A, B, esolver->evector[i], &evalue, esolver2);
      lis_esolver_work_destroy(esolver2); 
      esolver->evalue[i] = evalue;
      esolver->iter[i] = esolver2->iter[0];      
      esolver->resid[i] = esolver2->resid[0];

      if (i==0) 
	{
	  evalue0 = esolver->evalue[0];
	  iter0 = esolver2->iter[0];
	  resid0 = esolver2->resid[0];
	  if( output & LIS_EPRINT_MEM ) 
	    {
	      for (ic=0;ic<iter0+1;ic++)
		{
		  esolver->rhistory[ic] = esolver2->rhistory[ic]; 
		}
	    }
	  esolver->time = esolver2->time;	  
	  esolver->ptime += esolver2->ptime;
	  esolver->itime += esolver2->itime;
	  esolver->p_c_time += esolver2->p_c_time;
	  esolver->p_i_time += esolver2->p_i_time;
	}

      if (output) 
	{
	  lis_printf(comm,"Generalized Lanczos: mode number          = %D\n", i);
#ifdef _COMPLEX
	  lis_printf(comm,"Generalized Lanczos: eigenvalue           = (%e, %e)\n", (double)creal(esolver->evalue[i]), (double)cimag(esolver->evalue[i]));
#else	  
	  lis_printf(comm,"Generalized Lanczos: eigenvalue           = %e\n", (double)esolver->evalue[i]);
#endif
	  lis_printf(comm,"Generalized Lanczos: elapsed time         = %e sec.\n", esolver2->time);	  	  
	  lis_printf(comm,"Generalized Lanczos: number of iterations = %D\n",esolver2->iter[0]);
	  lis_printf(comm,"Generalized Lanczos: relative residual    = %e\n\n",(double)esolver2->resid[0]);
	}
    }
  esolver->evalue[0] = evalue0; 
  esolver->iter[0] = iter0;
  esolver->resid[0] = resid0;

  lis_vector_copy(esolver->evector[0], esolver->x);

  lis_esolver_destroy(esolver2); 

  lis_free(t); 
  lis_free(tq);
  lis_free(tr);
  
  lis_precon_destroy(precon);
  lis_solver_destroy(solver);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}



