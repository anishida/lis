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
#include <string.h>
#include <math.h>
#include <float.h>

#include "lislib.h"

/*
 * Native FEAST eigensolver foundation.
 *
 * FEAST 4.0 is used as an algorithmic reference.
 * This implementation has no external FEAST library dependency.
 *
 * Stage F0 provides registration, interval configuration and
 * parameter validation. The contour projector and shifted systems
 * are implemented in later stages.
 */

#define LIS_FEAST_INTERVAL_EMIN 0x01
#define LIS_FEAST_INTERVAL_EMAX 0x02

#define LIS_FEAST_PI \
        ((LIS_REAL)3.141592653589793238462643383279502884)


/*
 * Compute the n-point Gauss-Legendre rule on [-1,1].
 *
 * The implementation uses Newton refinement of the Legendre roots
 * and does not depend on a table of precomputed quadrature nodes.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_gauss_legendre"

LIS_INT lis_efeast_gauss_legendre(
        LIS_INT n,
        LIS_REAL *x,
        LIS_REAL *w)
{
        LIS_INT i;
        LIS_INT j;
        LIS_INT m;
        LIS_INT iter;

        LIS_REAL z;
        LIS_REAL zold;
        LIS_REAL p1;
        LIS_REAL p2;
        LIS_REAL p3;
        LIS_REAL pp;
        LIS_REAL wi;
        LIS_REAL tol;

        LIS_DEBUG_FUNC_IN;

        if( n<1 || x==NULL || w==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Gauss-Legendre quadrature requires n >= 1 and valid arrays\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        m = (n+1)/2;

        tol = (LIS_REAL)(64.0*DBL_EPSILON);

        for(i=0;i<m;i++)
        {
                /*
                 * Standard asymptotic initial estimate for the
                 * positive roots.
                 */
                z = cos(
                        LIS_FEAST_PI *
                        ((LIS_REAL)i + (LIS_REAL)0.75) /
                        ((LIS_REAL)n + (LIS_REAL)0.5));

                for(iter=0;iter<100;iter++)
                {
                        p1 = 1.0;
                        p2 = 0.0;

                        for(j=1;j<=n;j++)
                        {
                                p3 = p2;
                                p2 = p1;

                                p1 =
                                        (((LIS_REAL)(2*j-1))*z*p2 -
                                         ((LIS_REAL)(j-1))*p3) /
                                        (LIS_REAL)j;
                        }

                        pp =
                                ((LIS_REAL)n) *
                                (z*p1-p2) /
                                (z*z-(LIS_REAL)1.0);

                        zold = z;
                        z = zold-p1/pp;

                        if( fabs(z-zold) <=
                            tol*((LIS_REAL)1.0+fabs(z)) )
                        {
                                break;
                        }
                }

                if( iter==100 )
                {
                        LIS_SETERR1(
                                LIS_BREAKDOWN,
                                "Gauss-Legendre root %D did not converge\n",
                                i);

                        LIS_DEBUG_FUNC_OUT;
                        return LIS_BREAKDOWN;
                }

                /*
                 * Re-evaluate P_n and its derivative at the
                 * converged root before forming the weight.
                 */
                p1 = 1.0;
                p2 = 0.0;

                for(j=1;j<=n;j++)
                {
                        p3 = p2;
                        p2 = p1;

                        p1 =
                                (((LIS_REAL)(2*j-1))*z*p2 -
                                 ((LIS_REAL)(j-1))*p3) /
                                (LIS_REAL)j;
                }

                pp =
                        ((LIS_REAL)n) *
                        (z*p1-p2) /
                        (z*z-(LIS_REAL)1.0);

                wi =
                        ((LIS_REAL)2.0) /
                        (((LIS_REAL)1.0-z*z)*pp*pp);

                x[i] = -z;
                x[n-1-i] = z;

                w[i] = wi;
                w[n-1-i] = wi;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


/*
 * Generate the upper-half elliptic contour used by the native
 * real-symmetric FEAST path.
 *
 * The complex values are represented as separate real and imaginary
 * arrays so this geometry kernel is available in both real and
 * complex LIS builds.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_contour_gauss"

LIS_INT lis_efeast_contour_gauss(
        LIS_REAL emin,
        LIS_REAL emax,
        LIS_INT n,
        LIS_REAL ratio,
        LIS_REAL *z_real,
        LIS_REAL *z_imag,
        LIS_REAL *w_real,
        LIS_REAL *w_imag)
{
        LIS_REAL x[1024];
        LIS_REAL w[1024];

        LIS_REAL radius;
        LIS_REAL center;
        LIS_REAL theta;
        LIS_REAL ctheta;
        LIS_REAL stheta;

        LIS_INT i;
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        if( !(emin<emax) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour requires emin < emax\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( n<1 || n>1024 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour point count must be in [1,1024]\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !(ratio>0.0) || !isfinite((double)ratio) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour ratio must be finite and positive\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( z_real==NULL ||
            z_imag==NULL ||
            w_real==NULL ||
            w_imag==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST contour output arrays are undefined\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        err = lis_efeast_gauss_legendre(
                n,
                x,
                w);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        radius =
                (emax-emin) /
                (LIS_REAL)2.0;

        center = emin+radius;

        for(i=0;i<n;i++)
        {
                /*
                 * Map [-1,1] to the upper half ellipse:
                 *
                 * theta in [pi,0].
                 */
                theta =
                        -(LIS_FEAST_PI/(LIS_REAL)2.0)*x[i] +
                        LIS_FEAST_PI/(LIS_REAL)2.0;

                ctheta = cos(theta);
                stheta = sin(theta);

                z_real[i] =
                        center +
                        radius*ctheta;

                z_imag[i] =
                        radius*ratio*stheta;

                /*
                 * FEAST projector weight:
                 *
                 *   (1/4) * w_i * dz/dtheta
                 *
                 * after the half-contour parameter mapping.
                 */
                w_real[i] =
                        (LIS_REAL)0.25 *
                        w[i] *
                        radius *
                        ratio *
                        ctheta;

                w_imag[i] =
                        (LIS_REAL)0.25 *
                        w[i] *
                        radius *
                        stheta;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


/*
 * Apply the complex shifted operator
 *
 *     y = (z B - A) x
 *
 * for real A and B using split-complex vectors:
 *
 *     x = x_real + i x_imag
 *     z = z_real + i z_imag
 *
 * If B is NULL, the identity matrix is used.
 *
 * This representation is intentionally available in both real
 * and complex LIS builds.  It provides the foundation for the
 * real 2N block formulation required by applications that keep
 * the original matrices real.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_shifted_matvec_split"

LIS_INT lis_efeast_shifted_matvec_split(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_REAL z_real,
        LIS_REAL z_imag,
        LIS_VECTOR x_real,
        LIS_VECTOR x_imag,
        LIS_VECTOR y_real,
        LIS_VECTOR y_imag,
        LIS_VECTOR work_real,
        LIS_VECTOR work_imag)
{
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        if( A==NULL ||
            x_real==NULL ||
            x_imag==NULL ||
            y_real==NULL ||
            y_imag==NULL ||
            work_real==NULL ||
            work_imag==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted operator received an undefined argument\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !isfinite((double)z_real) ||
            !isfinite((double)z_imag) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted operator requires a finite complex shift\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        /*
         * y_real = -A*x_real
         */
        err = lis_matvec(
                A,
                x_real,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_scale(
                (LIS_SCALAR)-1.0,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        /*
         * y_imag = -A*x_imag
         */
        err = lis_matvec(
                A,
                x_imag,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_scale(
                (LIS_SCALAR)-1.0,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        /*
         * work_real = B*x_real
         * work_imag = B*x_imag
         *
         * B == NULL means B = I.
         */
        if( B!=NULL )
        {
                err = lis_matvec(
                        B,
                        x_real,
                        work_real);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }

                err = lis_matvec(
                        B,
                        x_imag,
                        work_imag);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }
        }
        else
        {
                err = lis_vector_copy(
                        x_real,
                        work_real);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }

                err = lis_vector_copy(
                        x_imag,
                        work_imag);

                if( err )
                {
                        LIS_DEBUG_FUNC_OUT;
                        return err;
                }
        }

        /*
         * Real part:
         *
         * y_real =
         *     -A*x_real
         *     + z_real*B*x_real
         *     - z_imag*B*x_imag
         */
        err = lis_vector_axpy(
                (LIS_SCALAR)z_real,
                work_real,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_axpy(
                (LIS_SCALAR)(-z_imag),
                work_imag,
                y_real);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        /*
         * Imaginary part:
         *
         * y_imag =
         *     -A*x_imag
         *     + z_imag*B*x_real
         *     + z_real*B*x_imag
         */
        err = lis_vector_axpy(
                (LIS_SCALAR)z_imag,
                work_real,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_axpy(
                (LIS_SCALAR)z_real,
                work_imag,
                y_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


/*
 * Matrix-free real 2N representation of the split-complex
 * shifted FEAST system.
 *
 * For
 *
 *     C = z B - A
 *
 * with z = z_real + i*z_imag, the split system is
 *
 *     [ z_real B - A   -z_imag B ] [ x_real ] = [ rhs_real ]
 *     [ z_imag B        z_real B-A ] [ x_imag ]   [ rhs_imag ]
 *
 * The original A and B matrices are retained without conversion.
 * Every shell matvec delegates to lis_efeast_shifted_matvec_split().
 */
typedef struct
{
        LIS_MATRIX A;
        LIS_MATRIX B;

        LIS_REAL z_real;
        LIS_REAL z_imag;

        LIS_INT n;

        LIS_VECTOR x_real;
        LIS_VECTOR x_imag;

        LIS_VECTOR y_real;
        LIS_VECTOR y_imag;

        LIS_VECTOR work_real;
        LIS_VECTOR work_imag;
} LIS_EFEAST_SHIFTED_SOLVE_CONTEXT;


static void lis_efeast_shifted_solve_context_destroy(
        LIS_EFEAST_SHIFTED_SOLVE_CONTEXT *ctx)
{
        if( ctx==NULL )
        {
                return;
        }

        if( ctx->work_imag!=NULL )
        {
                lis_vector_destroy(
                        ctx->work_imag);
        }

        if( ctx->work_real!=NULL )
        {
                lis_vector_destroy(
                        ctx->work_real);
        }

        if( ctx->y_imag!=NULL )
        {
                lis_vector_destroy(
                        ctx->y_imag);
        }

        if( ctx->y_real!=NULL )
        {
                lis_vector_destroy(
                        ctx->y_real);
        }

        if( ctx->x_imag!=NULL )
        {
                lis_vector_destroy(
                        ctx->x_imag);
        }

        if( ctx->x_real!=NULL )
        {
                lis_vector_destroy(
                        ctx->x_real);
        }

        memset(
                ctx,
                0,
                sizeof(*ctx));
}


static LIS_INT lis_efeast_shifted_solve_context_create(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_REAL z_real,
        LIS_REAL z_imag,
        LIS_EFEAST_SHIFTED_SOLVE_CONTEXT *ctx)
{
        LIS_INT err;

        if( A==NULL ||
            ctx==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        memset(
                ctx,
                0,
                sizeof(*ctx));

        ctx->A = A;
        ctx->B = B;

        ctx->z_real = z_real;
        ctx->z_imag = z_imag;

        ctx->n = A->n;

        err = lis_vector_duplicate(
                A,
                &ctx->x_real);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &ctx->x_imag);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &ctx->y_real);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &ctx->y_imag);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &ctx->work_real);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &ctx->work_imag);

        if( err ) goto fail;

        return LIS_SUCCESS;

fail:

        lis_efeast_shifted_solve_context_destroy(
                ctx);

        return err;
}


static LIS_INT lis_efeast_shifted_solve_matvec(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        LIS_EFEAST_SHIFTED_SOLVE_CONTEXT *ctx;

        LIS_INT i;
        LIS_INT err;

        ctx =
                (LIS_EFEAST_SHIFTED_SOLVE_CONTEXT *)
                user_data;

        if( ctx==NULL ||
            ctx->A==NULL ||
            ctx->n<=0 ||
            x==NULL ||
            y==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        for(i=0;i<ctx->n;i++)
        {
                ctx->x_real->value[i] =
                        x[i];

                ctx->x_imag->value[i] =
                        x[ctx->n+i];
        }

        err = lis_efeast_shifted_matvec_split(
                ctx->A,
                ctx->B,
                ctx->z_real,
                ctx->z_imag,
                ctx->x_real,
                ctx->x_imag,
                ctx->y_real,
                ctx->y_imag,
                ctx->work_real,
                ctx->work_imag);

        if( err )
        {
                return err;
        }

        for(i=0;i<ctx->n;i++)
        {
                y[i] =
                        ctx->y_real->value[i];

                y[ctx->n+i] =
                        ctx->y_imag->value[i];
        }

        return LIS_SUCCESS;
}


static LIS_INT lis_efeast_shifted_vector_compatible(
        LIS_MATRIX A,
        LIS_VECTOR v)
{
        if( A==NULL ||
            v==NULL )
        {
                return LIS_FALSE;
        }

        if( v->label!=LIS_LABEL_VECTOR )
        {
                return LIS_FALSE;
        }

        if( v->precision!=LIS_PRECISION_DEFAULT )
        {
                return LIS_FALSE;
        }

        if( v->n!=A->n ||
            v->gn!=A->gn ||
            v->my_rank!=A->my_rank ||
            v->nprocs!=A->nprocs )
        {
                return LIS_FALSE;
        }

        return LIS_TRUE;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast_shifted_solve_split"

LIS_INT lis_efeast_shifted_solve_split(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_REAL z_real,
        LIS_REAL z_imag,
        LIS_VECTOR rhs_real,
        LIS_VECTOR rhs_imag,
        LIS_VECTOR x_real,
        LIS_VECTOR x_imag,
        LIS_SOLVER solver)
{
        LIS_EFEAST_SHIFTED_SOLVE_CONTEXT ctx;

        LIS_MATRIX shifted;

        LIS_VECTOR packed_rhs;
        LIS_VECTOR packed_x;

        LIS_INT nsolver;
        LIS_INT i;
        LIS_INT err;
        LIS_INT solve_err;

        LIS_DEBUG_FUNC_IN;

        memset(
                &ctx,
                0,
                sizeof(ctx));

        shifted = NULL;

        packed_rhs = NULL;
        packed_x = NULL;

        if( A==NULL ||
            rhs_real==NULL ||
            rhs_imag==NULL ||
            x_real==NULL ||
            x_imag==NULL ||
            solver==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted solve received an undefined argument\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( A->n<=0 ||
            A->gn<=0 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted solve requires a non-empty matrix\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !isfinite((double)z_real) ||
            !isfinite((double)z_imag) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted solve requires a finite complex shift\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( B!=NULL &&
            (B->n!=A->n ||
             B->gn!=A->gn ||
             B->my_rank!=A->my_rank ||
             B->nprocs!=A->nprocs) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted solve requires compatible A and B layouts\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !lis_efeast_shifted_vector_compatible(
                    A,
                    rhs_real) ||
            !lis_efeast_shifted_vector_compatible(
                    A,
                    rhs_imag) ||
            !lis_efeast_shifted_vector_compatible(
                    A,
                    x_real) ||
            !lis_efeast_shifted_vector_compatible(
                    A,
                    x_imag) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST shifted solve requires compatible split vectors\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        /*
         * Stage F3-2 intentionally uses the matrix-free path only.
         *
         * No preconditioner, scaling, transpose operator or requested
         * storage conversion may force materialization of the 2N shell.
         */
        if( solver->options[LIS_OPTIONS_PRECON]
                    != LIS_PRECON_TYPE_NONE ||
            solver->options[LIS_OPTIONS_SCALE]
                    != LIS_SCALE_NONE ||
            solver->options[LIS_OPTIONS_USE_AT]
                    != LIS_FALSE ||
            solver->options[LIS_OPTIONS_STORAGE]
                    != 0 ||
            solver->options[LIS_OPTIONS_PRECISION]
                    != LIS_PRECISION_DOUBLE )
        {
                LIS_SETERR(
                        LIS_ERR_NOT_IMPLEMENTED,
                        "FEAST shifted shell solve requires "
                        "double precision, no scaling, no preconditioner "
                        "and no storage conversion\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_NOT_IMPLEMENTED;
        }

        nsolver =
                solver->options[
                        LIS_OPTIONS_SOLVER];

        if( nsolver!=LIS_SOLVER_GMRES &&
            nsolver!=LIS_SOLVER_FGMRES &&
            nsolver!=LIS_SOLVER_GCRODR &&
            nsolver!=LIS_SOLVER_FGCRODR )
        {
                LIS_SETERR(
                        LIS_ERR_NOT_IMPLEMENTED,
                        "FEAST shifted shell solve requires "
                        "GMRES, FGMRES, GCRO-DR or FGCRO-DR\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_NOT_IMPLEMENTED;
        }

        err =
                lis_efeast_shifted_solve_context_create(
                        A,
                        B,
                        z_real,
                        z_imag,
                        &ctx);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_matrix_create(
                A->comm,
                &shifted);

        if( err ) goto cleanup;

        err = lis_matrix_set_size(
                shifted,
                2*A->n,
                2*A->gn);

        if( err ) goto cleanup;

        err = lis_matrix_set_user(
                shifted,
                &ctx,
                lis_efeast_shifted_solve_matvec,
                NULL);

        if( err ) goto cleanup;

        err = lis_vector_duplicate(
                shifted,
                &packed_rhs);

        if( err ) goto cleanup;

        err = lis_vector_duplicate(
                shifted,
                &packed_x);

        if( err ) goto cleanup;

        for(i=0;i<A->n;i++)
        {
                packed_rhs->value[i] =
                        rhs_real->value[i];

                packed_rhs->value[A->n+i] =
                        rhs_imag->value[i];

                packed_x->value[i] =
                        x_real->value[i];

                packed_x->value[A->n+i] =
                        x_imag->value[i];
        }

        solve_err = lis_solve(
                shifted,
                packed_rhs,
                packed_x,
                solver);

        /*
         * lis_solve() keeps diagnostics and reusable work in solver,
         * but the temporary FEAST shell and packed vectors are owned
         * by this call. Remove borrowed pointers before destruction.
         */
        solver->A = NULL;
        solver->Ah = NULL;
        solver->b = NULL;
        solver->x = NULL;
        solver->xx = NULL;
        solver->d = NULL;
        solver->setup = LIS_FALSE;

        if( solve_err )
        {
                err = solve_err;
                goto cleanup;
        }

        for(i=0;i<A->n;i++)
        {
                x_real->value[i] =
                        packed_x->value[i];

                x_imag->value[i] =
                        packed_x->value[A->n+i];
        }

        err = LIS_SUCCESS;

cleanup:

        if( packed_x!=NULL )
        {
                lis_vector_destroy(
                        packed_x);
        }

        if( packed_rhs!=NULL )
        {
                lis_vector_destroy(
                        packed_rhs);
        }

        if( shifted!=NULL )
        {
                lis_matrix_destroy(
                        shifted);
        }

        lis_efeast_shifted_solve_context_destroy(
                &ctx);

        LIS_DEBUG_FUNC_OUT;
        return err;
}


/*
 * Apply the FEAST spectral projector to one split-complex vector.
 *
 * The generalized projector is
 *
 *     Q = (1 / (2*pi*i)) integral
 *             (z B - A)^(-1) B X dz .
 *
 * For the standard problem, B is NULL and the right-hand side is X.
 *
 * lis_efeast_contour_gauss() provides quadrature points and weights
 * for the upper half contour.  This routine evaluates both the upper
 * point z and the conjugate lower point explicitly.  Consequently the
 * implementation does not require a real matrix, a real right-hand
 * side, or the relation
 *
 *     Y(conj(z)) = conj(Y(z)).
 *
 * If
 *
 *     w = wr + i*wi
 *     Y = Yr + i*Yi,
 *
 * then
 *
 *     Re(wY) = wr*Yr - wi*Yi
 *     Im(wY) = wi*Yr + wr*Yi.
 *
 * The lower contour uses conj(w).
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_projector_apply_split"

LIS_INT lis_efeast_projector_apply_split(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_REAL emin,
        LIS_REAL emax,
        LIS_REAL ratio,
        LIS_INT npoints,
        LIS_VECTOR x_real,
        LIS_VECTOR x_imag,
        LIS_VECTOR q_real,
        LIS_VECTOR q_imag,
        LIS_SOLVER solver)
{
        LIS_REAL z_real[1024];
        LIS_REAL z_imag[1024];
        LIS_REAL w_real[1024];
        LIS_REAL w_imag[1024];

        LIS_VECTOR rhs_real;
        LIS_VECTOR rhs_imag;

        LIS_VECTOR y_real;
        LIS_VECTOR y_imag;

        LIS_INT i;
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        rhs_real = NULL;
        rhs_imag = NULL;

        y_real = NULL;
        y_imag = NULL;

        if( A==NULL ||
            x_real==NULL ||
            x_imag==NULL ||
            q_real==NULL ||
            q_imag==NULL ||
            solver==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projector received an undefined argument\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !lis_efeast_shifted_vector_compatible(
                    A,
                    x_real) ||
            !lis_efeast_shifted_vector_compatible(
                    A,
                    x_imag) ||
            !lis_efeast_shifted_vector_compatible(
                    A,
                    q_real) ||
            !lis_efeast_shifted_vector_compatible(
                    A,
                    q_imag) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projector requires compatible split vectors\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( B!=NULL &&
            (B->n!=A->n ||
             B->gn!=A->gn ||
             B->my_rank!=A->my_rank ||
             B->nprocs!=A->nprocs) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projector requires compatible A and B layouts\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        err = lis_efeast_contour_gauss(
                emin,
                emax,
                npoints,
                ratio,
                z_real,
                z_imag,
                w_real,
                w_imag);

        if( err )
        {
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

        err = lis_vector_duplicate(
                A,
                &rhs_real);

        if( err ) goto cleanup;

        err = lis_vector_duplicate(
                A,
                &rhs_imag);

        if( err ) goto cleanup;

        err = lis_vector_duplicate(
                A,
                &y_real);

        if( err ) goto cleanup;

        err = lis_vector_duplicate(
                A,
                &y_imag);

        if( err ) goto cleanup;

        /*
         * The generalized FEAST projector acts on B*X.
         */
        if( B!=NULL )
        {
                err = lis_matvec(
                        B,
                        x_real,
                        rhs_real);

                if( err ) goto cleanup;

                err = lis_matvec(
                        B,
                        x_imag,
                        rhs_imag);

                if( err ) goto cleanup;
        }
        else
        {
                err = lis_vector_copy(
                        x_real,
                        rhs_real);

                if( err ) goto cleanup;

                err = lis_vector_copy(
                        x_imag,
                        rhs_imag);

                if( err ) goto cleanup;
        }

        /*
         * The RHS has now been snapshotted, so output vectors may
         * safely alias the input vectors.
         */
        err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                q_real);

        if( err ) goto cleanup;

        err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                q_imag);

        if( err ) goto cleanup;

        for(i=0;i<npoints;i++)
        {
                /*
                 * Upper half:
                 *
                 *     w * Y(z)
                 */
                err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        y_real);

                if( err ) goto cleanup;

                err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        y_imag);

                if( err ) goto cleanup;

                err = lis_efeast_shifted_solve_split(
                        A,
                        B,
                        z_real[i],
                        z_imag[i],
                        rhs_real,
                        rhs_imag,
                        y_real,
                        y_imag,
                        solver);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)w_real[i],
                        y_real,
                        q_real);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)(-w_imag[i]),
                        y_imag,
                        q_real);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)w_imag[i],
                        y_real,
                        q_imag);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)w_real[i],
                        y_imag,
                        q_imag);

                if( err ) goto cleanup;

                /*
                 * Lower half:
                 *
                 *     conj(w) * Y(conj(z))
                 *
                 * Evaluate the solve explicitly instead of deriving it
                 * by conjugation.  This keeps the projector valid for
                 * complex A, B and X.
                 */
                err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        y_real);

                if( err ) goto cleanup;

                err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        y_imag);

                if( err ) goto cleanup;

                err = lis_efeast_shifted_solve_split(
                        A,
                        B,
                        z_real[i],
                        -z_imag[i],
                        rhs_real,
                        rhs_imag,
                        y_real,
                        y_imag,
                        solver);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)w_real[i],
                        y_real,
                        q_real);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)w_imag[i],
                        y_imag,
                        q_real);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)(-w_imag[i]),
                        y_real,
                        q_imag);

                if( err ) goto cleanup;

                err = lis_vector_axpy(
                        (LIS_SCALAR)w_real[i],
                        y_imag,
                        q_imag);

                if( err ) goto cleanup;
        }

        err = LIS_SUCCESS;

cleanup:

        if( y_imag!=NULL )
        {
                lis_vector_destroy(
                        y_imag);
        }

        if( y_real!=NULL )
        {
                lis_vector_destroy(
                        y_real);
        }

        if( rhs_imag!=NULL )
        {
                lis_vector_destroy(
                        rhs_imag);
        }

        if( rhs_real!=NULL )
        {
                lis_vector_destroy(
                        rhs_real);
        }

        LIS_DEBUG_FUNC_OUT;
        return err;
}


/*
 * Destroy an internally owned FEAST vector block.
 */
static void lis_efeast_vector_block_destroy(
        LIS_INT count,
        LIS_VECTOR block[])
{
        LIS_INT i;

        if( block==NULL )
        {
                return;
        }

        for(i=0;i<count;i++)
        {
                if( block[i]!=NULL )
                {
                        lis_vector_destroy(
                                block[i]);
                }
        }

        lis_free(
                block);
}


/*
 * Build a rank-revealing filtered FEAST block basis.
 *
 * Each input vector is independently passed through the contour
 * spectral projector.  The resulting raw block W is then compressed
 * and normalized by the common block rank-deflation core.
 *
 * Standard problem:
 *
 *     X -> P X -> W -> rank_deflate -> Q
 *
 * Generalized problem:
 *
 *     X -> P_B X -> W -> B-rank_deflate -> Q, BQ
 *
 * The rank-deflation helpers already perform normalization through
 * the Gram eigensystem.  No second Gram-Schmidt pass is required.
 *
 * In a complex build the split projector result is reconstructed as
 *
 *     W = q_real + i*q_imag.
 *
 * In a real build the full contour is conjugate symmetric and the
 * physical filtered basis is q_real.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_filtered_block_basis"

LIS_INT lis_efeast_filtered_block_basis(
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_REAL emin,
        LIS_REAL emax,
        LIS_REAL ratio,
        LIS_INT npoints,
        LIS_INT block_size,
        LIS_VECTOR X[],
        LIS_VECTOR Q[],
        LIS_VECTOR BQ[],
        LIS_REAL rank_tol,
        LIS_INT *active_rank,
        LIS_SOLVER solver)
{
        LIS_VECTOR *W;
        LIS_VECTOR *BW;

        LIS_VECTOR zero_imag;
        LIS_VECTOR q_real;
        LIS_VECTOR q_imag;

        LIS_SCALAR *gram;
        LIS_SCALAR *eval;
        LIS_SCALAR *evec;
        LIS_SCALAR *transform;

        size_t vector_bytes;
        size_t square_count;
        size_t square_bytes;
        size_t eval_bytes;

        LIS_INT i;
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        W = NULL;
        BW = NULL;

        zero_imag = NULL;
        q_real = NULL;
        q_imag = NULL;

        gram = NULL;
        eval = NULL;
        evec = NULL;
        transform = NULL;

        if( active_rank!=NULL )
        {
                *active_rank = 0;
        }

        if( A==NULL ||
            block_size<=0 ||
            X==NULL ||
            Q==NULL ||
            active_rank==NULL ||
            solver==NULL ||
            !isfinite((double)rank_tol) ||
            rank_tol<0.0 ||
            rank_tol>=1.0 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST filtered block received invalid arguments\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( B!=NULL && BQ==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Generalized FEAST filtered block requires BQ outputs\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( B==NULL && BQ!=NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Standard FEAST filtered block requires BQ to be NULL\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( B!=NULL &&
            (B->n!=A->n ||
             B->gn!=A->gn ||
             B->my_rank!=A->my_rank ||
             B->nprocs!=A->nprocs) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST filtered block requires compatible A and B layouts\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        for(i=0;i<block_size;i++)
        {
                if( X[i]==NULL ||
                    Q[i]==NULL ||
                    !lis_efeast_shifted_vector_compatible(
                            A,
                            X[i]) ||
                    !lis_efeast_shifted_vector_compatible(
                            A,
                            Q[i]) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST filtered block contains an incompatible vector\n");

                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                if( B!=NULL &&
                    (BQ[i]==NULL ||
                     !lis_efeast_shifted_vector_compatible(
                            A,
                            BQ[i])) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST filtered block contains an incompatible BQ vector\n");

                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }
        }

        vector_bytes =
                (size_t)block_size
                *
                sizeof(LIS_VECTOR);

        square_count =
                (size_t)block_size
                *
                (size_t)block_size;

        square_bytes =
                square_count
                *
                sizeof(LIS_SCALAR);

        eval_bytes =
                (size_t)block_size
                *
                sizeof(LIS_SCALAR);

        W = (LIS_VECTOR *)lis_malloc(
                vector_bytes,
                "lis_efeast_filtered_block_basis::W");

        if( W==NULL )
        {
                LIS_SETERR_MEM(
                        vector_bytes);

                LIS_DEBUG_FUNC_OUT;
                return LIS_OUT_OF_MEMORY;
        }

        memset(
                W,
                0,
                vector_bytes);

        if( B!=NULL )
        {
                BW = (LIS_VECTOR *)lis_malloc(
                        vector_bytes,
                        "lis_efeast_filtered_block_basis::BW");

                if( BW==NULL )
                {
                        lis_free(
                                W);

                        LIS_SETERR_MEM(
                                vector_bytes);

                        LIS_DEBUG_FUNC_OUT;
                        return LIS_OUT_OF_MEMORY;
                }

                memset(
                        BW,
                        0,
                        vector_bytes);
        }

        gram = (LIS_SCALAR *)lis_malloc(
                square_bytes,
                "lis_efeast_filtered_block_basis::gram");

        if( gram==NULL )
        {
                err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        square_bytes);
                goto cleanup;
        }

        eval = (LIS_SCALAR *)lis_malloc(
                eval_bytes,
                "lis_efeast_filtered_block_basis::eval");

        if( eval==NULL )
        {
                err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        eval_bytes);
                goto cleanup;
        }

        evec = (LIS_SCALAR *)lis_malloc(
                square_bytes,
                "lis_efeast_filtered_block_basis::evec");

        if( evec==NULL )
        {
                err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        square_bytes);
                goto cleanup;
        }

        transform = (LIS_SCALAR *)lis_malloc(
                square_bytes,
                "lis_efeast_filtered_block_basis::transform");

        if( transform==NULL )
        {
                err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        square_bytes);
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &zero_imag);

        if( err ) goto cleanup;

        err = lis_vector_duplicate(
                A,
                &q_real);

        if( err ) goto cleanup;

        err = lis_vector_duplicate(
                A,
                &q_imag);

        if( err ) goto cleanup;

        err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                zero_imag);

        if( err ) goto cleanup;

        for(i=0;i<block_size;i++)
        {
                err = lis_vector_duplicate(
                        A,
                        &W[i]);

                if( err ) goto cleanup;

                if( B!=NULL )
                {
                        err = lis_vector_duplicate(
                                A,
                                &BW[i]);

                        if( err ) goto cleanup;
                }
        }

        /*
         * Filter each candidate independently.
         *
         * x_imag=0 is sufficient even for a complex X because the
         * split representation is linear:
         *
         *     X = X + i*0.
         */
        for(i=0;i<block_size;i++)
        {
                err = lis_efeast_projector_apply_split(
                        A,
                        B,
                        emin,
                        emax,
                        ratio,
                        npoints,
                        X[i],
                        zero_imag,
                        q_real,
                        q_imag,
                        solver);

                if( err ) goto cleanup;

                err = lis_vector_copy(
                        q_real,
                        W[i]);

                if( err ) goto cleanup;

#ifdef _COMPLEX

                err = lis_vector_axpy(
                        (LIS_SCALAR)(1.0*I),
                        q_imag,
                        W[i]);

                if( err ) goto cleanup;

#endif
        }

        /*
         * Rank deflation is also the normalization step.
         *
         * The generalized helper additionally returns BQ and reuses
         * BW internally, avoiding a separate post-normalization B*Q.
         */
        if( B!=NULL )
        {
                err = lis_eigen_block_brank_deflate(
                        B,
                        block_size,
                        block_size,
                        W,
                        BW,
                        Q,
                        BQ,
                        rank_tol,
                        active_rank,
                        gram,
                        eval,
                        evec,
                        transform);
        }
        else
        {
                err = lis_eigen_block_rank_deflate(
                        block_size,
                        block_size,
                        W,
                        Q,
                        rank_tol,
                        active_rank,
                        gram,
                        eval,
                        evec,
                        transform);
        }

cleanup:

        if( q_imag!=NULL )
        {
                lis_vector_destroy(
                        q_imag);
        }

        if( q_real!=NULL )
        {
                lis_vector_destroy(
                        q_real);
        }

        if( zero_imag!=NULL )
        {
                lis_vector_destroy(
                        zero_imag);
        }

        lis_efeast_vector_block_destroy(
                block_size,
                BW);

        lis_efeast_vector_block_destroy(
                block_size,
                W);

        if( transform!=NULL )
        {
                lis_free(
                        transform);
        }

        if( evec!=NULL )
        {
                lis_free(
                        evec);
        }

        if( eval!=NULL )
        {
                lis_free(
                        eval);
        }

        if( gram!=NULL )
        {
                lis_free(
                        gram);
        }

        LIS_DEBUG_FUNC_OUT;
        return err;
}


/*
 * Compute the Rayleigh-Ritz eigensystem on an already normalized
 * FEAST filtered basis.
 *
 * For a standard problem Q is Euclidean-orthonormal.
 *
 * For a generalized problem F3-4 supplies a B-orthonormal basis:
 *
 *     Q^H B Q = I.
 *
 * Therefore both cases reduce to the ordinary projected Hermitian
 * eigensystem
 *
 *     H C = C Lambda,
 *
 * where
 *
 *     H = Q^H A Q.
 *
 * The physical Ritz vectors are reconstructed as
 *
 *     X = Q C.
 *
 * No projected B matrix is required here.
 */
#undef __FUNC__
#define __FUNC__ "lis_efeast_projected_rayleigh_ritz"

LIS_INT lis_efeast_projected_rayleigh_ritz(
        LIS_MATRIX A,
        LIS_INT active_rank,
        LIS_VECTOR Q[],
        LIS_SCALAR ritz_value[],
        LIS_VECTOR ritz_vector[],
        LIS_INT *projected_niter,
        LIS_REAL *projected_error)
{
        LIS_VECTOR *AQ;

        LIS_SCALAR *projected;
        LIS_SCALAR *coeff;

        LIS_SCALAR value;

        size_t vector_bytes;
        size_t square_count;
        size_t square_bytes;

        LIS_INT i;
        LIS_INT j;
        LIS_INT k;
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        AQ = NULL;

        projected = NULL;
        coeff = NULL;

        if( projected_niter!=NULL )
        {
                *projected_niter = 0;
        }

        if( projected_error!=NULL )
        {
                *projected_error = 0.0;
        }

        if( A==NULL ||
            active_rank<=0 ||
            Q==NULL ||
            ritz_value==NULL ||
            ritz_vector==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projected Rayleigh-Ritz received "
                        "invalid arguments\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        for(i=0;i<active_rank;i++)
        {
                if( Q[i]==NULL ||
                    ritz_vector[i]==NULL ||
                    !lis_efeast_shifted_vector_compatible(
                            A,
                            Q[i]) ||
                    !lis_efeast_shifted_vector_compatible(
                            A,
                            ritz_vector[i]) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST projected Rayleigh-Ritz contains "
                                "an incompatible vector\n");

                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                for(j=0;j<active_rank;j++)
                {
                        if( ritz_vector[i]==Q[j] )
                        {
                                LIS_SETERR(
                                        LIS_ERR_ILL_ARG,
                                        "FEAST Ritz reconstruction requires "
                                        "out-of-place vectors\n");

                                LIS_DEBUG_FUNC_OUT;
                                return LIS_ERR_ILL_ARG;
                        }
                }

                for(j=0;j<i;j++)
                {
                        if( ritz_vector[i]==ritz_vector[j] )
                        {
                                LIS_SETERR(
                                        LIS_ERR_ILL_ARG,
                                        "FEAST Ritz vectors must be distinct\n");

                                LIS_DEBUG_FUNC_OUT;
                                return LIS_ERR_ILL_ARG;
                        }
                }
        }

        vector_bytes =
                (size_t)active_rank
                *
                sizeof(LIS_VECTOR);

        square_count =
                (size_t)active_rank
                *
                (size_t)active_rank;

        square_bytes =
                square_count
                *
                sizeof(LIS_SCALAR);

        AQ = (LIS_VECTOR *)lis_malloc(
                vector_bytes,
                "lis_efeast_projected_rayleigh_ritz::AQ");

        if( AQ==NULL )
        {
                LIS_SETERR_MEM(
                        vector_bytes);

                LIS_DEBUG_FUNC_OUT;
                return LIS_OUT_OF_MEMORY;
        }

        memset(
                AQ,
                0,
                vector_bytes);

        projected = (LIS_SCALAR *)lis_malloc(
                square_bytes,
                "lis_efeast_projected_rayleigh_ritz::projected");

        if( projected==NULL )
        {
                err = LIS_OUT_OF_MEMORY;

                LIS_SETERR_MEM(
                        square_bytes);

                goto cleanup;
        }

        coeff = (LIS_SCALAR *)lis_malloc(
                square_bytes,
                "lis_efeast_projected_rayleigh_ritz::coeff");

        if( coeff==NULL )
        {
                err = LIS_OUT_OF_MEMORY;

                LIS_SETERR_MEM(
                        square_bytes);

                goto cleanup;
        }

        for(i=0;i<active_rank;i++)
        {
                err = lis_vector_duplicate(
                        A,
                        &AQ[i]);

                if( err ) goto cleanup;
        }

        /*
         * Use the block matvec path when the matrix supplies one.
         * Ordinary matrix formats retain the scalar fallback.
         */
        err = lis_matvec_block(
                A,
                active_rank,
                Q,
                AQ);

        if( err ) goto cleanup;

        /*
         * H = Q^H A Q.
         */
        err = lis_eigen_block_gram(
                active_rank,
                Q,
                active_rank,
                AQ,
                projected);

        if( err ) goto cleanup;

        /*
         * Remove only roundoff-level non-Hermitian skew introduced
         * by distributed or threaded dot-product reductions.
         *
         * The projected eigensolver intentionally accepts a
         * symmetric/Hermitian matrix, so enforce that contract
         * explicitly before entering the dense solve.
         */
        for(i=0;i<active_rank;i++)
        {
#ifdef _COMPLEX

                projected[i*active_rank+i] =
                        (LIS_SCALAR)creal(
                                projected[i*active_rank+i]);

#endif

                for(j=i+1;j<active_rank;j++)
                {
#ifdef _COMPLEX

                        value =
                                (LIS_SCALAR)0.5
                                *
                                (
                                        projected[i*active_rank+j]
                                        +
                                        conj(
                                                projected[
                                                        j*active_rank+i])
                                );

                        projected[i*active_rank+j] =
                                value;

                        projected[j*active_rank+i] =
                                conj(
                                        value);

#else

                        value =
                                (LIS_SCALAR)0.5
                                *
                                (
                                        projected[i*active_rank+j]
                                        +
                                        projected[j*active_rank+i]
                                );

                        projected[i*active_rank+j] =
                                value;

                        projected[j*active_rank+i] =
                                value;

#endif
                }
        }

        /*
         * The common dense helper overwrites the projected matrix
         * and returns projected eigenvectors as columns of coeff.
         */
        err = lis_eigen_projected_solve(
                active_rank,
                projected,
                ritz_value,
                coeff,
                projected_niter,
                projected_error);

        if( err ) goto cleanup;

        /*
         * ritz_vector[j] =
         *
         *     sum_i Q[i] * coeff[i*active_rank+j].
         */
        err = lis_eigen_ritz_reconstruct(
                active_rank,
                Q,
                active_rank,
                coeff,
                active_rank,
                ritz_vector);

        if( err ) goto cleanup;

        /*
         * Reject non-finite projected eigenvalues before returning
         * them to the outer FEAST iteration.
         */
        for(k=0;k<active_rank;k++)
        {
#ifdef _COMPLEX

                if( !isfinite((double)creal(ritz_value[k])) ||
                    !isfinite((double)cimag(ritz_value[k])) )

#else

                if( !isfinite((double)ritz_value[k]) )

#endif
                {
                        LIS_SETERR(
                                LIS_BREAKDOWN,
                                "FEAST projected Rayleigh-Ritz produced "
                                "a non-finite eigenvalue\n");

                        err = LIS_BREAKDOWN;
                        goto cleanup;
                }
        }

        err = LIS_SUCCESS;

cleanup:

        lis_efeast_vector_block_destroy(
                active_rank,
                AQ);

        if( coeff!=NULL )
        {
                lis_free(
                        coeff);
        }

        if( projected!=NULL )
        {
                lis_free(
                        projected);
        }

        LIS_DEBUG_FUNC_OUT;
        return err;
}


static LIS_INT lis_efeast_parse_real(
        const char *text,
        LIS_REAL *value)
{
        char *endptr;
        double parsed;

        if( text==NULL || value==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        endptr = NULL;
        parsed = strtod(text,&endptr);

        if( endptr==text ||
            endptr==NULL ||
            *endptr!='\0' ||
            !isfinite(parsed) )
        {
                return LIS_ERR_ILL_ARG;
        }

        *value = (LIS_REAL)parsed;

        return LIS_SUCCESS;
}


static LIS_INT lis_efeast_parse_int(
        const char *text,
        LIS_INT *value)
{
        char *endptr;
        long parsed;

        if( text==NULL || value==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        endptr = NULL;
        parsed = strtol(text,&endptr,10);

        if( endptr==text ||
            endptr==NULL ||
            *endptr!='\0' ||
            parsed<1 ||
            parsed>1024 )
        {
                return LIS_ERR_ILL_ARG;
        }

        *value = (LIS_INT)parsed;

        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast_set_option"

LIS_INT lis_efeast_set_option(
        char *arg1,
        char *arg2,
        LIS_ESOLVER esolver)
{
        LIS_REAL value;
        LIS_INT ivalue;
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL ||
            arg1==NULL ||
            arg2==NULL )
        {
                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( strcmp(arg1,"-feast_emin")==0 )
        {
                err = lis_efeast_parse_real(arg2,&value);

                if( err )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_emin must be a finite real value\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_emin = value;
                esolver->feast_interval_mask |=
                        LIS_FEAST_INTERVAL_EMIN;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( strcmp(arg1,"-feast_emax")==0 )
        {
                err = lis_efeast_parse_real(arg2,&value);

                if( err )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_emax must be a finite real value\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_emax = value;
                esolver->feast_interval_mask |=
                        LIS_FEAST_INTERVAL_EMAX;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( strcmp(arg1,"-feast_ratio")==0 )
        {
                err = lis_efeast_parse_real(arg2,&value);

                if( err ||
                    !(value>0.0) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_ratio must be finite and positive\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_ratio = value;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( strcmp(arg1,"-feast_points")==0 )
        {
                err = lis_efeast_parse_int(arg2,&ivalue);

                if( err )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "Parameter -feast_points must be an integer in [1,1024]\n");
                        LIS_DEBUG_FUNC_OUT;
                        return LIS_ERR_ILL_ARG;
                }

                esolver->feast_points = ivalue;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        LIS_SETERR1(
                LIS_ERR_ILL_ARG,
                "Unknown native FEAST option %s\n",
                arg1);

        LIS_DEBUG_FUNC_OUT;
        return LIS_ERR_ILL_ARG;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast_check_params"

LIS_INT lis_efeast_check_params(
        LIS_ESOLVER esolver)
{
        LIS_INT ss;

        LIS_DEBUG_FUNC_IN;

        if( esolver==NULL )
        {
                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

        if( ss<=0 )
        {
                LIS_SETERR1(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST requires positive subspace size -ss (got %D)\n",
                        ss);

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( esolver->feast_interval_mask !=
            (LIS_FEAST_INTERVAL_EMIN |
             LIS_FEAST_INTERVAL_EMAX) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST requires both -feast_emin and -feast_emax\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !(esolver->feast_emin <
              esolver->feast_emax) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST requires -feast_emin < -feast_emax\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( !(esolver->feast_ratio>0.0) ||
            !isfinite((double)esolver->feast_ratio) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST contour ratio must be finite and positive\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        if( esolver->feast_points<1 ||
            esolver->feast_points>1024 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST contour point count is invalid\n");

                LIS_DEBUG_FUNC_OUT;
                return LIS_ERR_ILL_ARG;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast_malloc_work"

LIS_INT lis_efeast_malloc_work(
        LIS_ESOLVER esolver)
{
        LIS_DEBUG_FUNC_IN;

        esolver->work = NULL;
        esolver->worklen = 0;

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_efeast"

LIS_INT lis_efeast(
        LIS_ESOLVER esolver)
{
        LIS_INT err;

        LIS_DEBUG_FUNC_IN;

        err = lis_efeast_check_params(esolver);

        if( err )
        {
                esolver->retcode = err;
                LIS_DEBUG_FUNC_OUT;
                return err;
        }

#ifndef _COMPLEX

        LIS_SETERR(
                LIS_ERR_NOT_IMPLEMENTED,
                "Native FEAST requires complex shifted systems; "
                "the real-build bridge is not implemented yet\n");

#else

        LIS_SETERR(
                LIS_ERR_NOT_IMPLEMENTED,
                "Native FEAST contour projector is not implemented in stage F0\n");

#endif

        esolver->retcode = LIS_ERR_NOT_IMPLEMENTED;

        LIS_DEBUG_FUNC_OUT;
        return LIS_ERR_NOT_IMPLEMENTED;
}
