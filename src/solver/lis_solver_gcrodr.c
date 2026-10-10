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

/*
 * GCRO-DR implementation for LIS.
 *
 * This is an original LIS-oriented implementation based on the published
 * GCRO-DR algorithm (Parks et al., 2006).
 *
 * Scope:
 *   - real/default precision
 *   - single RHS
 *   - right-preconditioned operator A M^{-1}
 *   - persistent recycle space across repeated lis_solve() calls
 *   - configurable recycle dimension and quality thresholds
 *   - cross-solve recycle quality measurement
 *   - KEEP / REFRESH / REJECT policy for changing operators/preconditioners
 *   - harmonic Ritz selection from the augmented GCRO relation
 *   - real representation of complex-conjugate Ritz pairs
 *
 * The recycle vectors solver->recycle_u live in the domain of the
 * right-preconditioned operator.  For each solve we refresh
 *
 *     ZU = M^{-1} U,   C = A ZU,
 *
 * and orthonormalize C while applying the same transformations to U and ZU.
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
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "lislib.h"

#define GCRODR_NWORK 4
#define GCRODR_TINY ((LIS_REAL)1.0e-14)
#define GCRODR_REORTH_TINY ((LIS_REAL)1.0e-12)
#define GCRODR_DEFAULT_RECYCLE_DIM ((LIS_INT)8)

#ifndef _COMPLEX

typedef struct
{
    LIS_REAL r;
    LIS_REAL i;
} LIS_GCRODR_CPX;

typedef struct
{
    LIS_REAL lr;
    LIS_REAL li;
    LIS_REAL score;
    LIS_INT pair;
} LIS_GCRODR_EIG;

static LIS_GCRODR_CPX lis_gcrodr_cpx(LIS_REAL r, LIS_REAL i)
{
    LIS_GCRODR_CPX z;
    z.r = r;
    z.i = i;
    return z;
}

static LIS_GCRODR_CPX lis_gcrodr_cadd(LIS_GCRODR_CPX a, LIS_GCRODR_CPX b)
{
    return lis_gcrodr_cpx(a.r+b.r,a.i+b.i);
}

static LIS_GCRODR_CPX lis_gcrodr_csub(LIS_GCRODR_CPX a, LIS_GCRODR_CPX b)
{
    return lis_gcrodr_cpx(a.r-b.r,a.i-b.i);
}

static LIS_GCRODR_CPX lis_gcrodr_cmul(LIS_GCRODR_CPX a, LIS_GCRODR_CPX b)
{
    return lis_gcrodr_cpx(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);
}

static LIS_REAL lis_gcrodr_cabs2(LIS_GCRODR_CPX a)
{
    return a.r*a.r+a.i*a.i;
}

static LIS_GCRODR_CPX lis_gcrodr_cdiv(LIS_GCRODR_CPX a, LIS_GCRODR_CPX b)
{
    LIS_REAL d = lis_gcrodr_cabs2(b);
    if( d<=GCRODR_TINY*GCRODR_TINY ) return lis_gcrodr_cpx(0.0,0.0);
    return lis_gcrodr_cpx((a.r*b.r+a.i*b.i)/d,
                          (a.i*b.r-a.r*b.i)/d);
}

/* Dense complex Gaussian elimination with partial pivoting.
 * Matrix is row-major and is overwritten.
 */
static LIS_INT lis_gcrodr_csolve(LIS_INT n, LIS_GCRODR_CPX *a,
                                 LIS_GCRODR_CPX *b, LIS_GCRODR_CPX *x)
{
    LIS_INT i,j,k,p;
    LIS_REAL vmax,v;
    LIS_GCRODR_CPX tmp,factor,sum;

    for(k=0;k<n;k++)
    {
        p = k;
        vmax = lis_gcrodr_cabs2(a[k*n+k]);
        for(i=k+1;i<n;i++)
        {
            v = lis_gcrodr_cabs2(a[i*n+k]);
            if( v>vmax )
            {
                vmax = v;
                p = i;
            }
        }
        if( vmax<=GCRODR_TINY*GCRODR_TINY ) return LIS_BREAKDOWN;

        if( p!=k )
        {
            for(j=k;j<n;j++)
            {
                tmp = a[k*n+j];
                a[k*n+j] = a[p*n+j];
                a[p*n+j] = tmp;
            }
            tmp = b[k]; b[k] = b[p]; b[p] = tmp;
        }

        for(i=k+1;i<n;i++)
        {
            factor = lis_gcrodr_cdiv(a[i*n+k],a[k*n+k]);
            a[i*n+k] = lis_gcrodr_cpx(0.0,0.0);
            for(j=k+1;j<n;j++)
            {
                a[i*n+j] = lis_gcrodr_csub(a[i*n+j],
                                           lis_gcrodr_cmul(factor,a[k*n+j]));
            }
            b[i] = lis_gcrodr_csub(b[i],lis_gcrodr_cmul(factor,b[k]));
        }
    }

    for(i=n-1;i>=0;i--)
    {
        sum = b[i];
        for(j=i+1;j<n;j++)
        {
            sum = lis_gcrodr_csub(sum,lis_gcrodr_cmul(a[i*n+j],x[j]));
        }
        if( lis_gcrodr_cabs2(a[i*n+i])<=GCRODR_TINY*GCRODR_TINY )
            return LIS_BREAKDOWN;
        x[i] = lis_gcrodr_cdiv(sum,a[i*n+i]);
    }
    return LIS_SUCCESS;
}

/* Inverse iteration for a (possibly complex) eigenvalue of a real dense T.
 * T uses LIS column-major storage.  A tiny deterministic shift makes the
 * inverse system nonsingular even when lambda is very accurate.
 */
static LIS_INT lis_gcrodr_inverse_eigenvector(
    LIS_INT n, const LIS_SCALAR *T,
    LIS_REAL lr, LIS_REAL li,
    LIS_REAL *vr, LIS_REAL *vi)
{
    LIS_GCRODR_CPX *a,*rhs,*x;
    LIS_GCRODR_CPX sigma;
    LIS_REAL delta,nrm,mag;
    LIS_INT i,j,it,err;

    a   = (LIS_GCRODR_CPX *)lis_malloc((size_t)n*n*sizeof(*a),
                                       "lis_gcrodr_inverse_eigenvector::a");
    rhs = (LIS_GCRODR_CPX *)lis_malloc((size_t)n*sizeof(*rhs),
                                       "lis_gcrodr_inverse_eigenvector::rhs");
    x   = (LIS_GCRODR_CPX *)lis_malloc((size_t)n*sizeof(*x),
                                       "lis_gcrodr_inverse_eigenvector::x");
    if( a==NULL || rhs==NULL || x==NULL )
    {
        lis_free(a); lis_free(rhs); lis_free(x);
        return LIS_ERR_OUT_OF_MEMORY;
    }

    mag = sqrt(lr*lr+li*li);
    delta = (LIS_REAL)1.0e-8*((LIS_REAL)1.0+mag);
    sigma.r = lr + delta;
    if( li==0.0 ) sigma.i = 0.0;
    else sigma.i = li + ((li>0.0)?(LIS_REAL)0.125*delta:-(LIS_REAL)0.125*delta);

    for(i=0;i<n;i++)
    {
        rhs[i].r = (LIS_REAL)1.0/(LIS_REAL)(i+1);
        rhs[i].i = (li!=0.0) ? (LIS_REAL)0.125/(LIS_REAL)(n-i) : 0.0;
    }

    err = LIS_SUCCESS;
    for(it=0;it<5;it++)
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                a[i*n+j].r = (LIS_REAL)T[i+j*n];
                a[i*n+j].i = 0.0;
            }
            a[i*n+i] = lis_gcrodr_csub(a[i*n+i],sigma);
            x[i] = lis_gcrodr_cpx(0.0,0.0);
        }

        err = lis_gcrodr_csolve(n,a,rhs,x);
        if( err ) break;

        nrm = 0.0;
        for(i=0;i<n;i++) nrm += lis_gcrodr_cabs2(x[i]);
        nrm = sqrt(nrm);
        if( nrm<=GCRODR_TINY || nrm!=nrm )
        {
            err = LIS_BREAKDOWN;
            break;
        }
        for(i=0;i<n;i++)
        {
            rhs[i].r = x[i].r/nrm;
            rhs[i].i = x[i].i/nrm;
        }
    }

    if( err==LIS_SUCCESS )
    {
        for(i=0;i<n;i++)
        {
            vr[i] = rhs[i].r;
            vi[i] = rhs[i].i;
        }
    }

    lis_free(a); lis_free(rhs); lis_free(x);
    return err;
}

static void lis_gcrodr_sort_eigs(LIS_GCRODR_EIG *e, LIS_INT n)
{
    LIS_INT i,j;
    LIS_GCRODR_EIG x;
    for(i=1;i<n;i++)
    {
        x = e[i];
        j = i-1;
        while( j>=0 && e[j].score < x.score )
        {
            e[j+1] = e[j];
            j--;
        }
        e[j+1] = x;
    }
}

static LIS_INT lis_gcrodr_extract_eigs(LIS_INT n, LIS_SCALAR *tqr,
                                       LIS_GCRODR_EIG *eig, LIS_INT *neig)
{
    LIS_INT i,m;
    LIS_REAL a,b,c,d,tr,det,disc,s,sub,scale;

    m = 0;
    i = 0;
    while( i<n )
    {
        scale = (LIS_REAL)1.0 + fabs((LIS_REAL)tqr[i+i*n]);
        sub = (i+1<n) ? fabs((LIS_REAL)tqr[(i+1)+i*n]) : 0.0;
        if( i+1<n && sub > (LIS_REAL)1.0e-10*scale )
        {
            a = (LIS_REAL)tqr[i+i*n];
            b = (LIS_REAL)tqr[i+(i+1)*n];
            c = (LIS_REAL)tqr[(i+1)+i*n];
            d = (LIS_REAL)tqr[(i+1)+(i+1)*n];
            tr = a+d;
            det = a*d-b*c;
            disc = tr*tr-(LIS_REAL)4.0*det;
            if( disc>=0.0 )
            {
                s = sqrt(disc);
                eig[m].lr = ((LIS_REAL)0.5)*(tr+s);
                eig[m].li = 0.0;
                eig[m].score = fabs(eig[m].lr);
                eig[m].pair = 0;
                m++;
                eig[m].lr = ((LIS_REAL)0.5)*(tr-s);
                eig[m].li = 0.0;
                eig[m].score = fabs(eig[m].lr);
                eig[m].pair = 0;
                m++;
            }
            else
            {
                eig[m].lr = ((LIS_REAL)0.5)*tr;
                eig[m].li = ((LIS_REAL)0.5)*sqrt(-disc);
                eig[m].score = sqrt(eig[m].lr*eig[m].lr +
                                    eig[m].li*eig[m].li);
                eig[m].pair = 1;
                m++;
            }
            i += 2;
        }
        else
        {
            eig[m].lr = (LIS_REAL)tqr[i+i*n];
            eig[m].li = 0.0;
            eig[m].score = fabs(eig[m].lr);
            eig[m].pair = 0;
            m++;
            i++;
        }
    }
    *neig = m;
    lis_gcrodr_sort_eigs(eig,m);
    return LIS_SUCCESS;
}

#endif /* !_COMPLEX */

/* ---------------------------------------------------------------------- */
/* persistent recycle state                                               */
/* ---------------------------------------------------------------------- */

LIS_INT lis_gcrodr_clear_recycle(LIS_SOLVER solver)
{
    LIS_INT i;
    if( solver==NULL ) return LIS_SUCCESS;

    if( solver->recycle_u )
    {
        for(i=0;i<solver->recycle_alloc;i++) lis_vector_destroy(solver->recycle_u[i]);
        lis_free(solver->recycle_u);
    }
    if( solver->recycle_zu )
    {
        for(i=0;i<solver->recycle_alloc;i++) lis_vector_destroy(solver->recycle_zu[i]);
        lis_free(solver->recycle_zu);
    }
    if( solver->recycle_c )
    {
        for(i=0;i<solver->recycle_alloc;i++) lis_vector_destroy(solver->recycle_c[i]);
        lis_free(solver->recycle_c);
    }

    solver->recycle_u = NULL;
    solver->recycle_zu = NULL;
    solver->recycle_c = NULL;
    solver->recycle_dim = 0;
    solver->recycle_alloc = 0;
    solver->recycle_valid = LIS_FALSE;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_gcrodr_recycle"
LIS_INT lis_solver_set_gcrodr_recycle(LIS_SOLVER solver, LIS_INT dim,
                                       LIS_REAL keep_quality,
                                       LIS_REAL reject_quality)
{
    if( solver==NULL )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"GCRO-DR solver is NULL\n");
        return LIS_ERR_ILL_ARG;
    }
    if( dim<0 )
    {
        LIS_SETERR1(LIS_ERR_ILL_ARG,
                    "GCRO-DR recycle dimension %D is negative\n",dim);
        return LIS_ERR_ILL_ARG;
    }
    if( reject_quality!=reject_quality || keep_quality!=keep_quality ||
        reject_quality<0.0 || reject_quality>1.0 ||
        keep_quality<0.0 || keep_quality>1.0 ||
        reject_quality>keep_quality )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "GCRO-DR quality thresholds must satisfy 0 <= reject <= keep <= 1\n");
        return LIS_ERR_ILL_ARG;
    }

    if( solver->recycle_target_dim!=dim )
    {
        lis_gcrodr_clear_recycle(solver);
    }

    solver->recycle_target_dim = dim;
    solver->recycle_keep_threshold = keep_quality;
    solver->recycle_reject_threshold = reject_quality;
    return LIS_SUCCESS;
}

LIS_INT lis_solver_get_gcrodr_recycle(LIS_SOLVER solver, LIS_INT *dim,
                                       LIS_INT *action,
                                       LIS_REAL *quality,
                                       LIS_REAL *quality_min)
{
    if( solver==NULL ) return LIS_ERR_ILL_ARG;
    if( dim ) *dim = solver->recycle_dim;
    if( action ) *action = solver->recycle_last_action;
    if( quality ) *quality = solver->recycle_quality;
    if( quality_min ) *quality_min = solver->recycle_quality_min;
    return LIS_SUCCESS;
}

static LIS_INT lis_gcrodr_recycle_compatible(LIS_SOLVER solver)
{
    LIS_INT i;
    if( solver->recycle_alloc<=0 || solver->recycle_u==NULL ||
        solver->recycle_zu==NULL || solver->recycle_c==NULL ) return LIS_FALSE;

    for(i=0;i<solver->recycle_alloc;i++)
    {
        if( !lis_solver_work_vector_compatible(solver->recycle_u[i],solver->A,
                                               LIS_PRECISION_DEFAULT) ||
            !lis_solver_work_vector_compatible(solver->recycle_zu[i],solver->A,
                                               LIS_PRECISION_DEFAULT) ||
            !lis_solver_work_vector_compatible(solver->recycle_c[i],solver->A,
                                               LIS_PRECISION_DEFAULT) )
            return LIS_FALSE;
    }
    return LIS_TRUE;
}

static LIS_INT lis_gcrodr_ensure_recycle(LIS_SOLVER solver, LIS_INT nvec)
{
    LIS_INT i,j,err;
    LIS_VECTOR *u,*zu,*c;

    if( nvec<=0 )
    {
        lis_gcrodr_clear_recycle(solver);
        return LIS_SUCCESS;
    }

    if( solver->recycle_alloc==nvec && lis_gcrodr_recycle_compatible(solver) )
        return LIS_SUCCESS;

    lis_gcrodr_clear_recycle(solver);

    u  = (LIS_VECTOR *)lis_malloc((size_t)nvec*sizeof(LIS_VECTOR),
                                  "lis_gcrodr_ensure_recycle::u");
    zu = (LIS_VECTOR *)lis_malloc((size_t)nvec*sizeof(LIS_VECTOR),
                                  "lis_gcrodr_ensure_recycle::zu");
    c  = (LIS_VECTOR *)lis_malloc((size_t)nvec*sizeof(LIS_VECTOR),
                                  "lis_gcrodr_ensure_recycle::c");
    if( u==NULL || zu==NULL || c==NULL )
    {
        lis_free(u); lis_free(zu); lis_free(c);
        return LIS_ERR_OUT_OF_MEMORY;
    }
    for(i=0;i<nvec;i++) { u[i]=NULL; zu[i]=NULL; c[i]=NULL; }

    for(i=0;i<nvec;i++)
    {
        err = lis_vector_duplicate(solver->A,&u[i]);
        if( err ) break;
        err = lis_vector_duplicate(solver->A,&zu[i]);
        if( err ) break;
        err = lis_vector_duplicate(solver->A,&c[i]);
        if( err ) break;
        lis_vector_set_all(0.0,u[i]);
        lis_vector_set_all(0.0,zu[i]);
        lis_vector_set_all(0.0,c[i]);
    }
    if( i<nvec )
    {
        for(j=0;j<nvec;j++)
        {
            lis_vector_destroy(u[j]);
            lis_vector_destroy(zu[j]);
            lis_vector_destroy(c[j]);
        }
        lis_free(u); lis_free(zu); lis_free(c);
        return err;
    }

    solver->recycle_u = u;
    solver->recycle_zu = zu;
    solver->recycle_c = c;
    solver->recycle_alloc = nvec;
    solver->recycle_dim = 0;
    solver->recycle_valid = LIS_FALSE;
    return LIS_SUCCESS;
}

/* Refresh ZU=M^-1 U and C=A ZU for the current A/preconditioner.  C is
 * orthonormalized, and exactly the same transformations are applied to U/ ZU.
 */
static LIS_INT lis_gcrodr_refresh_recycle(LIS_SOLVER solver)
{
    LIS_INT i,j,pass,keep,err;
    LIS_SCALAR alpha;
    LIS_REAL nrm;
    LIS_VECTOR ui,zui,ci;

    if( !solver->recycle_valid || solver->recycle_dim<=0 ) return LIS_SUCCESS;
    if( !lis_gcrodr_recycle_compatible(solver) )
    {
        solver->recycle_rejects++;
        solver->recycle_dim = 0;
        solver->recycle_valid = LIS_FALSE;
        return LIS_SUCCESS;
    }

    keep = 0;
    for(i=0;i<solver->recycle_dim;i++)
    {
        ui = solver->recycle_u[i];
        zui = solver->recycle_zu[i];
        ci = solver->recycle_c[i];

        err = lis_psolve(solver,ui,zui);
        if( err ) return err;
        lis_matvec(solver->A,zui,ci);

        for(pass=0;pass<2;pass++)
        {
            for(j=0;j<keep;j++)
            {
                lis_vector_dot(solver->recycle_c[j],ci,&alpha);
                lis_vector_axpy(-alpha,solver->recycle_c[j],ci);
                lis_vector_axpy(-alpha,solver->recycle_u[j],ui);
                lis_vector_axpy(-alpha,solver->recycle_zu[j],zui);
            }
        }

        lis_vector_nrm2(ci,&nrm);
        if( nrm<=GCRODR_REORTH_TINY || nrm!=nrm )
        {
            solver->recycle_rejects++;
            continue;
        }

        lis_vector_scale(1.0/nrm,ci);
        lis_vector_scale(1.0/nrm,ui);
        lis_vector_scale(1.0/nrm,zui);

        if( keep!=i )
        {
            lis_vector_copy(ui,solver->recycle_u[keep]);
            lis_vector_copy(zui,solver->recycle_zu[keep]);
            lis_vector_copy(ci,solver->recycle_c[keep]);
        }
        keep++;
    }

    solver->recycle_dim = keep;
    solver->recycle_valid = (keep>0) ? LIS_TRUE : LIS_FALSE;
    if( keep>0 ) solver->recycle_uses++;
    return LIS_SUCCESS;
}


/* Cross-solve quality test.
 *
 * The old C basis is orthonormal and represents the previous
 * right-preconditioned operator image.  For each old recycle direction u_i
 * compute
 *
 *     c_i^new = A_new M_new^{-1} u_i
 *
 * and measure how much of c_i^new lies in span(C_old):
 *
 *     q_i = || C_old^T c_i^new ||_2 / || c_i^new ||_2.
 *
 * This score is invariant to rotations inside the old recycle subspace and
 * therefore is more robust than comparing only c_i^old with c_i^new.
 *
 * KEEP    : every direction remains strongly aligned.
 * REFRESH : the recycle space remains useful, but the operator changed
 *           enough to require a full remap and re-QR.
 * REJECT  : too little of the old image survives; rebuild from the new
 *           Krylov sequence.
 *
 * scratch_z and scratch_c must provide at least recycle_dim compatible
 * vectors.  GCRO-DR passes unused Arnoldi work vectors before the new
 * residual is formed, so no additional large persistent allocation is
 * required.
 */
static LIS_INT lis_gcrodr_assess_refresh_recycle(
    LIS_SOLVER solver, LIS_VECTOR *scratch_z, LIS_VECTOR *scratch_c)
{
    LIS_INT i,j,good,need,err,action,oldk;
    LIS_SCALAR alpha;
    LIS_REAL nrm,proj2,qi,qsum,qmin;

    if( !solver->recycle_valid || solver->recycle_dim<=0 ) return LIS_SUCCESS;

    if( !lis_gcrodr_recycle_compatible(solver) )
    {
        solver->recycle_policy_rejects++;
        solver->recycle_rejects++;
        solver->recycle_last_action = LIS_GCRODR_RECYCLE_REJECT;
        solver->recycle_quality = 0.0;
        solver->recycle_quality_min = 0.0;
        solver->recycle_dim = 0;
        solver->recycle_valid = LIS_FALSE;
        return LIS_SUCCESS;
    }

    oldk = solver->recycle_dim;

    qsum = 0.0;
    qmin = 1.0;
    good = 0;

    /* First pass: build the new operator images without overwriting old C. */
    for(i=0;i<oldk;i++)
    {
        err = lis_psolve(solver,solver->recycle_u[i],scratch_z[i]);
        if( err )
        {
            return err;
        }
        lis_matvec(solver->A,scratch_z[i],scratch_c[i]);
        lis_vector_nrm2(scratch_c[i],&nrm);

        if( nrm<=GCRODR_REORTH_TINY || nrm!=nrm )
        {
            qi = 0.0;
        }
        else
        {
            proj2 = 0.0;
            for(j=0;j<oldk;j++)
            {
                lis_vector_dot(solver->recycle_c[j],scratch_c[i],&alpha);
                proj2 += (LIS_REAL)(fabs(alpha)*fabs(alpha));
            }
            if( proj2<0.0 ) proj2 = 0.0;
            qi = sqrt(proj2)/nrm;
            if( qi<0.0 ) qi = 0.0;
            if( qi>1.0 ) qi = 1.0;
        }

        qsum += qi;
        if( qi<qmin ) qmin = qi;
        if( qi>=solver->recycle_reject_threshold ) good++;
    }

    solver->recycle_quality = qsum/(LIS_REAL)oldk;
    solver->recycle_quality_min = qmin;

    need = _max((LIS_INT)1,(oldk+1)/2);

    if( qmin>=solver->recycle_keep_threshold )
    {
        action = LIS_GCRODR_RECYCLE_KEEP;
    }
    else if( solver->recycle_quality<solver->recycle_reject_threshold ||
             good<need )
    {
        action = LIS_GCRODR_RECYCLE_REJECT;
    }
    else
    {
        action = LIS_GCRODR_RECYCLE_REFRESH;
    }

    if( action==LIS_GCRODR_RECYCLE_REJECT )
    {
        solver->recycle_policy_rejects++;
        solver->recycle_rejects++;
        solver->recycle_last_action = action;
        solver->recycle_dim = 0;
        solver->recycle_valid = LIS_FALSE;
        return LIS_SUCCESS;
    }

    /*
     * Keep the policy assessment read-only with respect to U.  The refresh
     * path rebuilds
     *
     *     ZU = M^-1 U,   C = A ZU
     *
     * and applies a consistent QR transformation to U/ZU/C.  KEEP and
     * REFRESH both remap the complete accepted recycle space; REJECT is the
     * action that discards it.
     */
    err = lis_gcrodr_refresh_recycle(solver);
    if( err )
    {
        lis_gcrodr_clear_recycle(solver);
        return err;
    }

    if( !solver->recycle_valid || solver->recycle_dim<=0 )
    {
        solver->recycle_policy_rejects++;
        solver->recycle_last_action = LIS_GCRODR_RECYCLE_REJECT;
        solver->recycle_quality = 0.0;
        solver->recycle_quality_min = 0.0;
        return LIS_SUCCESS;
    }

    solver->recycle_cross_uses++;
    solver->recycle_last_action = action;

    if( action==LIS_GCRODR_RECYCLE_KEEP )
        solver->recycle_policy_keeps++;
    else
        solver->recycle_policy_refreshes++;

    return LIS_SUCCESS;
}

/* ---------------------------------------------------------------------- */
/* harmonic Ritz update                                                   */
/* ---------------------------------------------------------------------- */
#ifndef _COMPLEX

static LIS_INT lis_gcrodr_build_vector(LIS_SOLVER solver,
                                       LIS_VECTOR *v, LIS_INT jdim,
                                       LIS_INT oldk,
                                       const LIS_REAL *coef,
                                       LIS_VECTOR out)
{
    LIS_INT i;
    LIS_REAL nrm;
    lis_vector_set_all(0.0,out);
    for(i=0;i<oldk;i++)
        lis_vector_axpy((LIS_SCALAR)coef[i],solver->recycle_u[i],out);
    for(i=0;i<jdim;i++)
        lis_vector_axpy((LIS_SCALAR)coef[oldk+i],v[i],out);
    lis_vector_nrm2(out,&nrm);
    if( nrm<=GCRODR_TINY || nrm!=nrm ) return LIS_BREAKDOWN;
    lis_vector_scale(1.0/nrm,out);
    return LIS_SUCCESS;
}

static LIS_INT lis_gcrodr_update_recycle(LIS_SOLVER solver,
                                         LIS_VECTOR *v, LIS_INT jdim,
                                         const LIS_SCALAR *hraw, LIS_INT hld,
                                         const LIS_SCALAR *bcoef, LIS_INT bld)
{
    LIS_INT oldk,kmax,knew,p,q;
    LIS_INT i,j,l,r,c,err,qriter,neig,nsel,eidx;
    LIS_REAL qrerr;
    LIS_SCALAR *G,*F,*AP,*BP,*T,*TQR,*QW,*RW,*AINV;
    LIS_GCRODR_EIG *eig;
    LIS_REAL *vr,*vi;
    LIS_VECTOR *newu;

    kmax = _min(solver->recycle_target_dim, solver->options[LIS_OPTIONS_RESTART]-1);
    oldk = solver->recycle_valid ? solver->recycle_dim : 0;
    if( kmax<=0 || jdim<=0 ) return LIS_SUCCESS;

    p = oldk+jdim;
    q = oldk+jdim+1;
    if( p<=0 ) return LIS_SUCCESS;
    knew = _min(kmax,p);

    G    = (LIS_SCALAR *)lis_malloc((size_t)q*p*sizeof(LIS_SCALAR),"gcrodr::G");
    F    = (LIS_SCALAR *)lis_malloc((size_t)q*p*sizeof(LIS_SCALAR),"gcrodr::F");
    AP   = (LIS_SCALAR *)lis_malloc((size_t)p*p*sizeof(LIS_SCALAR),"gcrodr::AP");
    BP   = (LIS_SCALAR *)lis_malloc((size_t)p*p*sizeof(LIS_SCALAR),"gcrodr::BP");
    T    = (LIS_SCALAR *)lis_malloc((size_t)p*p*sizeof(LIS_SCALAR),"gcrodr::T");
    TQR  = (LIS_SCALAR *)lis_malloc((size_t)p*p*sizeof(LIS_SCALAR),"gcrodr::TQR");
    QW   = (LIS_SCALAR *)lis_malloc((size_t)p*p*sizeof(LIS_SCALAR),"gcrodr::QW");
    RW   = (LIS_SCALAR *)lis_malloc((size_t)p*p*sizeof(LIS_SCALAR),"gcrodr::RW");
    AINV = (LIS_SCALAR *)lis_malloc((size_t)p*p*sizeof(LIS_SCALAR),"gcrodr::AINV");
    eig  = (LIS_GCRODR_EIG *)lis_malloc((size_t)p*sizeof(*eig),"gcrodr::eig");
    vr   = (LIS_REAL *)lis_malloc((size_t)p*sizeof(LIS_REAL),"gcrodr::vr");
    vi   = (LIS_REAL *)lis_malloc((size_t)p*sizeof(LIS_REAL),"gcrodr::vi");
    newu = (LIS_VECTOR *)lis_malloc((size_t)knew*sizeof(LIS_VECTOR),"gcrodr::newu");
    if( newu!=NULL )
        for(i=0;i<knew;i++) newu[i]=NULL;

    if( G==NULL || F==NULL || AP==NULL || BP==NULL || T==NULL || TQR==NULL ||
        QW==NULL || RW==NULL || AINV==NULL || eig==NULL || vr==NULL || vi==NULL ||
        newu==NULL )
    {
        err = LIS_ERR_OUT_OF_MEMORY;
        goto cleanup;
    }

    memset(G,0,(size_t)q*p*sizeof(LIS_SCALAR));
    memset(F,0,(size_t)q*p*sizeof(LIS_SCALAR));

    /* Atilde [U,V] = [C,Vnext] G, Atilde = A M^{-1}. */
    for(i=0;i<oldk;i++) G[i+i*q] = 1.0;
    for(j=0;j<jdim;j++)
    {
        for(i=0;i<oldk;i++) G[i+(oldk+j)*q] = bcoef[i+j*bld];
        for(i=0;i<=jdim;i++) G[(oldk+i)+(oldk+j)*q] = hraw[i+j*hld];
    }

    /* F = [C,Vnext]^T [U,V].  C^T V is zero by projected Arnoldi;
     * Vnext^T V has [I;0] structure.  Only C^T U and Vnext^T U require dots.
     */
    for(c=0;c<oldk;c++)
    {
        for(r=0;r<oldk;r++)
            lis_vector_dot(solver->recycle_c[r],solver->recycle_u[c],
                           &F[r+c*q]);
        for(r=0;r<=jdim;r++)
            lis_vector_dot(v[r],solver->recycle_u[c],
                           &F[(oldk+r)+c*q]);
    }
    for(c=0;c<jdim;c++) F[(oldk+c)+(oldk+c)*q] = 1.0;

    /* AP=G^T G, BP=G^T F. */
    for(c=0;c<p;c++)
    {
        for(r=0;r<p;r++)
        {
            LIS_SCALAR ag=0.0,bg=0.0;
            for(l=0;l<q;l++)
            {
                ag += G[l+r*q]*G[l+c*q];
                bg += G[l+r*q]*F[l+c*q];
            }
            AP[r+c*p] = ag;
            BP[r+c*p] = bg;
        }
    }

    memcpy(AINV,AP,(size_t)p*p*sizeof(LIS_SCALAR));
    err = lis_array_ge(p,AINV);
    if( err ) goto reject;

    /* T = (G^T G)^-1 (G^T F), eig(T)=1/theta_harmonic. */
    for(c=0;c<p;c++)
    {
        for(r=0;r<p;r++)
        {
            LIS_SCALAR s=0.0;
            for(l=0;l<p;l++) s += AINV[r+l*p]*BP[l+c*p];
            T[r+c*p] = s;
        }
    }

    memcpy(TQR,T,(size_t)p*p*sizeof(LIS_SCALAR));
    err = lis_array_qr(p,TQR,QW,RW,&qriter,&qrerr);
    if( err ) goto reject;
    (void)qriter;
    (void)qrerr;

    err = lis_gcrodr_extract_eigs(p,TQR,eig,&neig);
    if( err ) goto reject;

    for(i=0;i<knew;i++)
    {
        err = lis_vector_duplicate(solver->A,&newu[i]);
        if( err ) goto cleanup;
    }

    nsel = 0;
    for(eidx=0;eidx<neig && nsel<knew;eidx++)
    {
        if( eig[eidx].pair && nsel+2>knew ) continue;
        err = lis_gcrodr_inverse_eigenvector(p,T,eig[eidx].lr,eig[eidx].li,vr,vi);
        if( err==LIS_ERR_OUT_OF_MEMORY ) goto cleanup;
        if( err ) continue;

        if( eig[eidx].pair )
        {
            err = lis_gcrodr_build_vector(solver,v,jdim,oldk,vr,newu[nsel]);
            if( err ) continue;
            err = lis_gcrodr_build_vector(solver,v,jdim,oldk,vi,newu[nsel+1]);
            if( err ) continue;
            nsel += 2;
        }
        else
        {
            err = lis_gcrodr_build_vector(solver,v,jdim,oldk,vr,newu[nsel]);
            if( err ) continue;
            nsel++;
        }
    }

    if( nsel<=0 ) goto reject;

    for(i=0;i<nsel;i++) lis_vector_copy(newu[i],solver->recycle_u[i]);
    solver->recycle_dim = nsel;
    solver->recycle_valid = LIS_TRUE;
    if( oldk>0 )
    {
        solver->recycle_updates++;
    }
    else
    {
        if( solver->recycle_builds==0 &&
            solver->recycle_last_action==LIS_GCRODR_RECYCLE_NONE )
            solver->recycle_last_action = LIS_GCRODR_RECYCLE_BUILD;
        solver->recycle_builds++;
    }

    err = lis_gcrodr_refresh_recycle(solver);
    if( err )
    {
        lis_gcrodr_clear_recycle(solver);
        goto cleanup;
    }
    goto cleanup;

reject:
    solver->recycle_rejects++;
    err = LIS_SUCCESS; /* recycling failure must not kill the linear solve */

cleanup:
    if( newu )
    {
        for(i=0;i<knew;i++) lis_vector_destroy(newu[i]);
    }
    lis_free(newu);
    lis_free(G); lis_free(F); lis_free(AP); lis_free(BP); lis_free(T);
    lis_free(TQR); lis_free(QW); lis_free(RW); lis_free(AINV);
    lis_free(eig); lis_free(vr); lis_free(vi);
    return err;
}

#endif /* !_COMPLEX */

/* ---------------------------------------------------------------------- */
/* public solver hooks                                                    */
/* ---------------------------------------------------------------------- */

#undef __FUNC__
#define __FUNC__ "lis_gcrodr_check_params"
LIS_INT lis_gcrodr_check_params(LIS_SOLVER solver)
{
    LIS_INT restart;
    LIS_DEBUG_FUNC_IN;

    restart = solver->options[LIS_OPTIONS_RESTART];
    if( restart<2 )
    {
        LIS_SETERR1(LIS_ERR_ILL_ARG,
                    "Parameter LIS_OPTIONS_RESTART(=%D) must be at least 2 for GCRO-DR\n",
                    restart);
        return LIS_ERR_ILL_ARG;
    }
    if( solver->recycle_target_dim<0 )
    {
        LIS_SETERR1(LIS_ERR_ILL_ARG,
                    "GCRO-DR recycle dimension %D is negative\n",
                    solver->recycle_target_dim);
        return LIS_ERR_ILL_ARG;
    }
    if( solver->recycle_reject_threshold!=solver->recycle_reject_threshold ||
        solver->recycle_keep_threshold!=solver->recycle_keep_threshold ||
        solver->recycle_reject_threshold<0.0 ||
        solver->recycle_reject_threshold>1.0 ||
        solver->recycle_keep_threshold<0.0 ||
        solver->recycle_keep_threshold>1.0 ||
        solver->recycle_reject_threshold>solver->recycle_keep_threshold )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "GCRO-DR quality thresholds must satisfy 0 <= reject <= keep <= 1\n");
        return LIS_ERR_ILL_ARG;
    }

#ifdef _COMPLEX
    LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
               "GCRO-DR currently supports real scalar builds only\n");
    return LIS_ERR_NOT_IMPLEMENTED;
#endif
    if( solver->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE )
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                   "GCRO-DR currently supports default/double precision only\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_gcrodr_malloc_work"
LIS_INT lis_gcrodr_malloc_work(LIS_SOLVER solver)
{
    LIS_VECTOR *work;
    LIS_INT i,j,restart,worklen,err;

    LIS_DEBUG_FUNC_IN;
    restart = solver->options[LIS_OPTIONS_RESTART];
    worklen = GCRODR_NWORK+(2*restart+1);
    work = (LIS_VECTOR *)lis_malloc((size_t)worklen*sizeof(LIS_VECTOR),
                                    "lis_gcrodr_malloc_work::work");
    if( work==NULL ) return LIS_ERR_OUT_OF_MEMORY;
    for(i=0;i<worklen;i++) work[i]=NULL;

    for(i=1;i<worklen;i++)
    {
        err = lis_vector_duplicate(solver->A,&work[i]);
        if( err ) break;
    }
    if( i<worklen )
    {
        for(j=1;j<i;j++) lis_vector_destroy(work[j]);
        lis_free(work);
        return err;
    }

    err = lis_vector_create(solver->A->comm,&work[0]);
    if( err )
    {
        for(j=1;j<worklen;j++) lis_vector_destroy(work[j]);
        lis_free(work);
        return err;
    }
    err = lis_vector_set_size(work[0],restart+1,0);
    if( err )
    {
        lis_vector_destroy(work[0]);
        for(j=1;j<worklen;j++) lis_vector_destroy(work[j]);
        lis_free(work);
        return err;
    }

    solver->worklen = worklen;
    solver->work = work;
    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_gcrodr"
LIS_INT lis_gcrodr(LIS_SOLVER solver)
{
#ifdef _COMPLEX
    (void)solver;
    return LIS_ERR_NOT_IMPLEMENTED;
#else
    LIS_Comm comm;
    LIS_MATRIX A;
    LIS_VECTOR b,x,s,*z,*v,wtmp;
    LIS_SCALAR *h,*hraw,*bcoef;
    LIS_SCALAR aa,bb,rr,a2,b2,t,alpha,gamma;
    LIS_REAL bnrm,nrm2,tol,min_nrm2,rnorm,wnorm;
    LIS_INT err,iter,maxiter,output,maxiter_noimp,noimp_count;
    LIS_INT i,j,k,m,kmax,kcur,cycle_steps,jdim,ii,i1,jj,hld,cs,sn,pass;
    LIS_INT lucky;
    double time,ptime;

    LIS_DEBUG_FUNC_IN;
    comm = LIS_COMM_WORLD;
    A = solver->A;
    b = solver->b;
    x = solver->x;
    maxiter = solver->options[LIS_OPTIONS_MAXITER];
    output = solver->options[LIS_OPTIONS_OUTPUT];
    m = solver->options[LIS_OPTIONS_RESTART];
    if( solver->recycle_target_dim<0 )
        solver->recycle_target_dim = GCRODR_DEFAULT_RECYCLE_DIM;
    kmax = _min(solver->recycle_target_dim, m-1);
    maxiter_noimp = solver->options[LIS_OPTIONS_MAXITER_NO_IMP];
    min_nrm2 = LIS_SCALAR_MAX;
    noimp_count = 0;
    ptime = 0.0;
    iter = 0;
    nrm2 = LIS_SCALAR_MAX;

    s = solver->work[0];
    wtmp = solver->work[1];
    z = &solver->work[2];
    v = &solver->work[m+2];
    hld = m+1;

    h = (LIS_SCALAR *)lis_malloc((size_t)(hld+1)*(hld+2)*sizeof(LIS_SCALAR),
                                 "lis_gcrodr::h");
    hraw = (LIS_SCALAR *)lis_malloc((size_t)(m+1)*m*sizeof(LIS_SCALAR),
                                    "lis_gcrodr::hraw");
    bcoef = (LIS_SCALAR *)lis_malloc((size_t)_max(1,kmax)*m*sizeof(LIS_SCALAR),
                                     "lis_gcrodr::bcoef");
    if( h==NULL || hraw==NULL || bcoef==NULL )
    {
        lis_free(h); lis_free(hraw); lis_free(bcoef);
        return LIS_ERR_OUT_OF_MEMORY;
    }
    cs = (m+1)*hld;
    sn = (m+2)*hld;

    err = lis_gcrodr_ensure_recycle(solver,kmax);
    if( err ) goto fail;

    /*
     * Assess cross-solve reuse before v[0] is occupied by the new residual.
     * z[] and v[] are unused scratch storage at this point.
     */
    if( solver->recycle_valid && solver->recycle_dim>0 )
    {
        err = lis_gcrodr_assess_refresh_recycle(solver,z,v);
        if( err ) goto fail;
    }

    /* Initial true residual in v[0]; bnrm is LIS' fixed residual scaling. */
    if( lis_solver_get_initial_residual(solver,NULL,NULL,v[0],&bnrm) )
    {
        err = LIS_SUCCESS;
        goto done;
    }
    tol = solver->tol;

    if( solver->recycle_valid && solver->recycle_dim>0 )
    {
        /* x <- x + ZU C^T r, r <- (I-CC^T)r */
        for(i=0;i<solver->recycle_dim;i++)
        {
            lis_vector_dot(solver->recycle_c[i],v[0],&alpha);
            lis_vector_axpy(alpha,solver->recycle_zu[i],x);
            lis_vector_axpy(-alpha,solver->recycle_c[i],v[0]);
        }
        lis_vector_nrm2(v[0],&rnorm);
        nrm2 = rnorm*bnrm;
        if( nrm2<=tol )
        {
            solver->retcode = LIS_SUCCESS;
            solver->iter = 0;
            solver->resid = nrm2;
            solver->ptime = ptime;
            err = LIS_SUCCESS;
            goto done;
        }
    }

    while( iter<maxiter )
    {
        kcur = solver->recycle_valid ? solver->recycle_dim : 0;
        memset(h,0,(size_t)(hld+1)*(hld+2)*sizeof(LIS_SCALAR));
        memset(hraw,0,(size_t)(m+1)*m*sizeof(LIS_SCALAR));
        memset(bcoef,0,(size_t)_max(1,kmax)*m*sizeof(LIS_SCALAR));

        lis_vector_nrm2(v[0],&rnorm);
        if( rnorm<=GCRODR_TINY || rnorm!=rnorm )
        {
            nrm2 = rnorm*bnrm;
            err = (nrm2<=tol) ? LIS_SUCCESS : LIS_BREAKDOWN;
            break;
        }
        lis_vector_scale(1.0/rnorm,v[0]);
        lis_vector_set_all(0.0,s);
        s->value[0] = rnorm;

        jdim = 0;
        lucky = LIS_FALSE;
        cycle_steps = m-kcur;
        if( cycle_steps<1 ) cycle_steps = 1;
        for(j=0;j<cycle_steps && iter<maxiter;j++)
        {
            iter++;

            time = lis_wtime();
            err = lis_psolve(solver,v[j],z[j]);
            ptime += lis_wtime()-time;
            if( err ) goto fail;

            lis_matvec(A,z[j],v[j+1]);

            /* Project A M^-1 v_j away from C, with a second MGS pass. */
            for(pass=0;pass<2;pass++)
            {
                for(i=0;i<kcur;i++)
                {
                    lis_vector_dot(solver->recycle_c[i],v[j+1],&t);
                    bcoef[i+j*_max(1,kmax)] += t;
                    lis_vector_axpy(-t,solver->recycle_c[i],v[j+1]);
                }
            }

            /* Arnoldi MGS + reorthogonalization. */
            for(pass=0;pass<2;pass++)
            {
                for(i=0;i<=j;i++)
                {
                    lis_vector_dot(v[i],v[j+1],&t);
                    hraw[i+j*hld] += t;
                    lis_vector_axpy(-t,v[i],v[j+1]);
                }
            }
            lis_vector_nrm2(v[j+1],&wnorm);
            hraw[(j+1)+j*hld] = wnorm;
            if( wnorm>GCRODR_TINY && wnorm==wnorm )
                lis_vector_scale(1.0/wnorm,v[j+1]);
            else
                lucky = LIS_TRUE;

            /* Copy this Hessenberg column to the QR work matrix. */
            for(i=0;i<=j+1;i++) h[i+j*hld] = hraw[i+j*hld];

            /* Apply previous Givens rotations. */
            for(k=1;k<=j;k++)
            {
                jj = k-1;
                t = h[jj+j*hld];
                aa = h[jj+cs]*t + h[jj+sn]*h[k+j*hld];
                bb = -h[jj+sn]*t + h[jj+cs]*h[k+j*hld];
                h[jj+j*hld] = aa;
                h[k+j*hld] = bb;
            }

            ii = j;
            i1 = j+1;
            aa = h[ii+j*hld];
            bb = h[i1+j*hld];
            a2 = aa*aa;
            b2 = bb*bb;
            rr = sqrt(a2+b2);
            if( rr==0.0 || rr!=rr )
            {
                if( lucky )
                {
                    jdim = j+1;
                    break;
                }
                err = LIS_BREAKDOWN;
                goto fail;
            }
            h[ii+cs] = aa/rr;
            h[ii+sn] = bb/rr;
            s->value[i1] = -h[ii+sn]*s->value[ii];
            s->value[ii] =  h[ii+cs]*s->value[ii];
            h[ii+j*hld] = h[ii+cs]*h[ii+j*hld] + h[ii+sn]*h[i1+j*hld];

            nrm2 = fabs(s->value[i1])*bnrm;
            if( output )
            {
                if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
                if( output & LIS_PRINT_OUT ) lis_print_rhistory(comm,iter,nrm2);
            }

            jdim = j+1;
            if( nrm2<=tol || lucky ) break;

            if( maxiter_noimp )
            {
                if( min_nrm2>nrm2 )
                {
                    min_nrm2=nrm2;
                    noimp_count=0;
                }
                else if( ++noimp_count>maxiter_noimp )
                {
                    err = LIS_MAXITER;
                    goto fail;
                }
            }
        }

        if( jdim<=0 )
        {
            err = LIS_BREAKDOWN;
            goto fail;
        }

        /* Back substitution in rotated H. */
        ii = jdim-1;
        if( h[ii+ii*hld]==0.0 || h[ii+ii*hld]!=h[ii+ii*hld] )
        {
            err = LIS_BREAKDOWN;
            goto fail;
        }
        s->value[ii] /= h[ii+ii*hld];
        for(k=1;k<=ii;k++)
        {
            jj = ii-k;
            t = s->value[jj];
            for(j=jj+1;j<=ii;j++) t -= h[jj+j*hld]*s->value[j];
            if( h[jj+jj*hld]==0.0 || h[jj+jj*hld]!=h[jj+jj*hld] )
            {
                err = LIS_BREAKDOWN;
                goto fail;
            }
            s->value[jj] = t/h[jj+jj*hld];
        }

        /* GCRO correction in physical x-space:
         *   delta x = M^-1(V y - U B y) = Z y - ZU(B y).
         */
        for(j=0;j<jdim;j++) lis_vector_axpy(s->value[j],z[j],x);
        for(i=0;i<kcur;i++)
        {
            gamma = 0.0;
            for(j=0;j<jdim;j++)
                gamma += bcoef[i+j*_max(1,kmax)]*s->value[j];
            lis_vector_axpy(-gamma,solver->recycle_zu[i],x);
        }

        /* Build/update recycle space before v[0] is overwritten by true r. */
        if( kmax>0 && !lucky )
        {
            err = lis_gcrodr_update_recycle(solver,v,jdim,hraw,hld,
                                            bcoef,_max(1,kmax));
            if( err ) goto fail;
        }

        /* True residual replacement at every restart boundary. */
        lis_matvec(A,x,wtmp);
        lis_vector_xpay(b,-1.0,wtmp);
        lis_vector_copy(wtmp,v[0]);

        /* New C may have changed after harmonic-Ritz update.  Reproject r and
         * apply the matching solution correction to preserve consistency.
         */
        if( solver->recycle_valid )
        {
            for(i=0;i<solver->recycle_dim;i++)
            {
                lis_vector_dot(solver->recycle_c[i],v[0],&alpha);
                lis_vector_axpy(alpha,solver->recycle_zu[i],x);
                lis_vector_axpy(-alpha,solver->recycle_c[i],v[0]);
            }
        }
        lis_vector_nrm2(v[0],&rnorm);
        nrm2 = rnorm*bnrm;

        if( nrm2!=nrm2 || fabs(nrm2)>LIS_SCALAR_MAX )
        {
            err = LIS_BREAKDOWN;
            goto fail;
        }
        if( output & LIS_PRINT_MEM ) solver->rhistory[iter] = nrm2;
        if( nrm2<=tol )
        {
            err = LIS_SUCCESS;
            break;
        }
    }

    if( iter>=maxiter && nrm2>tol ) err = LIS_MAXITER;
    solver->retcode = err;
    solver->iter = iter;
    solver->resid = nrm2;
    solver->ptime = ptime;
    goto done;

fail:
    solver->retcode = err;
    solver->iter = iter;
    solver->resid = nrm2;
    solver->ptime = ptime;

done:
    lis_free(h);
    lis_free(hraw);
    lis_free(bcoef);
    LIS_DEBUG_FUNC_OUT;
    return err;
#endif
}

#undef GCRODR_NWORK
#undef GCRODR_TINY
#undef GCRODR_REORTH_TINY
