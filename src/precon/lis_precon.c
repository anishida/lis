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
#include <math.h>
#include "lislib.h"

/************************************************
 * lis_precon_init
 * lis_precon_create
 * lis_precon_destroy
 * lis_psolve
 * lis_psolveh
 ************************************************/

LIS_PRECON_CREATE_XXX lis_precon_create_xxx[] = {
	lis_precon_create_none, lis_precon_create_jacobi, lis_precon_create_iluk,
	lis_precon_create_ssor, lis_precon_create_hybrid, lis_precon_create_is,
	lis_precon_create_sainv,lis_precon_create_saamg,  lis_precon_create_iluc,
	lis_precon_create_ilut, lis_precon_create_bjacobi
};

/*NEH support for extended "solve_kernel" workflow*/
LIS_PRECON_PSD_CREATE_XXX lis_precon_psd_create_xxx[] = {
	lis_precon_create_none, NULL,                        lis_precon_psd_create_iluk,
	NULL,                   NULL,                        NULL,
	NULL,                   lis_precon_psd_create_saamg, NULL,
	NULL,                   NULL
};

/*NEH support for extended "solve_kernel" workflow*/
LIS_PRECON_PSD_UPDATE_XXX lis_precon_psd_update_xxx[] = {
	NULL,                   NULL,                        lis_precon_psd_update_iluk,
	NULL,                   NULL,                        NULL,
	NULL,                   lis_precon_psd_update_saamg, NULL,
	NULL,                   NULL
};

LIS_PSOLVE_XXX lis_psolve_xxx[] = {
	lis_psolve_none,     lis_psolve_jacobi, lis_psolve_iluk_csr,
	lis_psolve_ssor,     lis_psolve_hybrid, lis_psolve_none,
	lis_psolve_sainv,    lis_psolve_saamg,  lis_psolve_iluc,
	lis_psolve_ilut_csr, lis_psolve_bjacobi, lis_psolve_adds
};

LIS_PSOLVEH_XXX lis_psolveh_xxx[] = {
	lis_psolveh_none,     lis_psolveh_jacobi, lis_psolveh_iluk_csr,
	lis_psolveh_ssor,     lis_psolveh_hybrid, lis_psolveh_none,
	lis_psolveh_sainv,    lis_psolveh_saamg,  lis_psolveh_iluc,
	lis_psolveh_ilut_csr, lis_psolveh_bjacobi, lis_psolveh_adds
};

LIS_PRECON_REGISTER *precon_register_top = NULL;

static LIS_PRECON_PSD_CREATE_XXX
        precon_register_psd_create[LIS_PRECON_REGISTER_MAX] = {NULL};

static LIS_PRECON_PSD_UPDATE_XXX
        precon_register_psd_update[LIS_PRECON_REGISTER_MAX] = {NULL};
LIS_INT	precon_register_type = LIS_PRECON_TYPE_USERDEF;

#if defined(USE_SAAMG)
	extern char *f_data_creation_ptr;
	extern char *f_v_cycle_ptr;
	extern char *f_clear_matrix_ptr;
#endif

LIS_REAL lis_precon_csr_row_scale(
    LIS_MATRIX A,
    LIS_INT row)
{
    LIS_INT j;
    LIS_REAL scale,t;

    scale = 0.0;

    if( A==NULL || row<0 || row>=A->n )
    {
        return scale;
    }

    if( A->is_splited &&
        A->D!=NULL &&
        A->L!=NULL &&
        A->U!=NULL )
    {
        scale = fabs(A->D->value[row]);

        for(j=A->L->ptr[row];j<A->L->ptr[row+1];j++)
        {

            t = fabs(A->L->value[j]);
            if( t>scale ) scale = t;
        }

        for(j=A->U->ptr[row];j<A->U->ptr[row+1];j++)
        {

            t = fabs(A->U->value[j]);
            if( t>scale ) scale = t;
        }
    }
    else
    {
        for(j=A->ptr[row];j<A->ptr[row+1];j++)
        {

            t = fabs(A->value[j]);
            if( t>scale ) scale = t;
        }
    }

    return scale;
}


#undef __FUNC__
#define __FUNC__ "lis_precon_check_ilu_pivot_tol"
LIS_INT lis_precon_check_ilu_pivot_tol(
    LIS_MATRIX A,
    LIS_REAL tol)
{
    LIS_INT i;
    LIS_REAL scale;

    if( tol<0.0 )
    {
        LIS_SETERR1(
            LIS_ERR_ILL_ARG,
            "ILU pivot tolerance (%e) is less than 0\n",
            (double)tol);
        return LIS_ERR_ILL_ARG;
    }

    if( tol==0.0 )
    {
        return LIS_SUCCESS;
    }

    for(i=0;i<A->n;i++)
    {
        scale = lis_precon_csr_row_scale(A,i);

        if( !(scale>0.0) )
        {
            LIS_SETERR1(
                LIS_BREAKDOWN,
                "zero row scale at row %D\n",
                i);
            return LIS_BREAKDOWN;
        }
    }

    return LIS_SUCCESS;
}


LIS_SCALAR lis_precon_regularize_ilu_pivot(
    LIS_SCALAR pivot,
    LIS_REAL row_scale,
    LIS_REAL tol)
{
    LIS_REAL mag,limit;

    if( tol<=0.0 )
    {
        return pivot;
    }

    limit = tol * row_scale;
    mag   = fabs(pivot);

    if( mag < limit )
    {
        if( mag>0.0 )
        {
            pivot *= (LIS_SCALAR)(limit/mag);
        }
        else
        {
            pivot = (LIS_SCALAR)limit;
        }
    }

    return pivot;
}


#undef __FUNC__
#define __FUNC__ "lis_precon_init"
LIS_INT lis_precon_init(LIS_PRECON precon)
{

	LIS_DEBUG_FUNC_IN;

	memset(precon,0,sizeof(struct LIS_PRECON_STRUCT));

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_set_user_data"
LIS_INT lis_precon_set_user_data(LIS_PRECON precon, void *user_data)
{
    LIS_DEBUG_FUNC_IN;

    if( precon==NULL )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"preconditioner is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    precon->user_data = user_data;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_get_user_data"
LIS_INT lis_precon_get_user_data(LIS_PRECON precon, void **user_data)
{
    LIS_DEBUG_FUNC_IN;

    if( precon==NULL )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"preconditioner is NULL\n");
        return LIS_ERR_ILL_ARG;
    }
    if( user_data==NULL )
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"user_data output pointer is NULL\n");
        return LIS_ERR_ILL_ARG;
    }

    *user_data = precon->user_data;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create"
LIS_INT lis_precon_create(LIS_SOLVER solver, LIS_PRECON *precon)
{
	LIS_INT	err;
	LIS_INT	precon_type;

	LIS_DEBUG_FUNC_IN;

	*precon     = NULL;
	precon_type = solver->options[LIS_OPTIONS_PRECON];

	*precon = (LIS_PRECON)lis_malloc( sizeof(struct LIS_PRECON_STRUCT),"lis_precon_create::precon" );
	if( NULL==*precon )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_PRECON_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_precon_init(*precon);
	(*precon)->precon_type = precon_type;

	if( precon_type>=LIS_PRECON_TYPE_USERDEF )
	{
		(*precon)->user_destroy =
		    precon_register_top[
		        precon_type-LIS_PRECON_TYPE_USERDEF].pdestroy;
		err = precon_register_top[precon_type-LIS_PRECON_TYPE_USERDEF].pcreate(solver,*precon);
	}
	else if( precon_type && solver->options[LIS_OPTIONS_ADDS] )
	{
		err = lis_precon_create_adds(solver,*precon);
		(*precon)->precon_type = LIS_PRECON_TYPE_ADDS;
	}
	else
	{
		err = lis_precon_create_xxx[precon_type](solver,*precon);
	}
	if( err )
	{
		lis_precon_destroy(*precon);
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_create"
LIS_INT lis_precon_psd_create(LIS_SOLVER solver, LIS_PRECON *precon)
{
	LIS_INT	err;
	LIS_INT	precon_type;

	LIS_DEBUG_FUNC_IN;

	*precon     = NULL;
	precon_type = solver->options[LIS_OPTIONS_PRECON];

	*precon = (LIS_PRECON)lis_malloc( sizeof(struct LIS_PRECON_STRUCT),"lis_precon_psd_create::precon" );
	if( NULL==*precon )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_PRECON_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_precon_init(*precon);
	(*precon)->precon_type = precon_type;

        if( precon_type>=LIS_PRECON_TYPE_USERDEF )
        {
                if( precon_register_top==NULL ||
                    precon_type>=precon_register_type )
                {
                        err = LIS_ERR_ILL_ARG;
                }
                else
                {
                        (*precon)->user_destroy =
                                precon_register_top[
                                    precon_type-LIS_PRECON_TYPE_USERDEF].pdestroy;

                        if( precon_register_psd_create[
                                precon_type-LIS_PRECON_TYPE_USERDEF] )
                        {
                                err = precon_register_psd_create[
                                        precon_type-LIS_PRECON_TYPE_USERDEF](
                                                solver,*precon);
                        }
                        else
                        {
                                err = LIS_ERR_NOT_IMPLEMENTED;
                        }
                }
        }
	else if( precon_type && solver->options[LIS_OPTIONS_ADDS] )
	{
/*        err = lis_precon_create_adds(solver,*precon);*/
/*        (*precon)->precon_type = LIS_PRECON_TYPE_ADDS;*/
        err = LIS_ERR_NOT_IMPLEMENTED;
	}
	else
	{
        switch(precon_type) {
            case LIS_PRECON_TYPE_JACOBI:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_SSOR:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_HYBRID:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_IS:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_SAI:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_ILUC:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_ILUT:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_BJACOBI:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            default:
                err = lis_precon_psd_create_xxx[precon_type](solver,*precon);
        }
	}
	if( err )
	{
		lis_precon_destroy(*precon);
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_precon_psd_update"
LIS_INT lis_precon_psd_update(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_INT	err;
	LIS_INT	precon_type;

	LIS_DEBUG_FUNC_IN;

    precon_type = precon->precon_type;

        if( precon_type>=LIS_PRECON_TYPE_USERDEF )
        {
                if( precon_register_top==NULL ||
                    precon_type>=precon_register_type )
                {
                        err = LIS_ERR_ILL_ARG;
                }
                else if( precon_register_psd_update[
                             precon_type-LIS_PRECON_TYPE_USERDEF] )
                {
                        err = precon_register_psd_update[
                                precon_type-LIS_PRECON_TYPE_USERDEF](
                                        solver,precon);
                }
                else
                {
                        err = LIS_ERR_NOT_IMPLEMENTED;
                }
        }
	else if( precon_type && solver->options[LIS_OPTIONS_ADDS] )
	{
/*        err = lis_precon_create_adds(solver,*precon);*/
/*        (*precon)->precon_type = LIS_PRECON_TYPE_ADDS;*/
        err = LIS_ERR_NOT_IMPLEMENTED;
	}
	else
	{
        switch(precon_type) {
            case LIS_PRECON_TYPE_NONE:
                /* its not that its not implemented, */
                /* but there is just nothing to do here . . . */
                err=0;
                break;
            case LIS_PRECON_TYPE_JACOBI:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_SSOR:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_HYBRID:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_IS:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_SAI:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_ILUC:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_ILUT:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            case LIS_PRECON_TYPE_BJACOBI:
                err = LIS_ERR_NOT_IMPLEMENTED;
                break;
            default:
                err = lis_precon_psd_update_xxx[precon_type](solver,precon);
        }
	}
	if( err )
	{
		lis_precon_destroy(precon);
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_none"
LIS_INT lis_precon_create_none(LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_DEBUG_FUNC_IN;

	precon->A       = solver->A;
	precon->is_copy = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_destroy"
LIS_INT lis_precon_destroy(LIS_PRECON precon)
{
	LIS_INT i, err = LIS_SUCCESS;

	LIS_DEBUG_FUNC_IN;

	if( precon )
	{
	    if( precon->user_destroy )
	    {
	        err = precon->user_destroy(precon);
	    }
		if( precon->is_copy ) lis_matrix_destroy(precon->A);
		lis_vector_destroy(precon->Pb);
		lis_vector_destroy(precon->D);
		lis_vector_destroy(precon->temp);
		lis_matrix_ilu_destroy(precon->L);
		lis_matrix_ilu_destroy(precon->U);
		lis_matrix_diag_destroy(precon->WD);
		if( precon->solver )
		{
			lis_vector_destroy(precon->solver->x);
			lis_precon_destroy(precon->solver->precon);
			lis_solver_destroy(precon->solver);
		}
#if defined(USE_SAAMG)
		lis_commtable_destroy(precon->commtable);
		if( precon->precon_type==LIS_PRECON_TYPE_SAAMG )
		{
		 ((void (*)(void *))f_clear_matrix_ptr)(&precon->level_num);
		}
#endif
		if( precon->work    )
		{
			for(i=0;i<precon->worklen;i++)
			{
				lis_vector_destroy(precon->work[i]);
			}
			lis_free(precon->work);
		}
		lis_free(precon);
	}

	LIS_DEBUG_FUNC_OUT;
	return err;
}


#undef __FUNC__
#define __FUNC__ "lis_psolve_dispatch"
LIS_INT lis_psolve_dispatch(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	LIS_INT precon_type;
	LIS_PSOLVE_XXX psolve;

	LIS_DEBUG_FUNC_IN;

	if( solver==NULL || solver->precon==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"preconditioner is not initialized\n");
		return LIS_ERR_ILL_ARG;
	}

	precon_type = solver->precon->precon_type;
	if( precon_type>=LIS_PRECON_TYPE_USERDEF )
	{
		if( precon_type>=precon_register_type || precon_register_top==NULL )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"invalid registered preconditioner type\n");
			return LIS_ERR_ILL_ARG;
		}
		psolve = precon_register_top[precon_type-LIS_PRECON_TYPE_USERDEF].psolve;
		if( psolve==NULL )
		{
			LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
			           "registered preconditioner has no psolve callback\n");
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		return psolve(solver,b,x);
	}

	if( precon_type<0 || precon_type>LIS_PRECON_TYPE_ADDS )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"invalid built-in preconditioner type\n");
		return LIS_ERR_ILL_ARG;
	}
	return lis_psolve_xxx[precon_type](solver,b,x);
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_dispatch"
LIS_INT lis_psolveh_dispatch(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	LIS_INT precon_type;
	LIS_PSOLVEH_XXX psolveh;

	LIS_DEBUG_FUNC_IN;

	if( solver==NULL || solver->precon==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"preconditioner is not initialized\n");
		return LIS_ERR_ILL_ARG;
	}

	precon_type = solver->precon->precon_type;
	if( precon_type>=LIS_PRECON_TYPE_USERDEF )
	{
		if( precon_type>=precon_register_type || precon_register_top==NULL )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"invalid registered preconditioner type\n");
			return LIS_ERR_ILL_ARG;
		}
		psolveh = precon_register_top[precon_type-LIS_PRECON_TYPE_USERDEF].psolveh;
		if( psolveh==NULL )
		{
			LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
			           "registered preconditioner has no psolveh callback\n");
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		return psolveh(solver,b,x);
	}

	if( precon_type<0 || precon_type>LIS_PRECON_TYPE_ADDS )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"invalid built-in preconditioner type\n");
		return LIS_ERR_ILL_ARG;
	}
	return lis_psolveh_xxx[precon_type](solver,b,x);
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_none"
LIS_INT lis_psolve_none(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	LIS_DEBUG_FUNC_IN;

	#ifndef USE_QUAD_PRECISION
		lis_vector_copy(b,x);
	#else
		if( solver->precision==LIS_PRECISION_DOUBLE )
		{
			lis_vector_copy(b,x);
		}
		else
		{
			lis_vector_copyex_mm(b,x);
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolveh_none"
LIS_INT lis_psolveh_none(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	LIS_DEBUG_FUNC_IN;

	#ifndef USE_QUAD_PRECISION
		lis_vector_copy(b,x);
	#else
		if( solver->precision==LIS_PRECISION_DOUBLE )
		{
			lis_vector_copy(b,x);
		}
		else
		{
			lis_vector_copyex_mm(b,x);
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_register"
LIS_INT lis_precon_register(char *name, LIS_PRECON_CREATE_XXX pcreate,
                            LIS_PSOLVE_XXX psolve, LIS_PSOLVEH_XXX psolveh)
{
	LIS_DEBUG_FUNC_IN;

	if( name==NULL || pcreate==NULL || psolve==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,
		           "registered preconditioner requires name, pcreate and psolve\n");
		return LIS_ERR_ILL_ARG;
	}

	if( precon_register_top==NULL )
	{
		precon_register_top = (LIS_PRECON_REGISTER *)lis_malloc(
			LIS_PRECON_REGISTER_MAX*sizeof(struct LIS_PRECON_REGISTER_STRUCT),
			"lis_precon_register::precon_register_top");
		if( precon_register_top==NULL )
		{
			LIS_SETERR_MEM(
				LIS_PRECON_REGISTER_MAX*sizeof(struct LIS_PRECON_REGISTER_STRUCT));
			return LIS_OUT_OF_MEMORY;
		}
	}

	if( precon_register_type-LIS_PRECON_TYPE_USERDEF==LIS_PRECON_REGISTER_MAX )
	{
		LIS_SETERR(LIS_FAILS,"lis_precon_resister is max\n");
		return LIS_FAILS;
	}

	precon_register_top[precon_register_type-LIS_PRECON_TYPE_USERDEF].pcreate = pcreate;
	precon_register_top[precon_register_type-LIS_PRECON_TYPE_USERDEF].psolve = psolve;
	precon_register_top[precon_register_type-LIS_PRECON_TYPE_USERDEF].psolveh = psolveh;
	precon_register_top[precon_register_type-LIS_PRECON_TYPE_USERDEF].pdestroy = NULL;
        precon_register_psd_create[
                precon_register_type-LIS_PRECON_TYPE_USERDEF] = NULL;
        precon_register_psd_update[
                precon_register_type-LIS_PRECON_TYPE_USERDEF] = NULL;
	precon_register_top[precon_register_type-LIS_PRECON_TYPE_USERDEF].precon_type =
		precon_register_type;
	strncpy(precon_register_top[precon_register_type-LIS_PRECON_TYPE_USERDEF].name,
	        name,LIS_PRECONNAME_MAX);
	precon_register_top[precon_register_type-LIS_PRECON_TYPE_USERDEF]
		.name[LIS_PRECONNAME_MAX] = '\0';
	precon_register_type++;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#undef __FUNC__
#define __FUNC__ "lis_precon_register_ex"
LIS_INT lis_precon_register_ex(char *name, LIS_PRECON_CREATE_XXX pcreate,
                               LIS_PSOLVE_XXX psolve, LIS_PSOLVEH_XXX psolveh,
                               LIS_PRECON_DESTROY_XXX pdestroy)
{
    LIS_INT err;

    LIS_DEBUG_FUNC_IN;

    err = lis_precon_register(name,pcreate,psolve,psolveh);
    if( err ) return err;

    precon_register_top[
        precon_register_type-LIS_PRECON_TYPE_USERDEF-1].pdestroy = pdestroy;

    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_register_psd"
LIS_INT lis_precon_register_psd(char *name,
                                LIS_PRECON_PSD_CREATE_XXX psd_create,
                                LIS_PRECON_PSD_UPDATE_XXX psd_update)
{
        LIS_INT i;
        LIS_INT nregister;

        LIS_DEBUG_FUNC_IN;

        if( name==NULL || psd_create==NULL || psd_update==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                           "PSD registration requires name, psd_create and psd_update\n");
                return LIS_ERR_ILL_ARG;
        }

        if( precon_register_top==NULL )
        {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                           "registered preconditioner not found\n");
                return LIS_ERR_ILL_ARG;
        }

        nregister = precon_register_type-LIS_PRECON_TYPE_USERDEF;

        for(i=0;i<nregister;i++)
        {
                if( strncmp(precon_register_top[i].name,
                            name,LIS_PRECONNAME_MAX)==0 )
                {
                        precon_register_psd_create[i] = psd_create;
                        precon_register_psd_update[i] = psd_update;

                        LIS_DEBUG_FUNC_OUT;
                        return LIS_SUCCESS;
                }
        }

        LIS_SETERR(LIS_ERR_ILL_ARG,
                   "registered preconditioner not found\n");
        return LIS_ERR_ILL_ARG;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_register_free"
LIS_INT lis_precon_register_free(void)
{

	LIS_DEBUG_FUNC_IN;

	if( precon_register_top )
	{
		lis_free(precon_register_top);
		precon_register_top = NULL;
	}
	precon_register_type = LIS_PRECON_TYPE_USERDEF;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
