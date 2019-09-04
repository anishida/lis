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
/*NEH deleted reference to "math.h" from lis.h, and placed here,*/
/*NEH for consistency . . .*/
#include <math.h>
#include <ctype.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_solver_init
 * lis_solver_create
 * lis_solver_destroy
 * lis_solver_work_destroy
 * lis_solver_set_option
 * lis_solver_get_option
 * lis_solve
 ************************************************/

#define LIS_SOLVERS_LEN			25
#define LIS_PRECON_TYPE_LEN		12


LIS_SOLVER_CHECK_PARAMS lis_solver_check_params[] = {
	NULL,
	lis_cg_check_params       , lis_bicg_check_params      , lis_cgs_check_params, 
	lis_bicgstab_check_params , lis_bicgstabl_check_params , lis_gpbicg_check_params, 
	lis_tfqmr_check_params    , lis_orthomin_check_params  , lis_gmres_check_params,
	lis_jacobi_check_params   , lis_gs_check_params        , lis_sor_check_params,
	lis_bicgsafe_check_params , lis_cr_check_params        , lis_bicr_check_params,
	lis_crs_check_params      , lis_bicrstab_check_params  , lis_gpbicr_check_params,
	lis_bicrsafe_check_params , lis_fgmres_check_params    , lis_idrs_check_params,
	lis_idr1_check_params     , lis_minres_check_params    , lis_cocg_check_params,
	lis_cocr_check_params
};

LIS_SOLVER_MALLOC_WORK lis_solver_malloc_work[] = {
	NULL,
	lis_cg_malloc_work       , lis_bicg_malloc_work      , lis_cgs_malloc_work, 
	lis_bicgstab_malloc_work , lis_bicgstabl_malloc_work , lis_gpbicg_malloc_work, 
	lis_tfqmr_malloc_work    , lis_orthomin_malloc_work  , lis_gmres_malloc_work,
	lis_jacobi_malloc_work   , lis_gs_malloc_work        , lis_sor_malloc_work,
	lis_bicgsafe_malloc_work , lis_cr_malloc_work        , lis_bicr_malloc_work,
	lis_crs_malloc_work      , lis_bicrstab_malloc_work  , lis_gpbicr_malloc_work,
	lis_bicrsafe_malloc_work , lis_fgmres_malloc_work    , lis_idrs_malloc_work,
	lis_idr1_malloc_work     , lis_minres_malloc_work    , lis_cocg_malloc_work,
        lis_cocr_malloc_work
};

LIS_SOLVER_EXECUTE lis_solver_execute[] = {
	NULL,
	lis_cg       , lis_bicg      , lis_cgs, 
	lis_bicgstab , lis_bicgstabl , lis_gpbicg, 
	lis_tfqmr    , lis_orthomin  , lis_gmres,
	lis_jacobi   , lis_gs        , lis_sor,
	lis_bicgsafe , lis_cr        , lis_bicr,
	lis_crs      , lis_bicrstab  , lis_gpbicr,
	lis_bicrsafe , lis_fgmres    , lis_idrs, 
	lis_idr1     , lis_minres    , lis_cocg,
	lis_cocr
};

LIS_SOLVER_EXECUTE lis_solver_execute_conv_cond[] = {
	NULL,
	lis_cg       , lis_bicg      , lis_cgs, 
	lis_bicgstab , lis_bicgstabl , lis_gpbicg, 
	NULL         , lis_orthomin  , NULL,
	NULL         , NULL          , NULL,
	lis_bicgsafe , lis_cr        , lis_bicr,
	lis_crs      , lis_bicrstab  , lis_gpbicr,
	lis_bicrsafe , NULL          , lis_idrs, 
	lis_idr1     , NULL          , lis_cocg,
	lis_cocr
};

#ifdef USE_QUAD_PRECISION
	LIS_SOLVER_EXECUTE lis_solver_execute_quad[] = {
		NULL,
		lis_cg_quad       , lis_bicg_quad      , lis_cgs_quad, 
		lis_bicgstab_quad , lis_bicgstabl_quad , lis_gpbicg_quad,
		lis_tfqmr_quad    , lis_orthomin_quad  , lis_gmres_quad,
		NULL              , NULL               , NULL,
		lis_bicgsafe_quad , lis_cr_quad        , lis_bicr_quad,
		lis_crs_quad      , lis_bicrstab_quad  , lis_gpbicr_quad,
		lis_bicrsafe_quad , lis_fgmres_quad    , NULL,
		NULL              , NULL               , NULL,
		NULL
	};
	LIS_SOLVER_EXECUTE lis_solver_execute_switch[] = {
		NULL,
		lis_cg_switch       , lis_bicg_switch , lis_cgs_switch, 
		lis_bicgstab_switch , NULL            , lis_gpbicg_switch,
		NULL                , NULL            , lis_gmres_switch,
		NULL                , NULL            , NULL,
		NULL                , NULL            , NULL,
		NULL                , NULL            , NULL,
		NULL                , NULL            , NULL,
		NULL                , NULL            , NULL,
		NULL
	};
	/*
	LIS_SOLVER_EXECUTE lis_solver_execute_periodic[] = {
		NULL,
		lis_cg_periodic       , lis_bicg_periodic , lis_cgs_periodic, 
		lis_bicgstab_periodic , NULL              , lis_gpbicg_periodic,
		NULL                  , NULL              , NULL,
		NULL                  , NULL              , NULL,
		NULL
	};
	*/
#endif

LIS_SOLVER_GET_RESIDUAL lis_solver_get_residual[] = {
	lis_solver_get_residual_nrm2_r, 
	lis_solver_get_residual_nrm2_r, 
	lis_solver_get_residual_nrm1_b
};

LIS_INT LIS_USE_AT_TYPE[] = {
	0,
	LIS_MATRIX_CSC,LIS_MATRIX_CSR
	};
#define LIS_SOLVER_OPTION_LEN		46
#define LIS_PRINT_LEN			4
#define LIS_SCALE_LEN			3
#define LIS_TRUEFALSE_LEN		2
#define LIS_PRECISION_LEN		3
#define LIS_STORAGE_LEN			11
#define LIS_CONV_COND_LEN		3

char *LIS_SOLVER_OPTNAME[] = {
	"-maxiter",           "-tol",           "-print",          "-scale",         "-ssor_omega",
	"-ilu_fill",          "-ilu_relax",     "-is_alpha",       "-is_level",      "-is_m",
	"-hybrid_maxiter",    "-hybrid_ell",    "-hybrid_restart", "-hybrid_tol",    "-hybrid_omega",
	"-hybrid_i",          "-sainv_drop",    "-ric2s_tau",      "-ric2s_sigma",   "-ric2s_gamma",
	"-restart",           "-ell",           "-omega",          "-i",             "-p",
	"-f",                 "-h",             "-ver",            "-hybrid_p",      "-initx_zeros",
	"-adds",              "-adds_iter",     "-f",              "-use_at",        "-switch_tol",
	"-switch_maxiter",    "-saamg_unsym",   "-iluc_drop",      "-iluc_gamma",    "-iluc_rate",
	"-storage",           "-storage_block", "-conv_cond",      "-tol_w",         "-saamg_theta",	"-irestart"
};

LIS_INT LIS_SOLVER_OPTACT[] = {
	LIS_OPTIONS_MAXITER          , LIS_PARAMS_RESID          , LIS_OPTIONS_OUTPUT        , LIS_OPTIONS_SCALE        , LIS_PARAMS_SSOR_OMEGA ,
	LIS_OPTIONS_FILL             , LIS_PARAMS_RELAX          , LIS_PARAMS_ALPHA          , LIS_OPTIONS_ISLEVEL      , LIS_OPTIONS_M     ,
	LIS_OPTIONS_PMAXITER         , LIS_OPTIONS_PELL          , LIS_OPTIONS_PRESTART      , LIS_PARAMS_PRESID        , LIS_PARAMS_POMEGA     ,
	LIS_OPTIONS_PSOLVER          , LIS_PARAMS_DROP           , LIS_PARAMS_TAU            , LIS_PARAMS_SIGMA         , LIS_PARAMS_GAMMA  ,
	LIS_OPTIONS_RESTART          , LIS_OPTIONS_ELL           , LIS_PARAMS_OMEGA          , LIS_OPTIONS_SOLVER       , LIS_OPTIONS_PRECON,
	LIS_OPTIONS_FILE             , LIS_OPTIONS_HELP          , LIS_OPTIONS_VER           , LIS_OPTIONS_PPRECON      , LIS_OPTIONS_INITGUESS_ZEROS,
	LIS_OPTIONS_ADDS             , LIS_OPTIONS_ADDS_ITER     , LIS_OPTIONS_PRECISION     , LIS_OPTIONS_USE_AT       , LIS_PARAMS_SWITCH_RESID,
	LIS_OPTIONS_SWITCH_MAXITER   , LIS_OPTIONS_SAAMG_UNSYM   , LIS_PARAMS_DROP           , LIS_PARAMS_GAMMA         , LIS_PARAMS_RATE, 
	LIS_OPTIONS_STORAGE          , LIS_OPTIONS_STORAGE_BLOCK , LIS_OPTIONS_CONV_COND     , LIS_PARAMS_RESID_WEIGHT  , LIS_PARAMS_SAAMG_THETA, LIS_OPTIONS_IDRS_RESTART
};

char *lis_solver_atoi[]    = {"cg", "bicg", "cgs", "bicgstab", "bicgstabl", "gpbicg", "tfqmr","orthomin", "gmres", "jacobi", "gs", "sor", "bicgsafe", "cr", "bicr", "crs", "bicrstab", "gpbicr", "bicrsafe", "fgmres", "idrs", "idr1", "minres", "cocg", "cocr"};
char *lis_precon_atoi[]    = {"none", "jacobi", "ilu", "ssor", "hybrid", "is", "sainv", "saamg", "iluc", "ilut", "bjacobi", ""};
char *lis_storage_atoi[]   = {"csr", "csc", "msr", "dia", "ell", "jad", "bsr", "bsc", "vbr", "coo", "dns"};
char *lis_print_atoi[]     = {"none", "mem", "out", "all"};
char *lis_scale_atoi[]     = {"none", "jacobi", "symm_diag"};
char *lis_truefalse_atoi[] = {"false", "true"};
char *lis_precision_atoi[] = {"double", "quad", "switch"};
char *lis_conv_cond_atoi[] = {"nrm2_r", "nrm2_b", "nrm1_b"};

char *lis_solvername[] = {"", "CG", "BiCG", "CGS", "BiCGSTAB", "BiCGSTAB(l)", "GPBiCG", "TFQMR", "Orthomin", "GMRES", "Jacobi",	"Gauss-Seidel", "SOR", "BiCGSafe", "CR", "BiCR", "CRS", "BiCRSTAB", "GPBiCR", "BiCRSafe", "FGMRES", "IDR(s)", "IDR(1)", "MINRES", "COCG", "COCR"};
char *lis_preconname[] = {"none", "Jacobi", "ILU", "SSOR", "Hybrid", "I+S", "SAINV", "SAAMG", "Crout ILU", "ILUT", "Block Jacobi"};

char *lis_returncode[] = {"LIS_SUCCESS", "LIS_ILL_OPTION", "LIS_BREAKDOWN", "LIS_OUT_OF_MEMORY", "LIS_MAXITER", "LIS_NOT_IMPLEMENTED", "LIS_ERR_FILE_IO"};
char *lis_precisionname[] = {"double", "quad", "switch"};
char *lis_storagename[]   = {"CSR", "CSC", "MSR", "DIA", "ELL", "JAD", "BSR", "BSC", "VBR", "COO", "DNS"};

LIS_VECTOR lis_solver_residual_history = NULL;

#undef __FUNC__
#define __FUNC__ "lis_solver_init"
LIS_INT lis_solver_init(LIS_SOLVER solver)
{

	LIS_DEBUG_FUNC_IN;

	solver->A        = NULL;
	solver->Ah       = NULL;
	solver->b        = NULL;
	solver->x        = NULL;
	solver->d        = NULL;
	solver->work     = NULL;
	solver->rhistory = NULL;
	solver->precon   = NULL;

	solver->worklen   = 0;
	solver->iter      = 0;
	solver->iter2     = 0;
	solver->resid     = 0;
	solver->time      = 0;
	solver->itime     = 0;
	solver->ptime     = 0;
	solver->precision = LIS_PRECISION_DOUBLE;

	solver->options[LIS_OPTIONS_SOLVER]               = LIS_SOLVER_BICG;
	solver->options[LIS_OPTIONS_PRECON]               = LIS_PRECON_TYPE_NONE;
	solver->options[LIS_OPTIONS_OUTPUT]               = LIS_FALSE;
	solver->options[LIS_OPTIONS_MAXITER]              = 1000;
	solver->options[LIS_OPTIONS_RESTART]              = 40;
	solver->options[LIS_OPTIONS_ELL]                  = 2;
	solver->options[LIS_OPTIONS_SCALE]                = LIS_SCALE_NONE;
	solver->options[LIS_OPTIONS_FILL]                 = 0;
	solver->options[LIS_OPTIONS_M]                    = 3;
	solver->options[LIS_OPTIONS_PSOLVER]              = LIS_SOLVER_SOR;
	solver->options[LIS_OPTIONS_PMAXITER]             = 25;
	solver->options[LIS_OPTIONS_PRESTART]             = 40;
	solver->options[LIS_OPTIONS_PELL]                 = 2;
	solver->options[LIS_OPTIONS_PPRECON]              = LIS_PRECON_TYPE_NONE;
	solver->options[LIS_OPTIONS_ISLEVEL]              = 1;
	solver->options[LIS_OPTIONS_INITGUESS_ZEROS]      = LIS_TRUE;
	solver->options[LIS_OPTIONS_ADDS]                 = LIS_FALSE;
	solver->options[LIS_OPTIONS_ADDS_ITER]            = 1;
	solver->options[LIS_OPTIONS_PRECISION]            = LIS_PRECISION_DOUBLE;
	solver->options[LIS_OPTIONS_USE_AT]               = LIS_FALSE;
	solver->options[LIS_OPTIONS_SWITCH_MAXITER]       = -1;
	solver->options[LIS_OPTIONS_SAAMG_UNSYM]          = LIS_FALSE;
	solver->options[LIS_OPTIONS_STORAGE]              = 0;
	solver->options[LIS_OPTIONS_STORAGE_BLOCK]        = 2;
	solver->options[LIS_OPTIONS_CONV_COND]            = 0;
	solver->options[LIS_OPTIONS_INIT_SHADOW_RESID]    = LIS_RESID;
	solver->options[LIS_OPTIONS_IDRS_RESTART]         = 2;

	solver->params[LIS_PARAMS_RESID        -LIS_OPTIONS_LEN] = 1.0e-12;
	solver->params[LIS_PARAMS_RESID_WEIGHT -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_OMEGA        -LIS_OPTIONS_LEN] = 1.9;
	solver->params[LIS_PARAMS_SSOR_OMEGA   -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_RELAX        -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_DROP         -LIS_OPTIONS_LEN] = 0.05;
	solver->params[LIS_PARAMS_ALPHA        -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_TAU          -LIS_OPTIONS_LEN] = 0.05;
	solver->params[LIS_PARAMS_SIGMA        -LIS_OPTIONS_LEN] = 2.0;
	solver->params[LIS_PARAMS_GAMMA        -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_PRESID       -LIS_OPTIONS_LEN] = 1.0e-3;
	solver->params[LIS_PARAMS_POMEGA       -LIS_OPTIONS_LEN] = 1.5;
	solver->params[LIS_PARAMS_SWITCH_RESID -LIS_OPTIONS_LEN] = 1.0e-12;
	solver->params[LIS_PARAMS_RATE         -LIS_OPTIONS_LEN] = 5.0;
	solver->params[LIS_PARAMS_SAAMG_THETA  -LIS_OPTIONS_LEN] = 0.05;

	/* reset solver->setup */
	solver->setup = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_create"
LIS_INT lis_solver_create(LIS_SOLVER *solver)
{
	LIS_DEBUG_FUNC_IN;

	*solver = NULL;

	*solver = (LIS_SOLVER)lis_malloc( sizeof(struct LIS_SOLVER_STRUCT),"lis_solver_create::solver" );
	if( NULL==*solver )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_SOLVER_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_solver_init(*solver);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_work_destroy"
LIS_INT lis_solver_work_destroy(LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( solver && solver->work )
	{
		for(i=0;i<solver->worklen;i++) lis_vector_destroy(solver->work[i]);
		lis_free(solver->work);
		solver->work    = NULL;
		solver->worklen = 0;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_destroy"
LIS_INT lis_solver_destroy(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;

	if( solver )
	{
		lis_solver_work_destroy(solver);
		lis_vector_destroy(solver->d);
		if( solver->Ah ) lis_matrix_destroy(solver->Ah);
		if( solver->rhistory ) lis_free(solver->rhistory);
		lis_free(solver);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_solver_set_matrix"
LIS_INT lis_solver_set_matrix(LIS_MATRIX A, LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;

	solver->A = A;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solve"
LIS_INT lis_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_SOLVER solver)
{
    LIS_INT	err;
	LIS_PRECON precon;

	LIS_DEBUG_FUNC_IN;

	solver->A = A;

	/* create preconditioner */

	if( solver->options[LIS_OPTIONS_PRECON] < 0 || solver->options[LIS_OPTIONS_PRECON] > LIS_PRECONNAME_MAX )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PRECON is %D (Set between 0 to %D)\n",solver->options[LIS_OPTIONS_PRECON], LIS_PRECONNAME_MAX);
		return LIS_ERR_ILL_ARG;
	}

	err = lis_precon_create(solver, &precon);
	if( err )
	{
		lis_solver_work_destroy(solver);
		solver->retcode = err;
		return err;
	}

	/* core kernel of lis_solve() */
	err = lis_solve_kernel(A, b, x, solver, precon);
	if( err )
	{
		lis_solver_work_destroy(solver);	  
		solver->retcode = err;
		return err;
	}

	lis_precon_destroy(precon);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solve_setup"
LIS_INT lis_solve_setup(LIS_MATRIX A, LIS_SOLVER solver)
{
        LIS_INT	err;
	LIS_PRECON precon;
	LIS_VECTOR b,x;

	LIS_DEBUG_FUNC_IN;

	lis_vector_duplicate(A,&b);
	lis_vector_duplicate(A,&x);

	/* setup solver for preconditioning */
	/* Do not call lis_solve_execute if solver->setup is true. 
	   See esolver/lis_esolver_cg.c, where only preconditioner is called. */
	solver->setup = LIS_TRUE;
	err = lis_solve(A,b,x,solver);
	if( err )
	  {
	    lis_solver_work_destroy(solver);
	    solver->retcode = err;
	    return err;
	  }

	lis_vector_destroy(b);
	lis_vector_destroy(x);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solve_kernel"
LIS_INT lis_solve_kernel(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_SOLVER solver, LIS_PRECON precon)
{
	LIS_Comm comm;  
	LIS_INT	nsolver, precon_type, maxiter;
	LIS_INT	err;
	LIS_REAL *rhistory;
	LIS_VECTOR xx;

	LIS_INT output;
	LIS_INT scale;
	LIS_INT conv_cond;
	LIS_INT precision,is_use_at,storage,block;
	LIS_INT i,n;
	double p_c_time, p_i_time,itime;
	LIS_REAL nrm2,tol,tol_w;
	LIS_VECTOR t;
	LIS_VECTOR bb;
	LIS_MATRIX AA,B;
	LIS_MATRIX Ah;
	char buf[64];

	LIS_VECTOR r,z;

	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	nsolver     = solver->options[LIS_OPTIONS_SOLVER];
	precon_type = solver->options[LIS_OPTIONS_PRECON];
	maxiter     = solver->options[LIS_OPTIONS_MAXITER];
	output      = solver->options[LIS_OPTIONS_OUTPUT];
	scale       = solver->options[LIS_OPTIONS_SCALE];
	precision   = solver->options[LIS_OPTIONS_PRECISION];
	is_use_at   = solver->options[LIS_OPTIONS_USE_AT];
	storage     = solver->options[LIS_OPTIONS_STORAGE];
	block       = solver->options[LIS_OPTIONS_STORAGE_BLOCK];
	conv_cond   = solver->options[LIS_OPTIONS_CONV_COND];
	tol         = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol_w       = solver->params[LIS_PARAMS_RESID_WEIGHT-LIS_OPTIONS_LEN];
	solver->precision = precision;

	if( nsolver < 1 || nsolver > LIS_SOLVERS_LEN )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SOLVER is %D (Set between 1 to %D)\n",nsolver, LIS_SOLVERS_LEN);
		return LIS_ERR_ILL_ARG;
	}
	if( precon_type < 0 || precon_type > precon_register_type )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PRECON is %D (Set between 0 to %D)\n",precon_type, precon_register_type-1);
		return LIS_ERR_ILL_ARG;
	}
	if( maxiter<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_MAXITER(=%D) is less than 0\n",maxiter);
		return LIS_ERR_ILL_ARG;
	}
	if( conv_cond>0 && lis_solver_execute_conv_cond[nsolver]==NULL )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Option conv_cond is not implemented for solver %s\n",lis_solvername[nsolver]);
		return LIS_ERR_ILL_ARG;
	}
	#ifdef USE_MPI
	if( precon_type == LIS_PRECON_TYPE_SAAMG  && solver->A->nprocs < 2)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter A->nprocs (=%D) is less than 2 (Set more than 1 when using parallel version of SAAMG)\n",solver->A->nprocs);
		return LIS_ERR_ILL_ARG;
	}
	#endif
	#ifdef USE_QUAD_PRECISION
		if( precision==LIS_PRECISION_QUAD && lis_solver_execute_quad[nsolver]==NULL )
		{
			LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Quad precision solver %s is not implemented\n",lis_solvername[nsolver]);
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		else if( precision==LIS_PRECISION_SWITCH && lis_solver_execute_switch[nsolver]==NULL )
		{
			LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Switch solver %s is not implemented\n",lis_solvername[nsolver]);
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		if( solver->options[LIS_OPTIONS_SWITCH_MAXITER]==-1 )
		{
			solver->options[LIS_OPTIONS_SWITCH_MAXITER] = maxiter;
		}
	#else		
		if( precision==LIS_PRECISION_QUAD )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"Quad precision is not enabled\n");
			return LIS_ERR_ILL_ARG;
		}
	#endif

	err = lis_solver_check_params[nsolver](solver);
	if( err )
	{
		solver->retcode = err;
		return err;
	}
	/* end parameter check */

	solver->A        = A;
	solver->b        = b;

	/* create initial vector */
	#ifndef USE_QUAD_PRECISION
		err = lis_vector_duplicate(A,&xx);
	#else
		if( precision==LIS_PRECISION_DOUBLE )
		{
			err = lis_vector_duplicate(A,&xx);
		}
		else
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,A,&xx);
		}
	#endif
	if( err )
	{
		solver->retcode = err;
		return err;
	}
	if( solver->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
	  if( output ) lis_printf(comm,"initial vector x      : all components set to 0\n");
		#ifndef USE_QUAD_PRECISION
	                lis_vector_set_all(0.0,xx);
		#else
			if( precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_set_all(0.0,xx);
			}
			else
			{
				lis_vector_set_allex_nm(0.0,xx);
			}
		#endif
	}
	else
	{
	  if( output ) lis_printf(comm,"initial vector x      : user defined\n"); 
		#ifndef USE_QUAD_PRECISION
			lis_vector_copy(x,xx);
		#else
			if( precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_copy(x,xx);
			}
			else
			{
				lis_vector_copyex_nm(x,xx);
			}
		#endif
	}

	/* create residual history vector */
	if( solver->rhistory ) lis_free(solver->rhistory);
	rhistory = (LIS_REAL *)lis_malloc((maxiter+2)*sizeof(LIS_REAL),"lis_solve::rhistory");
	if( rhistory==NULL )
	{
		LIS_SETERR_MEM((maxiter+2)*sizeof(LIS_SCALAR));
		lis_vector_destroy(xx);
		solver->retcode = err;
		return err;
	}
	rhistory[0] = 1.0;


	n       = A->n;
	t       = NULL;
	Ah      = NULL;


	p_c_time = lis_wtime();
	if( precon_type==LIS_PRECON_TYPE_IS )
	{
		if( solver->d==NULL )
		{
			err = lis_vector_duplicate(A,&solver->d);
			if( err )
			{
				return err;
			}
		}
		if( !A->is_scaled )
		{
			lis_matrix_scale(A,b,solver->d,LIS_SCALE_JACOBI);
		}
		else if( !b->is_scaled )
		{
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(i=0;i<n;i++)
			{
				b->value[i] = b->value[i]*solver->d->value[i];
			}
		}
		if( nsolver >= LIS_SOLVER_JACOBI && nsolver <= LIS_SOLVER_SOR )
		{
			solver->options[LIS_OPTIONS_ISLEVEL] = 0;
		}
	}
	else if( nsolver >= LIS_SOLVER_JACOBI && nsolver <= LIS_SOLVER_SOR && precon_type!=LIS_PRECON_TYPE_NONE )
	{
		if( solver->d==NULL )
		{
			err = lis_vector_duplicate(A,&solver->d);
			if( err )
			{
				return err;
			}
		}
		if( !A->is_scaled )
		{
			lis_matrix_scale(A,b,solver->d,LIS_SCALE_JACOBI);
		}
	}
	else if( scale )
	{
		if( storage==LIS_MATRIX_BSR && scale==LIS_SCALE_JACOBI )
		{
			if( A->matrix_type!=LIS_MATRIX_BSR )
			{
				err = lis_matrix_duplicate(A,&B);
				if( err ) return err;
				lis_matrix_set_blocksize(B,block,block,NULL,NULL);
				lis_matrix_set_type(B,storage);
				err = lis_matrix_convert(A,B);
				if( err ) return err;
				lis_matrix_storage_destroy(A);
				lis_matrix_DLU_destroy(A);
				lis_matrix_diag_destroy(A->WD);
				if( A->l2g_map ) lis_free( A->l2g_map );
				if( A->commtable ) lis_commtable_destroy( A->commtable );
				if( A->ranges ) lis_free( A->ranges );
				err = lis_matrix_copy_struct(B,A);
				if( err ) return err;
				lis_free(B);
			}
			err = lis_matrix_split(A);
			if( err ) return err;
			err = lis_matrix_diag_duplicate(A->D,&solver->WD);
			if( err ) return err;
			lis_matrix_diag_copy(A->D,solver->WD);
			lis_matrix_diag_inverse(solver->WD);
			lis_matrix_bscale_bsr(A,solver->WD);
			lis_vector_duplicate(A,&t);
			lis_matrix_diag_matvec(solver->WD,b,t);
			lis_vector_copy(t,b);
			lis_vector_destroy(t);
			t = NULL;
		}
		else
		{
			if( solver->d==NULL )
			{
				err = lis_vector_duplicate(A,&solver->d);
				if( err )
				{
					return err;
				}
			}
			if( scale==LIS_SCALE_JACOBI && nsolver==LIS_SOLVER_CG )
			{
				scale = LIS_SCALE_SYMM_DIAG;
			}
			if( !A->is_scaled )
			{
				lis_matrix_scale(A,b,solver->d,scale);
			}
			else if( !b->is_scaled )
			{
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(i=0;i<n;i++)
				{
					b->value[i] = b->value[i]*solver->d->value[i];
				}
			}
		}
	}

/*	precon_type = precon->precon_type;*/
	if( precon_type==LIS_PRECON_TYPE_IS )
	{
		if( nsolver < LIS_SOLVER_JACOBI || nsolver > LIS_SOLVER_SOR )
		{
			AA = solver->A;
			bb = solver->b;
		}
		else
		{
			AA = precon->A;
			bb = precon->Pb;
		}
	}
	else
	{
		AA = A;
		bb = b;
	}

	p_c_time = lis_wtime() - p_c_time;
	itime = lis_wtime();

	/* convert matrix */
	solver->A  = AA;
	solver->b  = bb;
	err = lis_matrix_convert_self(solver);
	if( err )
	{
		lis_vector_destroy(xx);
		lis_solver_work_destroy(solver);
		lis_free(rhistory);
		solver->retcode = err;
		return err;
	}
	block = solver->A->bnr;

	if( A->my_rank==0 )
	{
#ifdef _LONG__DOUBLE
	  if( output ) printf("precision             : long double\n");
#else
 	  if( output ) printf("precision             : %s\n", lis_precisionname[precision]); 
#endif
	  if( output ) printf("linear solver         : %s\n", lis_solvername[nsolver]); 
		switch( precon_type )
		{
		case LIS_PRECON_TYPE_ILU:
			i = solver->options[LIS_OPTIONS_FILL];
			if( A->matrix_type==LIS_MATRIX_BSR || A->matrix_type==LIS_MATRIX_VBR )
			{
#ifdef _LONG__LONG
			  if( output ) sprintf(buf,"Block %s(%lld)",lis_preconname[precon_type],i); 
#else
			  if( output ) sprintf(buf,"Block %s(%d)",lis_preconname[precon_type],i); 
#endif
			}
			else
			{
#ifdef _LONG__LONG
			  if( output ) sprintf(buf,"%s(%lld)",lis_preconname[precon_type],i); 
#else
			  if( output ) sprintf(buf,"%s(%d)",lis_preconname[precon_type],i); 
#endif
			}
			break;
		default:
		  if( output ) sprintf(buf,"%s",lis_preconname[precon_type]); 
			break;
		}
		if( solver->options[LIS_OPTIONS_ADDS] && precon_type )
		{
 		  if( output ) printf("preconditioner        : %s + Additive Schwarz\n", buf);
		}
		else
		{
		  if( output ) printf("preconditioner        : %s\n", buf); 
		}
	}
	switch(conv_cond)
	{
	case LIS_CONV_COND_NRM2_R:
	  if( output ) lis_printf(comm,"convergence condition : ||b-Ax||_2 <= %6.1e * ||b-Ax_0||_2\n", (double)tol); 		  
		break;		
	case LIS_CONV_COND_NRM2_B:
		lis_vector_nrm2(b,&nrm2);
		nrm2 = nrm2*tol;
		if( output ) lis_printf(comm,"convergence condition : ||b-Ax||_2 <= %6.1e*||b||_2 = %6.1e\n", (double)tol,(double)nrm2);
		break;
	case LIS_CONV_COND_NRM1_B:
		lis_vector_nrm1(b,&nrm2);
		nrm2 = nrm2*tol_w + tol;
		if( output ) lis_printf(comm,"convergence condition : ||b-Ax||_1 <= %6.1e*||b||_1 + %6.1e = %6.1e\n", (double)tol_w,(double)tol,(double)nrm2);
		break;
	}
	if( AA->matrix_type==LIS_MATRIX_BSR || AA->matrix_type==LIS_MATRIX_BSC )
	  {
	    if( output ) lis_printf(comm,"matrix storage format : %s(%D x %D)\n", lis_storagename[AA->matrix_type-1],block,block); 
	  }
	else
	  {
	    if( output ) lis_printf(comm,"matrix storage format : %s\n", lis_storagename[AA->matrix_type-1]); 
	  }

	/* create work vector */
	err = lis_solver_malloc_work[nsolver](solver); 
	if( err )
	{
		lis_vector_destroy(xx);
		lis_precon_destroy(precon);
		solver->retcode = err;
		return err;
	}
	if( nsolver==LIS_SOLVER_BICG && is_use_at )
	{
	  if( output ) lis_printf(comm,"Use Ah\n"); 
	  lis_matrix_duplicate(AA,&Ah);
	  lis_matrix_set_type(Ah,LIS_USE_AT_TYPE[AA->matrix_type]);
	  lis_matrix_convert(AA,Ah);
	  solver->Ah = Ah;
	}

	solver->x        = xx;
	solver->xx       = x;
	solver->precon   = precon;
	solver->rhistory = rhistory;

	/* Do not call lis_solve_execute if solver->setup is true. 
	   See esolver/lis_esolver_cg.c, where only preconditioner is called.
	   solver->setup is initialized in lis_solver_init, 
	   and reset by lis_solve_setup.
	*/
	
	if (!solver->setup)
	  {
#ifndef USE_QUAD_PRECISION
	    err = lis_solver_execute[nsolver](solver);
#else
	    if( precision==LIS_PRECISION_DOUBLE )
	      {
		err = lis_solver_execute[nsolver](solver);
	      }
	    else if( precision==LIS_PRECISION_QUAD )
	      {
		err = lis_solver_execute_quad[nsolver](solver);
	      }
	    else if( precision==LIS_PRECISION_SWITCH )
	      {
		err = lis_solver_execute_switch[nsolver](solver);
	      }
#endif
	    solver->retcode = err;
	  }

	if( scale==LIS_SCALE_SYMM_DIAG && precon_type!=LIS_PRECON_TYPE_IS)
	{
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(i=0;i<n;i++)
		{
			x->value[i] = xx->value[i]*solver->d->value[i];
		}
	}
	else
	{
		#ifndef USE_QUAD_PRECISION
			lis_vector_copy(xx,x);
		#else
			if( precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_copy(xx,x);
			}
			else
			{
				lis_vector_copyex_mn(xx,x);
			}
		#endif
	}
	itime = lis_wtime() - itime - solver->ptime;
	p_i_time = solver->ptime;
	solver->ptime = p_c_time + p_i_time;
	solver->p_c_time = p_c_time;
	solver->p_i_time = p_i_time;
	solver->time  = solver->ptime + itime;
	solver->itime = itime;
	lis_solver_work_destroy(solver);
	lis_vector_duplicate(A,&t);
	xx->precision = LIS_PRECISION_DEFAULT;
	lis_matvec(A,xx,t);
	lis_vector_xpay(b,-1.0,t);
	if( scale==LIS_SCALE_SYMM_DIAG && precon_type!=LIS_PRECON_TYPE_IS)
	{
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(i=0;i<n;i++)
		{
			t->value[i] = t->value[i]/solver->d->value[i];
		}
	}
	lis_vector_nrm2(t,&nrm2);

	/* solver->resid = nrm2; */
	if( err )
	  {
	    if( output ) lis_printf(comm,"linear solver status  : %s(code=%D)\n\n",lis_returncode[err],err); 
	  }
	else
	  {
	    if( output ) lis_printf(comm,"linear solver status  : normal end\n\n"); 
	  }

	if( precision==LIS_PRECISION_DOUBLE )
	{
		solver->iter2 = solver->iter;
	}
	else if( precision==LIS_PRECISION_QUAD )
	{
		solver->iter2 = 0;
	}


	lis_vector_destroy(t);

	/* lis_vector_destroy(d); */
	lis_vector_destroy(xx);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_initial_residual"
LIS_INT lis_solver_get_initial_residual(LIS_SOLVER solver, LIS_PRECON M, LIS_VECTOR t, LIS_VECTOR r, LIS_REAL *bnrm2)
{
	LIS_Comm comm;  
	LIS_INT	output,conv;
	#ifdef USE_QUAD_PRECISION
	LIS_INT	i;
	#endif
	LIS_MATRIX A;
	LIS_VECTOR x,xx,b,p;
	LIS_REAL nrm2;
	LIS_REAL tol,tol_w,tol_switch;
	
	LIS_DEBUG_FUNC_IN;

	comm = LIS_COMM_WORLD;

	A  = solver->A;
	b  = solver->b;
	x  = solver->x;
	xx = solver->x;
	output     = solver->options[LIS_OPTIONS_OUTPUT];
	conv       = solver->options[LIS_OPTIONS_CONV_COND];
	tol        = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol_w      = solver->params[LIS_PARAMS_RESID_WEIGHT-LIS_OPTIONS_LEN];
	tol_switch = solver->params[LIS_PARAMS_SWITCH_RESID-LIS_OPTIONS_LEN];


	/* initial residual */
	if( M==NULL )
	{
		p = r;
	}
	else
	{
		p = t;
	}

	if( !solver->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
		#ifndef USE_QUAD_PRECISION
			lis_matvec(A,x,p);           /* p = Ax    */
			lis_vector_xpay(b,-1,p);     /* p = b - p */
		#else
			if( solver->precision==LIS_PRECISION_DOUBLE )
			{
				lis_matvec(A,x,p);           /* p = Ax    */
				lis_vector_xpay(b,-1,p);     /* p = b - p */
			}
			else
			{
				lis_matvec(A,xx,p);          /* p = Ax    */
				lis_vector_xpay(b,-1,p);     /* p = b - p */
				
				#ifdef _OPENMP
				#pragma omp parallel for private(i)
				#endif
				for(i=0;i<A->n;i++)
				{
					p->value_lo[i] = 0.0;
				}
				
			}
		#endif
	}
	else
	{
		#ifndef USE_QUAD_PRECISION
			lis_vector_copy(b,p);
		#else
			if( solver->precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_copy(b,p);
			}
			else
			{
				lis_vector_copyex_nm(b,p);
			}
		#endif
	}

	switch(conv)
	{
	case LIS_CONV_COND_NRM2_R:
		lis_vector_nrm2(p,&nrm2);
		*bnrm2 = nrm2;
		solver->tol = tol;
		solver->tol_switch = tol_switch;
		break;
	case LIS_CONV_COND_NRM2_B:
		lis_vector_nrm2(p,&nrm2);
		lis_vector_nrm2(b,bnrm2);
		solver->tol = tol;
		solver->tol_switch = tol_switch;
		break;
	case LIS_CONV_COND_NRM1_B:
		lis_vector_nrm1(p,&nrm2);
		lis_vector_nrm1(b,bnrm2);
		solver->tol = *bnrm2*tol_w + tol;
		solver->tol_switch = *bnrm2*tol_w + tol_switch;
		break;
	}
	if( *bnrm2 == 0.0 )
	{
		*bnrm2 = 1.0;
	}
	else
	{
		*bnrm2 = 1.0 / *bnrm2;
	}
	solver->bnrm = *bnrm2;
	nrm2 = nrm2 * *bnrm2;

	if( output && (r->precision==LIS_PRECISION_QUAD && solver->precision!=LIS_PRECISION_SWITCH) )
	{
	  if( output & LIS_PRINT_MEM ) solver->rhistory[0] = nrm2;
	  if( output & LIS_PRINT_OUT ) lis_printf(comm,"iteration: %5d  relative residual = %e\n", 0, (double)nrm2); 
	}
	if( nrm2 <= fabs(solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN]) )
	{
		solver->retcode = LIS_SUCCESS;
		solver->iter    = 1;
		solver->resid   = nrm2; 
		LIS_DEBUG_FUNC_OUT;
		return LIS_FAILS;
	}

	if( M!=NULL )
	{
		/* r = M^-1 * p */
		lis_psolve(solver, p, r);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_optionC"
LIS_INT lis_solver_set_optionC(LIS_SOLVER solver)
{
	LIS_INT err;    
	LIS_ARGS p;

	LIS_DEBUG_FUNC_IN;

	p = cmd_args->next;
	while( p!=cmd_args )
	{
		err = lis_solver_set_option2(p->arg1,p->arg2,solver);
		if( err )
		  {
		    lis_solver_work_destroy(solver);
		    solver->retcode = err;
		    LIS_DEBUG_FUNC_OUT;		    
		    return err;
		  }
		p = p->next;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option"
LIS_INT lis_solver_set_option(char *text, LIS_SOLVER solver)
{
	LIS_INT err;  
	LIS_ARGS args,p;

	LIS_DEBUG_FUNC_IN;

	lis_text2args(text,&args);
	p = args->next;
	while( p!=args )
	{
		err = lis_solver_set_option2(p->arg1,p->arg2,solver);
		if( err )
		  {
		    lis_solver_work_destroy(solver);
		    solver->retcode = err;
		    LIS_DEBUG_FUNC_OUT;		    
		    return err;
		  }
		p = p->next;
	}
	lis_args_free(args);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option2"
LIS_INT lis_solver_set_option2(char* arg1, char *arg2, LIS_SOLVER solver)
{
	LIS_INT err;
	LIS_INT i;	
	double double_value;

	LIS_DEBUG_FUNC_IN;

	err = 0;

	for(i=0;i<LIS_SOLVER_OPTION_LEN;i++)
	{
		if( strcmp(arg1, LIS_SOLVER_OPTNAME[i])==0 )
		{
			switch( LIS_SOLVER_OPTACT[i] )
			{
			case LIS_OPTIONS_FILE:
				break;
			case LIS_OPTIONS_HELP:
				break;
			case LIS_OPTIONS_VER:
				break;
			case LIS_OPTIONS_SOLVER:
				err = lis_solver_set_option_solver(arg2,solver);
				break;
			case LIS_OPTIONS_PRECON:
				err = lis_solver_set_option_precon(arg2,solver);
				break;
			case LIS_OPTIONS_SCALE:
				err = lis_solver_set_option_scale(arg2,solver);
				break;
			case LIS_OPTIONS_OUTPUT:
				err = lis_solver_set_option_print(arg2,solver);
				break;
			case LIS_OPTIONS_PSOLVER:
				err = lis_solver_set_option_psolver(arg2,solver);
				break;
			case LIS_OPTIONS_PPRECON:
				err = lis_solver_set_option_pprecon(arg2,solver);
				break;
			case LIS_OPTIONS_INITGUESS_ZEROS:
				err = lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_INITGUESS_ZEROS,solver);
				break;
			case LIS_OPTIONS_ADDS:
				err = lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_ADDS,solver);
				break;
			case LIS_OPTIONS_PRECISION:
				err = lis_solver_set_option_precision(arg2,LIS_OPTIONS_PRECISION,solver);
				break;
			case LIS_OPTIONS_USE_AT:
				err = lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_USE_AT,solver);
				break;
			case LIS_OPTIONS_SAAMG_UNSYM:
				err = lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_SAAMG_UNSYM,solver);
				if (solver->options[LIS_OPTIONS_SAAMG_UNSYM])
				  {
				    solver->params[LIS_PARAMS_SAAMG_THETA  -LIS_OPTIONS_LEN] = 0.12;
				  }
				break;
			case LIS_OPTIONS_STORAGE:
				err = lis_solver_set_option_storage(arg2,solver);
				break;
			case LIS_OPTIONS_CONV_COND:
				err = lis_solver_set_option_conv_cond(arg2,solver);
				break;
			default:
				if( LIS_SOLVER_OPTACT[i] < LIS_OPTIONS_LEN )
				{
#ifdef _LONG__LONG
					sscanf(arg2, "%lld", &solver->options[LIS_SOLVER_OPTACT[i]]);
#else
					sscanf(arg2, "%d", &solver->options[LIS_SOLVER_OPTACT[i]]);
#endif
				}
				else
				{
					sscanf(arg2, "%lg", &double_value);
					solver->params[LIS_SOLVER_OPTACT[i]-LIS_OPTIONS_LEN] = double_value;
				}
				break;
			}
		}
		if( err )
		  {
		    lis_solver_work_destroy(solver);
		    solver->retcode = err;
		    LIS_DEBUG_FUNC_OUT;
		    return err;
		  }
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_solver"
LIS_INT lis_solver_set_option_solver(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_SOLVER]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_SOLVER]);
#endif
	}
	else 
	{
		for(i=0;i<LIS_SOLVER_LEN;i++)
		{
			if( strcmp(argv,lis_solver_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_SOLVER] = i+1;
				break;
			}
			else if( i==LIS_SOLVER_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SOLVER is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_psolver"
LIS_INT lis_solver_set_option_psolver(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_PSOLVER]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_PSOLVER]);
#endif
	}
	else
	{
		for(i=0;i<LIS_SOLVER_LEN;i++)
		{
			if( strcmp(argv,lis_solver_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_PSOLVER] = i+1;
				break;
			}
			else if( i==LIS_SOLVER_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PSOLVER is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_precon"
LIS_INT lis_solver_set_option_precon(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_PRECON]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_PRECON]);
#endif
	}
	else
	{
		for(i=0;i<LIS_PRECON_TYPE_LEN;i++)
		{
			if( strcmp(argv,lis_precon_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_PRECON] = i;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
			else if( i==LIS_PRECON_TYPE_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PRECON is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
		for(i=0;i<precon_register_type-LIS_PRECON_TYPE_USERDEF;i++)
		{
			if( strcmp(argv,precon_register_top[i].name)==0 )
			{
				solver->options[LIS_OPTIONS_PRECON] = i+LIS_PRECON_TYPE_USERDEF;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_pprecon"
LIS_INT lis_solver_set_option_pprecon(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_PPRECON]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_PPRECON]);
#endif
	}
	else
	{
		for(i=0;i<LIS_PRECON_TYPE_LEN;i++)
		{
			if( strcmp(argv,lis_precon_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_PPRECON] = i;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
			else if( i==LIS_PRECON_TYPE_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PPRECON is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
		for(i=0;i<precon_register_type-LIS_PRECON_TYPE_USERDEF;i++)
		{
			if( strcmp(argv,precon_register_top[i].name)==0 )
			{
				solver->options[LIS_OPTIONS_PPRECON] = i+LIS_PRECON_TYPE_USERDEF;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_print"
LIS_INT lis_solver_set_option_print(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='3' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_OUTPUT]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_OUTPUT]);
#endif
	}
	else
	{
		for(i=0;i<LIS_PRINT_LEN;i++)
		{
			if( strcmp(argv,lis_print_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_OUTPUT] = i;
				break;
			}
			else if( i==LIS_PRINT_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_OUTPUT is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_scale"
LIS_INT lis_solver_set_option_scale(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='2' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_SCALE]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_SCALE]);
#endif
	}
	else
	{
		for(i=0;i<LIS_SCALE_LEN;i++)
		{
			if( strcmp(argv,lis_scale_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_SCALE] = i;
				break;
			}
			else if( i==LIS_SCALE_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SCALE is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_truefalse"
LIS_INT lis_solver_set_option_truefalse(char *argv, LIS_INT opt, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='1' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[opt]);
#else
		sscanf(argv, "%d", &solver->options[opt]);
#endif
	}
	else
	{
		for(i=0;i<LIS_TRUEFALSE_LEN;i++)
		{
			if( strcmp(argv,lis_truefalse_atoi[i])==0 )
			{
				solver->options[opt] = i;
				break;
			}
			else if( i==LIS_TRUEFALSE_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_TRUEFALSE is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_precision"
LIS_INT lis_solver_set_option_precision(char *argv, LIS_INT opt, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='1' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[opt]);
#else
		sscanf(argv, "%d", &solver->options[opt]);
#endif
	}
	else
	{
		for(i=0;i<LIS_PRECISION_LEN;i++)
		{
			if( strcmp(argv,lis_precision_atoi[i])==0 )
			{
				solver->options[opt] = i;
				break;
			}
			else if( i==LIS_PRECISION_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PRECISION is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_storage"
LIS_INT lis_solver_set_option_storage(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_STORAGE]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_STORAGE]);
#endif
	}
	else
	{
		for(i=0;i<LIS_STORAGE_LEN;i++)
		{
			if( strcmp(argv,lis_storage_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_STORAGE] = i+1;
				break;
			}
			else if( i==LIS_STORAGE_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_STORAGE is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_conv_cond"
LIS_INT lis_solver_set_option_conv_cond(char *argv, LIS_SOLVER solver)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='3' )
	{
#ifdef _LONG__LONG
		sscanf(argv, "%lld", &solver->options[LIS_OPTIONS_CONV_COND]);
#else
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_CONV_COND]);
#endif
	}
	else
	{
		for(i=0;i<LIS_CONV_COND_LEN;i++)
		{
			if( strcmp(argv,lis_conv_cond_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_CONV_COND] = i;
				break;
			}
			else if( i==LIS_CONV_COND_LEN-1 )
			{
				LIS_SETERR(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_CONV_COND is not correct\n");
				LIS_DEBUG_FUNC_OUT;
				return LIS_ERR_ILL_ARG;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_iter"
LIS_INT lis_solver_get_iter(LIS_SOLVER solver, LIS_INT *iter)
{
	LIS_DEBUG_FUNC_IN;

	*iter = solver->iter;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_iterex"
LIS_INT lis_solver_get_iterex(LIS_SOLVER solver, LIS_INT *iter, LIS_INT *iter_double, LIS_INT *iter_quad)
{
	LIS_DEBUG_FUNC_IN;

	*iter = solver->iter;
	*iter_double = solver->iter2;
	*iter_quad = solver->iter - solver->iter2;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_time"
LIS_INT lis_solver_get_time(LIS_SOLVER solver, double *time)
{
	LIS_DEBUG_FUNC_IN;

	*time  = solver->time;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_timeex"
LIS_INT lis_solver_get_timeex(LIS_SOLVER solver, double *time, double *itime, double *ptime, double *p_c_time, double *p_i_time)
{
	LIS_DEBUG_FUNC_IN;

	*time  = solver->time;
	if( itime ) *itime = solver->itime;
	if( ptime ) *ptime = solver->ptime;
	if( p_c_time ) *p_c_time = solver->p_c_time;
	if( p_i_time ) *p_i_time = solver->p_i_time;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_residualnorm"
LIS_INT lis_solver_get_residualnorm(LIS_SOLVER solver, LIS_REAL *residual)
{
	LIS_DEBUG_FUNC_IN;

	*residual  = solver->resid;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_solver"
LIS_INT lis_solver_get_solver(LIS_SOLVER solver, LIS_INT *nsol)
{
	LIS_DEBUG_FUNC_IN;

	*nsol = solver->options[LIS_OPTIONS_SOLVER];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_precon"
LIS_INT lis_solver_get_precon(LIS_SOLVER solver, LIS_INT *precon_type)
{
	LIS_DEBUG_FUNC_IN;

	*precon_type = solver->options[LIS_OPTIONS_PRECON];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_status"
LIS_INT lis_solver_get_status(LIS_SOLVER solver, LIS_INT *status)
{
	LIS_DEBUG_FUNC_IN;

	*status = solver->retcode;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_rhistory"
LIS_INT lis_solver_get_rhistory(LIS_SOLVER solver, LIS_VECTOR v)
{
	LIS_INT	i,n,maxiter;

	maxiter = solver->iter+1;
        if( solver->retcode!=LIS_SUCCESS )
	  {
	    maxiter--;
	  }
	n = _min(v->n,maxiter);
	for(i=0;i<n;i++)
	{
		v->value[i] = solver->rhistory[i];
	}
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_solvername"
LIS_INT lis_solver_get_solvername(LIS_INT solver, char *solvername)
{
	LIS_DEBUG_FUNC_IN;

	if( solver < 1 || solver > LIS_SOLVERS_LEN )
	{
		solvername = NULL;
		return LIS_FAILS;
	}
	strcpy(solvername,lis_solvername[solver]);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_preconname"
LIS_INT lis_solver_get_preconname(LIS_INT precon_type, char *preconname)
{
	LIS_DEBUG_FUNC_IN;

	if( precon_type < 0 || precon_type > LIS_PRECON_TYPE_LEN-1 )
	{
		preconname = NULL;
		return LIS_FAILS;
	}
	strcpy(preconname,lis_preconname[precon_type]);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_residual_nrm2_r"
LIS_INT lis_solver_get_residual_nrm2_r(LIS_VECTOR r, LIS_SOLVER solver, LIS_REAL *res)
{
	LIS_DEBUG_FUNC_IN;

	lis_vector_nrm2(r,res);
	*res = *res * solver->bnrm;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_residual_nrm1_b"
LIS_INT lis_solver_get_residual_nrm1_b(LIS_VECTOR r, LIS_SOLVER solver, LIS_REAL *res)
{
	LIS_DEBUG_FUNC_IN;

	lis_vector_nrm1(r,res);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_shadowresidual"
LIS_INT lis_solver_set_shadowresidual(LIS_SOLVER solver, LIS_VECTOR r0, LIS_VECTOR rs0)
{
        unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	LIS_INT i,n,resid;

	LIS_DEBUG_FUNC_IN;

	resid = solver->options[LIS_OPTIONS_INIT_SHADOW_RESID];
	if( resid==LIS_RANDOM )
	{
		n     = solver->A->n;
		init_by_array(init, length);

		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		#endif
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				rs0->value[i] = genrand_real1();
			}
		}
		#ifdef USE_QUAD_PRECISION
		else
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				rs0->value[i]    = genrand_real1();
				rs0->value_lo[i] = 0.0;
			}
		}
		#endif
	}
	else
	{
		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		#endif
		{
			lis_vector_copy(r0,rs0);
			lis_vector_conjugate(rs0);
		}
		#ifdef USE_QUAD_PRECISION
		else
		{
			lis_vector_copyex_mm(r0,rs0);
		}
		#endif
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
