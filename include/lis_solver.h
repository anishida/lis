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


#ifndef __LIS_SOLVER_H__
#define __LIS_SOLVER_H__

typedef LIS_INT (*LIS_SOLVER_CHECK_PARAMS)(LIS_SOLVER solver);
typedef LIS_INT (*LIS_SOLVER_MALLOC_WORK)(LIS_SOLVER solver);
typedef LIS_INT (*LIS_SOLVER_EXECUTE)(LIS_SOLVER solver);
typedef LIS_INT (*LIS_SOLVER_GET_RESIDUAL)(LIS_VECTOR r, LIS_SOLVER solver, LIS_REAL *res);

#ifdef __cplusplus
extern "C"
{
#endif
	extern LIS_INT lis_solver_init(LIS_SOLVER solver);
	extern LIS_INT lis_solver_work_destroy(LIS_SOLVER solver);
	extern LIS_INT lis_solver_get_initial_residual(LIS_SOLVER solver, LIS_PRECON M, LIS_VECTOR t, LIS_VECTOR r, LIS_REAL *bnrm2);
	extern LIS_INT lis_solver_set_option2(char* arg1, char *arg2, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_solver(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_psolver(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_precon(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_pprecon(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_scale(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_print(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_truefalse(char *argv, LIS_INT opt, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_precision(char *argv, LIS_INT opt, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_storage(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_option_conv_cond(char *argv, LIS_SOLVER solver);
	extern LIS_INT lis_solver_get_residual_nrm2_r(LIS_VECTOR r, LIS_SOLVER solver, LIS_REAL *res);
	extern LIS_INT lis_solver_get_residual_nrm1_b(LIS_VECTOR r, LIS_SOLVER solver, LIS_REAL *res);
	extern LIS_INT lis_solver_set_shadowresidual(LIS_SOLVER solver, LIS_VECTOR r0, LIS_VECTOR rs0);
	extern LIS_SOLVER_GET_RESIDUAL lis_solver_get_residual[3];

/*******************/
/* CG              */
/*******************/

	extern LIS_INT lis_cg(LIS_SOLVER solver);
	extern LIS_INT lis_cg_quad(LIS_SOLVER solver);
	extern LIS_INT lis_cg_switch(LIS_SOLVER solver);
	extern LIS_INT lis_cg_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_cg_malloc_work(LIS_SOLVER solver);

/*******************/
/* BiCG            */
/*******************/

	extern LIS_INT lis_bicg(LIS_SOLVER solver);
	extern LIS_INT lis_bicg_quad(LIS_SOLVER solver);
	extern LIS_INT lis_bicg_switch(LIS_SOLVER solver);
	extern LIS_INT lis_bicg_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_bicg_malloc_work(LIS_SOLVER solver);

/*******************/
/* CGS             */
/*******************/

	extern LIS_INT lis_cgs(LIS_SOLVER solver);
	extern LIS_INT lis_cgs_quad(LIS_SOLVER solver);
	extern LIS_INT lis_cgs_switch(LIS_SOLVER solver);
	extern LIS_INT lis_cgs_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_cgs_malloc_work(LIS_SOLVER solver);

/*******************/
/* BiCGSTAB        */
/*******************/

	extern LIS_INT lis_bicgstab(LIS_SOLVER solver);
	extern LIS_INT lis_bicgstab_quad(LIS_SOLVER solver);
	extern LIS_INT lis_bicgstab_switch(LIS_SOLVER solver);
	extern LIS_INT lis_bicgstab_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_bicgstab_malloc_work(LIS_SOLVER solver);

/*******************/
/* BiCGSTAB(l)     */
/*******************/

	extern LIS_INT lis_bicgstabl(LIS_SOLVER solver);
	extern LIS_INT lis_bicgstabl_quad(LIS_SOLVER solver);
	extern LIS_INT lis_bicgstabl_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_bicgstabl_malloc_work(LIS_SOLVER solver);

/*******************/
/* GPBiCG          */
/*******************/

	extern LIS_INT lis_gpbicg(LIS_SOLVER solver);
	extern LIS_INT lis_gpbicg_quad(LIS_SOLVER solver);
	extern LIS_INT lis_gpbicg_switch(LIS_SOLVER solver);
	extern LIS_INT lis_gpbicg_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_gpbicg_malloc_work(LIS_SOLVER solver);

/*******************/
/* BiCGSafe        */
/*******************/

	extern LIS_INT lis_bicgsafe(LIS_SOLVER solver);
	extern LIS_INT lis_bicgsafe_quad(LIS_SOLVER solver);
	extern LIS_INT lis_bicgsafe_switch(LIS_SOLVER solver);
	extern LIS_INT lis_bicgsafe_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_bicgsafe_malloc_work(LIS_SOLVER solver);

/*******************/
/* TFQMR           */
/*******************/

	extern LIS_INT lis_tfqmr(LIS_SOLVER solver);
	extern LIS_INT lis_tfqmr_quad(LIS_SOLVER solver);
	extern LIS_INT lis_tfqmr_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_tfqmr_malloc_work(LIS_SOLVER solver);

/*******************/
/* Orthomin(m)     */
/*******************/

	extern LIS_INT lis_orthomin(LIS_SOLVER solver);
	extern LIS_INT lis_orthomin_quad(LIS_SOLVER solver);
	extern LIS_INT lis_orthomin_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_orthomin_malloc_work(LIS_SOLVER solver);

/*******************/
/* GMRES(m)        */
/*******************/

	extern LIS_INT lis_gmres(LIS_SOLVER solver);
	extern LIS_INT lis_gmres_quad(LIS_SOLVER solver);
	extern LIS_INT lis_gmres_switch(LIS_SOLVER solver);
	extern LIS_INT lis_gmres_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_gmres_malloc_work(LIS_SOLVER solver);

/*******************/
/* Jacobi          */
/*******************/

	extern LIS_INT lis_jacobi(LIS_SOLVER solver);
	extern LIS_INT lis_jacobi_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_jacobi_malloc_work(LIS_SOLVER solver);

/*******************/
/* Gauss-Seidel    */
/*******************/

	extern LIS_INT lis_gs(LIS_SOLVER solver);
	extern LIS_INT lis_gs_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_gs_malloc_work(LIS_SOLVER solver);

/*******************/
/* SOR             */
/*******************/

	extern LIS_INT lis_sor(LIS_SOLVER solver);
	extern LIS_INT lis_sor_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_sor_malloc_work(LIS_SOLVER solver);

/*******************/
/* CR              */
/*******************/

	extern LIS_INT lis_cr(LIS_SOLVER solver);
	extern LIS_INT lis_cr_quad(LIS_SOLVER solver);
/*	extern LIS_INT lis_cr_switch(LIS_SOLVER solver);*/
	extern LIS_INT lis_cr_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_cr_malloc_work(LIS_SOLVER solver);

/*******************/
/* BiCR            */
/*******************/

	extern LIS_INT lis_bicr(LIS_SOLVER solver);
	extern LIS_INT lis_bicr_quad(LIS_SOLVER solver);
/*	extern LIS_INT lis_bicr_switch(LIS_SOLVER solver);*/
	extern LIS_INT lis_bicr_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_bicr_malloc_work(LIS_SOLVER solver);

/*******************/
/* CRS             */
/*******************/

	extern LIS_INT lis_crs(LIS_SOLVER solver);
	extern LIS_INT lis_crs_quad(LIS_SOLVER solver);
/*	extern LIS_INT lis_crs_switch(LIS_SOLVER solver);*/
	extern LIS_INT lis_crs_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_crs_malloc_work(LIS_SOLVER solver);

/*******************/
/* BiCRSTAB        */
/*******************/

	extern LIS_INT lis_bicrstab(LIS_SOLVER solver);
	extern LIS_INT lis_bicrstab_quad(LIS_SOLVER solver);
/*	extern LIS_INT lis_bicrstab_switch(LIS_SOLVER solver);*/
	extern LIS_INT lis_bicrstab_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_bicrstab_malloc_work(LIS_SOLVER solver);

/*******************/
/* GPBiCR          */
/*******************/

	extern LIS_INT lis_gpbicr(LIS_SOLVER solver);
	extern LIS_INT lis_gpbicr_quad(LIS_SOLVER solver);
/*	extern LIS_INT lis_gpbicr_switch(LIS_SOLVER solver);*/
	extern LIS_INT lis_gpbicr_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_gpbicr_malloc_work(LIS_SOLVER solver);

/*******************/
/* BiCRSafe        */
/*******************/

	extern LIS_INT lis_bicrsafe(LIS_SOLVER solver);
	extern LIS_INT lis_bicrsafe_quad(LIS_SOLVER solver);
/*	extern LIS_INT lis_bicrsafe_switch(LIS_SOLVER solver);*/
	extern LIS_INT lis_bicrsafe_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_bicrsafe_malloc_work(LIS_SOLVER solver);

/*******************/
/* FGMRES(m)       */
/*******************/

	extern LIS_INT lis_fgmres(LIS_SOLVER solver);
	extern LIS_INT lis_fgmres_quad(LIS_SOLVER solver);
/*	extern LIS_INT lis_fgmres_switch(LIS_SOLVER solver);*/
	extern LIS_INT lis_fgmres_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_fgmres_malloc_work(LIS_SOLVER solver);

/*******************/
/* IDR(s)          */
/*******************/

	extern LIS_INT lis_idrs(LIS_SOLVER solver);
	extern LIS_INT lis_idrs_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_idrs_malloc_work(LIS_SOLVER solver);

/*******************/
/* IDR(1)          */
/*******************/

	extern LIS_INT lis_idr1(LIS_SOLVER solver);
	extern LIS_INT lis_idr1_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_idr1_malloc_work(LIS_SOLVER solver);

/*******************/
/* GMRES(m)        */
/*******************/

	extern LIS_INT lis_minres(LIS_SOLVER solver);
	extern LIS_INT lis_minres_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_minres_malloc_work(LIS_SOLVER solver);

/*******************/
/* COCG            */
/*******************/

	extern LIS_INT lis_cocg(LIS_SOLVER solver);
	extern LIS_INT lis_cocg_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_cocg_malloc_work(LIS_SOLVER solver);

/*******************/
/* COCR            */
/*******************/

	extern LIS_INT lis_cocr(LIS_SOLVER solver);
	extern LIS_INT lis_cocr_check_params(LIS_SOLVER solver);
	extern LIS_INT lis_cocr_malloc_work(LIS_SOLVER solver);

  
#ifdef __cplusplus
}
#endif
#endif
