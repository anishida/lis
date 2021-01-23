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


#ifndef __LIS_FORTRAN_H__
#define __LIS_FORTRAN_H__

#ifdef _WIN64
typedef long long LIS_MATRIX_F;
typedef long long LIS_VECTOR_F;
typedef long long LIS_SOLVER_F;
typedef long long LIS_PRECON_F;
typedef long long LIS_ESOLVER_F;
typedef long long LIS_SCALAR_F;
#else
typedef long LIS_MATRIX_F;
typedef long LIS_VECTOR_F;
typedef long LIS_SOLVER_F;
typedef long LIS_PRECON_F;
typedef long LIS_ESOLVER_F;
typedef long LIS_SCALAR_F;
#endif

#ifdef USE_MPI
	#include <mpi.h>
	typedef MPI_Fint LIS_Comm_f;
#else
	typedef LIS_INT	LIS_Comm_f;
#endif

#ifdef _WIN64
#define LIS_V2P(a) (*(long long *)(a))
#define LIS_P2V(a) ((long long)(a))
#else
#define LIS_V2P(a) (*(long *)(a))
#define LIS_P2V(a) ((long)(a))
#endif

extern LIS_Comm_f lis_comm_world_f;

/**************/
/* ARRAY      */
/**************/

#define lis_array_swap_f F77_FUNC_(lis_array_swap_f, LIS_ARRAY_SWAP_F) 
#define lis_array_copy_f F77_FUNC_(lis_array_copy_f, LIS_ARRAY_COPY_F) 
#define lis_array_axpy_f F77_FUNC_(lis_array_axpy_f, LIS_ARRAY_AXPY_F) 
#define lis_array_xpay_f F77_FUNC_(lis_array_xpay_f, LIS_ARRAY_XPAY_F) 
#define lis_array_axpyz_f F77_FUNC_(lis_array_axpyz_f, LIS_ARRAY_AXPYZ_F) 
#define lis_array_scale_f F77_FUNC_(lis_array_scale_f, LIS_ARRAY_SCALE_F) 
#define lis_array_pmul_f F77_FUNC_(lis_array_pmul_f, LIS_ARRAY_PMUL_F) 
#define lis_array_pdiv_f F77_FUNC_(lis_array_pdiv_f, LIS_ARRAY_PDIV_F) 
#define lis_array_set_all_f F77_FUNC_(lis_array_set_all_f, LIS_ARRAY_SET_ALL_F) 
#define lis_array_abs_f F77_FUNC_(lis_array_abs_f, LIS_ARRAY_ABS_F) 
#define lis_array_reciprocal_f F77_FUNC_(lis_array_reciprocal_f, LIS_ARRAY_RECIPROCAL_F)
#define lis_array_conjugate_f F77_FUNC_(lis_array_conjugate_f, LIS_ARRAY_CONJUGATE_F)
#define lis_array_shift_f F77_FUNC_(lis_array_shift_f, LIS_ARRAY_SHIFT_F) 
#define lis_array_dot_f F77_FUNC_(lis_array_dot_f, LIS_ARRAY_DOT_F)
#define lis_array_nhdot_f F77_FUNC_(lis_array_nhdot_f, LIS_ARRAY_NHDOT_F) 
#define lis_array_dot2_f F77_FUNC_(lis_array_dot2_f, LIS_ARRAY_DOT2_F) 
#define lis_array_nrm2_f F77_FUNC_(lis_array_nrm2_f, LIS_ARRAY_NRM2_F) 
#define lis_array_nrm1_f F77_FUNC_(lis_array_nrm1_f, LIS_ARRAY_NRM1_F) 
#define lis_array_nrmi_f F77_FUNC_(lis_array_nrmi_f, LIS_ARRAY_NRMI_F) 
#define lis_array_sum_f F77_FUNC_(lis_array_sum_f, LIS_ARRAY_SUM_F) 
#define lis_array_matvec_f F77_FUNC_(lis_array_matvec_f, LIS_ARRAY_MATVEC_F) 
#define lis_array_matvech_f F77_FUNC_(lis_array_matvech_f, LIS_ARRAY_MATVECH_F)   
#define lis_array_matvec2_f F77_FUNC_(lis_array_matvec2_f, LIS_ARRAY_matvec2_F) 
#define lis_array_matinv_f F77_FUNC_(lis_array_matinv_f, LIS_ARRAY_MATINV_F) 
#define lis_array_LUdecomp_f F77_FUNC_(lis_array_LUdecomp_f, LIS_ARRAY_LUDECOMP_F) 
#define lis_array_invGauss_f F77_FUNC_(lis_array_invGauss_f, LIS_ARRAY_INVGAUSS_F) 
#define lis_array_invvec_f F77_FUNC_(lis_array_invvec_f, LIS_ARRAY_INVVEC_F) 
#define lis_array_invvect_f F77_FUNC_(lis_array_invvect_f, LIS_ARRAY_INVVECT_F) 
#define lis_array_matmat_f F77_FUNC_(lis_array_matmat_f, LIS_ARRAY_MATMAT_F) 
#define lis_array_matmat2_f F77_FUNC_(lis_array_matmat2_f, LIS_ARRAY_MATMAT2_F) 
#define lis_array_solve_f F77_FUNC_(lis_array_solve_f, LIS_ARRAY_SOLVE_F) 
#define lis_array_cgs_f F77_FUNC_(lis_array_cgs_f, LIS_ARRAY_CGS_F) 
#define lis_array_mgs_f F77_FUNC_(lis_array_mgs_f, LIS_ARRAY_MGS_F) 
#define lis_array_qr_f F77_FUNC_(lis_array_qr_f, LIS_ARRAY_QR_F) 
#define lis_array_power_f F77_FUNC_(lis_array_power_f, LIS_ARRAY_POWER_F) 

/**************/
/* ESOLVER     */
/**************/

#define lis_esolver_create_f F77_FUNC_(lis_esolver_create_f, LIS_ESOLVER_CREATE_F) 
#define lis_esolver_destroy_f F77_FUNC_(lis_esolver_destroy_f, LIS_ESOLVER_DESTROY_F) 
#define lis_esolve_f F77_FUNC_(lis_esolve_f, LIS_ESOLVE_F)
#define lis_gesolve_f F77_FUNC_(lis_gesolve_f, LIS_GESOLVE_F)   
#define lis_esolver_set_option_f F77_FUNC_(lis_esolver_set_option_f, LIS_ESOLVER_SET_OPTION_F) 
#define lis_esolver_get_iter_f F77_FUNC_(lis_esolver_get_iter_f, LIS_ESOLVER_GET_ITER_F) 
#define lis_esolver_get_iterex_f F77_FUNC_(lis_esolver_get_iterex_f, LIS_ESOLVER_GET_ITEREX_F) 
#define lis_esolver_get_time_f F77_FUNC_(lis_esolver_get_time_f, LIS_ESOLVER_GET_TIME_F) 
#define lis_esolver_get_timeex_f F77_FUNC_(lis_esolver_get_timeex_f, LIS_ESOLVER_GET_TIMEEX_F) 
#define lis_esolver_get_residualnorm_f F77_FUNC_(lis_esolver_get_residualnorm_f, LIS_ESOLVER_GET_RESIDUALNORM_F) 
#define lis_esolver_get_rhistory_f F77_FUNC_(lis_esolver_get_rhistory_f, LIS_ESOLVER_GET_RHISTORY_F) 
#define lis_esolver_get_evalues_f F77_FUNC_(lis_esolver_get_evalues_f, LIS_ESOLVER_GET_EVALUES_F)
#define lis_esolver_get_evectors_f F77_FUNC_(lis_esolver_get_evectors_f, LIS_ESOLVER_GET_EVECTORS_F) 
#define lis_esolver_get_residualnorms_f F77_FUNC_(lis_esolver_get_residualnorms_f, LIS_ESOLVER_GET_RESIDUALNORMS_F) 
#define lis_esolver_get_iters_f F77_FUNC_(lis_esolver_get_iters_f, LIS_ESOLVER_GET_ITERS_F) 
#define lis_esolver_get_specific_evalue_f F77_FUNC_(lis_esolver_get_specific_evalue_f, LIS_ESOLVER_GET_SPECIFIC_EVALUE_F)
#define lis_esolver_get_specific_evector_f F77_FUNC_(lis_esolver_get_specific_evector_f, LIS_ESOLVER_GET_SPECIFIC_EVECTOR_F) 
#define lis_esolver_get_specific_residualnorm_f F77_FUNC_(lis_esolver_get_specific_residualnorm_f, LIS_ESOLVER_GET_SPECIFIC_RESIDUALNORM_F) 
#define lis_esolver_get_specific_iter_f F77_FUNC_(lis_esolver_get_specific_iter_f, LIS_ESOLVER_GET_SPECIFIC_ITER_F) 
#define lis_esolver_get_esolver_f F77_FUNC_(lis_esolver_get_esolver_f, LIS_ESOLVER_GET_ESOLVER_F) 
#define lis_esolver_get_esolvername_f F77_FUNC_(lis_esolver_get_esolvername_f, LIS_ESOLVER_GET_ESOLVERNAME_F) 
#define lis_esolver_set_optionC_f F77_FUNC_(lis_esolver_set_optionc_f, LIS_ESOLVER_SET_OPTIONC_F) 

/**************/
/* MATRIX     */
/**************/

#define lis_matrix_create_f F77_FUNC_(lis_matrix_create_f, LIS_MATRIX_CREATE_F) 
#define lis_matrix_duplicate_f F77_FUNC_(lis_matrix_duplicate_f, LIS_MATRIX_DUPLICATE_F) 
#define lis_matrix_destroy_f F77_FUNC_(lis_matrix_destroy_f, LIS_MATRIX_DESTROY_F) 
#define lis_matrix_set_size_f F77_FUNC_(lis_matrix_set_size_f, LIS_MATRIX_SET_SIZE_F) 
#define lis_matrix_get_size_f F77_FUNC_(lis_matrix_get_size_f, LIS_MATRIX_GET_SIZE_F) 
#define lis_matrix_get_range_f F77_FUNC_(lis_matrix_get_range_f, LIS_MATRIX_GET_RANGE_F) 
#define lis_matrix_get_nnz_f F77_FUNC_(lis_matrix_get_nnz_f, LIS_MATRIX_GET_NNZ_F) 
#define lis_matrix_set_value_f F77_FUNC_(lis_matrix_set_value_f, LIS_MATRIX_SET_VALUE_F) 
/*NEH support for extended "solve_kernel" workflow*/
#define lis_matrix_psd_set_value_f F77_FUNC_(lis_matrix_psd_set_value_f, LIS_MATRIX_PSD_SET_VALUE_F) 
#define lis_matrix_get_type_f F77_FUNC_(lis_matrix_get_type_f, LIS_MATRIX_GET_TYPE_F) 
#define lis_matrix_set_type_f F77_FUNC_(lis_matrix_set_type_f, LIS_MATRIX_SET_TYPE_F)
#define lis_matrix_set_csr_f F77_FUNC_(lis_matrix_set_csr_f, LIS_MATRIX_SET_CSR_F)
#define lis_matrix_set_csc_f F77_FUNC_(lis_matrix_set_csc_f, LIS_MATRIX_SET_CSC_F)
#define lis_matrix_set_msr_f F77_FUNC_(lis_matrix_set_msr_f, LIS_MATRIX_SET_MSR_F)
#define lis_matrix_set_dia_f F77_FUNC_(lis_matrix_set_dia_f, LIS_MATRIX_SET_DIA_F)
#define lis_matrix_set_ell_f F77_FUNC_(lis_matrix_set_ell_f, LIS_MATRIX_SET_ELL_F)
#define lis_matrix_set_jad_f F77_FUNC_(lis_matrix_set_jad_f, LIS_MATRIX_SET_JAD_F)
#define lis_matrix_set_bsr_f F77_FUNC_(lis_matrix_set_bsr_f, LIS_MATRIX_SET_BSR_F)
#define lis_matrix_set_bsc_f F77_FUNC_(lis_matrix_set_bsc_f, LIS_MATRIX_SET_BSC_F)
#define lis_matrix_set_coo_f F77_FUNC_(lis_matrix_set_coo_f, LIS_MATRIX_SET_COO_F)
#define lis_matrix_set_dns_f F77_FUNC_(lis_matrix_set_dns_f, LIS_MATRIX_SET_DNS_F)
#define lis_matrix_set_vbr_f F77_FUNC_(lis_matrix_set_vbr_f, LIS_MATRIX_SET_VBR_F)
#define lis_matrix_assemble_f F77_FUNC_(lis_matrix_assemble_f, LIS_MATRIX_ASSEMBLE_F) 
#define lis_matrix_is_assembled_f F77_FUNC_(lis_matrix_is_assembled_f, LIS_MATRIX_IS_ASSEMBLED_F) 
#define lis_matrix_malloc_f F77_FUNC_(lis_matrix_malloc_f, LIS_MATRIX_MALLOC_F)
#define lis_matrix_malloc_csr_f F77_FUNC_(lis_matrix_malloc_csr_f, LIS_MATRIX_MALLOC_CSR_F)
#define lis_matrix_malloc_csc_f F77_FUNC_(lis_matrix_malloc_csc_f, LIS_MATRIX_MALLOC_CSC_F)
#define lis_matrix_malloc_bsr_f F77_FUNC_(lis_matrix_malloc_bsr_f, LIS_MATRIX_MALLOC_BSR_F)
#define lis_matrix_malloc_msr_f F77_FUNC_(lis_matrix_malloc_msr_f, LIS_MATRIX_MALLOC_MSR_F)
#define lis_matrix_malloc_ell_f F77_FUNC_(lis_matrix_malloc_ell_f, LIS_MATRIX_MALLOC_ELL_F)
#define lis_matrix_malloc_jad_f F77_FUNC_(lis_matrix_malloc_jad_f, LIS_MATRIX_MALLOC_JAD_F)
#define lis_matrix_malloc_dia_f F77_FUNC_(lis_matrix_malloc_dia_f, LIS_MATRIX_MALLOC_DIA_F)
#define lis_matrix_malloc_bsc_f F77_FUNC_(lis_matrix_malloc_bsc_f, LIS_MATRIX_MALLOC_BSC_F)
#define lis_matrix_malloc_vbr_f F77_FUNC_(lis_matrix_malloc_vbr_f, LIS_MATRIX_MALLOC_VBR_F)
#define lis_matrix_malloc_coo_f F77_FUNC_(lis_matrix_malloc_coo_f, LIS_MATRIX_MALLOC_COO_F)
#define lis_matrix_malloc_dns_f F77_FUNC_(lis_matrix_malloc_dns_f, LIS_MATRIX_MALLOC_DNS_F)
#define lis_matrix_convert_f F77_FUNC_(lis_matrix_convert_f, LIS_MATRIX_CONVERT_F)
#define lis_matrix_copy_f F77_FUNC_(lis_matrix_copy_f, LIS_MATRIX_COPY_F)
#define lis_matrix_axpy_f F77_FUNC_(lis_matrix_axpy_f, LIS_MATRIX_AXPY_F)
#define lis_matrix_xpay_f F77_FUNC_(lis_matrix_xpay_f, LIS_MATRIX_XPAY_F)
#define lis_matrix_axpyz_f F77_FUNC_(lis_matrix_axpyz_f, LIS_MATRIX_AXPYZ_F)  
#define lis_matrix_scale_f F77_FUNC_(lis_matrix_scale_f, LIS_MATRIX_SCALE_F)
#define lis_matrix_psd_reset_scale_f F77_FUNC_(lis_matrix_psd_reset_scale_f, LIS_MATRIX_PSD_RESET_SCALE_F)
#define lis_matrix_get_diagonal_f F77_FUNC_(lis_matrix_get_diagonal_f, LIS_MATRIX_GET_DIAGONAL_F)
#define lis_matrix_shift_diagonal_f F77_FUNC_(lis_matrix_shift_diagonal_f, LIS_MATRIX_SHIFT_DIAGONAL_F)
#define lis_matrix_shift_matrix_f F77_FUNC_(lis_matrix_shift_matrix_f, LIS_MATRIX_SHIFT_GENERAL_F)
#define lis_matrix_set_blocksize_f F77_FUNC_(lis_matrix_set_blocksize_f, LIS_MATRIX_SET_BLOCKSIZE_F)
#define lis_matrix_unset_f F77_FUNC_(lis_matrix_unset_f, LIS_MATRIX_UNSET_F)

/**************/
/* MATVEC     */
/**************/

#define lis_matvec_f F77_FUNC_(lis_matvec_f, LIS_MATVEC_F) 
#define lis_matvech_f F77_FUNC_(lis_matvech_f, LIS_MATVECH_F)   

/**************/
/* PRECON     */
/**************/

#define lis_precon_create_f F77_FUNC_(lis_precon_create_f, LIS_PRECON_CREATE_F) 
/*NEH support for extended "solve_kernel" workflow*/
#define lis_precon_psd_create_f F77_FUNC_(lis_precon_psd_create_f, LIS_PRECON_PSD_CREATE_F) 
/*NEH support for extended "solve_kernel" workflow*/
#define lis_precon_psd_update_f F77_FUNC_(lis_precon_psd_update_f, LIS_PRECON_PSD_UPDATE_F) 
#define lis_precon_destroy_f F77_FUNC_(lis_precon_destroy_f, LIS_PRECON_DESTROY_F) 

/**************/
/* SOLVER     */
/**************/

#define lis_solver_create_f F77_FUNC_(lis_solver_create_f, LIS_SOLVER_CREATE_F) 
#define lis_solver_destroy_f F77_FUNC_(lis_solver_destroy_f, LIS_SOLVER_DESTROY_F) 
/*NEH support for extended "solve_kernel" workflow*/
#define lis_solver_set_matrix_f F77_FUNC_(lis_solver_set_matrix_f, LIS_SOLVER_SET_MATRIX_F) 
#define lis_solve_f F77_FUNC_(lis_solve_f, LIS_SOLVE_F)
#define lis_solve_setup_f F77_FUNC_(lis_solve_setup_f, LIS_SOLVE_SETUP_F) 
#define lis_solve_kernel_f F77_FUNC_(lis_solve_kernel_f, LIS_SOLVE_KERNEL_F)
#define lis_solver_get_status_f F77_FUNC_(lis_solver_get_status_f, LIS_SOLVER_GET_STATUS_F) 
#define lis_solver_set_option_f F77_FUNC_(lis_solver_set_option_f, LIS_SOLVER_SET_OPTION_F) 
#define lis_solver_get_iter_f F77_FUNC_(lis_solver_get_iter_f, LIS_SOLVER_GET_ITER_F) 
#define lis_solver_get_iterex_f F77_FUNC_(lis_solver_get_iterex_f, LIS_SOLVER_GET_ITEREX_F) 
#define lis_solver_get_time_f F77_FUNC_(lis_solver_get_time_f, LIS_SOLVER_GET_TIME_F) 
#define lis_solver_get_timeex_f F77_FUNC_(lis_solver_get_timeex_f, LIS_SOLVER_GET_TIMEEX_F) 
#define lis_solver_get_residualnorm_f F77_FUNC_(lis_solver_get_residualnorm_f, LIS_SOLVER_GET_RESIDUALNORM_F) 
#define lis_solver_get_rhistory_f F77_FUNC_(lis_solver_get_rhistory_f, LIS_SOLVER_GET_RHISTORY_F) 
#define lis_solver_get_solver_f F77_FUNC_(lis_solver_get_solver_f, LIS_SOLVER_GET_SOLVER_F) 
#define lis_solver_get_precon_f F77_FUNC_(lis_solver_get_precon_f, LIS_SOLVER_GET_PRECON_F) 
#define lis_solver_get_solvername_f F77_FUNC_(lis_solver_get_solvername_f, LIS_SOLVER_GET_SOLVERNAME_F) 
#define lis_solver_get_preconname_f F77_FUNC_(lis_solver_get_preconname_f, LIS_SOLVER_GET_PRECONNAME_F) 
#define lis_solver_set_optionC_f F77_FUNC_(lis_solver_set_optionc_f, LIS_SOLVER_SET_OPTIONC_F) 

/**************/
/* SYSTEM     */
/**************/

#define lis_set_comm_world_f F77_FUNC_(lis_set_comm_world_f, LIS_SET_COMM_WORLD_F)
#define lis_finitialize_f F77_FUNC_(lis_finitialize_f, LIS_FINITIALIZE_F) 
#define lis_finalize_f F77_FUNC_(lis_finalize_f, LIS_FINALIZE_F) 
#define lis_wtime_f F77_FUNC_(lis_wtime_f , LIS_WTIME_F)
#define CHKERR_f F77_FUNC_(chkerr_f, CHKERR_F) 
#define lis_set_argv_begin_f F77_FUNC_(lis_set_argv_begin_f, LIS_SET_ARGV_BEGIN_F) 
#define lis_set_argv_f F77_FUNC_(lis_set_argv_f, LIS_SET_ARGV_F) 
#define lis_set_argv_end_f F77_FUNC_(lis_set_argv_end_f, LIS_SET_ARGV_END_F) 
#define lis_arg2args_f F77_FUNC_(lis_arg2args_f, LIS_ARG2ARGS_F) 
#define lis_input_f F77_FUNC_(lis_input_f, LIS_INPUT_F)
#define lis_input_matrix_f F77_FUNC_(lis_input_matrix_f, LIS_INPUT_MATRIX_F)
#define lis_input_vector_f F77_FUNC_(lis_input_vector_f, LIS_INPUT_VECTOR_F)
#define lis_output_f F77_FUNC_(lis_output_f, LIS_OUTPUT_F) 
#define lis_output_matrix_f F77_FUNC_(lis_output_matrix_f, LIS_OUTPUT_MATRIX_F) 
#define lis_output_vector_f F77_FUNC_(lis_output_vector_f, LIS_OUTPUT_VECTOR_F) 
#define lis_solver_output_rhistory_f F77_FUNC_(lis_solver_output_rhistory_f, LIS_SOLVER_OUTPUT_RHISTORY_F) 
#define lis_esolver_output_rhistory_f F77_FUNC_(lis_esolver_output_rhistory_f, LIS_ESOLVER_OUTPUT_RHISTORY_F) 

/**************/
/* VECTOR     */
/**************/

#define lis_vector_create_f F77_FUNC_(lis_vector_create_f, LIS_VECTOR_CREATE_F) 
#define lis_vector_duplicate_f F77_FUNC_(lis_vector_duplicate_f, LIS_VECTOR_DUPLICATE_F) 
#define lis_vector_destroy_f F77_FUNC_(lis_vector_destroy_f, LIS_VECTOR_DESTROY_F) 
#define lis_vector_set_size_f F77_FUNC_(lis_vector_set_size_f, LIS_VECTOR_SET_SIZE_F) 
/*NEH support for extended "solve_kernel" workflow*/
#define lis_vector_psd_reset_scale_f F77_FUNC_(lis_vector_psd_reset_scale_f, LIS_VECTOR_PSD_RESET_SCALE_F) 
#define lis_vector_get_size_f F77_FUNC_(lis_vector_get_size_f, LIS_VECTOR_GET_SIZE_F) 
#define lis_vector_get_range_f F77_FUNC_(lis_vector_get_range_f, LIS_VECTOR_GET_RANGE_F) 
#define lis_vector_set_value_f F77_FUNC_(lis_vector_set_value_f, LIS_VECTOR_SET_VALUE_F) 
#define lis_vector_set_values_f F77_FUNC_(lis_vector_set_values_f, LIS_VECTOR_SET_VALUES_F) 
#define lis_vector_set_values2_f F77_FUNC_(lis_vector_set_values2_f, LIS_VECTOR_SET_VALUES2_F) 
#define lis_vector_get_value_f F77_FUNC_(lis_vector_get_value_f, LIS_VECTOR_GET_VALUE_F) 
#define lis_vector_get_values_f F77_FUNC_(lis_vector_get_values_f, LIS_VECTOR_GET_VALUES_F) 
#define lis_vector_scatter_f F77_FUNC_(lis_vector_scatter_f, LIS_VECTOR_SCATTER_F) 
#define lis_vector_gather_f F77_FUNC_(lis_vector_gather_f, LIS_VECTOR_GATHER_F) 
#define lis_vector_print_f F77_FUNC_(lis_vector_print_f, LIS_VECTOR_PRINT_F) 
#define lis_vector_swap_f F77_FUNC_(lis_vector_swap_f, LIS_VECTOR_SWAP_F) 
#define lis_vector_copy_f F77_FUNC_(lis_vector_copy_f, LIS_VECTOR_COPY_F) 
#define lis_vector_axpy_f F77_FUNC_(lis_vector_axpy_f, LIS_VECTOR_AXPY_F) 
#define lis_vector_xpay_f F77_FUNC_(lis_vector_xpay_f, LIS_VECTOR_XPAY_F) 
#define lis_vector_axpyz_f F77_FUNC_(lis_vector_axpyz_f, LIS_VECTOR_AXPYZ_F) 
#define lis_vector_scale_f F77_FUNC_(lis_vector_scale_f, LIS_VECTOR_SCALE_F) 
#define lis_vector_pmul_f F77_FUNC_(lis_vector_pmul_f, LIS_VECTOR_PMUL_F) 
#define lis_vector_pdiv_f F77_FUNC_(lis_vector_pdiv_f, LIS_VECTOR_PDIV_F) 
#define lis_vector_set_all_f F77_FUNC_(lis_vector_set_all_f, LIS_VECTOR_SET_ALL_F) 
#define lis_vector_abs_f F77_FUNC_(lis_vector_abs_f, LIS_VECTOR_ABS_F) 
#define lis_vector_reciprocal_f F77_FUNC_(lis_vector_reciprocal_f, LIS_VECTOR_RECIPROCAL_F)
#define lis_vector_conjugate_f F77_FUNC_(lis_vector_conjugate_f, LIS_VECTOR_CONJUGATE_F) 
#define lis_vector_shift_f F77_FUNC_(lis_vector_shift_f, LIS_VECTOR_SHIFT_F) 
#define lis_vector_dot_f F77_FUNC_(lis_vector_dot_f, LIS_VECTOR_DOT_F)
#define lis_vector_nhdot_f F77_FUNC_(lis_vector_nhdot_f, LIS_VECTOR_NHDOT_F) 
#define lis_vector_nrm1_f F77_FUNC_(lis_vector_nrm1_f, LIS_VECTOR_NRM1_F) 
#define lis_vector_nrm2_f F77_FUNC_(lis_vector_nrm2_f, LIS_VECTOR_NRM2_F) 
#define lis_vector_nrmi_f F77_FUNC_(lis_vector_nrmi_f, LIS_VECTOR_NRMI_F) 
#define lis_vector_sum_f F77_FUNC_(lis_vector_sum_f, LIS_VECTOR_SUM_F) 
#define lis_vector_is_null_f F77_FUNC_(lis_vector_is_null_f, LIS_VECTOR_IS_NULL_F) 

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef __cplusplus
}
#endif

#endif
