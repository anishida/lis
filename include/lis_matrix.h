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


#ifndef __LIS_MATRIX_H__
#define __LIS_MATRIX_H__


#define LIS_MATRIX_CSR_STR "csr"
#define LIS_MATRIX_CSC_STR "csc"
#define LIS_MATRIX_MSR_STR "msr"
#define LIS_MATRIX_DIA_STR "dia"
#define LIS_MATRIX_ELL_STR "ell"
#define LIS_MATRIX_JAD_STR "jad"
#define LIS_MATRIX_BSR_STR "bsr"
#define LIS_MATRIX_BSC_STR "bsc"
#define LIS_MATRIX_VBR_STR "vbr"
#define LIS_MATRIX_DNS_STR "dns"
#define LIS_MATRIX_COO_STR "coo"
#define LIS_MATRIX_TJD_STR "tjd"

#define LIS_MATRIX_CHECK_ALL 0
#define LIS_MATRIX_CHECK_SIZE 1
#define LIS_MATRIX_CHECK_NULL 2
#define LIS_MATRIX_CHECK_TYPE 3
#define LIS_MATRIX_CHECK_NOT_ASSEMBLED 4
#define LIS_MATRIX_CHECK_SET 5

#define LIS_MATRIX_OPTION_CALL_BY 0

#define LIS_CALL_BY_REFERENCE 0
#define LIS_CALL_BY_VALUE 1

#define LIS_MATRIX_W_ANNZ 10

#ifdef __cplusplus
extern "C"
{
#endif
	extern LIS_INT lis_matrix_init(LIS_MATRIX *Amat);
	extern LIS_INT lis_matrix_check(LIS_MATRIX A, LIS_INT level);
	extern LIS_INT lis_matrix_storage_destroy(LIS_MATRIX Amat);
	extern LIS_INT lis_matrix_set_destroyflag(LIS_MATRIX A, LIS_INT flag);
	extern LIS_INT lis_matrix_get_destroyflag(LIS_MATRIX A, LIS_INT *flag);
	extern LIS_INT lis_matrix_LU_create(LIS_MATRIX A);
	extern LIS_INT lis_matrix_DLU_destroy(LIS_MATRIX Amat);
	extern LIS_INT lis_matrix_unset(LIS_MATRIX A);
	extern LIS_INT lis_matrix_copy_struct(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_self(LIS_SOLVER solver);

/*******************/
/* Operations      */
/*******************/

	extern LIS_INT lis_matrix_split(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_create(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_update(LIS_MATRIX A);
	extern LIS_INT lis_matrix_merge(LIS_MATRIX A);
	extern LIS_INT lis_matrix_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT flag);
	extern LIS_INT lis_matrix_copyDLU(LIS_MATRIX Ain, LIS_MATRIX_DIAG *D, LIS_MATRIX *L, LIS_MATRIX *U);
	extern LIS_INT lis_matrix_axpy(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B);
	extern LIS_INT lis_matrix_xpay(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B);
	extern LIS_INT lis_matrix_axpyz(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B, LIS_MATRIX C);  
	extern LIS_INT lis_matrix_axpy_dns(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B);
	extern LIS_INT lis_matrix_xpay_dns(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B);
	extern LIS_INT lis_matrix_axpyz_dns(LIS_SCALAR alpha, LIS_MATRIX A, LIS_MATRIX B, LIS_MATRIX C);  
        extern LIS_INT lis_matrix_shift_diagonal(LIS_MATRIX A, LIS_SCALAR sigma);
        extern LIS_INT lis_matrix_shift_diagonal_csr(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_csc(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_msr(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_dia(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_ell(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_jad(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_bsr(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_bsc(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_dns(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_coo(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_diagonal_vbr(LIS_MATRIX A, LIS_SCALAR sigma);
	extern LIS_INT lis_matrix_shift_matrix(LIS_MATRIX A, LIS_MATRIX B, LIS_SCALAR sigma);  

/*******************/
/* Diagonal Matrix */
/*******************/

	extern LIS_INT lis_matrix_diag_init(LIS_MATRIX_DIAG *D);
	extern LIS_INT lis_matrix_diag_check(LIS_MATRIX_DIAG D, LIS_INT level);
	extern LIS_INT lis_matrix_diag_create(LIS_INT local_n, LIS_INT global_n, LIS_Comm comm, LIS_MATRIX_DIAG *D);
	extern LIS_INT lis_matrix_diag_destroy(LIS_MATRIX_DIAG D);
	extern LIS_INT lis_matrix_diag_duplicate(LIS_MATRIX_DIAG Din, LIS_MATRIX_DIAG *Dout);
	extern LIS_INT lis_matrix_diag_duplicateM(LIS_MATRIX Ain, LIS_MATRIX_DIAG *Dout);
	extern LIS_INT lis_matrix_diag_get_range(LIS_MATRIX_DIAG D, LIS_INT *is, LIS_INT *ie);
	extern LIS_INT lis_matrix_diag_get_size(LIS_MATRIX_DIAG D, LIS_INT *local_n, LIS_INT *global_n);
	extern LIS_INT lis_matrix_diag_set_blocksize(LIS_MATRIX_DIAG D, LIS_INT bn, LIS_INT *bns);
	extern LIS_INT lis_matrix_diag_copy(LIS_MATRIX_DIAG X, LIS_MATRIX_DIAG Y);
	extern LIS_INT lis_matrix_diag_scale(LIS_SCALAR alpha, LIS_MATRIX_DIAG D);
	extern LIS_INT lis_matrix_diag_inverse(LIS_MATRIX_DIAG D);
	extern LIS_INT lis_matrix_diag_prLIS_INT(LIS_MATRIX_DIAG D);
	extern LIS_INT lis_matrix_diag_mallocM(LIS_MATRIX A, LIS_SCALAR **diag);
	extern LIS_INT lis_matrix_diag_matvec(LIS_MATRIX_DIAG D, LIS_VECTOR X, LIS_VECTOR Y);
	extern LIS_INT lis_matrix_diag_matvech(LIS_MATRIX_DIAG D, LIS_VECTOR X, LIS_VECTOR Y);

/*******************/
/* CSR             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_csr(LIS_INT nnzl, LIS_INT nnzu, LIS_SCALAR *diag, LIS_INT *lptr, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uptr, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_csr(LIS_INT n, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_csr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_copyDLU_csr(LIS_MATRIX Ain, LIS_MATRIX_DIAG *D, LIS_MATRIX *L, LIS_MATRIX *U);
	extern LIS_INT lis_matrix_get_diagonal_csr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_csr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_csr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_sort_csr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_normf_csr(LIS_MATRIX A, LIS_SCALAR *nrm);
	extern LIS_INT lis_matrix_transpose_csr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_split_csr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_create_csr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_update_csr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_merge_csr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_solve_csr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_csr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_split2_csr(LIS_MATRIX A);

/*******************/
/* CSC             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_csc(LIS_INT nnzl, LIS_INT nnzu, LIS_SCALAR *diag, LIS_INT *lptr, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uptr, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_csc(LIS_INT n, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_csc(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_get_diagonal_csc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_csc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_csc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_sort_csc(LIS_MATRIX A);
	extern LIS_INT lis_matrix_normf_csc(LIS_MATRIX A, LIS_SCALAR *nrm);
	extern LIS_INT lis_matrix_transpose_csc(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_split_csc(LIS_MATRIX A);
	extern LIS_INT lis_matrix_merge_csc(LIS_MATRIX A);
	extern LIS_INT lis_matrix_solve_csc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_csc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_convert_csr2csc(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_csc2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* BSR             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_bsr(LIS_INT bnr, LIS_INT bnc, LIS_INT lbnnz, LIS_INT ubnnz, LIS_MATRIX_DIAG D, LIS_INT *lbptr, LIS_INT *lbindex, LIS_SCALAR *lvalue, LIS_INT *ubptr, LIS_INT *ubindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_bsr(LIS_INT n, LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_bsr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_csr2bsr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_bsr2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_get_diagonal_bsr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_bsr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_bsr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_bscale_bsr(LIS_MATRIX A, LIS_MATRIX_DIAG D);
	extern LIS_INT lis_matrix_split_bsr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_merge_bsr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_solve_bsr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_bsr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_sort_bsr(LIS_MATRIX A);

/*******************/
/* MSR             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_msr(LIS_INT lnnz, LIS_INT unnz, LIS_INT lndz, LIS_INT undz, LIS_SCALAR *diag, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_msr(LIS_INT n, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_msr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_merge_msr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_msr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_get_diagonal_msr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_msr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_msr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_solve_msr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_msr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_convert_csr2msr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_msr2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* ELL             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_ell(LIS_INT lmaxnzr, LIS_INT umaxnzr, LIS_SCALAR *diag, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_ell(LIS_INT n, LIS_INT maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_ell(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_merge_ell(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_ell(LIS_MATRIX A);
	extern LIS_INT lis_matrix_get_diagonal_ell(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_ell(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_ell(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_solve_ell(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_ell(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_convert_csr2ell(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_ell2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* JAD             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_jad(LIS_INT lnnz, LIS_INT unnz, LIS_INT lmaxnzr, LIS_INT umaxnzr, LIS_SCALAR *diag, LIS_INT *lperm, LIS_INT *lptr, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uperm, LIS_INT *uptr, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_jad(LIS_INT n, LIS_INT maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_perm, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_jad(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_merge_jad(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_jad(LIS_MATRIX A);
	extern LIS_INT lis_matrix_get_diagonal_jad(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_jad(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_jad(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_solve_jad(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_jad(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_convert_csr2jad(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_jad2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* DIA             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_dia(LIS_INT lnnd, LIS_INT unnd, LIS_SCALAR *diag, LIS_INT *lindex, LIS_SCALAR *lvalue, LIS_INT *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_dia(LIS_INT n, LIS_INT nnd, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_dia(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_merge_dia(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_dia(LIS_MATRIX A);
	extern LIS_INT lis_matrix_get_diagonal_dia(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_dia(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_dia(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_solve_dia(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_dia(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_convert_csr2dia(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_dia2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* BSC             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_bsc(LIS_INT bnr, LIS_INT bnc, LIS_INT lbnnz, LIS_INT ubnnz, LIS_MATRIX_DIAG D, LIS_INT *lbptr, LIS_INT *lbindex, LIS_SCALAR *lvalue, LIS_INT *ubptr, LIS_INT *ubindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_bsc(LIS_INT n, LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *o_ptr, LIS_INT *o_index, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_bsc(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_csc2bsc(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_bsc2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_get_diagonal_bsc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_bsc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_bsc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_split_bsc(LIS_MATRIX A);
	extern LIS_INT lis_matrix_merge_bsc(LIS_MATRIX A);
	extern LIS_INT lis_matrix_solve_bsc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_bsc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);

/*******************/
/* VBR             */
/*******************/

	extern LIS_INT lis_matrix_elements_copy_vbr(LIS_INT n, LIS_INT nr, LIS_INT nc, LIS_INT bnnz, LIS_INT *row, LIS_INT *col, LIS_INT *ptr, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_INT *o_row, LIS_INT *o_col, LIS_INT *o_ptr, LIS_INT *o_bptr, LIS_INT *o_bindex, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_vbr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_csr2vbr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_vbr2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_get_diagonal_vbr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_vbr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_vbr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_split_vbr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_merge_vbr(LIS_MATRIX A);
	extern LIS_INT lis_matrix_solve_vbr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_vbr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);

/*******************/
/* DNS             */
/*******************/

	extern LIS_INT lis_matrix_setDLU_dns(LIS_SCALAR *diag, LIS_SCALAR *lvalue, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern LIS_INT lis_matrix_elements_copy_dns(LIS_INT n, LIS_INT gn, LIS_SCALAR *value, LIS_SCALAR *o_value);
	extern LIS_INT lis_matrix_copy_dns(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_merge_dns(LIS_MATRIX A);
	extern LIS_INT lis_matrix_split_dns(LIS_MATRIX A);
	extern LIS_INT lis_matrix_get_diagonal_dns(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_dns(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_dns(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_solve_dns(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_solveh_dns(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, LIS_INT flag);
	extern LIS_INT lis_matrix_convert_csr2dns(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_dns2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* COO             */
/*******************/

	extern LIS_INT lis_matrix_copy_coo(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_get_diagonal_coo(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_coo(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_scale_symm_coo(LIS_MATRIX A, LIS_SCALAR d[]);
	extern LIS_INT lis_matrix_sort_coo(LIS_MATRIX A);
	extern LIS_INT lis_matrix_normf_coo(LIS_MATRIX A, LIS_SCALAR *nrm);
	extern LIS_INT lis_matrix_transpose_coo(LIS_MATRIX Ain, LIS_MATRIX *Aout);
	extern LIS_INT lis_matrix_split_coo(LIS_MATRIX A);
	extern LIS_INT lis_matrix_merge_coo(LIS_MATRIX A);
	extern LIS_INT lis_matrix_convert_csr2coo(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_coo2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* RCO             */
/*******************/

	extern LIS_INT lis_matrix_create_rco(LIS_INT local_n, LIS_INT global_n, LIS_Comm comm, LIS_INT annz, LIS_INT *nnz, LIS_MATRIX *Amat);
	extern LIS_INT lis_matrix_malloc_rco(LIS_INT n, LIS_INT nnz[], LIS_INT **row, LIS_INT ***index, LIS_SCALAR ***value);
	extern LIS_INT lis_matrix_realloc_rco(LIS_INT row, LIS_INT nnz, LIS_INT ***index, LIS_SCALAR ***value);
	extern LIS_INT lis_matrix_convert_rco2csr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_rco2bsr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_convert_rco2csc(LIS_MATRIX Ain, LIS_MATRIX Aout);

/*******************/
/* ILU             */
/*******************/

	extern LIS_INT lis_matrix_ilu_create(LIS_INT n, LIS_INT bs, LIS_MATRIX_ILU *A);
	extern LIS_INT lis_matrix_ilu_setCR(LIS_MATRIX_ILU A);
	extern LIS_INT lis_matrix_ilu_setVR(LIS_MATRIX_ILU A);
	extern LIS_INT lis_matrix_ilu_destroy(LIS_MATRIX_ILU A);
	extern LIS_INT lis_matrix_ilu_premalloc(LIS_INT nnzrow, LIS_MATRIX_ILU A);
	extern LIS_INT lis_matrix_ilu_realloc(LIS_INT row, LIS_INT nnz, LIS_MATRIX_ILU A);
	extern LIS_INT lis_matvec_ilu(LIS_MATRIX A, LIS_MATRIX_ILU LU, LIS_VECTOR X, LIS_VECTOR Y);
	extern LIS_INT lis_matvech_ilu(LIS_MATRIX A, LIS_MATRIX_ILU LU, LIS_VECTOR X, LIS_VECTOR Y);


#ifdef __cplusplus
}
#endif

#endif
