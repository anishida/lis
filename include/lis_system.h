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


#ifndef __LIS_SYSTEM_H__
#define __LIS_SYSTEM_H__

#include <stdarg.h>


#define LIS_SETERR(code,mess) lis_error(__FILE__,__FUNC__,__LINE__,code,mess)
#define LIS_SETERR1(code,mess,a1) lis_error(__FILE__,__FUNC__,__LINE__,code,mess,a1)
#define LIS_SETERR2(code,mess,a1,a2) lis_error(__FILE__,__FUNC__,__LINE__,code,mess,a1,a2)
#define LIS_SETERR3(code,mess,a1,a2,a3)	lis_error(__FILE__,__FUNC__,__LINE__,code,mess,a1,a2,a3)
#define LIS_SETERR4(code,mess,a1,a2,a3,a4) lis_error(__FILE__,__FUNC__,__LINE__,code,mess,a1,a2,a3,a4)
#define LIS_SETERR_MEM(sz) lis_error(__FILE__,__FUNC__,__LINE__,LIS_ERR_OUT_OF_MEMORY,"malloc size = %D\n",sz)
#define LIS_SETERR_IMP lis_error(__FILE__,__FUNC__,__LINE__,LIS_ERR_NOT_IMPLEMENTED,"not implemented\n")
#define LIS_SETERR_FIO lis_error(__FILE__,__FUNC__,__LINE__,LIS_ERR_FILE_IO,"file i/o error\n")


typedef struct LIS_ARGS_STRUCT
{
  struct LIS_ARGS_STRUCT *next, *prev;
  char *arg1;
  char *arg2;
} *LIS_ARGS;

typedef struct LIS_HASH_STRUCT
{
	struct LIS_HASH_STRUCT *next;
	LIS_INT	index;
	LIS_INT	value;
} *LIS_HASH;

typedef struct LIS_HASH_STRUCT **LIS_HASHTABLE;



#ifdef __cplusplus
extern "C"
{
#endif
	extern LIS_ARGS cmd_args;
	extern LIS_SCALAR *lis_vec_tmp;
	extern int lis_mpi_initialized;
#ifdef USE_MPI
	extern MPI_Op LIS_MPI_MSUM;
	extern MPI_Datatype LIS_MPI_MSCALAR;
#endif
	extern void lis_memory_init(void);
	extern void lis_free_mat(LIS_MATRIX A);
	extern LIS_INT lis_text2args(char *text, LIS_ARGS *args);
	extern LIS_INT lis_arg2args(LIS_INT argc, char *argv[], LIS_ARGS *args);
	extern LIS_INT lis_args_free(LIS_ARGS args);
	extern void lis_debug_set_comm(LIS_Comm comm);
	extern LIS_INT lis_replace(char *buf, const char *str1, const char *str2);
	extern LIS_INT lis_error(const char *file, const char *func, const LIS_INT line, const LIS_INT code, const char *mess, ...);
	extern LIS_INT lis_print_rhistory(LIS_Comm comm,LIS_INT iter, LIS_REAL resid);
	extern void lis_free_all(void);
	extern void lis_sort_i(LIS_INT is, LIS_INT ie, LIS_INT *i1);
	extern void lis_sort_ii(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2);
	extern void lis_sort_id(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_SCALAR *d1);
	extern void lis_sort_di(LIS_INT is, LIS_INT ie, LIS_SCALAR *d1, LIS_INT *i1);
  	extern void lis_sort_d(LIS_INT is, LIS_INT ie, LIS_SCALAR *d1);
	extern void lis_sort_dd(LIS_INT is, LIS_INT ie, LIS_SCALAR *d1, LIS_VECTOR *d2);
	extern void lis_sort_iid(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2, LIS_SCALAR *d1);
	extern void lis_sort_iiid(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2, LIS_INT *i3, LIS_SCALAR *d1);
	extern void lis_sortr_ii(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2);
	extern void lis_sort_jad(LIS_INT is, LIS_INT ie, LIS_INT maxnzr, LIS_INT *i1, LIS_INT *i2);
	extern void lis_sort_id_block(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_SCALAR *d1, LIS_INT bs);
	extern LIS_INT lis_bswap_int(LIS_INT n, LIS_INT *buf);
	extern LIS_INT lis_bswap_double(LIS_INT n, double *buf);
	extern LIS_INT lis_bswap_scalar(LIS_INT n, LIS_SCALAR *buf);	  
	extern LIS_INT lis_bswap_size_t(LIS_INT n, size_t *buf);
	extern LIS_INT lis_ranges_create(LIS_Comm comm, LIS_INT *local_n, LIS_INT *global_n, LIS_INT **ranges, LIS_INT *is, LIS_INT *ie, LIS_INT *nprocs, LIS_INT *my_rank);
	extern LIS_INT lis_hashtable_create(LIS_HASHTABLE *hashtable);
	extern LIS_INT lis_hashtable_destroy(LIS_HASHTABLE hashtable);
	extern LIS_INT lis_hashtable_clear(LIS_HASHTABLE hashtable);
	extern LIS_HASH lis_hashtable_search(LIS_HASHTABLE hashtable, LIS_INT index);
	extern LIS_INT lis_hashtable_set_value(LIS_HASHTABLE hashtable, LIS_INT index, LIS_INT value);
	extern LIS_INT lis_hashtable_get_value(LIS_HASHTABLE hashtable, LIS_INT index);

#ifdef __cplusplus
}
#endif

#endif
