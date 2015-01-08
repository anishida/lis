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


#ifndef __LIS_VECTOR_H__
#define __LIS_VECTOR_H__

#define LIS_VECTOR_CHECK_ALL 0
#define LIS_VECTOR_CHECK_NULL 1
#define LIS_VECTOR_CHECK_ASSEMBLE 2



#ifdef __cplusplus
extern "C"
{
#endif

	extern LIS_INT lis_vector_init(LIS_VECTOR *vec);
	extern LIS_INT lis_vector_createex(LIS_INT precision, LIS_Comm comm, LIS_VECTOR *vec);
	extern LIS_INT lis_vector_reuse(LIS_VECTOR *vec);
	extern LIS_INT lis_vector_duplicateex(LIS_INT precision, void *A, LIS_VECTOR *vout);
	extern LIS_INT lis_vector_unset(LIS_VECTOR vec);
	extern LIS_INT lis_vector_set(LIS_VECTOR vec, LIS_SCALAR *value);
	extern LIS_INT lis_vector_set_destroyflag(LIS_VECTOR v, LIS_INT flag);
	extern LIS_INT lis_vector_get_destroyflag(LIS_VECTOR v, LIS_INT *flag);
	extern LIS_INT lis_vector_check(LIS_VECTOR v, LIS_INT level);

/*******************/
/* Operations      */
/*******************/

	extern LIS_INT lis_vector_axpyex_mmm(LIS_QUAD_PTR alpha, LIS_VECTOR vx, LIS_VECTOR vy);
	extern LIS_INT lis_vector_xpayex_mmm(LIS_VECTOR vx, LIS_QUAD_PTR alpha, LIS_VECTOR vy);
	extern LIS_INT lis_vector_axpyzex_mmmm(LIS_QUAD_PTR alpha, LIS_VECTOR vx, LIS_VECTOR vy, LIS_VECTOR vz);
	extern LIS_INT lis_vector_copyex_mm(LIS_VECTOR vx, LIS_VECTOR vy);
	extern LIS_INT lis_vector_copyex_mn(LIS_VECTOR vx, LIS_VECTOR vy);
	extern LIS_INT lis_vector_copyex_nm(LIS_VECTOR vx, LIS_VECTOR vy);
	extern LIS_INT lis_vector_scaleex_mm(LIS_QUAD_PTR alpha, LIS_VECTOR vx);
	extern LIS_INT lis_vector_scaleex_nm(LIS_SCALAR alpha, LIS_VECTOR vx);
	extern LIS_INT lis_vector_set_allex_nm(LIS_SCALAR alpha, LIS_VECTOR vx);
	extern LIS_INT lis_vector_nrm2ex_mm(LIS_VECTOR vx, LIS_QUAD_PTR *value);
	extern LIS_INT lis_vector_dotex_mmm(LIS_VECTOR vx, LIS_VECTOR vy, LIS_QUAD_PTR *value);
	extern LIS_INT lis_vector_reciprocalex_m(LIS_VECTOR vx);

#ifdef __cplusplus
}
#endif

#endif
