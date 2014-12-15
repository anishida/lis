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


#ifndef __LIS_MPI_H__
#define __LIS_MPI_H__




#ifdef __cplusplus
extern "C"
{
#endif

	extern LIS_INT lis_commtable_create(LIS_MATRIX A);
	extern void lis_commtable_destroy(LIS_COMMTABLE table);
	extern LIS_INT lis_commtable_duplicate(LIS_COMMTABLE tin, LIS_COMMTABLE *tout);
	extern LIS_INT lis_commtable_duplicateM(LIS_MATRIX Ain, LIS_MATRIX *Aout);
	extern LIS_INT lis_send_recv(LIS_COMMTABLE commtable, LIS_SCALAR x[]);
	extern LIS_INT lis_reduce(LIS_COMMTABLE commtable, LIS_SCALAR x[]);
	extern LIS_INT lis_send_recv_mp(LIS_COMMTABLE commtable, LIS_VECTOR X);
	extern LIS_INT lis_reduce_mp(LIS_COMMTABLE commtable, LIS_VECTOR X);
	extern LIS_INT lis_matrix_g2l(LIS_MATRIX A);
		extern LIS_INT lis_matrix_g2l_csr(LIS_MATRIX A);
		extern LIS_INT lis_matrix_g2l_rco(LIS_MATRIX A);

#ifdef __cplusplus
}
#endif

#endif
