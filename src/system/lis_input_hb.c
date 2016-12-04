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
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_input_hb
 * lis_input_hb_get_fmt
 * lis_input_hb_csr
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_input_hb"
LIS_INT lis_input_hb(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file)
{
	LIS_INT	err;
	LIS_INT	matrix_type;
	LIS_MATRIX B;

	LIS_DEBUG_FUNC_IN;

	matrix_type = A->matrix_type;

	err = lis_input_hb_csr(A,b,x,file);
	if( err ) return err;

	if( matrix_type!=LIS_MATRIX_CSR && matrix_type!=LIS_MATRIX_CSC )
	{
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,matrix_type);
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
		if( A->matrix_type==LIS_MATRIX_JAD )
		{
			A->work = (LIS_SCALAR *)lis_malloc(A->n*sizeof(LIS_SCALAR),"lis_input_hb::A->work");
			if( A->work==NULL )
			{
				LIS_SETERR_MEM(A->n*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
		}
	}


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_hb_get_fmt"
LIS_INT lis_input_hb_get_fmt(char *buf, LIS_INT size, LIS_INT *iter, LIS_INT *width)
{
	char tmp[64];
	char *p, *s, *t;

	LIS_DEBUG_FUNC_IN;

	strncpy(tmp, buf, size); tmp[size] = '\0';
	for(p=tmp;*p!='\0';p++)     *p = (char)tolower(*p);
	p = strchr(tmp, '(');
	if( p!=NULL )
	{
		s = p+1;
		p = strchr(s, ')'); *p = '\0';
		p = strchr(s, 'i');
		if( p==NULL )
		{
			p = strchr(s, 'e');
			if( p==NULL )
			{
				p = strchr(s, 'd');
				if( p==NULL ) return LIS_FAILS;
			}
			t = strchr(s, '.'); *t = '\0';
		}
		*p = '\0'; 
		*iter  = atoi(s);
		*width = atoi(p+1);
	}
	else
	{
		*iter  = 0;
		*width = 0;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_hb_csr"
LIS_INT lis_input_hb_csr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file)
{
	char buf[BUFSIZE];
	char title[128], key[128], mtx[64], dat[128];
	char *p;
	char MXTYPE_F,MXTYPE_S,MXTYPE_T;
	char RHSTYP_F,RHSTYP_S,RHSTYP_T;
	LIS_INT TOTCRD,PTRCRD,INDCRD,VALCRD,RHSCRD;
	LIS_INT	NROW,NCOL,NNZERO,NELTVL;
	LIS_INT	NRHS,NRHSIX;
	LIS_INT	iptr,iind,ival,irhs;
	LIS_INT	wptr,wind,wval,wrhs;
	LIS_INT	i,k,j,my_rank;
	int int_my_rank;
	LIS_INT	err;
	LIS_INT	n,is,ie;
	LIS_INT	*ptr, *index;
	LIS_INT	matrix_type;
	LIS_SCALAR *value;
	LIS_REAL re,im;
	LIS_MATRIX B;

	#ifdef USE_MPI
		MPI_Comm_rank(A->comm,&int_my_rank);
		my_rank = int_my_rank;
	#else
		my_rank = 0;
	#endif

	matrix_type = A->matrix_type;

	/* Line 1 */
	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	strncpy(title, buf    ,72); title[72] = '\0';
	strncpy(key  ,&buf[72], 8); key[8]    = '\0';
	/*
	printf("title: %s\n",title);
	printf("key  : %s\n",key);
	*/

	/* Line 2 */
	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	RHSCRD = 0;
#ifdef _LONG__LONG
	if( sscanf(buf, "%14lld%14lld%14lld%14lld%14lld", &TOTCRD, &PTRCRD, &INDCRD, &VALCRD, &RHSCRD) < 4 )
#else
	if( sscanf(buf, "%14d%14d%14d%14d%14d", &TOTCRD, &PTRCRD, &INDCRD, &VALCRD, &RHSCRD) < 4 )
#endif
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	/*
#ifdef _LONG__LONG
	printf("%14lld%14lld%14lld%14lld%14lld\n",TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD);
#else
	printf("%14d%14d%14d%14d%14d\n",TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD);
#endif
	*/

	/* Line 3 */
	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
#ifdef _LONG__LONG
	if( sscanf(buf, "%s %lld %lld %lld %lld", mtx, &NROW, &NCOL, &NNZERO, &NELTVL) != 5 )
#else
	if( sscanf(buf, "%s %d %d %d %d", mtx, &NROW, &NCOL, &NNZERO, &NELTVL) != 5 )
#endif
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	for(p=mtx;*p!='\0';p++)     *p = (char)tolower(*p);
	MXTYPE_F = mtx[0];
	MXTYPE_S = mtx[1];
	MXTYPE_T = mtx[2];
#ifdef _COMPLEX	
	if( mtx[0]!='r' && mtx[0]!='c' )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not real or complex\n");
		return LIS_ERR_FILE_IO;
	}
#else
	if( mtx[0]!='r' )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not real\n");
		return LIS_ERR_FILE_IO;
	}
#endif	
	if( mtx[1]!='u' )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not unsymmetric\n");
		return LIS_ERR_FILE_IO;
	}
	if( mtx[2]!='a' )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not assembled\n");
		return LIS_ERR_FILE_IO;
	}
	if( NROW!=NCOL )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"matrix is not square\n");
		return LIS_ERR_FILE_IO;
	}
#ifdef _LONG__LONG
	if( my_rank==0 ) printf("matrix size = %lld x %lld (%lld nonzero entries)\n\n",NROW,NCOL,NNZERO);
#else
	if( my_rank==0 ) printf("matrix size = %d x %d (%d nonzero entries)\n\n",NROW,NCOL,NNZERO);
#endif

	/* Line 4 */
	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	lis_input_hb_get_fmt( buf    ,16,&iptr,&wptr);
	lis_input_hb_get_fmt(&buf[16],16,&iind,&wind);
	lis_input_hb_get_fmt(&buf[32],20,&ival,&wval);
	lis_input_hb_get_fmt(&buf[52],20,&irhs,&wrhs);
	/*
#ifdef _LONG__LONG
	printf("%lld %lld %lld %lld\n",iptr,iind,ival,irhs);
	printf("%lld %lld %lld %lld\n",wptr,wind,wval,wrhs);
#else
	printf("%d %d %d %d\n",iptr,iind,ival,irhs);
	printf("%d %d %d %d\n",wptr,wind,wval,wrhs);
#endif
	*/

	/* Line 5 */
	if( RHSCRD!=0 )
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
#ifdef _LONG__LONG
		sscanf(buf, "%s %lld %lld", mtx, &NRHS, &NRHSIX);
#else
		sscanf(buf, "%s %d %d", mtx, &NRHS, &NRHSIX);
#endif
/*
#ifdef _LONG__LONG
		if( sscanf(buf, "%s %lld %lld", mtx, &NRHS, &NRHSIX) != 3 )
#else
		if( sscanf(buf, "%s %d %d", mtx, &NRHS, &NRHSIX) != 3 )
#endif
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
*/
		for(p=mtx;*p!='\0';p++)     *p = (char)tolower(*p);
		RHSTYP_F = mtx[0];
		RHSTYP_S = mtx[1];
		RHSTYP_T = mtx[2];
		/*
#ifdef _LONG__LONG
		printf("%c%c%c %lld %lld\n",RHSTYP_F, RHSTYP_S, RHSTYP_T, NRHS, NRHSIX);
#else
		printf("%c%c%c %d %d\n",RHSTYP_F, RHSTYP_S, RHSTYP_T, NRHS, NRHSIX);
#endif
		*/
	}

	err = lis_matrix_set_size(A,0,NROW);
	if( err )
	{
		return err;
	}
	n = A->n;
	lis_matrix_get_range(A,&is,&ie);
	err = lis_matrix_malloc_csr(n,NNZERO,&ptr,&index,&value);
	if( err )
	{
		return err;
	}

	/* read data */
	k = 0;
	for( i=0; i<PTRCRD; i++ )
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
		p = buf;
		for(j=0;j<iptr&&k<n+1;j++)
		{
			strncpy(dat, p, wptr); dat[wptr] = '\0';
			ptr[k] = atoi(dat) - 1;
			p += wptr;
			k++;
		}
	}

	k = 0;
	for( i=0; i<INDCRD; i++ )
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
		p = buf;
		for(j=0;j<iind&&k<NNZERO;j++)
		{
			strncpy(dat, p, wind); dat[wind] = '\0';
			index[k] = atoi(dat) - 1;
			p += wind;
			k++;
		}
	}

	k = 0;
	for( i=0; i<VALCRD; i++ )
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
		p = buf;
		for(j=0;j<ival&&k<NNZERO;j++)
		{
			if (MXTYPE_F == 'r')
			  {
			    strncpy(dat, p, wval); dat[wval] = '\0';
			    value[k] = atof(dat);
			    p += wval;
			    k++;
			  }
#ifdef _COMPLEX
			else if (MXTYPE_F == 'c')
			  {
			    strncpy(dat, p, wval); dat[wval] = '\0'; 
			    re = atof(dat);
			    p += wval;
			    strncpy(dat, p, wval); dat[wval] = '\0'; 
			    im = atof(dat);
			    p += wval;			    
			    value[k] = re + im * _Complex_I;
			    k++;
			    j++;
			  }
#endif
			else
			  {
			    LIS_SETERR_FIO;
			    return LIS_ERR_FILE_IO;
			  }
		}
	}

	if( RHSCRD>0 )
	{
		/*
		k = 0;
		for( i=0; i<RHSCRD; i++ )
		{
			if( fgets(buf, BUFSIZE, file) == NULL )
			{
				LIS_SETERR_FIO;
				return LIS_ERR_FILE_IO;
			}
			p = buf;
			for(j=0;j<ival&&k<NNZERO;j++)
			{
				strncpy(dat, p, wval); dat[wval] = '\0';
				value[k] = atof(dat);
				p += wval;
				printf("%e ",value[k]);
				k++;
			}
			printf("\n");
		}
		*/
	}
	err = lis_matrix_set_csc(NNZERO,ptr,index,value,A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_assemble(A);
	if( err ) return err;

	if( matrix_type!=LIS_MATRIX_CSC )
	{
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CSR);
		err = lis_matrix_convert_csc2csr(A,B);
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

	return LIS_SUCCESS;
}

