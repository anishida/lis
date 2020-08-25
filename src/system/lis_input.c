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
#include <ctype.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_input_option
 * lis_input
 * lis_input_matrix
 * lis_input_vector
 * lis_input_vector_mm
 * lis_input_vector_plain
 * lis_input_vector_lis
 * lis_input_vector_lis_ascii
 * lis_fscan_scalar
 ************************************************/

LIS_INT lis_input_option(LIS_MATRIX_INPUT_OPTION *option, char *path);

#undef __FUNC__
#define __FUNC__ "lis_input"
LIS_INT lis_input(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, char *filename)
{
	LIS_INT	err;
	LIS_INT	fileformat;
	char buf[256],banner[128];
	FILE *file;

	LIS_DEBUG_FUNC_IN;

	err = lis_matrix_check(A,LIS_MATRIX_CHECK_NULL);
	if( err ) return err;
	if( b!=NULL && x!=NULL )
	{
		err = lis_vector_check(b,LIS_VECTOR_CHECK_NULL);
		if( err ) return err;
		err = lis_vector_check(x,LIS_VECTOR_CHECK_NULL);
		if( err ) return err;
	}

	if( filename==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"filname is NULL\n");
		return LIS_ERR_ILL_ARG;
	}
	file = fopen(filename, "r");
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n",filename);
		return LIS_ERR_FILE_IO;
	}

	/* file format check */
	if( fgets(buf, 256, file) == NULL )
	{
		fclose(file);
		return LIS_ERR_FILE_IO;
	}
	sscanf(buf, "%s", banner);
	if( strncmp(banner, MM_BANNER, strlen(MM_BANNER)) == 0)
	{
		fileformat = LIS_FMT_MM;
	}
/*	else if( strncmp(banner, LISBanner, strlen(LISBanner)) == 0)
	{
		fileformat = LIS_FMT_LIS;
	}
	else if( strncmp(banner, ITBLBanner, strlen(ITBLBanner)) == 0)
	{
		fileformat = LIS_FMT_ITBL;
	}
*/
	else
	{
		fileformat = LIS_FMT_HB;
	}
	rewind(file);

/*
	if( fileformat==LIS_FMT_FREE )
	{
		fclose(file);
		err = lis_input_option(&option, filename);
		if( err ) return err;
		file = fopen(option.filename, "r");
		if( file==NULL )
		{
			LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n",filename);
			return LIS_ERR_FILE_IO;
		}
	}
*/

	switch( fileformat )
	{
	case LIS_FMT_MM:
		err = lis_input_mm(A,b,x,file);
		break;
	case LIS_FMT_HB:
		err = lis_input_hb(A,b,x,file);
		break;
/*
	case LIS_FMT_ITBL:
		err = lis_input_mmm(A,b,x,file,comm,matrix_type,bnr,bnc,row,col);
		break;
	case LIS_FMT_LIS:
		err = lis_input_lis(A,b,x,filename,file,comm,matrix_type,bnr,bnc,row,col);
		break;
	case LIS_FMT_FREE:
		err = lis_input_free(A,b,x,option,file,comm,matrix_type,bnr,bnc,row,col);
		break;
*/
	default:
		fclose(file);
		return err;
	}
	fclose(file);
#ifdef USE_MPI
	MPI_Barrier(A->comm);
#endif

	LIS_DEBUG_FUNC_OUT;
	return err;
}


#undef __FUNC__
#define __FUNC__ "lis_input_matrix"
LIS_INT lis_input_matrix(LIS_MATRIX A, char *filename)
{
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_input(A,NULL,NULL,filename);
	if( err ) return err;

	LIS_DEBUG_FUNC_OUT;
	return err;
}


#undef __FUNC__
#define __FUNC__ "lis_input_vector"
LIS_INT lis_input_vector(LIS_VECTOR v, char *filename)
{
	LIS_INT	fileformat;
	char buf[256],banner[128];
	LIS_INT err;
	FILE *file;
	LIS_Comm comm;

	comm = v->comm;
	if( filename==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"filname is NULL\n");
		return LIS_ERR_ILL_ARG;
	}
	file = fopen(filename, "r");
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n",filename);
		return LIS_ERR_FILE_IO;
	}

	if( fgets(buf, 256, file) == NULL )
	{
		fclose(file);
		return LIS_ERR_FILE_IO;
	}
	sscanf(buf, "%s", banner);
	if( strncmp(banner, MM_BANNER, strlen(MM_BANNER)) == 0)
	{
		fileformat = LIS_FMT_MM;
	}
	else if( strncmp(banner, LISBanner, strlen(LISBanner)) == 0)
	{
		fileformat = LIS_FMT_LIS;
	}
	else
	{
		fileformat = LIS_FMT_PLAIN;
	}
	rewind(file);

	switch( fileformat )
	{
	case LIS_FMT_MM:
		err = lis_input_vector_mm(v,file);
		break;
	case LIS_FMT_LIS:
		err = lis_input_vector_lis(v,filename,file);
		break;
	case LIS_FMT_PLAIN:
		err = lis_input_vector_plain(v,file);
		break;
	}
	fclose(file);
#ifdef USE_MPI
	MPI_Barrier(comm);
#endif
	return err;
}


#undef __FUNC__
#define __FUNC__ "lis_input_vector_mm"
LIS_INT lis_input_vector_mm(LIS_VECTOR v, FILE *file)
{
	char buf[BUFSIZE];
	char banner[64], mtx[64], fmt[64], dtype[64], dstruct[64];
	char *p;
	LIS_INT i;
	LIS_INT	err;
	LIS_INT	n,is,ie;
	LIS_INT	idx;
	LIS_INT mmtype;	
	double val;
	double re,im;	


	/* check banner */
	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	sscanf(buf, "%s %s %s %s %s", banner, mtx, fmt, dtype, dstruct);

	for(p=mtx;*p!='\0';p++)     *p = (char)tolower(*p);
	for(p=fmt;*p!='\0';p++)     *p = (char)tolower(*p);
	for(p=dtype;*p!='\0';p++)   *p = (char)tolower(*p);
	for(p=dstruct;*p!='\0';p++) *p = (char)tolower(*p);

	if( strncmp(banner, MM_BANNER, strlen(MM_BANNER))!=0 || strncmp(mtx, MM_VEC, strlen(MM_VEC))!=0 )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not Matrix Market banner\n");
		return LIS_ERR_FILE_IO;
	}
	if( strncmp(fmt, MM_FMT_COO, strlen(MM_FMT_COO))!=0 && strncmp(fmt, MM_FMT_DNS, strlen(MM_FMT_DNS))!=0 )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not Matrix Market format\n");
		return LIS_ERR_FILE_IO;
	}
	if( strncmp(dtype, MM_TYPE_REAL, strlen(MM_TYPE_REAL))==0 )
	{
		mmtype = MM_REAL;
	}
#ifdef _COMPLEX	
	else if( strncmp(dtype, MM_TYPE_COMPLEX, strlen(MM_TYPE_COMPLEX))==0 )
	{
		mmtype = MM_COMPLEX;
	}
	else
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not real or complex\n");
		return LIS_ERR_FILE_IO;
	}
#else
	else
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not real\n");
		return LIS_ERR_FILE_IO;
	}
#endif	

	if( strncmp(dstruct, MM_TYPE_GENERAL, strlen(MM_TYPE_GENERAL))!=0 )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not general\n");
		return LIS_ERR_FILE_IO;
	}

	/* check size */		
	do
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
	}while( buf[0]=='%' );
#ifdef _LONG__LONG
	if( sscanf(buf, "%lld", &n) != 1 )
#else
	if( sscanf(buf, "%d", &n) != 1 )
#endif
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}

	/* read data */
	err = lis_vector_set_size(v,0,n);
	if( err )
	{
		return err;
	}
	lis_vector_get_range(v,&is,&ie);

	for(i=0;i<n;i++)
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
#ifdef _COMPLEX		
#ifdef _LONG__LONG
		if( mmtype==MM_REAL && sscanf(buf, "%lld %lg", &idx, &re) != 2 )
#else
		if( mmtype==MM_REAL && sscanf(buf, "%d %lg", &idx, &re) != 2 )
#endif
#else
#ifdef _LONG__LONG
		if( mmtype==MM_REAL && sscanf(buf, "%lld %lg", &idx, &val) != 2 )
#else
		if( mmtype==MM_REAL && sscanf(buf, "%d %lg", &idx, &val) != 2 )
#endif
#endif		  
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
#ifdef _LONG__LONG
		if( mmtype==MM_COMPLEX && sscanf(buf, "%lld %lg %lg", &idx, &re, &im) != 3 )
#else
		  if( mmtype==MM_COMPLEX && sscanf(buf, "%d %lg %lg", &idx, &re, &im) != 3 )
#endif
		  
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
		idx--;
		if( idx>=is && idx<ie )
		{
			if( mmtype==MM_REAL )
			{
#ifdef _COMPLEX	  
				v->value[idx-is] = re;
#else
				v->value[idx-is] = val;
#endif
			}
#ifdef _COMPLEX	  
			else
			{
				v->value[idx-is] = re + im * _Complex_I;
			}	  
#endif
		}
	}
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_vector_plain"
LIS_INT lis_input_vector_plain(LIS_VECTOR v, FILE *file)
{
	char buf[BUFSIZE];
	LIS_INT i;
	LIS_INT	err;
	LIS_INT	n,is,ie;
	double val;


	/* check size */
	n = 0;
	do
	{
		err = fscanf(file, "%lg", &val);
		if( err==1 ) n++;
	}while( err==1 );
	rewind(file);

	/* read data */
	err = lis_vector_set_size(v,0,n);
	if( err )
	{
		return err;
	}
	lis_vector_get_range(v,&is,&ie);

	for(i=0;i<n;i++)
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
		if( i>=is && i<ie )
		{
			if( sscanf(buf, "%lg", &val) != 1 )
			{
				LIS_SETERR_FIO;
				return LIS_ERR_FILE_IO;
			}
			v->value[i-is] = val;
		}
	}
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_vector_lis"
LIS_INT lis_input_vector_lis(LIS_VECTOR v, char *filename, FILE *file)
{
	LIS_INT err;
	char buf[BUFSIZE],banner[128],mode[128],mattype[128];
	LIS_INT in_mode;

	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	buf[10] = '\0';
	sscanf(buf, "%s %s %s", banner, mode, mattype);
	if( strncmp(banner, LISBanner, strlen(LISBanner)) != 0)
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"not lis file format\n");
		return LIS_ERR_FILE_IO;
	}

	in_mode = LIS_FMT_LIS_ASCII;
	if( mode[0]=='B' || mode[0]=='L' )
	{
		fclose(file);
		file = fopen(filename, "rb");
		if( file==NULL )
		{
			LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", filename);
			return LIS_ERR_FILE_IO;
		}
		err = fread(buf, sizeof(char), 10, file);
		if( err )
		  {
		    return err;
		  }
		in_mode = 1;
		in_mode = *(char *)&in_mode;
		if( (in_mode==LIS_BINARY_BIG && mode[0]=='L') || (in_mode==LIS_BINARY_LITTLE && mode[0]=='B') )
		{
			in_mode = LIS_TRUE;
		}
		else
		{
			in_mode = LIS_FALSE;
		}
	}

	if( strncmp(mattype, "vec", 3) == 0 )
	{
		if( in_mode==LIS_FMT_LIS_ASCII )
		{
			lis_input_vector_lis_ascii(v,file);
		}
		else
		{
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
		}
	}
	else
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"not lis file format\n");
		return LIS_ERR_FILE_IO;
	}

	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_vector_lis_ascii"
LIS_INT lis_input_vector_lis_ascii(LIS_VECTOR v, FILE *file)
{
	LIS_INT nprocs,my_rank,pe,n;  
	int int_nprocs,int_my_rank;
	LIS_INT	err;
	LIS_INT	ibuf[10];
	char cbuf[BUFSIZE];
	char c;
	double val;
	LIS_Comm comm;

	comm = v->comm;
#ifdef USE_MPI
	MPI_Comm_size(comm,&int_nprocs);
	MPI_Comm_rank(comm,&int_my_rank);
	nprocs = int_nprocs;
	my_rank = int_my_rank;
#else
	nprocs  = 1;
	my_rank = 0;
#endif

	if( fgets(cbuf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
#ifdef _LONG__LONG
	if( sscanf(cbuf, "%lld",&ibuf[0])!=1 )
#else
	if( sscanf(cbuf, "%d",&ibuf[0])!=1 )
#endif
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	if( nprocs!=ibuf[0] )
	{
		LIS_SETERR2(LIS_ERR_FILE_IO,"The number of PE=(%D) is different (in file PE=%D)\n",nprocs,ibuf[0]);
		return LIS_ERR_FILE_IO;
	}
	pe=-1;
	do
	{
		if( fgets(cbuf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
		if( cbuf[0]=='#' )
		{
#ifdef _LONG__LONG
			if( sscanf(cbuf, "%c %lld %lld",&c, &pe, &ibuf[1])!=3 )
#else
			if( sscanf(cbuf, "%c %d %d",&c, &pe, &ibuf[1])!=3 )
#endif
			{
				LIS_SETERR_FIO;
				return LIS_ERR_FILE_IO;
			}
		}
	}while(pe!=my_rank);

	n   = ibuf[1];
	err = lis_vector_set_size(v,0,n);
	if( err )
	{
		return err;
	}

	err = lis_fscan_scalar(n,file,v->value);
	if( err )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}

	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_fscan_scalar"
LIS_INT lis_fscan_scalar(LIS_INT n, FILE *file, LIS_SCALAR val[])
{
	LIS_INT err;
	LIS_INT i;
	double double_value;

	i=0;
	while( i<n )
	{
		err = fscanf(file, "%lg", &double_value);
		val[i++] = double_value;
		if( err )
		  {
		    return err;
		  }
	}
	return LIS_SUCCESS;
}
