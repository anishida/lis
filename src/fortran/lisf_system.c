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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_initialize_f
 * lis_finalize_f
 ************************************************/

#ifdef USE_FORTRAN

char **f_argv_tmp;
LIS_INT	f_argc_tmp;
LIS_Comm_f lis_comm_world_f = 0;

#undef __FUNC__
#define __FUNC__ "lis_check_comm_f"
void lis_set_comm_world_f(LIS_Comm_f *comm)
{
	LIS_DEBUG_FUNC_IN;

	lis_comm_world_f = *comm;

#ifdef _DEBUG
	  printf("c_comm = %d f_comm = %d\n",LIS_COMM_WORLD,*comm);
#endif

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_finitialize_f"
void lis_finitialize_f(LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_initialize(&f_argc_tmp,&f_argv_tmp);
	if( *ierr )	return;	
	/*
	#ifdef USE_MPI
		MPI_Init(&f_argc_tmp, &f_argv_tmp);
	#endif
	*/

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_finalize_f"
void lis_finalize_f(LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_finalize();
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_wtime_f"
double lis_wtime_f(void)
{
	LIS_DEBUG_FUNC_IN;

	return lis_wtime();

	LIS_DEBUG_FUNC_OUT;
}

#undef __FUNC__
#define __FUNC__ "CHKERR_f"
void CHKERR_f(LIS_INT *err)
{
	LIS_DEBUG_FUNC_IN;

	if(*err)
	{
		lis_finalize();
		exit(*err);
	}

	LIS_DEBUG_FUNC_OUT;
}


#undef __FUNC__
#define __FUNC__ "lis_set_argv_begin_f"
void lis_set_argv_begin_f(LIS_INT *argc, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	f_argc_tmp = *argc+1;
	f_argv_tmp = (char **)lis_malloc(f_argc_tmp*sizeof(char *),"lis_set_argv_begin_f::f_argv_tmp");
	if( f_argv_tmp==NULL )
	{
		LIS_SETERR_MEM(f_argc_tmp*sizeof(char *));
		*ierr = LIS_OUT_OF_MEMORY;
		return;
	}
	*ierr = LIS_SUCCESS;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_set_argv_f"
void lis_set_argv_f(LIS_INT *no, char *argv, LIS_INT *ierr, LIS_INT len)
{
	LIS_INT i;
	char *p;

	LIS_DEBUG_FUNC_IN;

	i = *no;
	f_argv_tmp[i] = (char *)lis_malloc((len+1)*sizeof(char),"lis_set_argv_f::f_argv_tmp");
	if( f_argv_tmp[i]==NULL )
	{
		LIS_SETERR_MEM((len+1)*sizeof(char));
		*ierr = LIS_OUT_OF_MEMORY;
		return;
	}
	memset(f_argv_tmp[i],0x20,(len+1)*sizeof(char));
	strncpy(f_argv_tmp[i],argv,len);
	p = &f_argv_tmp[i][len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';
	*ierr = LIS_SUCCESS;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_set_argv_end_f"
void lis_set_argv_end_f(LIS_INT *ierr)
{
	LIS_INT i;

	LIS_DEBUG_FUNC_IN;

	for(i=0;i<f_argc_tmp;i++)
	{
		lis_free(f_argv_tmp[i]);
	}
	lis_free(f_argv_tmp);
	f_argv_tmp = NULL;
	f_argc_tmp = 0;
	*ierr = LIS_SUCCESS;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_arg2args_f"
void lis_arg2args_f(LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_arg2args(f_argc_tmp,f_argv_tmp,&cmd_args);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_input_f"
void lis_input_f(LIS_MATRIX_F *A, LIS_VECTOR_F *b, LIS_VECTOR_F *x, char *filename, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,filename,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_input((LIS_MATRIX)LIS_V2P(A),(LIS_VECTOR)LIS_V2P(b),(LIS_VECTOR)LIS_V2P(x),buf);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_input_matrix_f"
void lis_input_matrix_f(LIS_MATRIX_F *A, char *filename, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,filename,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_input_matrix((LIS_MATRIX)LIS_V2P(A),buf);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_input_vector_f"
void lis_input_vector_f(LIS_VECTOR_F *v, char *filename, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,filename,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_input_vector((LIS_VECTOR)LIS_V2P(v),buf);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_output_f"
void lis_output_f(LIS_MATRIX_F *A, LIS_VECTOR_F *b, LIS_VECTOR_F *x, LIS_INT *mode, char *path, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,path,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_output((LIS_MATRIX)LIS_V2P(A),(LIS_VECTOR)LIS_V2P(b),(LIS_VECTOR)LIS_V2P(x),*mode,buf);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_output_matrix_f"
void lis_output_matrix_f(LIS_MATRIX_F *A, LIS_INT *mode, char *path, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,path,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_output_matrix((LIS_MATRIX)LIS_V2P(A),*mode,buf);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_output_vector_f"
void lis_output_vector_f(LIS_VECTOR_F *v, LIS_INT *format, char *filename, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,filename,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_output_vector((LIS_VECTOR)LIS_V2P(v),*format,buf);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_output_rhistory_f"
void lis_solver_output_rhistory_f(LIS_SOLVER_F *solver, char *filename, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,filename,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_solver_output_rhistory((LIS_SOLVER)LIS_V2P(solver),buf);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_output_rhistory_f"
void lis_esolver_output_rhistory_f(LIS_ESOLVER_F *esolver, char *filename, LIS_INT *ierr, LIS_INT len)
{
	char buf[1024];
	char *p;

	LIS_DEBUG_FUNC_IN;

	memset(buf,0x20,sizeof(buf));
	strncpy(buf,filename,len);
	p = &buf[len];
	if( len>0 )
	{
		while( *p==' ' ) p--;
		p++;
	}
	*p = '\0';

	*ierr = lis_esolver_output_rhistory((LIS_ESOLVER)LIS_V2P(esolver),buf);
	if( *ierr )	return;	

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif

