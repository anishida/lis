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
#include <string.h>
#include <stdarg.h>
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_debug_trace_func
 * lis_printf
 * lis_error
 * lis_print_rhistory
 ************************************************/

char *LIS_ERR_MESS[] = {
	"ILL_ARG",
	"BREAKDOWN",
	"OUT_OF_MEMORY",
	"MAXITER",
	"NOT_IMPLEMENTED",
	"FILE_IO_ERROR"
};

LIS_Comm lis_debug_comm = LIS_COMM_WORLD;

void lis_debug_set_comm(LIS_Comm comm)
{
	lis_debug_comm = LIS_COMM_WORLD;
}

LIS_INT lis_debug_trace_func(LIS_INT flag, char *func)
{
	static int lis_tr_func = 0;
	char buf[1024];
	int argc=0;

	#ifdef USE_MPI
		MPI_Initialized(&lis_mpi_initialized);
		if (!lis_mpi_initialized) MPI_Init(&argc, NULL);
	#endif
	if( flag )
	{
		sprintf(buf, "%%%ds : %%s\n",lis_tr_func+3);
		lis_printf(lis_debug_comm, buf, "IN ", func);
		lis_tr_func++;
	}
	else
	{
		lis_tr_func--;
		sprintf(buf, "%%%ds : %%s\n",lis_tr_func+3);
		lis_printf(lis_debug_comm, buf, "OUT", func);
	}
	#ifdef USE_MPI
		if( strcmp(func,"main")==0 && !flag )
		{
			MPI_Finalize();
		}
	#endif
	return LIS_SUCCESS;
}

LIS_INT lis_replace(char *buf, const char *str1, const char *str2)
{
  char tmp[1024+1];
  char *p;

  while ((p = strstr(buf,str1)) != NULL)
    {
      *p = '\0';
      p += strlen(str1);
      strcpy(tmp,p);
      strcat(buf,str2);
      strcat(buf,tmp);
    }
  return LIS_SUCCESS;
}

LIS_INT lis_vprintf(const char *mess, va_list vvlist)
{
	char str[1024];

	strcpy(str,mess);
	#ifdef _LONG__LONG
	lis_replace(str,"%D","%lld");
	#else
	lis_replace(str,"%D","%d");
	#endif
	vprintf(str,vvlist);
	return LIS_SUCCESS;
}

LIS_INT lis_printf(LIS_Comm comm, const char *mess, ...)
{
	va_list vvlist;
	int my_rank;
	char str[1024];

	#ifdef USE_MPI
        #ifdef _DEBUG
	MPI_Barrier(comm);
	#endif
	MPI_Comm_rank(comm,&my_rank);
	#else
	my_rank = 0;
	#endif

	if( my_rank==0 )
	{
		va_start(vvlist,mess);
		lis_vprintf(mess,vvlist);
		va_end(vvlist);
	}
	#ifdef USE_MPI
        #ifdef _DEBUG
	MPI_Barrier(comm);
	#endif
	#endif
	return LIS_SUCCESS;
}

LIS_INT lis_error(const char *file, const char *func, const LIS_INT line, const LIS_INT code, const char *mess, ...)
{
	va_list vvlist;
	int my_rank;
	char str[1024];

	#ifdef USE_MPI
	MPI_Comm_rank(lis_debug_comm,&my_rank);
	#else
	my_rank = 0;
	#endif

	if( my_rank==0 )
	{
		va_start(vvlist, mess);
		lis_printf(lis_debug_comm,"%s(%D) : %s : error %s :",file,line,func,LIS_ERR_MESS[code-LIS_ERR_ILL_ARG]);
		if( mess ) lis_vprintf(mess,vvlist);
		va_end(vvlist);
	}
	return LIS_SUCCESS;
}

void CHKERR(LIS_INT ierr)
{
	if(ierr)
	{
		lis_finalize();
		exit((int)ierr);
	}
}


LIS_INT lis_print_rhistory(LIS_Comm comm, LIS_INT iter, LIS_REAL resid)
{

#ifdef _LONG__LONG	
	lis_printf(comm,"iteration: %5lld  relative residual = %e\n", iter, (double)resid);
#else
	lis_printf(comm,"iteration: %5d  relative residual = %E\n", iter, (double)resid);
#endif

	return LIS_SUCCESS;
}



