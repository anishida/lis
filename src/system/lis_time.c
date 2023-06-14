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
#include <time.h>

#if defined(WIN32) || defined(_WIN32)
	#include <windows.h>
#else
	#include <sys/time.h>
#endif

#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif

/************************************************
 * lis_wtime
 * lis_date
 ************************************************/

#ifdef USE_MPI
double lis_wtime(void)
{
	double etime;
	
	etime = MPI_Wtime();

	return etime;
}
#else
#ifdef _OPENMP
double lis_wtime(void)
{
	double etime;
	
	etime = omp_get_wtime();

	return etime;
}
#else
#ifdef WIN32
static LARGE_INTEGER t = {0, 0};
double lis_wtime(void)
{

  LARGE_INTEGER l;
  if (t.QuadPart == 0)
    QueryPerformanceFrequency(&t);
  QueryPerformanceCounter(&l);
  return (double)l.QuadPart / (double)t.QuadPart;
}
#else
#ifdef HAVE_CLOCK_GETTIME
double lis_wtime(void)
{
	double etime;
	struct timespec ts;

	clock_gettime(CLOCK_REALTIME, &ts);
	etime = tv.tv_sec + (double)tv.tv_usec*1.0e-9;
	return etime;
}
#else
double lis_wtime(void)
{
	double etime;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	etime = tv.tv_sec + (double)tv.tv_usec*1.0e-6;
	return etime;
}
#endif
#endif
#endif
#endif

#ifdef WIN32
void lis_date(char *date)
{
	SYSTEMTIME now;
	GetLocalTime(&now);
	sprintf(date, "%04d%02d%02d%02d%02d%02d",
        now.wYear, now.wMonth, now.wDay,
        now.wHour, now.wMinute, now.wSecond);
}
#else
void lis_date(char *date)
{
	struct timeval tv;
	struct tm *tp;

	gettimeofday(&tv,NULL);
	tp = localtime(&tv.tv_sec);
	sprintf(date, "%04d%02d%02d%02d%02d%02d",
        tp->tm_year + 1900, tp->tm_mon + 1, tp->tm_mday,
        tp->tm_hour, tp->tm_min, tp->tm_sec);
}
#endif
