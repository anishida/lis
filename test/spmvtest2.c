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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

char *lis_storagename2[]   = {"CSR", "CSC", "MSR", "DIA", "ELL", "JAD", "BSR", "BSC", "VBR", "COO", "DNS"};

#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(int argc, char* argv[])
{
  LIS_Comm comm;
  LIS_MATRIX A,A0;
  LIS_VECTOR b,x;
  LIS_SCALAR *value;
  int nprocs,my_rank;
  LIS_INT nthreads,maxthreads;
  LIS_INT gn,nnz;
  LIS_INT i,ii,j,jj,ctr,np;
  LIS_INT m,n,nn;
  LIS_INT is,ie;
  LIS_INT err,iter,matrix_type,storage,ss,se;
  LIS_INT *ptr,*index;
  double time,time2,nnzs,nnzap,nnzt;
  LIS_REAL val;
  double commtime,comptime,flops;

  LIS_DEBUG_FUNC_IN;

  lis_initialize(&argc, &argv);

  comm = LIS_COMM_WORLD;
    
#ifdef USE_MPI
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&my_rank);
#else
  nprocs  = 1;
  my_rank = 0;
#endif

  if( argc < 4 )
    {
      lis_printf(comm,"Usage: %s m n iter [matrix_type]\n", argv[0]);
      CHKERR(1);      
    }

  m  = atoi(argv[1]);
  n  = atoi(argv[2]);
  iter = atoi(argv[3]);
  if (argv[4] == NULL) {
    storage = 0;
  }
  else {
    storage = atoi(argv[4]);
  }

  if( iter<=0 )
    {
      lis_printf(comm,"iter=%D <= 0\n",iter);
      CHKERR(1);
    }
  if( n<=0 || m<=0 )
    {
      lis_printf(comm,"m=%D <=0 or n=%D <=0\n",m,n);
      CHKERR(1);
    }
  if( storage<0 || storage>11 )
    {
      lis_printf(comm,"matrix_type=%D < 0 or matrix_type=%D > 11\n",storage,storage);
      CHKERR(1);
    }

  lis_printf(comm,"\n");
  lis_printf(comm,"number of processes = %d\n",nprocs);

#ifdef _OPENMP
  nthreads = omp_get_num_procs();
  maxthreads = omp_get_max_threads();
  lis_printf(comm,"max number of threads = %d\n", nthreads);
  lis_printf(comm,"number of threads = %d\n", maxthreads);
#else
  nthreads = 1;
  maxthreads = 1;
#endif

  /* create matrix and vectors */
  nn = m*n;
  err = lis_matrix_create(comm,&A0);
  err = lis_matrix_set_size(A0,0,nn);
  CHKERR(err);

  ptr   = (LIS_INT *)malloc((A0->n+1)*sizeof(LIS_INT));
  if( ptr==NULL ) CHKERR(1);
  index = (LIS_INT *)malloc(5*A0->n*sizeof(LIS_INT));
  if( index==NULL ) CHKERR(1);
  value = (LIS_SCALAR *)malloc(5*A0->n*sizeof(LIS_SCALAR));
  if( value==NULL ) CHKERR(1);

  lis_matrix_get_range(A0,&is,&ie);
  ctr = 0;
  for(ii=is;ii<ie;ii++)
    {
      i = ii/n;
      j = ii - i*n;
      if( i>0 )   { jj = ii - n; index[ctr] = jj; value[ctr++] = -1.0;}
      if( i<m-1 ) { jj = ii + n; index[ctr] = jj; value[ctr++] = -1.0;}
      if( j>0 )   { jj = ii - 1; index[ctr] = jj; value[ctr++] = -1.0;}
      if( j<n-1 ) { jj = ii + 1; index[ctr] = jj; value[ctr++] = -1.0;}
      index[ctr] = ii; value[ctr++] = 4.0;
      ptr[ii-is+1] = ctr;
    }
  ptr[0] = 0;
  err = lis_matrix_set_csr(ptr[ie-is],ptr,index,value,A0);
  CHKERR(err);
  err = lis_matrix_assemble(A0);
  CHKERR(err);

  n   = A0->n;
  gn  = A0->gn;
  nnz = A0->nnz;
  np  = A0->np-n;

#ifdef USE_MPI
  MPI_Allreduce(&nnz,&i,1,LIS_MPI_INT,MPI_SUM,A0->comm);
  nnzap = (double)i / (double)nprocs;
  nnzt  = ((double)nnz -nnzap)*((double)nnz -nnzap);
  nnz   = i;
  MPI_Allreduce(&nnzt,&nnzs,1,MPI_DOUBLE,MPI_SUM,A0->comm);
  nnzs  = (nnzs / (double)nprocs)/nnzap;
  MPI_Allreduce(&np,&i,1,LIS_MPI_INT,MPI_SUM,A0->comm);
  np = i;
#endif

  lis_printf(comm,"matrix size = %D x %D (%D nonzero entries)\n",gn,gn,nnz);
  lis_printf(comm,"number of iterations = %D\n\n",iter);

  err = lis_vector_duplicate(A0,&x);
  if( err ) CHKERR(err);
  err = lis_vector_duplicate(A0,&b);
  if( err ) CHKERR(err);

  lis_matrix_get_range(A0,&is,&ie);
  for(i=0;i<n;i++)
    {
      err = lis_vector_set_value(LIS_INS_VALUE,i+is,1.0,x);
    }
  for(i=0;i<n;i++)
    {
      lis_sort_id(A0->ptr[i],A0->ptr[i+1]-1,A0->index,A0->value);
    }
		
  /* 
     MPI version of VBR is not implemented.
     DNS is also excluded to reduce memory usage.
  */

  if (storage==0) 
    {
      ss = 1;
      se = 11;
    }
  else
    {
      ss = storage;
      se = storage+1;
    }
	
  for (matrix_type=ss;matrix_type<se;matrix_type++)
    {
      if ( nprocs>1 && matrix_type==9 ) continue;
      lis_matrix_duplicate(A0,&A); 
      lis_matrix_set_type(A,matrix_type);
      err = lis_matrix_convert(A0,A);
      if( err ) CHKERR(err);
		    
      comptime = 0.0;
      commtime = 0.0;

      for(i=0;i<iter;i++)
	{
#ifdef USE_MPI
	  MPI_Barrier(A->comm);
	  time = lis_wtime();
	  lis_send_recv(A->commtable,x->value);
	  commtime += lis_wtime() - time;
#endif
	  time2 = lis_wtime();
	  lis_matvec(A,x,b);
	  comptime += lis_wtime() - time2;
	}
      lis_vector_nrm2(b,&val);

      flops = 2.0*nnz*iter*1.0e-6 / comptime;
#ifdef USE_MPI
      lis_printf(comm,"matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f %%, 2-norm = %e\n",(int)matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100,(double)val);
#else
      printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, 2-norm = %e\n",(int)matrix_type,lis_storagename2[matrix_type-1],comptime,flops,(double)val);
#endif
      lis_matrix_destroy(A);
    }
  
  lis_matrix_destroy(A0);
  lis_vector_destroy(b);
  lis_vector_destroy(x);

  lis_finalize();

  LIS_DEBUG_FUNC_OUT;

  return 0;
}


