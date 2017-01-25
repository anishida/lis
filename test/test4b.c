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
#include "lis.h"
LIS_INT main(int argc, char* argv[])
{
    LIS_Comm comm;  
    LIS_MATRIX A;
    LIS_VECTOR b,x,u;
    LIS_SOLVER solver;
    int nprocs,my_rank;
    LIS_INT err,l,m,n,nn,i,j,k,ii,jj,gn,is,ie,step,s,iter;
    LIS_INT *nnz;
    double time,timea,timeb,timec,timed,timee,time0,time1;    

    lis_initialize(&argc, &argv);
    comm = LIS_COMM_WORLD;

#ifdef USE_MPI
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&my_rank);
#else
    my_rank = 0;
#endif

    if( argc < 5 )
      {
	lis_printf(comm,"Usage: %s l m n step [options]\n", argv[0]);
	CHKERR(1);	  	    
      }

    l  = atoi(argv[1]);
    m  = atoi(argv[2]);
    n  = atoi(argv[3]);
    step = atoi(argv[4]);

    if( l<=0 || m<=0 || n<=0 )
      {
	lis_printf(comm,"l=%D <=0, m=%D <=0 or n=%D <=0\n",l,m,n);
	CHKERR(1);
      } 

    lis_printf(comm,"\n");
    lis_printf(comm,"number of processes = %d\n",nprocs);

#ifdef _OPENMP
    lis_printf(comm,"max number of threads = %d\n",omp_get_num_procs());
    lis_printf(comm,"number of threads = %d\n",omp_get_max_threads());
#endif
		
    /* create matrix and vectors */
    nn = l*m*n;    
    time0 = lis_wtime();    
    lis_matrix_create(comm,&A); 
    err = lis_matrix_set_size(A,0,nn); CHKERR(err);
    time1 = lis_wtime();
    /*
    err = lis_matrix_malloc(A,7*nn,nnz); CHKERR(err);    
    */
    timea = lis_wtime() - time1;
    lis_vector_duplicate(A,&u);
    lis_vector_duplicate(A,&b);
    lis_vector_duplicate(A,&x);
    lis_solver_create(&solver);
    lis_solver_set_option("-print mem",solver);
    lis_solver_set_optionC(solver);
    lis_matrix_get_size(A,&nn,&gn);
    lis_matrix_get_range(A,&is,&ie);

    timeb = 0;
    timec = 0;
    timed = 0;
    for(s=0;s<step;s++)
      {
	lis_printf(comm,"\n");
	lis_printf(comm,"step = %D\n",s);
	lis_printf(comm,"\n");

	time1 = lis_wtime();	
	for(ii=is;ii<ie;ii++)
	{
	  i = ii/(m*n);
	  j = (ii - i*m*n)/n;
	  k = ii - i*m*n - j*n;
	  if( i>0 )   lis_matrix_set_value(LIS_INS_VALUE,ii,ii-m*n,-1.0,A);
	  if( i<l-1 ) lis_matrix_set_value(LIS_INS_VALUE,ii,ii+m*n,-1.0,A);
	  if( j>0 ) lis_matrix_set_value(LIS_INS_VALUE,ii,ii-n,-1.0,A); 
	  if( j<m-1 ) lis_matrix_set_value(LIS_INS_VALUE,ii,ii+n,-1.0,A);
	  if( k>0 ) lis_matrix_set_value(LIS_INS_VALUE,ii,ii-1,-1.0,A);
	  if( k<n-1 ) lis_matrix_set_value(LIS_INS_VALUE,ii,ii+1,-1.0,A);
	  lis_matrix_set_value(LIS_INS_VALUE,ii,ii,6.0,A);
	}
	timeb += lis_wtime() - time1;
	lis_matrix_set_type(A,LIS_MATRIX_CSR);
	time1 = lis_wtime();
	lis_matrix_assemble(A);
	timec += lis_wtime() - time1;
	lis_vector_set_all(1.0,u);
	lis_matvec(A,u,b);
	time1 = lis_wtime();
	lis_solve(A,b,x,solver);
	timed += lis_wtime() - time1;
	lis_solver_get_iter(solver,&iter);

	lis_printf(comm,"number of iterations = %D\n",iter);
	lis_printf(comm,"\n");
	
	lis_matrix_unset(A);

      }
    
    time1 = lis_wtime();    
    lis_matrix_destroy(A);
    timee = lis_wtime() - time1;    
    lis_vector_destroy(b);
    lis_vector_destroy(x);
    lis_vector_destroy(u);
    lis_solver_destroy(solver);
    time = lis_wtime() - time0;    

    lis_printf(comm,"\n");
    lis_printf(comm,"elapsed time           = %e sec.\n", time);
    lis_printf(comm,"  lis_matrix_malloc    = %e sec.\n", timea);    
    lis_printf(comm,"  lis_matrix_set_value = %e sec.\n", timeb);
    lis_printf(comm,"  lis_matrix_assemble  = %e sec.\n", timec);
    lis_printf(comm,"  lis_solve            = %e sec.\n", timed);
    lis_printf(comm,"  lis_matrix_destroy   = %e sec.\n", timee);
    
    lis_finalize();
    return 0;
}
