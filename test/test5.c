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
#include "lis.h"


#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(int argc, char* argv[])
{
	LIS_Comm comm;  
	LIS_MATRIX A;
	LIS_VECTOR x,b,u;
	LIS_SOLVER solver;
	LIS_INT	k,n,gn,ii,jj;
	LIS_INT	is,ie;
	int nprocs,my_rank;
	LIS_INT	nsol;
	LIS_INT	err,iter,iter_double,iter_quad;
	double time,itime,ptime,p_c_time,p_i_time;
	LIS_REAL resid;
	char solvername[128];
	LIS_INT	*ptr,*index;
	LIS_SCALAR *value,gamma;


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

	if( argc < 3 )
	{
	  lis_printf(comm,"Usage: %s n gamma [options]\n", argv[0]);
	  CHKERR(1);	  	    	    	  
	}

	gn  = atoi(argv[1]);
	gamma  = atof(argv[2]);
	if( gn<=0 )
	{
	  lis_printf(comm,"n=%D <=0 \n",gn);
	  CHKERR(1);
	}

	lis_printf(comm,"\n");
	lis_printf(comm,"number of processes = %d\n",nprocs);

#ifdef _OPENMP
	lis_printf(comm,"max number of threads = %d\n",omp_get_num_procs());
	lis_printf(comm,"number of threads = %d\n",omp_get_max_threads());
#endif
		
	lis_printf(comm,"n = %D, gamma = %f\n\n",gn,(double)gamma);
		
	/* create matrix and vectors */
	err = lis_matrix_create(comm,&A); CHKERR(err);
	err = lis_matrix_set_size(A,0,gn); CHKERR(err);
	err = lis_matrix_get_size(A,&n,&gn); CHKERR(err);
	err = lis_matrix_malloc_csr(n,3*n,&ptr,&index,&value); CHKERR(err);
	err = lis_matrix_get_range(A,&is,&ie); CHKERR(err);

	k = 0;
	ptr[0] = 0;
	for(ii=is;ii<ie;ii++)
	{
		if( ii>1 )    { jj = ii - 2; index[k] = jj; value[k++] = gamma;}
		if( ii<gn-1 ) { jj = ii + 1; index[k] = jj; value[k++] = 1.0;}
		index[k] = ii; value[k++] = 2.0;
		ptr[ii-is+1] = k;
	}
	err = lis_matrix_set_csr(ptr[ie-is],ptr,index,value,A); CHKERR(err);
	err = lis_matrix_assemble(A); CHKERR(err);

	err = lis_vector_duplicate(A,&u); CHKERR(err);
	err = lis_vector_duplicate(u,&b); CHKERR(err);
	err = lis_vector_duplicate(u,&x); CHKERR(err);

	err = lis_vector_set_all(1.0,u);
	lis_matvec(A,u,b);

	err = lis_solver_create(&solver); CHKERR(err);
	lis_solver_set_option("-print mem",solver);
	err = lis_solver_set_optionC(solver);
	CHKERR(err);	

	err = lis_solve(A,b,x,solver); CHKERR(err);
	lis_solver_get_iterex(solver,&iter,&iter_double,&iter_quad);
	lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
	lis_solver_get_residualnorm(solver,&resid);
	lis_solver_get_solver(solver,&nsol);
	lis_solver_get_solvername(nsol,solvername);

	lis_printf(comm,"%s: number of iterations = %D\n",solvername,iter);
#ifndef _LONG__DOUBLE
	lis_printf(comm,"%s:   double             = %D\n",solvername,iter_double);
	lis_printf(comm,"%s:   quad               = %D\n",solvername,iter_quad);		
#endif
	lis_printf(comm,"%s: elapsed time         = %e sec.\n",solvername,time);
	lis_printf(comm,"%s:   preconditioner     = %e sec.\n",solvername, ptime);
	lis_printf(comm,"%s:     matrix creation  = %e sec.\n",solvername, p_c_time);
	lis_printf(comm,"%s:   linear solver      = %e sec.\n",solvername, itime);
	lis_printf(comm,"%s: relative residual    = %e\n\n",solvername,(double)resid);
	
	lis_solver_destroy(solver);
	lis_matrix_destroy(A);
	lis_vector_destroy(b);
	lis_vector_destroy(x);
	lis_vector_destroy(u);

	lis_finalize();

	LIS_DEBUG_FUNC_OUT;

	return 0;
}


