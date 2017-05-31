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
#include <string.h>
#include "lis.h"

#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(int argc, char* argv[])
{
    LIS_Comm comm;
    LIS_INT err,nnz,is,ie;
    LIS_INT l,m,n,nn,i,j,k,ii,jj,ctr;
    int nprocs,mtype,my_rank;
    LIS_INT nesol;
    LIS_MATRIX A,A0;
    LIS_VECTOR x;
    LIS_SCALAR evalue0;
    LIS_ESOLVER esolver;
    LIS_REAL residual;
    LIS_INT iter;
    double time;
    double itime,ptime,p_c_time,p_i_time;
    LIS_INT *ptr,*index;
    LIS_SCALAR *value;
    char esolvername[128];
    
    LIS_DEBUG_FUNC_IN;
    
    lis_initialize(&argc, &argv);

    comm = LIS_COMM_WORLD;

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#else
    nprocs  = 1;
    my_rank = 0;
#endif
    
    if( argc < 7 )
      {
	lis_printf(comm,"Usage: %s l m n matrix_type evector_filename rhistory_filename [options]\n", argv[0]);
	CHKERR(1);		
      }

    lis_printf(comm,"\n");
    lis_printf(comm,"number of processes = %d\n",nprocs);

#ifdef _OPENMP
    lis_printf(comm,"max number of threads = %d\n",omp_get_num_procs());
    lis_printf(comm,"number of threads = %d\n",omp_get_max_threads());
#endif
		
    /* generate matrix */
    l  = atoi(argv[1]);
    m  = atoi(argv[2]);
    n  = atoi(argv[3]);
    mtype  = atoi(argv[4]);
    if( l<=0 || m<=0 || n<=0 )
      {
	lis_printf(comm,"l=%D <=0, m=%D <=0 or n=%D <=0\n",l,m,n);
	CHKERR(1);
      }
    nn = l*m*n;
    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    err = lis_matrix_set_size(A,0,nn);
    CHKERR(err);

    ptr   = (LIS_INT *)malloc((A->n+1)*sizeof(LIS_INT));
    if( ptr==NULL ) CHKERR(1);
    index = (LIS_INT *)malloc(7*A->n*sizeof(LIS_INT));
    if( index==NULL ) CHKERR(1);
    value = (LIS_SCALAR *)malloc(7*A->n*sizeof(LIS_SCALAR));
    if( value==NULL ) CHKERR(1);

    lis_matrix_get_range(A,&is,&ie);
    ctr = 0;
    for(ii=is;ii<ie;ii++)
      {
	i = ii/(m*n);
	j = (ii - i*m*n)/n;
	k = ii - i*m*n - j*n;
	if( i>0 )   { jj = ii - m*n; index[ctr] = jj; value[ctr++] = -1.0;}
	if( i<l-1 ) { jj = ii + m*n; index[ctr] = jj; value[ctr++] = -1.0;}
	if( j>0 )   { jj = ii - n; index[ctr] = jj; value[ctr++] = -1.0;}
	if( j<m-1 ) { jj = ii + n; index[ctr] = jj; value[ctr++] = -1.0;}
	if( k>0 )   { jj = ii - 1; index[ctr] = jj; value[ctr++] = -1.0;}
	if( k<n-1 ) { jj = ii + 1; index[ctr] = jj; value[ctr++] = -1.0;}
	index[ctr] = ii; value[ctr++] = 6.0;
	ptr[ii-is+1] = ctr;
      }
    ptr[0] = 0;
    err = lis_matrix_set_csr(ptr[ie-is],ptr,index,value,A);
    CHKERR(err);
    err = lis_matrix_assemble(A);
    CHKERR(err);

    nnz = A->nnz;
#ifdef USE_MPI
    MPI_Allreduce(&nnz,&i,1,LIS_MPI_INT,MPI_SUM,A->comm);
    nnz   = i;
#endif
    lis_printf(comm,"matrix size = %D x %D (%D nonzero entries)\n\n",nn,nn,nnz);

    err = lis_matrix_duplicate(A,&A0);
    CHKERR(err);
    lis_matrix_set_type(A0,mtype);
    err = lis_matrix_convert(A,A0);
    CHKERR(err);
    lis_matrix_destroy(A);
    A = A0;
    lis_vector_duplicate(A,&x);

    err = lis_esolver_create(&esolver);
    CHKERR(err);

    lis_esolver_set_option("-eprint mem",esolver);
    err = lis_esolver_set_optionC(esolver);
    CHKERR(err);    
    err = lis_esolve(A, x, &evalue0, esolver);
    CHKERR(err);    
    lis_esolver_get_esolver(esolver,&nesol);
    lis_esolver_get_esolvername(nesol,esolvername);
    lis_esolver_get_residualnorm(esolver, &residual);
    lis_esolver_get_iter(esolver, &iter);
    lis_esolver_get_timeex(esolver,&time,&itime,&ptime,&p_c_time,&p_i_time);

    lis_printf(comm,"%s: mode number          = %d\n", esolvername, 0);
#ifdef _COMPLEX      
    lis_printf(comm,"%s: eigenvalue           = (%e, %e)\n", esolvername, (double)creal(evalue0), (double)cimag(evalue0));
#else
    lis_printf(comm,"%s: eigenvalue           = %e\n", esolvername, (double)evalue0);
#endif
    lis_printf(comm,"%s: number of iterations = %D\n",esolvername, iter);
    lis_printf(comm,"%s: elapsed time         = %e sec.\n", esolvername, time);
    lis_printf(comm,"%s:   preconditioner     = %e sec.\n", esolvername, ptime);
    lis_printf(comm,"%s:     matrix creation  = %e sec.\n", esolvername, p_c_time);
    lis_printf(comm,"%s:   linear solver      = %e sec.\n", esolvername, itime);
    lis_printf(comm,"%s: relative residual    = %e\n\n",esolvername, (double)residual);

    /* write eigenvector */
    lis_output_vector(x,LIS_FMT_MM,argv[5]);

    /* write residual history */
    lis_esolver_output_rhistory(esolver, argv[6]);

    CHKERR(err);

    lis_esolver_destroy(esolver);
    lis_matrix_destroy(A);
    lis_matrix_destroy(A0);
    lis_vector_destroy(x);

    lis_finalize();

    LIS_DEBUG_FUNC_OUT;

    return 0;
}

 
