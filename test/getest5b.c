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
#include "lis.h"

#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(int argc, char* argv[])
{
  LIS_Comm comm;
  LIS_INT err;
  int nprocs,my_rank;
  LIS_INT nesol;
  LIS_MATRIX A,B;
  LIS_VECTOR x,y;
  LIS_SCALAR evalue0;
  LIS_ESOLVER esolver;
  double time;
  char esolvername[128];

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
      lis_printf(comm,"Usage: %s matrix_a_filename matrix_b_filename evalues_filename [options]\n", argv[0]);
      CHKERR(1);
    }

  lis_printf(comm,"\n");
  lis_printf(comm,"number of processes = %d\n",nprocs);

#ifdef _OPENMP
  lis_printf(comm,"max number of threads = %d\n",omp_get_num_procs());
  lis_printf(comm,"number of threads = %d\n",omp_get_max_threads());
#endif
		
  /* create matrix and vectors */

  lis_matrix_create(comm,&A);
  lis_matrix_create(comm,&B);
  lis_printf(comm,"\nmatrix A:\n");
  lis_input_matrix(A,argv[1]);
  lis_printf(comm,"matrix B:\n");
  lis_input_matrix(B,argv[2]);
  lis_vector_duplicate(A,&x);
  lis_esolver_create(&esolver);
  lis_esolver_set_option("-e gli -ss 1 -eprint mem -rval true",esolver);
  err = lis_esolver_set_optionC(esolver);
  CHKERR(err);  
  err = lis_gesolve(A,B,x,&evalue0,esolver);
  CHKERR(err);  
  lis_esolver_get_esolver(esolver,&nesol);
  lis_esolver_get_esolvername(nesol,esolvername);

  lis_vector_create(comm,&y);
  lis_esolver_get_evalues(esolver,y);

  /* write eigenvalues */
  lis_output_vector(y,LIS_FMT_MM,argv[3]);

  lis_esolver_destroy(esolver);
  lis_matrix_destroy(A);
  lis_matrix_destroy(B);
  lis_vector_destroy(x);
  lis_vector_destroy(y);

  lis_finalize();

  LIS_DEBUG_FUNC_OUT;

  return 0;
}

 
