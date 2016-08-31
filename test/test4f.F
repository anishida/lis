C   Copyright (C) 2005 The Scalable Software Infrastructure Project. 
C   All rights reserved.
C
C   Redistribution and use in source and binary forms, with or without
C   modification, are permitted provided that the following conditions
C   are met:
C   1. Redistributions of source code must retain the above copyright
C      notice, this list of conditions and the following disclaimer.
C   2. Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C   3. Neither the name of the project nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C
C   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE
C   PROJECT ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
C   BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
C   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
C   THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT BE LIABLE FOR ANY
C   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
C   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
C   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
C   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
C   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
C   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

      implicit none
      
#include "lisf.h"

      LIS_INTEGER my_rank,nprocs
      LIS_INTEGER matrix_type,comm_world
      LIS_INTEGER omp_get_num_procs,omp_get_max_threads
      LIS_INTEGER i,n,gn,ln,is,ie,iter,ierr
      LIS_MATRIX A
      LIS_VECTOR b,x,u
      LIS_SOLVER solver
      
      call lis_initialize(ierr)

#ifdef USE_MPI
      call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
#else
      nprocs  = 1
      my_rank = 0
#endif

      n = 12
      ln = 0
      matrix_type = LIS_MATRIX_CSR
      comm_world = LIS_COMM_WORLD

      if (my_rank .eq. 0) then
         write(*,*) ''
         write(*,*) 'number of processes = ',nprocs
#ifdef _OPENMP
         write(*,*) 'max number of threads = ',omp_get_num_procs()
         write(*,*) 'number of threads = ', omp_get_max_threads()
#endif
      endif

      call lis_matrix_create(comm_world,A,ierr)
      call lis_matrix_set_size(A,ln,n,ierr)
      call lis_matrix_get_size(A,n,gn,ierr)
      call lis_matrix_get_range(A,is,ie,ierr)
      do i=is,ie-1
        if( i>1  ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0d0,
     .                                        A,ierr)
        if( i<gn ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0d0,
     .                                        A,ierr)
        call lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0d0,A,ierr)   
      enddo
      call lis_matrix_set_type(A,matrix_type,ierr)
      call lis_matrix_assemble(A,ierr)
      call lis_vector_duplicate(A,u,ierr)
      call lis_vector_duplicate(A,b,ierr)
      call lis_vector_duplicate(A,x,ierr)
      call lis_vector_set_all(1.0d0,u,ierr)
      call lis_matvec(A,u,b,ierr)
      
      call lis_solver_create(solver,ierr)
      call lis_solver_set_option("-print mem",solver,ierr);
      call lis_solver_set_optionC(solver,ierr)
      call lis_solve(A,b,x,solver,ierr)
      call lis_solver_get_iter(solver,iter,ierr)
!      write(*,*) 'number of iterations = ', iter
      call lis_vector_print(x,ierr)
      call lis_matrix_destroy(A,ierr)
      call lis_vector_destroy(b,ierr)
      call lis_vector_destroy(x,ierr)
      call lis_vector_destroy(u,ierr)
      call lis_solver_destroy(solver,ierr)
      call lis_finalize(ierr)

      stop
      end
      
