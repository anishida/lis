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

      LIS_MATRIX A,A0
      LIS_VECTOR b,x,u
      LIS_SOLVER solver
      LIS_INTEGER ierr
      integer*4 my_rank,nprocs
      LIS_INTEGER matrix_type,comm
      LIS_INTEGER omp_get_num_procs,omp_get_max_threads
      LIS_INTEGER n,gn,i
      LIS_INTEGER nsol,iter,rhs,iter_double,iter_quad
      real*8 time,itime,ptime,p_c_time,p_i_time
      LIS_REAL resid
      character*256 fname,solname,resname,argc
      character*20 solvername
      integer*4 iargc

      call lis_initialize(ierr)

      comm = LIS_COMM_WORLD

#ifdef USE_MPI
      call MPI_Comm_size(comm,nprocs,ierr)
      call MPI_Comm_rank(comm,my_rank,ierr)
#else
      nprocs  = 1
      my_rank = 0
#endif

      matrix_type = LIS_MATRIX_CSR

      i = iargc()
      if( i.lt.4 ) then
         if( my_rank.eq.0 ) then
            write(*,'(a)') 'Usage: test1f matrix_filename rhs_setting ',
     .              'solution_filename rhistory_filename [options]'
            call lis_finalize(ierr)
         endif
         stop
      endif
      call getarg(2,argc)
      if( argc.eq.'0' ) then
        rhs = 0
      elseif( argc.eq.'1' ) then
        rhs = 1
      elseif( argc.eq.'2' ) then
        rhs = 2
      else
        rhs = -1
      endif
      
      if (my_rank .eq. 0) then
         write(*,'(a)') ''
         write(*,'(a,i0)') 'number of processes = ',nprocs
#ifdef _OPENMP
         write(*,'(a,i0)') 'max number of threads = ',
     .        omp_get_num_procs()
         write(*,'(a,i0)') 'number of threads = ', omp_get_max_threads()
#endif
      endif

! read matrix and vectors from file 
      call getarg(1,fname)
      call lis_matrix_create(comm,A,ierr)
      call CHKERR(ierr)
      call lis_vector_create(comm,b,ierr)
      call CHKERR(ierr)
      call lis_vector_create(comm,x,ierr)
      call CHKERR(ierr)
      call lis_matrix_set_type(A,matrix_type,ierr)
      call lis_input(A,b,x,fname,ierr)
      call CHKERR(ierr);

      call lis_vector_duplicate(A,u,ierr)
      
      call lis_matrix_get_size(A,n,gn,ierr)

      call CHKERR(ierr)
      call lis_vector_is_null(b,ierr)
      if( ierr.eq.LIS_TRUE ) then
        call lis_vector_destroy(b,ierr)
        call lis_vector_duplicate(A,b,ierr)
        call CHKERR(ierr)
        if( rhs.eq.0 ) then
          call lis_finalize(ierr)
          stop
        elseif( rhs.eq.1 ) then
#ifdef COMPLEX
#ifdef LONG__DOUBLE      
           call lis_vector_set_all((1.0q0,0.0q0),b,ierr)
#else
           call lis_vector_set_all((1.0d0,0.0d0),b,ierr)
#endif      
#else
#ifdef LONG__DOUBLE      
           call lis_vector_set_all(1.0q0,b,ierr)
#else
           call lis_vector_set_all(1.0d0,b,ierr)
#endif      
#endif      
        else
#ifdef COMPLEX
#ifdef LONG__DOUBLE      
           call lis_vector_set_all((1.0q0,0.0q0),u,ierr)
#else
           call lis_vector_set_all((1.0d0,0.0d0),u,ierr)
#endif      
#else
#ifdef LONG__DOUBLE      
           call lis_vector_set_all(1.0q0,u,ierr)
#else
           call lis_vector_set_all(1.0d0,u,ierr)
#endif      
#endif      
           call lis_matvec(A,u,b,ierr)
        endif
      endif
      if( rhs.eq.-1 ) then
        call getarg(2,fname)
        call lis_input_vector(b,fname,ierr)
        call CHKERR(ierr)
      endif

      call lis_vector_is_null(x,ierr)
      if( ierr.eq.LIS_TRUE ) then
        call lis_vector_destroy(x,ierr)
        call lis_vector_duplicate(u,x,ierr)
        call CHKERR(ierr)
      endif

      call lis_solver_create(solver,ierr)
      call CHKERR(ierr)
      call lis_solver_set_option('-print mem',solver,ierr)
      call lis_solver_set_optionC(solver,ierr)
      call CHKERR(ierr)      

      call lis_solve(A,b,x,solver,ierr)

      call CHKERR(ierr)

	call lis_solver_get_iterex(solver,iter,iter_double,iter_quad,
     .                           ierr)
	call lis_solver_get_timeex(solver,time,itime,ptime,
     .                           p_c_time,p_i_time,ierr)
	call lis_solver_get_residualnorm(solver,resid,ierr)
      call lis_solver_get_solver(solver,nsol,ierr)
      call lis_solver_get_solvername(nsol,solvername,ierr)

      if( my_rank.eq.0 ) then
        write(*,'(a,a,i0)') solvername,': number of iterations = ',iter
#ifndef LONG__DOUBLE
        write(*,'(a,a,i0)') solvername,':   double             = ',
     .       iter_double
        write(*,'(a,a,i0)') solvername,':   quad               = ',
     .       iter_quad
#endif
        write(*,'(a,a,e14.7e2,a)') solvername,
     .       ': elapsed time         = ',time,' sec.'
        write(*,'(a,a,e14.7e2,a)') solvername,
     .       ':   preconditioner     = ',ptime,' sec.'
        write(*,'(a,a,e14.7e2,a)') solvername,
     .       ':     matrix creation  = ',p_c_time,' sec.'
        write(*,'(a,a,e14.7e2,a)') solvername,
     .       ':   linear solver      = ',itime,' sec.'
        write(*,'(a,a,e14.7e2)') solvername,': relative residual    = ',
     .       resid
        write(*,'(a)') ''
      endif

! write solution 
      call getarg(3,solname)
      call lis_output_vector(x,LIS_FMT_MM,solname,ierr);

! write residual 
      call getarg(4,resname)
      call lis_solver_output_rhistory(solver, resname, ierr)

      call lis_solver_destroy(solver,ierr)
      call lis_vector_destroy(u,ierr)
      call lis_vector_destroy(x,ierr)
      call lis_vector_destroy(b,ierr)
      call lis_matrix_destroy(A,ierr)

      call lis_finalize(ierr)

      stop
      end
      
