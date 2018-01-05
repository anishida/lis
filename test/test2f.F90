!   Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!   1. Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!   3. Neither the name of the project nor the names of its contributors 
!      may be used to endorse or promote products derived from this software 
!      without specific prior written permission.
!
!   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
!   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
!   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!   POSSIBILITY OF SUCH DAMAGE.

      implicit none
      
#include "lisf.h"

      LIS_MATRIX :: A,A0
      LIS_VECTOR :: b,x,u
      LIS_SOLVER :: solver
      LIS_INTEGER :: ia,ierr
      integer*4 :: nprocs,my_rank
      LIS_INTEGER :: matrix_type,comm
      LIS_INTEGER :: omp_get_num_procs,omp_get_max_threads
      LIS_INTEGER :: m,n,nn,nnz,innz
      LIS_INTEGER :: i,j,ii,jj,ctr
      LIS_INTEGER :: is,ie
      LIS_INTEGER,allocatable :: ptr(:),index(:)
      LIS_SCALAR,allocatable :: value(:)
      LIS_SCALAR :: one=1.0d0
      LIS_INTEGER :: nsol,iter,iter_double,iter_quad
      real*8 :: time,itime,ptime,p_c_time,p_i_time
      LIS_REAL :: resid
      character*256 :: solname,resname
      character*20 :: m_char,n_char,mt_char,solvername
      integer*4 :: iargc
      
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

      ia = iargc()
      if( ia.lt.5 ) then
         if( my_rank.eq.0 ) then
          write(*,'(a)') 'Usage: test2f m n matrix_type solution_filename residual_filename [options]'
          call lis_finalize(ierr)
        endif
        stop
      endif

      if ( my_rank .eq. 0 ) then
         write(*,'(a)') ''
         write(*,'(a,i0)') 'number of processes = ',nprocs
#ifdef _OPENMP
         write(*,'(a,i0)') 'max number of threads = ',omp_get_num_procs()
         write(*,'(a,i0)') 'number of threads = ', omp_get_max_threads()
#endif
      endif

! read matrix and vectors from file 
      call getarg(1,m_char)
      read(m_char,*) m
      call getarg(2,n_char)
      read(n_char,*) n
      call getarg(3,mt_char)
      read(mt_char,*) matrix_type

      nn = m * n
      call lis_matrix_create(comm,A,ierr)
      call CHKERR(ierr)
      call lis_matrix_set_size(A,0,nn,ierr)
      call CHKERR(ierr)

      allocate(ptr(0:nn))
      allocate(index(0:5*nn-1))
      allocate(value(0:5*nn-1))

      call lis_matrix_get_range(A,is,ie,ierr)
      ctr = 0
#ifdef COMPLEX      
      do ii=is-1,ie-2
         i = ii/m
         j = ii - i*m
         if ( i .GT. 0 ) then
            jj = ii - m
            index(ctr) = jj
            value(ctr) = (-1.0d0,0.0d0)
            ctr = ctr + 1
         end if
         if ( i .LT. n-1 ) then
            jj = ii + m
            index(ctr) = jj
            value(ctr) = (-1.0d0,0.0d0)
            ctr = ctr + 1
         end if
         if ( j .GT. 0 ) then
            jj = ii - 1
            index(ctr) = jj
            value(ctr) = (-1.0d0,0.0d0)
            ctr = ctr + 1
         end if
         if ( j .LT. m-1 ) then
            jj = ii + 1
            index(ctr) = jj
            value(ctr) = (-1.0d0,0.0d0)
            ctr = ctr + 1
         end if
         index(ctr) = ii
         value(ctr) = (4.0d0,0.0d0)
         ctr = ctr + 1
         ptr(ii-(is-1)+1) = ctr
      end do
#else
      do ii=is-1,ie-2
         i = ii/m
         j = ii - i*m
         if ( i .GT. 0 ) then
            jj = ii - m
            index(ctr) = jj
            value(ctr) = -1.0d0
            ctr = ctr + 1
         end if
         if ( i .LT. n-1 ) then
            jj = ii + m
            index(ctr) = jj
            value(ctr) = -1.0d0
            ctr = ctr + 1
         end if
         if ( j .GT. 0 ) then
            jj = ii - 1
            index(ctr) = jj
            value(ctr) = -1.0d0
            ctr = ctr + 1
         end if
         if ( j .LT. m-1 ) then
            jj = ii + 1
            index(ctr) = jj
            value(ctr) = -1.0d0
            ctr = ctr + 1
         end if
         index(ctr) = ii
         value(ctr) = 4.0d0
         ctr = ctr + 1
         ptr(ii-(is-1)+1) = ctr
      end do
#endif      
      ptr(0) = 0
      call lis_matrix_set_csr(ptr(ie-is),ptr,index,value,A,ierr)
      call CHKERR(ierr)
      call lis_matrix_assemble(A,ierr)
      call CHKERR(ierr)
      call lis_matrix_get_nnz(A,nnz,ierr)

#ifdef USE_MPI
      call MPI_Allreduce(nnz,innz,1,LIS_MPI_INTEGER,MPI_SUM,comm,ierr)
      nnz = innz
#endif

      if ( my_rank .EQ. 0 ) then
         write(*,'(a,i0,a,i0,a,i0,a)') 'matrix size = ', nn, ' x ', nn, ' (', nnz, ' nonzero entries)'
         write(*,'(a)')
      endif

      call lis_matrix_duplicate(A,A0,ierr)
      call CHKERR(ierr)
      call lis_matrix_set_type(A0,matrix_type,ierr)
      call lis_matrix_convert(A,A0,ierr)
      call CHKERR(ierr)
      call lis_matrix_destroy(A,ierr)
      A = A0
      
!      call lis_vector_create(comm,u,ierr)      
!      call lis_vector_set_size(u,0,nn,ierr)        

      call lis_vector_duplicate(A,u,ierr)
      call CHKERR(ierr)
      call lis_vector_duplicate(A,b,ierr)
      call CHKERR(ierr)
      call lis_vector_duplicate(A,x,ierr)
      call CHKERR(ierr)

      call lis_vector_set_all(one,u,ierr)      
      call lis_matvec(A,u,b,ierr)

      call lis_solver_create(solver,ierr)
      call CHKERR(ierr)
      call lis_solver_set_option('-print mem',solver,ierr)
      call lis_solver_set_optionC(solver,ierr)
      call CHKERR(ierr)      
      
      call lis_solve(A,b,x,solver,ierr)

      call CHKERR(ierr)

      call lis_solver_get_iterex(solver,iter,iter_double,iter_quad,ierr)
      call lis_solver_get_timeex(solver,time,itime,ptime,p_c_time,p_i_time,ierr)
      call lis_solver_get_residualnorm(solver,resid,ierr)
      call lis_solver_get_solver(solver,nsol,ierr)
      call lis_solver_get_solvername(nsol,solvername,ierr)

      if( my_rank.eq.0 ) then
        write(*,'(a,a,i0)') solvername,': number of iterations = ', iter
#ifndef LONG__DOUBLE
        write(*,'(a,a,i0)') solvername,':   double             = ', iter_double
        write(*,'(a,a,i0)') solvername,':   quad               = ', iter_quad
#endif
        write(*,'(a,a,e14.7e2,a)') solvername,': elapsed time         = ', time, ' sec.'
        write(*,'(a,a,e14.7e2,a)') solvername,':   preconditioner     = ', ptime, ' sec.'
        write(*,'(a,a,e14.7e2,a)') solvername,':     matrix creation  = ', p_c_time, ' sec.'
        write(*,'(a,a,e14.7e2,a)') solvername,':   linear solver      = ', itime, ' sec.'
        write(*,'(a,a,e14.7e2)') solvername,': relative residual    = ', resid
        write(*,'(a)') ''        
      endif

! write solution 
      call getarg(4,solname)
      call lis_output_vector(x,LIS_FMT_MM,solname,ierr);

! write residual history
      call getarg(5,resname)
      call lis_solver_output_rhistory(solver, resname, ierr)

      call lis_solver_destroy(solver,ierr)
      call lis_matrix_destroy(A,ierr)
      call lis_vector_destroy(u,ierr)
      call lis_vector_destroy(x,ierr)
      call lis_vector_destroy(b,ierr)

      call lis_finalize(ierr)

      stop
      end
      
