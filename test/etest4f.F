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
 
      LIS_INTEGER i,n,gn,ln,nnz,is,ie,iter,ierr
      LIS_MATRIX A
      LIS_VECTOR x
      LIS_SCALAR evalue0
      LIS_ESOLVER esolver
      integer*4 my_rank,nprocs,comm
      LIS_INTEGER matrix_type
      LIS_INTEGER omp_get_num_procs,omp_get_max_threads      
      LIS_INTEGER nsol,iter_double,iter_quad
      real*8 time,itime,ptime,p_c_time,p_i_time
      LIS_REAL resid
      character*20 nchar,esolvername
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

      call lis_matrix_create(comm,A,ierr)
      i = iargc()

      if( i.lt.1 ) then
        if( my_rank.eq.0 ) then
          write(*,'(a)') 'etest4f n [options]'
          call lis_finalize(ierr)
        endif
        stop
      endif

      if (my_rank .eq. 0) then
         write(*,'(a)') ''
         write(*,'(a,i0)') 'number of processes = ',nprocs
      endif

#ifdef _OPENMP
      write(*,'(a,i0)') 'max number of threads = ',omp_get_num_procs()
      write(*,'(a,i0)') 'number of threads = ', omp_get_max_threads()
#endif

      call getarg(1,nchar)
      read (nchar,*) n

      ln = 0

      call lis_matrix_create(comm,A,ierr)
      call lis_matrix_set_size(A,ln,n,ierr)
      call lis_matrix_get_size(A,n,gn,ierr)
      call lis_matrix_get_range(A,is,ie,ierr)
#ifdef COMPLEX
#ifdef LONG__DOUBLE      
      do i=is,ie-1
        if( i>1  ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
     .                                        (-1.0q0,0.0q0),A,ierr)
        if( i<gn ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
     .                                        (-1.0q0,0.0q0),A,ierr)
        call lis_matrix_set_value(LIS_INS_VALUE,i,i,(2.0q0,0.0q0),
     .                                        A,ierr)
      enddo
#else
      do i=is,ie-1
        if( i>1  ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
     .                                        (-1.0d0,0.0d0),A,ierr)
        if( i<gn ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
     .                                        (-1.0d0,0.0d0),A,ierr)
        call lis_matrix_set_value(LIS_INS_VALUE,i,i,(2.0d0,0.0d0),
     .                                        A,ierr)
      enddo
#endif      
#else
#ifdef LONG__DOUBLE      
      do i=is,ie-1
        if( i>1  ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0q0,
     .                                        A,ierr)
        if( i<gn ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0q0,
     .                                        A,ierr)
        call lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0q0,A,ierr)
      enddo
#else
      do i=is,ie-1
        if( i>1  ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0d0,
     .                                        A,ierr)
        if( i<gn ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0d0,
     .                                        A,ierr)
        call lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0d0,A,ierr)
      enddo
#endif      
#endif      
      call lis_matrix_set_type(A,matrix_type,ierr)
      call lis_matrix_assemble(A,ierr)
      call lis_matrix_get_nnz(A,nnz,ierr)      
      
      if ( my_rank .EQ. 0 ) then
         write(*,'(a,i0,a,i0,a,i0,a)') 'matrix size = ', n, ' x ', n,
     .        ' (', nnz, ' nonzero entries)'
         write(*,'(a)')
      endif

      call lis_vector_duplicate(A,x,ierr)
#ifdef COMPLEX
#ifdef LONG__DOUBLE      
      call lis_vector_set_all((1.0q0,0.0q0),x,ierr)
#else
      call lis_vector_set_all((1.0d0,0.0d0),x,ierr)
#endif      
#else
#ifdef LONG__DOUBLE      
      call lis_vector_set_all(1.0q0,x,ierr)
#else
      call lis_vector_set_all(1.0d0,x,ierr)
#endif      
#endif      
      
      call lis_esolver_create(esolver,ierr)
      call lis_esolver_set_option('-eprint mem',esolver,ierr)
      call lis_esolver_set_optionC(esolver,ierr)
      call CHKERR(ierr)      
      call lis_esolve(A,x,evalue0,esolver,ierr)
      call CHKERR(ierr)
      call lis_esolver_get_iterex(esolver,iter,iter_double,iter_quad,
     .     ierr)
      call lis_esolver_get_timeex(esolver,time,itime,ptime,
     .     p_c_time,p_i_time,ierr)
      call lis_esolver_get_residualnorm(esolver,resid,ierr)
      call lis_esolver_get_esolver(esolver,nsol,ierr)
      call lis_esolver_get_esolvername(nsol,esolvername,ierr)

      if( my_rank.eq.0 ) then
        write(*,'(a,a,i0)') esolvername,': mode number          = ',0
#ifdef COMPLEX         
        write(*,'(a,a,a,e14.7e2,a,e14.7e2,a)') esolvername,
     .       ': eigenvalue           = ','(',real(evalue0),', ',
     .       imag(evalue0),')'
#else
        write(*,'(a,a,e14.7e2)') esolvername,': eigenvalue           = '
     .       ,evalue0
#endif        
        write(*,'(a,a,i0)') esolvername,': number of iterations = ',iter
        write(*,'(a,a,e14.7e2,a)') esolvername,
     .       ': elapsed time         = ',time,' sec.'
        write(*,'(a,a,e14.7e2,a)') esolvername,
     .       ':   preconditioner     = ',ptime,' sec.'
        write(*,'(a,a,e14.7e2,a)') esolvername,
     .       ':     matrix creation  = ',p_c_time,' sec.'
        write(*,'(a,a,e14.7e2,a)') esolvername,
     .       ':   linear solver      = ',itime,' sec.'
        write(*,'(a,a,e14.7e2)') esolvername,': relative residual    = '
     .       ,resid
        write(*,'(a)') ''
      endif

      call lis_matrix_destroy(A,ierr)
      call lis_vector_destroy(x,ierr)
      call lis_esolver_destroy(esolver,ierr)
      call lis_finalize(ierr)

      stop
      end
      
