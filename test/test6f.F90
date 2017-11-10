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

      LIS_Comm comm
      LIS_SCALAR,allocatable :: a(:),x(:),b(:),u(:),w(:)
      LIS_INTEGER :: m,n,nn
      LIS_INTEGER :: i,j,ii,jj,nnz
      LIS_INTEGER :: ia,ierr
      real*8 :: time,time0,lis_wtime
      LIS_REAL :: resid_r,resid_b
      character*20 :: m_char,n_char
      integer*4 :: iargc

      call lis_initialize(ierr) 

      comm = LIS_COMM_WORLD

      ia = iargc()
      if( ia.lt.2 ) then
         write(*,'(a)') 'Usage: test6f m n'
         call lis_finalize(ierr)
         stop
      end if

! read m and n

      call getarg(1,m_char)
      read(m_char,*) m
      call getarg(2,n_char)
      read(n_char,*) n

      nn = m * n

      allocate(a(0:nn*nn-1))
      allocate(b(0:nn-1))
      allocate(x(0:nn-1))
      allocate(u(0:nn-1))
      allocate(w(0:nn*nn-1))

! define two-dimensional Laplacian

#ifdef COMPLEX
#ifdef LONG__DOUBLE      
      call lis_array_set_all(nn*nn,(0.0q0,0.0q0),a,ierr)
#else
      call lis_array_set_all(nn*nn,(0.0d0,0.0d0),a,ierr)
#endif      
#else
#ifdef LONG__DOUBLE            
      call lis_array_set_all(nn*nn,0.0q0,a,ierr)
#else
      call lis_array_set_all(nn*nn,0.0d0,a,ierr)
#endif      
#endif      

      nnz = 0
#ifdef COMPLEX      
      do ii=0,nn-1
         i = ii/m
         j = ii - i*m
         if (i .GT. 0) then
            jj = ii - m
            a(ii + nn * jj) = (-1.0d0,0.0d0)
            nnz = nnz + 1
         end if
         if (i .LT. n-1) then
            jj = ii + m
            a(ii + nn * jj) = (-1.0d0,0.0d0)
            nnz = nnz + 1
         end if
         if (j .GT. 0) then
            jj = ii - 1
            a(ii + nn * jj) = (-1.0d0,0.0d0)
            nnz = nnz + 1
         end if
         if (j .LT. m-1) then
            jj = ii + 1
            a(ii + nn * jj) = (-1.0d0,0.0d0)
            nnz = nnz + 1
         end if
         jj = ii
         a(ii + nn * jj) = (4.0d0,0.0d0)
         nnz = nnz + 1
      end do
#else
      do ii=0,nn-1
         i = ii/m
         j = ii - i*m
         if (i .GT. 0) then
            jj = ii - m
            a(ii + nn * jj) = -1.0d0
            nnz = nnz + 1
         end if
         if (i .LT. n-1) then
            jj = ii + m
            a(ii + nn * jj) = -1.0d0
            nnz = nnz + 1
         end if
         if (j .GT. 0) then
            jj = ii - 1
            a(ii + nn * jj) = -1.0d0
            nnz = nnz + 1
         end if
         if (j .LT. m-1) then
            jj = ii + 1
            a(ii + nn * jj) = -1.0d0
            nnz = nnz + 1
         end if
         jj = ii
         a(ii + nn * jj) = 4.0d0
         nnz = nnz + 1
      end do
#endif      
      write(*,'(a,i0,a,i0,a,i0,a)') 'matrix size = ', nn, ' x ', nn, ' (', nnz, ' nonzero entries)'
      write(*,'(a)')

#ifdef COMPLEX
#ifdef LONG__DOUBLE
      call lis_array_set_all(nn,(1.0q0,0.0q0),u,ierr)
#else
      call lis_array_set_all(nn,(1.0d0,0.0d0),u,ierr)
#endif      
#else
#ifdef LONG__DOUBLE      
      call lis_array_set_all(nn,1.0q0,u,ierr)
#else
      call lis_array_set_all(nn,1.0d0,u,ierr)
#endif      
#endif      
      call lis_array_matvec(nn,a,u,b,LIS_INS_VALUE,ierr)

! solve linear system

      time0 = lis_wtime()
      call lis_array_solve(nn,a,b,x,w,ierr)
      time = lis_wtime() - time0

#ifdef COMPLEX
#ifdef LONG__DOUBLE
      call lis_array_xpay(nn,x,(-1.0q0,0.0q0),u,ierr)
#else      
      call lis_array_xpay(nn,x,(-1.0d0,0.0d0),u,ierr)
#endif      
#else
#ifdef LONG__DOUBLE      
      call lis_array_xpay(nn,x,-1.0q0,u,ierr)
#else
      call lis_array_xpay(nn,x,-1.0d0,u,ierr)
#endif
#endif      
      call lis_array_nrm2(nn,u,resid_r,ierr)
      call lis_array_nrm2(nn,b,resid_b,ierr)

      write(*,'(a,e14.7e2,a)') 'Direct: elapsed time         = ',time, ' sec.'
      write(*,'(a,e14.7e2,a)') 'Direct:   linear solver      = ',time, ' sec.'
      write(*,'(a,e14.7e2)') 'Direct: relative residual    = ',resid_r/resid_b
      write(*,'(a)') ''
      
      deallocate(a,b,x,u,w)
      
      call lis_finalize(ierr)

      stop
      end
      
