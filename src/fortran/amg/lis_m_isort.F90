!C   Copyright (C) The Scalable Software Infrastructure Project. All rights reserved.
!C
!C   Redistribution and use in source and binary forms, with or without
!C   modification, are permitted provided that the following conditions are met:
!C   1. Redistributions of source code must retain the above copyright
!C      notice, this list of conditions and the following disclaimer.
!C   2. Redistributions in binary form must reproduce the above copyright
!C      notice, this list of conditions and the following disclaimer in the
!C      documentation and/or other materials provided with the distribution.
!C   3. Neither the name of the project nor the names of its contributors 
!C      may be used to endorse or promote products derived from this software 
!C      without specific prior written permission.
!C
!C   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
!C   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!C   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!C   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
!C   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!C   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!C   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!C   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!C   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!C   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!C   POSSIBILITY OF SUCH DAMAGE.

#ifdef LONGLONG
#define LIS_MPI_INTEGER MPI_INTEGER8
#else
#define LIS_MPI_INTEGER MPI_INTEGER
#endif

!C   ************************************************
!C   * MODULE isort
!C     CONTAINS
!C   * SUBROUTINE qsorti
!C   ************************************************

module isort
contains
  subroutine  qsorti(ord,N,A)
    implicit none
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    integer(kind=kint) :: N
    integer(kind=kint),DIMENSION(N) ::ORD
    integer(kind=kint),DIMENSION(N) ::A

    integer(kind=kint) :: i
    DO i=1,N
       ord(i)=i
    END DO
    call qsort_rec(ord,N,A,1,N)
    contains

      recursive subroutine qsort_rec(ord,N,A,l,h)
        implicit none
#ifdef LONGLONG
        include 'precision_longlong.inc'
#else
        include 'precision.inc'
#endif
        integer(kind=kint) :: N,l,h
        integer(kind=kint),DIMENSION(N) ::ORD
        integer(kind=kint),DIMENSION(N) ::A

        integer(kind=kint) :: ii,jj,i,j,k,p,temp,t

        if(l>=h) return

        p=(l+h)/2
        temp=ord(p);ord(p)=ord(l);ord(l)=temp;
        t=A(temp);i=l;j=h;

        DO
           DO ii=i+1,h
              k=ord(ii)
              if(A(k)>=t) exit
           END DO
           i=ii
           DO jj=j,l+1,-1
              k=ord(jj)
              if(A(k)<t) exit
           END DO
           j=jj

           if(i>j) exit
           temp=ord(i);ord(i)=ord(j);ord(j)=temp;
        END DO
        temp=ord(l);ord(l)=ord(j);ord(j)=temp;

        call qsort_rec(ord,N,A,l,j-1)
        call qsort_rec(ord,N,A,j+1,h)
      end subroutine qsort_rec
      
  end subroutine qsorti
end module isort









