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

#ifndef ZERO_ORIGIN
#define ZERO_ORIGIN 0
#endif

#ifdef LONGLONG
#define LIS_MPI_INTEGER MPI_INTEGER8
#else
#define LIS_MPI_INTEGER MPI_INTEGER
#endif

!C   ************************************************
!C   * MODULE solver_SR2
!C     CONTAINS
!C   * subroutine  SOLVER_SEND_RECV2I
!C   * subroutine  SOLVER_SEND_RECV2
!C   ************************************************

module solver_SR2
  contains
  subroutine  SOLVER_SEND_RECV2I                                                &
       &                ( N, S, R, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
       &                                        STACK_EXPORT, NOD_EXPORT,       &
       &                  WS, WR, X, SOLVER_COMM, my_rank)
    
    implicit none
    include  'mpif.h'
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif

    ! ......................................................................

    integer(kind=kint )                , intent(in)   ::  N
    integer(kind=kint )                , intent(in)   ::  S
    integer(kind=kint )                , intent(in)   ::  R
    integer(kind=kint )                , intent(in)   ::  NEIBPETOT
    integer(kind=kint ) :: NEIBPE      (1:NEIBPETOT)
    integer(kind=kint ) :: STACK_IMPORT(0:NEIBPETOT)
    integer(kind=kint ) :: NOD_IMPORT  (1:R)
    integer(kind=kint ) :: STACK_EXPORT(0:NEIBPETOT)
    integer(kind=kint ) :: NOD_EXPORT  (1:S)
    integer(kind=kint ), intent(inout):: WS(1:S)
    integer(kind=kint ), intent(inout):: WR(1:R)
    !C- dimension(R)  
    INTEGER(kind=kint ), intent(inout):: X(1:N)
    !C- dimension(N)  
    integer(kind=kint ), intent(in)   :: SOLVER_COMM
    integer(kind=kint ), intent(in)   :: my_rank

    integer(kind=kint ), dimension(:,:), save, allocatable :: sta1
    integer(kind=kint ), dimension(:,:), save, allocatable :: sta2
    integer(kind=kint ), dimension(:  ), save, allocatable :: req1
    integer(kind=kint ), dimension(:  ), save, allocatable :: req2  
    
    integer(kind=kint) :: neib, istart, inum, k, ierr
    integer(kind=kint), save :: NFLAG = 0
    integer(kind=kint), save :: EARLIER_NEIBPETOT = 0


    !C
    !C-- INIT.
    if (NFLAG.eq.0) then
       allocate (sta1(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (sta2(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (req1(NEIBPETOT))
       allocate (req2(NEIBPETOT))
       NFLAG = 1
       EARLIER_NEIBPETOT = NEIBPETOT
    else if(EARLIER_NEIBPETOT < NEIBPETOT) then
       deallocate(sta1, sta2, req1, req2)
       allocate (sta1(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (sta2(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (req1(NEIBPETOT))
       allocate (req2(NEIBPETOT))
       EARLIER_NEIBPETOT = NEIBPETOT
    endif

    
    !C
    !C-- SEND
    do neib= 1, NEIBPETOT
       istart= STACK_EXPORT(neib-1)
       inum  = STACK_EXPORT(neib  ) - istart
       do k= istart+1, istart+inum
          WS(k)= X(NOD_EXPORT(k)+ZERO_ORIGIN)

       enddo
       call MPI_ISEND (WS(istart+1), inum, LIS_MPI_INTEGER,       &
            &                  NEIBPE(neib), 13, SOLVER_COMM,  &
            &                  req1(neib), ierr)
    enddo

    !C
    !C-- RECEIVE
    do neib= 1, NEIBPETOT
       istart= STACK_IMPORT(neib-1)
       inum  = STACK_IMPORT(neib  ) - istart
       call MPI_IRECV (WR(istart+1), inum, LIS_MPI_INTEGER,        &
            &                  NEIBPE(neib), 13, SOLVER_COMM,  &
            &                  req2(neib), ierr)
    enddo
    call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
    
    do neib= 1, NEIBPETOT
       istart= STACK_IMPORT(neib-1)
       inum  = STACK_IMPORT(neib  ) - istart
       do k= istart+1, istart+inum
          X(NOD_IMPORT(k)+ZERO_ORIGIN)= WR(k)

       enddo
    enddo
    
    call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)

  end subroutine solver_send_recv2I


  subroutine  SOLVER_SEND_RECV2                                                 &
       &                ( N, S, R, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
       &                                        STACK_EXPORT, NOD_EXPORT,       &
       &                  WS, WR, X, SOLVER_COMM, my_rank)
    implicit none
    include  'mpif.h'
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif

    ! ......................................................................

    integer(kind=kint), intent(in)   ::  N
    integer(kind=kint), intent(in)   ::  S
    integer(kind=kint), intent(in)   ::  R
    integer(kind=kint), intent(in)   ::  NEIBPETOT

    integer(kind=kint) :: NEIBPE      (1:NEIBPETOT)
    integer(kind=kint) :: STACK_IMPORT(0:NEIBPETOT)
    integer(kind=kint) :: NOD_IMPORT  (1:R)
    integer(kind=kint) :: STACK_EXPORT(0:NEIBPETOT)
    integer(kind=kint) :: NOD_EXPORT  (1:S)
    real   (kind=kreal), intent(inout):: WS(:)
    !C--  WS(S)
    real   (kind=kreal), intent(inout):: WR(:)
    !C--  WR(R)
    real   (kind=kreal), intent(inout):: X(1:N)
    !C-- X(N)

    integer(kind=kint), intent(in)   :: SOLVER_COMM
    integer(kind=kint), intent(in)   :: my_rank

    ! ....................................................................
    integer(kind=kint), dimension(:,:), save, allocatable :: sta1
    integer(kind=kint), dimension(:,:), save, allocatable :: sta2
    integer(kind=kint), dimension(:  ), save, allocatable :: req1
    integer(kind=kint), dimension(:  ), save, allocatable :: req2  
    
    integer(kind=kint) :: neib, istart, inum, k, ierr,j,l
    integer(kind=kint), save :: NFLAG = 0
    integer(kind=kint), save :: EARLIER_NEIBPETOT = 0
    logical :: flag1, flag2


#ifdef DEBUG
    write(*,*) EARLIER_NEIBPETOT,NEIBPETOT, my_rank,"send_recv_2"    
#endif    
    !C
    !C-- INIT.
    if (NFLAG.eq.0) then
       allocate (sta1(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (sta2(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (req1(NEIBPETOT))
       allocate (req2(NEIBPETOT))
       NFLAG = 1
       EARLIER_NEIBPETOT = NEIBPETOT
    else if(EARLIER_NEIBPETOT < NEIBPETOT) then
       deallocate(sta1, sta2, req1, req2)
       allocate (sta1(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (sta2(MPI_STATUS_SIZE, NEIBPETOT))
       allocate (req1(NEIBPETOT))
       allocate (req2(NEIBPETOT))
       EARLIER_NEIBPETOT = NEIBPETOT
    endif


    !C
    !C-- SEND
    do neib = 1, NEIBPETOT
       istart = STACK_EXPORT(neib - 1)
       inum  = STACK_EXPORT(neib) - istart
       
       do k = istart + 1, istart + inum
          WS(k) = X(NOD_EXPORT(k)+ZERO_ORIGIN)
       enddo
       call MPI_ISEND (WS(istart + 1), inum, MPI_DOUBLE_PRECISION,     &
            &          NEIBPE(neib), 0, SOLVER_COMM, req1(neib), ierr)
       
    enddo

    !C
    !C-- RECEIVE

    do neib = 1, NEIBPETOT
       istart = STACK_IMPORT(neib - 1)
       inum   = STACK_IMPORT(neib) - istart
       CALL MPI_IRECV (WR(istart + 1), inum, MPI_DOUBLE_PRECISION,     &
            &          NEIBPE(neib), 0, SOLVER_COMM, req2(neib), ierr)
    enddo
    
    call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
#ifdef DEBUG
    call MPI_TESTALL (NEIBPETOT, req2, flag2, sta2, ierr)
#endif
    
    do neib = 1, NEIBPETOT
       istart = STACK_IMPORT(neib - 1)
       inum  = STACK_IMPORT(neib) - istart
       do k = istart + 1, istart + inum
          X(NOD_IMPORT(k)+ZERO_ORIGIN) = WR(k)
       enddo
    enddo
    
    call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)
#ifdef DEBUG
    call MPI_TESTALL (NEIBPETOT, req1, flag1, sta1, ierr)
#endif



#ifdef DEBUG
    if(.not.(flag1.and.flag2).and.NEIBPETOT>0) then
       write(*,*) "error in send_recv",flag1,flag2
       stop "hehe"
    end if

    write(*,*) "exit send_recv_2"    
#endif    

  end subroutine solver_send_recv2
    
end module solver_SR2
