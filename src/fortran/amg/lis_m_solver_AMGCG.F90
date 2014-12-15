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
!C   * MODULE solver_AMGCG
!C     CONTAINS
!C   * SUBROUTINE clear_matrix
!C   * SUBROUTINE clear_matrix_ssi_amg
!C   * SUBROUTINE v_cycle_ssi_amg
!C   * SUBROUTINE v_cycle
!C   * SUBROUTINE ll_slv
!C   * SUBROUTINE sgs
!C   * SUBROUTINE matrix_counting
!C   * SUBROUTINE matrix_arrange
!C   * SUBROUTINE AMGCG
!C   ************************************************

MODULE solver_AMGCG
CONTAINS

  SUBROUTINE clear_matrix(LEVEL_NUM)
    use data_structure_for_AMG
    IMPLICIT NONE
    INTEGER(kind=kint) :: i, LEVEL_NUM
    NULLIFY(       HIERARCHICAL_DATA(1) % INU,      &
            &      HIERARCHICAL_DATA(1) % INL,      &
            &      HIERARCHICAL_DATA(1) % IAU,      &
            &      HIERARCHICAL_DATA(1) % IAL,      &
            &      HIERARCHICAL_DATA(1) % AU,       &
            &      HIERARCHICAL_DATA(1) % AL,       &
            &      HIERARCHICAL_DATA(1) % D,        &
            &      HIERARCHICAL_DATA(1) % X,        &
            &      HIERARCHICAL_DATA(1) % B)


    DO i=2,LEVEL_NUM
       deallocate( HIERARCHICAL_DATA(i) % INU,      &
            &      HIERARCHICAL_DATA(i) % INL,      &
            &      HIERARCHICAL_DATA(i) % IAU,      &
            &      HIERARCHICAL_DATA(i) % IAL,      &
            &      HIERARCHICAL_DATA(i) % AU,       &
            &      HIERARCHICAL_DATA(i) % AL,       &
            &      HIERARCHICAL_DATA(i) % D,        &
            &      HIERARCHICAL_DATA(i) % X,        &
            &      HIERARCHICAL_DATA(i) % B,        &
            &      HIERARCHICAL_DATA(i) % R % IN,   &
            &      HIERARCHICAL_DATA(i) % R % CN,   &
            &      HIERARCHICAL_DATA(i) % R % V,    &
            &      HIERARCHICAL_DATA(i) % P % IN,   &
            &      HIERARCHICAL_DATA(i) % P % CN,   &
            &      HIERARCHICAL_DATA(i) % P % V     )
       deallocate( HIERARCHICAL_DATA(i) % P,        &
            &      HIERARCHICAL_DATA(i) % R)

       if(HIERARCHICAL_DATA(i) % COMM_TABLE % NEIBPETOT > 0) then
          deallocate(HIERARCHICAL_DATA(i) % COMM_TABLE % NEIBPE,       &
               &     HIERARCHICAL_DATA(i) % COMM_TABLE % STACK_IMPORT, & 
               &     HIERARCHICAL_DATA(i) % COMM_TABLE % STACK_EXPORT, & 
               &     HIERARCHICAL_DATA(i) % COMM_TABLE % NOD_IMPORT,   &
               &     HIERARCHICAL_DATA(i) % COMM_TABLE % NOD_EXPORT)
       end if
       if(HIERARCHICAL_DATA(i) % INT_LVL_TABLE % NEIBPETOT > 0) then
          deallocate(HIERARCHICAL_DATA(i) % INT_LVL_TABLE % NEIBPE, &
               &     HIERARCHICAL_DATA(i) % INT_LVL_TABLE % STACK_IMPORT, & 
               &     HIERARCHICAL_DATA(i) % INT_LVL_TABLE % STACK_EXPORT, & 
               &     HIERARCHICAL_DATA(i) % INT_LVL_TABLE % NOD_IMPORT,   &
               &     HIERARCHICAL_DATA(i) % INT_LVL_TABLE % NOD_EXPORT)
       end if

    END DO

  END SUBROUTINE clear_matrix

  subroutine clear_matrix_ssi_amg(LEVEL_NUM)
    IMPLICIT NONE
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    integer(kind=kint) :: LEVEL_NUM
    
    call clear_matrix(LEVEL_NUM)
    
  end subroutine clear_matrix_ssi_amg

  !C B(NP), X(NP), TEMP(NP), WS(WSIZE), WR(WSIZE)
  SUBROUTINE v_cycle_ssi_amg(B, X, TEMP, LEVEL_NUM, SOLVER_COMM, WS, WR, NP, WSIZE)
    IMPLICIT NONE
    include 'mpif.h'
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    INTEGER(kind=kint), intent(in) :: WSIZE, NP
    REAL   (kind=kreal) :: B(1:NP), X(1:NP), TEMP(1:NP), WS(1:WSIZE), WR(1:WSIZE)
    INTEGER(kind=kint), intent(in) :: LEVEL_NUM, SOLVER_COMM
    INTEGER(kind=kint) :: ierr, my_rank, NPROCS

    CALL MPI_COMM_SIZE (SOLVER_COMM, NPROCS, ierr) 
    CALL MPI_COMM_RANK (SOLVER_COMM, my_rank, ierr) 

    X(1:NP) = 0.0d0
    call v_cycle(B(1:NP), X(1:NP), LEVEL_NUM, SOLVER_COMM, my_rank, WS(1:WSIZE), WR(1:WSIZE), Temp(1:NP), NPROCS, WSIZE, NP)
    
  END SUBROUTINE v_cycle_ssi_amg

  
  SUBROUTINE v_cycle(problem_B, problem_X, LEVEL_NUM, SOLVER_COMM, my_rank, WS, WR, Temp, &
       &             NPROCS, dim, problem_NP)
    USE solver_SR2
    USE data_structure_for_AMG
    USE solver_Gnumbering
    USE count_time_mod
    
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER(kind=kint), INTENT(in)    :: NPROCS
    INTEGER(kind=kint), INTENT(in)    :: dim, problem_NP
    INTEGER(kind=kint), INTENT(in)    :: LEVEL_NUM
    INTEGER(kind=kint)                :: N,NP,NPL,NPU
    REAL   (kind=kreal), DIMENSION(:), pointer:: D
    !C DIMENSION(N)
    REAL   (kind=kreal), target::  problem_B(1:problem_NP)
    !C DIMENSION(N)
    REAL   (kind=kreal), target::  problem_X(1:problem_NP)   
    !C DIMENSION(N)

    REAL   (kind=kreal), DIMENSION(:), pointer::  B
    !C DIMENSION(N)
    REAL   (kind=kreal), DIMENSION(:), pointer::  X
    !C DIMENSION(N)

    REAL   (kind=kreal), DIMENSION(:), pointer::  AU
    REAL   (kind=kreal), DIMENSION(:), pointer::  AL
    
    INTEGER(kind=kint), DIMENSION(:), pointer ::  INU
    !C DIMENSION(0:N)
    INTEGER(kind=kint), DIMENSION(:), pointer ::  IAU
    !C DIMENSION(NPU)
    INTEGER(kind=kint), DIMENSION(:), pointer ::  INL
    !C DIMENSION(0:N)
    INTEGER(kind=kint), DIMENSION(:), pointer ::  IAL
    
    INTEGER(kind=kint), intent(in) :: SOLVER_COMM
    INTEGER(kind=kint), intent(in) :: my_rank
    INTEGER(kind=kint) :: ierr
    
    REAL   (kind=kreal), DIMENSION(:), pointer :: coarser_B
    REAL   (kind=kreal), DIMENSION(:), pointer :: coarser_X
    REAL   (kind=kreal)  :: Temp(1:problem_NP)
    REAL   (kind=kreal)  :: WS(1:dim), WR(1:dim)
    
    INTEGER(kind=kint ) :: nth_lev,i,inod,j,ieL,isL,is,ie,coarser_N,isU,ieU,k,coarser_NP
    INTEGER(kind=kint ) :: coarser_NEIBPETOT
    REAL   (kind=kreal) :: R_val,B_val,w,X_val,R_norm,GR_norm
    
    type(INTER_LEVEL_OPERATOR), pointer :: R, P

    INTEGER(kind=kint) :: NEIBPETOT, iNEIBPETOT
    INTEGER(kind=kint), pointer :: NEIBPE(:)
    INTEGER(kind=kint), pointer :: STACK_IMPORT(:),STACK_EXPORT(:)
    INTEGER(kind=kint), pointer :: NOD_IMPORT(:),NOD_EXPORT(:)

    
#ifdef CHOLESCKY
    INTEGER(kind=kint) :: recvcounts(0:NPROCS-1)
#endif

    logical :: print_residual_on_each_level, print_residual_on_level_one
    print_residual_on_each_level = .false.
    print_residual_on_level_one = .false.

    HIERARCHICAL_DATA(1) % B => problem_B
    HIERARCHICAL_DATA(1) % X => problem_X

    
    DO nth_lev=1,LEVEL_NUM-1
       N  = HIERARCHICAL_DATA(nth_lev) % N
       NP = HIERARCHICAL_DATA(nth_lev) % NP
       NPL= HIERARCHICAL_DATA(nth_lev) % NPL
       NPU= HIERARCHICAL_DATA(nth_lev) % NPU
       INU=>HIERARCHICAL_DATA(nth_lev) % INU 
       INL=>HIERARCHICAL_DATA(nth_lev) % INL 
       X => HIERARCHICAL_DATA(nth_lev) % X
       B => HIERARCHICAL_DATA(nth_lev) % B
       D => HIERARCHICAL_DATA(nth_lev) % D
       IAU=>HIERARCHICAL_DATA(nth_lev) % IAU
       IAL=>HIERARCHICAL_DATA(nth_lev) % IAL
       AU =>HIERARCHICAL_DATA(nth_lev) % AU
       AL =>HIERARCHICAL_DATA(nth_lev) % AL
       NEIBPETOT    =  HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NEIBPETOT
       STACK_IMPORT => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % STACK_IMPORT
       STACK_EXPORT => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % STACK_EXPORT
       NOD_IMPORT   => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NOD_IMPORT
       NOD_EXPORT   => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NOD_EXPORT
       NEIBPE       => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NEIBPE

       
       NEIBPETOT = HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NEIBPETOT
       
       coarser_B => HIERARCHICAL_DATA(nth_lev+1) % B
       coarser_X => HIERARCHICAL_DATA(nth_lev+1) % X
       
#ifdef CHOLESCKY
       if(nth_lev==LEVEL_NUM-1) then
          coarser_B => HIERARCHICAL_DATA(nth_lev+1) % X
       end if
#endif
       !C calculate residual
       if(print_residual_on_each_level) Then
          DO j= 1, N
             R_val= B(j) - D(j)*X(j)
             isU= INU(j-1)+1
             ieU= INU(j)
             DO i= isU, ieU
                inod=IAU(i)+ZERO_ORIGIN
                R_val= R_val - AU(i) * X(inod)
             END DO
             isL= INL(j-1)+1
             ieL= INL(j)
             DO i= isL, ieL
                inod= IAL(i)+ZERO_ORIGIN
                R_val= R_val - AL(i) * X(inod)
             END DO
             Temp(j)=R_val
          END DO

          R_norm=0.0
          DO i=1,N
             R_norm = R_norm + Temp(i)**2
          END DO
          CALL MPI_REDUCE(R_norm,GR_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,SOLVER_COMM,ierr)
          if(my_rank==0) then
             GR_norm = dsqrt(GR_norm)
             write(*,*) nth_lev, 'pre ',log10(GR_norm)
          end if
       end if
       CALL count_time(1,SOLVER_COMM,my_rank,8)


       CALL sgs(N,NP,NPL,NPU,D(1:NP),AL(1:NPL),INL(0:N),IAL(1:NPL),AU(1:NPU),INU(0:N),IAU(1:NPU),B(1:NP),&
            & X(1:NP),1,nth_lev,SOLVER_COMM,my_rank,WS(1:dim),WR(1:dim),dim, nth_lev)

       CALL count_time(2,SOLVER_COMM,my_rank,8)
       !C barrier in ssor ##############

       !C calculate residual
       DO j= 1, N
          R_val= B(j) - D(j)*X(j)
          isU= INU(j-1)+1
          ieU= INU(j)
          DO i= isU, ieU
             inod  = IAU(i) + ZERO_ORIGIN
             R_val = R_val - AU(i) * X(inod)
          END DO
          isL= INL(j-1)+1
          ieL= INL(j)
          DO i= isL, ieL
             inod  = IAL(i) + ZERO_ORIGIN
             R_val = R_val - AL(i) * X(inod)
          END DO
          Temp(j)=R_val
       END DO


       IF((nth_lev == 1 .and. print_residual_on_level_one).or. print_residual_on_each_level) THEN
          R_norm=0.0
          DO i=1,N
             R_norm = R_norm + Temp(i)**2
          END DO
          CALL MPI_REDUCE(R_norm,GR_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,SOLVER_COMM,ierr)
          if(my_rank==0) then
             GR_norm = dsqrt(GR_norm)
             write(*,*) nth_lev,'post',log10(GR_norm)
          end if
       END IF
       
       
       !C restrict the residual vector
!!$       R => HIERARCHICAL_DATA(nth_lev+1) % R
       coarser_N= HIERARCHICAL_DATA(nth_lev+1) % N
       coarser_NP=HIERARCHICAL_DATA(nth_lev+1) % NP
       coarser_NEIBPETOT=HIERARCHICAL_DATA(nth_lev+1) % COMM_TABLE % NEIBPETOT
       

       DO j= 1, coarser_NP
          is= HIERARCHICAL_DATA(nth_lev+1) % R % IN(j-1)+1
          ie= HIERARCHICAL_DATA(nth_lev+1) % R % IN(j)
          B_val=0.0
          DO i=is, ie
             inod = HIERARCHICAL_DATA(nth_lev+1) % R % CN(i)
             B_val= B_val + HIERARCHICAL_DATA(nth_lev+1) % R % V(i) * Temp(inod)
          END DO
          coarser_B(j)= B_val
       END DO

       iNEIBPETOT = HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % NEIBPETOT

       CALL SOLVER_SR_REDUCE(coarser_NP, dim,&
            & HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_EXPORT(iNEIBPETOT),    &
            & HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_IMPORT(iNEIBPETOT),    &
            & iNEIBPETOT,                                                                 &
            & HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % NEIBPE(1:iNEIBPETOT),        &
            & HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_IMPORT(0:iNEIBPETOT),  &
            & HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % NOD_IMPORT                   &
            & (1:HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_IMPORT(iNEIBPETOT)),&
            & HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_EXPORT(0:iNEIBPETOT),  &
            & HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % NOD_EXPORT                   &
            & (1:HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_EXPORT(iNEIBPETOT)),&
            & WS(1:dim), WR(1:dim), coarser_B(1:coarser_NP), SOLVER_COMM, my_rank, 1)



       DO j=1,coarser_NP
          coarser_X(j)= 0.0
       END DO
    END DO

    N = HIERARCHICAL_DATA(LEVEL_NUM) % N
    NP = HIERARCHICAL_DATA(LEVEL_NUM) % NP
    NPL = HIERARCHICAL_DATA(LEVEL_NUM) % NPL
    NPU = HIERARCHICAL_DATA(LEVEL_NUM) % NPU
    NEIBPETOT    =  HIERARCHICAL_DATA(LEVEL_NUM) % COMM_TABLE % NEIBPETOT 
    NEIBPE       => HIERARCHICAL_DATA(LEVEL_NUM) % COMM_TABLE % NEIBPE
    STACK_IMPORT => HIERARCHICAL_DATA(LEVEL_NUM) % COMM_TABLE % STACK_IMPORT
    STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NUM) % COMM_TABLE % STACK_EXPORT
    NOD_IMPORT   => HIERARCHICAL_DATA(LEVEL_NUM) % COMM_TABLE % NOD_IMPORT
    NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NUM) % COMM_TABLE % NOD_EXPORT

#ifdef PSPASES
    X => HIERARCHICAL_DATA(LEVEL_NUM) % X
    B => HIERARCHICAL_DATA(LEVEL_NUM) % B

#elif defined(CHOLESCKY)
    X   => HIERARCHICAL_DATA(LEVEL_NUM) % X
    B   => HIERARCHICAL_DATA(LEVEL_NUM) % B
    AL  => HIERARCHICAL_DATA(LEVEL_NUM) % AL
    INL => HIERARCHICAL_DATA(LEVEL_NUM) % INL 
#else
    X   => HIERARCHICAL_DATA(LEVEL_NUM) % X
    B   => HIERARCHICAL_DATA(LEVEL_NUM) % B
    AL  => HIERARCHICAL_DATA(LEVEL_NUM) % AL
    INL => HIERARCHICAL_DATA(LEVEL_NUM) % INL 
    INU => HIERARCHICAL_DATA(LEVEL_NUM) % INU 
    D   => HIERARCHICAL_DATA(LEVEL_NUM) % D
    IAU => HIERARCHICAL_DATA(LEVEL_NUM) % IAU
    IAL => HIERARCHICAL_DATA(LEVEL_NUM) % IAL
    AU  => HIERARCHICAL_DATA(LEVEL_NUM) % AU
    
    CALL count_time(1,SOLVER_COMM,my_rank,8)


    CALL sgs(N,NP,NPL,NPU,D(1:NP),AL(1:NPL),INL(0:N),IAL(1:NPL),AU(1:NPU),INU(0:N),IAU(1:NPU),B(1:NP),&
         & X(1:NP),30,LEVEL_NUM,SOLVER_COMM,my_rank,WS(1:dim),WR(1:dim),dim, LEVEL_NUM)

    CALL count_time(2,SOLVER_COMM,my_rank,8)

    !C barrier in ssor ########
#endif

#ifdef PSPASES
    CALL DPSPACET(rowdistbx,nrhs,coarsest_B,ldb,coarsest_X,ldx,ioptions,pspcommf,SOLVER_COMM)
!!$    CALL SOLVER_SEND_RECV2(NP,dim,dim,NEIBPETOT,NEIBPE,STACK_IMPORT,&
!!$         & NOD_IMPORT,STACK_EXPORT, NOD_EXPORT, WS, WR, coarsest_X, &
!!$         & SOLVER_COMM, my_rank)

    
#elif defined(CHOLESCKY)
    !C--X has data of each PE's B
    X(1:N) = B(1:N)
    DO i = 1, NPROCS
       recvcounts(i-1) = INL(i) - INL(i - 1)
    END DO
    CALL MPI_ALLGATHERV(X, N, MPI_DOUBLE_PRECISION, B, recvcounts(0), INL(0), &
         & MPI_DOUBLE_PRECISION, SOLVER_COMM, ierr)
    
    CALL ll_slv(AL, B, INL(NPROCS))

    i = INL(my_rank) + 1
    j = INL(my_rank + 1)
    DO k = i, j
       X(k - i + 1) = B(k)
    END DO

    if(NEIBPETOT > 0) then
       CALL SOLVER_SEND_RECV2(NP, STACK_EXPORT(NEIBPETOT), STACK_IMPORT(NEIBPETOT), &
            & NEIBPETOT, NEIBPE(1:NEIBPETOT), STACK_IMPORT(0:NEIBPETOT), &
            & NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)), STACK_EXPORT(0:NEIBPETOT), &
            & NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)), &
            & WS, WR, X(1:NP), SOLVER_COMM, my_rank)
    end if

    !C-- end: X has data of each PE's B
#endif

    !C calculate residual
    if(print_residual_on_each_level)Then
       DO j= 1, N
!!$          CB is used in the case LU
!!$          R_val= CB(j) - D(j)*X(j)

          R_val = B(j)-D(j)*X(j)
          isU= INU(j-1)+1
          ieU= INU(j)
          DO i= isU, ieU
             inod=IAU(i)+ZERO_ORIGIN
             R_val= R_val - AU(i) * X(inod)
          END DO
          isL= INL(j-1)+1
          ieL= INL(j)
          DO i= isL, ieL
             inod= IAL(i)+ZERO_ORIGIN
             R_val= R_val - AL(i) * X(inod)
          END DO
          TEMP(j)=R_val
       END DO

       
       R_norm=0.0
       DO i=1,N
          R_norm = R_norm + TEMP(i)**2
       END DO


       
       CALL MPI_REDUCE(R_norm,GR_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,SOLVER_COMM,ierr)
       if(my_rank==0) then
          GR_norm = dsqrt(GR_norm)
          write(*,*) nth_lev,'DIRECT',log10(GR_norm)
       end if
    end if

    DO nth_lev=LEVEL_NUM-1,1,-1
       N= HIERARCHICAL_DATA(nth_lev) % N
       NP=HIERARCHICAL_DATA(nth_lev) % NP
       NPL=HIERARCHICAL_DATA(nth_lev) % NPL
       NPU=HIERARCHICAL_DATA(nth_lev) % NPU

       INU          => HIERARCHICAL_DATA(nth_lev) % INU
       INL          => HIERARCHICAL_DATA(nth_lev) % INL
       X            => HIERARCHICAL_DATA(nth_lev) % X
       B            => HIERARCHICAL_DATA(nth_lev) % B
       D            => HIERARCHICAL_DATA(nth_lev) % D
       IAU          => HIERARCHICAL_DATA(nth_lev) % IAU
       IAL          => HIERARCHICAL_DATA(nth_lev) % IAL
       AU           => HIERARCHICAL_DATA(nth_lev) % AU
       AL           => HIERARCHICAL_DATA(nth_lev) % AL
       STACK_IMPORT => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % STACK_IMPORT
       STACK_EXPORT => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % STACK_EXPORT
       NOD_IMPORT   => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NOD_IMPORT
       NOD_EXPORT   => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NOD_EXPORT
       NEIBPE       => HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NEIBPE
       NEIBPETOT    =  HIERARCHICAL_DATA(nth_lev) % COMM_TABLE % NEIBPETOT
       coarser_X    => HIERARCHICAL_DATA(nth_lev+1) % X
!!$       P => HIERARCHICAL_DATA(nth_lev+1) % P

       IF(NEIBPETOT > 0) THEN
          iNEIBPETOT = HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % NEIBPETOT 
          CALL SOLVER_SEND_RECV2 &
               & (HIERARCHICAL_DATA(nth_lev+1) % NP, &
               &  HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_EXPORT(iNEIBPETOT), &
               &  HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_IMPORT(iNEIBPETOT), &
               &  iNEIBPETOT, &
               &  HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % NEIBPE(1:iNEIBPETOT),     &
               &  HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % &
               &   STACK_IMPORT(0:iNEIBPETOT), &
               &  HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % &
               &   NOD_IMPORT(1:HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_IMPORT(iNEIBPETOT)), &
               &  HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_EXPORT(0:iNEIBPETOT),              &
               &  HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % &
               &   NOD_EXPORT(1:HIERARCHICAL_DATA(nth_lev+1) % INT_LVL_TABLE % STACK_EXPORT(iNEIBPETOT)), &
               &  WS(1:dim), WR(1:dim), coarser_X(1:HIERARCHICAL_DATA(nth_lev+1) % NP), SOLVER_COMM,      &
               &  my_rank)
       end IF

       DO j= 1, N
          is= HIERARCHICAL_DATA(nth_lev+1) % P % IN(j-1)+1
          ie= HIERARCHICAL_DATA(nth_lev+1) % P % IN(j)
          X_val=0.0
          DO i=is, ie
             inod = HIERARCHICAL_DATA(nth_lev+1) % P % CN(i)
             X_val= X_val+HIERARCHICAL_DATA(nth_lev+1) % P % V(i) * coarser_X(inod)
          END DO
          X(j)=X(j)+X_val
       END DO

       if(NEIBPETOT > 0) then       
          CALL SOLVER_SEND_RECV2(NP,STACK_EXPORT(NEIBPETOT),STACK_IMPORT(NEIBPETOT),&
               & NEIBPETOT,NEIBPE(1:NEIBPETOT),STACK_IMPORT(0:NEIBPETOT),&
               & NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)),STACK_EXPORT(0:NEIBPETOT), &
               & NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)), &
               & WS(1:dim), WR(1:dim), X(1:NP), SOLVER_COMM, my_rank)
       end if
       
       !C calculate residual
       if(print_residual_on_each_level) Then
          DO j= 1, N
             R_val= B(j) - D(j)*X(j)
             isU= INU(j-1)+1
             ieU= INU(j)
             DO i= isU, ieU
                inod=IAU(i)+ZERO_ORIGIN
                R_val= R_val - AU(i) * X(inod)
             END DO
             isL= INL(j-1)+1
             ieL= INL(j)
             DO i= isL, ieL
                inod= IAL(i)+ZERO_ORIGIN
                R_val= R_val - AL(i) * X(inod)
             END DO
             TEMP(j)=R_val
          END DO

          R_norm=0.0
          DO i=1,N
             R_norm = R_norm + TEMP(i)**2
          END DO

          CALL MPI_REDUCE(R_norm,GR_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,SOLVER_COMM,ierr)
          if(my_rank==0) then
             GR_norm = dsqrt(GR_norm)
             write(*,*) nth_lev, 'pre ',log10(GR_norm)
          end if
       end if
       
       CALL count_time(1,SOLVER_COMM,my_rank,8)

       CALL sgs(N,NP,NPL,NPU,D(1:NP),AL(1:NPL),INL(0:N),IAL(1:NPL),AU(1:NPU),INU(0:N),&
            & IAU(1:NPU),B(1:NP),X(1:NP),1,nth_lev,SOLVER_COMM,my_rank,WS(1:dim),WR(1:dim),dim, nth_lev)
       
       CALL count_time(2,SOLVER_COMM,my_rank,8)
       
       !C barrier in ssor ########

       !C calculate residual
       if(print_residual_on_each_level) Then
          DO j= 1, N
             R_val= B(j) - D(j)*X(j)
             isU= INU(j-1)+1
             ieU= INU(j)
             DO i= isU, ieU
                inod=IAU(i)+ZERO_ORIGIN
                R_val= R_val - AU(i) * X(inod)
             END DO
             isL= INL(j-1)+1
             ieL= INL(j)
             DO i= isL, ieL
                inod= IAL(i)+ZERO_ORIGIN
                R_val= R_val - AL(i) * X(inod)
             END DO
             TEMP(j)=R_val
          END DO

          R_norm=0.0
          DO i=1,N
             R_norm = R_norm + TEMP(i)**2
          END DO
          if(NEIBPETOT > 0) then
             CALL SOLVER_SEND_RECV2(NP,STACK_EXPORT(NEIBPETOT),STACK_IMPORT(NEIBPETOT), &
                  & NEIBPETOT,NEIBPE(1:NEIBPETOT),STACK_IMPORT(0:NEIBPETOT),  &
                  & NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)), STACK_EXPORT(0:NEIBPETOT), &
                  & NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)),&
                  & WS(1:dim), WR(1:dim), Temp(1:NP), SOLVER_COMM, my_rank)
          end if
          
          CALL MPI_REDUCE(R_norm,GR_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,SOLVER_COMM,ierr)
          if(my_rank==0) then
             write(*,*) nth_lev, 'post ',log10(GR_norm)/2
          end if
       end if
    END DO
  END SUBROUTINE v_cycle

  SUBROUTINE ll_slv(lower_mat, b, N)
    IMPLICIT NONE
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    REAL(kind=kreal),   intent(in) :: lower_mat(:)
    REAL(kind=kreal),   intent(out):: b(:)
    INTEGER(kind=kint), intent(in) :: N
    INTEGER(kind=kint)             ::i,j,k,cj,ci
    
    DO i = 1, N
       DO j = 1, i - 1
          cj = N * (j - 1) - j * (j - 1) / 2
          b(i) = b(i) - b(j) * lower_mat(i + cj)
       END DO
       b(i) = b(i) * lower_mat(i + N * (i - 1) - i * (i - 1) / 2)
    END DO
    DO i = N, 1, -1
       ci = N * (i - 1) - i * (i - 1) / 2
       DO j = i + 1, N
          b(i) = b(i) - b(j) * lower_mat(j + ci)
       END DO
       b(i) = b(i) * lower_mat(i + ci)
    END DO

  END SUBROUTINE ll_slv


  SUBROUTINE sgs(N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, B, X, count, LEVEL_NO, SOLVER_COMM, &
       &         my_rank, WS, WR, dim, level)
    USE solver_SR2
    USE data_structure_for_AMG
    USE solver_Gnumbering
    USE count_time_mod

    IMPLICIT NONE

    INTEGER(kind=kint), INTENT(in)   :: N, NP, NPL, NPU, count, dim, level

    REAL   (kind=kreal) :: D(1:NP), B(1:NP), X(1:NP), AU(1:NPU), AL(1:NPL)

    INTEGER(kind=kint) :: INU(0:N), INL(0:N), IAU(1:NPU), IAL(1:NPL)

    INTEGER(kind=kint), intent(in) :: SOLVER_COMM
    INTEGER(kind=kint), intent(in) :: my_rank
    INTEGER(kind=kint)             :: ierr
    INTEGER(kind=kint), intent(in) :: LEVEL_NO

    INTEGER(kind=kint), pointer :: NEIBPE(:)
    INTEGER(kind=kint), pointer :: STACK_IMPORT(:), STACK_EXPORT(:)
    INTEGER(kind=kint), pointer :: NOD_EXPORT(:), NOD_IMPORT(:)
    INTEGER(kind=kint) :: NEIBPETOT
    
    !C WS(dim), WR(dim)
    REAL   (kind=kreal), intent(out) :: WS(1:dim), WR(1:dim)
    
    INTEGER(kind=kint)  :: cnt, isL, ieL, isU, ieU, i, j, inod
    REAL   (kind=kreal) :: R_val, R_norm, GR_norm, v





    NEIBPETOT    =  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPETOT 
    NEIBPE       => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPE
    STACK_IMPORT => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT
    STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT
    NOD_IMPORT   => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_IMPORT
    NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_EXPORT
    
    DO cnt = 1, count
       call count_time(1, solver_comm, my_rank, 10)

       DO i = 1, N
          v = B(i)
          
          isL = INL(i - 1) + 1
          ieL = INL(i)
          DO j = isL, ieL
             inod = IAL(j) + ZERO_ORIGIN
             v = v - AL(j) * X(inod)
          END DO

          isU = INU(i - 1) + 1
          ieU = INU(i)
          DO j = isU, ieU
             inod = IAU(j) + ZERO_ORIGIN
             v = v - AU(j) * X(inod)
          END DO

#ifdef DEBUG
          if(D(i) == 0) then
             write(*,*) my_rank, i,"diagonal part equals 0!"
             stop "0 in diagonal part in sgs()"
          end if
#endif
          
          X(i) = v / D(i)
       END DO
       call count_time(2, solver_comm, my_rank, 10)       

       IF(NEIBPETOT > 0) then
          CALL SOLVER_SEND_RECV2(NP, STACK_EXPORT(NEIBPETOT), STACK_IMPORT(NEIBPETOT),  &
               &  NEIBPETOT, NEIBPE(1:NEIBPETOT), STACK_IMPORT(0:NEIBPETOT), &
               &  NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)), STACK_EXPORT(0:NEIBPETOT), &
               &  NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)), &
               &  WS(1:dim), WR(1:dim), X(1:NP), SOLVER_COMM, my_rank)
       end IF

       call count_time(1, solver_comm, my_rank, 10)       
       DO i = N, 1, -1
          v = B(i)

          isL = INL(i - 1) + 1
          ieL = INL(i)
          DO j = isL, ieL
             inod = IAL(j) + ZERO_ORIGIN
             v = v - AL(j) * X(inod)
          END DO
          
          isU = INU(i - 1) + 1
          ieU = INU(i)
          DO j = isU, ieU
             inod = IAU(j) + ZERO_ORIGIN
             v = v - AU(j) * X(inod)
          END DO

#ifdef DEBUG
          if(D(i) == 0) then
             write(*,*) my_rank, i,"diagonal part equals 0!"
             stop "0 in diagonal part in sgs()"
          end if
#endif

          X(i) = v / D(i)
       END DO
       CALL count_time(2, solver_comm, my_rank, 10)       
       
       if(NEIBPETOT > 0) then
          CALL SOLVER_SEND_RECV2(NP, STACK_EXPORT(NEIBPETOT), STACK_IMPORT(NEIBPETOT),  &
               &  NEIBPETOT, NEIBPE(1:NEIBPETOT), STACK_IMPORT(0:NEIBPETOT), &
               &  NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)), STACK_EXPORT(0:NEIBPETOT), &
               &  NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)), &
               &  WS(1:dim), WR(1:dim), X(1:NP), SOLVER_COMM, my_rank)
       end if
       
    END DO

  END SUBROUTINE sgs

  subroutine matrix_counting(LEVEL_NUM, SOLVER_COMM, my_rank)
    USE data_structure_for_AMG
    USE hash_mod
    IMPLICIT NONE
    
    include 'mpif.h'
    integer(kind=kint), intent(in) :: LEVEL_NUM, SOLVER_COMM, my_rank
    integer(kind=kint) :: i, j, nonzero_l, nonzero_g, level, n, ierr, np
    REAL(kind=kreal), allocatable :: comm_rate(:)
    REAL(kind=kreal) :: rate, crate_l, crate_g
    
    allocate(comm_rate(LEVEL_NUM))
    DO level = 1, LEVEL_NUM
       n  = HIERARCHICAL_DATA(level) % N
       np = HIERARCHICAL_DATA(level) % NP
       !C in the case of n == 0
       if(np == 0) then
          crate_l = 1
       else
          crate_l = REAL(n) / np
       end if
       
       CALL MPI_REDUCE(crate_l, crate_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, SOLVER_COMM, ierr)       
       comm_rate(level) = crate_g
    END DO
    
    nonzero_l=0
    
#ifndef CHOLESCKY

    DO level = 1, LEVEL_NUM
       n = HIERARCHICAL_DATA(level) % N
       if(n > 0) then
          nonzero_l = nonzero_l + n + HIERARCHICAL_DATA(level) % INU(n)
          nonzero_l = nonzero_l + HIERARCHICAL_DATA(level) % INL(n)
       end if
    END DO
       
#endif

#ifdef CHOLESCKY
    DO level = 1, LEVEL_NUM - 1
       n = HIERARCHICAL_DATA(level) % N
       nonzero_l = nonzero_l + n + HIERARCHICAL_DATA(level) % INU(n)
       nonzero_l = nonzero_l + HIERARCHICAL_DATA(level) % INL(n)
    END DO
#endif

    nonzero_g = 0
    CALL MPI_REDUCE(nonzero_l, nonzero_g, 1, LIS_MPI_INTEGER, MPI_SUM, 0, SOLVER_COMM, ierr)
    i = nonzero_g

#ifdef CHOLESCKY
    n = HIERARCHICAL_DATA(LEVEL_NUM) % N
    i = HIERARCHICAL_DATA(LEVEL_NUM) % NPL
    nonzero_g = nonzero_g + i * 2 - n
#endif


    n = HIERARCHICAL_DATA(1) % N
    
    !C in the case of n == 0
    nonzero_l = 0 

    if(n > 0) then
       nonzero_l = n + HIERARCHICAL_DATA(1) % INU(n) + HIERARCHICAL_DATA(1) % INL(n)
    end if

    CALL MPI_REDUCE(nonzero_l, nonzero_g, 1, LIS_MPI_INTEGER, MPI_SUM, 0, SOLVER_COMM, ierr)
    j = nonzero_g
    if(my_rank == 0)then
       rate = REAL(i, kind = kreal) / j
       write(*,*) "nonzeros in all level:", i, "rate of level:", rate
       write(*,*) "worst rate of own nodes of vector"
       write(*,*) comm_rate
    end if

  end subroutine matrix_counting

  !C Symmetrize the matrix especially the rows: N+1..NP
  SUBROUTINE matrix_arrange (N,NP,D,AL,INL,IAL,AU,INU,IAU)

    IMPLICIT NONE
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif

    INTEGER(kind=kint ), intent(in)    :: N,NP
    REAL   (kind=kreal), intent(inout),DIMENSION(:) :: D
    REAl   (kind=kreal), intent(inout),DIMENSION(:) :: AU,AL
    INTEGER(kind=kint ), intent(in),DIMENSION(:) :: IAU,IAL
    INTEGER(kind=kint ), intent(in),DIMENSION(0:) :: INU,INL

    INTEGER(kind=kint ) :: i,j,isU,ieU,isL,ieL,k,column
    logical :: flag


    DO j=N+1,NP
       isL = INL(j-1)+1
       ieL = INL(j)
       DO i=isL,ieL
          column = IAL(i)+ZERO_ORIGIN
          if(column>N) then
             AL(i)=0
          else if(column>0 .and. column<=N) then
             isU = INU(column-1)+1
             ieU = INU(column)
             flag = .false.
             DO k=isU,ieU
                if(IAU(k)+ZERO_ORIGIN==j)then
                   AL(i)=AU(k)
                   flag= .true.
                   exit
                end if
             end DO
             if(.not.flag) then
                AL(i) =0.0
             end if
          end if
       END DO

       isU = INU(j-1)+1
       ieU = INU(j)
       DO i=isU,ieU
          AU(i)= 0.0
       END DO
    END DO
  END SUBROUTINE matrix_arrange


  !C
  !C*** CG
  !C
  SUBROUTINE AMGCG    (N, NP,  NPL, NPU, D,  AL, INL, IAL, AU, INU, IAU,  &
       & B, X, RESID, ITER, ERROR, my_rank, NEIBPETOT, NEIBPE,   &
       & STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, SOLVER_COMM, &
       & NPROCS)
    
    USE solver_SR2
    USE data_structure_for_AMG
    USE solver_Gnumbering
    USE count_time_mod
    USE data_creation_AMGCG
    IMPLICIT NONE
    include 'mpif.h'
    
    INTEGER(kind=kint ), intent(in)  :: N,NP,NPL,NPU
    REAL   (kind=kreal), INTENT(inout)::  RESID
    INTEGER(kind=kint ), INTENT(inout)::  ITER
    INTEGER(kind=kint ), INTENT(inout)::  ERROR
    INTEGER(kind=kint ), INTENT(in   )::  my_rank

    INTEGER(kind=kint ), INTENT(in)   :: SOLVER_COMM
    INTEGER(kind=kint ), INTENT(in)   :: NPROCS
    
    REAL   (kind=kreal) ::  D(:)
    REAL   (kind=kreal) ::  B(:)

    !C NP
    REAL   (kind=kreal) ::  X(:)
    REAL   (kind=kreal) ::  AU(:)
    REAL   (kind=kreal) ::  AL(:)

    INTEGER(kind=kint ) ::  INU(0:)
    INTEGER(kind=kint ) ::  IAU(:)
    INTEGER(kind=kint ) ::  INL(0:)
    INTEGER(kind=kint ) ::  IAL(:)
    
    INTEGER(kind=kint ) :: NEIBPE      (:)
    INTEGER(kind=kint ) :: STACK_IMPORT(0:)
    INTEGER(kind=kint ) :: NOD_IMPORT  (:)

    INTEGER(kind=kint ) :: STACK_EXPORT(0:)
    INTEGER(kind=kint ) :: NOD_EXPORT  (:)

    REAL   (kind=kreal), ALLOCATABLE :: WS(:), WR(:)
    REAL   (kind=kreal), ALLOCATABLE :: WW(:,:)
    REAL   (kind=kreal), ALLOCATABLE :: Temp(:)

    REAL   (kind=kreal), ALLOCATABLE, SAVE :: DD(:)
    REAL   (kind=kreal), ALLOCATABLE, SAVE :: SCALE(:)

    INTEGER(kind=kint ) :: P, Q, R, Z, MAXIT, IFLAG=0
    REAL   (kind=kreal) :: TOL, W, SS

    INTEGER(kind=kint)  :: LEVEL_NUM,cnt,NEIBPETOT,ISU,IEU,ISL,IEL,WSIZE,WSIZE2
    INTEGER(kind=kint)  :: INOD, IERR, i,j,k,l,m
    REAL   (kind=kreal) :: WVAL,BNRM20,BNRM2,RHO0,RHO,BETA,RHO1,C10,C1,ALPHA,DNRM20
    REAL   (kind=kreal) :: DNRM2,X1,X2, theta
    REAL   (kind=kreal) :: STARTTIME,ENDTIME,RTIME,STARTTIME2,ENDTIME2,RTIME2

    !C-- INIT. 
    ERROR= 0


    ALLOCATE (WW(NP,3))

    !C- decision of the allocation size for WS, WR
    WSIZE = 0
    if(NEIBPETOT > 0) then
       if( STACK_IMPORT(NEIBPETOT) > STACK_EXPORT(NEIBPETOT) ) then
          WSIZE = STACK_IMPORT(NEIBPETOT)
       else
          WSIZE = STACK_EXPORT(NEIBPETOT)
       end if
    end if
    if( WSIZE < NP ) WSIZE = NP
    
    ALLOCATE ( WS(WSIZE) )
    ALLOCATE ( WR(WSIZE) )
    ALLOCATE ( Temp(NP) )
    
    R = 1
    Z = 2
    Q = 2
    P = 3

    MAXIT  = ITER
    TOL   = RESID           



!!$    call  matrix_arrange (N,NP,D,AL,INL,IAL,AU,INU,IAU)
    if (NEIBPETOT > 0) then
       !C diagonal part is modified to have correct values in overlapped region.
       call SOLVER_SEND_RECV2                                                     &
            &   ( NP,STACK_EXPORT(NEIBPETOT),STACK_IMPORT(NEIBPETOT),NEIBPETOT,   &
            &     NEIBPE, STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, WS, &
            &     WR, D , SOLVER_COMM, my_rank)

    end if
    
    
    STARTTIME2=MPI_WTIME()    
    
    CALL data_creation_ssi_amg(N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,    &
         & NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, &
         & STACK_IMPORT(NEIBPETOT), STACK_EXPORT(NEIBPETOT), SOLVER_COMM,         &
         & LEVEL_NUM, WSIZE2,theta)
    
    if(WSIZE2 > WSIZE) then
       deallocate(WS,WR)
       allocate(WS(WSIZE2),WR(WSIZE2))
       WSIZE = WSIZE2
    END if

    ENDTIME2=MPI_WTIME()
    
    
    STARTTIME = MPI_WTIME()

    !C
    !C +-----------------------+
    !C | {r0}= {b} - [A]{xini} |
    !C +-----------------------+
    !C===
    !C-- BEGIN calculation
    do j= 1, N
       WVAL= B(j) - D(j) * X(j)
       isU= INU(j-1) + 1
       ieU= INU(j  ) 

       do i= isU, ieU
          inod= IAU(i)+ZERO_ORIGIN
          WVAL= WVAL - AU(i) * X(inod)
       enddo

       isL= INL(j-1) + 1
       ieL= INL(j  ) 
       do i= isL, ieL
          inod= IAL(i)+ZERO_ORIGIN
          WVAL= WVAL - AL(i) * X(inod)
       enddo
       WW(j,R)= WVAL
    enddo

    BNRM20 = 0.d0
    do i= 1, N
       BNRM20 = BNRM20 + B(i)**2
    enddo

    call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
         &                    MPI_SUM, SOLVER_COMM, ierr)

    if (BNRM2.eq.0.d0) BNRM2= 1.d0
    ITER = 0

!C===
    do iter= 1, MAXIT
       call MPI_BARRIER  (SOLVER_COMM, ierr)
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===
       
       do i = 1, NP
          WW(i,Z) = 0.0
       end do

       call count_time(1, SOLVER_COMM, my_rank, 9)


       CALL v_cycle_ssi_amg(WW(:,R), WW(:,Z), TEMP, LEVEL_NUM, SOLVER_COMM, WS, WR, NP, WSIZE)


       call count_time(2,SOLVER_COMM,my_rank,9)

       if(NEIBPETOT > 0) then
          call SOLVER_SEND_RECV2                                                      &
               &   ( NP, WSIZE, WSIZE, NEIBPETOT,                                     &
               &     NEIBPE, STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, WS,  &
               &     WR, WW(:,Z), SOLVER_COMM, my_rank)
       END if

!C===

!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
       RHO0= 0.0

       do i= 1, N
          RHO0= RHO0 + WW(i,R)*WW(i,Z)
       enddo


       call MPI_allREDUCE (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
       &                    MPI_SUM, SOLVER_COMM, ierr)

!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
       if ( ITER.eq.1 ) then
          do i= 1, N
             WW(i,P)= WW(i,Z)
          enddo
       else
          BETA= RHO / RHO1
          do i= 1, N
             WW(i,P)= WW(i,Z) + BETA*WW(i,P)
          enddo
       endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===        

!C
!C-- INTERFACE data EXCHANGE
       if(NEIBPETOT > 0) then
          call SOLVER_SEND_RECV2                                                     &
               &   ( NP, WSIZE, WSIZE, NEIBPETOT,                                    &
               &     NEIBPE, STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, WS, &
               &     WR, WW(:,P), SOLVER_COMM, my_rank)
       end if


!C
!C-- BEGIN calculation

       do j= 1, N
          WVAL= D(j) * WW(j,P)
          
          isU= INU(j-1) + 1
          ieU= INU(j  ) 
          
          do i= isU, ieU
             inod= IAU(i)+ZERO_ORIGIN
             WVAL= WVAL + AU(i) * WW(inod,P)
          enddo

          isL= INL(j-1) + 1
          ieL= INL(j  ) 
          
          do i= isL, ieL
             inod= IAL(i)+ZERO_ORIGIN
             WVAL= WVAL + AL(i) * WW(inod,P)
          enddo
          WW(j,Q)= WVAL
       enddo
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
       C10= 0.d0
       
       do i= 1, N
          C10= C10 + WW(i,P)*WW(i,Q)
       enddo

       call MPI_allREDUCE (C10, C1, 1, MPI_DOUBLE_PRECISION,             &
            &                    MPI_SUM, SOLVER_COMM, ierr)

       ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===
      DNRM20= 0.d0

      X1= 0.0d0
      X2= 0.0d0
      
      do i= 1, N
         X(i)  = X (i)   + ALPHA * WW(i,P)
        WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
      enddo

      DNRM20 = 0.0
      do i= 1, N
         DNRM20= DNRM20 + WW(i,R)**2
      enddo

      call MPI_allREDUCE (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, SOLVER_COMM, ierr)

      RESID= dsqrt(DNRM2/BNRM2)


#ifdef PRINT_REZ
      if (my_rank.eq.0) write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)
      if ( RESID.le.TOL   ) then
         if(my_rank==0)then
            ENDTIME=MPI_WTIME()
            RTIME=ENDTIME-STARTTIME
            RTIME2=ENDTIME2-STARTTIME2
            write (*, '(i5, 1pe16.6)') ITER, RESID
            write (*, '("*** datacreation & iteration", 1pe16.6," sec",&
                 &   1pe16.6," sec." )') RTIME2,RTIME
         end if
         exit
      end if
#endif

      if ( ITER .eq.MAXIT ) ERROR= -100
      
      RHO1 = RHO                                                             
   end do
!C===


#ifdef PRINT_INFO
   call matrix_counting(LEVEL_NUM,SOLVER_COMM,my_rank)
   call count_time(4,solver_comm,my_rank,1)
#endif
     
     call MPI_BARRIER  (SOLVER_COMM,ierr)
     
     DEALLOCATE (WW)
     DEALLOCATE (WR)
     DEALLOCATE (WS)
     DEALLOCATE (Temp)


     call clear_matrix_ssi_amg(LEVEL_NUM)


#ifdef PSPASES
     CALL PSPACEC(pspcommf,0)
#endif
     
   END SUBROUTINE AMGCG
 END MODULE  solver_AMGCG
