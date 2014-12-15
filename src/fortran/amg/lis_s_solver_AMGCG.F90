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

!C   ************************************************
!C   * MODULE solver_AMGCG
!C     CONTAINS
!C   * SUBROUTINE v_cycle
!C   * SUBROUTINE sgs
!C   * SUBROUTINE matrix_counting
!C   * SUBROUTINE matrix_arrange
!C   * SUBROUTINE AMGCG
!C   * SUBROUTINE clear_matrix
!C   ************************************************

MODULE solver_AMGCG
CONTAINS
  
  SUBROUTINE v_cycle(N, problem_B, problem_X, LEVEL_NUM, Temp)
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint),  INTENT(in)    :: LEVEL_NUM
    INTEGER(kind=kint)                 :: N, NP
    INTEGER(kind=kint)                 :: TMP_INT
    REAL   (kind=kreal), DIMENSION(:), pointer :: D
    !C DIMENSION(N)
    REAL   (kind=kreal), DIMENSION(N), target ::  problem_B
    !C DIMENSION(N)
    REAL   (kind=kreal), DIMENSION(N), target::  problem_X
    !C DIMENSION(N)

    REAL   (kind=kreal), DIMENSION(: ), pointer::  B
    !C DIMENSION(N)
    REAL   (kind=kreal), DIMENSION(: ), pointer::  X   
    !C DIMENSION(N)

    REAL   (kind=kreal), DIMENSION(:), pointer::  AU
    REAL   (kind=kreal), DIMENSION(:), pointer::  AL

    INTEGER(kind=kint ), DIMENSION(:), pointer ::  INU
    !C DIMENSION(0:N)
    INTEGER(kind=kint ), DIMENSION(:), pointer ::  IAU
    !C DIMENSION(NPU)
    INTEGER(kind=kint ), DIMENSION(:), pointer ::  INL
    !C DIMENSION(0:N)
    INTEGER(kind=kint ), DIMENSION(:), pointer ::  IAL

    REAL   (kind=kreal), DIMENSION(:), pointer :: coarser_B
    REAL   (kind=kreal), DIMENSION(:), pointer :: coarser_X
    REAL   (kind=kreal), DIMENSION(N) :: Temp

    INTEGER(kind=kint ) :: nth_lev,i,inod,j,ieL,isL,is,ie,coarser_N,isU,ieU,k,coarser_NP
    INTEGER(kind=kint ) :: coarser_NEIBPETOT
    REAL   (kind=kreal) :: R_val,B_val,w,X_val,R_norm,GR_norm


    type(INTER_LEVEL_OPERATOR),pointer :: R, P



    HIERARCHICAL_DATA(1) % B =>problem_B
    HIERARCHICAL_DATA(1) % X =>problem_X
    
    do i = 1, N
       problem_X(i) = 0.0D0
    END do

    DO nth_lev=1,LEVEL_NUM-1
       N= HIERARCHICAL_DATA(nth_lev) % N
       NP=HIERARCHICAL_DATA(nth_lev) % NP
       INU=>HIERARCHICAL_DATA(nth_lev) % INU 
       INL=>HIERARCHICAL_DATA(nth_lev) % INL 
       X => HIERARCHICAL_DATA(nth_lev) % X
       B => HIERARCHICAL_DATA(nth_lev) % B
       D => HIERARCHICAL_DATA(nth_lev) % D
       IAU=>HIERARCHICAL_DATA(nth_lev) % IAU
       IAL=>HIERARCHICAL_DATA(nth_lev) % IAL
       AU =>HIERARCHICAL_DATA(nth_lev) % AU
       AL =>HIERARCHICAL_DATA(nth_lev) % AL
       
       coarser_B => HIERARCHICAL_DATA(nth_lev+1) % B
       coarser_X => HIERARCHICAL_DATA(nth_lev+1) % X

       TMP_INT=1
       CALL sgs(N,NP,D,AL,INL,IAL,AU,INU,IAU,B,X,TMP_INT,nth_lev)

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

       !C restrict the residual vector
       R => HIERARCHICAL_DATA(nth_lev+1) % R
       coarser_N= HIERARCHICAL_DATA(nth_lev+1) % N
       coarser_NP=HIERARCHICAL_DATA(nth_lev+1) % NP

       DO j= 1, coarser_NP
          is= R % IN(j-1)+1
          ie= R % IN(j)
          B_val=0.0
          DO i=is, ie
             inod = R % CN(i)
             B_val= B_val + R % V(i) * Temp(inod)
          END DO
          coarser_B(j)= B_val
       END DO

       DO j=1,coarser_NP
          coarser_X(j)= 0.0
       END DO
    END DO

    N=HIERARCHICAL_DATA(LEVEL_NUM) % N
    NP=HIERARCHICAL_DATA(LEVEL_NUM) % NP

    X => HIERARCHICAL_DATA(LEVEL_NUM) % X
    B => HIERARCHICAL_DATA(LEVEL_NUM) % B
    AL=> HIERARCHICAL_DATA(LEVEL_NUM) % AL
    INL=>HIERARCHICAL_DATA(LEVEL_NUM) % INL 

    INU => HIERARCHICAL_DATA(LEVEL_NUM) % INU 
    D   => HIERARCHICAL_DATA(LEVEL_NUM) % D
    IAU => HIERARCHICAL_DATA(LEVEL_NUM) % IAU
    IAL => HIERARCHICAL_DATA(LEVEL_NUM) % IAL
    AU  => HIERARCHICAL_DATA(LEVEL_NUM) % AU
    TMP_INT=30
    CALL sgs(N,NP,D,AL,INL,IAL,AU,INU,IAU,B,X,TMP_INT,LEVEL_NUM)


    DO nth_lev=LEVEL_NUM-1,1,-1
       N= HIERARCHICAL_DATA(nth_lev) % N
       NP=HIERARCHICAL_DATA(nth_lev) % NP

       INU          => HIERARCHICAL_DATA(nth_lev) % INU
       INL          => HIERARCHICAL_DATA(nth_lev) % INL
       X            => HIERARCHICAL_DATA(nth_lev) % X
       B            => HIERARCHICAL_DATA(nth_lev) % B
       D            => HIERARCHICAL_DATA(nth_lev) % D
       IAU          => HIERARCHICAL_DATA(nth_lev) % IAU
       IAL          => HIERARCHICAL_DATA(nth_lev) % IAL
       AU           => HIERARCHICAL_DATA(nth_lev) % AU
       AL           => HIERARCHICAL_DATA(nth_lev) % AL
       coarser_X    => HIERARCHICAL_DATA(nth_lev+1) % X
       P => HIERARCHICAL_DATA(nth_lev+1) % P


       DO j= 1, N
          is= P % IN(j-1)+1
          ie= P % IN(j)
          X_val=0.0
          DO i=is, ie
             inod = P % CN(i)
             X_val= X_val+P % V(i) * coarser_X(inod)
          END DO
          X(j)=X(j)+X_val
       END DO

       TMP_INT=1
       CALL sgs(N,NP,D,AL,INL,IAL,AU,INU,IAU,B,X,TMP_INT,nth_lev)
    END DO
  END SUBROUTINE v_cycle


  SUBROUTINE sgs(N, NP, D, AL, INL, IAL, AU, INU, IAU, B, X, count, LEVEL_NO)

    USE data_structure_for_AMG


    IMPLICIT NONE

    INTEGER(kind=kint),INTENT(in )  :: N,NP,count

    REAL   (kind=kreal), DIMENSION(:)   ::  D
    REAL   (kind=kreal), DIMENSION(:)   ::  B
    REAL   (kind=kreal), DIMENSION(:)   ::  X

    REAL   (kind=kreal), DIMENSION(:)   ::  AU
    REAL   (kind=kreal), DIMENSION(:)   ::  AL

    INTEGER(kind=kint ), DIMENSION(0:)  ::  INU
    INTEGER(kind=kint ), DIMENSION(:)   ::  IAU
    INTEGER(kind=kint ), DIMENSION(0:)  ::  INL
    INTEGER(kind=kint ), DIMENSION(:)   ::  IAL

    INTEGER(kind=kint), intent(in) :: LEVEL_NO

!!$    INTEGER(kind=kint), DIMENSION(:),pointer :: NEIBPE
!!$    INTEGER(kind=kint), DIMENSION(:),pointer :: STACK_IMPORT,STACK_EXPORT
!!$    INTEGER(kind=kint), DIMENSION(:),pointer :: NOD_EXPORT,NOD_IMPORT
!!$    INTEGER(kind=kint)                       :: NEIBPETOT
!!$



    INTEGER(kind=kint)  :: cnt,isL,ieL,isU,ieU,i,j,inod
    REAL   (kind=kreal) :: R_val,R_norm,GR_norm,v

!!$    NEIBPETOT    =  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPETOT 
!!$    NEIBPE       => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPE
!!$    STACK_IMPORT => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT
!!$    STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT
!!$    NOD_IMPORT   => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_IMPORT
!!$    NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_EXPORT


    DO cnt=1,count

       DO i=1,N
          v = B(i)

          isL=INL(i-1)+1
          ieL=INL(i)
          DO j=isL,ieL
             inod=IAL(j)+ZERO_ORIGIN
             v = v - AL(j)*X(inod)
          END DO

          isU=INU(i-1)+1
          ieU=INU(i)
          DO j=isU,ieU
             inod=IAU(j)+ZERO_ORIGIN
             v = v - AU(j)*X(inod)
          END DO

          X(i) = v / D(i)
       END DO
       DO i=N,1,-1
          v = B(i)

          isL=INL(i-1)+1
          ieL=INL(i)
          DO j=isL,ieL
             inod=IAL(j)+ZERO_ORIGIN
             v = v - AL(j)*X(inod)
          END DO

          isU=INU(i-1)+1
          ieU=INU(i)
          DO j=isU,ieU
             inod=IAU(j)+ZERO_ORIGIN
             v = v - AU(j)*X(inod)
          END DO

          X(i) = v / D(i)
       END DO
    END DO
  END SUBROUTINE sgs

  SUBROUTINE matrix_counting(LEVEL_NUM,SOLVER_COMM,my_rank)
    USE data_structure_for_AMG
    IMPLICIT NONE

    integer(kind=kint) ,intent(in) :: LEVEL_NUM,SOLVER_COMM,my_rank
    integer(kind=kint) :: i,j,nonzero_l,nonzero_g,level,n,np
    REAL,DIMENSION(:),allocatable::comm_rate
    REAL:: rate,crate_l,crate_g

    allocate(comm_rate(LEVEL_NUM))
    DO level=1,LEVEL_NUM
       n= HIERARCHICAL_DATA(level) % N
       np=HIERARCHICAL_DATA(level) % NP
       !C in the case of n == 0
       if(np == 0) then
          crate_l=1
       else
          crate_l=REAL(n)/np
       end if
       crate_g=0

       crate_g=crate_l
       comm_rate(level) = crate_g
    END DO

    nonzero_l=0

    DO level=1,LEVEL_NUM
       n= HIERARCHICAL_DATA(level) % N
       nonzero_l=nonzero_l+1 + HIERARCHICAL_DATA(level)%INU(n)
       nonzero_l=nonzero_l + HIERARCHICAL_DATA(level)%INL(n)
    END DO

    nonzero_g=0

    nonzero_g=nonzero_l
    i=nonzero_g


    n=HIERARCHICAL_DATA(1)%N

    !C in the case of n == 0
    nonzero_l = 0 
    if(n>0) then
       nonzero_l=1+HIERARCHICAL_DATA(1)%INU(n)+HIERARCHICAL_DATA(1)%INL(n)
    end if


    nonzero_g = 0

    nonzero_g=nonzero_l
    j=nonzero_g

    if(my_rank==0)then
       rate=REAL(i)/j
       write(*,*) "nonzeros in all level:", i, "rate of level:",rate
       write(*,*) "worst rate of own nodes of vector"
       write(*,*) comm_rate
    end if
  END SUBROUTINE matrix_counting

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
  SUBROUTINE AMGCG    (N, NP,  NPL, NPU,                                  &
       &                  D,  AL, INL, IAL, AU, INU, IAU,                 &
       &                  B,  X, resid,iter)     

    USE data_structure_for_AMG
    USE data_creation_AMGCG
    IMPLICIT NONE
    INTEGER(kind=kint ), intent(in)  :: N,NP,NPL,NPU
    REAL   (kind=kreal)              ::  RESID
    REAL   (kind=kreal)              ::  SIGMA_DIAG
    REAL   (kind=kreal)              ::  SIGMA
    INTEGER(kind=kint )              ::  ITER
    INTEGER(kind=kint )              ::  ERROR
    INTEGER(kind=kint )              ::  my_rank

    INTEGER(kind=kint )              :: NPROCS,SOLVER_COMM

    REAL   (kind=kreal), DIMENSION(: ),  pointer ::  D
    REAL   (kind=kreal), DIMENSION(: ),  pointer ::  B
    REAL   (kind=kreal), DIMENSION(NP )  ::  X

    REAL   (kind=kreal), DIMENSION(:), pointer ::  AU
    REAL   (kind=kreal), DIMENSION(:), pointer ::  AL

    INTEGER(kind=kint ), DIMENSION(:), pointer ::  INU
    INTEGER(kind=kint ), DIMENSION(:), pointer ::  IAU
    INTEGER(kind=kint ), DIMENSION(:), pointer ::  INL
    INTEGER(kind=kint ), DIMENSION(:), pointer ::  IAL


    REAL   (kind=kreal), DIMENSION(:,:), ALLOCATABLE :: WW
    REAL   (kind=kreal), DIMENSION(:),   ALLOCATABLE :: Temp

    REAL   (kind=kreal), DIMENSION(:), ALLOCATABLE, SAVE :: DD
    REAL   (kind=kreal), DIMENSION(:), ALLOCATABLE, SAVE :: SCALE

    INTEGER(kind=kint ) :: P, Q, R, Z, MAXIT, IFLAG=0
    REAL   (kind=kreal) :: TOL, W, SS

    INTEGER(kind=kint)  :: LEVEL_NUM,cnt,NEIBPETOT,I,J,ISU,IEU,ISL,IEL,WSIZE
    INTEGER(kind=kint)  :: INOD
    REAL   (kind=kreal) :: WVAL,BNRM20,BNRM2,RHO0,RHO,BETA,RHO1,C10,C1,ALPHA,DNRM20
    REAL   (kind=kreal) :: DNRM2,X1,X2
    REAL   (kind=kreal) :: STARTTIME,ENDTIME,RTIME,STARTTIME2,ENDTIME2,RTIME2

    REAL   (kind=kreal) :: theta
    !C-- theta means the critierion for strong connection.
    !C-- if theta=0 then all nonzero elements are recognized as 
    !C-- strong connection.

    !C-- INIT. 
    ERROR= 0
    NPROCS=0
    my_rank=0

    ALLOCATE (WW(NP,3))


    ALLOCATE ( Temp(NP) )

    R = 1
    Z = 2
    Q = 2
    P = 3

    MAXIT  = ITER
    TOL   = RESID           

    CALL data_creation(NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, LEVEL_NUM, theta)
    

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

    BNRM2=BNRM20

    if (BNRM2.eq.0.d0) BNRM2= 1.d0
    ITER = 0

    !C===
    do iter= 1, MAXIT
       !C
       !C************************************************* Conjugate Gradient Iteration

       !C
       !C +----------------+
       !C | {z}= [Minv]{r} |
       !C +----------------+
       !C===
       do i=1,NP
          WW(i,Z)=0.0
       end do
       CALL v_cycle(N, WW(:,R), WW(:,Z), LEVEL_NUM, Temp)
!!$       WW(:,Z)     =  WW(:,R)
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


       RHO=RHO0

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

       C1=C10

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


       DNRM2=DNRM20
       RESID= dsqrt(DNRM2/BNRM2)


!!$#ifdef PRINT_REZ
       if (my_rank.eq.0) write (*, 1000) ITER, RESID
1000   format (i5, 1pe16.6)
!!$#endif
       if ( RESID.le.TOL   ) then
          exit
       end if
       if ( ITER .eq.MAXIT ) ERROR= -100

       RHO1 = RHO                                                             
    end do
    !C===
    
    
    call clear_matrix(LEVEL_NUM)
    
    
    !C count non zeros of hierarchy of matrixes
    DEALLOCATE (WW)
    DEALLOCATE (Temp)

  END SUBROUTINE        AMGCG

  subroutine clear_matrix(LEVEL_NUM)
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
            &      HIERARCHICAL_DATA(1) % B         )

    DO i=2,LEVEL_NUM
       deallocate(HIERARCHICAL_DATA(i) % INU,      &
            &     HIERARCHICAL_DATA(i) % INL,      &
            &     HIERARCHICAL_DATA(i) % IAU,      &
            &     HIERARCHICAL_DATA(i) % IAL,      &
            &     HIERARCHICAL_DATA(i) % AU,       &
            &     HIERARCHICAL_DATA(i) % AL,       &
            &     HIERARCHICAL_DATA(i) % D,        &
            &     HIERARCHICAL_DATA(i) % X,        &
            &     HIERARCHICAL_DATA(i) % B,        &
            &     HIERARCHICAL_DATA(i) % R % IN,   &
            &     HIERARCHICAL_DATA(i) % R % CN,   &
            &     HIERARCHICAL_DATA(i) % R % V,    &
            &     HIERARCHICAL_DATA(i) % P % IN,   &
            &     HIERARCHICAL_DATA(i) % P % CN,   &
            &     HIERARCHICAL_DATA(i) % P % V     )
       deallocate(HIERARCHICAL_DATA(i) % P,        &
            &     HIERARCHICAL_DATA(i) % R         )
    END DO


  end subroutine clear_matrix

END MODULE  solver_AMGCG


