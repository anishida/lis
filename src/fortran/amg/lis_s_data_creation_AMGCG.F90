!C   Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.
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
!C   * MODULE data_creation_AMGCG
!C     CONTAINS
!C   * SUBROUTINE data_creation_unsym
!C   * SUBROUTINE data_creation
!C   * SUBROUTINE RAP2
!C   * SUBROUTINE make_rap
!C   * SUBROUTINE smooth_aggregate
!C   * SUBROUTINE neighbors
!C   * SUBROUTINE neighbors_unsym
!C   ************************************************

MODULE data_creation_AMGCG
CONTAINS
  SUBROUTINE data_creation_unsym(NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, LEVEL_NUM, theta)
    USE isort
    USE data_structure_for_AMG
    USE aggregate_mod

    IMPLICIT NONE
    INTEGER(kind=kint), INTENT(in)  :: NP, NPL, NPU
    INTEGER(kind=kint), INTENT(out) :: LEVEL_NUM
    REAL   (kind=kreal),INTENT(inout),TARGET :: D(NP )
    REAL   (kind=kreal),INTENT(inout),TARGET :: AU(NPU)
    REAL   (kind=kreal),INTENT(inout),TARGET :: AL(NPL)
    INTEGER(kind=kint ),INTENT(in),   TARGET :: IAU(NPU)
    INTEGER(kind=kint ),INTENT(in),   TARGET :: IAL(NPL)
    INTEGER(kind=kint ),INTENT(in),   TARGET :: INU(0:NP), INL(0:NP)

    INTEGER(kind=kint) :: NPROCS

    REAL   (kind=kreal), ALLOCATABLE :: tilde_D(:)
    INTEGER(kind=kint ), ALLOCATABLE :: NI(:),PNI(:)
    INTEGER(kind=kint ), ALLOCATABLE :: node_index(:)

    INTEGER(kind=kint)  :: coarser_level_size,level,finer_level_size,count
    INTEGER(kind=kint)  :: finer_level_NP,i,j,k
    INTEGER(kind=kint)  :: finer_level_NPL,finer_level_NPU
    
    INTEGER(kind=kint), POINTER ::in_aggregates_result(:)
    INTEGER(kind=kint), POINTER :: aggregates_result(:)
    INTEGER(kind=kint) :: in_aggregates_result_size

    TYPE(row_node),     POINTER :: Temp_P(:,:)
    INTEGER(kind=kint), POINTER :: Temp_P_row_size(:)

    LOGICAL            :: finish_flag, global_finish_flag
    
    INTEGER(kind=kint) :: local_aggre_size,col
    TYPE(INTER_LEVEL_OPERATOR),POINTER :: P,R

    REAL   (kind=kreal) :: theta
    !C-- theta means the critierion for strong connection.
    !C-- if theta=0 then all nonzero elements are recognized as 
    !C-- strong connection.

    INTEGER(kind=kint ) :: RATE_OF_SPACE
    INTEGER(kind=kint ) :: TEMP_COLSIZE, NONZEROS_PER_ROW
    INTEGER(kind=kint ) :: TEMP_RAP_WORKSIZE
    !C SYM_FLAG == 1: symmetric mat, 0: unsymmetric mat
    INTEGER(kind=kint) :: SYM_FLAG = 1
    NPROCS=1

    HIERARCHICAL_DATA(1) % N   = NP
    HIERARCHICAL_DATA(1) % NP  = NP
    HIERARCHICAL_DATA(1) % NPL = NPL
    HIERARCHICAL_DATA(1) % NPU = NPU
    HIERARCHICAL_DATA(1) % D   => D
    HIERARCHICAL_DATA(1) % INL => INL 
    HIERARCHICAL_DATA(1) % INU => INU 
    HIERARCHICAL_DATA(1) % IAL => IAL
    HIERARCHICAL_DATA(1) % IAU => IAU
    HIERARCHICAL_DATA(1) % AL  => AL
    HIERARCHICAL_DATA(1) % AU  => AU

    !C-- upper bound of the number of nonzeros per row 
    !C-- on next coarser level
    !C-- The allocated memory size is 
    !C-- unknowns*12*NONZEROS_PER_ROW byte.
    NONZEROS_PER_ROW=100

    !C-- upper bound of the number of nonzeros per row 
    !C-- in prolongation matrix. This value is smaller 
    !C-- than NONZEROS_PER_ROW
    !C-- The allocated memory size is 
    !C-- unknowns*12*RATE_OF_SPACE byte.
    RATE_OF_SPACE=30
    
    !C-- work area of the calculation of RAP
    !C-- this value must be bigger than or equal to 1.
    !C-- The allocated memory size is 
    !C-- unknowns*8*TEMP_RAP_WORKSIZE byte.
    TEMP_RAP_WORKSIZE=40
    
    
    
    
    DO level = 2,MAX_LEVEL_SIZE
       
       finer_level_size = HIERARCHICAL_DATA(level - 1) % N
       finer_level_NP   = HIERARCHICAL_DATA(level - 1) % NP
       finer_level_NPL  = HIERARCHICAL_DATA(level - 1) % NPL
       finer_level_NPU  = HIERARCHICAL_DATA(level - 1) % NPU       


       ALLOCATE(NI(finer_level_NPL + finer_level_NPU))
       ALLOCATE(PNI(0:2*finer_level_NP))
       ALLOCATE(node_index(finer_level_NP))
       ALLOCATE(tilde_D   (finer_level_NP))

       if(SYM_FLAG == 1) then
          CALL neighbors(PNI, NI, theta, level, node_index, tilde_D)
       else
          CALL neighbors_unsym(PNI, NI, theta, level, node_index, tilde_D)
       end if

       CALL indpdt_agrgt(PNI, NI, level, node_index, in_aggregates_result,   &
            & aggregates_result, in_aggregates_result_size)

       DEALLOCATE(node_index)
       
       coarser_level_size = in_aggregates_result_size
       global_finish_flag = (coarser_level_size < MIN_NODE_SIZE)
       

       !C temporary work area for inter-level matrix P is allocated
       !C (RATE_OF_SPACE+0.5)*finer_level_NP*8 byte
       RATE_OF_SPACE = RATE_OF_SPACE + 80 * (level - 2)
       IF(RATE_OF_SPACE > coarser_level_size) RATE_OF_SPACE = coarser_level_size
       ALLOCATE(Temp_P(RATE_OF_SPACE, finer_level_NP))
       ALLOCATE(Temp_P_row_size(finer_level_NP))

       CALL smooth_aggregate(PNI, NI, REAL(0.6666,KIND=KREAL), level,      &
            & in_aggregates_result_size, in_aggregates_result,             &
            & aggregates_result, Temp_P, Temp_P_row_size, tilde_D,         &
            & RATE_OF_SPACE)

       CALL agrgt_dealloc (in_aggregates_result, aggregates_result)
       DEALLOCATE(NI)
       DEALLOCATE(PNI)
       DEALLOCATE(tilde_D)

       !C for single CPU case
       local_aggre_size = in_aggregates_result_size

       !--C check Temp_P
       DO i=1,finer_level_NP
          IF(Temp_P_row_size(i) > RATE_OF_SPACE) THEN
             WRITE(*,*) "RATE_OF_SPACE should be enlarged "
             STOP 'The allocation error : Temp_P in make_rap'
          END IF
       END DO

       !C-- Temp_P is recorded in HIERARCHICAL_DATA()..
       ALLOCATE(HIERARCHICAL_DATA(level) % P)
       P => HIERARCHICAL_DATA(level) % P
       P % ROW_SIZE = HIERARCHICAL_DATA(level - 1) % NP
       k = P % ROW_SIZE
       ALLOCATE(P % IN(0:k)) 
       P % IN(0) = 0
       DO i = 1, k
          P % IN(i) = P % IN(i - 1) + Temp_P_row_size(i)
       END DO
       i = P % IN(k)

       ALLOCATE(HIERARCHICAL_DATA(level) % P % CN(i))
       ALLOCATE(HIERARCHICAL_DATA(level) % P % V(i))

       k = 1
       DO i = 1, P % ROW_SIZE
          DO j = 1, Temp_P_row_size(i)
             P % CN(k) = Temp_P(j,i) % column
             P % V(k)  = Temp_P(j,i) % value
             k = k + 1
          END DO
       END DO

       DEALLOCATE(Temp_P)
       DEALLOCATE(Temp_P_row_size)
       !C-- end:Temp_P is recorded in HIERARCHICAL_DATA()..
       

       !C-- make R : 1..finer_level_NP -> 1..local_aggre_size
       ALLOCATE(HIERARCHICAL_DATA(level) % R)
       R => HIERARCHICAL_DATA(level) % R
       R % ROW_SIZE = local_aggre_size
       ALLOCATE(R % IN(0:local_aggre_size))
       R % IN(0:local_aggre_size) = 0

       count = 0
       DO i = 1,finer_level_NP
          DO j = P % IN(i - 1) + 1, P % IN(i)
             col = P % CN(j)
             count = count + 1
             IF(col < local_aggre_size) R % IN(col + 1) = R % IN(col + 1) + 1
          END DO
       END DO

       ALLOCATE(R % CN(count))
       ALLOCATE(R % V(count))

       DO i = 2, R % ROW_SIZE
          R % IN(i) = R % IN(i) + R % IN(i - 1)
       END DO

       DO i = 1, finer_level_NP
          DO j = P % IN(i - 1) + 1, P % IN(i)
             col = P % CN(j)
             k = R % IN(col) + 1
             R % IN(col) = k

             R % CN(k) = i
             R %  V(k) = P % V(j)
          END DO
       END DO
       !C--end: make R


       !C-- allocate temporaly domain and RAP().
       !C   TEMP_COLSIZE*coarser_level_size*12+coarser_level byte
       TEMP_COLSIZE = NONZEROS_PER_ROW + 600 * (level - 1)
       IF( TEMP_COLSIZE > coarser_level_size) TEMP_COLSIZE = coarser_level_size
       
       CALL make_rap(level, NPROCS, coarser_level_size, .FALSE.,    &
            &        local_aggre_size, TEMP_COLSIZE, TEMP_RAP_WORKSIZE, SYM_FLAG)
       
       IF(global_finish_flag)   EXIT

!!$#ifdef ANISO
!!$       theta = theta * 0.95
!!$#else
       theta = theta * 0.5       
!!$#endif
       
    END DO
    
    LEVEL_NUM = level

  END SUBROUTINE data_creation_unsym

  SUBROUTINE data_creation(NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, LEVEL_NUM, theta)
    USE isort
    USE data_structure_for_AMG
    USE aggregate_mod

    IMPLICIT NONE
    INTEGER(kind=kint), INTENT(in)  :: NP, NPL, NPU
    INTEGER(kind=kint), INTENT(out) :: LEVEL_NUM
    REAL   (kind=kreal),INTENT(inout),TARGET :: D(NP )
    REAL   (kind=kreal),INTENT(inout),TARGET :: AU(NPU)
    REAL   (kind=kreal),INTENT(inout),TARGET :: AL(NPL)
    INTEGER(kind=kint ),INTENT(in),   TARGET :: IAU(NPU)
    INTEGER(kind=kint ),INTENT(in),   TARGET :: IAL(NPL)
    INTEGER(kind=kint ),INTENT(in),   TARGET :: INU(0:NP), INL(0:NP)

    INTEGER(kind=kint) :: NPROCS

    REAL   (kind=kreal), ALLOCATABLE :: tilde_D(:)
    INTEGER(kind=kint ), ALLOCATABLE :: NI(:),PNI(:)
    INTEGER(kind=kint ), ALLOCATABLE :: node_index(:)

    INTEGER(kind=kint)  :: coarser_level_size,level,finer_level_size,count
    INTEGER(kind=kint)  :: finer_level_NP,i,j,k
    INTEGER(kind=kint)  :: finer_level_NPL,finer_level_NPU
    
    INTEGER(kind=kint), POINTER ::in_aggregates_result(:)
    INTEGER(kind=kint), POINTER :: aggregates_result(:)
    INTEGER(kind=kint) :: in_aggregates_result_size

    TYPE(row_node),     POINTER :: Temp_P(:,:)
    INTEGER(kind=kint), POINTER :: Temp_P_row_size(:)

    LOGICAL            :: finish_flag, global_finish_flag
    
    INTEGER(kind=kint) :: local_aggre_size,col
    TYPE(INTER_LEVEL_OPERATOR),POINTER :: P,R

    REAL   (kind=kreal) :: theta
    !C-- theta means the critierion for strong connection.
    !C-- if theta=0 then all nonzero elements are recognized as 
    !C-- strong connection.

    INTEGER(kind=kint ) :: RATE_OF_SPACE
    INTEGER(kind=kint ) :: TEMP_COLSIZE, NONZEROS_PER_ROW
    INTEGER(kind=kint ) :: TEMP_RAP_WORKSIZE
    !C SYM_FLAG == 1: symmetric mat, 0: unsymmetric mat
    INTEGER(kind=kint) :: SYM_FLAG = 0
    NPROCS=1

    HIERARCHICAL_DATA(1) % N   = NP
    HIERARCHICAL_DATA(1) % NP  = NP
    HIERARCHICAL_DATA(1) % NPL = NPL
    HIERARCHICAL_DATA(1) % NPU = NPU
    HIERARCHICAL_DATA(1) % D   => D
    HIERARCHICAL_DATA(1) % INL => INL 
    HIERARCHICAL_DATA(1) % INU => INU 
    HIERARCHICAL_DATA(1) % IAL => IAL
    HIERARCHICAL_DATA(1) % IAU => IAU
    HIERARCHICAL_DATA(1) % AL  => AL
    HIERARCHICAL_DATA(1) % AU  => AU

    !C-- upper bound of the number of nonzeros per row 
    !C-- on next coarser level
    !C-- The allocated memory size is 
    !C-- unknowns*12*NONZEROS_PER_ROW byte.
    NONZEROS_PER_ROW=100

    !C-- upper bound of the number of nonzeros per row 
    !C-- in prolongation matrix. This value is smaller 
    !C-- than NONZEROS_PER_ROW
    !C-- The allocated memory size is 
    !C-- unknowns*12*RATE_OF_SPACE byte.
    RATE_OF_SPACE=30
    
    !C-- work area of the calculation of RAP
    !C-- this value must be bigger than or equal to 1.
    !C-- The allocated memory size is 
    !C-- unknowns*8*TEMP_RAP_WORKSIZE byte.
    TEMP_RAP_WORKSIZE=40
    
!!$!!$#ifdef ANISO
!!$!!$    theta = 0.11
!!$!!$#else
!!$    theta = 0.05
!!$!!$#endif
    
    
    
    DO level = 2,MAX_LEVEL_SIZE
       
       finer_level_size = HIERARCHICAL_DATA(level - 1) % N
       finer_level_NP   = HIERARCHICAL_DATA(level - 1) % NP
       finer_level_NPL  = HIERARCHICAL_DATA(level - 1) % NPL
       finer_level_NPU  = HIERARCHICAL_DATA(level - 1) % NPU       


       ALLOCATE(NI(finer_level_NPL + finer_level_NPU))
       ALLOCATE(PNI(0:2*finer_level_NP))
       ALLOCATE(node_index(finer_level_NP))
       ALLOCATE(tilde_D   (finer_level_NP))

       if(SYM_FLAG == 1) then
          CALL neighbors(PNI, NI, theta, level, node_index, tilde_D)
       else
          CALL neighbors_unsym(PNI, NI, theta, level, node_index, tilde_D)
       end if

       CALL indpdt_agrgt(PNI, NI, level, node_index, in_aggregates_result,   &
            & aggregates_result, in_aggregates_result_size)

       DEALLOCATE(node_index)
       
       coarser_level_size = in_aggregates_result_size
       global_finish_flag = (coarser_level_size < MIN_NODE_SIZE)
       

       !C temporary work area for inter-level matrix P is allocated
       !C (RATE_OF_SPACE+0.5)*finer_level_NP*8 byte
       RATE_OF_SPACE = RATE_OF_SPACE + 80 * (level - 2)
       IF(RATE_OF_SPACE > coarser_level_size) RATE_OF_SPACE = coarser_level_size
       ALLOCATE(Temp_P(RATE_OF_SPACE, finer_level_NP))
       ALLOCATE(Temp_P_row_size(finer_level_NP))

       CALL smooth_aggregate(PNI, NI, REAL(0.6666,KIND=KREAL), level,      &
            & in_aggregates_result_size, in_aggregates_result,             &
            & aggregates_result, Temp_P, Temp_P_row_size, tilde_D,         &
            & RATE_OF_SPACE)

       CALL agrgt_dealloc (in_aggregates_result, aggregates_result)
       DEALLOCATE(NI)
       DEALLOCATE(PNI)
       DEALLOCATE(tilde_D)

       !C for single CPU case
       local_aggre_size = in_aggregates_result_size

       !--C check Temp_P
       DO i=1,finer_level_NP
          IF(Temp_P_row_size(i) > RATE_OF_SPACE) THEN
             WRITE(*,*) "RATE_OF_SPACE should be enlarged "
             STOP 'The allocation error : Temp_P in make_rap'
          END IF
       END DO

       !C-- Temp_P is recorded in HIERARCHICAL_DATA()..
       ALLOCATE(HIERARCHICAL_DATA(level) % P)
       P => HIERARCHICAL_DATA(level) % P
       P % ROW_SIZE = HIERARCHICAL_DATA(level - 1) % NP
       k = P % ROW_SIZE
       ALLOCATE(P % IN(0:k)) 
       P % IN(0) = 0
       DO i = 1, k
          P % IN(i) = P % IN(i - 1) + Temp_P_row_size(i)
       END DO
       i = P % IN(k)

       ALLOCATE(HIERARCHICAL_DATA(level) % P % CN(i))
       ALLOCATE(HIERARCHICAL_DATA(level) % P % V(i))

       k = 1
       DO i = 1, P % ROW_SIZE
          DO j = 1, Temp_P_row_size(i)
             P % CN(k) = Temp_P(j,i) % column
             P % V(k)  = Temp_P(j,i) % value
             k = k + 1
          END DO
       END DO

       DEALLOCATE(Temp_P)
       DEALLOCATE(Temp_P_row_size)
       !C-- end:Temp_P is recorded in HIERARCHICAL_DATA()..
       

       !C-- make R : 1..finer_level_NP -> 1..local_aggre_size
       ALLOCATE(HIERARCHICAL_DATA(level) % R)
       R => HIERARCHICAL_DATA(level) % R
       R % ROW_SIZE = local_aggre_size
       ALLOCATE(R % IN(0:local_aggre_size))
       R % IN(0:local_aggre_size) = 0

       count = 0
       DO i = 1,finer_level_NP
          DO j = P % IN(i - 1) + 1, P % IN(i)
             col = P % CN(j)
             count = count + 1
             IF(col < local_aggre_size) R % IN(col + 1) = R % IN(col + 1) + 1
          END DO
       END DO

       ALLOCATE(R % CN(count))
       ALLOCATE(R % V(count))

       DO i = 2, R % ROW_SIZE
          R % IN(i) = R % IN(i) + R % IN(i - 1)
       END DO

       DO i = 1, finer_level_NP
          DO j = P % IN(i - 1) + 1, P % IN(i)
             col = P % CN(j)
             k = R % IN(col) + 1
             R % IN(col) = k

             R % CN(k) = i
             R %  V(k) = P % V(j)
          END DO
       END DO
       !C--end: make R


       !C-- allocate temporaly domain and RAP().
       !C   TEMP_COLSIZE*coarser_level_size*12+coarser_level byte
       TEMP_COLSIZE = NONZEROS_PER_ROW + 600 * (level - 1)
       IF( TEMP_COLSIZE > coarser_level_size) TEMP_COLSIZE = coarser_level_size
       
       CALL make_rap(level, NPROCS, coarser_level_size, .FALSE.,    &
            &        local_aggre_size, TEMP_COLSIZE, TEMP_RAP_WORKSIZE, SYM_FLAG)
       
!C       WRITE(*,*) '############',level,'############'

       IF(global_finish_flag)   EXIT

!!$#ifdef ANISO
!!$       theta = theta * 0.95
!!$#else
       theta = theta * 0.5       
!!$#endif
       
    END DO
    
    LEVEL_NUM = level

  END SUBROUTINE data_creation

  
  SUBROUTINE RAP2(Temp_N, Temp_CN, Temp_V, hash_size, N, ownaggre_size, LEVEL_NO, &
       & local_aggre_size, Temp_SIZE, SYM_FLAG)

    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint), INTENT(in)   :: hash_size, N, ownaggre_size, LEVEL_NO, local_aggre_size
    INTEGER(kind=kint), INTENT(inout):: Temp_N(:)
    INTEGER(kind=kint), INTENT(inout):: Temp_CN(:,:)
    TYPE(INTER_LEVEL_OPERATOR), POINTER :: P, R
    REAL(kind=kreal), INTENT(inout) :: Temp_V(:,:)

    INTEGER(kind=kint), INTENT(in) :: Temp_SIZE, SYM_FLAG

    !C pointers for matrix
    REAL   (kind=kreal), POINTER ::  D(:), AU(:), AL(:)
    INTEGER(kind=kint ), POINTER ::  INU(:), IAU(:), INL(:), IAL(:)
    
    INTEGER(kind=kint ), ALLOCATABLE :: Temp_agr(:), Temp_sz(:)
    INTEGER(kind=kint ), ALLOCATABLE :: Temp_nd(:)
    INTEGER(kind=kint ), ALLOCATABLE :: Temp_int(:,:)
    REAL   (kind=kreal), ALLOCATABLE :: Temp(:,:), Temp_1_col(:)

    INTEGER(kind=kint) :: i, j, k, l, m, rowins, rowine, colins, coline, row, col, row2
    INTEGER(kind=kint) :: is, ie
    INTEGER(kind=kint) :: js, je, jin, rowin, colin, count, hash_no, NP
    LOGICAL :: cancel_flag, finish_flag, flag

    INTEGER(kind=kint ) :: in1, in2, P_col, T_col, p_width, t_width, cancel_sz, nd_sz
    REAL   (kind=kreal) :: v

    NP   =  HIERARCHICAL_DATA(LEVEL_NO - 1) % NP 
    D    => HIERARCHICAL_DATA(LEVEL_NO - 1) % D 
    INL  => HIERARCHICAL_DATA(LEVEL_NO - 1) % INL  
    INU  => HIERARCHICAL_DATA(LEVEL_NO - 1) % INU 
    IAL  => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAL
    IAU  => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAU
    AL   => HIERARCHICAL_DATA(LEVEL_NO - 1) % AL
    AU   => HIERARCHICAL_DATA(LEVEL_NO - 1) % AU
    P    => HIERARCHICAL_DATA(LEVEL_NO    ) % P
    R    => HIERARCHICAL_DATA(LEVEL_NO    ) % R

    !C Calculation of RAP. This uses R information. 
    !C R % IN must be rewrite @ the latter part.
    !C So Temp_R_IN is used for R % IN(:).
    !C ------
    ALLOCATE(Temp_sz (N), Temp_int(Temp_SIZE,N), Temp(Temp_SIZE,N),          &  
         &   Temp_agr(local_aggre_size), Temp_nd(local_aggre_size),          &
         &   Temp_1_col(local_aggre_size))
       
    Temp_1_col = 0    
    
    Temp_agr = 0
    !C-- RA
    finish_flag = .TRUE.
    DO WHILE(finish_flag)
       finish_flag = .FALSE.
       Temp_sz = 0
       
       DO col = 1, local_aggre_size
          IF(Temp_agr(col) > 0) CYCLE
          cancel_flag = .FALSE.
          
          DO i = R % IN(col - 1) + 1, R % IN(col) 
             row = R % CN(i)
             IF(row <= N .AND. Temp_sz(row) >= Temp_SIZE) THEN
                cancel_flag = .TRUE.
                EXIT
             END IF
          END DO
          IF(.NOT.cancel_flag) THEN
             DO i = R % IN(col - 1) + 1, R % IN(col) 
                row = R % CN(i)
                DO j = INL(row - 1) + 1, INL(row)
                   row2 = IAL( j )+ZERO_ORIGIN
                   IF(row2 <= N .AND. Temp_sz(row2) >= Temp_SIZE) THEN
                      cancel_flag = .TRUE.
                      EXIT
                   END IF
                END DO
                IF(cancel_flag) EXIT
                DO j = INU(row - 1) + 1, INU(row)
                   row2 = IAU( j )+ZERO_ORIGIN
                   IF(row2 <= N .AND. Temp_sz(row2) >= Temp_SIZE) THEN
                      cancel_flag = .TRUE.
                      EXIT
                   END IF
                END DO
                IF(cancel_flag) EXIT
             END DO
          END IF

          finish_flag = finish_flag .OR. cancel_flag
          IF(cancel_flag) CYCLE
          
          Temp_agr(col) = 1
          
          DO i = R % IN(col - 1) + 1, R % IN(col) 
             row = R % CN(i)
             v   = R % V(i)
             IF(row <= N) THEN
                k = Temp_sz( row )

                flag = (k==0)
                IF(k>0) THEN
                   IF(Temp_int(k,row) /= col) THEN
                      flag=.TRUE.
                   END IF
                END IF
                IF(flag) THEN                
!!$                IF( k == 0 .OR. Temp_int(k,row) /= col ) THEN                
                   k = k+1
                   Temp_sz ( row ) = k
                   Temp_int(k,row) = col
                   Temp    (k,row) = v * D(row)
                ELSE 
                   Temp(k,row) = Temp(k,row) + v * D(row)
                END IF
             END IF

             DO j = INL(row - 1) +1, INL(row)
                row2 = IAL( j )+ZERO_ORIGIN
                !C --
                IF(row2 > N) THEN
                   CYCLE
                END IF
                !C --

                k = Temp_sz(row2)

                flag = (k==0)
                IF(k > 0) THEN
                   IF(Temp_int(k,row2) /= col) THEN
                      flag = .TRUE.
                   END IF
                END IF
                IF(flag) THEN                
!!$                IF( k == 0 .OR. (k>0 .and. Temp_int(k,row2)/=col) ) THEN
                   k = k + 1
                   Temp_sz(row2) = k
                   Temp_int(k, row2) = col
                   Temp    (k, row2) = v * AL(j)
                ELSE
                   Temp(k, row2) = Temp(k, row2) + v * AL(j)
                END IF
             END DO

             DO j = INU(row - 1) + 1, INU(row)
                row2 = IAU(j)+ZERO_ORIGIN
                !C --
                IF(row2 > N) CYCLE
                !C --
                   
                k = Temp_sz(row2)

                flag = (k==0)
                IF(k > 0) THEN
                   IF(Temp_int(k, row2) /= col) THEN
                      flag = .TRUE.
                   END IF
                END IF
                IF(flag) THEN                
!!$                IF( k==0 .OR. (k>0 .and. Temp_int(k,row2)/=col) ) THEN
                   k = k + 1
                   Temp_sz(row2) = k
                   Temp_int(k, row2) = col
                   Temp    (k, row2) = v * AU(j)
                ELSE
                   Temp(k, row2) = Temp(k, row2) + v * AU(j)
                END IF
             END DO
          END DO
       END DO
       !C-- end: RA


       !C-- temp * P
       DO col = 1, R % ROW_SIZE
          nd_sz = 0
              
          js = R % IN(col - 1) + 1
          je = R % IN(col)
          DO j = js, je
             k = R % CN(j)
             v = R % V(j)

             DO l = 1, Temp_sz(k)
                row = Temp_int(l, k)

                IF(SYM_FLAG /= 1 .OR. col >= row .OR. row > ownaggre_size) THEN
                   
                   DO m = 1, nd_sz
                      IF(Temp_nd(m) == row) EXIT
                   END DO
                   IF(m > nd_sz) THEN
                      nd_sz = nd_sz + 1
                      Temp_nd(nd_sz) = row
                   END IF
                   Temp_1_col(row) = Temp_1_col(row) + v * Temp(l, k)
                END IF
             END DO
          END DO

          !C-- input
          DO j = 1, nd_sz
             row = Temp_nd(j)
             hash_no = MOD(col, hash_size)+1
             DO count = hash_no, hash_size
                IF(Temp_CN(count, row) == 0) EXIT
             END DO
             IF(count > hash_size) THEN
                DO count = 1, hash_no - 1
                   IF(Temp_CN(count,row) == 0) EXIT
                END DO
                IF(count >= hash_no) THEN
                   WRITE(*,*) temp_nd(1:nd_sz)
                   WRITE(*,*) count,row, hash_no,hash_size,Temp_CN(count,row)
                   STOP "hash overflow : RAP2()"
                END IF
             END IF
             
             !C assign the value
             Temp_N(row) = Temp_N(row) + 1
             Temp_CN(count,row) = col
             Temp_V(count,row)  = Temp_1_col(row)

             Temp_1_col(row) = 0
          END DO
          !C-- end: input
          
       END DO
       !C-- end:temp * P
    END DO
    
    DEALLOCATE(Temp_nd, Temp_sz, Temp_int, Temp, Temp_agr, Temp_1_col)
    
  END SUBROUTINE RAP2

  SUBROUTINE make_rap(LEVEL_NO, NPROCS, coarser_level_size,        &
       &              global_finish_flag, local_aggre_size, TEMP_COLSIZE,   &
       &              TEMP_RAP_WORKSIZE, SYM_FLAG)
    
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint), INTENT(in) :: LEVEL_NO
    LOGICAL,            INTENT(in) :: global_finish_flag
    INTEGER(kind=kint), INTENT(in) :: coarser_level_size
    INTEGER(kind=kint), INTENT(in) :: SYM_FLAG
    
    TYPE(INTER_LEVEL_OPERATOR), POINTER :: R, P

    !C temporary space of the result of RAP
    INTEGER(kind=kint ), INTENT(in)                 :: Temp_COLSIZE
    REAL   (kind=kreal), ALLOCATABLE :: Temp_V(:,:)
    INTEGER(kind=kint ), ALLOCATABLE :: Temp_CN(:,:)
    INTEGER(kind=kint ), ALLOCATABLE :: Temp_N(:)
    !C temporary workspace of RAP
    INTEGER(kind=kint), INTENT(in) :: TEMP_RAP_WORKSIZE
    
    !C pointers for matrix
    REAL   (kind=kreal), POINTER  :: newD(:), newAU(:), newAL(:)
    INTEGER(kind=kint ), POINTER  :: newINU(:), newIAU(:), newINL(:), newIAL(:)

    !C communication data structures
    INTEGER(kind=kint )                   , INTENT(in) :: NPROCS
    INTEGER(kind=kint),INTENT(inout) :: local_aggre_size

    INTEGER(kind=kint) :: i, j, k, l, m, g, N, NP, row, colin, col, n_lower, n_upper
    INTEGER(kind=kint) :: rin, trin, global_CN, hash_no, count, index

    INTEGER(kind=kint) :: newNP, newN, newNPU, newNPL

    global_CN = coarser_level_size
    N  = HIERARCHICAL_DATA(LEVEL_NO - 1) % N 
    NP = HIERARCHICAL_DATA(LEVEL_NO - 1) % NP 

    
    ALLOCATE(Temp_N (local_aggre_size))
    ALLOCATE(Temp_V (TEMP_COLSIZE, local_aggre_size))
    ALLOCATE(Temp_CN(TEMP_COLSIZE, local_aggre_size))
    Temp_N (1:local_aggre_size) = 0
    Temp_CN(1:TEMP_COLSIZE, 1:local_aggre_size) = 0
    
    CALL RAP2(Temp_N, Temp_CN, Temp_V, TEMP_COLSIZE, N, &
         & coarser_level_size, LEVEL_NO, local_aggre_size, TEMP_RAP_WORKSIZE, SYM_FLAG)

    !C-- end:allocate temporaly domain and RAP().
    
    DO i = 1, local_aggre_size
       IF(Temp_N(i) > TEMP_COLSIZE) STOP 'TEMP_COLSIZE should be enlarged: RAP'
    END DO

    !C-- allocate the location and copy Temp matrix to correct location
    
    ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % INU(0:local_aggre_size))
    ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % INL(0:local_aggre_size))
    newINL => HIERARCHICAL_DATA(LEVEL_NO) % INL 
    newINU => HIERARCHICAL_DATA(LEVEL_NO) % INU 

    newINL = 0
    newINU = 0

    if(SYM_FLAG /= 1) then
       !C-UnSymmetric case
       newNPU = 0; newNPL = 0

       !C-- negligible nonzeros are cleared
!!$       DO row = 1, coarser_level_size
!!$          DO colin = 1, TEMP_COLSIZE
!!$             col = Temp_CN(colin, row)
!!$             IF(col > 0 .AND. EPS > abs(Temp_V(colin,row))) THEN
!!$                TEMP_CN(colin, row) = 0
!!$                Temp_N(row) = Temp_N(row) - 1 
!!$                Temp_V(colin, row) = 0
!!$             END IF
!!$          END DO
!!$       END DO

       !C-- structure is symmetrized. Upper part is considered first.
       DO row = 1, coarser_level_size
          DO colin = 1, TEMP_COLSIZE
             col = Temp_CN(colin, row)
             IF(col > 0) THEN
                IF(col < row) THEN
                   hash_no = MOD(row, TEMP_COLSIZE) + 1
                   DO count = hash_no, TEMP_COLSIZE
                      IF(Temp_CN(count, col) == row .or. Temp_CN(count, col) == 0) EXIT
                   END DO
                   IF(count > TEMP_COLSIZE) THEN
                      DO count = 1, hash_no - 1
                         IF(Temp_CN(count, col) == row .or. Temp_CN(count, col) == 0) EXIT
                      END DO
                      IF(count >= hash_no) stop 'error in RAP2()'
                   END IF

                   if(Temp_CN(count, col) == 0) then
                      !C-- the case of modification of the structure
                      Temp_N(col) = Temp_N(col) + 1
                      Temp_CN(count, col) = row
                      Temp_V(count, col) = 0.d0
                   end if
                END IF
             END IF
          END DO
       END DO

       !C- lower part is considered.
       DO row = 1, coarser_level_size
          DO colin = 1, TEMP_COLSIZE
             col = Temp_CN(colin, row)
             IF(col > row) THEN
                hash_no = MOD(row, TEMP_COLSIZE) + 1
                DO count = hash_no, TEMP_COLSIZE
                   IF(Temp_CN(count, col) == row .or. Temp_CN(count, col) == 0) EXIT
                END DO
                IF(count > TEMP_COLSIZE) THEN
                   DO count = 1, hash_no - 1
                      IF(Temp_CN(count, col) == row .or. Temp_CN(count, col) == 0) EXIT
                   END DO
                   IF(count >= hash_no) stop 'error in RAP2()'
                END IF
                
                if(Temp_CN(count, col) == 0) then
                   !C--  modification of the structure is needed.
                   Temp_N(col) = Temp_N(col) + 1
                   Temp_CN(count, col) = row
                   Temp_V(count, col) = 0.d0
                end if

             END IF
          END DO
       END DO

       !C-- after symmetrization of the nonzero structure, allocating arrays
       DO row = 1, coarser_level_size
          DO colin = 1, TEMP_COLSIZE
             col = Temp_CN(colin, row)
             IF(col > 0) THEN
                IF(col > row) THEN
                   newNPU = newNPU + 1
                   newINU(row) = newINU(row) + 1
                ELSE IF(col < row) THEN
                   newNPL = newNPL + 1
                   newINL(row) = newINL(row) + 1                
                END IF
             END IF
          END DO
       END DO
       DO i = 1, local_aggre_size
          newINU(i) = newINU(i-1) + newINU(i)
          newINL(i) = newINL(i-1) + newINL(i)
       END DO

       HIERARCHICAL_DATA(LEVEL_NO) % NPU = newNPU
       HIERARCHICAL_DATA(LEVEL_NO) % NPL = newNPL

       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % IAL(newNPL))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % IAU(newNPU))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % AL(newNPL))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % AU(newNPU))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % D(local_aggre_size))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % X(local_aggre_size))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % B(local_aggre_size))

       newD   => HIERARCHICAL_DATA(LEVEL_NO) % D 
       newIAL => HIERARCHICAL_DATA(LEVEL_NO) % IAL
       newIAU => HIERARCHICAL_DATA(LEVEL_NO) % IAU
       newAL  => HIERARCHICAL_DATA(LEVEL_NO) % AL
       newAU  => HIERARCHICAL_DATA(LEVEL_NO) % AU

       n_lower = 0
       n_upper = 0

       DO row = 1, coarser_level_size
          n_upper = newINU(row - 1)
          DO colin = 1, TEMP_COLSIZE
             col = Temp_CN(colin, row)
             IF(col > 0) THEN
                IF(row == col) THEN 
                   newD(row) = Temp_V(colin, row)
                ELSE IF(row < col) THEN
                   !C upper part
                   n_upper = n_upper + 1
                   newAU(n_upper) = Temp_V(colin, row)
                   newIAU(n_upper) = col - ZERO_ORIGIN
                ELSE IF(row > col) THEN
                   n_lower = n_lower + 1
                   newAL(n_lower) = Temp_V(colin, row)
                   newIAL(n_lower) =col - ZERO_ORIGIN               
                END IF
             END IF
          END DO
       END DO

    else
       !C-Symmetric case

       !C count upper nodes: symmetricity is assumed
       !C col is shifted for the assignment afterwards 
       count = 0
       DO row = 1, coarser_level_size
          DO colin = 1, TEMP_COLSIZE
             col = Temp_CN(colin, row)
             IF(col > row .AND. EPS < abs(Temp_V(colin, row))) THEN
                count = count + 1
                newINU(row) = newINU(row) + 1
                IF(col < local_aggre_size) newINL(col + 1) = newINL(col + 1) + 1
             END IF
          END DO
       END DO

       DO i = 1, local_aggre_size
          newINU(i) = newINU(i-1) + newINU(i)
          newINL(i) = newINL(i-1) + newINL(i)
       END DO

       HIERARCHICAL_DATA(LEVEL_NO) % NPU = count
       HIERARCHICAL_DATA(LEVEL_NO) % NPL = count

       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % IAL(count))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % IAU(count))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % AL(count))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % AU(count))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % D(local_aggre_size))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % X(local_aggre_size))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % B(local_aggre_size))

       newD   => HIERARCHICAL_DATA(LEVEL_NO) % D 
       newIAL => HIERARCHICAL_DATA(LEVEL_NO) % IAL
       newIAU => HIERARCHICAL_DATA(LEVEL_NO) % IAU
       newAL  => HIERARCHICAL_DATA(LEVEL_NO) % AL
       newAU  => HIERARCHICAL_DATA(LEVEL_NO) % AU

       n_lower = 0
       n_upper = 0

       DO row = 1, coarser_level_size
          n_upper = newINU(row - 1)
          DO colin = 1, TEMP_COLSIZE
             col = Temp_CN(colin, row)
             IF(col > 0) THEN
                IF(row == col) THEN 
                   newD(row) = Temp_V(colin, row)
                ELSE IF(row < col .AND. EPS < abs(Temp_V(colin, row))) THEN
                   !C upper part
                   n_upper = n_upper + 1
                   newAU(n_upper) = Temp_V(colin, row)
                   newIAU(n_upper) = col - ZERO_ORIGIN

                   !C lower part
                   n_lower = newINL(col)
                   n_lower = n_lower + 1
                   newINL(col) = n_lower
                   newAL(n_lower) = Temp_V(colin, row)
                   newIAL(n_lower) = row - ZERO_ORIGIN               
                END IF
             END IF
          END DO
       END DO
    end if
    
    DEALLOCATE(Temp_V)
    DEALLOCATE(Temp_CN)
    DEALLOCATE(Temp_N)

    HIERARCHICAL_DATA(LEVEL_NO) % N = coarser_level_size
    HIERARCHICAL_DATA(LEVEL_NO) % NP = local_aggre_size
    !C--end: allocate the location and copy Temp matrix to correct location
    
  END SUBROUTINE make_rap

  
  SUBROUTINE smooth_aggregate (PNI, NI, omega, LEVEL_NO, in_aggregates_result_size, &
       & in_aggregates_result, aggregates_result, Temp_P, Temp_P_row_size, tilde_D, &
       & Temp_P_SIZE)
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint ), INTENT(in) :: PNI(0:)
    INTEGER(kind=kint ), INTENT(in) :: NI(:)
    REAL   (kind=kreal), INTENT(in) :: omega
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO
    INTEGER(kind=kint ), INTENT(in) :: in_aggregates_result_size
    REAL   (kind=kreal), INTENT(in) :: tilde_D(:)
    INTEGER(kind=kint ), INTENT(in) :: Temp_P_SIZE
    
    INTEGER(kind=kint) :: in_aggregates_result(0:)
    INTEGER(kind=kint) :: aggregates_result(:)
    TYPE(row_node)     :: Temp_P(:,:)
    INTEGER(kind=kint) :: Temp_P_row_size(:)

    REAL   (kind=kreal),ALLOCATABLE :: Temp(:,:)

    REAL   (kind=kreal), POINTER :: D  (:)
    REAL   (kind=kreal), POINTER :: AU (:)
    REAL   (kind=kreal), POINTER :: AL (:)
    INTEGER(kind=kint ), POINTER :: INU(:)
    INTEGER(kind=kint ), POINTER :: IAU(:)
    INTEGER(kind=kint ), POINTER :: INL(:)
    INTEGER(kind=kint ), POINTER :: IAL(:)
    
    REAL   (kind=kreal) :: omega2
    INTEGER(kind=kint ) :: i, j, k, l, row, index, isL, ieL, isU, ieU, is, ie, NP, N
    INTEGER(kind=kint ) :: connected_node, connected_node_start, connected_node_end

    INTEGER(kind=kint)  :: row_no, row_size, A_row, A_col, P_col
    INTEGER(kind=kint), ALLOCATABLE :: Temp_int(:)

    NP = HIERARCHICAL_DATA(LEVEL_NO - 1) % NP
    N  = HIERARCHICAL_DATA(LEVEL_NO - 1) % N 

    ALLOCATE(Temp_int(NP))

    D   => HIERARCHICAL_DATA(LEVEL_NO - 1) % D 
    INL => HIERARCHICAL_DATA(LEVEL_NO - 1) % INL 
    INU => HIERARCHICAL_DATA(LEVEL_NO - 1) % INU 
    IAL => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AU


    DO i = 1, NP
       Temp_P_row_size(i) = 0
    END DO
    
    
    !C-- Temp_int() has aggregate NO for each node.
    Temp_int = 0
    DO k = 1, in_aggregates_result_size
       is = in_aggregates_result(k - 1) + 1
       ie = in_aggregates_result(k)
       DO i = is, ie
          row_no = aggregates_result(i)
          Temp_int(row_no) = k
       END DO
    END DO
    
    !C-- Temp_P's space is calculated for smoothed aggregations
    !C   and input column and value( initial values 0 )
    DO A_row = 1, NP
       !C  element's  value of the aggregate is 1 in array P      
!!$#ifdef SMOOTH_TILDE_A             
       isL = PNI(A_row-1) + 1
       ieL = PNI(A_row)
!!$#else
!!$       isL = INL(A_row-1) + 1
!!$       ieL = INL(A_row)
!!$#endif
       DO i = isL, ieL
!!$#ifdef SMOOTH_TILDE_A             
          l = NI(i)
!!$#else
!!$          l = i
!!$#endif
          A_col = IAL(l)+ZERO_ORIGIN
          P_col = Temp_int(A_col)
          IF(P_col == 0) CYCLE
          
          DO k = 1, Temp_P_row_size(A_row) 
             IF( Temp_P(k,A_row) % column == P_col ) EXIT
          END DO
          IF(k > Temp_P_row_size(A_row)) THEN
             IF(k <= Temp_P_SIZE) THEN
                Temp_P_row_size(A_row)    = k
                Temp_P(k, A_row) % column = P_col
                Temp_P(k, A_row) % value  = 0.0
             ELSE
                STOP 'Temp_P should be enlarged : smooth_aggregate'
             END IF
          END IF
       END DO
!!$#ifdef SMOOTH_TILDE_A             
       isU = PNI(A_row + NP -1) + 1
       ieU = PNI(A_row + NP)
!!$#else
!!$       isU = INU(A_row - 1) + 1
!!$       ieU = INU(A_row)
!!$#endif
       DO i = isU,ieU
!!$#ifdef SMOOTH_TILDE_A             
          l = NI(i)
!!$#else
!!$          l = i
!!$#endif
          A_col = IAU(l)+ZERO_ORIGIN
          P_col = Temp_int(A_col)
          IF(P_col == 0) CYCLE

          DO k=1, Temp_P_row_size(A_row)
             IF( Temp_P(k, A_row) % column == P_col) EXIT
          END DO
          IF(k > Temp_P_row_size(A_row)) THEN
             IF(k <= Temp_P_SIZE) THEN                
                Temp_P_row_size(A_row) = k
                Temp_P(k, A_row) % column = P_col
                Temp_P(k, A_row) % value  = 0.0
             ELSE
                STOP 'Temp_P should be enlarged : smooth_aggregate'                                      
             END IF
          END IF
       END DO
       
       !C- diagonal part
       A_col = A_row
       P_col = Temp_int(A_col)
       IF(P_col == 0) CYCLE
       
       DO k = 1, Temp_P_row_size(A_row)
          IF(Temp_P(k,A_row) % column == P_col) EXIT
       END DO
       IF(k > Temp_P_row_size(A_row)) THEN
          IF(k <= Temp_P_SIZE) THEN                
             Temp_P_row_size(A_row) = k
             Temp_P(k, A_row) % column = P_col
             Temp_P(k, A_row) % value  = 0.0
          ELSE
             STOP 'Temp_P should be enlarged : smooth_aggregate'                                      
          END IF
       END IF
    END DO

    !C-- smooth aggregates
    DO A_row = 1, NP
       row_size = Temp_P_row_size(A_row)
       IF(row_size == 0) CYCLE

       !C  element's  value of the aggregate is 1 in array P      
!!$#ifdef SMOOTH_TILDE_A             
       isL = PNI(A_row - 1) + 1
       ieL = PNI(A_row)
!!$#else
!!$       isL = INL(A_row - 1) + 1
!!$       ieL = INL(A_row)
!!$#endif
       DO i = isL, ieL
!!$#ifdef SMOOTH_TILDE_A             
          l = NI(i)
!!$#else
!!$          l = i
!!$#endif
          A_col = IAL(l)+ZERO_ORIGIN
          P_col = Temp_int(A_col)
          IF(P_col == 0) CYCLE
          
          DO k = 1, Temp_P_row_size(A_row) 
             IF( Temp_P(k,A_row) % column == P_col ) EXIT
          END DO
          IF(k > Temp_P_row_size(A_row)) THEN
             STOP 'unexpected error : smooth_aggregate()'
          END IF
          Temp_P(k, A_row) % value = Temp_P(k, A_row) % value - AL(l)
       END DO
       

!!$#ifdef SMOOTH_TILDE_A             
       isU = PNI(A_row + NP - 1) + 1
       ieU = PNI(A_row + NP)
!!$#else
!!$       isU = INU(A_row - 1) + 1
!!$       ieU = INU(A_row)
!!$#endif
       DO i = isU, ieU
!!$#ifdef SMOOTH_TILDE_A             
          l = NI(i)
!!$#else
!!$          l = i
!!$#endif
          A_col = IAU(l)+ZERO_ORIGIN
          P_col = Temp_int(A_col)
          IF(P_col == 0) CYCLE

          DO k = 1, Temp_P_row_size(A_row)
             IF( Temp_P(k,A_row) % column == P_col) EXIT
          END DO
          IF(k > Temp_P_row_size(A_row)) THEN
             STOP 'unexpected error : smooth_aggregate()'
          END IF
          
          Temp_P(k, A_row) % value = Temp_P(k, A_row) % value - AU(l)
       END DO
    END DO
    
    !C-- diagonal part
    DO A_row = 1, NP
       DO i= 1, Temp_P_row_size(A_row)
!!$#ifdef SMOOTH_TILDE_A                            
          Temp_P(i, A_row) % value = Temp_P(i, A_row) % value * omega / &
               & tilde_D(A_row)
!!$#else
!!$          Temp_P(i, A_row) % value = Temp_P(i, A_row) % value * omega / &
!!$               & D(A_row)
!!$#endif
       END DO
    END DO

    omega2 = 1. - omega
    DO row_no = 1, N
       P_col = Temp_int(row_no)
       IF(P_col == 0) CYCLE
       DO j = 1, Temp_P_row_size(row_no)
          IF(P_col == Temp_P(j, row_no) % column) EXIT
       END DO
       IF(j > Temp_P_row_size(row_no)) STOP 'unexpected error : smooth_aggregate()'
       Temp_P(j, row_no) % value = Temp_P(j, row_no) % value + omega2
    END DO
    !C-- end: smooth aggregates
    
    DEALLOCATE(Temp_int)
  END SUBROUTINE smooth_aggregate
  
  
  SUBROUTINE neighbors(PNI, NI, theta, LEVEL_NO, node_index, tilde_D)
    USE data_structure_for_AMG
    IMPLICIT NONE
    
    INTEGER(kind=kint ) :: N, NP
    REAL   (kind=kreal), INTENT(in):: theta
    
    INTEGER(kind=kint ), INTENT(out) :: NI  (:)
    INTEGER(kind=kint ), INTENT(out) :: PNI(0:)
    INTEGER(kind=kint ), INTENT(out) :: node_index(:)
    
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO
    REAL   (kind=kreal), POINTER :: D  (:)
    REAL   (kind=kreal), POINTER :: AU (:)
    REAL   (kind=kreal), POINTER :: AL (:)
    INTEGER(kind=kint ), POINTER :: INU(:)
    INTEGER(kind=kint ), POINTER :: IAU(:)
    INTEGER(kind=kint ), POINTER :: INL(:)
    INTEGER(kind=kint ), POINTER :: IAL(:)

    REAL(kind=kreal)  :: tilde_D(:)

    INTEGER(kind=kint) :: i,j,isU,ieU,isL,ieL,inod,nni,nni_before,size,count
    REAL   (kind=kreal):: t

    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP     
    N  = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    
    D   => HIERARCHICAL_DATA(LEVEL_NO-1) % D 
    INL => HIERARCHICAL_DATA(LEVEL_NO-1) % INL 
    INU => HIERARCHICAL_DATA(LEVEL_NO-1) % INU 
    IAL => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO-1) % AU
    
    PNI(0) = 0
    nni    = 0

    DO i = 1,NP
       tilde_D(i) = D(i)
       node_index(i) = 0
    END DO
    
    DO j = 1,NP
       isL = INL(j-1)+1
       ieL = INL(j  )

       DO i= isL, ieL
          inod= IAL(i)+ZERO_ORIGIN
          t= AL(i)

          IF(t*t > abs(D(j)*D(inod)*theta*theta) .AND. D(j)*D(inod)*t<0) THEN
             nni=nni+1
             NI(nni)=i
          ELSE
!!$             if(t*D(j)<0) then 
             tilde_D(j)=tilde_D(j)+t
!!$             end if
          ENDIF
       ENDDO
       PNI(j)=nni
    END DO

    DO j = 1,NP
       isU = INU(j-1)+1
       ieU = INU(j  )
       
       DO i = isU, ieU
          inod= IAU(i)+ZERO_ORIGIN
          t=AU(i)
          IF(t*t > abs(D(j)*D(inod)*theta*theta) .AND. D(j)*D(inod)*t<0) THEN
             nni=nni+1
             NI(nni)=i
          ELSE
!!$             if(t*D(j)<0) then 
             tilde_D(j)=tilde_D(j)+t
!!$             end if
          END IF
       ENDDO
       PNI(NP+j)=nni
    END DO

    DO j=1,N
       size=PNI(j)-PNI(j-1)+PNI(NP+j)-PNI(NP+j-1)
       IF(size==0 ) THEN
          node_index(j)=-1
       END IF
    ENDDO

  END SUBROUTINE neighbors



  !C this subroutine considers the case of unsymmetric problem matrix.
  !C But symmetricity of nonzero structure of the matrix is assumed.
  SUBROUTINE neighbors_unsym(PNI, NI, theta, LEVEL_NO, node_index, tilde_D)
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint ) :: N,NP
    REAL   (kind=kreal), INTENT(in):: theta
    
    INTEGER(kind=kint ), INTENT(out) :: NI  (:)
    INTEGER(kind=kint ), INTENT(out) :: PNI(0:)
    INTEGER(kind=kint ), INTENT(out) :: node_index(:)
    
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO
    REAL   (kind=kreal), POINTER :: D  (:)
    REAL   (kind=kreal), POINTER :: AU (:)
    REAL   (kind=kreal), POINTER :: AL (:)
    INTEGER(kind=kint ), POINTER :: INU(:)
    INTEGER(kind=kint ), POINTER :: IAU(:)
    INTEGER(kind=kint ), POINTER :: INL(:)
    INTEGER(kind=kint ), POINTER :: IAL(:)

    REAL(kind=kreal)  :: tilde_D(:)

    INTEGER(kind=kint), ALLOCATABLE :: check_U(:)

    INTEGER(kind=kint) :: i,j,k,isU,ieU,isL,ieL,inod,nni,nni_before,size,count
    INTEGER(kind=kint) :: NPU,NPL
    REAL   (kind=kreal):: tl,tu

    NP   = HIERARCHICAL_DATA(LEVEL_NO-1) % NP     
    N    = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    NPU  = HIERARCHICAL_DATA(LEVEL_NO-1) % NPU 
    NPL  = HIERARCHICAL_DATA(LEVEL_NO-1) % NPL 
    
    D   => HIERARCHICAL_DATA(LEVEL_NO-1) % D 
    INL => HIERARCHICAL_DATA(LEVEL_NO-1) % INL 
    INU => HIERARCHICAL_DATA(LEVEL_NO-1) % INU 
    IAL => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO-1) % AU
    
    PNI(0) = 0
    nni    = 0
    
    ALLOCATE(check_U(NPU))
    check_U=0
    
    DO i = 1,NP
       tilde_D(i) = D(i)
       node_index(i) = 0
    END DO
    
    DO j = 1,NP
       isL = INL(j-1)+1
       ieL = INL(j  )

       DO i = isL, ieL
          inod = IAL(i)+ZERO_ORIGIN
          tl = AL(i)
          
          !C determining the cordinate position of upper matrix
          isU = INU(inod - 1) + 1; ieU = INU(inod)
          DO k= isU, ieU
             IF(IAU(k)+ZERO_ORIGIN == j) EXIT
          END DO
          IF(k > ieU) STOP "error in neighbor()"
          tu = AU(k)
          
          IF(tl * tl > abs(D(j) * D(inod) * theta * theta) .AND. D(j) * D(inod) * tl < 0 .OR. &
               & tu * tu > abs(D(j) * D(inod) * theta * theta) .AND. D(j) * D(inod) * tu < 0) THEN
             check_U(k) = 1
             nni = nni + 1
             NI(nni) = i
          ELSE
             tilde_D(j) = tilde_D(j) + tl
          ENDIF
       ENDDO
       PNI(j) = nni
    END DO
    
    DO j = 1, NP
       isU = INU(j - 1) + 1; ieU = INU(j)
       DO i = isU, ieU
          IF(check_U(i) == 1) THEN
             nni = nni + 1
             NI(nni) = i
          ELSE
             tilde_D(j) = tilde_D(j) + AU(i)
          END IF
       ENDDO
       PNI(NP + j) = nni
    ENDDO
    
    DEALLOCATE(check_U)

    DO j = 1, N
       size = PNI(j) - PNI(j - 1) + PNI(NP + j) - PNI(NP + j - 1)
       IF(size == 0) node_index(j) = -1
    END DO

  END SUBROUTINE neighbors_unsym
END MODULE  data_creation_AMGCG
