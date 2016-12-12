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

#define SMOOTH_TILDE_A

#ifndef ZERO_ORIGIN
#define ZERO_ORIGIN 0
#endif

#ifdef LONG__LONG
#define LIS_MPI_INTEGER MPI_INTEGER8
#else
#define LIS_MPI_INTEGER MPI_INTEGER
#endif

!C   ************************************************
!C   * MODULE data_creation_AMGCG
!C     CONTAINS
!C   * SUBROUTINE data_creation_ssi_amg
!C   * SUBROUTINE data_creation_unsym_ssi_amg
!C   * SUBROUTINE data_creation
!C   * SUBROUTINE repeated_nds_tbl_crtn
!C   * SUBROUTINE make_PE_lists
!C   * SUBROUTINE RAP
!C   * SUBROUTINE RAP_unsym
!C   * SUBROUTINE RAP2
!C   * SUBROUTINE DEBUG_BARRIER_PRINT
!C   * SUBROUTINE make_rap
!C   * SUBROUTINE make_PSPSE_MAT
!C   * SUBROUTINE left_looking_chorescky
!C   * SUBROUTINE smooth_aggregate
!C   * SUBROUTINE smooth_aggregate_old
!C   * SUBROUTINE smooth_aggregate_unsym
!C   * SUBROUTINE lower_ext_matrix_crtn
!C   * SUBROUTINE neighbors
!C   * SUBROUTINE neighbors_unsym
!C   ************************************************  

MODULE data_creation_AMGCG
CONTAINS

  SUBROUTINE data_creation_ssi_amg(N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,  &
       & NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,     &
       & NOD_IMPORT_SIZE, NOD_EXPORT_SIZE, SOLVER_COMM, LEVEL_NUM, WSIZE, theta)
    USE data_structure_for_AMG
    IMPLICIT NONE
    include 'mpif.h'
    
    INTEGER(kind=kint), intent(out)  :: WSIZE, LEVEL_NUM
    INTEGER(kind=kint), intent(in)   :: N, NP, NPL, NPU, NOD_IMPORT_SIZE, NOD_EXPORT_SIZE
    !C NP
    REAL   (kind=kreal) :: D(NP)
    !C NPU
    REAL   (kind=kreal) :: AU(NPU)
    !C NPL
    REAL   (kind=kreal) :: AL(NPL)
    
    !C 0:NP, 0:NP
    INTEGER(kind=kint ) :: INU(0:N), INL(0:N)
    
    !C NPU, NPL
    INTEGER(kind=kint ) :: IAU(1:NPU), IAL(1:NPL)
    
    INTEGER(kind=kint)  :: my_rank
    INTEGER(kind=kint), INTENT(in) :: NEIBPETOT
    INTEGER(kind=kint), INTENT(in) :: SOLVER_COMM
    !C-     NEIBPE      (NEIBPETOT) 
    !C-     STACK_IMPORT(0:NEIBPETOT)
    !C-     NOD_IMPORT  (:)
    !C-     STACK_EXPORT(0:NEIBPETOT)
    !C-     NOD_EXPORT  (:)
    INTEGER(kind=kint ) :: NEIBPE      (1:NEIBPETOT) 
    INTEGER(kind=kint ) :: STACK_IMPORT(0:NEIBPETOT)
    INTEGER(kind=kint ) :: NOD_IMPORT  (1:NOD_IMPORT_SIZE)
    INTEGER(kind=kint ) :: STACK_EXPORT(0:NEIBPETOT)
    INTEGER(kind=kint ) :: NOD_EXPORT  (1:NOD_EXPORT_SIZE)

    CHARACTER(len=20)   :: PRECOND
    INTEGER(kind=kint) :: ierr
    INTEGER(kind=kint) :: i,j,k,l,m

    REAL   (kind=kreal) :: theta 
    !C-- theta means the critierion for strong connection.
    !C-- if theta=0 then all nonzero elements are recognized as
    !C-- strong connection.
    
    CALL MPI_COMM_RANK (SOLVER_COMM, my_rank, ierr )    
    PRECOND = "SHA2"
!!$    PRECOND = "IND1"
!!$    PRECOND = "MIX"


    call data_creation(N, NP, NPL, NPU, D(1:NP), AL(1:NPL), INL(0:N), IAL(1:NPL), &
         & AU(1:NPU), INU(0:N), IAU(1:NPU), LEVEL_NUM,   &
         & my_rank, NEIBPETOT, NEIBPE(1:NEIBPETOT), STACK_IMPORT(0:NEIBPETOT), &
         & NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)), STACK_EXPORT(0:NEIBPETOT),      &
         & NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)), SOLVER_COMM, PRECOND, 1, theta)


    
    !C- WSIZE

    WSIZE = 0
    !C level 1
    if(NEIBPETOT > 0) then
       if( STACK_IMPORT(NEIBPETOT) > STACK_EXPORT(NEIBPETOT) ) then
          WSIZE = STACK_IMPORT(NEIBPETOT)
       else
          WSIZE = STACK_EXPORT(NEIBPETOT)
       end if
    end if
    if( WSIZE < NP ) WSIZE = NP

    !C level 2..
    if(NEIBPETOT > 0) then
       m = WSIZE
       DO i = 2, LEVEL_NUM
          j = HIERARCHICAL_DATA(i) % COMM_TABLE % NEIBPETOT 
          k = HIERARCHICAL_DATA(i) % COMM_TABLE % STACK_IMPORT(j)
          l = HIERARCHICAL_DATA(i) % COMM_TABLE % STACK_EXPORT(j)
          m = max(k,l,m)

          j = HIERARCHICAL_DATA(i) % INT_LVL_TABLE % NEIBPETOT 
          k = HIERARCHICAL_DATA(i) % INT_LVL_TABLE % STACK_IMPORT(j)
          l = HIERARCHICAL_DATA(i) % INT_LVL_TABLE % STACK_EXPORT(j)
          m = max(k,l,m)
       END DO
       if(m > WSIZE) then
          !C-- WS, WR must be reallocated
          WSIZE = m
       end if
    end if

  END SUBROUTINE data_creation_ssi_amg

  SUBROUTINE data_creation_unsym_ssi_amg(N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,  &
       & NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,     &
       & NOD_IMPORT_SIZE,  NOD_EXPORT_SIZE, SOLVER_COMM, LEVEL_NUM, WSIZE, theta)
    USE data_structure_for_AMG
    IMPLICIT NONE
    include 'mpif.h'
    
    INTEGER(kind=kint), intent(out)  :: WSIZE, LEVEL_NUM
    INTEGER(kind=kint), intent(in)   :: N, NP, NPL, NPU, NOD_IMPORT_SIZE, NOD_EXPORT_SIZE
    !C NP
    REAL   (kind=kreal) :: D(NP)
    !C NPU
    REAL   (kind=kreal) :: AU(NPU)
    !C NPL
    REAL   (kind=kreal) :: AL(NPL)
    
    !C 0:NP, 0:NP
    INTEGER(kind=kint ) :: INU(0:N), INL(0:N)
    !C NPU, NPL
    INTEGER(kind=kint ) :: IAU(NPU), IAL(NPL)
    
    INTEGER(kind=kint)  :: my_rank
    INTEGER(kind=kint), INTENT(in) :: NEIBPETOT
    INTEGER, INTENT(in) :: SOLVER_COMM
    !C-     NEIBPE      (NEIBPETOT) 
    !C-     STACK_IMPORT(0:NEIBPETOT)
    !C-     NOD_IMPORT  (:)
    !C-     STACK_EXPORT(0:NEIBPETOT)
    !C-     NOD_EXPORT  (:)
    INTEGER(kind=kint ) :: NEIBPE      (NEIBPETOT) 
    INTEGER(kind=kint ) :: STACK_IMPORT(0:NEIBPETOT)
    INTEGER(kind=kint ) :: NOD_IMPORT  (NOD_IMPORT_SIZE)
    INTEGER(kind=kint ) :: STACK_EXPORT(0:NEIBPETOT)
    INTEGER(kind=kint ) :: NOD_EXPORT  (NOD_EXPORT_SIZE)

    CHARACTER(len=20)   :: PRECOND
    INTEGER(kind=kint) :: ierr
    INTEGER(kind=kint) :: i,j,k,l,m

    REAL   (kind=kreal) :: theta 
    !C-- theta means the critierion for strong connection.
    !C-- if theta=0 then all nonzero elements are recognized as
    !C-- strong connection.
    

    CALL MPI_COMM_RANK (SOLVER_COMM, my_rank, ierr )    


    PRECOND = "SHA2"

    call data_creation(N, NP, NPL, NPU, D, AL, INL, IAL, &
         & AU, INU, IAU, LEVEL_NUM, &
         & my_rank, NEIBPETOT, NEIBPE, STACK_IMPORT(0:NEIBPETOT), &
         & NOD_IMPORT, STACK_EXPORT(0:NEIBPETOT), &
         & NOD_EXPORT, SOLVER_COMM, PRECOND, 0, theta)

    !C- WSIZE

    WSIZE = 0
    !C level 1
    if(NEIBPETOT > 0) then
       if( STACK_IMPORT(NEIBPETOT) > STACK_EXPORT(NEIBPETOT) ) then
          WSIZE = STACK_IMPORT(NEIBPETOT)
       else
          WSIZE = STACK_EXPORT(NEIBPETOT)
       end if
    end if
    if( WSIZE < NP ) WSIZE = NP

    !C level 2..
    if(NEIBPETOT > 0) then
       m = WSIZE
       DO i = 2, LEVEL_NUM
          j = HIERARCHICAL_DATA(i) % COMM_TABLE % NEIBPETOT 
          k = HIERARCHICAL_DATA(i) % COMM_TABLE % STACK_IMPORT(j)
          l = HIERARCHICAL_DATA(i) % COMM_TABLE % STACK_EXPORT(j)
          m = max(k,l,m)

          j = HIERARCHICAL_DATA(i) % INT_LVL_TABLE % NEIBPETOT 
          k = HIERARCHICAL_DATA(i) % INT_LVL_TABLE % STACK_IMPORT(j)
          l = HIERARCHICAL_DATA(i) % INT_LVL_TABLE % STACK_EXPORT(j)
          m = max(k,l,m)
       END DO
       if(m > WSIZE) then
          !C-- WS, WR must be reallocated
          WSIZE = m
       end if
    end if

  END SUBROUTINE data_creation_unsym_ssi_amg


!==========================================================================
  SUBROUTINE data_creation(N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, LEVEL_NUM, &
       & my_rank, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, STACK_EXPORT,          &
       & NOD_EXPORT, SOLVER_COMM, PRECOND, symmetric_mat, theta)
    USE isort
    USE solver_SR2
    USE data_structure_for_AMG
    USE solver_Gnumbering
    USE aggregate_mod
    USE count_time_mod

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER(kind=kint), INTENT(in) :: N, NP, NPL, NPU
    INTEGER(kind=kint), INTENT(out):: LEVEL_NUM
    INTEGER(kind=kint) :: NPROCS

    !C NP
    REAL   (kind=kreal), TARGET :: D(1:NP)
    REAL   (kind=kreal), ALLOCATABLE :: tilde_D(:)
    !C NPU
    REAL   (kind=kreal), TARGET :: AU(1:NPU)
    !C NPL
    REAL   (kind=kreal), TARGET :: AL(1:NPL)
    
    !C 0:NP, 0:NP
    INTEGER(kind=kint ), TARGET :: INU(0:N), INL(0:N)
    !C NPU, NPL
    INTEGER(kind=kint ), TARGET :: IAU(1:NPU), IAL(1:NPL)
    INTEGER(kind=kint ), ALLOCATABLE :: NI(:), PNI(:)

    !C symmetric_mat == 1 means symmetric matrix
    INTEGER(kind=kint ), intent(in) :: symmetric_mat
    
    INTEGER(kind=kint), INTENT(in) :: my_rank
    INTEGER(kind=kint), INTENT(in) :: NEIBPETOT
    INTEGER(kind=kint), INTENT(in) :: SOLVER_COMM
    CHARACTER(len=20), INTENT(in) :: PRECOND

    !C-     NEIBPE      (NEIBPETOT) 
    !C-     STACK_IMPORT(0:NEIBPETOT)
    !C-     NOD_IMPORT  (:)
    !C-     STACK_EXPORT(0:NEIBPETOT)
    !C-     NOD_EXPORT  (:)
    INTEGER(kind=kint ), TARGET :: NEIBPE      (1:NEIBPETOT) 
    INTEGER(kind=kint ), TARGET :: STACK_IMPORT(0:NEIBPETOT)
    INTEGER(kind=kint ), TARGET :: NOD_IMPORT  (:)
    INTEGER(kind=kint ), TARGET :: STACK_EXPORT(0:NEIBPETOT)
    INTEGER(kind=kint ), TARGET :: NOD_EXPORT  (:)
    INTEGER(kind=kint ) :: ierr

    INTEGER(kind=kint ), ALLOCATABLE :: sorted_fNEIBPE(:)
    INTEGER(kind=kint )              :: finerNEIBPETOT
    INTEGER(kind=kint ), POINTER :: finerNEIBPE      (:) 
    INTEGER(kind=kint ), POINTER :: finerSTACK_IMPORT(:)
    INTEGER(kind=kint ), POINTER :: finerNOD_IMPORT  (:)
    INTEGER(kind=kint ), POINTER :: finerSTACK_EXPORT(:)

    INTEGER(kind=kint ), POINTER :: finerNOD_EXPORT(:)
    INTEGER(kind=kint ), ALLOCATABLE :: node_index(:)

    INTEGER(kind=kint) :: coarser_level_size, level, finer_level_size, count
    INTEGER(kind=kint) :: c1, c2, c3, c4, k, finer_level_NP, i, j, agmtd
    INTEGER(kind=kint) :: finer_level_NPL, finer_level_NPU

    INTEGER(kind=kint ) :: PE_num_size
    INTEGER(kind=kint ), POINTER :: PE_nums(:)
    INTEGER(kind=kint ), POINTER :: nodes_for_each_PE(:)
    INTEGER(kind=kint ), POINTER :: in_nodes_for_each_PE(:)

    INTEGER(kind=kint ), POINTER :: PE_list(:,:)
    INTEGER(kind=kint )          :: PE_list_size 

    REAL   (kind=kreal) :: theta 
    !C-- theta means the critierion for strong connection.
    !C-- if theta=0 then all nonzero elements are recognized as
    !C-- strong connection.
    
    INTEGER(kind=kint ) :: TEMP_RAP_WORKSIZE, NONZEROS_PER_ROW, TEMP_COLSIZE
    REAL   (kind=kreal) :: RATE_OF_SPACE

    INTEGER(kind=kint), POINTER :: in_aggregates_result(:)
    INTEGER(kind=kint), POINTER :: aggregates_result(:)
    INTEGER(kind=kint) :: in_aggregates_result_size

    TYPE(row_node),     POINTER :: Temp_P(:,:)
    INTEGER(kind=kint), POINTER :: Temp_P_row_size(:)

    REAL(kind=kreal),   ALLOCATABLE  :: WS(:), WR(:)
    INTEGER(kind=kint), ALLOCATABLE  :: WSI(:), WRI(:)

    LOGICAL            :: finish_flag, global_finish_flag, local_flag, global_flag
    LOGICAL            :: border_flag
    INTEGER(kind=kint) :: global_CN
    
    !C - begin: aggregate_table_array
    INTEGER(kind=kint) :: aggregate_table_size
    !C two arguments: communication table, each PE
    INTEGER(kind=kint), POINTER :: aggregate_table_array(:,:)
    !C argument: the number of aggregates
    INTEGER(kind=kint), POINTER :: aggregate_number_in_table(:)
    !C - end: aggregate_table_array
    INTEGER(kind=kint), ALLOCATABLE :: GIN_aggregate(:)

    !C - begin: hash_table
    INTEGER(kind=kint) :: HASH_TABLE_SIZE
    INTEGER(kind=kint), POINTER :: global_local_hash_table(:,:)
    !C - end: hash_table
    INTEGER(kind=kint) :: local_aggre_size,col
    TYPE(INTER_LEVEL_OPERATOR),POINTER :: P,R
    INTEGER(kind=kint), ALLOCATABLE :: Temp_R_IN(:)


    !C debug for lower_ext_matrix_crtn
    INTEGER(kind=kint ) :: EX_NPL
    INTEGER(kind=kint ), POINTER :: EX_INL (:), EX_IAL (:)
    REAL   (kind=kreal), POINTER :: EX_AL  (:)


    REAL(kind=kreal) :: dump_jacobi_weight
    INTEGER(kind=kint) :: max_neib_pe_size

    INTEGER(kind=kint) :: nfepe_size

    !C debug
    CHARACTER(len= 4 )  :: penum

    EX_INL=>NULL()
    EX_IAL=>NULL()
    EX_AL=>NULL()


    CALL MPI_COMM_SIZE(SOLVER_COMM, NPROCS, IERR)

    
    !C-- upper bound of the number of the external aggregates
    !C-- on next coarser level
    !C-- The allocated memory size is 8*HASH_TABLE_SIZE byte.
    HASH_TABLE_SIZE = 25849
!!$    HASH_TABLE_SIZE = 65029

    !C-- upper bound of the number of nonzeros per row
    !C-- on next coarser level
    !C-- The allocated memory size is
    !C-- coarse level unknowns*12*NONZEROS_PER_ROW byte.
!!$    NONZEROS_PER_ROW = 903
!!$    NONZEROS_PER_ROW = 401
    NONZEROS_PER_ROW = 401

    
    !C-- upper bound of the number of nonzeros per row
    !C-- in prolongation matrix. This value is smaller
    !C-- than NONZEROS_PER_ROW
    !C-- The allocated memory size is
    !C-- unknowns*12*RATE_OF_SPACE byte.
    RATE_OF_SPACE = 19
!!$    RATE_OF_SPACE = 5
    
!!$    !C-- work area of the calculation of RAP
!!$    !C-- this value must be bigger than or equal to 1.
!!$    !C-- The allocated memory size is
!!$    !C-- unknowns*8*TEMP_RAP_WORKSIZE byte.
!!$    TEMP_RAP_WORKSIZE = 5

!!$#ifdef ANISO
!!$    theta = 0.12
!!$#else
!!$    !C hoho
!!$    theta = 0.05
!!$!!$    theta = 0.2
!!$#endif
  
!!$    dump_jacobi_weight = 0.88888
    dump_jacobi_weight = 0.66666

    !C-- max number of external aggregates
    !C-- The allocated memory size is
    !C-- aggregate_table_size*2*4*max_neib_pe_size
    aggregate_table_size = 12000
!!$    aggregate_table_size = 7000

    max_neib_pe_size = NPROCS
!!$    max_neib_pe_size = 5
    if(NPROCS < max_neib_pe_size ) then
       max_neib_pe_size = NPROCS
    end if
       
    if(symmetric_mat == 1) then
       SYMMETRIC_FLAG = .true.
    else
       SYMMETRIC_FLAG = .false.
    end if

    agmtd = 0 
    IF (       ((PRECOND(1:1).EQ.'I').OR.(PRECOND(1:1).EQ.'i')) .AND.       &
         &     ((PRECOND(2:2).EQ.'N').OR.(PRECOND(2:2).EQ.'n')) .AND.       &
         &     ((PRECOND(3:3).EQ.'D').OR.(PRECOND(3:3).EQ.'d')) .AND.       &
         &     (PRECOND(4:4) .EQ.'1') ) THEN
       agmtd = 1
    ELSE IF (  ((PRECOND(1:1).EQ.'I').OR.(PRECOND(1:1).EQ.'i')) .AND.       &
         &     ((PRECOND(2:2).EQ.'N').OR.(PRECOND(2:2).EQ.'n')) .AND.       &
         &     ((PRECOND(3:3).EQ.'D').OR.(PRECOND(3:3).EQ.'d')) .AND.       &
         &     (PRECOND(4:4) .EQ.'2') ) THEN
       agmtd = 2
    ELSE IF (  ((PRECOND(1:1).EQ.'S').OR.(PRECOND(1:1).EQ.'s')) .AND.       &
         &     ((PRECOND(2:2).EQ.'H').OR.(PRECOND(2:2).EQ.'h')) .AND.       &
         &     ((PRECOND(3:3).EQ.'A').OR.(PRECOND(3:3).EQ.'a')) .AND.       &
         &     (PRECOND(4:4) .EQ.'1') ) THEN
       agmtd = 3
    ELSE IF (  ((PRECOND(1:1).EQ.'S').OR.(PRECOND(1:1).EQ.'s')) .AND.       &
         &     ((PRECOND(2:2).EQ.'H').OR.(PRECOND(2:2).EQ.'h')) .AND.       &
         &     ((PRECOND(3:3).EQ.'A').OR.(PRECOND(3:3).EQ.'a')) .AND.       &
         &     (PRECOND(4:4).EQ.'2') ) THEN
       agmtd = 4
    ELSE IF (  ((PRECOND(1:1).EQ.'S').OR.(PRECOND(1:1).EQ.'s')) .AND.       &
         &     ((PRECOND(2:2).EQ.'H').OR.(PRECOND(2:2).EQ.'h')) .AND.       &
         &     ((PRECOND(3:3).EQ.'A').OR.(PRECOND(3:3).EQ.'a')) .AND.       &
         &     (PRECOND(4:4).EQ.'3') ) THEN
       agmtd = 5
    ELSE IF (  ((PRECOND(1:1).EQ.'M').OR.(PRECOND(1:1).EQ.'m')) .AND.       &
         &     ((PRECOND(2:2).EQ.'I').OR.(PRECOND(2:2).EQ.'i')) .AND.       &
         &     ((PRECOND(3:3).EQ.'X').OR.(PRECOND(3:3).EQ.'x')) ) THEN
       agmtd = 6
    ENDIF
    IF(agmtd == 0)THEN 
       IF(my_rank == 0) THEN
          WRITE(*,*) "Select correct aggregate strategy"
          WRITE(*,*) "Aggregation strategies are ind1,ind2,sha1,sha2,sha3,and mix."
       END IF
       STOP 
    END IF


    HIERARCHICAL_DATA(1) % N   = N
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

    HIERARCHICAL_DATA(1) % COMM_TABLE % NEIBPETOT = NEIBPETOT
    if(NEIBPETOT > 0) then
       HIERARCHICAL_DATA(1) % COMM_TABLE % NEIBPE => NEIBPE
       HIERARCHICAL_DATA(1) % COMM_TABLE % STACK_IMPORT => STACK_IMPORT
       HIERARCHICAL_DATA(1) % COMM_TABLE % STACK_EXPORT => STACK_EXPORT
       HIERARCHICAL_DATA(1) % COMM_TABLE % NOD_IMPORT => NOD_IMPORT
       HIERARCHICAL_DATA(1) % COMM_TABLE % NOD_EXPORT => NOD_EXPORT
    end if



    !C PE_list(*,1): PE number
    !C PE_list(*,2): NEIBPE number|
    !C PE_list(*,3): PE_num number|these two has 0 if no correspondings.
    ALLOCATE(PE_list(NPROCS, 3))
    

    DO level = 2, MAX_LEVEL_SIZE
       finer_level_size = HIERARCHICAL_DATA(level - 1) % N
       finer_level_NP   = HIERARCHICAL_DATA(level - 1) % NP
       finer_level_NPL  = HIERARCHICAL_DATA(level - 1) % NPL
       finer_level_NPU  = HIERARCHICAL_DATA(level - 1) % NPU       
       
       finerNEIBPETOT    =   HIERARCHICAL_DATA(level - 1) % COMM_TABLE % NEIBPETOT
       finerNEIBPE       =>  HIERARCHICAL_DATA(level - 1) % COMM_TABLE % NEIBPE
       finerSTACK_IMPORT =>  HIERARCHICAL_DATA(level - 1) % COMM_TABLE % STACK_IMPORT
       finerSTACK_EXPORT =>  HIERARCHICAL_DATA(level - 1) % COMM_TABLE % STACK_EXPORT
       finerNOD_IMPORT   =>  HIERARCHICAL_DATA(level - 1) % COMM_TABLE % NOD_IMPORT
       finerNOD_EXPORT   =>  HIERARCHICAL_DATA(level - 1) % COMM_TABLE % NOD_EXPORT
       
       
       ALLOCATE(sorted_fNEIBPE(finerNEIBPETOT))
       IF(NPROCS > 1 .AND. N>0) THEN
          CALL QSORTI(sorted_fNEIBPE, finerNEIBPETOT, finerNEIBPE)
       END IF


#ifdef DEBUG       
       allocate(Temp_R_IN(finer_level_size))
       Temp_R_IN=0

       if(NEIBPETOT>0) then
          do i = 1, finerSTACK_EXPORT(finerNEIBPETOT)
             Temp_R_IN(finerNOD_EXPORT(i)+ZERO_ORIGIN)=Temp_R_IN(finerNOD_EXPORT(i)+ZERO_ORIGIN)+1       
          end do
          do i = 1, finer_level_size
             if(Temp_R_IN(i) >= NPROCS) then
                write(*,*) Temp_R_IN(i), i, "Error in export table!!"
                stop
             end if
          end do
       end if
       deallocate(Temp_R_IN)
#endif       

       CALL DEBUG_BARRIER_PRINT("## repeated_nds_tbl_crtn ##", SOLVER_COMM, my_rank)
       
       IF(NPROCS > 1) THEN
           IF (N>0) THEN
               CALL repeated_nds_tbl_crtn(level, sorted_fNEIBPE, SOLVER_COMM, my_rank, &
               & PE_num_size, PE_nums, nodes_for_each_PE, in_nodes_for_each_PE, NPROCS)
           ELSE
               PE_num_size=0
           END IF
       END IF

       CALL DEBUG_BARRIER_PRINT("## lower_ext_matrix_crtn ##", SOLVER_COMM, my_rank)
       
       if(.not.SYMMETRIC_FLAG) then
            IF (N>0) THEN
                CALL lower_ext_matrix_crtn(EX_NPL, EX_INL, EX_IAL, EX_AL, level, my_rank, &
                                           SOLVER_COMM, NPROCS)
            ELSE
                EX_NPL=0
            END IF
        else
            EX_NPL = 0
        end if
          
       CALL DEBUG_BARRIER_PRINT("## neighbors             ##", SOLVER_COMM, my_rank)
       
       ALLOCATE(NI(finer_level_NPL + finer_level_NPU + EX_NPL))
       ALLOCATE(PNI(0:2 * finer_level_NP))
       ALLOCATE(node_index(finer_level_NP))
       ALLOCATE(tilde_D   (finer_level_NP))

       CALL count_time(1, solver_comm, my_rank, 5)                     

       if(SYMMETRIC_FLAG) then
          CALL neighbors(2 * finer_level_NP, finer_level_NPL + finer_level_NPU + EX_NPL, &
               & finer_level_NP, PNI, NI, theta, level, node_index, tilde_D)
       else 
          CALL neighbors_unsym(EX_INL, EX_IAL, EX_AL, PNI, NI, theta, level, node_index, tilde_D)
       end if


       
       CALL count_time(2, solver_comm, my_rank, 5)                     


       CALL DEBUG_BARRIER_PRINT("## exchange node_index   ##", SOLVER_COMM, my_rank)


       IF(NPROCS > 1 .and. finerNEIBPETOT > 0) THEN

          i = finerSTACK_EXPORT(finerNEIBPETOT)
          j = finerSTACK_IMPORT(finerNEIBPETOT)
          ALLOCATE(WSI(i))
          ALLOCATE(WRI(j))

          

          CALL SOLVER_SEND_RECV2I                                     &
               &   (finer_level_NP, i, j, finerNEIBPETOT,             &
               &    finerNEIBPE,                    &
               &    finerSTACK_IMPORT,              &
               &    finerNOD_IMPORT,                             &
               &    finerSTACK_EXPORT,              &
               &    finerNOD_EXPORT, WSI, WRI,         &
               &    node_index,                     &
               &    SOLVER_COMM, my_rank)
          
          DEALLOCATE(WSI)
          DEALLOCATE(WRI)
       END IF

       CALL DEBUG_BARRIER_PRINT("## exchange tilde_D      ##", SOLVER_COMM, my_rank)






       IF(NPROCS > 1 .and. finerNEIBPETOT > 0) THEN
          ALLOCATE (WS(i))
          ALLOCATE (WR(j))
          CALL SOLVER_SEND_RECV2                                &
               &   (finer_level_NP, i, j, finerNEIBPETOT,       &
               &    finerNEIBPE,              &
               &    finerSTACK_IMPORT,        &
               &    finerNOD_IMPORT,                       &
               &    finerSTACK_EXPORT,        &
               &    finerNOD_EXPORT, WS, WR,     &
               &    tilde_D,                  &
               &    SOLVER_COMM, my_rank)
          DEALLOCATE(WS)
          DEALLOCATE(WR)
       END IF




       !C-- switch aggregation strategy
       !C single case takes up independent aggregation
       IF(NPROCS == 1) agmtd = 1

       CALL count_time(1, solver_comm, my_rank, 6)       
       IF(agmtd == 1 .OR. agmtd == 2) THEN
          
          CALL DEBUG_BARRIER_PRINT("## independent aggregate ##", SOLVER_COMM, my_rank)
          
          border_flag = .TRUE.
          IF(agmtd == 1) border_flag = .FALSE.
          
          CALL indpdt_agrgt (PNI, NI, level, node_index, in_aggregates_result, &
               & aggregates_result, in_aggregates_result_size, border_flag)
          
          CALL make_PE_lists (level, sorted_fNEIBPE, PE_num_size, PE_nums, PE_list_size, &
               &              PE_list, NPROCS)
          DEALLOCATE(sorted_fNeibPE)
          ALLOCATE(aggregate_table_array(aggregate_table_size * 2, max_neib_pe_size))
          ALLOCATE(aggregate_number_in_table(max_neib_pe_size))
          ALLOCATE(GIN_aggregate(0:NPROCS))

          !C--initialize aggregate_table and set GIN_aggregate(:) for indpdt_agrgt
          aggregate_number_in_table = 0
          CALL MPI_ALLGATHER(in_aggregates_result_size, 1, LIS_MPI_INTEGER, &
               & GIN_aggregate(1), 1, LIS_MPI_INTEGER, SOLVER_COMM, ierr)
          GIN_aggregate(0) = 0
          DO i = 2, NPROCS
             GIN_aggregate(i) = GIN_aggregate(i - 1) + GIN_aggregate(i)
          END DO
          !C--end:initialize aggregate_table and set GIN_aggregate(:) for indpdt_agrgt


       ELSE IF(agmtd == 3 .OR. agmtd == 4 .OR. agmtd == 5) THEN
          
          
          CALL DEBUG_BARRIER_PRINT("## border_aggregate      ##", SOLVER_COMM, my_rank)


          CALL border_aggregate (PNI(0:2*finer_level_NP), &
               &                 NI(1:finer_level_NPL+finer_level_NPU+EX_NPL), level, &
               &                 node_index, sorted_fNEIBPE(1:finerNEIBPETOT),        &
               &                 my_rank, PE_num_size, PE_nums, nodes_for_each_PE,    &
               &                 in_nodes_for_each_PE, PE_list_size, PE_list, NPROCS, &
               &                 SOLVER_COMM)
          DEALLOCATE(sorted_fNeibPE)



          CALL DEBUG_BARRIER_PRINT("## aggregate             ##", SOLVER_COMM, my_rank)
          

          ALLOCATE(aggregate_table_array(aggregate_table_size * 2, max_neib_pe_size))
          ALLOCATE(aggregate_number_in_table(max_neib_pe_size))
          ALLOCATE(GIN_aggregate(0:NPROCS))

          CALL aggregate(PNI(0:2*finer_level_NP), &
               &          NI(1:finer_level_NPL+finer_level_NPU+EX_NPL), level,       &
               &          node_index(1:finer_level_NP), in_aggregates_result,        &
               &          aggregates_result, in_aggregates_result_size, my_rank,     &
               &          aggregate_table_size,                                      &
               &          aggregate_table_array(1:aggregate_table_size*2, 1:max_neib_pe_size), &
               &          aggregate_number_in_table(1:max_neib_pe_size), SOLVER_COMM,NPROCS,   &
               &          GIN_aggregate(0:NPROCS), PE_num_size, PE_nums,             &
               &          nodes_for_each_PE, &
               &          in_nodes_for_each_PE, PE_list_size, PE_list, agmtd - 2)
          
       ELSE IF(agmtd == 6) THEN
          local_flag = finer_level_size > SWITCH_SIZE
          CALL MPI_ALLREDUCE(local_flag, global_flag, 1, MPI_LOGICAL, MPI_LAND, SOLVER_COMM, ierr)
          IF(global_flag) THEN

             CALL DEBUG_BARRIER_PRINT("## independent aggregate ##", SOLVER_COMM, my_rank)
             
             CALL indpdt_agrgt(PNI, NI, level, node_index, in_aggregates_result, &
                  & aggregates_result, in_aggregates_result_size, .TRUE.)


             CALL make_PE_lists(level, sorted_fNEIBPE, PE_num_size, PE_nums, PE_list_size, &
                  &              PE_list, NPROCS)
             DEALLOCATE(sorted_fNeibPE)
             ALLOCATE(aggregate_table_array(aggregate_table_size * 2, max_neib_pe_size))
             ALLOCATE(aggregate_number_in_table(max_neib_pe_size))
             ALLOCATE(GIN_aggregate(0:NPROCS))

             !C--initialize aggregate_table and set GIN_aggregate(:) for indpdt_agrgt
             aggregate_number_in_table = 0
             CALL MPI_ALLGATHER(in_aggregates_result_size, 1, LIS_MPI_INTEGER, &
                  & GIN_aggregate(1), 1, LIS_MPI_INTEGER, SOLVER_COMM, ierr)
             GIN_aggregate(0) = 0
             DO i = 2, NPROCS
                GIN_aggregate(i) = GIN_aggregate(i - 1) + GIN_aggregate(i)
             END DO
             !C--end:initialize aggregate_table and set GIN_aggregate(:) for indpdt_agrgt

          ELSE 

             CALL DEBUG_BARRIER_PRINT("## border_aggregate      ##", SOLVER_COMM, my_rank)

             CALL border_aggregate (PNI, NI, level, node_index, sorted_fNEIBPE,        &
                  &                 my_rank, PE_num_size, PE_nums, nodes_for_each_PE,  &
                  &                 in_nodes_for_each_PE, PE_list_size, PE_list,NPROCS,&
                  &                 SOLVER_COMM)
             DEALLOCATE(sorted_fNeibPE)

             CALL DEBUG_BARRIER_PRINT("## aggregate             ##", SOLVER_COMM, my_rank)


             ALLOCATE(aggregate_table_array(aggregate_table_size * 2, max_neib_pe_size))
             ALLOCATE(aggregate_number_in_table(max_neib_pe_size))
             ALLOCATE(GIN_aggregate(0:NPROCS))

             CALL aggregate(PNI, NI, level, node_index, in_aggregates_result,           &
                  &         aggregates_result, in_aggregates_result_size, my_rank,      &
                  &         aggregate_table_size, aggregate_table_array,                &
                  &         aggregate_number_in_table, SOLVER_COMM, NPROCS,             &
                  &         GIN_aggregate, PE_num_size, PE_nums, nodes_for_each_PE,     &
                  &         in_nodes_for_each_PE, PE_list_size, PE_list, 2)

          END IF
          
       END IF

       CALL count_time(2, solver_comm, my_rank, 6)       

       DEALLOCATE(node_index)

       coarser_level_size = GIN_aggregate(my_rank + 1) - GIN_aggregate(my_rank)

       
       finish_flag = (coarser_level_size < MIN_NODE_SIZE) .OR. (level >= MAX_LEVEL_SIZE)       
       global_CN = 0


       CALL MPI_ALLREDUCE(coarser_level_size, global_CN, 1, LIS_MPI_INTEGER, MPI_SUM, SOLVER_COMM, ierr)


       CALL MPI_ALLREDUCE(finish_flag, global_finish_flag, 1, MPI_LOGICAL, MPI_LAND, SOLVER_COMM, ierr)


       RATE_OF_SPACE = RATE_OF_SPACE + 79 * (level - 2)
!!$       RATE_OF_SPACE = RATE_OF_SPACE + 40 * (level - 2)

       IF(RATE_OF_SPACE > global_CN) RATE_OF_SPACE = global_CN
       
       ALLOCATE(Temp_P(INT(RATE_OF_SPACE), finer_level_NP))
       ALLOCATE(Temp_P_row_size(finer_level_NP))


       CALL DEBUG_BARRIER_PRINT("## smooth_aggregate      ##", SOLVER_COMM, my_rank)



       CALL count_time(1, solver_comm, my_rank, 7)              
       if(SYMMETRIC_FLAG) then
          CALL smooth_aggregate_old (PNI(0:2*finer_level_NP),                                 &
               & NI(1:finer_level_NPL+finer_level_NPU+EX_NPL),                                &
               & dump_jacobi_weight, level,                                             &
               & in_aggregates_result_size, in_aggregates_result(0:in_aggregates_result_size),&
               & aggregates_result(1:in_aggregates_result(in_aggregates_result_size)),        &
               & Temp_P(1:INT(RATE_OF_SPACE),1:finer_level_NP),                               &
               & Temp_P_row_size(1:finer_level_NP), tilde_D(1:finer_level_NP),                &
               & INT(RATE_OF_SPACE),finer_level_NP, finer_level_NPL+finer_level_NPU+EX_NPL,   &
               & in_aggregates_result(in_aggregates_result_size))
       else
          
          CALL smooth_aggregate_unsym (PNI(0:2*finer_level_NP), &
               & NI(1:finer_level_NPL+finer_level_NPU+EX_NPL), &
               & dump_jacobi_weight, level, &
               & in_aggregates_result_size, in_aggregates_result(0:in_aggregates_result_size),&
               & aggregates_result, Temp_P, Temp_P_row_size, tilde_D,            &
               & INT(RATE_OF_SPACE), EX_INL, EX_IAL, EX_AL)
       end if
       

       CALL count_time(2, solver_comm, my_rank, 7)              
       DEALLOCATE(NI)
       DEALLOCATE(PNI)
       DEALLOCATE(tilde_D)
       DEALLOCATE(aggregates_result) 
       DEALLOCATE(in_aggregates_result)
       
       CALL mpi_barrier(SOLVER_COMM, ierr)

       
       CALL DEBUG_BARRIER_PRINT("## sr_aggregates         ##", SOLVER_COMM, my_rank)       

       ALLOCATE(global_local_hash_table(HASH_TABLE_SIZE, 2))
       
       IF(NPROCS > 1) THEN
          IF (N>0) THEN
          CALL count_time(1, solver_comm, my_rank, 2)       
          nfepe_size = in_nodes_for_each_PE(PE_num_size)          
          !C own aggregates is exchanged           
          CALL solver_SR_aggregates &
               &    ( PE_num_size, PE_nums(1:PE_num_size), nfepe_size, nodes_for_each_PE(1:nfepe_size),       &
               &      in_nodes_for_each_PE(0:PE_num_size), &
               &      SOLVER_COMM,  INT(RATE_OF_SPACE),finer_level_NP,&
               &      Temp_P(1:INT(RATE_OF_SPACE), 1:finer_level_NP), Temp_P_row_size(1:finer_level_NP), &
               &      my_rank, GIN_aggregate(0:NPROCS), NPROCS,                     &
               &      HASH_TABLE_SIZE, global_local_hash_table(1:HASH_TABLE_SIZE, 1:2), aggregate_table_size, &
               &      aggregate_table_array(1:aggregate_table_size*2, 1:max_neib_pe_size), &
               &      aggregate_number_in_table(1:max_neib_pe_size), &
               &      local_aggre_size, PE_list(1:NPROCS, 1:3), PE_list_size, level) 

          CALL count_time(2,solver_comm,my_rank,2)       
          DEALLOCATE(PE_nums)
          DEALLOCATE(nodes_for_each_PE)
          DEALLOCATE(in_nodes_for_each_PE)
          ELSE
              local_aggre_size=0
          END IF
       ELSE
          !C for single CPU case
          local_aggre_size = in_aggregates_result_size
       END IF
       
       !--C check Temp_P
       DO i = 1, finer_level_NP
          IF(Temp_P_row_size(i) > RATE_OF_SPACE) THEN
             WRITE(*,*) "RATE_OF_SPACE should be enlarged! "
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
             P % CN(k) = Temp_P(j, i) % column
             P % V(k)  = Temp_P(j, i) % value
             k = k + 1
          END DO
       END DO


       DEALLOCATE(Temp_P)
       DEALLOCATE(Temp_P_row_size)
       !C-- end:Temp_P is recorded in HIERARCHICAL_DATA()..


#ifdef DEBUG
       DO i=1,PE_list_size
          IF(aggregate_number_in_table(i)==0) THEN
             WRITE(*,*)  "0 size table", PE_list(i,1)," in " ,my_rank
          END IF
       END DO
#endif

       !C-- make R : 1..finer_level_NP -> 1..local_aggre_size
       !C   after RAP calculation,
       !C   later chanremake R:  1..finer_level_size -> 1..local_aggre_size
       ALLOCATE( HIERARCHICAL_DATA(level) % R )
       R => HIERARCHICAL_DATA(level) % R
       P => HIERARCHICAL_DATA(level) % P
       R % ROW_SIZE = local_aggre_size
       ALLOCATE( Temp_R_IN(0:local_aggre_size) )
       Temp_R_IN(0:local_aggre_size) = 0
       
       count = 0
       DO i = 1, finer_level_NP
          DO j = P % IN(i - 1) + 1, P % IN(i)
             col = P % CN(j)
             count = count + 1
             IF(col < local_aggre_size) Temp_R_IN(col + 1) = Temp_R_IN(col + 1) + 1
          END DO
       END DO
       
       ALLOCATE(R % CN(count))
       ALLOCATE(R % V(count))
    
       DO i = 2, R % ROW_SIZE
          Temp_R_IN(i) = Temp_R_IN(i) + Temp_R_IN(i - 1)
       END DO

       DO i = 1, finer_level_NP
          DO j = P % IN(i - 1) + 1, P % IN(i)
             col = P % CN(j)
             Temp_R_IN(col) = Temp_R_IN(col) + 1
             k = Temp_R_IN(col)
             R % CN(k) = i
             R %  V(k) = P % V(j)
          END DO
       END DO
       !C--end: make R

       CALL DEBUG_BARRIER_PRINT("## make_rap              ##", SOLVER_COMM, my_rank)       

       !C-- allocate temporal and calculate RAP.
!!$       TEMP_COLSIZE = NONZEROS_PER_ROW + 601 * (level - 2)
       TEMP_COLSIZE = NONZEROS_PER_ROW + 501 * (level - 2)
!!$       TEMP_COLSIZE = NONZEROS_PER_ROW
       IF(TEMP_COLSIZE > global_CN) TEMP_COLSIZE = global_CN
       

#if defined(CHOLESCKY)||defined(PSPASES)
       CALL make_rap(level, my_rank, SOLVER_COMM, NPROCS,                             &
            & coarser_level_size, global_finish_flag, GIN_aggregate(0:NPROCS),        &
            & PE_list_size, PE_list(1:NPROCS, 1:3), aggregate_table_size,             &
            & aggregate_table_array(1:aggregate_table_size*2,1:max_neib_pe_size),               &
            & aggregate_number_in_table(1:max_neib_pe_size), HASH_TABLE_SIZE,                   &
            & global_local_hash_table(1:HASH_TABLE_SIZE,1:2), local_aggre_size,       &
            & Temp_R_IN(0:local_aggre_size), &
            & TEMP_COLSIZE, TEMP_RAP_WORKSIZE)
#else
       finish_flag = .FALSE.
       CALL make_rap(level, my_rank, SOLVER_COMM, NPROCS,                             &
            & coarser_level_size, finish_flag, GIN_aggregate(0:NPROCS),               &
            & PE_list_size, PE_list(1:NPROCS,1:3), aggregate_table_size,                            &
            & aggregate_table_array(1:aggregate_table_size*2,1:max_neib_pe_size),               &
            & aggregate_number_in_table(1:max_neib_pe_size), HASH_TABLE_SIZE,                   &
            & global_local_hash_table(1:HASH_TABLE_SIZE,1:2), local_aggre_size,       &
            & Temp_R_IN(0:local_aggre_size), TEMP_COLSIZE, TEMP_RAP_WORKSIZE)
#endif

       DEALLOCATE(Temp_R_IN)
       DEALLOCATE(global_local_hash_table)
       DEALLOCATE(aggregate_table_array)
       DEALLOCATE(aggregate_number_in_table)
       DEALLOCATE(GIN_aggregate)

       

#ifdef PRINT_DATA_CREATE
       WRITE(*,*) '##', HIERARCHICAL_DATA(level)%N, HIERARCHICAL_DATA(level)%NP,"@",level, my_rank, '##'
       IF(my_rank==0)THEN
          WRITE(*,*) '############',level,'############'
          
          DO i = 1, HIERARCHICAL_DATA(level)%N
             write(*,*) "rowrowrowrow",i,"rowrowrowrow"
             DO j = HIERARCHICAL_DATA(level)%INU(i-1)+1, HIERARCHICAL_DATA(level)%INU(i)
                write(*,*) HIERARCHICAL_DATA(level)%IAU(j),HIERARCHICAL_DATA(level)%AU(j)
             END DO
             write(*,*) i, HIERARCHICAL_DATA(level)%D(i)
             DO j = HIERARCHICAL_DATA(level)%INL(i-1)+1, HIERARCHICAL_DATA(level)%INL(i)
                write(*,*) HIERARCHICAL_DATA(level)%IAL(j),HIERARCHICAL_DATA(level)%AL(j)
             END DO
          END DO
       END IF
#endif
       
       if(.not.SYMMETRIC_FLAG) then       
          IF (ASSOCIATED(EX_INL)) DEALLOCATE(EX_INL)
          IF (ASSOCIATED(EX_IAL)) DEALLOCATE(EX_IAL)
          IF (ASSOCIATED(EX_AL)) DEALLOCATE(EX_AL)
       end if
       
       IF(global_finish_flag)   EXIT
       
       if(max_neib_pe_size < 6*finerNEIBPETOT) then
          max_neib_pe_size = 6*finerNEIBPETOT
       end if
       if(max_neib_pe_size > NPROCS) max_neib_pe_size = NPROCS
       

    END DO
    DEALLOCATE(PE_list)
    
    LEVEL_NUM = level

    
  END SUBROUTINE data_creation

  
!==========================================================================
  SUBROUTINE repeated_nds_tbl_crtn(LEVEL_NO,sorted_NEIBPE,SOLVER_COMM,my_rank, &
       &       PE_num_size,PE_nums,nodes_for_each_PE,in_nodes_for_each_PE,NPROCS)
    USE data_structure_for_AMG
    USE solver_Gnumbering
    USE count_time_mod
    IMPLICIT NONE

    INTEGER(kind=kint ) :: N,NEIBPETOT
    
    INTEGER(kind=kint ), INTENT(in) :: sorted_NEIBPE(:)
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO,NPROCS
    INTEGER(kind=kint ), INTENT(in) :: SOLVER_COMM
    INTEGER(kind=kint ), INTENT(in) :: my_rank

    INTEGER(kind=kint ), POINTER    :: NEIBPE      (:)
    INTEGER(kind=kint ), POINTER    :: NOD_IMPORT  (:), STACK_IMPORT(:)
    INTEGER(kind=kint ), POINTER    :: NOD_EXPORT  (:), STACK_EXPORT(:)
    
    INTEGER(kind=kint ), ALLOCATABLE:: export_node_list(:,:)
    INTEGER(kind=kint ), ALLOCATABLE:: repeated_nodes_at_boundary(:,:)

    INTEGER(kind=kint ) :: i,j,k,l,js,je,count
    INTEGER(kind=kint ) :: current_node,repeated_nodes_tot

    INTEGER(kind=kint ) :: REPEATED_NODE_SIZE
    INTEGER(kind=kint ), ALLOCATABLE   :: export_repeated(:,:)
    INTEGER(kind=kint ), ALLOCATABLE   :: in_export_repeated(:)
    INTEGER(kind=kint ), ALLOCATABLE   :: repeated_nodes_send_buf(:)

    INTEGER(kind=kint ),INTENT(out) :: PE_num_size
    INTEGER(kind=kint ),  POINTER :: PE_nums(:)
    INTEGER(kind=kint ),  POINTER :: nodes_for_each_PE(:)
    INTEGER(kind=kint ),  POINTER :: in_nodes_for_each_PE(:)

    N            = HIERARCHICAL_DATA(LEVEL_NO-1) % N
    NEIBPETOT    =  HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPETOT     
    if(NEIBPETOT > 0) then
       NEIBPE       => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPE
       STACK_IMPORT => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_IMPORT
       STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_EXPORT
       i = STACK_IMPORT(NEIBPETOT)
       NOD_IMPORT   => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_IMPORT
       i = STACK_EXPORT(NEIBPETOT)
       NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_EXPORT
    else
       if (my_rank==0) then
          write(*,*) "domain decomposition failed in SAAMG"
          write(*,*) "0 size NEIBPETOT at LEVEL_NO:", LEVEL_NO
          stop 'error in repeated_nds_tbl_crtn()'
       end if
    end if

    
    
    k=STACK_EXPORT(NEIBPETOT)
    ALLOCATE(export_node_list(3,k))
    
    !C export_node_list 1: node 2: PE NO 3: next
    DO i=1,k
       export_node_list(3,i)=0
    END DO
    ALLOCATE(repeated_nodes_at_boundary(3,N))

    !C repeated_nodes_at_boundary 1:TAIL 2:SIZE,3:HEAD
    DO i=1,N
       repeated_nodes_at_boundary(1,i)=0
       repeated_nodes_at_boundary(2,i)=0
    END DO


    DO i = 1, NEIBPETOT
       js = STACK_EXPORT(i - 1) + 1
       je = STACK_EXPORT(i)
       DO j = js, je
          k = NOD_EXPORT(j)+ZERO_ORIGIN

          repeated_nodes_at_boundary(2, k) = repeated_nodes_at_boundary(2, k) + 1
          current_node = repeated_nodes_at_boundary(1, k)
          IF(current_node == 0) THEN
             repeated_nodes_at_boundary(3, k) = j
          ELSE
             export_node_list(3, current_node) = j
          END IF
          repeated_nodes_at_boundary(1, k) = j
          !C node is the number from  STACK_EXPORT's associated PE js
          export_node_list(1, j) = j - js + 1
          export_node_list(2, j) = i
       END DO
    END DO


    !C-- REPEATED_NODE_SIZE is max size of send table
    REPEATED_NODE_SIZE = 0
    DO i = 1, NEIBPETOT
       j = STACK_EXPORT(i)-STACK_EXPORT(i-1)
       IF(REPEATED_NODE_SIZE < j) REPEATED_NODE_SIZE = j
    END DO
    !C--end: REPEATED_NODE_SIZE is calculated.

    ALLOCATE(export_repeated(REPEATED_NODE_SIZE, NEIBPETOT))
    ALLOCATE(in_export_repeated(0:NEIBPETOT))
    
    repeated_nodes_tot = 0
    in_export_repeated(0:NEIBPETOT) = 0

    DO i = 1, N
       l = repeated_nodes_at_boundary(2, i)
       IF(l >= 2) THEN
          repeated_nodes_tot = repeated_nodes_tot + l * (l + 1)
          !C l:number of PEs  *  l+1:node_NO,PE_list_size,PE_list  
          current_node = repeated_nodes_at_boundary(3, i)
          DO WHILE(current_node > 0)
             k = export_node_list(2, current_node)
             in_export_repeated(k) = in_export_repeated(k) + 1
             j = in_export_repeated(k)
             !C node_number
             export_repeated(j,k) = i
             current_node=export_node_list(3, current_node)
          END DO
       END IF
    END DO

    !C No nodes are repeated.-> assign 0 
    DO i = 1, NEIBPETOT
       IF(in_export_repeated(i) == 0)THEN
          repeated_nodes_tot = repeated_nodes_tot + 1
       END IF
    END DO
    
    ALLOCATE(repeated_nodes_send_buf(repeated_nodes_tot))
    count = 0
    DO i = 1, NEIBPETOT
       IF(in_export_repeated(i) == 0) THEN
          !C No nodes are repeated for that PE.
          count = count + 1
          repeated_nodes_send_buf(count) = 0
       ELSE
          DO j = 1, in_export_repeated(i)
             k = export_repeated(j, i)
             current_node = repeated_nodes_at_boundary(3, k)
             count = count + 1
             repeated_nodes_send_buf(count) = repeated_nodes_at_boundary(2, k)
             
             count = count + 1
             !C l is the place for node number in NOD_EXPORT.
             l = count
             DO WHILE(current_node > 0)
                !C PE in NEIBPE(i) is in the stack_export.
                !C export_node_list:: 1:node  2:pe  3:next
                IF(export_node_list(2, current_node) == i) THEN
                   repeated_nodes_send_buf(l) = export_node_list(1, current_node) 
                ELSE
                   count = count + 1
                   repeated_nodes_send_buf(count) = NEIBPE(export_node_list(2, current_node))
#ifdef DEBUG
                   IF(NEIBPE(export_node_list(2, current_node)) > NPROCS) THEN
                      STOP "repeated_nodes creation send buffer"
                   END IF
#endif
                END IF
                current_node = export_node_list(3, current_node)
             END DO
          END DO
       END IF
       !C in_export_repeated is changed to the index for each PE in 
       !C repeated_nodes_send_buf.
       in_export_repeated(i)=count
    END DO
    in_export_repeated(0)=0
    
    !C deallocate repeated_nodes_at_boundary
    DEALLOCATE(repeated_nodes_at_boundary)
    DEALLOCATE(export_node_list)
    DEALLOCATE(export_repeated)
    
    CALL count_time(1,solver_comm,my_rank,1)    
    !C repeated nodes information is exchanged
    !C PE_num_size,PE_nums,in_nodes_for_each_PE,nodes_for_each_PE is created


    CALL solver_SR_repeated_nodes &
         &    ( repeated_nodes_tot,     &
         &      NEIBPETOT, NEIBPE(1:NEIBPETOT), STACK_IMPORT(0:NEIBPETOT),&
         &      NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)), STACK_EXPORT(0:NEIBPETOT), &
         &      NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)), &
         &      in_export_repeated(0:NEIBPETOT),&
         &      repeated_nodes_send_buf(1:repeated_nodes_tot),   &
         &      PE_num_size,PE_nums,nodes_for_each_PE,in_nodes_for_each_PE,SOLVER_COMM, &
         &      sorted_NEIBPE(1:NEIBPETOT), my_rank,NPROCS) 

    CALL count_time(2,solver_comm,my_rank,1)    
    
    DEALLOCATE(repeated_nodes_send_buf)
    DEALLOCATE(in_export_repeated)

  CONTAINS
    SUBROUTINE  solver_SR_repeated_nodes &
         &     ( SIZE, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,STACK_EXPORT,&
         &       NOD_EXPORT, in_export_repeated,repeated_nodes_send_buf,                &
         &       PE_num_size,PE_nums,nodes_for_each_PE,in_nodes_for_each_PE,SOLVER_COMM,&
         &       sorted_recv_order_rank,my_rank,NPROCS)

      USE isort
      IMPLICIT NONE
      INCLUDE  'mpif.h'
#include "precision.inc"          

      INTEGER(kind=kint ), INTENT(in) :: SIZE,NEIBPETOT,NPROCS
      INTEGER(kind=kint ), INTENT(in) :: NEIBPE(:)
      INTEGER(kind=kint ), INTENT(in) :: STACK_IMPORT(0:), STACK_EXPORT(0:)
      INTEGER(kind=kint ), INTENT(in) :: NOD_IMPORT(:), NOD_EXPORT(:)
      INTEGER(kind=kint ), INTENT(in) :: sorted_recv_order_rank(:)

      INTEGER(kind=kint ), INTENT(in)  :: in_export_repeated(0:)
      INTEGER(kind=kint ), INTENT(in)  :: repeated_nodes_send_buf(SIZE)
      INTEGER(kind=kint ), INTENT(in)  :: SOLVER_COMM
      INTEGER(kind=kint ), INTENT(in)  :: my_rank

      INTEGER(kind=kint ),  ALLOCATABLE :: status(:)
      INTEGER(kind=kint ),  ALLOCATABLE :: sta1(:,:)
      INTEGER(kind=kint ),  ALLOCATABLE :: sta2(:,:)
      INTEGER(kind=kint ),  ALLOCATABLE :: req1(:)
      INTEGER(kind=kint ),  ALLOCATABLE :: req2(:)

      INTEGER(kind=kint ) :: neib,istart,inum,message_size,buf_size,count,w,s,e,hn,hc,hs
      INTEGER(kind=kint ) :: i,j,k,l,m,counter_PE_nums,counter_nodes_for_each_PE,base
      INTEGER(kind=kint ) :: PE_list_size,src_PE,PE_NO,Local_node_NO
      INTEGER(kind=kint ) :: ierr

      INTEGER(kind=kint), ALLOCATABLE :: buf(:)
      INTEGER(kind=kint), ALLOCATABLE :: cnt_nd_for_echPE(:)
      INTEGER(kind=kint), ALLOCATABLE :: in_recv_order_rank(:)

      INTEGER(kind=kint ), INTENT(out) :: PE_num_size
      INTEGER(kind=kint ), POINTER :: PE_nums(:)
      INTEGER(kind=kint ), POINTER :: nodes_for_each_PE(:)
      INTEGER(kind=kint ), POINTER :: in_nodes_for_each_PE(:)


      ALLOCATE (status(MPI_STATUS_SIZE))
      ALLOCATE (sta1(MPI_STATUS_SIZE,NEIBPETOT))
      ALLOCATE (sta2(MPI_STATUS_SIZE,NEIBPETOT))
      ALLOCATE (req1(NEIBPETOT))
      ALLOCATE (req2(NEIBPETOT))      
      ALLOCATE (in_recv_order_rank(0:NEIBPETOT))

      !C-- message size is communicated
      DO neib= 1, NEIBPETOT
         istart= in_export_repeated(neib-1)
         inum  = in_export_repeated(neib  ) - istart
         CALL MPI_ISEND (inum, 1, LIS_MPI_INTEGER, NEIBPE(neib), 0, &
              &          SOLVER_COMM, req1(neib), ierr)
      END DO
      DO neib= 1, NEIBPETOT      
         CALL MPI_IRECV(in_recv_order_rank(neib),1,LIS_MPI_INTEGER, &
              &         NEIBPE(neib), 0, SOLVER_COMM, req2(neib), ierr)
      END DO
      CALL MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
      CALL MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)

      in_recv_order_rank(0) = 0
      DO i=1,NEIBPETOT
         in_recv_order_rank(i)=in_recv_order_rank(i)+in_recv_order_rank(i-1)
      END DO

      ALLOCATE( buf( in_recv_order_rank(NEIBPETOT) ) )
      buf = 0
      !C-- end: message size is communicated


      !C-- repeated node and PE numbers are communicated
      DO neib= 1, NEIBPETOT
         istart= in_export_repeated(neib-1)
         inum  = in_export_repeated(neib  ) - istart
         CALL MPI_ISEND (repeated_nodes_send_buf(istart+1), inum, LIS_MPI_INTEGER,     &
              &                NEIBPE(neib), 0, SOLVER_COMM,                       &
              &                req1(neib), ierr)
      END DO
      DO neib= 1, NEIBPETOT      
         istart=in_recv_order_rank(neib-1)
         inum=in_recv_order_rank(neib)-istart
         CALL MPI_IRECV(buf(istart+1),inum,LIS_MPI_INTEGER,NEIBPE(neib),0,SOLVER_COMM, &
              &         req2(neib),ierr)
      END DO

      CALL MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
      CALL MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)
      !C--end: repeated node and PE numbers are communicated

      !C-- determine node sized for all PEs on cnt_nd_for_echPE
      ALLOCATE(cnt_nd_for_echPE(0:NPROCS-1))
      cnt_nd_for_echPE=0

      DO i = 1, NEIBPETOT
         s = in_recv_order_rank(i - 1) + 1
         e = in_recv_order_rank(i)
         k = s
         IF(buf(k) > 0) THEN
            DO WHILE(k < e)
               l = buf(k)
               DO j = k + 2, k + l
                  PE_NO = buf(j)
#ifdef DEBUG
                  IF(PE_NO >= NPROCS) THEN 
                     WRITE(*,*) "PE_NO:", PE_NO, s,e ,j-k,l, my_rank, NEIBPE(i) 
                     STOP "PE_NO in repeated tables creation"
                  END IF
#endif
                  cnt_nd_for_echPE(PE_NO) = cnt_nd_for_echPE(PE_NO) + 1
               END DO
               k = k + l + 1
            END DO
         END IF
      END DO
      !C--end: determine node size for all PEs on cnt_nd_for_echPE


      !C--PE_num_size and in_nodes_for_each_PE are created from cnt_nd_for_echPE
      !C--cnt_nd_for_echPE is used for transpose PE_nums in lateral part.
      PE_num_size=0
      DO i=0,NPROCS-1
         IF(cnt_nd_for_echPE(i)>0) PE_num_size=PE_num_size+1
      END DO
      ALLOCATE(PE_nums(PE_num_size))
      ALLOCATE(in_nodes_for_each_PE(0:PE_num_size+1))

      in_nodes_for_each_PE(0)=0
      in_nodes_for_each_PE(1)=0
      j=0
      DO i=0,NPROCS-1
         IF(cnt_nd_for_echPE(i)>0) THEN
            j=j+1
            PE_nums(j)=i
            in_nodes_for_each_PE(j+1)=in_nodes_for_each_PE(j)+cnt_nd_for_echPE(i)

            !C PE_nums index can be retrieved from PE number 
            cnt_nd_for_echPE(i)=j
         END IF
      END DO
      !C--end: PE_num_size and in_nodes_for_each_PE are created from cnt_nd_for_echPE    


      !C--nodes_for_each_PE is determined
      ALLOCATE(nodes_for_each_PE(in_nodes_for_each_PE(PE_num_size+1)))
      DO i=1,NEIBPETOT
         j=sorted_recv_order_rank(i)                
         base=STACK_IMPORT(j-1)
         s=in_recv_order_rank(j-1)+1
         e=in_recv_order_rank(j)
         k=s
         IF(buf(k)>0) THEN
            DO WHILE(k < e)
               PE_list_size=buf(k)
               Local_node_NO=NOD_IMPORT(buf(k+1)+base)+ZERO_ORIGIN
               DO j=k+2,k+PE_list_size
                  PE_NO=buf(j)
                  l=cnt_nd_for_echPE(PE_NO)

                  m=in_nodes_for_each_PE(l)+1
                  in_nodes_for_each_PE(l)=m
                  nodes_for_each_PE(m)=Local_node_NO-ZERO_ORIGIN

               END DO
               k=k+PE_list_size+1
            END DO
         END IF
      END DO
      DEALLOCATE(cnt_nd_for_echPE)
      !C--end: nodes_for_each_PE is determined

      DEALLOCATE (in_recv_order_rank)
      DEALLOCATE (buf)

      DEALLOCATE (status)
      DEALLOCATE (sta1)
      DEALLOCATE (sta2)
      DEALLOCATE (req1)
      DEALLOCATE (req2)
    END SUBROUTINE solver_SR_repeated_nodes

  END SUBROUTINE repeated_nds_tbl_crtn
  
  
  SUBROUTINE make_PE_lists(LEVEL_NO, sorted_NEIBPE, PE_num_size, PE_nums, &
       &                   PE_list_size, PE_list, NPROCS)
    USE data_structure_for_AMG
    USE solver_Gnumbering
    IMPLICIT NONE
    INCLUDE  'mpif.h'
    INTEGER(kind=kint), INTENT(in)    :: LEVEL_NO, NPROCS
    INTEGER(kind=kint), INTENT(in)    :: sorted_NEIBPE(:)
    INTEGER(kind=kint), INTENT(in)    :: PE_num_size
    INTEGER(kind=kint ), POINTER :: PE_nums(:)
    INTEGER(kind=kint ), INTENT(out) :: PE_list_size
    INTEGER(kind=kint ), INTENT(out) :: PE_list(NPROCS, 3)
    
    INTEGER(kind=kint )          :: NEIBPETOT
    INTEGER(kind=kint ), POINTER :: NEIBPE(:)
    
    INTEGER(kind=kint) :: i, j, k, l, pe, node
    INTEGER(kind=kint) :: index, Neib_index, Pe_num_index
    
    NEIBPE    => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NEIBPE
    NEIBPETOT =  HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NEIBPETOT     
    
    
    !C PE_list(*,1): PE number
    !C PE_list(*,2): NEIBPE number|
    !C PE_list(*,3): PE_num number|these two has 0 if no correspondings.
    !C-- making PE_list..
    index = 1
    Neib_index = 1
    Pe_num_index = 1

    DO WHILE((Neib_index <= NEIBPETOT) .AND. (Pe_num_index <= PE_num_size))
       i = sorted_NEIBPE(Neib_index)
       IF(PE_nums(Pe_num_index) == NEIBPE(i)) THEN
          PE_list(index, 1) = PE_nums(Pe_num_index)
          PE_list(index, 2) = i
          PE_list(index, 3) = Pe_num_index
          Neib_index = Neib_index + 1
          PE_num_index = Pe_num_index + 1
       ELSE IF(PE_nums(Pe_num_index) < NEIBPE(i)) THEN
          PE_list(index, 1) = PE_nums(Pe_num_index)
          PE_list(index, 2) = 0
          PE_list(index, 3) = Pe_num_index
          PE_num_index = PE_num_index + 1
       ELSE      
          PE_list(index, 1) = NEIBPE(i)
          PE_list(index, 2) = i
          PE_list(index, 3) = 0
          Neib_index = Neib_index + 1
       END IF
       index = index + 1
    END DO
    
    IF(Neib_index > NEIBPETOT) THEN
       DO i = PE_num_index, PE_num_size
          PE_list(index, 1) = PE_nums(i)
          PE_list(index, 2) = 0
          PE_list(index, 3) = i
          index = index + 1
       END DO
    ELSE
       DO i = Neib_index, NEIBPETOT
          j = sorted_NEIBPE(i)
          PE_list(index, 1) = NEIBPE(j)
          PE_list(index, 2) = j
          PE_list(index, 3) = 0
          index = index + 1
       END DO
    END IF
    PE_list_size = index - 1
    !C-- end: making PE_list..
  END SUBROUTINE make_PE_lists
  
  
!==========================================================================
  SUBROUTINE RAP(Temp_N, Temp_CN, Temp_V, hash_size, s, ownaggre_size, LEVEL_NO, local_aggre_size, my_rank)
    
    USE data_structure_for_AMG
    USE hash_mod
    IMPLICIT NONE

    INTEGER(kind=kint),INTENT(in)   :: hash_size,s,ownaggre_size,LEVEL_NO,local_aggre_size,my_rank
    INTEGER(kind=kint),INTENT(inout):: Temp_N(1:local_aggre_size)
    INTEGER(kind=kint),INTENT(inout):: Temp_CN(1:hash_size,1:local_aggre_size)
    TYPE(INTER_LEVEL_OPERATOR),POINTER :: P
    REAL(kind=kreal),INTENT(inout):: Temp_V(1:hash_size,1:local_aggre_size)

    !C pointers for matrix
    REAL   (kind=kreal),POINTER ::  D(:)
    REAL   (kind=kreal),POINTER ::  AU(:)
    REAL   (kind=kreal),POINTER ::  AL(:)
    INTEGER(kind=kint ),POINTER ::  INU(:)
    INTEGER(kind=kint ),POINTER ::  IAU(:)
    INTEGER(kind=kint ),POINTER ::  INL(:)
    INTEGER(kind=kint ),POINTER ::  IAL(:)

    INTEGER(kind=kint) ::i,j,k,rowins,rowine,colins,coline,row,col
    INTEGER(kind=kint) ::js,je,jin,rowin,colin,count,hash_no
    LOGICAL :: flag

    D    => HIERARCHICAL_DATA(LEVEL_NO-1) % D
    INL  => HIERARCHICAL_DATA(LEVEL_NO-1) % INL
    INU  => HIERARCHICAL_DATA(LEVEL_NO-1) % INU
    IAL  => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU  => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL   => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU   => HIERARCHICAL_DATA(LEVEL_NO-1) % AU
    P    => HIERARCHICAL_DATA(LEVEL_NO  ) % P



    !C Calculation of RAP. 
    !C ------
    DO i = 1, s
       rowins = P % IN(i - 1) + 1
       rowine = P % IN(i)

       !C A's diagonal part
       DO rowin = rowins, rowine
          row = P % CN(rowin)
          colins = rowins
          coline = rowine
          DO colin = colins, coline
             col = P % CN(colin)

             IF(row <= ownaggre_size .AND. col < row) CYCLE

             !C search (row col) in new matrix
             CALL hash(hash_size, Temp_CN(1:hash_size, row), col, count, "NONZEROS_PER_ROW")
             
             !C assign the value
             IF(Temp_CN(count, row) == col) THEN
                Temp_V(count, row) = Temp_V(count,row) + P % V(rowin) * D(i) * P % V(colin)
             ELSE
                Temp_N(row) = Temp_N(row) + 1
                Temp_CN(count, row) = col
                Temp_V(count, row) = P % V(rowin) * D(i) * P % V(colin)
             END IF
          END DO
       END DO

       !C A's lower part
       js = INL(i - 1) + 1
       je = INL(i)
       DO jin = js, je
          j = IAL(jin)+ZERO_ORIGIN

          DO rowin = rowins, rowine
             row = P % CN(rowin)
             colins = P % IN(j - 1) + 1
             coline = P % IN(j)
             DO colin = colins, coline
                col = P % CN(colin)
                IF(row <= ownaggre_size .AND. col < row) CYCLE


                !C search (row col) in new matrix
                CALL hash(hash_size, Temp_CN(1:hash_size, row), col, count, "NONZEROS_PER_ROW")

                !C assign the value
                IF(Temp_CN(count, row) == col) THEN
                   Temp_V(count, row) = Temp_V(count, row) + P % V(rowin) * AL(jin) * P % V(colin)
                ELSE
                   Temp_N(row) = Temp_N(row) + 1
                   Temp_CN(count, row) = col
                   Temp_V(count, row) = P % V(rowin) * AL(jin) * P % V(colin)
                END IF
             END DO
          END DO
       END DO

       !C A's upper part      
       js = INU(i - 1) + 1
       je = INU(i)
       DO jin = js, je
          j = IAU(jin)+ZERO_ORIGIN

          DO rowin = rowins, rowine
             row = P % CN(rowin)

             colins = P % IN(j - 1) + 1
             coline = P % IN(j)
             DO colin = colins, coline
                col = P % CN(colin)

                IF(row <= ownaggre_size .AND. col < row) CYCLE

                !C search (row col) in new matrix
                CALL hash(hash_size, Temp_CN(1:hash_size, row), col, count, "NONZEROS_PER_ROW")

                !C assign the value
                IF(Temp_CN(count, row) == col) THEN
                   Temp_V(count, row) = Temp_V(count, row) + P % V(rowin) * AU(jin) * P % V(colin)
                ELSE
                   Temp_N(row) = Temp_N(row) + 1
                   Temp_CN(count, row) = col
                   Temp_V(count, row) = P % V(rowin) * AU(jin) * P % V(colin)
                END IF
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE RAP


!==========================================================================
  SUBROUTINE RAP_unsym(Temp_N, Temp_CN, Temp_V, hash_size, s, ownaggre_size, LEVEL_NO, local_aggre_size)
    
    USE data_structure_for_AMG
    USE hash_mod
    IMPLICIT NONE

    INTEGER(kind=kint), INTENT(in)   :: hash_size,s,ownaggre_size,LEVEL_NO,local_aggre_size
    INTEGER(kind=kint), INTENT(inout):: Temp_N(1:local_aggre_size)
    INTEGER(kind=kint), INTENT(inout):: Temp_CN(1:hash_size,1:local_aggre_size)
    TYPE(INTER_LEVEL_OPERATOR), POINTER :: P
    REAL(kind=kreal), INTENT(inout) :: Temp_V(1:hash_size,1:local_aggre_size)
    
    !C pointers for matrix
    REAL   (kind=kreal), POINTER ::  D(:)
    REAL   (kind=kreal), POINTER ::  AU(:)
    REAL   (kind=kreal), POINTER ::  AL(:)
    INTEGER(kind=kint ), POINTER ::  INU(:)
    INTEGER(kind=kint ), POINTER ::  IAU(:)
    INTEGER(kind=kint ), POINTER ::  INL(:)
    INTEGER(kind=kint ), POINTER ::  IAL(:)

    INTEGER(kind=kint) ::i,j,k,rowins,rowine,colins,coline,row,col
    INTEGER(kind=kint) ::js,je,jin,rowin,colin,count,hash_no
    LOGICAL :: flag
    
    D    => HIERARCHICAL_DATA(LEVEL_NO-1) % D
    INL  => HIERARCHICAL_DATA(LEVEL_NO-1) % INL
    INU  => HIERARCHICAL_DATA(LEVEL_NO-1) % INU
    IAL  => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU  => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL   => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU   => HIERARCHICAL_DATA(LEVEL_NO-1) % AU
    P    => HIERARCHICAL_DATA(LEVEL_NO  ) % P

    !C Calculation of RAP. 
    !C ------
    DO i = 1, s
       rowins = P % IN(i - 1) + 1
       rowine = P % IN(i)

       !C A's diagonal part
       DO rowin = rowins, rowine
          row = P % CN(rowin)
          colins = rowins
          coline = rowine
          DO colin = colins, coline
             col = P % CN(colin)
!!$             IF(row <= ownaggre_size) CYCLE

             !C search (row col) in new matrix
             CALL hash(hash_size, Temp_CN(:, row), col, count, "NONZEROS_PER_ROW")
             
             !C assign the value
             IF(Temp_CN(count, row) == col) THEN
                Temp_V(count, row) = Temp_V(count,row) + P % V(rowin) * D(i) * P % V(colin)
             ELSE
                Temp_N(row) = Temp_N(row) + 1
                Temp_CN(count, row) = col
                Temp_V(count, row) = P % V(rowin) * D(i) * P % V(colin)
             END IF
          END DO
       END DO

       !C A's lower part
       js = INL(i - 1) + 1
       je = INL(i)
       DO jin = js, je
          j = IAL(jin)+ZERO_ORIGIN

          DO rowin = rowins, rowine
             row = P % CN(rowin)
             colins = P % IN(j - 1) + 1
             coline = P % IN(j)
             DO colin = colins, coline
                col = P % CN(colin)
!!$                IF(row <= ownaggre_size) CYCLE

                !C search (row col) in new matrix
                CALL hash(hash_size, Temp_CN(:, row), col, count, "NONZEROS_PER_ROW")

                !C assign the value
                IF(Temp_CN(count, row) == col) THEN
                   Temp_V(count, row) = Temp_V(count, row) + P % V(rowin) * AL(jin) * P % V(colin)
                ELSE
                   Temp_N(row) = Temp_N(row) + 1
                   Temp_CN(count, row) = col
                   Temp_V(count, row) = P % V(rowin) * AL(jin) * P % V(colin)
                END IF
             END DO
          END DO
       END DO

       !C A's upper part      
       js = INU(i - 1) + 1
       je = INU(i)
       DO jin = js, je
          j = IAU(jin)+ZERO_ORIGIN

          DO rowin = rowins, rowine
             row = P % CN(rowin)

             colins = P % IN(j - 1) + 1
             coline = P % IN(j)
             DO colin = colins, coline
                col = P % CN(colin)

!!$                IF(row <= ownaggre_size) CYCLE

                !C search (row col) in new matrix
                CALL hash(hash_size, Temp_CN(:, row), col, count, "NONZEROS_PER_ROW")
!!$                !C debug
!!$                if(P % V(rowin) * AU(jin) * P % V(colin) > 10D+10) write(*,*) "here!!",P % V(rowin),AU(jin),P % V(colin)
                !C assign the value
                IF(Temp_CN(count, row) == col) THEN
                   Temp_V(count, row) = Temp_V(count, row) + P % V(rowin) * AU(jin) * P % V(colin)
                ELSE
                   Temp_N(row) = Temp_N(row) + 1
                   Temp_CN(count, row) = col
                   Temp_V(count, row) = P % V(rowin) * AU(jin) * P % V(colin)
                END IF
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE RAP_UNSYM


  
  !C R(AP) and this routine needs the rows from N+1 to NP in L part.
  SUBROUTINE RAP2(Temp_N, Temp_CN, Temp_V, hash_size, N, ownaggre_size, LEVEL_NO, &
       & local_aggre_size, Temp_R_IN, Temp_SIZE)

    USE data_structure_for_AMG
    USE hash_mod
    IMPLICIT NONE

    INTEGER(kind=kint),INTENT(in) :: hash_size, N, ownaggre_size, LEVEL_NO, local_aggre_size
    INTEGER(kind=kint),INTENT(in) :: Temp_SIZE
    INTEGER(kind=kint),INTENT(in)   :: Temp_R_IN(0:)
    INTEGER(kind=kint),INTENT(inout):: Temp_N(local_aggre_size)
    INTEGER(kind=kint),INTENT(inout):: Temp_CN(hash_size,local_aggre_size)
    TYPE(INTER_LEVEL_OPERATOR),POINTER :: P, R
    REAL(kind=kreal),INTENT(inout):: Temp_V(hash_size,local_aggre_size)

    !C pointers for matrix
    REAL   (kind=kreal), POINTER ::  D(:)
    REAL   (kind=kreal), POINTER ::  AU(:)
    REAL   (kind=kreal), POINTER ::  AL(:)
    INTEGER(kind=kint ), POINTER ::  INU(:)
    INTEGER(kind=kint ), POINTER ::  IAU(:)
    INTEGER(kind=kint ), POINTER ::  INL(:)
    INTEGER(kind=kint ), POINTER ::  IAL(:)
    
    INTEGER(kind=kint ),ALLOCATABLE :: Temp_agr(:),Temp_sz(:),Temp_nd(:)
    INTEGER(kind=kint ),ALLOCATABLE :: Temp_int(:,:)
    REAL   (kind=kreal),ALLOCATABLE :: Temp(:,:),Temp_1_col(:)

    INTEGER(kind=kint) :: i,j,k,l,m,rowins,rowine,colins,coline,row,col,row2
    INTEGER(kind=kint) :: do_whl_cnt,is,ie
    INTEGER(kind=kint) :: js,je,jin,rowin,colin,count,hash_no,NP
    LOGICAL :: cancel_flag,finish_flag

    INTEGER(kind=kint ) :: in1,in2,P_col,T_col,p_width,t_width,cancel_sz,nd_sz
    REAL   (kind=kreal) :: v

    NP   =  HIERARCHICAL_DATA(LEVEL_NO-1) % NP 
    D    => HIERARCHICAL_DATA(LEVEL_NO-1) % D 
    INL  => HIERARCHICAL_DATA(LEVEL_NO-1) % INL
    INU  => HIERARCHICAL_DATA(LEVEL_NO-1) % INU
    IAL  => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU  => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL   => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU   => HIERARCHICAL_DATA(LEVEL_NO-1) % AU
    P    => HIERARCHICAL_DATA(LEVEL_NO  ) % P
    R    => HIERARCHICAL_DATA(LEVEL_NO  ) % R

    !C Calculation of RAP. This uses R information. 
    !C R % IN must be rewrite @ the latter part.
    !C So Temp_R_IN is used for R % IN(:).
    !C ------

    ALLOCATE(Temp_sz   (N)               )
    ALLOCATE(Temp_int  (Temp_SIZE, N)    )
    ALLOCATE(Temp      (Temp_SIZE, N)    )
    ALLOCATE(Temp_agr  (local_aggre_size))
    ALLOCATE(Temp_1_col(local_aggre_size))
    ALLOCATE(Temp_nd   (local_aggre_size))
       
    Temp_1_col = 0    
    
    Temp_agr = 0
    !C-- RA
    finish_flag = .TRUE.
    do_whl_cnt = 0    
    DO WHILE(finish_flag)
       finish_flag = .FALSE.
       do_whl_cnt = do_whl_cnt + 1
       Temp_sz = 0
       
       DO col = 1, local_aggre_size
          IF(Temp_agr(col) > 0) CYCLE
          cancel_flag = .FALSE.
          
          DO i = Temp_R_IN(col - 1) + 1, Temp_R_IN(col) 
             row = R % CN(i)
             IF(row <= N .AND. Temp_sz(row) >= Temp_SIZE) THEN
                cancel_flag = .TRUE.
                EXIT
             END IF
          END DO
          IF(.NOT.cancel_flag) THEN
             DO i = Temp_R_IN(col - 1) + 1, Temp_R_IN(col) 
                row = R % CN(i)
                DO j = INL(row - 1) + 1, INL(row)
                   row2 = IAL( j ) + ZERO_ORIGIN
                   IF(row2 <= N .AND. Temp_sz(row2) >= Temp_SIZE) THEN
                      cancel_flag = .TRUE.
                      EXIT
                   END IF
                END DO
                IF(cancel_flag) EXIT
                DO j = INU(row - 1) + 1, INU(row)
                   row2 = IAU(j)+ZERO_ORIGIN
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



          DO i = Temp_R_IN(col - 1) + 1, Temp_R_IN(col) 
             row = R % CN(i)
             v   = R % V(i)
             IF(row <= N) THEN
                k = Temp_sz(row)
                IF(k == 0 .OR. Temp_int(k, row) /= col) THEN
                   k = k + 1
                   Temp_sz (row)  = k
                   Temp_int(k, row) = col
                   Temp    (k, row) = v * D(row)
                ELSE 
                   Temp(k, row) = Temp(k, row) + v * D(row)
                END IF
             END IF
                             

             DO j = INL(row - 1) + 1, INL(row)
                row2 = IAL(j)+ZERO_ORIGIN
                !C --
                IF( row2 > N) THEN
                   CYCLE
                END IF
                !C --

                k = Temp_sz(row2)
                
                IF(k == 0 .OR. Temp_int(k, row2) /= col) THEN
                   k = k + 1
                   Temp_sz(row2) = k
                   Temp_int(k, row2) = col
                   Temp    (k, row2) = v * AL(j)
                ELSE
                   Temp(k, row2) = Temp(k, row2) + v * AL(j)
                END IF
             END DO

             DO j = INU(row - 1) + 1, INU(row)
                row2 = IAU(j) + ZERO_ORIGIN
                !C --
                IF(row2 > N) CYCLE
                !C --
                   
                k = Temp_sz(row2)
                IF(k == 0 .OR. Temp_int(k, row2) /= col) THEN
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
              
          js = Temp_R_IN(col - 1) + 1
          je = Temp_R_IN(col)
          DO j = js, je
             k = R % CN(j)
             v = R % V(j)
             !C--
             IF(k > N) CYCLE
             !C--
             DO l = 1,Temp_sz(k)
                row = Temp_int(l,k)
                IF(col >= row .OR. row > ownaggre_size) THEN
                   
                   DO m=1, nd_sz
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
             
             CALL hash_zero(hash_size, Temp_CN(1:hash_size,row), col, count, "NONZEROS_PER_ROW")
             
             !C assign the value
             Temp_N(row) = Temp_N(row) + 1
             Temp_CN(count, row) = col
             Temp_V(count, row)  = Temp_1_col(row)

             Temp_1_col(row) = 0
          END DO
          !C-- end: input
          
       END DO
       !C-- end:temp * P

    END DO
    
    DEALLOCATE( Temp_nd    )
    DEALLOCATE( Temp_sz    )
    DEALLOCATE( Temp_int   )
    DEALLOCATE( Temp       )
    DEALLOCATE( Temp_agr   )
    DEALLOCATE( Temp_1_col )
  END SUBROUTINE RAP2


  SUBROUTINE  DEBUG_BARRIER_PRINT(string, comm, my_rank)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
#include "precision.inc"    

    CHARACTER(len=*), intent(in) :: string
    INTEGER(kind=kint), intent(in) :: comm, my_rank
    INTEGER(kind=kint) :: ierr

#ifdef DEBUG    
    CALL MPI_BARRIER(comm, ierr)
    IF (my_rank == 0) THEN
       WRITE(*,*) string
    END IF
    CALL MPI_Barrier(comm, ierr)
#endif    

  END SUBROUTINE DEBUG_BARRIER_PRINT
  
  
!==========================================================================
  SUBROUTINE make_rap(LEVEL_NO, my_rank, SOLVER_COMM, NPROCS,                      &
       & coarser_level_size, global_finish_flag, GIN_aggregate,                    &
       & PE_list_size, PE_list, aggregate_table_size, aggregate_table_array,       &
       & aggregate_number_in_table, HASH_TABLE_SIZE, global_local_hash_table,      &
       & local_aggre_size, Temp_R_IN, TEMP_COLSIZE, TEMP_RAP_WORKSIZE)

    USE data_structure_for_AMG
    USE solver_SR2
    USE solver_Gnumbering
    USE count_time_mod
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER(kind=kint), INTENT(in) :: LEVEL_NO
    LOGICAL,            INTENT(in) :: global_finish_flag
    INTEGER(kind=kint), INTENT(in) :: coarser_level_size
    INTEGER(kind=kint), INTENT(in) :: TEMP_COLSIZE, TEMP_RAP_WORKSIZE
    
    TYPE(INTER_LEVEL_OPERATOR), POINTER :: R, P

    !C temporary space of RAP
    REAL   (kind=kreal), ALLOCATABLE :: Temp_V(:, :)
    INTEGER(kind=kint ), ALLOCATABLE :: Temp_CN(:, :)
    INTEGER(kind=kint ), ALLOCATABLE :: Temp_N(:)

    INTEGER(kind=kint),INTENT(inout) :: local_aggre_size
    INTEGER(kind=kint ) :: Temp_R_IN(0:local_aggre_size)    
    
    !C pointers for matrix
    REAL   (kind=kreal), POINTER ::  newD(:), newAU(:), newAL(:)
    INTEGER(kind=kint ), POINTER ::  newINU(:), newIAU(:), newINL(:), newIAL(:)

    !C communication data structures
    INTEGER(kind=kint ), INTENT(in) :: NPROCS, SOLVER_COMM
    INTEGER(kind=kint ), INTENT(in) :: my_rank

    INTEGER(kind=kint )          :: NEIBPETOT
    INTEGER(kind=kint ), POINTER :: NEIBPE    (:)
    INTEGER(kind=kint ), POINTER :: NOD_IMPORT(:), STACK_IMPORT(:)
    INTEGER(kind=kint ), POINTER :: NOD_EXPORT(:), STACK_EXPORT(:)
    INTEGER(kind=kint ) :: ierr
    
    !C GIN_aggregate is used as rowdistbx on toplevel.
    INTEGER(kind=kint ), INTENT(inout) :: GIN_aggregate(0:NPROCS)

    INTEGER(kind=kint) :: i, j, k, l, m, g, N, NP, row, colin, col, n_lower, n_upper, count, index
    INTEGER(kind=kint) :: rin, trin, l_count, u_count

    REAL   (kind=kreal), ALLOCATABLE :: WS(:), WR(:)


    INTEGER(kind=kint ), INTENT(inout) :: PE_list(1:NPROCS,1:3)
    INTEGER(kind=kint ), INTENT(inout) :: PE_list_size 

    !C - begin: hash_table
    INTEGER(kind=kint), INTENT(in)    :: HASH_TABLE_SIZE
    INTEGER(kind=kint) :: global_local_hash_table(1:HASH_TABLE_SIZE, 1:2)
    !C - end: hash_table
    
    INTEGER(kind=kint), INTENT(in)    :: aggregate_table_size
    !C- aggregate_table_array(aggregate_table_size*2,max_neib_pe_size)
    !C- aggregate_number_in_table(max_neib_pe_size)
!!$    INTEGER(kind=kint) :: aggregate_table_array(1:aggregate_table_size * 2, 1:NPROCS)
!!$    INTEGER(kind=kint) :: aggregate_number_in_table(1:NPROCS)
    INTEGER(kind=kint) :: aggregate_table_array(:, :)
    INTEGER(kind=kint) :: aggregate_number_in_table(:)
    
    
    INTEGER(kind=kint) :: newNEIBPETOT, newNP, newN
    INTEGER(kind=kint), POINTER :: newNOD_IMPORT(:), newNOD_EXPORT(:)
    INTEGER(kind=kint), POINTER :: newSTACK_IMPORT(:), newSTACK_EXPORT(:)
    INTEGER(kind=kint), POINTER :: newNEIBPE(:)
    
    INTEGER(kind=kint) :: gl_crsest_sz,gbase, displs(0:NPROCS), recvcounts(0:NPROCS)
    INTEGER(kind=kint),ALLOCATABLE :: local_global_vector(:)
    REAL   (kind=kreal),ALLOCATABLE   :: send_buf(:)
    
    INTEGER(kind=kint),ALLOCATABLE :: export_global_buffer(:)
    INTEGER(kind=kint) :: idx_PE_INT,idx_NOD_INT,idx_PE,idx_NOD
    INTEGER(kind=kint) :: newNEIBPETOT_INT
    INTEGER(kind=kint), POINTER :: newNOD_IMPORT_INT(:), newNOD_EXPORT_INT(:)
    INTEGER(kind=kint), POINTER :: newSTACK_IMPORT_INT(:), newSTACK_EXPORT_INT(:)
    INTEGER(kind=kint), POINTER :: newNEIBPE_INT(:)
    
    INTEGER(kind=kint), ALLOCATABLE :: export_int_lvl_buffer(:)
    INTEGER(kind=kint) :: m1, m2
    INTEGER(kind=kint), ALLOCATABLE :: comm_recv_node(:)

#ifdef CHECK_make_export_table
    CHARACTER(len=6) :: fn
    INTEGER(kind=kint), ALLOCATABLE :: Temp(:)
#endif

    N  = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP 

    NEIBPETOT    =  HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NEIBPETOT         
    if(NEIBPETOT > 0) then
       NEIBPE       => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NEIBPE
       STACK_IMPORT => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % STACK_IMPORT
       NOD_IMPORT   => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NOD_IMPORT
       STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % STACK_EXPORT
       NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NOD_EXPORT
    end if


    CALL DEBUG_BARRIER_PRINT("### allocation and RAP  ###", SOLVER_COMM, my_rank)

    ALLOCATE(Temp_N(local_aggre_size))
    ALLOCATE(Temp_V(TEMP_COLSIZE, local_aggre_size))
    ALLOCATE(Temp_CN(TEMP_COLSIZE, local_aggre_size))

    Temp_N (1:local_aggre_size) = 0
    Temp_CN(1:TEMP_COLSIZE, 1:local_aggre_size) = 0
    CALL count_time(1, solver_comm, my_rank, 4)       


    if(SYMMETRIC_FLAG) then
       call RAP(Temp_N(1:local_aggre_size), Temp_CN(1:TEMP_COLSIZE,1:local_aggre_size),&
            & Temp_V(1:TEMP_COLSIZE,1:local_aggre_size), TEMP_COLSIZE, N, &
            & coarser_level_size, LEVEL_NO, local_aggre_size, my_rank)

!!$       call RAP(Temp_N, Temp_CN,&
!!$            & Temp_V, TEMP_COLSIZE, N, &
!!$            & coarser_level_size, LEVEL_NO, local_aggre_size, my_rank)
       
!!$       CALL RAP2(Temp_N, Temp_CN, Temp_V, TEMP_COLSIZE, N, coarser_level_size, &
!!$            & LEVEL_NO, local_aggre_size, Temp_R_IN, TEMP_RAP_WORKSIZE)


    else
       IF (N>0) THEN
       call RAP_unsym(Temp_N(1:local_aggre_size), Temp_CN(1:TEMP_COLSIZE,1:local_aggre_size), &
            & Temp_V(1:TEMP_COLSIZE,1:local_aggre_size), TEMP_COLSIZE, N, &
            & coarser_level_size, LEVEL_NO, local_aggre_size)
       END IF
    end if


    
    CALL count_time(2, solver_comm, my_rank, 4)       
    !C-- end:allocate temporaly domain and RAP().



    CALL DEBUG_BARRIER_PRINT("### sr_elements         ###", SOLVER_COMM, my_rank)

    !C-- send the not-own aggregates and receive own aggregates in matrix(Temp_CN,_N,_V)
    CALL count_time(1, solver_comm, my_rank, 3)       
    IF(NPROCS > 1) THEN
       i = local_aggre_size
       CALL solver_SR_elements &
            & (SOLVER_COMM, my_rank, GIN_aggregate(0:NPROCS), NPROCS, HASH_TABLE_SIZE,         &
            & global_local_hash_table(1:HASH_TABLE_SIZE, 1:2), aggregate_table_size, &
!!$            & aggregate_table_array(1:aggregate_table_size*2, 1:max_neib_pe_size), &
            & aggregate_table_array(:, :), &
!!$            & aggregate_number_in_table(1:max_neib_pe_size), local_aggre_size, PE_list(1:NPROCS,1:3), &
            & aggregate_number_in_table(:), local_aggre_size, PE_list(1:NPROCS,1:3), &
            & PE_list_size, Temp_CN, Temp_N, &
            & Temp_V, TEMP_COLSIZE, i)
    END IF



    CALL count_time(2,solver_comm,my_rank,3)       
    !C-- end:send the not-own aggregates and receive own aggregates in matrix(Temp_CN,_N,_V)
#ifdef ALLOCATION_CHECK
    DO i = 1, coarser_level_size
       IF(Temp_N(i) > TEMP_COLSIZE) THEN
          STOP 'Temp_col_size should be enlarged: sr_elements'
       END IF
    END DO
#endif

    !C-- allocate the location and copy Temp matrix to correct location
    IF(global_finish_flag) THEN
#ifdef CHOLESCKY       
       !C== LU factoryzation and assign to AL
       !C symmetricity is assumed.
       !C INL has node size for each PEs (equal to GIN_aggregate)
       !C NPL is the size of triangular matrix's elements

       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % INL(0:NPROCS))
       newINL => HIERARCHICAL_DATA(LEVEL_NO) % INL
       newINL =  GIN_aggregate
       gl_crsest_sz = newINL(NPROCS)
       HIERARCHICAL_DATA(LEVEL_NO) % NPL = (1 + gl_crsest_sz) * gl_crsest_sz / 2
       
       !C local_aggre_size is meant to be NP on coarser level

       !C--making local_global_vector
       ALLOCATE(local_global_vector(local_aggre_size))
       local_global_vector = 0
       gbase = GIN_aggregate(my_rank)
       DO i = 1, coarser_level_size
          local_global_vector(i) = gbase + i
       END DO
       DO i = 1, PE_list_size
          DO j = 1, aggregate_number_in_table(i)
             !C local
             k = aggregate_table_array(j * 2, i)
             !C global
             local_global_vector(k) = aggregate_table_array(j * 2 - 1, i)
          END DO
       END DO
       !C--end: making local_global_vector       

       !C-- displs and recvcounts,send_buf  are created
       !C send_buf has data in lower triangular format.
       displs(0) = 0
       DO i = 1, NPROCS
          j = newINL(i)
          displs(i) = (gl_crsest_sz * 2 - j + 1) * j / 2
          recvcounts(i - 1) = displs(i) - displs(i - 1)
       END DO
       ALLOCATE(send_buf(recvcounts(my_rank) + 1))
       send_buf = 0

       l = gl_crsest_sz - gbase
       DO row = 1, coarser_level_size
          k = l * (row - 1) - row * (row - 1) / 2
          DO colin = 1, TEMP_COLSIZE
             col = Temp_CN(colin, row)
             IF(col > 0) THEN
                i = local_global_vector(col) - gbase
                IF(row <= i) THEN
                   send_buf(k + i) = Temp_V(colin, row)
                END IF
             END IF
          END DO
       END DO
       !C--end: displs and recvcounts,send_buf  are created
       
       !C-- exchange the data of matrix
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % AL( HIERARCHICAL_DATA(LEVEL_NO) % NPL ))
       newAL => HIERARCHICAL_DATA(LEVEL_NO) % AL
       newAL = 0
       CALL MPI_ALLGATHERV(send_buf, recvcounts(my_rank), MPI_DOUBLE_PRECISION, &
            & newAL, recvcounts(0), displs(0), MPI_DOUBLE_PRECISION, SOLVER_COMM, ierr)

       DEALLOCATE(send_buf)
       !C--end: exchange the data of matrix       

       CALL left_looking_chorescky(newAL, gl_crsest_sz)
       
       !C-- X is array of local NP size vector. B is global vector on coarsest level.
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % X(local_aggre_size))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % B(gl_crsest_sz))
       !C--
       
       !C== end: LU factoryzation and assign to AL
#elif defined(PSPASES)
       CALL mpi_barrier(solver_comm,ierr)
       IF(my_rank==0) WRITE(*,*) "in"
       CALL make_PSPSE_MAT(Temp_CN,Temp_N,Temp_V,TEMP_COLSIZE,GIN_aggregate,&
            & aggregate_table_array,aggregate_number_in_table,PE_list_size,PE_list,&
            & local_aggre_size,coarser_level_size,GIN_aggregate(my_rank),SOLVER_COMM)
       
       CALL mpi_barrier(solver_comm,ierr)
       IF(my_rank==0) WRITE(*,*) "out"
       
       ALLOCATE(rowdistbx(0:NPROCS))
       rowdistbx = GIN_aggregate
       HIERARCHICAL_DATA(LEVEL_NO) % X=>coarsest_X(:,1)
       HIERARCHICAL_DATA(LEVEL_NO) % B=>coarsest_B(:,1)
#endif

    ELSE
       
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % INU(0:local_aggre_size))
       ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % INL(0:local_aggre_size))
       newINL => HIERARCHICAL_DATA(LEVEL_NO) % INL
       newINU => HIERARCHICAL_DATA(LEVEL_NO) % INU

       newINL(0:local_aggre_size) = 0
       newINU(0:local_aggre_size) = 0


       if(SYMMETRIC_FLAG) then

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
             newINU(i) = newINU(i - 1) + newINU(i)
             newINL(i) = newINL(i - 1) + newINL(i)
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

          !C this array stores recv nodes.
          allocate(comm_recv_node(local_aggre_size - coarser_level_size))
          comm_recv_node = 0

          DO row = 1, coarser_level_size
             n_upper = newINU(row - 1)
             DO colin = 1, TEMP_COLSIZE
                col = Temp_CN(colin, row)
                IF(col > 0) THEN
                   IF(row == col) THEN 
                      newD(row) = Temp_V(colin, row)

                   ELSE IF(row < col .AND. EPS < abs(Temp_V(colin,row))) THEN
                      !C upper part
                      n_upper = n_upper + 1
                      newAU(n_upper) = Temp_V(colin, row)
                      newIAU(n_upper) = col - ZERO_ORIGIN
                      !C record for recv table
                      if(col > coarser_level_size) then
                         comm_recv_node(col - coarser_level_size) = 1
                      end if
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

       else
          !C unsymmetric case
          DO row = 1, coarser_level_size
             DO colin = 1, TEMP_COLSIZE
                col = Temp_CN(colin, row)
                IF(col > 0) then
                   IF(col > row .AND. EPS < abs(Temp_V(colin, row))) THEN
                      newINU(row) = newINU(row) + 1
                   ELSE IF(col < row .AND. EPS < abs(Temp_V(colin, row))) THEN
                      newINL(row) = newINL(row) + 1
                   END IF
                end IF
             END DO
          END DO



          DO i = 1, local_aggre_size
             newINU(i) = newINU(i - 1) + newINU(i)
             newINL(i) = newINL(i - 1) + newINL(i)
          END DO
          u_count = newINU(local_aggre_size)
          l_count = newINL(local_aggre_size)

          HIERARCHICAL_DATA(LEVEL_NO) % NPU = u_count
          HIERARCHICAL_DATA(LEVEL_NO) % NPL = l_count

          ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % IAL(l_count))
          ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % IAU(u_count))
          ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % AL (l_count))
          ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % AU (u_count))
          ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % D(local_aggre_size))
          ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % X(local_aggre_size))
          ALLOCATE(HIERARCHICAL_DATA(LEVEL_NO) % B(local_aggre_size))

          newD  => HIERARCHICAL_DATA(LEVEL_NO) % D
          newIAL=> HIERARCHICAL_DATA(LEVEL_NO) % IAL
          newIAU=> HIERARCHICAL_DATA(LEVEL_NO) % IAU
          newAL => HIERARCHICAL_DATA(LEVEL_NO) % AL
          newAU => HIERARCHICAL_DATA(LEVEL_NO) % AU

          !C this array stores recv nodes.
          allocate(comm_recv_node(local_aggre_size - coarser_level_size))
          comm_recv_node = 0
          
          DO row = 1, coarser_level_size
             n_upper = newINU(row - 1)
             n_lower = newINL(row - 1)
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
                      !C record for recv table
                      if(col > coarser_level_size) then
                         comm_recv_node(col - coarser_level_size) = 1
                      end if
                   ELSE IF(row > col .AND. EPS < abs(Temp_V(colin, row))) THEN
                      !C lower part
                      n_lower = n_lower + 1
                      newAL(n_lower) = Temp_V(colin, row)
                      newIAL(n_lower) = col - ZERO_ORIGIN
                   END IF
                END IF
             END DO
          END DO
       end IF
    END IF

    DEALLOCATE(Temp_N)
    DEALLOCATE(Temp_V)
    DEALLOCATE(Temp_CN)
    !C--end: allocate the location and copy Temp matrix to correct location

    !C-- remake R
    R => HIERARCHICAL_DATA(LEVEL_NO) % R
    P => HIERARCHICAL_DATA(LEVEL_NO) % P
    ALLOCATE(R % IN(0:local_aggre_size))
    if(local_aggre_size > 0) then
       R % IN = 0
    end if
    DO i = 1, N
       DO j = P % IN(i - 1) + 1, P % IN(i)
          col = P % CN(j)
          R % IN(col) = R % IN(col) + 1
       END DO
    END DO

    DO i = 2, local_aggre_size
       R % IN(i) = R % IN(i) + R % IN(i - 1)
    END DO

    DO i = 2, R % ROW_SIZE
       rin = R % IN(i - 1)
       trin = Temp_R_IN(i - 1)
       j = trin - rin
       k = R % IN(i) - rin

       DO l = trin + 1, trin + k
          R % V (l - j) = R %  V(l)
          R % CN(l - j) = R % CN(l)
       END DO
    END DO

    R % ROW_SIZE = local_aggre_size    
    !C-- end:remake R

    CALL DEBUG_BARRIER_PRINT("### make communication tables ###", SOLVER_COMM, my_rank)
    !C-- communication table creation

    j = 0
    DO i = 1, PE_list_size
       j = j + aggregate_number_in_table(i)
    END DO


    newN = coarser_level_size
    newNP = local_aggre_size
    HIERARCHICAL_DATA(LEVEL_NO) % N = newN
    HIERARCHICAL_DATA(LEVEL_NO) % NP = newNP

    ALLOCATE(export_global_buffer(j))    
    allocate(newNOD_IMPORT(j))
    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_IMPORT => newNOD_IMPORT
    allocate(newSTACK_IMPORT(0:NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT => newSTACK_IMPORT
    allocate(newSTACK_EXPORT(0:NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT => newSTACK_EXPORT
    allocate(newNEIBPE(NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPE => newNEIBPE



    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(0) = 0
    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT(0) = 0 


    ALLOCATE(export_int_lvl_buffer(j))    
    allocate(newNOD_IMPORT_INT(j))
    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NOD_IMPORT => newNOD_IMPORT_INT
    allocate(newSTACK_IMPORT_INT(0:NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_IMPORT => newSTACK_IMPORT_INT
    allocate(newSTACK_EXPORT_INT(0:NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_EXPORT => newSTACK_EXPORT_INT
    allocate(newNEIBPE_INT(NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPE => newNEIBPE_INT

    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_EXPORT(0) = 0
    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_IMPORT(0) = 0

#if defined(CHOLESCKY)||defined(PSPASES)
    !C- only PSPASES needs no communication tables on the coarsest level
    IF(global_finish_flag) THEN
       !C-- check for comm_table and select important nodes
       idx_NOD = 0
       idx_PE  = 0
       idx_NOD_INT = 0
       idx_PE_INT  = 0
       DO i = 1, PE_list_size
          k = aggregate_number_in_table(i)
          DO j = 1, k
             l = aggregate_table_array(j * 2, i)
             g = aggregate_table_array(j * 2 - 1, i)
             m2 = R % IN(l) - R % IN(l - 1)

             !C we ignore coefficients among external nodes,so m1 /=0 means 
             !C l is connected to own node.
             IF(m2 > 0) THEN
                idx_NOD_INT = idx_NOD_INT + 1
                export_int_lvl_buffer(idx_NOD_INT) = g
                newNOD_IMPORT_INT(idx_NOD_INT)     = l-ZERO_ORIGIN
             END IF
          END DO
          IF(idx_NOD_INT > newSTACK_IMPORT_INT(idx_PE_INT)) THEN
             idx_PE_INT = idx_PE_INT + 1
             newNEIBPE_INT(idx_PE_INT) = PE_list(i, 1)
             newSTACK_IMPORT_INT(idx_PE_INT) = idx_NOD_INT
          END IF
       END DO


       newNEIBPETOT_INT= idx_PE_INT
       newNEIBPETOT    = idx_PE
       HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPETOT = newNEIBPETOT_INT
       !C--end: check for comm_table and select important nodes

       CALL make_export_NOD &
            & (SOLVER_COMM, export_int_lvl_buffer, GIN_aggregate, my_rank, NPROCS, &
            &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPETOT,            &
            &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPE,               &
            &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NOD_EXPORT,           &
            &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_IMPORT,         &
            &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_EXPORT)

       newNEIBPETOT  =  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPETOT 


       DEALLOCATE(export_int_lvl_buffer)
       DEALLOCATE(export_global_buffer)
       !C-- end:communication table creation

    ELSE
#endif

    !C-- check for comm_table and select important nodes
    idx_NOD = 0
    idx_PE = 0
    idx_NOD_INT = 0
    idx_PE_INT = 0
    DO i = 1, PE_list_size
       k = aggregate_number_in_table(i)
       DO j = 1, k
          l = aggregate_table_array(j * 2, i)
          g = aggregate_table_array(j * 2 - 1, i)

          m1 = comm_recv_node(l - coarser_level_size)

          m2 = R % IN(l) - R % IN(l - 1)

          
          !C we ignore coefficients among external nodes, so m1 /=0 means 
          !C l is connected to own node.
          IF(m1 > 0) THEN
             idx_NOD = idx_NOD + 1

             export_global_buffer(idx_NOD) = g

             HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_IMPORT(idx_NOD) = l-ZERO_ORIGIN

          END IF

          IF(m2 > 0) THEN
             idx_NOD_INT = idx_NOD_INT + 1

             export_int_lvl_buffer(idx_NOD_INT) = g

             HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NOD_IMPORT(idx_NOD_INT) = l-ZERO_ORIGIN

          END IF
       END DO
       IF(idx_NOD > HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(idx_PE)) THEN

          idx_PE = idx_PE + 1
          HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPE(idx_PE) = PE_list(i, 1)
          HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(idx_PE) = idx_NOD       
       END IF
       IF(idx_NOD_INT > HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_IMPORT(idx_PE_INT)) THEN

          idx_PE_INT = idx_PE_INT + 1
          HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPE(idx_PE_INT) = PE_list(i, 1)
          HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_IMPORT(idx_PE_INT) = idx_NOD_INT
       END IF
    End DO

    newNEIBPETOT_INT = idx_PE_INT
    newNEIBPETOT     = idx_PE

    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE    % NEIBPETOT = newNEIBPETOT
    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPETOT = newNEIBPETOT_INT
    !C--end: check for comm_table and select important nodes

    CALL make_export_NOD &
         & (SOLVER_COMM, export_int_lvl_buffer,GIN_aggregate,my_rank,NPROCS,     &
         &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPETOT,             &
         &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NEIBPE(1:NPROCS),      &
         &  newNOD_EXPORT_INT,            &
         &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_IMPORT(0:NPROCS),&
         &  HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % STACK_EXPORT(0:NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % INT_LVL_TABLE % NOD_EXPORT => newNOD_EXPORT_INT

    CALL mpi_barrier(solver_comm,ierr)

    CALL make_export_NOD &
         & (SOLVER_COMM, export_global_buffer,GIN_aggregate,my_rank,NPROCS,   &
         &  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPETOT,             &
         &  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPE(1:NPROCS),      &
         &  newNOD_EXPORT,            &
         &  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(0:NPROCS),&
         &  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT(0:NPROCS))
    HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_EXPORT => newNOD_EXPORT

    newNEIBPETOT  =  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPETOT 
    m1 = HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT(newNEIBPETOT)


    DEALLOCATE(comm_recv_node)
    DEALLOCATE(export_int_lvl_buffer)
    DEALLOCATE(export_global_buffer)
    !C-- end:communication table creation

#if defined(CHOLESCKY)||defined(PSPASES)
    endif
#endif


    !C diagonal part is to have correct values in overlapped region.
    IF(.NOT. global_finish_flag) THEN


       if(NPROCS>1 .and. newNEIBPETOT>0) then
          ALLOCATE ( WS(HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT(newNEIBPETOT)) )
          ALLOCATE ( WR(HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(newNEIBPETOT)) )

          CALL DEBUG_BARRIER_PRINT("newD change", SOLVER_COMM, my_rank)
          m1 =  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT(newNEIBPETOT)
          m2 =  HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(newNEIBPETOT)
          CALL SOLVER_SEND_RECV2                                                           &
               &   ( newNP,&
               & HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT(newNEIBPETOT),  &
               & HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(newNEIBPETOT),  &
               & newNEIBPETOT, &
               & HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NEIBPE(1:newNEIBPETOT),       &
               & HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_IMPORT(0:newNEIBPETOT), &
               & HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_IMPORT(1:m2),   &
               & HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % STACK_EXPORT(0:newNEIBPETOT), &
               & HIERARCHICAL_DATA(LEVEL_NO) % COMM_TABLE % NOD_EXPORT(1:m1),   &
               & WS, WR, newD(1:newNP), SOLVER_COMM, my_rank)

          DEALLOCATE(WS)
          DEALLOCATE(WR)
       end if
    END IF
    

  CONTAINS

    SUBROUTINE solver_SR_elements &
         &    (SOLVER_COMM, my_rank, GIN_aggregate, NPROCS, HASH_TABLE_SIZE, &
         &     global_local_hash_table, aggregate_table_size,                &
         &     aggregate_table_array, aggregate_number_in_table,             &
         &     local_aggre_size, PE_list, PE_list_size, Temp_CN, Temp_N,     &
         &     Temp_V, TEMP_COLSIZE, temp_2_size) 
      USE hash_mod
      USE data_structure_for_AMG
      IMPLICIT NONE
      INCLUDE  'mpif.h'

      INTEGER(kind=kint ), INTENT(in) :: NPROCS
      !C--  DIMENSION(0:NPROCS)
      INTEGER(kind=kint ), INTENT(in) :: GIN_aggregate(0:NPROCS)
      INTEGER(kind=kint ), INTENT(inout) :: PE_list_size
      !C-- DIMENSION(NPROCS,3)
      INTEGER(kind=kint ), INTENT(inout) :: PE_list(1:NPROCS, 1:3)

      INTEGER(kind=kint ), INTENT(in) :: SOLVER_COMM
      INTEGER(kind=kint ), INTENT(in) :: my_rank
      INTEGER(kind=kint ), ALLOCATABLE :: req1(:), req2(:)
      INTEGER(kind=kint ), ALLOCATABLE :: sta1(:,:), sta2(:,:)
      INTEGER(kind=kint ) :: ierr

      INTEGER(kind=kint ), INTENT(inout) :: local_aggre_size

      INTEGER(kind=kint),INTENT(in):: HASH_TABLE_SIZE
      !C DIMENSION(HASH_TABLE_SIZE, 2)
      INTEGER(kind=kint) :: global_local_hash_table(1:HASH_TABLE_SIZE, 1:2)

      INTEGER(kind=kint), INTENT(in) :: aggregate_table_size
      !C two arguments: communication table,each PE,(aggregate_table_size*2,max_neib_pe_size)
!!$      INTEGER(kind=kint) :: aggregate_table_array(1:aggregate_table_size*2, 1:NPROCS)
      INTEGER(kind=kint) :: aggregate_table_array(:, :)
      !C argument: the number of aggregates,(NPROCS)
!!$      INTEGER(kind=kint) :: aggregate_number_in_table(1:NPROCS)
      INTEGER(kind=kint) :: aggregate_number_in_table(:)

      INTEGER(kind=kint), INTENT(in)    :: TEMP_COLSIZE, temp_2_size
      INTEGER(kind=kint), INTENT(inout) :: Temp_N(1:temp_2_size)
      REAL(kind=kreal),   INTENT(inout) :: Temp_V(1:TEMP_COLSIZE, 1:temp_2_size)
      INTEGER(kind=kint), INTENT(inout) :: Temp_CN(1:TEMP_COLSIZE, 1:temp_2_size)

      INTEGER(kind=kint) :: rs_in_buf_sz, sd_buf_sz, sd_in_buf_sz
      INTEGER(kind=kint) :: sd_in_buf_index, sd_buf_index, rs_buf_sz, index_record
      INTEGER(kind=kint), ALLOCATABLE :: rs_pes(:)
      INTEGER(kind=kint), ALLOCATABLE :: in_sd_in_buf(:), in_sd_buf(:)
      INTEGER(kind=kint), ALLOCATABLE :: sd_in_buf(:), sd_row_glo_buf(:)
      INTEGER(kind=kint), ALLOCATABLE :: sd_col_glo_buf(:), rs_col_glo_buf(:)
      INTEGER(kind=kint), ALLOCATABLE :: rs_in_buf(:), rs_row_glo_buf(:)
      REAL   (kind=kreal),ALLOCATABLE :: sd_values_buf(:), rs_values_buf(:)
      
      INTEGER(kind=kint),ALLOCATABLE :: in_recv_buf(:)
      
      INTEGER(kind=kint) :: dest, istart, inum
      INTEGER(kind=kint), ALLOCATABLE :: LG_vector(:)

      INTEGER(kind=kint) :: i, j, k, l, g, ownnodes, buf_index, buf_size, message_size, aggregate_number
      INTEGER(kind=kint) :: in_buf_size, ins, ine, n, global_base
      INTEGER(kind=kint) :: row, col, global_bound, hash_no, hash_pointer, GLOBAL, LOCAL, m, rank
      INTEGER(kind=kint) :: count, n_pes

      INTEGER(kind=kint), ALLOCATABLE  :: pes_s(:), pes_r(:)


      INTEGER(kind=kint), ALLOCATABLE :: in_buf(:), row_glo_buf(:), col_glo_buf(:)
      PARAMETER(GLOBAL = 1, LOCAL = 2)

      ALLOCATE(pes_s(3 * NPROCS), pes_r(3 * NPROCS))




      !C-- LG_vector creation
      ALLOCATE(LG_vector(local_aggre_size))
      global_base = GIN_aggregate(my_rank)
      global_bound = GIN_aggregate(my_rank + 1)
      ownnodes = global_bound - global_base

#ifdef DEBUG
      LG_vector = 0
#endif

      !C create LG_vector
      DO i = 1, ownnodes
         LG_vector(i) = i + global_base
      END DO

      DO i = 1, PE_list_size
         DO j = 1, aggregate_number_in_table(i)
            k = 2 * j
            l = aggregate_table_array(k, i)
            LG_vector(l) = aggregate_table_array(k - 1, i)
         END DO
      END DO

#ifdef DEBUG
      DO i = 1, local_aggre_size
         IF(LG_vector(i) == 0) WRITE(*,*) "Node", i , " is only for inter level communication."
      END DO
#endif
      !C-- end: LG_vector creation


      !C-- send buffer creation
      sd_in_buf_sz = 0
      sd_buf_sz    = 0
      DO i = 1, PE_list_size
         in_buf_size = aggregate_number_in_table(i)
         sd_in_buf_sz = sd_in_buf_sz + in_buf_size
         DO j = 1, in_buf_size
            l = aggregate_table_array(j * 2, i)
            sd_buf_sz = sd_buf_sz + Temp_N(l)
         END DO
      END DO


      ALLOCATE(sd_in_buf     (sd_in_buf_sz + 1))
      ALLOCATE(sd_row_glo_buf(sd_in_buf_sz + 1))
      ALLOCATE(sd_col_glo_buf(sd_buf_sz + 1))
      ALLOCATE(sd_values_buf (sd_buf_sz + 1))

      ALLOCATE(in_sd_in_buf(0:PE_list_size))
      ALLOCATE(in_sd_buf   (0:PE_list_size))
      in_sd_in_buf(0) = 0
      in_sd_buf(0)    = 0

      pes_s = 0
      sd_in_buf_index = 0; sd_buf_index = 0
      DO i = 1, PE_list_size
         in_buf_size = aggregate_number_in_table(i)
         dest = PE_list(i, 1)

#ifdef DEBUG
         IF(dest < 0 .OR. dest >= NPROCS) STOP "unexpected error! in sr_elements"
#endif 

         !C create index in in_buf of col_glo_buf  and values_buf
         DO j = 1, in_buf_size
            l = aggregate_table_array(j * 2, i)
            k = j + sd_in_buf_index
            sd_in_buf(k)      = Temp_N(l)
            sd_row_glo_buf(k) = LG_vector(l)
         END DO
         sd_in_buf_index = sd_in_buf_index + in_buf_size

         buf_index = 0
         DO j = 1, in_buf_size
            l = aggregate_table_array(j * 2, i)
            DO k = 1, TEMP_COLSIZE
               IF(Temp_CN(k ,l) > 0) THEN
                  buf_index = buf_index + 1
                  m = buf_index + sd_buf_index
                  sd_col_glo_buf(m) = LG_vector(Temp_CN(k, l))
                  sd_values_buf(m)  = Temp_V(k, l)
               END IF
            END DO
         END DO
         sd_buf_index = sd_buf_index + buf_index

         pes_s((dest + 1) * 3)  = buf_index
         pes_s((dest + 1) * 3 - 1)= in_buf_size
         pes_s((dest + 1) * 3 - 2)= 1
         
         in_sd_in_buf(i) = sd_in_buf_index
         in_sd_buf(i)    = sd_buf_index
      END DO

      DEALLOCATE(LG_vector)
      !C-- end: send buffer creation


      !C-- send recv
      !C-- exchange PEno,in_buf_size,buf_size
      pes_r = 0

      CALL MPI_ALLTOALL(pes_s, 3, LIS_MPI_INTEGER, pes_r, 3, LIS_MPI_INTEGER, &
           &            SOLVER_COMM, ierr)
      !C-- end: exchange in_buf_size, buf_size

      
      !C--    input from alltoall data
      n_pes=0;rs_in_buf_sz=0;rs_buf_sz=0
      DO i=1,NPROCS
         n_pes        = n_pes        + pes_r(i * 3 - 2)
         rs_in_buf_sz = rs_in_buf_sz + pes_r(i * 3 - 1)
         rs_buf_sz    = rs_buf_sz    + pes_r(i * 3)
      END DO
      
      ALLOCATE(rs_in_buf     (rs_in_buf_sz + 1))
      ALLOCATE(rs_row_glo_buf(rs_in_buf_sz + 1))
      ALLOCATE(rs_col_glo_buf(rs_buf_sz + 1)   )
      ALLOCATE(rs_values_buf (rs_buf_sz + 1)   )
      
      ALLOCATE(rs_pes(n_pes))
      ALLOCATE(in_recv_buf(0:n_pes))
      
      j = 0
      in_recv_buf(0) = 0
      DO i = 1, NPROCS
         IF(pes_r(i * 3 - 2) == 1) THEN
            j = j + 1
            rs_pes(j) = i - 1
            in_recv_buf(j) = in_recv_buf(j - 1) + pes_r(i * 3 - 1)
         END IF
      END DO
      !C--  alltoall data  input
      
      ALLOCATE(req1(PE_list_size))
      ALLOCATE(sta1(MPI_STATUS_SIZE, PE_list_size))
      ALLOCATE(req2(n_pes))
      ALLOCATE(sta2(MPI_STATUS_SIZE, n_pes))

      
      !C-- exchange 1: each row size of Temp_CN _N _V
      DO i = 1, PE_list_size
         dest  = PE_list(i, 1)
         istart= in_sd_in_buf(i-1)
         inum  = in_sd_in_buf(i)-istart
         CALL MPI_ISEND(sd_in_buf(istart+1), inum, LIS_MPI_INTEGER, dest, 0, SOLVER_COMM, &
              &         req1(i), ierr)
      END DO
      
      DO i = 1, n_pes
         istart = in_recv_buf(i - 1)
         inum   = in_recv_buf(i) - istart
         dest   = rs_pes(i) 
         CALL MPI_IRECV(rs_in_buf(istart + 1), inum, LIS_MPI_INTEGER, dest, 0, SOLVER_COMM, &
              &         req2(i), ierr)
      END DO

      CALL MPI_WAITALL (n_pes,        req2, sta2, ierr)    
      CALL MPI_WAITALL (PE_list_size, req1, sta1, ierr)    
      !C-- end: exchange 1

      !C-- exchange 2: row number is exchanged
      DO i = 1, PE_list_size
         dest   = PE_list(i, 1)
         istart = in_sd_in_buf(i - 1)
         inum   = in_sd_in_buf(i) - istart
         CALL MPI_ISEND(sd_row_glo_buf(istart + 1), inum, LIS_MPI_INTEGER, dest, 0, SOLVER_COMM, &
              &         req1(i), ierr)
      END DO

      DO i = 1, n_pes
         istart = in_recv_buf(i - 1)
         inum   = in_recv_buf(i) - istart
         dest   = rs_pes(i)

         CALL MPI_IRECV(rs_row_glo_buf(istart + 1), inum, LIS_MPI_INTEGER, dest , 0, &
              &         SOLVER_COMM, req2(i), ierr)
      END DO

      CALL MPI_WAITALL (n_pes,        req2, sta2, ierr)
      CALL MPI_WAITALL (PE_list_size, req1, sta1, ierr)    
      !C-- end: exchange 2

      !C-- change in_recv_buf(:)
      in_recv_buf(0) = 0
      DO i = 1, n_pes
         j = rs_pes(i) + 1
         in_recv_buf(i) = in_recv_buf(i - 1) + pes_r(j * 3)
      END DO
      !C-- end:change in_recv_buf(:)


      !C-- exchange 3: column number is exchanged
      DO i = 1, PE_list_size
         dest  = PE_list(i, 1)
         istart= in_sd_buf(i - 1)
         inum  = in_sd_buf(i) - istart
         CALL MPI_ISEND(sd_col_glo_buf(istart + 1), inum, LIS_MPI_INTEGER, dest, 0, SOLVER_COMM, &
              &         req1(i), ierr)
      END DO

      DO i = 1, n_pes
         dest = rs_pes(i)
         istart = in_recv_buf(i - 1)
         inum   = in_recv_buf(i) - istart
         CALL MPI_IRECV(rs_col_glo_buf(istart+1), inum, LIS_MPI_INTEGER, dest, 0, SOLVER_COMM, &
              &         req2(i), ierr)
      END DO
      
      CALL MPI_WAITALL (n_pes,        req2, sta2, ierr)          
      CALL MPI_WAITALL (PE_list_size, req1, sta1, ierr)          
      !C-- end: exchange 3


      !C-- exchange 4: value of Temp_N _CN _V is exchanged
      DO i=1,PE_list_size
         dest   = PE_list(i,1)
         istart = in_sd_buf(i - 1)
         inum   = in_sd_buf(i) - istart
         CALL MPI_ISEND(sd_values_buf(istart + 1), inum, MPI_DOUBLE_PRECISION, dest, 0, SOLVER_COMM, &
              &         req1(i), ierr)
      END DO

      DO i = 1, n_pes
         istart = in_recv_buf(i - 1)
         inum   = in_recv_buf(i) - istart
         dest   = rs_pes(i)

         CALL MPI_IRECV(rs_values_buf(istart + 1), inum, MPI_DOUBLE_PRECISION, dest, 0, SOLVER_COMM, &
              &         req2(i), ierr)
      END DO

      CALL MPI_WAITALL (n_pes,        req2, sta2, ierr)
      CALL MPI_WAITALL (PE_list_size, req1, sta1, ierr)    
      !C-- end: exchange 4


      DEALLOCATE(in_sd_in_buf)
      DEALLOCATE(in_sd_buf)
      DEALLOCATE(in_recv_buf)
      DEALLOCATE(sd_in_buf)
      DEALLOCATE(sd_row_glo_buf)
      DEALLOCATE(sd_col_glo_buf)
      DEALLOCATE(sd_values_buf)
      DEALLOCATE(req1)
      DEALLOCATE(sta1)
      DEALLOCATE(req2)
      DEALLOCATE(sta2)
      !C-- end: send recv


      index_record = 0
      DO j = 1, rs_in_buf_sz
         ins = index_record + 1
         ine = index_record + rs_in_buf(j)
         index_record = ine

         row = rs_row_glo_buf(j)-global_base
#ifdef DEBUG
         IF(row>ownnodes) STOP "unexpected error !! in sr_elements"
#endif
         DO k = ins,ine

            !C input data into associated location with hash function
            g = rs_col_glo_buf(k)
            !C  g in own aggregates
            !C  g not in own aggregates
            !C    g not in hashtable
            !C    g in hashtable
            IF(global_base < g .AND. g <= global_bound) THEN
               !C  g in own aggregates
               col = g - global_base

#ifdef DEBUG                  
               IF(col < 0) STOP "not own-aggregate is send in sr_elements"
#endif                  
               
               if(SYMMETRIC_FLAG .and. (col < row)) CYCLE
               !C search (row col) in new matrix

               CALL hash(TEMP_COLSIZE, Temp_CN(1:TEMP_COLSIZE,row), col, count, "NONZEROS_PER_ROW")               

               !C assign the value
               IF(Temp_CN(count,row)==col) THEN
                  Temp_V(count,row)=Temp_V(count,row)+rs_values_buf(k)
               ELSE
                  Temp_N(row)=Temp_N(row)+1
                  Temp_CN(count,row)=col
                  Temp_V(count,row)=rs_values_buf(k)
               END IF

            ELSE
               !C  g in neigbor aggregates
               CALL hash(HASH_TABLE_SIZE, global_local_hash_table(1:HASH_TABLE_SIZE, GLOBAL), g, hash_pointer, "HASH_TABLE_SIZE")

               l = global_local_hash_table(hash_pointer, LOCAL)
               IF(l == 0) THEN
                  !C g not in hashtable
                  !C record the new node
                  !C into hash_table
                  local_aggre_size = local_aggre_size + 1
                  global_local_hash_table(hash_pointer, GLOBAL) = g
                  global_local_hash_table(hash_pointer, LOCAL)  = local_aggre_size
                  l = local_aggre_size

                  !C into aggregate_table_array
                  DO m = 1, PE_list_size
                     rank = PE_list(m, 1)
                     IF(GIN_aggregate(rank) < g .AND. g <= GIN_aggregate(rank + 1)) THEN
                        EXIT
                     END IF
                  END DO

                  IF(m > PE_list_size) THEN
                     !C add new entry in PE_list in case of relayed node
                     DO n = 1, NPROCS
                        IF(GIN_aggregate(n - 1) < g .AND. g <= GIN_aggregate(n)) THEN
                           EXIT
                        END IF
                     END DO

                     IF(n > NPROCS) THEN
                        WRITE(*,*) "unknown node:solver_sr_elements: ",g
                        STOP 
                     END IF

                     PE_list_size = PE_list_size + 1
                     m = PE_list_size
                     PE_list(m,1) = n - 1
                     PE_list(m,2) = 0
                     PE_list(m,3) = 0
                     aggregate_number_in_table(m) = 0
                  END IF
                  
                  aggregate_number = aggregate_number_in_table(m) + 1
                  aggregate_number_in_table(m) = aggregate_number

#ifdef ALLOCATION_CHECK
                  IF(aggregate_table_size < aggregate_number) THEN
                     WRITE(*,*) "enlarge aggregate_table_size"
                     STOP 'error in solver_sr_element():'
                  END IF
#endif                  

                  aggregate_table_array(aggregate_number * 2 - 1, m) = g
                  aggregate_table_array(aggregate_number * 2, m) = l

                  col = l
                  !C search (row col) in new matrix
                  CALL hash(TEMP_COLSIZE, Temp_CN(1:TEMP_COLSIZE, row), col, count, "NONZEROS_PER_ROW")

                  !C assign the value
                  IF(Temp_CN(count, row) == col) THEN
                     STOP "solver_sr_element(): unsupported error!"
                  ELSE
                     Temp_N(row) = Temp_N(row) + 1
                     Temp_V(count,row) = rs_values_buf(k)
                     Temp_CN(count,row) = col
                  END IF

               ELSE
                  !C  g in hashtable
                  !C the case in which col is found in tables
                  col = l

                  !C search (row col) in new matrix
                  CALL hash(TEMP_COLSIZE, Temp_CN(1:TEMP_COLSIZE, row), col, count, "NONZEROS_PER_ROW")

                  !C assign the value
                  IF(Temp_CN(count, row) == col) THEN
                     Temp_V(count, row) = Temp_V(count, row) + rs_values_buf(k)
                  ELSE
                     Temp_N(row) = Temp_N(row) + 1
                     Temp_CN(count, row) = col
                     Temp_V(count, row) = rs_values_buf(k) 
                  END IF
               END IF
            END IF
         END DO
      END DO


      DEALLOCATE(rs_row_glo_buf)       
      DEALLOCATE(rs_in_buf)
      DEALLOCATE(rs_col_glo_buf)
      DEALLOCATE(rs_values_buf)
      DEALLOCATE(pes_s, pes_r)


#ifdef DEBUG
      DO i=1,PE_list_size
         IF(aggregate_number_in_table(i) > aggregate_table_size) THEN
            WRITE(*,*) "aggregate_table_array should be enlarged"
            STOP "aggregate_table_array:solver_SR_elements:solver_Gnum.f90:"
         END IF
      END DO
#endif

    END SUBROUTINE solver_SR_elements

    SUBROUTINE make_export_NOD (SOLVER_COMM, export_global_buffer,            &
         &     GIN_aggregate, my_rank, NPROCS, NEIBPETOT, NEIBPE, NOD_EXPORT, &
         &     STACK_IMPORT, STACK_EXPORT)

      IMPLICIT NONE
      INCLUDE  'mpif.h'
      INTEGER(kind=kint), INTENT(in)     :: NPROCS
      INTEGER(kind=kint), INTENT(inout)  :: NEIBPETOT
      INTEGER(kind=kint), INTENT(inout)  :: NEIBPE(NPROCS)
      INTEGER(kind=kint), INTENT(inout)  :: STACK_IMPORT(0:), STACK_EXPORT(0:)
      INTEGER(kind=kint), POINTER        :: NOD_EXPORT(:)
      
      INTEGER(kind=kint),  INTENT(in)    :: export_global_buffer(:)
      INTEGER(kind=kint), ALLOCATABLE    :: pes(:), in_recv_buf(:), recv_buf(:)
      
      INTEGER(kind=kint) :: pes_s(NPROCS), pes_r(NPROCS)

      INTEGER(kind=kint), INTENT(in)     :: GIN_aggregate(0:)

      INTEGER(kind=kint ) :: ierr
      INTEGER(kind=kint ), INTENT(in) :: SOLVER_COMM
      INTEGER(kind=kint ), INTENT(in) :: my_rank
      INTEGER(kind=kint), ALLOCATABLE :: req1(:), req2(:)
      INTEGER(kind=kint), ALLOCATABLE :: sta1(:,:), sta2(:,:)

      INTEGER(kind=kint) :: i, j, k, dest, inum, istart, recv_buf_size, message_size, l
      INTEGER(kind=kint) :: GLOBAL, LOCAL, base_index, n_pes, neib, from_pe, msg_sz
      PARAMETER(GLOBAL=1,LOCAL=2)
      LOGICAL :: hash_flag


      neib = NEIBPETOT


      !C-- determine number of recv
      pes_s = -1
      DO i = 1, neib
         dest = NEIBPE(i)
         pes_s(dest + 1)     = STACK_IMPORT(i) - STACK_IMPORT(i - 1)

#ifdef DEBUG
         IF(pes_s(dest+1)==0) STOP "0 size tables in make_export_NOD()"
#endif

      END DO
      pes_r = -1


      CALL MPI_ALLTOALL(pes_s, 1, LIS_MPI_INTEGER, pes_r, 1, LIS_MPI_INTEGER, &
           &            SOLVER_COMM, ierr)

      !C-- determine number of recv
      !C RECEIVE
      n_pes = 0
      ALLOCATE(in_recv_buf(0:NPROCS))
      ALLOCATE(pes(NPROCS))
      in_recv_buf(0) = 0
      DO i = 0, NPROCS - 1
         from_pe = i
         msg_sz  = pes_r(i + 1)
         IF(msg_sz >= 0) THEN
            n_pes              = n_pes + 1
            pes(n_pes)         = from_pe
            in_recv_buf(n_pes) = in_recv_buf(n_pes - 1) + msg_sz
         END IF
      END DO
      
      ALLOCATE(req1(neib))    
      ALLOCATE(sta1(MPI_STATUS_SIZE, neib))
      ALLOCATE(req2(n_pes))
      ALLOCATE(sta2(MPI_STATUS_SIZE, n_pes))    
      recv_buf_size = in_recv_buf(n_pes)
      ALLOCATE(recv_buf(recv_buf_size))
      ALLOCATE(NOD_EXPORT(recv_buf_size))
      DO i=1, n_pes
         istart = in_recv_buf(i - 1)
         inum   = in_recv_buf(i) - istart
         dest   = pes(i)
         CALL MPI_IRECV(recv_buf(istart + 1), inum, LIS_MPI_INTEGER, dest, 11, &
              &         SOLVER_COMM, req2(i), ierr)
      END DO

      !C SEND
      DO i = 1, neib
         istart = STACK_IMPORT(i - 1)
         inum   = STACK_IMPORT(i) - istart
         dest = NEIBPE(i)

         CALL MPI_ISEND(export_global_buffer(istart + 1), inum, LIS_MPI_INTEGER, dest, 11, & 
              &         SOLVER_COMM, req1(i), ierr)
      END DO

      IF(n_pes > 0) CALL MPI_WAITALL (n_pes, req2, sta2, ierr)
      IF(neib  > 0) CALL MPI_WAITALL (neib,  req1, sta1, ierr)

      base_index = GIN_aggregate(my_rank)    
      STACK_EXPORT(0) = 0
      DO i = 1, neib
         dest = NEIBPE(i)

         DO j = 1, n_pes
            if(pes(j) == dest) exit
         END DO

         IF(j > n_pes) THEN
            !C there is no nodes to send.
            STACK_EXPORT(i) = STACK_EXPORT(i - 1)
         ELSE
            STACK_EXPORT(i) = STACK_EXPORT(i - 1) + in_recv_buf(j) - in_recv_buf(j - 1)
            l = STACK_EXPORT(i - 1) - in_recv_buf(j - 1)
            DO k = in_recv_buf(j - 1) + 1, in_recv_buf(j)
               NOD_EXPORT(l + k) = recv_buf(k) - base_index - ZERO_ORIGIN

#ifdef DEBUG
               IF(pes(j)==-1) STOP "error!!! in make_export_table"
               IF(NOD_EXPORT(k+l)+ZERO_ORIGIN<0 .OR.NOD_EXPORT(k+l)+ZERO_ORIGIN>GIN_aggregate(my_rank+1)-base_index) THEN
                  WRITE(*,*) "unexpected error in make_export_nod():", NOD_EXPORT(k+l)+ZERO_ORIGIN,my_rank
                  STOP
               END IF
#endif            
               
            END DO
            !C mark which means un-needed 
            pes(j) = -1
         END IF
      END DO

      DO i = 1, n_pes
         IF(pes(i) /= -1) THEN
            dest = pes(i)
            NEIBPETOT = NEIBPETOT + 1
            NEIBPE(NEIBPETOT) = dest
            STACK_IMPORT(NEIBPETOT) = STACK_IMPORT(NEIBPETOT - 1)
            STACK_EXPORT(NEIBPETOT) = STACK_EXPORT(NEIBPETOT - 1) + in_recv_buf(i) - in_recv_buf(i - 1)
            l = STACK_EXPORT(NEIBPETOT - 1) - in_recv_buf(i - 1)
            DO k = in_recv_buf(i - 1) + 1, in_recv_buf(i)
               NOD_EXPORT(k + l) = recv_buf(k) - base_index-ZERO_ORIGIN
               
#ifdef DEBUG
               IF(NOD_EXPORT(k+l)+ZERO_ORIGIN<0 .OR. NOD_EXPORT(k+l)+ZERO_ORIGIN>GIN_aggregate(my_rank+1)-base_index) THEN
                  WRITE(*,*) "unexpected error in make_export_nod():" ,NOD_EXPORT(k+l)+ZERO_ORIGIN,my_rank
                  STOP
               END IF
#endif            
            END DO
         END IF
      END DO

      DEALLOCATE(pes)    
      DEALLOCATE(in_recv_buf)    
      DEALLOCATE(recv_buf)    

      DEALLOCATE(req1)
      DEALLOCATE(req2)
      DEALLOCATE(sta2)
      DEALLOCATE(sta1)
    END SUBROUTINE make_export_NOD
    
    
  END SUBROUTINE make_rap

#ifdef PSPASES
  SUBROUTINE make_PSPSE_MAT(Temp_CN, Temp_N, Temp_V, TEMP_COLSIZE, GIN_aggregate, &
       & aggregate_table_array, aggregate_number_in_table, PE_list_size, PE_list,&
       & NP, N, gbase, SOLVER_COMM)
    USE data_structure_for_AMG
    IMPLICIT NONE
    
    !C gbase is the last aggregate number on the proc of my_rank-1.
    INTEGER(kind=kint),                 INTENT(in) :: PE_list_size,NP
    INTEGER(kind=kint),                 INTENT(in) :: N,gbase,SOLVER_COMM
    INTEGER(kind=kint), DIMENSION(:,:), INTENT(in) :: aggregate_table_array,PE_list
    INTEGER(kind=kint), DIMENSION(:),   INTENT(in) :: aggregate_number_in_table

    INTEGER(kind=kint),                INTENT(in)  :: TEMP_COLSIZE
    INTEGER(kind=kint), DIMENSION(:,:),INTENT(in)  :: Temp_CN
    REAL(kind=kreal),   DIMENSION(:,:),INTENT(in)  :: Temp_V
    INTEGER(kind=kint), DIMENSION(:),  INTENT(in)  :: Temp_N
    INTEGER(kind=kint), DIMENSION(0:),  INTENT(in) :: GIN_aggregate
    
    INTEGER(kind=kint), ALLOCATABLE,DIMENSION(:)   :: local_global_vector,ainds
    REAL(kind=kreal),ALLOCATABLE, DIMENSION(:)     :: avals
    INTEGER(kind=kint), ALLOCATABLE,DIMENSION(:,:) :: aptrs
    
    INTEGER(kind=kint) :: i,j,k,col,count,l,m,ierr

    ALLOCATE(local_global_vector(NP))
    
    local_global_vector=0

    !C PSPESE's index starts from 0
    DO i=1,N
       local_global_vector(i)=gbase+i-1
    END DO

    DO i=1,PE_list_size
       DO j=1,aggregate_number_in_table(i)
          !C local
          k=aggregate_table_array(j*2,i)
          !C global
          local_global_vector(k)=aggregate_table_array(j*2-1,i)-1
       END DO
    END DO


    !C TempCN,Temp_N,Temp_V don't have lower matrix
    j = 0
    DO i=1,N
       j = j+Temp_N(i)
    END DO
    j=2*j-N
    ALLOCATE(ainds(j))
    ALLOCATE(avals(j))
    ALLOCATE(aptrs(2,N))

    aptrs=0
    
    DO i=1,N
       DO j=1,TEMP_COLSIZE
          col=Temp_CN(j,i)
!!$          if(col>i)then
          IF(col>i  .AND. EPS_COARSEST<abs(Temp_V(j,i)))THEN
!!$          if(col>i  .and. -EPS_COARSEST>Temp_V(j,i))then
             !C(k,i)
             !C upper part
             aptrs(2,i)=aptrs(2,i)+1
             !C lower part
             IF(col<=N) aptrs(2,col)=aptrs(2,col)+1
          END IF
       END DO
       !C diagonal part
       aptrs(2,i)=aptrs(2,i)+1
    END DO

    aptrs(1,1)=1
    DO i=2,N
       aptrs(1,i)=aptrs(1,i-1)+aptrs(2,i-1)
    END DO

    aptrs(2,:)=0

    DO i=1,N
       DO j=1,TEMP_COLSIZE
          col=Temp_CN(j,i)
          IF(col>i.AND. EPS_COARSEST<abs(Temp_V(j,i)))THEN
!!$          if(col>i.and. -EPS_COARSEST>Temp_V(j,i))then
             k=aptrs(1,i)+aptrs(2,i)
             aptrs(2,i)=aptrs(2,i)+1
             ainds(k)=local_global_vector(col)
             avals(k)=Temp_V(j,i)

             IF(col<=N) THEN
                k=aptrs(1,col)+aptrs(2,col)
                aptrs(2,col)=aptrs(2,col)+1
                ainds(k)=local_global_vector(i)
                avals(k)=Temp_V(j,i)
             END IF
          ELSE IF(col==i)THEN
             k=aptrs(1,i)+aptrs(2,i)
             aptrs(2,i)=aptrs(2,i)+1
             ainds(k)=local_global_vector(i)
             avals(k)=Temp_V(j,i)
          END IF
       END DO
    END DO

    !C block size
    ioptions(1)=64
    !C symmetric
    ioptions(2)=0
    !C sort on each row
    ioptions(3)=1
    ioptions(4)=0
    CALL DPSPACEF(GIN_aggregate,aptrs,ainds,avals,ioptions,doptions,&
         & pspcommf,SOLVER_COMM)
    
    !C parametar setting for DPSPACET
    ioptions=0
    ioptions(1) = 1
    ldb=NP
    ldx=NP
    nrhs=1
    DEALLOCATE(ainds)
    DEALLOCATE(avals)
    DEALLOCATE(aptrs)
    DEALLOCATE(local_global_vector)
    
    ALLOCATE(coarsest_X(NP,1))
    ALLOCATE(coarsest_B(NP,1))

  END SUBROUTINE make_PSPSE_MAT
#endif
  
  SUBROUTINE left_looking_chorescky(lower_mat,N)
    IMPLICIT NONE
#include "precision.inc"        
    
    REAL(kind=kreal), INTENT(inout)   :: lower_mat(:)
    INTEGER(kind=kint), INTENT(in)    :: N
    INTEGER(kind=kint)                :: i,j,k,rj,rk
    
    REAL(kind=kreal)                  :: v
    !C a(r,c) => lower_mat(r-c+1+ (2N-c+2)*(c-1)/2)
    !C lower_mat(r+N*(c-1)-c(c-1)/2)
    DO j=1,N
       rj=N*(j-1)-j*(j-1)/2
       DO k=1,j-1
          rk=N*(k-1)-k*(k-1)/2
          v=lower_mat(j+rk)
          DO i=j,N
             !C aij=aij - aik*akj
             lower_mat(i+rj)=lower_mat(i+rj)-lower_mat(i+rk)*v
          END DO
       END DO
       
       v=SQRT(lower_mat(j+rj))
       lower_mat(j+rj)=1/v
       DO k=j+1,N
          lower_mat(k+rj)=lower_mat(k+rj)/v
       END DO
    END DO
  END SUBROUTINE left_looking_chorescky
  
  

  !C symmetric and NP times NP matrix
  SUBROUTINE smooth_aggregate (PNI, NI, omega, LEVEL_NO, in_aggregates_result_size, &
       & in_aggregates_result, aggregates_result, Temp_P, Temp_P_row_size, tilde_D, &
       & RATE_OF_SPACE_i)
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint ), INTENT(in) :: PNI(0:)
    INTEGER(kind=kint ), INTENT(in) :: NI(:)
    REAL   (kind=kreal), INTENT(in) :: omega
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO, RATE_OF_SPACE_i
    INTEGER(kind=kint ), INTENT(in) :: in_aggregates_result_size
    REAL   (kind=kreal), INTENT(in) :: tilde_D(:)
    
    INTEGER(kind=kint) :: in_aggregates_result(0:)
    INTEGER(kind=kint) :: aggregates_result(:)
    TYPE(row_node)     :: Temp_P(:,:)
    INTEGER(kind=kint) :: Temp_P_row_size(:)

    REAL   (kind=kreal), ALLOCATABLE :: Temp(:,:)

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

    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP
    N  = HIERARCHICAL_DATA(LEVEL_NO-1) % N 

    ALLOCATE(Temp_int(NP))

    D   => HIERARCHICAL_DATA(LEVEL_NO-1) % D
    INL => HIERARCHICAL_DATA(LEVEL_NO-1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO-1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO-1) % AU

    Temp_P_row_size(1:NP) = 0


    
    !C-- Temp_int() has aggregate NO for each node.
    Temp_int = 0
    DO k = 1, in_aggregates_result_size
       is = in_aggregates_result(k-1)+1
       ie = in_aggregates_result(k)
       DO i = is,ie
          row_no = aggregates_result(i)
          Temp_int(row_no) = k
       END DO
    END DO

    !C-- Temp_P's space is calculated for smoothed aggregations
    !C   and input column and value( initial values 0 )
    DO k = 1, in_aggregates_result_size

       is = in_aggregates_result(k - 1) + 1
       ie = in_aggregates_result(k)
       DO i = is,ie
          row_no   = aggregates_result(i)
!!$          j        = 1
!!$          row_size = Temp_P_row_size(row_no)
!!$          DO WHILE(j <= row_size .AND. Temp_P(j, row_no) % column /= k)
!!$             j = j + 1
!!$          END DO

          row_size = Temp_P_row_size(row_no)
          DO j = 1, row_size
             if(Temp_P(j, row_no) % column == k) exit
          end DO
          IF(j > row_size) THEN
             Temp_P_row_size(row_no)   = j

#ifdef ALLOCATION_CHECK
             if(j > RATE_OF_SPACE_i) then
                write(*,*) "enlarge RATE_OF_SPACE"
                stop 'error in smooth_aggregate():'
             end if
#endif    

             Temp_P(j, row_no) % column = k
             Temp_P(j, row_no) % value  = 0.0

          END IF

#ifdef SMOOTH_TILDE_A
          connected_node_start = PNI(2*row_no-2)+1
          connected_node_end   = PNI(2*row_no-1)
#else
          connected_node_start = INL(row_no-1)+1
          connected_node_end   = INL(row_no)          
#endif
          DO connected_node = connected_node_start,connected_node_end
#ifdef SMOOTH_TILDE_A
             A_col = IAL(NI(connected_node)) + ZERO_ORIGIN 
#else
             A_col = IAL(connected_node) + ZERO_ORIGIN
#endif
!!$             j = 1
!!$             row_size = Temp_P_row_size(A_col)
!!$             DO WHILE(j <= row_size .AND. Temp_P(j, A_col) % column /= k)
!!$                j = j+1
!!$             END DO

             row_size = Temp_P_row_size(A_col)
             DO j = 1, row_size
                if(Temp_P(j, A_col) % column == k) exit
             END DO
             IF(j > row_size) THEN
                Temp_P_row_size(A_col) = j

#ifdef ALLOCATION_CHECK
                if(j > RATE_OF_SPACE_i) then
                   write(*,*) "enlarge RATE_OF_SPACE"
                   stop 'error in smooth_aggregate():'
                end if
#endif    
                Temp_P(j, A_col) % column = k
                Temp_P(j, A_col) % value  = 0.0

             END IF
          END DO

#ifdef SMOOTH_TILDE_A
          connected_node_start = PNI(2*row_no-1)+1
          connected_node_end = PNI(2*row_no)
#else
          connected_node_start = INU(row_no-1)+1
          connected_node_end = INU(row_no)          
#endif
          DO connected_node = connected_node_start, connected_node_end
#ifdef SMOOTH_TILDE_A             
             A_col = IAU(NI(connected_node)) + ZERO_ORIGIN
#else
             A_col = IAU(connected_node) + ZERO_ORIGIN
#endif
!!$             j = 1
!!$             row_size = Temp_P_row_size(A_col)
!!$             DO WHILE(j <= row_size .AND. Temp_P(j,A_col) % column /= k)
!!$                j = j+1
!!$             END DO

             row_size = Temp_P_row_size(A_col)
             DO j = 1, row_size
                if(Temp_P(j, A_col) % column == k) exit
             END DO
             IF(j > row_size) THEN
                Temp_P_row_size(A_col) = j

#ifdef ALLOCATION_CHECK
                if(j > RATE_OF_SPACE_i) then
                   write(*,*) "enlarge RATE_OF_SPACE"
                   stop 'error in smooth_aggregate():'
                end if
#endif    

                Temp_P(j, A_col) % column = k
                Temp_P(j, A_col) % value  = 0.0
                
             END IF
          END DO
       END DO
    END DO
    
    !C-- smooth aggregates
    DO A_row = 1, NP
       row_size = Temp_P_row_size(A_row)
       IF(row_size == 0) CYCLE

       !C  element's  value of the aggregate is 1 in array P      
#ifdef SMOOTH_TILDE_A             
       isL = PNI(2 * A_row - 2) + 1
       ieL = PNI(2 * A_row - 1)
#else
       isL = INL(A_row - 1) + 1
       ieL = INL(A_row)
#endif
       DO i = isL, ieL
#ifdef SMOOTH_TILDE_A             
          l = NI(i)
#else
          l = i
#endif
          A_col = IAL(l) + ZERO_ORIGIN
          P_col = Temp_int(A_col)
          IF(P_col == 0) CYCLE
          

!!$          k = 1
!!$          DO WHILE(Temp_P(k,A_row) % column /= P_col)
!!$             k = k + 1
!!$          END DO

          DO k = 1, Temp_P_row_size(A_row)
             if(Temp_P(k, A_row) % column == P_col) exit
          end DO

#ifdef DEBUG
          IF(k > row_size) then
             write(*,*) A_row,"1"
             STOP "unexpected error in smooth_aggregate()"
          end IF
#endif
          
          Temp_P(k,A_row) % value = Temp_P(k,A_row) % value - AL( l )
       END DO
#ifdef SMOOTH_TILDE_A             
       isU = PNI(2*A_row-1)+1
       ieU = PNI(2*A_row)
#else
       isU = INU(A_row-1)+1
       ieU = INU(A_row)
#endif

       DO i = isU,ieU
#ifdef SMOOTH_TILDE_A             
          l = NI(i)
#else
          l = i
#endif
          A_col = IAU(l) + ZERO_ORIGIN
          P_col = Temp_int(A_col)
          IF(P_col == 0) CYCLE

!!$          k = 1
!!$          DO WHILE(Temp_P(k,A_row) % column /= P_col)
!!$             k = k + 1
!!$          END DO
          DO k = 1, Temp_P_row_size(A_row)
             if(Temp_P(k, A_row) % column == P_col) exit
          end DO
          
#ifdef DEBUG
          IF(k > row_size) then
             write(*,*) A_row, row_size, k,Temp_P(1,A_row) % column ,P_col,"2"
             STOP "unexpected error in smooth_aggregate()"
          end IF
#endif
          Temp_P(k,A_row) % value = Temp_P(k,A_row) % value - AU( l )
       END DO
    END DO
    
    !C-- diagonal part
    DO A_row = 1, NP
       DO i = 1, Temp_P_row_size(A_row)
#ifdef SMOOTH_TILDE_A                            
          Temp_P(i, A_row) % value = Temp_P(i, A_row) % value * omega / &
               & tilde_D(A_row)
#else
          Temp_P(i, A_row) % value = Temp_P(i, A_row) % value * omega / &
               & D(A_row)
#endif
       END DO
    END DO


    omega2 = 1.-omega
    DO row_no = 1,N
       P_col = Temp_int(row_no)
       IF(P_col == 0) CYCLE



!!$       j = 1
!!$       DO WHILE(P_col /= Temp_P(j,row_no) % column)
!!$          j = j + 1
!!$       END DO

       DO j = 1, Temp_P_row_size(row_no)
          if(Temp_P(j, row_no) % column == P_col) exit
       end DO
       

#ifdef DEBUG
       IF(j > Temp_P_row_size(row_no)) then
          write(*,*)row_no,"3"
          STOP "unexpected error in smooth_aggreggate()"
       end IF
#endif
       Temp_P(j, row_no) % value = Temp_P(j, row_no) % value + omega2
    END DO
    !C-- end: smooth aggregates
    
    DEALLOCATE(Temp_int)
  END SUBROUTINE smooth_aggregate


  !C assumes symmetric matrix N times NP
  SUBROUTINE smooth_aggregate_old (PNI, NI, omega, LEVEL_NO, in_aggregates_result_size, &
       & in_aggregates_result, aggregates_result, Temp_P, Temp_P_row_size, tilde_D,     &
       & RATE_OF_SPACE_i, flnp, ni_sz, aggre_sz)
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint ), INTENT(in) :: flnp,ni_sz,aggre_sz
    INTEGER(kind=kint ), INTENT(in) :: PNI(0:2*flnp)
    INTEGER(kind=kint ), INTENT(in) :: NI(1:ni_sz)
    REAL   (kind=kreal), INTENT(in) :: omega
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO, RATE_OF_SPACE_i
    INTEGER(kind=kint ), INTENT(in) :: in_aggregates_result_size
    REAL   (kind=kreal), INTENT(in) ::  tilde_D(1:flnp)
    
    INTEGER(kind=kint) :: in_aggregates_result(0:in_aggregates_result_size)
    INTEGER(kind=kint) :: aggregates_result(1:aggre_sz)
    TYPE(row_node)     :: Temp_P(1:RATE_OF_SPACE_i,1:flnp)
    INTEGER(kind=kint) :: Temp_P_row_size(1:flnp)

    REAL   (kind=kreal),POINTER, DIMENSION(:) ::  D
    REAL   (kind=kreal),POINTER, DIMENSION(:) ::  AU
    REAL   (kind=kreal),POINTER, DIMENSION(:) ::  AL
    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  INU
    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  IAU
    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  INL
    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  IAL
    
    REAL   (kind=kreal) :: omega2, v
    INTEGER(kind=kint ) :: s, j, i, k, l, row, index, is, ie, NP
    INTEGER(kind=kint ) :: connected_node, connected_node_start, connected_node_end
    INTEGER(kind=kint ) :: P_CN, count, node_num


    
    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP
    s = HIERARCHICAL_DATA(LEVEL_NO-1) % N 

    D   => HIERARCHICAL_DATA(LEVEL_NO-1) % D
    INL => HIERARCHICAL_DATA(LEVEL_NO-1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO-1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO-1) % AU

    !C R=>HIERARCHICAL_DATA(LEVEL_NO) % R    
    
    Temp_P_row_size(1:NP)=0

    !    write(*,*) s,NP
    
    !C Dumped yacobi method
    DO k = 1, in_aggregates_result_size

       !C  element's  value of the aggregate is 1 in array P      
       is = in_aggregates_result(k - 1) + 1
       ie = in_aggregates_result(k)
       DO i = is, ie
          P_CN = aggregates_result(i)

          !C node in the aggregate
          DO count = 1, Temp_P_row_size(P_CN)
             if(Temp_P(count, P_CN) % column == k) exit
          END DO
          if(count > Temp_P_row_size(P_CN)) then
             Temp_P_row_size(P_CN) = count
             Temp_P(count, P_CN) % value = 0.
             Temp_P(count, P_CN) % column = k
          end if

          !C lower part
#ifdef SMOOTH_TILDE_A          
          connected_node_start = PNI(2 * P_CN - 2) + 1
          connected_node_end   = PNI(2 * P_CN - 1)
#else
          connected_node_start = INL(P_CN - 1) + 1
          connected_node_end   = INL(P_CN)
#endif
          DO connected_node = connected_node_start, connected_node_end


#ifdef SMOOTH_TILDE_A             
             node_num = IAL(NI(connected_node)) + ZERO_ORIGIN
             v = AL(NI(connected_node))
#else
             node_num = IAL(connected_node) + ZERO_ORIGIN
             v = AL(connected_node)
#endif
             DO count = 1, Temp_P_row_size(node_num)
                if(Temp_P(count, node_num) % column == k) exit
             END DO
             IF(count > Temp_P_row_size(node_num)) THEN
                Temp_P_row_size(node_num) = count
                Temp_P(count, node_num) % value  = -v
                Temp_P(count, node_num) % column = k
             else
                Temp_P(count, node_num) % value  = Temp_P(count, node_num) % value - v
             END IF
                             
          END DO

          !C upper part
#ifdef SMOOTH_TILDE_A          
          connected_node_start = PNI(2 * P_CN - 1) + 1
          connected_node_end = PNI(2 * P_CN)
#else
          connected_node_start = INU(P_CN - 1) + 1
          connected_node_end = INU(P_CN)
#endif
          
          DO connected_node = connected_node_start, connected_node_end

#ifdef SMOOTH_TILDE_A
             node_num = IAU(NI(connected_node)) + ZERO_ORIGIN
             v = AU(NI(connected_node))
#else
             node_num = IAU(connected_node) + ZERO_ORIGIN
             v = AU(connected_node)
#endif
             DO count = 1, Temp_P_row_size(node_num)
                if(Temp_P(count, node_num) % column == k) exit
             END DO
             if(count > Temp_P_row_size(node_num)) then
                Temp_P_row_size(node_num) = count
                Temp_P(count, node_num) % value  = -v
                Temp_P(count, node_num) % column = k
             else
                Temp_P(count, node_num) % value = Temp_P(count, node_num) % value - v
             END IF
          END DO
       END DO
    END DO


    !C-- diagonal part
    DO j = 1, NP
       DO i = 1, Temp_P_row_size(j)
#ifdef SMOOTH_TILDE_A                            
          Temp_P(i, j) % value = Temp_P(i, j) % value * omega / &
               & tilde_D(j)
#else
          Temp_P(i, j) % value = Temp_P(i, j) % value * omega / &
               & D(j)
#endif
       END DO
    END DO

    omega2 = 1. - omega
    DO k = 1, in_aggregates_result_size
       
       !C  element's  value of the aggregate is 1 in array P      
       is = in_aggregates_result(k - 1) + 1
       ie = in_aggregates_result(k)

       DO i = is, ie
          P_CN = aggregates_result(i)
          DO count = 1, Temp_P_row_size(P_CN)
             if(Temp_P(count, P_CN) % column == k) exit
          END DO
          if(count > Temp_P_row_size(P_CN)) then
             stop 'unexpected error in smooth_aggregate()'
          else
             Temp_P(count, P_CN) % value = Temp_P(count, P_CN) % value + omega2
          end if
       END DO
    END DO
    
  END SUBROUTINE smooth_aggregate_old

  !C This routine uses symmetricity structure of PNI in order to treat 
  !C unstructured problem matrix.
  !C This routine uses PNI, NI and no option other than SMOOTH_TILDE_A
  SUBROUTINE smooth_aggregate_unsym (PNI, NI, omega, LEVEL_NO, in_aggregates_result_size, &
       & in_aggregates_result, aggregates_result, Temp_P, Temp_P_row_size, tilde_D,       &
       & RATE_OF_SPACE_i, EX_INL, EX_IAL, EX_AL)
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint ), INTENT(in) :: PNI(0:)
    INTEGER(kind=kint ), INTENT(in) :: NI(:)
    REAL   (kind=kreal), INTENT(in) :: omega
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO, RATE_OF_SPACE_i
    INTEGER(kind=kint ), INTENT(in) :: in_aggregates_result_size
    REAL   (kind=kreal), INTENT(in) ::  tilde_D(:)
    INTEGER(kind=kint ), INTENT(in) :: EX_INL (0:), EX_IAL (:)
    REAL   (kind=kreal), INTENT(in) :: EX_AL  (:)
    
    
    INTEGER(kind=kint) :: in_aggregates_result(0:)
    INTEGER(kind=kint) :: aggregates_result(:)
    TYPE(row_node)     :: Temp_P(:,:)
    INTEGER(kind=kint) :: Temp_P_row_size(:)
    INTEGER(kind=kint), allocatable :: Temp_agr_num(:)

    REAL   (kind=kreal), POINTER ::  D(:)
    REAL   (kind=kreal), POINTER ::  AU(:)
    REAL   (kind=kreal), POINTER ::  AL(:)
    INTEGER(kind=kint ), POINTER ::  INU(:)
    INTEGER(kind=kint ), POINTER ::  IAU(:)
    INTEGER(kind=kint ), POINTER ::  INL(:)
    INTEGER(kind=kint ), POINTER ::  IAL(:)
    
    REAL   (kind=kreal) :: omega2, v
    INTEGER(kind=kint ) :: N, i, j, k, l, row, index, is, ie, NP, c
    INTEGER(kind=kint ) :: connected_node, connected_node_start, connected_node_end
    INTEGER(kind=kint ) :: P_CN, count, node_num


    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP
    N  = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    
    D   => HIERARCHICAL_DATA(LEVEL_NO - 1) % D 
    INL => HIERARCHICAL_DATA(LEVEL_NO - 1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO - 1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AU
    
    !C R=>HIERARCHICAL_DATA(LEVEL_NO) % R    

    !C assign the aggregate number in Temp_agr_num
    allocate(Temp_agr_num(N))
    Temp_agr_num = 0
    DO i = 1, in_aggregates_result_size
       DO j = in_aggregates_result(i - 1) + 1, in_aggregates_result(i)
          k = aggregates_result(j)
          Temp_agr_num(k) = i
       END DO
    END DO

    Temp_P_row_size(1:NP) = 0
    
    !    write(*,*) s,NP

    !C determine nonzero structure which utilizes the symmetricity of SMOOTH_TILDE_A
    DO k = 1, in_aggregates_result_size

       !C  element's  value of the aggregate is 1 in array P      
       is = in_aggregates_result(k - 1) + 1
       ie = in_aggregates_result(k)
       
       DO i = is, ie
          P_CN = aggregates_result(i)
          !C node in the aggregate
          DO count = 1, Temp_P_row_size(P_CN)
             if(Temp_P(count, P_CN) % column == k) exit
          END DO
          if(count > Temp_P_row_size(P_CN)) then
             Temp_P_row_size(P_CN) = count
             Temp_P(count, P_CN) % column = k
             Temp_P(count, P_CN) % value = 0.
          end if
          !C lower part
          connected_node_start = PNI(2 * P_CN - 2) + 1
          connected_node_end   = PNI(2 * P_CN - 1)
          DO connected_node = connected_node_start, connected_node_end
             node_num = IAL(NI(connected_node)) + ZERO_ORIGIN
             DO count = 1, Temp_P_row_size(node_num)
                if(Temp_P(count, node_num) % column == k) exit
             END DO
             IF(count > Temp_P_row_size(node_num)) THEN
                Temp_P_row_size(node_num) = count
                Temp_P(count, node_num) % column = k
                Temp_P(count, node_num) % value = 0.
             END IF
          END DO
          !C upper part
          connected_node_start = PNI(2 * P_CN - 1) + 1
          connected_node_end = PNI(2 * P_CN)
          DO connected_node = connected_node_start, connected_node_end
             node_num = IAU(NI(connected_node)) + ZERO_ORIGIN
             DO count = 1, Temp_P_row_size(node_num)
                if(Temp_P(count, node_num) % column == k) exit
             END DO
             if(count > Temp_P_row_size(node_num)) then
                Temp_P_row_size(node_num) = count
                Temp_P(count, node_num) % column = k
                Temp_P(count, node_num) % value = 0.
             END IF
          END DO
       END DO
    END DO
    
    !C Dumped yacobi method
    DO i = 1, N
       !C lower part
       DO j = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
          v = AL(NI(j))
          k = Temp_agr_num(IAL(NI(j))+ZERO_ORIGIN)

          DO count = 1, Temp_P_row_size(i)
             if(Temp_P(count, i) % column == k) exit
          END DO
          IF(count > Temp_P_row_size(i)) THEN
             stop '0:error in smooth_aggregate_unsym'
          else
             Temp_P(count, i) % value  = Temp_P(count, i) % value - v
          END IF
       END DO

       !C upper part
       DO j = PNI(2 * i - 1) + 1, PNI(2 * i)
          v = AU(NI(j))
          c = IAU(NI(j)) + ZERO_ORIGIN
          if(c <= N) then
             k = Temp_agr_num(c)

             DO count = 1, Temp_P_row_size(i)
                if(Temp_P(count, i) % column == k) exit
             END DO
             IF(count > Temp_P_row_size(i)) THEN
                stop '1:error in smooth_aggregate_unsym'
             else
                Temp_P(count, i) % value  = Temp_P(count, i) % value - v
             END IF
          end if
       END DO
    END DO

    !C external lower part
    DO i = N + 1, NP
       DO j = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
          v = EX_AL(NI(j))
          c = EX_IAL(NI(j))+ZERO_ORIGIN
          if(c <= N) then
             k = Temp_agr_num(c)

             DO count = 1, Temp_P_row_size(i)
                if(Temp_P(count, i) % column == k) exit
             END DO

             IF(count > Temp_P_row_size(i)) THEN
                stop '2:error in smooth_aggregate_unsymb'
             else
                Temp_P(count, i) % value  = Temp_P(count, i) % value - v
             END IF
          END if
       end DO
    END DO

    !C-- diagonal part
    DO j = 1, NP
       DO i = 1, Temp_P_row_size(j)
          Temp_P(i, j) % value = Temp_P(i, j) % value * omega / &
               & tilde_D(j)
       END DO
    END DO

    omega2 = 1. - omega
    DO k = 1, in_aggregates_result_size
       
       !C  element's  value of the aggregate is 1 in array P      
       is = in_aggregates_result(k - 1) + 1
       ie = in_aggregates_result(k)

       DO i = is, ie
          P_CN = aggregates_result(i)
          DO count = 1, Temp_P_row_size(P_CN)
             if(Temp_P(count, P_CN) % column == k) exit
          END DO
          if(count > Temp_P_row_size(P_CN)) then
             stop 'unexpected error in smooth_aggregate()'
          else
             Temp_P(count, P_CN) % value = Temp_P(count, P_CN) % value + omega2
          end if
       END DO
    END DO

    deallocate(Temp_agr_num)
  END SUBROUTINE smooth_aggregate_unsym

  !C create the lower external part of the matrix. 
  !C Row numbers of that part are from N+1 to NP on COMM tables.
  SUBROUTINE lower_ext_matrix_crtn(EX_NPL, EX_INL, EX_IAL, EX_AL, LEVEL_NO, my_rank, SOLVER_COMM, PETOT)
    USE data_structure_for_AMG
    USE solver_SR2
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO, PETOT

    !C These three arrays are allocated
    INTEGER(kind=kint ) :: EX_NPL
    INTEGER(kind=kint ), POINTER :: EX_INL (:), EX_IAL (:)
    REAL   (kind=kreal), POINTER :: EX_AL  (:)

    INTEGER(kind=kint ) :: N, NP
    REAL   (kind=kreal), POINTER :: D  (:)
    REAL   (kind=kreal), POINTER :: AU (:)
    REAL   (kind=kreal), POINTER :: AL (:)
    INTEGER(kind=kint ), POINTER :: INU(:)
    INTEGER(kind=kint ), POINTER :: IAU(:)
    INTEGER(kind=kint ), POINTER :: INL(:)
    INTEGER(kind=kint ), POINTER :: IAL(:)
    !C communication data structures
    INTEGER(kind=kint ), INTENT(in) :: SOLVER_COMM
    INTEGER(kind=kint ), INTENT(in) :: my_rank

    INTEGER(kind=kint )          :: NEIBPETOT
    INTEGER(kind=kint ), POINTER :: NEIBPE    (:)
    INTEGER(kind=kint ), POINTER :: NOD_IMPORT(:), STACK_IMPORT(:)
    INTEGER(kind=kint ), POINTER :: NOD_EXPORT(:), STACK_EXPORT(:)
    

    INTEGER(kind=kint ), allocatable :: Temp_N_NP(:), in_send_buf(:), send_buf(:)
    REAL   (kind=kreal), allocatable :: send_buf_dbl(:)
!!$    INTEGER(kind=kint ), allocatable :: localglobal_vec(:), index(:)
    INTEGER(kind=kint ) :: i, j, k, l, m, r, c, send_buf_size
    INTEGER(kind=kint ) :: send_buf_base
!!$    INTEGER(kind=kint ), ALLOCATABLE  :: WSI(:), WRI(:)

    INTEGER(kind=kint ), allocatable :: in_recv_buf(:), recv_buf(:)
    REAL   (kind=kreal), allocatable :: recv_buf_dbl(:)
    INTEGER(kind=kint ) :: recv_buf_size, recv_buf_base

    INTEGER(kind=kint ) :: ierr, istart, inum
    INTEGER(kind=kint ), allocatable :: req1(:),req2(:),sta1(:,:),sta2(:,:)
    INTEGER(kind=kint ), allocatable :: req3(:),req4(:),sta3(:,:),sta4(:,:)

    NP = HIERARCHICAL_DATA(LEVEL_NO - 1) % NP     
    N  = HIERARCHICAL_DATA(LEVEL_NO - 1) % N 
    
    D   => HIERARCHICAL_DATA(LEVEL_NO - 1) % D
    INL => HIERARCHICAL_DATA(LEVEL_NO - 1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO - 1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AU
    NEIBPETOT    =  HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NEIBPETOT
    if(NEIBPETOT > 0) then
       NEIBPE       => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NEIBPE
       STACK_IMPORT => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % STACK_IMPORT
       NOD_IMPORT   => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NOD_IMPORT
       STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % STACK_EXPORT
       NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NOD_EXPORT
    end if

    
    

    !C- communication for determining global number of the nodes from N to NP only on COMM table.
!!$    allocate(localglobal_vec(NP), index(0:PETOT))
!!$
!!$    CALL MPI_ALLGATHER(N, 1, LIS_MPI_INTEGER, index(1), 1, LIS_MPI_INTEGER, SOLVER_COMM, ierr)
!!$
!!$    index(0) = 0
!!$    DO i = 1, PETOT
!!$       index(i) = index(i) + index(i - 1)
!!$    END DO
!!$    
!!$    DO i = 1, N
!!$       localglobal_vec(i) = i + index(my_rank)
!!$    END DO
!!$    
!!$    IF(PETOT > 1) THEN
!!$       ALLOCATE(WSI(STACK_EXPORT(NEIBPETOT) + 1))
!!$       ALLOCATE(WRI(STACK_IMPORT(NEIBPETOT) + 1))
!!$       CALL SOLVER_SEND_RECV2I                                            &
!!$            &   (NP, STACK_EXPORT(NEIBPETOT) + 1,                         &
!!$            &    STACK_IMPORT(NEIBPETOT) + 1, NEIBPETOT,                  &
!!$            &    NEIBPE, STACK_IMPORT, NOD_IMPORT,                        &
!!$            &    STACK_EXPORT, NOD_EXPORT, WSI, WRI, localglobal_vec,     &
!!$            &    SOLVER_COMM, my_rank)
!!$       DEALLOCATE(WSI)
!!$       DEALLOCATE(WRI)
!!$    END IF
    !C- communication for determining global number of the nodes from N to NP only on COMM table.
    
    !C- count up the communication size for each PE
       
    allocate(Temp_N_NP(NP - N), in_send_buf(STACK_EXPORT(NEIBPETOT)))
    Temp_N_NP = 0    
    in_send_buf = 0
    DO i = 1, NEIBPETOT
       
       !C numbering the order in RECV table
       DO j = STACK_IMPORT(i - 1) + 1, STACK_IMPORT(i)
          r = NOD_IMPORT(j)+ZERO_ORIGIN
          Temp_N_NP(r - N) = j - STACK_IMPORT(i - 1)
       END DO

       !C count up the nonzeros for each row in the transposed upper external matrix.
       DO j = STACK_EXPORT(i - 1) + 1, STACK_EXPORT(i)
          r = NOD_EXPORT(j)+ZERO_ORIGIN
          DO k = INU(r - 1) + 1, INU(r)
             c = IAU(k) + ZERO_ORIGIN
             if(c > N) then
                l = Temp_N_NP(c - N)
                if(l > 0) then
                   !C count up l
                   in_send_buf(j) = in_send_buf(j) + 1
                end if
             END if
          END DO
       END DO
       
       !C reset the Temp_N_NP
       DO j = STACK_IMPORT(i - 1) + 1, STACK_IMPORT(i)
          r = NOD_IMPORT(j)+ZERO_ORIGIN
          Temp_N_NP(r - N) = 0
       END DO
    END DO
    !C- count up the communication size for each PE


    !C- create the send_buf and send_buf_dbl
    !C indexing for each PE
    DO i = 1, NEIBPETOT
       DO j = STACK_EXPORT(i - 1) + 2, STACK_EXPORT(i)
          in_send_buf(j) = in_send_buf(j - 1) + in_send_buf(j)
       END DO
    END DO
    
    !C send_buf_size
    send_buf_size = 0
    DO i = 1, NEIBPETOT
       send_buf_size = send_buf_size + in_send_buf(STACK_EXPORT(i))
    END DO

    allocate(send_buf(send_buf_size), send_buf_dbl(send_buf_size))
    send_buf_base = 0
    DO i = 1, NEIBPETOT

       !C numbering the order in RECV table
       DO j = STACK_IMPORT(i - 1) + 1, STACK_IMPORT(i)
          r = NOD_IMPORT(j)+ZERO_ORIGIN
          Temp_N_NP(r - N) = j - STACK_IMPORT(i - 1)
       END DO

       !C shift the in_send_buf array
       DO j = STACK_EXPORT(i), STACK_EXPORT(i - 1) + 2, -1
          in_send_buf(j) = in_send_buf(j - 1)
       END DO
       in_send_buf(STACK_EXPORT(i - 1) + 1) = 0
       
       !C assign the nonzeros for each row in the upper external matrix.
       DO j = STACK_EXPORT(i - 1) + 1, STACK_EXPORT(i)
          r = NOD_EXPORT(j)+ZERO_ORIGIN
          DO k = INU(r - 1) + 1, INU(r)
             c = IAU(k) + ZERO_ORIGIN
             if(c > N) then
                l = Temp_N_NP(c - N)
                if(l > 0) then
                   m = in_send_buf(j) + 1
                   in_send_buf(j) = m
                   send_buf    (send_buf_base + m) = l
                   send_buf_dbl(send_buf_base + m) = AU(k)
                end if
             END if
          END DO
       END DO
       
       if(STACK_EXPORT(i) > STACK_EXPORT(i - 1)) then
          send_buf_base = send_buf_base + in_send_buf(STACK_EXPORT(i))
       end if
       
       !C reset the Temp_N_NP
       DO j = STACK_IMPORT(i - 1) + 1, STACK_IMPORT(i)
          r = NOD_IMPORT(j)+ZERO_ORIGIN
          Temp_N_NP(r - N) = 0
       END DO
    END DO
    !C- create the send_buf and send_buf_dbl
    
    !C   in_send_buf(1:STACK_EXPORT(NEIBPETOT)) stores each row's size of 
    !C the external upper part.
    !C   send_buf(1:send_buf_size) stores the numbers in IMPORT TABLE 
    !C in order to show the column number on the other PE.
    !C   send_buf_dbl(1:send_buf_size) stores the value of the nonzero element

    
    !C- sendrecv in_send_buf
    allocate(in_recv_buf(STACK_IMPORT(NEIBPETOT)))
    allocate(sta1(MPI_STATUS_SIZE, NEIBPETOT), sta2(MPI_STATUS_SIZE, NEIBPETOT))
    allocate(sta3(MPI_STATUS_SIZE, NEIBPETOT), sta4(MPI_STATUS_SIZE, NEIBPETOT))
    allocate(req1(NEIBPETOT), req2(NEIBPETOT))
    allocate(req3(NEIBPETOT), req4(NEIBPETOT))

    DO i = 1, NEIBPETOT
       istart = STACK_EXPORT(i - 1)
       inum = STACK_EXPORT(i) - istart
       CALL MPI_ISEND(in_send_buf(istart + 1), inum, LIS_MPI_INTEGER, NEIBPE(i), 0, &
            & SOLVER_COMM, req1(i), ierr)
    END DO
    DO i = 1, NEIBPETOT
       istart = STACK_IMPORT(i - 1)
       inum = STACK_IMPORT(i) - istart
       CALL MPI_IRECV(in_recv_buf(istart + 1), inum, LIS_MPI_INTEGER, NEIBPE(i), 0, &
            & SOLVER_COMM, req2(i), ierr)
    END DO
    call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
    call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)
    !C- sendrecv in_send_buf

    !C- sendrecv send_buf, send_buf_dbl
    recv_buf_size = 0
    DO i = 1, NEIBPETOT
       recv_buf_size = recv_buf_size + in_recv_buf(STACK_IMPORT(i))
    END DO
    allocate(recv_buf(recv_buf_size), recv_buf_dbl(recv_buf_size))
    
    send_buf_base = 0
    DO i = 1, NEIBPETOT
       istart = send_buf_base
       inum = in_send_buf(STACK_EXPORT(i))
       call MPI_ISEND(send_buf(istart + 1),     inum, LIS_MPI_INTEGER,          &
            & NEIBPE(i), 0, SOLVER_COMM, req1(i), ierr)
       call MPI_ISEND(send_buf_dbl(istart + 1), inum, MPI_DOUBLE_PRECISION, &
            & NEIBPE(i), 0, SOLVER_COMM, req2(i), ierr)
       send_buf_base = send_buf_base + inum
    END DO

    recv_buf_base = 0
    DO i = 1, NEIBPETOT
       istart = recv_buf_base
       inum = in_recv_buf(STACK_IMPORT(i))
       call MPI_IRECV(recv_buf(istart + 1),     inum, LIS_MPI_INTEGER,          &
            & NEIBPE(i), 0, SOLVER_COMM, req3(i), ierr)
       call MPI_IRECV(recv_buf_dbl(istart + 1), inum, MPI_DOUBLE_PRECISION, &
            & NEIBPE(i), 0, SOLVER_COMM, req4(i), ierr)
       recv_buf_base = recv_buf_base + inum
    END DO
    call MPI_WAITALL (NEIBPETOT, req3, sta3, ierr)
    call MPI_WAITALL (NEIBPETOT, req4, sta4, ierr)
    call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
    call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)
    
    deallocate(in_send_buf, send_buf, send_buf_dbl, Temp_N_NP)
    !C- sendrecv send_buf, send_buf_dbl

    !C- input the external 
    allocate(EX_INL(0:NP - N))
    EX_INL(0:NP-N) = 0

    DO i = 1, NEIBPETOT
       !C front node in IMPORT table from each PE
       j = STACK_IMPORT(i - 1) + 1
       r = NOD_IMPORT(j)+ZERO_ORIGIN - N
       EX_INL(r) = in_recv_buf(j)
       !C the other nodes in IMPORT table from each PE
       DO j = STACK_IMPORT(i - 1) + 2, STACK_IMPORT(i)
          r = NOD_IMPORT(j)+ZERO_ORIGIN - N
          EX_INL(r) = in_recv_buf(j) - in_recv_buf(j - 1)
       END DO
    END DO

    !C 1. EX_INL is indexed
    DO i = 1, NP - N
       EX_INL(i) = EX_INL(i) + EX_INL(i - 1)
    END DO
    
    !C 2. EX_IAL, EX_AL is allocated
    EX_NPL = EX_INL(NP - N)
    allocate(EX_IAL(EX_NPL), EX_AL(EX_NPL))


    
    !C 3. EX_INL is shifted
    DO i = NP - N, 1, -1
       EX_INL(i) = EX_INL(i - 1)
    END DO
    
    !C 4. assigning the data in EX_IAL, EX_AL with counting up in EX_INL
    recv_buf_base = 0
    DO i = 1, NEIBPETOT
       !C front node
       j = STACK_IMPORT(i - 1) + 1
       r = NOD_IMPORT(j)+ZERO_ORIGIN - N
       l = EX_INL(r)
       DO k = 1, in_recv_buf(j)
          EX_IAL(l + k) = NOD_EXPORT(recv_buf(recv_buf_base + k) + STACK_EXPORT(i - 1))
          EX_AL (l + k) = recv_buf_dbl(recv_buf_base + k)
       END DO
       EX_INL(r) = EX_INL(r) + in_recv_buf(j)

       !C rest nodes
       DO j = STACK_IMPORT(i - 1) + 2, STACK_IMPORT(i)
          r = NOD_IMPORT(j)+ZERO_ORIGIN - N
          l = EX_INL(r)
          DO k = in_recv_buf(j - 1) + 1, in_recv_buf(j)
             EX_IAL(l + k - in_recv_buf(j - 1)) = NOD_EXPORT(recv_buf(recv_buf_base + k) + STACK_EXPORT(i - 1))
             EX_AL (l + k - in_recv_buf(j - 1)) = recv_buf_dbl(recv_buf_base + k)
          END DO
          EX_INL(r) = EX_INL(r) + in_recv_buf(j) - in_recv_buf(j - 1)
       END DO
       recv_buf_base = recv_buf_base + in_recv_buf(STACK_IMPORT(i))
    END DO
    
    deallocate(in_recv_buf, sta1, sta2, sta3, sta4, req1, req2, req3, req4)
    deallocate(recv_buf, recv_buf_dbl)

  END SUBROUTINE lower_ext_matrix_crtn

  
  SUBROUTINE neighbors(S1, S2, NP, PNI, NI, theta, LEVEL_NO, node_index, tilde_D)
    USE data_structure_for_AMG
    IMPLICIT NONE
    
    INTEGER(kind=kint ) :: N, NP
    REAL   (kind=kreal), INTENT(in):: theta
    INTEGER(kind=kint ), INTENT(in) :: S1, S2
    INTEGER(kind=kint ), INTENT(out) :: PNI(0:S1)
    INTEGER(kind=kint ), INTENT(out) :: NI(1:S2)
    INTEGER(kind=kint ), INTENT(out) :: node_index(1:NP)
    
    INTEGER(kind=kint ), INTENT(in) :: LEVEL_NO
    REAL   (kind=kreal), POINTER :: D  (:)
    REAL   (kind=kreal), POINTER :: AU (:)
    REAL   (kind=kreal), POINTER :: AL (:)
    INTEGER(kind=kint ), POINTER :: INU(:)
    INTEGER(kind=kint ), POINTER :: IAU(:)
    INTEGER(kind=kint ), POINTER :: INL(:)
    INTEGER(kind=kint ), POINTER :: IAL(:)

    REAL(kind=kreal)  :: tilde_D(1:NP)

    INTEGER(kind=kint)  :: i,j,isU,ieU,isL,ieL,inod,nni,nni_before,size,count
    INTEGER(kind=kint)  :: NPL, NPU
    REAL   (kind=kreal) :: t

    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP     
    N  = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    
    NPL = HIERARCHICAL_DATA(LEVEL_NO-1) % NPL
    NPU = HIERARCHICAL_DATA(LEVEL_NO-1) % NPU
    
    D   => HIERARCHICAL_DATA(LEVEL_NO - 1) % D
    INL => HIERARCHICAL_DATA(LEVEL_NO - 1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO - 1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AU
    

    PNI(0) = 0

    nni    = 0

    DO i = 1, NP
       tilde_D(i) = D(i)
       node_index(i) = 0
    END DO

    


!!$    DO j = 1, NP
    DO j = 1, N

       isU = INU(j-1)+1
       ieU = INU(j  )
       isL = INL(j-1)+1
       ieL = INL(j  )

       nni_before=nni
       DO i= isL, ieL

          inod= IAL(i)+ZERO_ORIGIN
          t= AL(i)
          IF(t*t > abs(D(j)*D(inod)*theta*theta) .AND. D(j)*D(inod)*t<0) THEN
             nni=nni+1
             NI(nni)=i
          ELSE
!!$             if(t*D(j)<0) then 
             tilde_D(j) = tilde_D(j) + t
!!$             end if
          ENDIF
       ENDDO

       PNI(2*j-1)=nni
       DO i = isU, ieU
          inod = IAU(i) + ZERO_ORIGIN
          t = AU(i)
          IF(t*t > abs(D(j)*D(inod)*theta*theta) .AND. D(j)*D(inod)*t<0) THEN
             nni = nni + 1
             NI(nni) = i
          ELSE
!!$             if(t*D(j)<0) then 
             tilde_D(j) = tilde_D(j) + t
!!$             end if
          END IF
       ENDDO
       PNI(2*j) = nni

       size = nni - nni_before
       IF(size == 0 .AND. j <= N) node_index(j) = -1
    ENDDO

  END SUBROUTINE neighbors
  
  
  !C this subroutine considers the case of unsymmetric problem matrix.
  !C But symmetricity of nonzero structure of the matrix is assumed.
  SUBROUTINE neighbors_unsym(EX_INL, EX_IAL, EX_AL, PNI, NI, theta, LEVEL_NO,  &
       & node_index, tilde_D)
    USE data_structure_for_AMG
    IMPLICIT NONE

    INTEGER(kind=kint ) :: N, NP
    REAL   (kind=kreal), INTENT(in) :: theta
    INTEGER(kind=kint ), INTENT(in) :: EX_INL (0:), EX_IAL (:)
    REAL   (kind=kreal), INTENT(in) :: EX_AL  (:)

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

    INTEGER(kind=kint) :: i, j, k, isU, ieU, isL, ieL, inod, nni
    REAL   (kind=kreal):: tl, tu

    NP   = HIERARCHICAL_DATA(LEVEL_NO-1) % NP     
    N    = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    
    D   => HIERARCHICAL_DATA(LEVEL_NO-1) % D
    INL => HIERARCHICAL_DATA(LEVEL_NO-1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO-1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO-1) % AU
    
    PNI(0) = 0
    nni    = 0
    
    DO i = 1, NP
       tilde_D(i) = D(i)
       node_index(i) = 0
    END DO
    
    DO j = 1, N
       isL = INL(j - 1) + 1; ieL = INL(j)
       DO i = isL, ieL
          inod = IAL(i) + ZERO_ORIGIN
          tl = AL(i)
          !C determining the cordinate position of upper matrix
          
          isU = INU(inod - 1) + 1; ieU = INU(inod)
          call cordinating_upper_lower_mat(NI, nni, D, theta, tilde_D, IAU, AU, isU, ieU, &
               & i, j, tl, inod)
       ENDDO
       PNI(2 * j - 1) = nni
       
       isU = INU(j - 1) + 1; ieU = INU(j)
       DO i = isU, ieU
          inod = IAU(i) + ZERO_ORIGIN
          tu = AU(i)

          !C determining the cordinate position of lower matrix          
          if(inod > N) then
             !C This part uses external lower part matrix
             isL = EX_INL(inod - N - 1) + 1; ieL = EX_INL(inod - N)
             call cordinating_upper_lower_mat(NI, nni, D, theta, tilde_D, EX_IAL, EX_AL, &
                  & isL, ieL, i, j, tu, inod)
          ELSE
             !C determining the cordinate position of lower matrix
             isL = INL(inod - 1) + 1; ieL = INL(inod)
             call cordinating_upper_lower_mat(NI, nni, D, theta, tilde_D, IAL, AL, isL, &
                  & ieL, i, j, tu, inod)
          ENDIF
       ENDDO
       PNI(2 * j) = nni
    END DO

    !C This part is for external lower submatrix
    DO j = N + 1, NP
       isL = EX_INL(j - N - 1) + 1
       ieL = EX_INL(j - N)
       DO i = isL, ieL
          inod = EX_IAL(i)+ZERO_ORIGIN
          tl = EX_AL(i)
          !C determining the cordinate position of upper matrix
          isU = INU(inod - 1) + 1; ieU = INU(inod)
          call cordinating_upper_lower_mat(NI, nni, D, theta, tilde_D, IAU, AU,  &
               & isU, ieU, i, j, tl, inod)
       ENDDO
       PNI(2 * j - 1) = nni
       PNI(2 * j) = nni
    END DO
    
    DO j = 1, N
       IF(PNI(2 * j) - PNI(2 * j - 2) == 0) node_index(j) = -1
    END DO
  
  contains
    !C search the cordinating nonzero element in matrix(cord_*A) and
    !C record chosen nodes in NI(:).
    subroutine cordinating_upper_lower_mat(NI, nni, D, theta, tilde_D, cord_IA, &
         & cord_A, isc, iec, i, r, tl, inod)
      implicit none
      
      INTEGER(kind=kint ), intent(in)    :: i, r, inod, isc, iec
      INTEGER(kind=kint ), intent(inout) :: nni,  NI(:)
      INTEGER(kind=kint ), intent(in)    :: cord_IA(:)
      REAL   (kind=kreal), intent(in)    :: cord_A(:), D(:), theta, tl
      REAL   (kind=kreal), intent(inout) :: tilde_D(:)

      REAL   (kind=kreal) :: tc
      INTEGER(kind=kint ) :: k 
      
      DO k = isc, iec
         IF(cord_IA(k)+ZERO_ORIGIN == r) EXIT
      END DO
      IF(k > iec) then
         tilde_D(r) = tilde_D(r) + tl
      ELSE
         tc = cord_A(k)
         IF(tl * tl > abs(D(r) * D(inod) * theta * theta) .AND. D(r) * D(inod) * tl < 0 .OR. &
              & tc * tc > abs(D(r) * D(inod) * theta * theta) .AND. D(r) * D(inod) * tc < 0) THEN
            nni = nni + 1
            NI(nni) = i
         ELSE
            tilde_D(r) = tilde_D(r) + tl
         ENDIF
      END IF
      
    END subroutine cordinating_upper_lower_mat
  END SUBROUTINE neighbors_unsym
  
END MODULE  data_creation_AMGCG
