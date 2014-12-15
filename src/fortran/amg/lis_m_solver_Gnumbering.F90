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
!C   * MODULE solver_Gnumbering
!C     CONTAINS
!C   * SUBROUTINE SOLVER_SR_REDUCE
!C   * SUBROUTINE check_Temp_P
!C   * SUBROUTINE solver_SR_aggregates
!C   ************************************************

MODULE solver_Gnumbering
CONTAINS
  SUBROUTINE  SOLVER_SR_REDUCE                                                  &
       &                ( N, dim, S, R, NEIBPETOT, NEIBPE, STACK_IMPORT,        &
       &                  NOD_IMPORT, STACK_EXPORT, NOD_EXPORT,&
       &                  WS, WR, X, SOLVER_COMM,my_rank,time_kind)
    
    IMPLICIT NONE
    INCLUDE  'mpif.h'
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif

    INTEGER(kind=kint), INTENT(in) ::  N, dim
    INTEGER(kind=kint), INTENT(in) ::  S
    INTEGER(kind=kint), INTENT(in) ::  R
    INTEGER(kind=kint), INTENT(in) ::  NEIBPETOT
    INTEGER(kind=kint) :: NEIBPE      (1:NEIBPETOT)
    INTEGER(kind=kint) :: STACK_IMPORT(0:NEIBPETOT)
    INTEGER(kind=kint) :: NOD_IMPORT  (1:R)
    INTEGER(kind=kint) :: STACK_EXPORT(0:NEIBPETOT)
    INTEGER(kind=kint) :: NOD_EXPORT  (1:S)
    REAL   (kind=kreal), INTENT(inout):: WS(1:dim)
    REAL   (kind=kreal), INTENT(inout):: WR(1:dim)

    REAL   (kind=kreal), INTENT(inout):: X(1:N)
    INTEGER(kind=kint), INTENT(in) :: SOLVER_COMM
    INTEGER(kind=kint), INTENT(in) :: my_rank

    INTEGER(kind=kint), DIMENSION(:,:), SAVE, ALLOCATABLE :: sta1
    INTEGER(kind=kint), DIMENSION(:,:), SAVE, ALLOCATABLE :: sta2
    INTEGER(kind=kint), DIMENSION(:  ), SAVE, ALLOCATABLE :: req1
    INTEGER(kind=kint), DIMENSION(:  ), SAVE, ALLOCATABLE :: req2  
    
    INTEGER(kind=kint) :: time_kind
    
    INTEGER(kind=kint) :: neib,istart,inum,k,ierr,i
    INTEGER(kind=kint), SAVE :: NFLAG=0
    INTEGER(kind=kint), SAVE :: EARLIER_NEIBPETOT=0

    !C
    !C-- INIT.
    IF (NFLAG.EQ.0) THEN
       ALLOCATE(sta1(MPI_STATUS_SIZE,NEIBPETOT))
       ALLOCATE(sta2(MPI_STATUS_SIZE,NEIBPETOT))
       ALLOCATE(req1(NEIBPETOT))
       ALLOCATE(req2(NEIBPETOT))
       NFLAG= 1
       EARLIER_NEIBPETOT=NEIBPETOT
    ELSE IF(EARLIER_NEIBPETOT<NEIBPETOT)THEN
       DEALLOCATE(sta1,sta2)
       DEALLOCATE(req1,req2)
       ALLOCATE(sta1(MPI_STATUS_SIZE,NEIBPETOT))
       ALLOCATE(sta2(MPI_STATUS_SIZE,NEIBPETOT))
       ALLOCATE(req1(NEIBPETOT))
       ALLOCATE(req2(NEIBPETOT))
       EARLIER_NEIBPETOT=NEIBPETOT
    ENDIF

    
    !C
    !C-- SEND
    
    DO neib= 1, NEIBPETOT
       istart= STACK_IMPORT(neib-1)
       inum  = STACK_IMPORT(neib  ) - istart
       
       DO k= istart+1, istart+inum
          WR(k)= X(NOD_IMPORT(k)+ZERO_ORIGIN)
       ENDDO

#ifdef DEBUG
       IF(istart+1>R) WRITE(*,*) "allocation error in send_recv",istart+1,R
#endif
       CALL MPI_ISEND (WR(istart+1), inum, MPI_DOUBLE_PRECISION,       &
            &                  NEIBPE(neib), 0, SOLVER_COMM,           &
            &                  req1(neib), ierr)
    ENDDO

    !C
    !C-- RECEIVE

    DO neib= 1, NEIBPETOT
       istart= STACK_EXPORT(neib-1)
       inum  = STACK_EXPORT(neib  ) - istart

#ifdef DEBUG
       IF(istart+1>S) WRITE(*,*) "allocation error in send_recv",istart+1,S
#endif

       CALL MPI_IRECV (WS(istart+1), inum, MPI_DOUBLE_PRECISION,       &
            &                  NEIBPE(neib), 0, SOLVER_COMM,           &
            &                  req2(neib), ierr)
       
    ENDDO

    IF(NEIBPETOT>0) CALL MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
    
    DO neib= 1, NEIBPETOT
       istart= STACK_EXPORT(neib-1)
       inum  = STACK_EXPORT(neib  ) - istart
       DO k= istart+1, istart+inum
          i = NOD_EXPORT(k)+ZERO_ORIGIN
          X(i)= X(i)+WS(k)
       ENDDO
    ENDDO

    IF(NEIBPETOT>0) CALL MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)
    
  END SUBROUTINE SOLVER_SR_REDUCE


  SUBROUTINE check_Temp_P(Temp_P, Temp_P_row_size, my_rank, SOLVER_COMM, LEVEL_NO, GIN_aggregate, &
       &     aggregate_table_size, aggregate_table_array, aggregate_number_in_table, &
       &     PE_list_size, local_aggre_size, RATE_OF_SPACE)
    USE data_structure_for_AMG
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    
    INTEGER(kind=kint), INTENT(in) :: LEVEL_NO, SOLVER_COMM, my_rank, RATE_OF_SPACE
    
    TYPE(row_node), INTENT(IN)      :: Temp_P(:,:)
    INTEGER(kind=kint ), INTENT(IN) :: Temp_P_row_size(:)
    
    INTEGER(kind=kint) :: N,NP,NEIBPETOT
    INTEGER(kind=kint), POINTER :: NEIBPE(:), STACK_EXPORT(:), STACK_IMPORT(:)
    INTEGER(kind=kint), POINTER :: NOD_EXPORT(:), NOD_IMPORT(:)

    INTEGER(kind=kint), ALLOCATABLE :: recv_buf_i(:)
    REAL(kind=kreal),  ALLOCATABLE :: recv_buf_v(:)

    INTEGER(kind=kint), ALLOCATABLE :: send_buf_i(:)
    REAL(kind=kreal), ALLOCATABLE :: send_buf_v(:)

    INTEGER(kind=kint) :: i,j,k,l,haba,istart,inum,neib,ierr,m,col,ownnodes,global_base
    REAL(kind=kreal) :: v1,v2

    INTEGER(kind=kint ), ALLOCATABLE :: sta1(:,:)
    INTEGER(kind=kint ), ALLOCATABLE :: sta2(:,:)
    INTEGER(kind=kint ), ALLOCATABLE :: req1(:  )
    INTEGER(kind=kint ), ALLOCATABLE :: req2(:  )  
    
    INTEGER(kind=kint), INTENT(in) :: GIN_aggregate(0:)
    INTEGER(kind=kint), INTENT(in) :: PE_list_size,local_aggre_size
    INTEGER(kind=kint), INTENT(in) :: aggregate_table_size
    !C two arguments: communication table,each PE
    INTEGER(kind=kint), DIMENSION(:,:), INTENT(inout) :: aggregate_table_array
    !C argument: the number of aggregates
    INTEGER(kind=kint), DIMENSION(:), INTENT(inout) :: aggregate_number_in_table

    INTEGER(kind=kint), ALLOCATABLE, DIMENSION(:) :: LG_vector(:)


    NP=HIERARCHICAL_DATA(LEVEL_NO-1) % NP
    N =HIERARCHICAL_DATA(LEVEL_NO-1) % N
    NEIBPETOT     =  HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPETOT
    NEIBPE        => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPE
    STACK_IMPORT  => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_IMPORT
    NOD_IMPORT    => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_IMPORT
    STACK_EXPORT  => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_EXPORT
    NOD_EXPORT    => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_EXPORT


    ALLOCATE (sta1(MPI_STATUS_SIZE,NEIBPETOT))
    ALLOCATE (sta2(MPI_STATUS_SIZE,NEIBPETOT))
    ALLOCATE (req1(NEIBPETOT))
    ALLOCATE (req2(NEIBPETOT))
    haba=INT(rate_of_space)
    ALLOCATE(recv_buf_i(haba*STACK_IMPORT(NEIBPETOT)))
    ALLOCATE(recv_buf_v(haba*STACK_IMPORT(NEIBPETOT)))
    ALLOCATE(send_buf_i(haba*STACK_EXPORT(NEIBPETOT)))
    ALLOCATE(send_buf_v(haba*STACK_EXPORT(NEIBPETOT)))


    ownnodes = GIN_aggregate(my_rank+1)-GIN_aggregate(my_rank)
    global_base = GIN_aggregate(my_rank)

    ALLOCATE(LG_vector(local_aggre_size))
    LG_vector=0
    DO i=1,ownnodes
       LG_vector(i)=i+global_base
    END DO
    
    DO i=1,PE_list_size
       DO j=1,aggregate_number_in_table(i)
          k=2*j
          l=aggregate_table_array(k,i)
          LG_vector(l)=aggregate_table_array(k-1,i)
       END DO
    END DO

    send_buf_i=-1
    send_buf_v=-1

    DO i=1,NEIBPETOT
       DO j=STACK_EXPORT(i-1)+1,STACK_EXPORT(i)
          k=NOD_EXPORT(j)

          IF(Temp_P_row_size(k)>haba) STOP "Temp_P_row_size in check_temp_P"

          DO l=1,Temp_P_row_size(k)
             send_buf_i(haba*(j-1)+l)=LG_vector( Temp_P(l,k) % column )
             send_buf_v(haba*(j-1)+l)=Temp_P(l,k) % value
          END DO
       END DO
    END DO

    DO neib=1,NEIBPETOT
       istart= stack_export(neib-1)
       inum  = stack_export(neib  ) - istart
       CALL MPI_ISEND(send_buf_i(haba*istart+1), inum*haba, LIS_MPI_INTEGER, &
            &         NEIBPE(neib),0, SOLVER_COMM, req1(neib), ierr )
    END DO
    DO neib=1,NEIBPETOT
       istart= stack_import(neib-1)
       inum  = stack_import(neib  ) - istart
       CALL MPI_IRECV(recv_buf_i(haba*istart+1), inum*haba, LIS_MPI_INTEGER, &
            &         NEIBPE(neib),0, SOLVER_COMM, req2(neib), ierr )
    END DO
    CALL MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)    
    CALL MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)


    DO neib=1,NEIBPETOT
       istart= stack_export(neib-1)
       inum  = stack_export(neib  ) - istart
       CALL MPI_ISEND(send_buf_v(haba*istart+1), inum*haba, MPI_DOUBLE_PRECISION, &
            &         NEIBPE(neib),0, SOLVER_COMM, req1(neib), ierr )
    END DO
    DO neib=1,NEIBPETOT
       istart= stack_import(neib-1)
       inum  = stack_import(neib  ) - istart
       CALL MPI_IRECV(recv_buf_v(haba*istart+1), inum*haba, MPI_DOUBLE_PRECISION, &
            &         NEIBPE(neib),0, SOLVER_COMM, req2(neib), ierr )
    END DO
    CALL MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)    
    CALL MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)

    DEALLOCATE(send_buf_i)
    DEALLOCATE(send_buf_v)


    DO i=1,NEIBPETOT
       DO j=STACK_IMPORT(i-1)+1,STACK_IMPORT(i)
          k=NOD_IMPORT(j)+ZERO_ORIGIN
          DO m=1,haba
             col = recv_buf_i(haba*(j-1)+m)
             IF(col>0) THEN
                DO l=1,Temp_P_row_size(k)
                   IF(LG_vector(Temp_P(l,k) % column) == col) THEN
                      v1=Temp_P(l,k)%value;v2=recv_buf_v(haba*(j-1)+m)
                      IF(dabs(v1-v2)>EPS) THEN
                         WRITE(*,*) my_rank,":",k,Temp_P(l,k) % column,v1,v2
                      END IF
                      EXIT
                   END IF
                END DO
                IF(l>Temp_P_row_size(k)) WRITE(*,*) "nasi.in Temp_P:",k,col,recv_buf_v(haba*(j-1)+m)
             END IF
          END DO
       END DO
    END DO

    DEALLOCATE(recv_buf_i)
    DEALLOCATE(recv_buf_v)

    DEALLOCATE (sta1)
    DEALLOCATE (sta2)
    DEALLOCATE (req1)
    DEALLOCATE (req2)


  END SUBROUTINE check_Temp_P
    



!C  side effect: Temp_P,Temp_P_row_size,global_local_hash_table
!C              aggregate_table_array,     
!C This routine exchanges the information of the aggregates on EXPORT table, IMPORT table,
!C and Repeated Table. In Smoothe_aggregate(), elements of aggregates must be only on these
!C three tables. In order to guaranty that condition, PNI and NI choose the symmetric part of 
!C the problem matrix.
  SUBROUTINE solver_SR_aggregates &
       &    ( PE_num_size, PE_nums, nfepe_size, nodes_for_each_PE, in_nodes_for_each_PE, SOLVER_COMM,  &
       &      temp_P_1,temp_P_2,Temp_P, Temp_P_row_size, my_rank, GIN_aggregate, NPROCS,                     &
       &      hash_table_size, global_local_hash_table, aggregate_table_size,              &
       &      aggregate_table_array, aggregate_number_in_table, local_aggre_size, PE_list, &
       &      PE_list_size, LEVEL_NO) 
    USE hash_mod
    USE data_structure_for_AMG
    IMPLICIT NONE
    INCLUDE  'mpif.h'

    INTEGER(kind=kint ), INTENT(IN) :: LEVEL_NO,NPROCS, nfepe_size 
    INTEGER(kind=kint ), INTENT(IN) :: temp_P_1,temp_P_2
    
    INTEGER(kind=kint)              :: NEIBPETOT
    INTEGER(kind=kint),POINTER      :: NEIBPE(:),STACK_IMPORT(:),STACK_EXPORT(:)
    INTEGER(kind=kint),POINTER      :: NOD_IMPORT(:),NOD_EXPORT(:)
    
    INTEGER(kind=kint ), INTENT(IN) :: PE_num_size
!!$    INTEGER(kind=kint ), INTENT(IN) :: PE_nums(PE_num_size)
    INTEGER(kind=kint ), INTENT(IN) :: PE_nums(1:PE_num_size)
    INTEGER(kind=kint ), INTENT(IN) :: nodes_for_each_PE(1:nfepe_size)
!!$    INTEGER(kind=kint ), INTENT(IN) :: in_nodes_for_each_PE(0:PE_num_size)
    INTEGER(kind=kint ), INTENT(IN) :: in_nodes_for_each_PE(0:PE_num_size)
    TYPE(row_node), INTENT(INOUT)      :: Temp_P(1:temp_P_1,1:temp_P_2)
    INTEGER(kind=kint ), INTENT(INOUT) :: Temp_P_row_size(1:temp_P_2)

    INTEGER(kind=kint ), INTENT(inout) :: PE_list_size
    INTEGER(kind=kint ), INTENT(inout) :: PE_list(1:NPROCS,1:3)
    

    INTEGER(kind=kint ), INTENT(in) :: GIN_aggregate(0:NPROCS)
    
    INTEGER(kind=kint ), INTENT(in) :: SOLVER_COMM
    INTEGER(kind=kint ), INTENT(in) :: my_rank
    INTEGER(kind=kint ), ALLOCATABLE :: req1(:),req2(:)
    INTEGER(kind=kint ), ALLOCATABLE :: sta1(:,:),sta2(:,:)
    INTEGER(kind=kint) :: ierr

    INTEGER(kind=kint ), ALLOCATABLE :: in_send_buf_for_PEs(:),in_send_buf_for_PEs_dbl(:)

    INTEGER(kind=kint ), INTENT(out) :: local_aggre_size

    INTEGER(kind=kint ) :: PE_list_alloc_size,istart,inum

    INTEGER(kind=kint) :: Neib_index,Pe_num_index,pack_index,pack_index_dbl,dest,  &
         &                src,message_size
    INTEGER(kind=kint) :: i,j,k,l,l1,l2,m,n,g,GLOBAL,LOCAL,hash_no,index,h,own_aggre_size,pls

    INTEGER(kind=kint), ALLOCATABLE :: LG_vector(:)

    INTEGER(kind=kint) :: recv_buf_size,recv_buf_size_dbl
    INTEGER(kind=kint) :: send_buf_size,wait_count,send_buf_size_dbl
    INTEGER(kind=kint) :: position,position_dbl
    PARAMETER(GLOBAL=1,LOCAL=2)    
    INTEGER(kind=kint ),ALLOCATABLE :: send_buf    (:),recv_buf    (:)
    REAL   (kind=kreal),ALLOCATABLE  :: send_buf_dbl(:),recv_buf_dbl(:)


    INTEGER(kind=kint),INTENT(in):: HASH_TABLE_SIZE
    !C DIMENSION(hash_table_size,2)
    INTEGER(kind=kint), INTENT(out) :: global_local_hash_table(1:HASH_TABLE_SIZE, 1:2)

    INTEGER(kind=kint),INTENT(in) :: aggregate_table_size
!!$    !C two arguments: communication table,each PE
!!$    INTEGER(kind=kint), INTENT(inout) :: aggregate_table_array(1:aggregate_table_size*2, 1:NPROCS)
!!$    !C argument: the number of aggregates
!!$    INTEGER(kind=kint), INTENT(inout) :: aggregate_number_in_table(1:NPROCS)
    !C two arguments: communication table,each PE
    INTEGER(kind=kint), INTENT(inout) :: aggregate_table_array(:, :)
    !C argument: the number of aggregates
    INTEGER(kind=kint), INTENT(inout) :: aggregate_number_in_table(:)

    INTEGER(kind=kint), DIMENSION(-1:2*NPROCS)::s_ary, r_ary

    INTEGER(kind=kint) :: NP
#ifdef DEBUG
    INTEGER(kind=kint),ALLOCATABLE :: deb1(:),deb2(:)
#endif


    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP
    
    NEIBPETOT     =  HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPETOT
    NEIBPE        => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPE
    STACK_IMPORT  => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_IMPORT
    NOD_IMPORT    => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_IMPORT
    STACK_EXPORT  => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_EXPORT
    NOD_EXPORT    => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_EXPORT


    global_local_hash_table = 0


    !C-- count up the external aggregates and set local_aggre_size
    j=0
    DO i=1,PE_list_size
       j=j+aggregate_number_in_table(i)
    END DO

    ALLOCATE(LG_vector(j))
    own_aggre_size   = GIN_aggregate(my_rank+1)-GIN_aggregate(my_rank)
    local_aggre_size = own_aggre_size+j
    !C-- end:count up the external aggregates and set local_aggre_size


    !C-- assign aggregate_table_array into hash_table, and create LG_vector(:)    
    LG_vector=0

    DO i = 1, PE_list_size
       j = aggregate_number_in_table(i)
       DO k = 1, j
          g = aggregate_table_array(2 * k - 1, i)
          l = aggregate_table_array(2 * k, i)
          LG_vector(l - own_aggre_size) = g
          
          
          CALL hash_zero(HASH_TABLE_SIZE, global_local_hash_table(:, GLOBAL), g,  &
               &         h, "HASH_TABLE_SIZE")
          global_local_hash_table(h, GLOBAL) = g
          global_local_hash_table(h, LOCAL ) = l
       END DO
    END DO
    !C-- end:assign aggregate_table_array into hash_table    

    
    !C-- count up  Send bufsize and make in_send_buf_for_PEs
    send_buf_size    = 0
    send_buf_size_dbl= 0
    ALLOCATE(in_send_buf_for_PEs    (0:PE_list_size))
    ALLOCATE(in_send_buf_for_PEs_dbl(0:PE_list_size))
    in_send_buf_for_PEs    (0)=0
    in_send_buf_for_PEs_dbl(0)=0
    DO i=1,PE_list_size
       dest         = PE_list(i,1)       
       Neib_index   = PE_list(i,2)
       PE_num_index = PE_list(i,3)
       IF(Neib_index > 0) THEN
          j=STACK_IMPORT(NEIBPETOT)
          CALL count_NOD(STACK_IMPORT(0:NEIBPETOT),NOD_IMPORT(1:j),               &
               &        Temp_P_row_size(1:NP),SOLVER_COMM,                  &
               &        send_buf_size,send_buf_size_dbl,Neib_index,               &
               &        j,NP,NEIBPETOT)
          j=STACK_EXPORT(NEIBPETOT)
          CALL count_NOD(STACK_EXPORT(0:NEIBPETOT),NOD_EXPORT(1:j),               &
               &        Temp_P_row_size(1:NP),SOLVER_COMM,                  &
               &        send_buf_size,send_buf_size_dbl,Neib_index,               &
               &        j,NP,NEIBPETOT)

       END IF
       IF(PE_num_index > 0) THEN
          !C   PE_num_size,PE_nums,in_nodes_for_each_PE,nodes_for_each_PE
          j = in_nodes_for_each_PE(PE_num_size)
          CALL count_NOD(in_nodes_for_each_PE(0:PE_num_size), &
               &        nodes_for_each_PE(1:j),Temp_P_row_size(1:NP),&
               &        SOLVER_COMM,send_buf_size,send_buf_size_dbl,PE_num_index,    &
               &        j, NP, PE_num_size)

       END IF
       in_send_buf_for_PEs    (i) = send_buf_size
       in_send_buf_for_PEs_dbl(i) = send_buf_size_dbl

    END DO

    ALLOCATE(send_buf    (send_buf_size+1    ))
    ALLOCATE(send_buf_dbl(send_buf_size_dbl+1))
    !C-- end: count up  Send bufsize and make in_send_buf_for_PEs
    

    !C-- sendrecv  in_send_buf_for_PEs
    ALLOCATE(req1(PE_list_size))
    ALLOCATE(req2(PE_list_size))
    ALLOCATE(sta2(MPI_STATUS_SIZE,PE_list_size))
    ALLOCATE(sta1(MPI_STATUS_SIZE,PE_list_size))

    r_ary=0
    DO i=1,PE_list_size
       src=PE_list(i,1)
       CALL MPI_IRECV(r_ary(i*2-1),2,LIS_MPI_INTEGER,src ,0,SOLVER_COMM,&
            &         req2(i),ierr)
    END DO

    s_ary=0
    DO i=1,PE_list_size
       s_ary(i*2-1)=in_send_buf_for_PEs    (i)-in_send_buf_for_PEs    (i-1)
       s_ary(i*2)  =in_send_buf_for_PEs_dbl(i)-in_send_buf_for_PEs_dbl(i-1)
    END DO

    DO i=1,PE_list_size
       dest=PE_list(i,1)
       CALL MPI_ISEND(s_ary(2*i-1),2,LIS_MPI_INTEGER,dest,0,SOLVER_COMM,&
            &         req1(i),ierr)
    END DO
    !C-- end: sendrecv in_send_buf_for_PEs


    !C-- Create Send buffer for neighboring PEs
    !C PE_list,in_send_buf_for_PEs has data from 1 to PE_list_size.
    !C position is the next space for packing.
    DO i=1,PE_list_size
       dest=PE_list(i,1)
       Neib_index=PE_list(i,2)
       PE_num_index=PE_list(i,3)
       pack_index    = in_send_buf_for_PEs    (i-1)
       pack_index_dbl= in_send_buf_for_PEs_dbl(i-1)

       CALL create_buf(dest, Neib_index, PE_num_index, NEIBPETOT, NEIBPE,      &
            &  STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, PE_num_size,&
            &  PE_nums, nodes_for_each_PE, in_nodes_for_each_PE, Temp_P,       &
            &  Temp_P_row_size, send_buf, send_buf_size, send_buf_dbl,         &
            &  send_buf_size_dbl, my_rank, GIN_aggregate(my_rank),             &
            &  SOLVER_COMM, pack_index, pack_index_dbl, own_aggre_size,        &
            &  LG_vector)

    END DO
    DEALLOCATE(LG_vector)
    !C-- end:Create Send buffer for neighboring PEs

    !C-- wait sendrecv of in_send_buf_for_PEs
    !C   r_ary is index of recv_buf(:)
    CALL MPI_WAITALL (PE_list_size, req2, sta1, ierr)
    CALL MPI_WAITALL (PE_list_size, req1, sta2, ierr)

    r_ary(0)=0;r_ary(-1)=0
    
    DO i=1,PE_list_size
       r_ary(2*i-1)=r_ary(2*i-3) +r_ary(2*i-1)
       r_ary(2*i)  =r_ary(2*i-2) +r_ary(2*i)
    END DO
    !C-- end:  wait sendrecv of in_send_buf_for_PEs
    
    !C-- communication the aggregates on border table of PE_list and input.
    !C   in this part PE_list(:) means neighboring pes within two path.
    pls=PE_list_size
    recv_buf_size    =r_ary(pls*2-1)
    recv_buf_size_dbl=r_ary(pls*2)
    ALLOCATE(recv_buf(recv_buf_size+1))
    ALLOCATE(recv_buf_dbl(recv_buf_size_dbl+1))

    DO i=1,pls
       src=PE_list(i,1)
       istart = r_ary(2*i-3)
       inum   = r_ary(2*i-1)-istart

       CALL MPI_IRECV(recv_buf(istart+1),inum,LIS_MPI_INTEGER,src,30,SOLVER_COMM, &
            &         req2(i),ierr)
    END DO
    

    
    DO i=1,pls
       dest=PE_list(i,1)
       istart = in_send_buf_for_PEs(i-1)
       inum   = in_send_buf_for_PEs(i)-istart
       CALL MPI_ISEND(send_buf(istart+1), inum, LIS_MPI_INTEGER,dest,30,SOLVER_COMM,&
            &         req1(i),ierr)
    END DO

    CALL MPI_WAITALL(pls, req2, sta2, ierr)
    CALL MPI_WAITALL(pls, req1, sta1, ierr)
    
    DO i=1,pls
       src=PE_list(i,1)
       istart = r_ary(2*i-2)
       inum   = r_ary(2*i)-istart

       CALL MPI_IRECV(recv_buf_dbl(istart+1),inum,MPI_DOUBLE_PRECISION,src,31,SOLVER_COMM, &
            &         req2(i),ierr)
    END DO

    DO i=1,pls
       dest=PE_list(i,1)
       istart = in_send_buf_for_PEs_dbl(i-1)
       inum   = in_send_buf_for_PEs_dbl(i)-istart
       CALL MPI_ISEND(send_buf_dbl(istart+1),inum, MPI_DOUBLE_PRECISION,dest,31,SOLVER_COMM,&
            &         req1(i),ierr)
    END DO

    CALL MPI_WAITALL(pls, req2, sta2, ierr)
    CALL MPI_WAITALL(pls, req1, sta1, ierr)


    !C UNPACK should be done in accordance with  PE_list()
    !C PE_list_size is changed
    DO i=1,pls
       src=PE_list(i,1)
       Neib_index=PE_list(i,2)
       PE_num_index=PE_list(i,3)
       j=r_ary(2*i-3)+1; k=r_ary(2*i-1)
       m=r_ary(2*i-2)+1; n=r_ary(2*i)
       l1=k-j+1;l2=n-m+1

       position=0
       position_dbl=0


       IF(Neib_index /= 0) THEN
          !C receive NOD_EXPORT and NOD_IMPORT
          CALL input_buf(Neib_index, Temp_P, Temp_P_row_size, recv_buf(j:k), recv_buf_dbl(m:n), &
               & hash_table_size, global_local_hash_table, position, position_dbl,              &
               & STACK_EXPORT, NOD_EXPORT, aggregate_table_size, aggregate_table_array,         &
               & aggregate_number_in_table, local_aggre_size, SOLVER_COMM, my_rank, l1, l2,     &
               & GIN_aggregate, PE_list_size, PE_list, NPROCS)


          CALL input_buf(Neib_index, Temp_P, Temp_P_row_size, recv_buf(j:k), recv_buf_dbl(m:n), &
               & HASH_TABLE_SIZE, global_local_hash_table, position, position_dbl,              &
               & STACK_IMPORT, NOD_IMPORT, aggregate_table_size, aggregate_table_array,         &
               & aggregate_number_in_table, local_aggre_size, SOLVER_COMM, my_rank, l1, l2,     &
               & GIN_aggregate, PE_list_size, PE_list, NPROCS)

       END IF
       IF(PE_num_index /= 0 ) THEN
          !C receive (in_nodes_for_each_PE, nodes_for_each_PE)
          CALL input_buf(PE_num_index, Temp_P, Temp_P_row_size, recv_buf(j:k),                &
               & recv_buf_dbl(m:n), hash_table_size, global_local_hash_table, position,       &
               & position_dbl, in_nodes_for_each_PE, nodes_for_each_PE, aggregate_table_size, &
               & aggregate_table_array, aggregate_number_in_table, local_aggre_size,          &
               & SOLVER_COMM, my_rank, l1, l2, GIN_aggregate, PE_list_size, PE_list, NPROCS)

       END IF
    END DO



    DEALLOCATE(send_buf)     
    DEALLOCATE(send_buf_dbl)     
    DEALLOCATE(in_send_buf_for_PEs)
    DEALLOCATE(in_send_buf_for_PEs_dbl)
    DEALLOCATE(recv_buf)
    DEALLOCATE(recv_buf_dbl)
    DEALLOCATE(req1)
    DEALLOCATE(req2)
    DEALLOCATE(sta1)
    DEALLOCATE(sta2)
    
    !C-- communication the aggregates on border table of PE_list and input.



#ifdef DEBUG
    DO i = 1, PE_list_size
       DO j = 1, aggregate_number_in_table(i)
          k = 2 * j
          l = aggregate_table_array(k - 1, i)
          dest  = PE_list(i, 1)
          IF(l <= GIN_aggregate(dest) .OR. l > GIN_aggregate(dest + 1)) THEN
             WRITE(*,*)"error in sr_aggregates:", l, GIN_aggregate(dest), GIN_aggregate(dest + 1), dest
             STOP "aggregate_table_array is corrupted : sr_aggregate "
          END IF
       END DO
    END DO
#endif


  CONTAINS

    SUBROUTINE count_NOD(STACK_PORT, NOD_PORT, Temp_P_row_size, solver_comm, &
         &              send_buf_size, send_buf_size_dbl, Neib_index, NPT_sz, &
         &              TP_sz, NEIBPETOT)
      
      USE data_structure_for_AMG
      IMPLICIT NONE
      INCLUDE  'mpif.h'
      INTEGER(kind=kint), INTENT(in) :: NEIBPETOT
      INTEGER(kind=kint), INTENT(in)    :: Neib_index,NPT_sz,TP_sz
      INTEGER(kind=kint), INTENT(in) :: STACK_PORT(0:NEIBPETOT)
      INTEGER(kind=kint), INTENT(in)  :: NOD_PORT(1:NPT_sz),Temp_P_row_size(TP_sz)
      INTEGER(kind=kint), INTENT(inout) :: send_buf_size,send_buf_size_dbl
      INTEGER(kind=kint), INTENT(in)                :: SOLVER_COMM

      INTEGER(kind=kint) :: ierr,in_buf_size,membersize
      INTEGER(kind=kint) :: i,node_index_base,buf_size,node

      !C in_buf(0:in_buf_size)
      !C aggregate_num_buf(buf_size)
      !C value_buf(buf_size)
      node_index_base = STACK_PORT(Neib_index - 1)
      in_buf_size = STACK_PORT(Neib_index) - node_index_base



      
      buf_size = 0
      DO i = 1, in_buf_size
         node = NOD_PORT(node_index_base + i)+ZERO_ORIGIN
         buf_size = buf_size + Temp_P_row_size(node)
      END DO

      
      !C first LIS_MPI_INTEGER is the size of in_buf
      !C zero:in_buf(0) is not sended.
      
      send_buf_size = send_buf_size + 1 + 1 + in_buf_size + buf_size
      send_buf_size_dbl = send_buf_size_dbl + buf_size

    END SUBROUTINE count_NOD
    
    

!C side effect: index,Temp_P,Temp_P_row_size,global_local_hash_table,
!C              position,aggregate_table_size,aggregate_table_each_PE,
!C              aggregate_number,local_aggre_size

    SUBROUTINE input_buf (index, Temp_P, Temp_P_row_size, recv_buf, recv_buf_dbl,    &
         &                hash_table_size, global_local_hash_table, position,        &
         &                position_dbl, STACK_PORT, NOD_PORT, aggregate_table_size,  &
         &                aggregate_table_array, aggregate_number_in_table,          &
         &                local_aggre_size, SOLVER_COMM, my_rank, l1, l2,            &
         &                GIN_aggregate, PE_list_size, PE_list, NPROCS)

      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER(kind=kint), INTENT(in) :: my_rank, NPROCS, l1, l2
      INTEGER(kind=kint), INTENT(inout) :: PE_list_size
      INTEGER(kind=kint), INTENT(inout) :: PE_list(:,:)

      !C max number of local aggregate number
      INTEGER(kind=kint), INTENT(inout) :: local_aggre_size
      !C comunication table is one dimensional array of aggregate_table_size*2
      INTEGER(kind=kint), INTENT(in) :: aggregate_table_size
      !C communication table: global,local
!!$      INTEGER(kind=kint),  INTENT(inout) :: aggregate_table_array(aggregate_table_size * 2, NPROCS)
      INTEGER(kind=kint),  INTENT(inout) :: aggregate_table_array(:, :)
      !C argument: the number of aggregates
!!$      INTEGER(kind=kint), INTENT(inout) :: aggregate_number_in_table(NPROCS)
      INTEGER(kind=kint), INTENT(inout) :: aggregate_number_in_table(:)

      INTEGER(kind=kint), INTENT(inout) :: position,position_dbl    
      INTEGER(kind=kint), DIMENSION(0:),INTENT(in) :: STACK_PORT
      INTEGER(kind=kint), DIMENSION(:), INTENT(in) :: NOD_PORT
      INTEGER(kind=kint), INTENT(in) :: index

      INTEGER(kind=kint), INTENT(in) :: GIN_aggregate(0:)

      INTEGER(kind=kint) :: GLOBAL,LOCAL
      PARAMETER(GLOBAL=1,LOCAL=2)

      TYPE(row_node), DIMENSION(:,:), INTENT(inout) :: Temp_P
      INTEGER(kind=kint ), DIMENSION(:) ,INTENT(inout):: Temp_P_row_size
      INTEGER(kind=kint ), TARGET,INTENT(in) :: recv_buf(l1)
!!$      INTEGER(kind=kint ), target,intent(in) :: recv_buf(:)
!!$      INTEGER(kind=kint ), pointer :: recv_buf(:)
      REAL   (kind=kreal), TARGET,INTENT(in) :: recv_buf_dbl(l2)
!!$      REAL   (kind=kreal), target :: recv_buf_dbl(:)
!!$      REAL   (kind=kreal), pointer :: recv_buf_dbl(:)
      INTEGER(kind=kint) , INTENT(in) :: hash_table_size
      INTEGER(kind=kint) , DIMENSION(:,:),INTENT(inout) :: global_local_hash_table
      INTEGER(kind=kint) ,INTENT(in) :: SOLVER_COMM

      INTEGER(kind=kint ) :: ierr

      LOGICAL :: hash_flag,Temp_flag
      INTEGER(kind=kint) :: i,j,k,l,m,Neib_index,PE_num_index,buf_size,geta_NOD_PORT
      INTEGER(kind=kint) :: ins,ine,g_ag_no,local_ag_no,hash_no,row,base,bound,pe
      REAL(kind=kreal) :: value
      INTEGER(kind=kint) :: buf_node_size
      INTEGER(kind=kint), POINTER :: buf_index(:)
      INTEGER(kind=kint), POINTER :: buf_of_global_aggregate_number(:)
      REAL(kind=kreal),   POINTER :: buf_of_value(:)

      base =GIN_aggregate(my_rank)
      bound=GIN_aggregate(my_rank+1)

      !C buf_node_size
      position=position+1
      buf_node_size=recv_buf(position)
      

      !C buf_index
      buf_size=0
      IF(buf_node_size>0) THEN
         !C  buf_index starts from 0.
         buf_index => recv_buf((position+1):(position+buf_node_size+1))
         buf_size  =  buf_index(buf_node_size+1)
      END IF
      position=position+buf_node_size+1

      !C buf_of_global_aggregate_number,and buf_of_value
      IF(buf_size>0) THEN

         buf_of_global_aggregate_number=>recv_buf((position+1):(position+buf_size))
         buf_of_value                  =>recv_buf_dbl((position_dbl+1):(position_dbl+buf_size))
         position    = position     + buf_size
         position_dbl= position_dbl + buf_size
         
      END IF
      
      geta_NOD_PORT=STACK_PORT(index-1)
      DO i=1,buf_node_size
         !C buf_index starts from 1
         ins=buf_index(i+1-1)+1
         ine=buf_index(i+1)


         DO j=ins,ine
            g_ag_no = buf_of_global_aggregate_number(j)
            value = buf_of_value(j)

            !C g_ag_no in own aggregates
            !C g_ag_no not in own aggregates
            !C   g_ag_no not in hash_table
            !C   g_ag_no in hash_table

            IF(base < g_ag_no .AND. g_ag_no <= bound)  THEN
               local_ag_no = g_ag_no - base
            ELSE
               !C-- insert g_ag_no to hash_table
               !C hash_table should be initialized to (0,0).
               CALL hash(HASH_TABLE_SIZE, global_local_hash_table(:, GLOBAL), g_ag_no, k, &
                    & "HASH_TABLE_SIZE")               
               local_ag_no = global_local_hash_table(k,LOCAL)
            END IF

            IF(local_ag_no == 0) THEN
               !C local_ag_no not in hashtable
               !C record the new aggregate
               local_aggre_size = local_aggre_size+1
               global_local_hash_table(k, GLOBAL) = g_ag_no
               global_local_hash_table(k, LOCAL)  = local_aggre_size
               local_ag_no = local_aggre_size
#ifdef DEBUG
               IF(g_ag_no < 1 .OR. g_ag_no > GIN_aggregate(NPROCS)) THEN
                  WRITE(*,*) "unexpected error!! in input_buf()", g_ag_no,GIN_aggregate(NPROCS)
                  STOP
               END IF
#endif               
               !C -- begin: insert aggregate number into aggregate_table_array
               DO k = 1, PE_list_size
                  l = PE_list(k, 1)
                  IF(GIN_aggregate(l) < g_ag_no .AND. g_ag_no <= GIN_aggregate(l + 1)) EXIT
               END DO
               IF(k > PE_list_size)THEN 
                  DO pe = 0, NPROCS - 1
                     IF(GIN_aggregate(pe) < g_ag_no .AND. g_ag_no <= GIN_aggregate(pe + 1)) EXIT
                  END DO
#ifdef DEBUG                
                  IF(pe >= NPROCS) STOP "unexpected error!! in input_buf()"
#endif                  
                  k = PE_list_size + 1
                  PE_list_size = k
                  PE_list(k,1) = pe
                  PE_list(k,2) = 0
                  PE_list(k,3) = 0
                  aggregate_number_in_table(k) = 0
               END IF
               
               aggregate_number_in_table(k) = aggregate_number_in_table(k) + 1
               m = aggregate_number_in_table(k)
               
#ifdef ALLOCATION_CHECK
               IF(aggregate_table_size < m) THEN
                  WRITE(*,*) "enlarge aggregate_table_size"
                  STOP 'allocation error in input_buf():'
               END IF
#endif                   

               aggregate_table_array(m * 2 - 1, k) = g_ag_no
               aggregate_table_array(m * 2, k) = local_ag_no
               !C -- end: insert aggregat number into aggregate_table_array

               !C-- insert node to Temp_P
               !C row is node number, column is aggregate number sent by other PE
               row = NOD_PORT(i + geta_NOD_PORT)+ZERO_ORIGIN
               
               k = Temp_P_row_size(row) + 1
               Temp_P_row_size(row) = k
               
               Temp_P(k,row) % column = local_ag_no
               Temp_P(k,row) % value  = value
               !C-- end: insert node to Temp_P
            ELSE
               !C-- insert node to Temp_P : insertion by local_ag_no
               row = NOD_PORT(i + geta_NOD_PORT)+ZERO_ORIGIN
               k   = Temp_P_row_size(row)

               l=1
               DO WHILE( l <= k .AND. Temp_P(l,row) % column < local_ag_no) 
                  l = l+1
               END DO
               
               IF( l <= k .AND. Temp_P(l,row) % column == local_ag_no) THEN
                  Temp_P(l,row) % value = Temp_P(l,row) % value + value
               ELSE
                  !C insertion by column
                  Temp_P_row_size(row) = k+1
                  DO m = k,l,-1
                     Temp_P(m+1,row) % column = Temp_P(m,row) % column
                     Temp_P(m+1,row) % value  = Temp_P(m,row) % value 
                  END DO
                  Temp_P(l,row) % column = local_ag_no
                  Temp_P(l,row) % value  = value
               END IF
               !C-- end: insert node to Temp_P
            END IF
            !C-- end:insert g_ag_no to hash_table
         END DO
      END DO
    END SUBROUTINE input_buf


    SUBROUTINE create_buf(dest,Neib_index,PE_num_index,NEIBPETOT, NEIBPE,   &
         &  STACK_IMPORT, NOD_IMPORT, STACK_EXPORT, NOD_EXPORT, PE_num_size,&
         &  PE_nums,nodes_for_each_PE,in_nodes_for_each_PE,Temp_P,          &
         &  Temp_P_row_size,send_buf,send_buf_size,send_buf_dbl,            &
         &  send_buf_size_dbl,my_rank,GIN_rank,SOLVER_COMM,                 &
         &  pack_index,pack_index_dbl,own_aggre_size,LG_vector)

      USE data_structure_for_AMG
      IMPLICIT NONE
      INCLUDE  'mpif.h'

      INTEGER(kind=kint), INTENT(in) :: dest, Neib_index, PE_num_index, NEIBPETOT
      INTEGER(kind=kint), INTENT(in) :: my_rank,pack_index,pack_index_dbl,own_aggre_size
      INTEGER(kind=kint), INTENT(in) :: SOLVER_COMM
      INTEGER(kind=kint), INTENT(in) :: send_buf_size,send_buf_size_dbl
      INTEGER(kind=kint)             :: ierr

      INTEGER(kind=kint ), INTENT(out):: send_buf(:)
      REAL   (kind=kreal), INTENT(out):: send_buf_dbl(:)

      INTEGER(kind=kint ), DIMENSION(:), INTENT(in) :: NEIBPE
      INTEGER(kind=kint ), DIMENSION(0:),INTENT(in) :: STACK_IMPORT,STACK_EXPORT
      INTEGER(kind=kint ), DIMENSION(:), INTENT(in) :: NOD_IMPORT,NOD_EXPORT
      INTEGER(kind=kint ), DIMENSION(:), INTENT(in) :: LG_vector

      INTEGER(kind=kint ), INTENT(in) :: PE_num_size
      INTEGER(kind=kint ), DIMENSION(:), INTENT(in) :: PE_nums
      INTEGER(kind=kint ), DIMENSION(:), INTENT(in) :: nodes_for_each_PE
      INTEGER(kind=kint ), DIMENSION(0:), INTENT(in) :: in_nodes_for_each_PE
      TYPE(row_node), DIMENSION(:,:) :: Temp_P
      INTEGER(kind=kint ), DIMENSION(:), INTENT(in)  :: Temp_P_row_size
      INTEGER(kind=kint ), INTENT(in)    :: GIN_rank
      INTEGER(kind=kint )                :: position, position_dbl

      !C Neib_index != 0 
      position=0;position_dbl=0
      IF(Neib_index > 0) THEN
         CALL pack_NOD(STACK_IMPORT,NOD_IMPORT,Temp_P_row_size,Temp_P,SOLVER_COMM,&
              &        GIN_rank,send_buf,send_buf_size,send_buf_dbl,              &
              &        send_buf_size_dbl,position,Neib_index,pack_index,          &
              &        position_dbl,pack_index_dbl,own_aggre_size,LG_vector)

         CALL pack_NOD(STACK_EXPORT,NOD_EXPORT,Temp_P_row_size,Temp_P,SOLVER_COMM,&
              &        GIN_rank,send_buf,send_buf_size,send_buf_dbl,              &
              &        send_buf_size_dbl,position,Neib_index,pack_index,          &
              &        position_dbl,pack_index_dbl,own_aggre_size,LG_vector)
      END IF
      IF(PE_num_index > 0) THEN
         !C   PE_num_size,PE_nums,in_nodes_for_each_PE,nodes_for_each_PE
         CALL pack_NOD(in_nodes_for_each_PE,nodes_for_each_PE,Temp_P_row_size,Temp_P,&
              &        SOLVER_COMM,GIN_rank,send_buf,send_buf_size,send_buf_dbl,     &
              &        send_buf_size_dbl, position,PE_num_index, pack_index,         &
              &        position_dbl,pack_index_dbl,own_aggre_size,LG_vector)
      END IF
    END SUBROUTINE create_buf


    SUBROUTINE pack_NOD(STACK_PORT,NOD_PORT,Temp_P_row_size,Temp_P,solver_comm, &
         &              GIN_rank,send_buf,send_buf_size,send_buf_dbl,           &
         &              send_buf_size_dbl,position,Neib_index,pack_index,       &
         &              position_dbl,pack_index_dbl,own_aggre_size,LG_vector)

      USE data_structure_for_AMG
      IMPLICIT NONE
      INCLUDE  'mpif.h'
      TYPE(row_node)    :: Temp_P(:,:)
      INTEGER(kind=kint ), INTENT(in)    :: own_aggre_size
      INTEGER(kind=kint ), INTENT(in)    :: STACK_PORT(0:)
      INTEGER(kind=kint ), INTENT(in)    :: NOD_PORT(:),Temp_P_row_size(:)
      INTEGER(kind=kint ), INTENT(in)    :: LG_vector(:)

      INTEGER(kind=kint ), INTENT(in)    :: send_buf_size,send_buf_size_dbl,GIN_rank
      INTEGER(kind=kint ), INTENT(in)    :: Neib_index,pack_index,pack_index_dbl
      INTEGER(kind=kint ), INTENT(inout) :: position, position_dbl
      INTEGER(kind=kint ), INTENT(inout) :: send_buf(:)
      REAL   (kind=kreal), INTENT(inout) :: send_buf_dbl(:)
      INTEGER(kind=kint ), INTENT(in)    :: SOLVER_COMM

      INTEGER(kind=kint ) :: ierr,in_buf_size,tprsz
      INTEGER(kind=kint ) :: i,j,k,l,node_index_base,node,count,pack_size
      
      node_index_base = STACK_PORT(Neib_index-1)
      in_buf_size = STACK_PORT(Neib_index)-node_index_base

      position =position+1
      send_buf(pack_index+position)=in_buf_size

      position =position+1
      send_buf(pack_index+position)=0

      DO i=1,in_buf_size
         node = NOD_PORT(node_index_base+i)+ZERO_ORIGIN
         send_buf(pack_index+position+i) = send_buf(pack_index+position+i-1) + &
              &                            Temp_P_row_size(node)
      END DO
      position = position+in_buf_size


      count = 0
      DO i=1,in_buf_size
         node = NOD_PORT(node_index_base+i)+ZERO_ORIGIN
         TPrsz=Temp_P_row_size(node)
         DO j=1,TPrsz
            send_buf_dbl( pack_index_dbl+position_dbl+j ) = Temp_P(j,node) % value

            k = Temp_P(j,node) % column
            IF(k<=own_aggre_size) THEN
               send_buf(pack_index+position+j) = k + GIN_rank
            ELSE 
               send_buf(pack_index+position+j) = LG_vector(k-own_aggre_size) 
#ifdef DEBUG
               IF(LG_vector(k-own_aggre_size)==0) THEN
                  WRITE(*,*) "unexpected error in LG_vector in pack_NOD()"
                  STOP
               END IF
#endif
            END IF
         END DO
         position     = position     + tprsz
         position_dbl = position_dbl + tprsz
      END DO
    END SUBROUTINE pack_NOD

  END SUBROUTINE solver_SR_aggregates



END MODULE solver_Gnumbering

