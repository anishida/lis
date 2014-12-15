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
!C   * MODULE aggregate_mod
!C     CONTAINS
!C   * SUBROUTINE indpdt_agrgt
!C   * SUBROUTINE border_aggregate
!C   * SUBROUTINE aggregate
!C   * SUBROUTINE assign_aggregate_table
!C   ************************************************

MODULE aggregate_mod
CONTAINS
  SUBROUTINE indpdt_agrgt ( PNI,NI,LEVEL_NO,node_index, &
       & in_aggregates_result,aggregates_result, in_aggregates_result_size, &
       & border_flag)
    USE data_structure_for_AMG
    USE queue_mod
    IMPLICIT NONE
    
    INTEGER(kind=kint), INTENT(in) :: NI(:) 
    INTEGER(kind=kint), INTENT(in) :: PNI(0:)
    INTEGER(kind=kint), INTENT(in) :: LEVEL_NO
    INTEGER(kind=kint), INTENT(inout) :: node_index(:)
    INTEGER(kind=kint), pointer :: in_aggregates_result(:)
    INTEGER(kind=kint), pointer :: aggregates_result(:)
    INTEGER(kind=kint), INTENT(out) :: in_aggregates_result_size
    LOGICAL, INTENT(in) :: border_flag

    !C    TYPE(INTER_LEVEL_OPERATOR),POINTER :: P
    INTEGER(kind=kint), ALLOCATABLE :: Temp(:,:)
    
    INTEGER(kind=kint), ALLOCATABLE :: P_queue(:,:)
    INTEGER(kind=kint), ALLOCATABLE :: aggregate_flags(:)

    INTEGER(kind=kint) :: head,tail,queue_size
    
    REAL   (kind=kreal), POINTER :: D(:), AU(:), AL(:)
    INTEGER(kind=kint ), POINTER :: INU(:), IAU(:), INL(:), IAL(:)
    
    INTEGER(kind=kint ), POINTER :: NOD_EXPORT(:)
    INTEGER(kind=kint ), POINTER :: STACK_EXPORT(:)
    INTEGER(kind=kint ) :: NEIBPETOT
    INTEGER(kind=kint ), ALLOCATABLE :: node_record(:)

    INTEGER(kind=kint) :: i, s, s_u, k, f, aggr, i_CN, size_of_aggregates, count, l, j, m, front, back
    INTEGER(kind=kint) :: node_record_counter, node, score, insert_point, max_score

    REAL   (kind=kreal) :: r
    INTEGER(kind=kint), parameter :: aggregate_number=1
    INTEGER(kind=kint), parameter :: aggregate_index=2
    logical :: flag

    s = HIERARCHICAL_DATA(LEVEL_NO - 1) % N 

    ALLOCATE(Temp(0:s, 2))
    ALLOCATE(P_queue(2, s))
    ALLOCATE(aggregate_flags(s))
    ALLOCATE(node_record(s))
    
    D   => HIERARCHICAL_DATA(LEVEL_NO - 1) % D
    INL => HIERARCHICAL_DATA(LEVEL_NO - 1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO - 1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO - 1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO - 1) % AU

    s_u = INU(s)
    
    NEIBPETOT =     HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NEIBPETOT
    if(NEIBPETOT > 0) then
       STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % STACK_EXPORT
       NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NO - 1) % COMM_TABLE % NOD_EXPORT
    end if


    !C aggregates is stored in P temporarily 
    !C    P=>HIERARCHICAL_DATA(LEVEL_NO) % P    
    
    !C Temp =0
    Temp(1:s, aggregate_number) = 0
    Temp(0:s, aggregate_index) = 0


    DO i = 1, s
       if(node_index(i) == 0) node_index(i) = 1
    END DO
    
    !C initialize Priority queue
    !C P_queue has score and node index. node_index has the index of P_queue
    !C -1 in node_index : selected as part of aggregate
    !C  0 in node_index : not selected as part of aggregate, but it's neighbor is aggregate.
    
    P_queue(1:2, 1:s) = 0
    aggregate_flags(1:s) = 0
    
    
    !C interface part is scored, and enqueued
    !C Temp(,aggregate_number) is used as buffer
    if(NEIBPETOT > 0 .and. border_flag) then
       DO count = 1, STACK_EXPORT(NEIBPETOT)
          j = NOD_EXPORT(count)+ZERO_ORIGIN
          DO k = PNI(2 * j - 2) + 1,PNI(2 * j - 1)
             l = IAL(NI(k))+ZERO_ORIGIN
             Temp(l, aggregate_number) = 1
          ENDDO
          DO k = PNI(2 * j - 1) + 1, PNI(2 * j)
             if(IAU(NI(k))+ZERO_ORIGIN > s) cycle
             l = IAU(NI(k))+ZERO_ORIGIN
             Temp(l, aggregate_number) = 1
          ENDDO
       END DO
       !C interface node is set to 0.
       DO count = 1, STACK_EXPORT(NEIBPETOT)
          j = NOD_EXPORT(count)+ZERO_ORIGIN
          Temp(j, aggregate_number) = 0
       END DO
    end if
    queue_size = 0
    DO count = 1, s
       if(Temp(count, aggregate_number) == 1) then
          
          !C inserting: score is 10,position is queue_size
          queue_size = queue_size + 1
          node_index(count) = queue_size
          P_queue(1, queue_size) = 10
          P_queue(2, queue_size) = count
       END if
    END DO
    
    !C enqueue nodes whose score is 1
    DO count = 1, s
       if(Temp(count, aggregate_number) == 1)then
          Temp(count, aggregate_number) = 0
          !C this case is enqueued
       else if(node_index(count) == 1)then
          queue_size = queue_size + 1
          node_index(count) = queue_size
          P_queue(1, queue_size) = 1
          P_queue(2, queue_size) = count
       end if
    END DO
    
    size_of_aggregates = 0

    if(s > 0) then
       DO while(P_queue(1,1) > 0)

          i = P_queue(2,1)
          size_of_aggregates = size_of_aggregates + 1

          !C----layer1 start
          DO k = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
             j=IAL(NI(k))+ZERO_ORIGIN

             Temp(j, aggregate_number) = size_of_aggregates

             if(node_index(j) > 0) then
                !C dequeue node j
                call dequeue(j, P_queue, node_index, queue_size)
             end if
             node_index(j) = -1
             aggregate_flags(j) = i

          END DO
          DO k = PNI(2 * i - 1) + 1, PNI(2 * i)
             if(IAU(NI(k))+ZERO_ORIGIN > s) cycle
             j = IAU(NI(k))+ZERO_ORIGIN
             Temp(j, aggregate_number) = size_of_aggregates

             if(node_index(j) > 0) then
                !C dequeue node j
                call dequeue(j, P_queue, node_index, queue_size)
             end if
             node_index(j) = -1
             aggregate_flags(j) = i
          END DO

          Temp(i , aggregate_number) = size_of_aggregates

          !C dequeue node i
          call dequeue(i, P_queue, node_index, queue_size)

          node_index(i) = -1
          aggregate_flags(i) = i
          !C----layer1 end

          !C----layer2 and 3 
          call dequeue_scoreup_Layer23(i, s, queue_size, P_queue,  & 
               & node_index, aggregate_flags, node_record, NI, PNI, IAL, IAU)

       END DO
    endif

!!$       !C greedy case
!!$       !C== creating independent aggregate
!!$       size_of_aggregates=0
!!$       DO i=s,1,-1
!!$          IF(node_index(i)==0) THEN
!!$             l=0
!!$             DO k=PNI(2*i-2)+1,PNI(2*i-1)
!!$                j=IAL(NI(k))
!!$                l=l+node_index(j)
!!$             END DO
!!$             DO k=PNI(2*i-1)+1,PNI(2*i)
!!$                j=IAU(NI(k))
!!$                l=l+node_index(j)
!!$             END DO
!!$             IF(l==0) THEN
!!$                size_of_aggregates=size_of_aggregates+1
!!$
!!$                DO k=PNI(2*i-2)+1,PNI(2*i-1)
!!$                   j=IAL(NI(k))
!!$                   Temp(j, aggregate_number) = size_of_aggregates
!!$                   node_index(j)=-1
!!$                END DO
!!$                DO k=PNI(2*i-1)+1,PNI(2*i)
!!$                   j=IAU(NI(k))
!!$                   Temp(j, aggregate_number) = size_of_aggregates
!!$                   node_index(j)=-1
!!$                END DO
!!$
!!$                Temp(i, aggregate_number) = size_of_aggregates
!!$                node_index(i)=-1
!!$
!!$             END IF
!!$          END IF
!!$       END DO
!!$       write(*,*) size_of_aggregates

    
    deallocate(node_record)
    
    aggregate_check: DO count = 1, s
       if(node_index(count) == 0) then
          lower_part: DO k = PNI(2 * count - 2) + 1, PNI(2 * count - 1)
             !C aggr is the aggregate number of next to count
             aggr = Temp(IAL( NI(k) )+ZERO_ORIGIN, aggregate_number)
             IF(aggr > 0) Then
                Temp(count, aggregate_number) = aggr
                exit lower_part
             END IF
          END DO lower_part

          IF(Temp(count, aggregate_number) > 0) Then
             cycle aggregate_check
          END IF
          
          upper_part: DO k = PNI(2 * count - 1) + 1, PNI(2 * count)
             !C aggr is the aggregate number of next to i
             if(IAU(NI(k))+ZERO_ORIGIN > s) cycle
             aggr = Temp(IAU(NI(k))+ZERO_ORIGIN, aggregate_number)
             
             IF(aggr > 0) Then
                Temp(count, aggregate_number) = aggr
                exit upper_part
             END IF
          END DO upper_part
          IF(Temp(count, aggregate_number) == 0) Then
             write(*,*) "hehehe", node_index(count)
          END IF
       end if
    END DO aggregate_check

    DEALLOCATE(P_queue)
  
    in_aggregates_result_size = size_of_aggregates
    ALLOCATE(in_aggregates_result(0:size_of_aggregates))
    ALLOCATE(aggregates_result(s))

    in_aggregates_result(0) = 0

    DO i = 1, s
       aggr = Temp(i, aggregate_number)
       if(aggr > 0) Temp(aggr, aggregate_index) = Temp(aggr, aggregate_index) + 1
    END DO

    !C size_of_aggregates has the size of aggregates
    !C P % IN(aggr)+1 points the place for new data of aggr
    DO i = 1, size_of_aggregates
       Temp(i, aggregate_index) = Temp(i, aggregate_index) + Temp(i - 1, aggregate_index)
       in_aggregates_result(i) = Temp(i - 1, aggregate_index)
    END DO

    !C unneeded nodes belong to aggregate NO 0.
    DO i = 1, s
       aggr = Temp(i, aggregate_number)
       if(aggr > 0) then
          k = in_aggregates_result(aggr) + 1
          in_aggregates_result(aggr) = k
          aggregates_result(k) = i
       end if
    END DO
#ifdef PRINT_DATA_CREATE
    write(*,*) size_of_aggregates, in_aggregates_result(size_of_aggregates)
#endif

    i_CN = Temp(size_of_aggregates, aggregate_index)

    DEALLOCATE(Temp)

  END SUBROUTINE indpdt_agrgt


  SUBROUTINE border_aggregate (PNI, NI, LEVEL_NO, node_index, sorted_NEIBPE, my_rank,     &
       &                      PE_num_size,PE_nums,nodes_for_each_PE,in_nodes_for_each_PE, &
       &                      PE_list_size,PE_list,NPROCS,SOLVER_COMM)
    USE data_structure_for_AMG
    USE solver_Gnumbering
    USE solver_SR2
    IMPLICIT NONE
    INCLUDE  'mpif.h'
    INTEGER(kind=kint), INTENT(in)    :: NI(:)
    INTEGER(kind=kint), INTENT(in)    :: PNI(0:)
    INTEGER(kind=kint), INTENT(in)    :: LEVEL_NO
    INTEGER(kind=kint), INTENT(inout) :: node_index(:)
    INTEGER(kind=kint), INTENT(in)    :: sorted_NEIBPE(:)
    INTEGER(kind=kint), INTENT(in)    :: my_rank
    INTEGER(kind=kint), INTENT(in)    :: SOLVER_COMM
    INTEGER(kind=kint), INTENT(in)     :: PE_num_size
    INTEGER(kind=kint), POINTER :: PE_nums(:)
    INTEGER(kind=kint), POINTER :: nodes_for_each_PE(:)
    INTEGER(kind=kint), POINTER :: in_nodes_for_each_PE(:)
    INTEGER(kind=kint), INTENT(out) :: PE_list_size
    INTEGER(kind=kint), INTENT(out) :: PE_list(:,:)
    INTEGER(kind=kint), INTENT(in)  :: NPROCS
    
    INTEGER(kind=kint)          :: NEIBPETOT
    INTEGER(kind=kint), POINTER :: NEIBPE      (:)
    INTEGER(kind=kint), POINTER :: NOD_IMPORT  (:),STACK_IMPORT(:)
    INTEGER(kind=kint), POINTER :: NOD_EXPORT  (:),STACK_EXPORT(:)
    INTEGER(kind=kint)          :: ierr
    INTEGER(kind=kint), POINTER :: IAU (:)
    INTEGER(kind=kint), POINTER :: IAL (:)

    
    INTEGER(kind=kint) :: i,j,k,l,pe,node,NP,N,sni,pli,nfei,counter
    INTEGER(kind=kint) :: des_pe,nd,nds_sz_for_each_PE,in_s,sz,rbd
    INTEGER(kind=kint) :: index,Neib_index,Pe_num_index,sbi,nd_base
    INTEGER(kind=kint) :: PNI_L,PNI_H,PNI_M,sb_base,nni,sz_in,nd_in
    INTEGER(kind=kint) :: smaller_PEs

    INTEGER(kind=kint) :: inum,istart,recv_count
    INTEGER(kind=kint), ALLOCATABLE :: change_vector(:)
    INTEGER(kind=kint), ALLOCATABLE :: inverse_vector(:)
    INTEGER(kind=kint)              :: inverse_vector_size
    INTEGER(kind=kint)              :: send_buf_size
    INTEGER(kind=kint), ALLOCATABLE :: send_buf(:)
    INTEGER(kind=kint), ALLOCATABLE :: in_send_buf(:)
    INTEGER(kind=kint)              :: recv_buf_size
    INTEGER(kind=kint), ALLOCATABLE :: recv_buf(:)
    INTEGER(kind=kint), ALLOCATABLE :: recv_cnt_ary(:)
    INTEGER(kind=kint), ALLOCATABLE :: recv_buf2(:) 
    INTEGER(kind=kint), ALLOCATABLE :: recv_from_PE(:)
    
    INTEGER(kind=kint), ALLOCATABLE :: req1(:),req2(:),sta1(:,:),sta2(:,:)
    INTEGER(kind=kint), ALLOCATABLE :: PE_in_charge(:)
    INTEGER(kind=kint), ALLOCATABLE :: WSI(:),WRI(:)
    LOGICAL :: tf

    NEIBPETOT    =  HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPETOT     
    if(NEIBPETOT > 0) then
       NEIBPE       => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPE
       STACK_IMPORT => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_IMPORT
       NOD_IMPORT   => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_IMPORT
       STACK_EXPORT => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_EXPORT
       NOD_EXPORT   => HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_EXPORT
    end if
    NP           =  HIERARCHICAL_DATA(LEVEL_NO-1) % NP
    N            =  HIERARCHICAL_DATA(LEVEL_NO-1) % N
    IAU          => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    IAL          => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    
    !C PE_list(*,1): PE number
    !C PE_list(*,2): NEIBPE number|
    !C PE_list(*,3): PE_num number|these two has 0 if no correspondings.
    !C-- making PE_list..
    index=1
    Neib_index=1
    Pe_num_index=1

    DO WHILE((Neib_index<=NEIBPETOT) .AND. (Pe_num_index<=PE_num_size))
       i=sorted_NEIBPE(Neib_index)
       IF(PE_nums(Pe_num_index)==NEIBPE(i))THEN
          PE_list(index,1) = PE_nums(Pe_num_index)
          PE_list(index,2) = i
          PE_list(index,3) = Pe_num_index
          Neib_index = Neib_index+1
          PE_num_index = Pe_num_index+1
       ELSE IF(PE_nums(Pe_num_index)<NEIBPE(i))THEN
          PE_list(index,1) = PE_nums(Pe_num_index)
          PE_list(index,2) = 0
          PE_list(index,3) = Pe_num_index
          PE_num_index = PE_num_index+1
       ELSE      
          PE_list(index,1) = NEIBPE(i)
          PE_list(index,2) = i
          PE_list(index,3) = 0
          Neib_index = Neib_index+1
       END IF
       index=index+1
    END DO
    
    IF(Neib_index>NEIBPETOT)THEN
       DO i=PE_num_index,PE_num_size
          PE_list(index,1)=PE_nums(i)
          PE_list(index,2)=0
          PE_list(index,3)=i
          index=index+1
       END DO
    ELSE
       DO i=Neib_index,NEIBPETOT
          j=sorted_NEIBPE(i)
          PE_list(index,1)=NEIBPE(j)
          PE_list(index,2)=j
          PE_list(index,3)=0
          index=index+1
       END DO
    END IF
    PE_list_size = index - 1
    !C-- end: making PE_list..

    !C-- PE_in_charge has its PE in charge
    ALLOCATE(PE_in_charge(NP))

    DO i = 1, N
       PE_in_charge(i) = my_rank
    END DO
    DO i = NEIBPETOT, 1, -1
       j = sorted_NEIBPE(i)
       pe = NEIBPE(j)
       IF(pe > my_rank) pe = my_rank
       DO k = STACK_EXPORT(j - 1) + 1, STACK_EXPORT(j)
          node = NOD_EXPORT(k)+ZERO_ORIGIN
          PE_in_charge(node) = pe
       END DO
    END DO
    
    ALLOCATE(WSI(STACK_EXPORT(NEIBPETOT)))
    ALLOCATE(WRI(STACK_IMPORT(NEIBPETOT)))
    CALL SOLVER_SEND_RECV2I(NP, STACK_EXPORT(NEIBPETOT),                 &
         &       STACK_IMPORT(NEIBPETOT), NEIBPETOT, NEIBPE(1:NEIBPETOT),&
         &       STACK_IMPORT(0:NEIBPETOT), &
         &       NOD_IMPORT(1:STACK_IMPORT(NEIBPETOT)), &
	 &       STACK_EXPORT(0:NEIBPETOT), &
         &       NOD_EXPORT(1:STACK_EXPORT(NEIBPETOT)), &
         &       WSI(1:STACK_EXPORT(NEIBPETOT)), &
         &       WRI(1:STACK_IMPORT(NEIBPETOT)), &
	 &       PE_in_charge(1:NP), SOLVER_COMM, my_rank)

    DEALLOCATE(WSI)
    DEALLOCATE(WRI)

    !C-- end: PE_in_charge has its PE in charge

    
    !C-- information of neighbors in COMM tables is exchanged
    !C   * change_vector has serial number of nod_import,
    !C    node_export, nodes_for_each_PE. irrelevant node has 0.
    ALLOCATE(change_vector(NP))
    send_buf_size = NEIBPETOT + STACK_EXPORT(NEIBPETOT) + PNI(2 * N)
    ALLOCATE(send_buf(send_buf_size))
    ALLOCATE(in_send_buf(0:NEIBPETOT))
    
    sbi = 0
    Pli = 1

    !C PE_list's PE number is sorted. so the loop is execed 
    !C on PE sorted order.
    DO i = 1, NEIBPETOT
       sni = sorted_NEIBPE(i)
       des_pe = NEIBPE(sni)
       DO WHILE(PE_list(pli, 1) /= des_pe)
          pli = pli + 1
       END DO
       IF(pli > PE_list_size) STOP "unexpected error in border_aggregates"
       nfei = PE_list(pli, 3)
       
       !C== making change_vector
       !C counter is used to create serial number
       !C change_vector: 0 or serial number on extended comm table
       change_vector = 0
       counter = 0
       DO j = STACK_IMPORT(sni - 1) + 1, STACK_IMPORT(sni)
          counter = counter + 1
#ifdef DEBUG
          IF(change_vector(NOD_IMPORT(j)+ZERO_ORIGIN) /= 0) & 
          & STOP "error in change_vector of border_aggregate()"
#endif
          change_vector(NOD_IMPORT(j)+ZERO_ORIGIN) = counter
       END DO

       !C nds_sz_for_each_PE: size of nodes to send from my_rank to des_pe
       nds_sz_for_each_PE = 0
       DO j = STACK_EXPORT(sni - 1) + 1, STACK_EXPORT(sni)
          counter = counter + 1
          nd = NOD_EXPORT(j)+ZERO_ORIGIN
#ifdef DEBUG
          IF(change_vector(NOD_EXPORT(j)+ZERO_ORIGIN) /= 0) &
          & STOP "error in change_vector of border_aggregate()"
#endif

          change_vector(nd) = counter
          IF(PE_in_charge(nd)==des_pe) nds_sz_for_each_PE=nds_sz_for_each_PE+1
       END DO
       IF(nfei > 0) THEN
          DO j = in_nodes_for_each_PE(nfei - 1) + 1, in_nodes_for_each_PE(nfei)
             counter = counter + 1
#ifdef DEBUG
             IF(change_vector(nodes_for_each_PE(j)+ZERO_ORIGIN)/=0)  &
             & STOP "error in change_vector of border_aggregate()"
#endif
             change_vector(nodes_for_each_PE(j)+ZERO_ORIGIN) = counter
          END DO
       END IF
       !C== end:making change_vector       

       !C== buffering each_node's neighbors for each PE
       !C   sbi: send buffer index
       sbi = sbi + 1

       send_buf(sbi) = nds_sz_for_each_PE
       !C sb_base means to be index of send_buf
       sb_base = sbi + nds_sz_for_each_PE
       DO j = STACK_EXPORT(sni - 1) + 1, STACK_EXPORT(sni)
          k = NOD_EXPORT(j)+ZERO_ORIGIN
          IF(PE_in_charge(k) == des_pe) THEN
             PNI_H = PNI(2 * k)
             PNI_M = PNI(2 * k - 1)
             PNI_L = PNI(2 * k - 2)
             !C neighboring info is thrown to the buffer
             sbi = sbi + 1
             send_buf(sbi) = PNI_H - PNI_L
             DO l = PNI_L + 1, PNI_M
                nni = NI(l)
                send_buf(sb_base + l - PNI_L) = change_vector(IAL(nni)+ZERO_ORIGIN)
             END DO
             sb_base = sb_base + PNI_M - PNI_L
             DO l = PNI_M + 1, PNI_H
                nni = NI(l)
                send_buf(sb_base + l - PNI_M) = change_vector(IAU(nni)+ZERO_ORIGIN)
             END DO
             sb_base = sb_base + PNI_H - PNI_M
          END IF
       END DO
       sbi = sb_base
       in_send_buf(i) = sbi
       !C== end:buffering each_node's neighbors for each PE
    END DO
    in_send_buf(0) = 0
    DEALLOCATE(change_vector)
    
    !C== count up recv_count, and allcate recv_buf
    !C maybe recv_buf_size is too small..
    !C recv_count is the number of receive times

    recv_count = 0
    DO i = 1, NEIBPETOT
       pe = NEIBPE(i)
       IF(pe > my_rank) recv_count = recv_count + 1
    END DO


    ALLOCATE(recv_from_PE(NEIBPETOT))
    ALLOCATE(recv_cnt_ary(0:recv_count))
    ALLOCATE(req1(PE_list_size))
    ALLOCATE(sta1(MPI_STATUS_SIZE,PE_list_size))
    ALLOCATE(req2(PE_list_size))
    ALLOCATE(sta2(MPI_STATUS_SIZE,PE_list_size))
    
    k = 0
    l = 0
    recv_from_PE = -1
    DO i = 1, NEIBPETOT
       sni = sorted_NEIBPE(i)
       des_pe = NEIBPE(sni)
       IF(des_pe > my_rank) THEN
          k = k + 1
          recv_from_PE(sni) = k
          CALL MPI_IRECV(recv_cnt_ary(k),1,LIS_MPI_INTEGER,des_pe,0,&
               &        SOLVER_COMM,req2(k),ierr)
       ELSE IF(des_pe<my_rank)THEN
          l=l+1
          istart=in_send_buf(i-1)
          inum=in_send_buf(i)-istart
          CALL MPI_ISEND (inum,1,LIS_MPI_INTEGER,des_pe,0,&
               &    SOLVER_COMM,req1(l),ierr)
       ELSE
          WRITE(*,*) "error!",des_pe,my_rank
          STOP
       END IF
    END DO
    IF(l > 0) CALL MPI_WAITALL(l, req1, sta1, ierr)
    IF(k > 0) CALL MPI_WAITALL(k, req2, sta2, ierr)
    recv_cnt_ary(0) = 0
    DO i = 1, recv_count
       recv_cnt_ary(i) = recv_cnt_ary(i - 1) + recv_cnt_ary(i)
    END DO
    recv_buf_size = recv_cnt_ary(recv_count)
    ALLOCATE(recv_buf(recv_buf_size))
    !C== end:count up recv_count, and allcate recv_buf    



    k = 0
    l = 0
    DO i = 1, NEIBPETOT
       sni = sorted_NEIBPE(i)
       des_pe = NEIBPE(sni)
       IF(des_pe > my_rank) THEN
          k = k + 1
          istart = recv_cnt_ary(k - 1)
          inum = recv_cnt_ary(k) - istart
          CALL MPI_IRECV(recv_buf(istart+1),inum,LIS_MPI_INTEGER,des_pe,0,&
               &        SOLVER_COMM,req2(k),ierr)
       ELSE
          l = l + 1
          istart = in_send_buf(i - 1)
          inum = in_send_buf(i) - istart
          CALL MPI_ISEND (send_buf(istart+1),inum,LIS_MPI_INTEGER,des_pe,0,&
               &    SOLVER_COMM,req1(l),ierr)
       END IF
    END DO
    IF(k > 0) CALL MPI_WAITALL(k, req2, sta2, ierr)
    IF(l > 0) CALL MPI_WAITALL(l, req1, sta1, ierr)
    DEALLOCATE(send_buf)
    DEALLOCATE(in_send_buf)
    !C-- end: infomation of neighbors in COMM tables is exchanged



    
    !C-- exec aggregation in COMM tables
    pli=1
    DO i=1,NEIBPETOT
       sni = sorted_NEIBPE(i)
       des_PE=NEIBPE(sni)
       IF(my_rank<des_PE) THEN
          DO WHILE(PE_list(pli,1)/=des_pe)
             pli=pli+1
          END DO
          IF(pli>PE_list_size) STOP "unexpected error in border_aggregates"
          nfei=PE_list(pli,3)
          
          !C== inverse_vector creation
          inverse_vector_size=STACK_EXPORT(NEIBPETOT)+STACK_IMPORT(NEIBPETOT)
          IF(PE_num_size/=0) inverse_vector_size=inverse_vector_size+        &
               & in_nodes_for_each_PE(PE_num_size)

          ALLOCATE(inverse_vector(inverse_vector_size))
          k=0
          DO j=STACK_EXPORT(sni-1)+1,STACK_EXPORT(sni)
             k=k+1
             inverse_vector(k)=NOD_EXPORT(j)+ZERO_ORIGIN
          END DO
          DO j=STACK_IMPORT(sni-1)+1,STACK_IMPORT(sni)
             k=k+1
             inverse_vector(k)=NOD_IMPORT(j)+ZERO_ORIGIN
          END DO
          IF(nfei>0) THEN
             DO j=in_nodes_for_each_PE(nfei-1)+1,in_nodes_for_each_PE(nfei)
                k=k+1
                inverse_vector(k)=nodes_for_each_PE(j)+ZERO_ORIGIN
             END DO
          END IF
          !C== end inverse_vector creation


          
          
          !C== aggregation on border
          DO j=STACK_EXPORT(sni-1)+1,STACK_EXPORT(sni)
             tf=.TRUE.
             node=NOD_EXPORT(j)+ZERO_ORIGIN
             IF(PE_in_charge(node)==my_rank .AND. node_index(node)==0) THEN
                l=node_index(node)
                DO k=PNI(2*node-2)+1,PNI(2*node-1)
                   index=IAL(NI(k))+ZERO_ORIGIN
                   l=l+node_index(index)
                   tf=tf.AND.(my_rank<=PE_in_charge(index))
                END DO
                DO k=PNI(2*node-1)+1,PNI(2*node)
                   index=IAU(NI(k))+ZERO_ORIGIN
                   l=l+node_index(index)
                   tf=tf.AND.(my_rank<=PE_in_charge(index))
                END DO
                IF(l==0 .AND. tf) THEN
                   node_index(node)=-2
                   DO k=PNI(2*node-2)+1,PNI(2*node-1)
                      node_index(IAL(NI(k))+ZERO_ORIGIN)=-1
                   END DO
                   DO k=PNI(2*node-1)+1,PNI(2*node)
                      node_index(IAU(NI(k))+ZERO_ORIGIN)=-1
                   END DO
                END IF
             END IF
          END DO


          
          k=recv_from_PE(sni)
          IF(k>0) THEN
             nds_sz_for_each_PE=recv_buf(recv_cnt_ary(k-1)+1)
             nd_base=nds_sz_for_each_PE+1+recv_cnt_ary(k-1)
!!$             recv_buf(1,k)=0

             nd_in=nd_base
             sz_in=1+recv_cnt_ary(k-1)
             DO j=STACK_IMPORT(sni-1)+1,STACK_IMPORT(sni)
                node=NOD_IMPORT(j)+ZERO_ORIGIN
                IF(PE_in_charge(node)==my_rank) THEN
                   in_s=nd_in+1
                   sz_in=sz_in+1
                   sz=recv_buf(sz_in)
                   nd_in=nd_in+sz
                   
                   tf=.TRUE.
                   l=node_index(node)
                   DO index=in_s,nd_in
                      rbd=recv_buf(index)
                      IF(rbd>0) THEN
                         nd=inverse_vector(rbd)
                         l=l+node_index(nd)
                         tf=tf.AND.(my_rank<=PE_in_charge(nd))
                      END IF
                   END DO
                   IF(l==0 .AND. tf) THEN
                      node_index(node)=-2
                      DO index=in_s,nd_in
                         rbd=recv_buf(index)
                         IF(rbd>0) node_index(inverse_vector(rbd)) = -1
                      END DO
                   END IF
                END IF
             END DO
          END IF
          !C== end: aggregation on border
          DEALLOCATE(inverse_vector)
       END IF
    END DO
    DEALLOCATE(recv_buf)
    DEALLOCATE(recv_cnt_ary)
    DEALLOCATE(recv_from_PE)
    
    !C-- end: exec aggregation in COMM tables




    !C== the result of aggregation is communicated
    k=0
    l=0
    send_buf_size=STACK_IMPORT(NEIBPETOT)+1
    ALLOCATE(send_buf(send_buf_size))
    recv_buf_size=STACK_EXPORT(NEIBPETOT)+1
    ALLOCATE(recv_buf2(recv_buf_size))

    DO i=1,NEIBPETOT
       des_PE=NEIBPE(i)
       IF(des_pe>my_rank) THEN
          l=l+1
          istart=STACK_IMPORT(i-1)
          inum=STACK_IMPORT(i)-istart
          DO j=istart+1,istart+inum
             send_buf(j)=node_index(NOD_IMPORT(j)+ZERO_ORIGIN)
          END DO
          !C MPI_ISEND
          CALL MPI_ISEND(send_buf(istart+1), inum, LIS_MPI_INTEGER, &
               &  des_PE,0,SOLVER_COMM,req1(l),ierr)
       ELSE
          k=k+1
          istart=STACK_EXPORT(i-1)
          inum=STACK_EXPORT(i)-istart
          !C MPI_IRECV
          CALL MPI_IRECV(recv_buf2(istart+1),inum,LIS_MPI_INTEGER,  &
               &  des_PE, 0, SOLVER_COMM, req2(k), ierr)
       END IF
    END DO
    IF(k>0) CALL MPI_WAITALL(k,req2,sta2,ierr)
    IF(l>0) CALL MPI_WAITALL(l,req1,sta1,ierr)
    DEALLOCATE(send_buf)    



    !C-- assign by large to small number order of PE 
    DO i=NEIBPETOT,1,-1
       sni = sorted_NEIBPE(i)
       IF(NEIBPE(sni)<my_rank) THEN
          istart = STACK_EXPORT(sni-1)
          inum   = STACK_EXPORT(sni)-istart
          DO j=istart+1, istart+inum
             k=NOD_EXPORT(j)+ZERO_ORIGIN
             node_index(k)=recv_buf2(j)
          END DO
       END IF
    END DO
    DEALLOCATE(recv_buf2)
    !C== end: the result of aggregation is communicated



    
    !C-- overlooked nodes in export table by dependency of PEs are checked.
    !C   PE_list's have information  in PE number order.
    
    ALLOCATE(in_send_buf(0:PE_list_size))
    in_send_buf(0)=0
    counter=0
    smaller_PEs=0
    DO i=1,PE_list_size
       des_pe=PE_list(i,1)
       Neib_index  =PE_list(i,2)
       Pe_num_index=PE_list(i,3)
       IF(Neib_index>0) THEN
          counter=counter+STACK_EXPORT(Neib_index)-STACK_EXPORT(Neib_index-1)&
               &         +STACK_IMPORT(Neib_index)-STACK_IMPORT(Neib_index-1)
       END IF
       IF(Pe_num_index>0) THEN
          counter=counter+in_nodes_for_each_PE(Pe_num_index) &
               &         -in_nodes_for_each_PE(Pe_num_index-1)
       END IF
       in_send_buf(i)=counter
       IF(des_pe<my_rank) smaller_PEs=smaller_PEs+1
    END DO
    
    !C-- receive the message from smaller PEs and assign node_index
    !C   assgining to the node_index in decendant order.
    ALLOCATE(recv_buf2(in_send_buf(smaller_PEs)+1))
    DO i=1,smaller_PEs
       des_PE=PE_list(i,1)
       istart=in_send_buf(i-1)
       inum=in_send_buf(i)-istart
#ifdef DEBUG
       IF(inum==0) WRITE(*,*) "unexpected error in border_aggregate:2"
#endif

       CALL MPI_IRECV(recv_buf2(istart+1),inum,LIS_MPI_INTEGER,des_PE,0,&
            &         SOLVER_COMM,req2(i),ierr)
    END DO
    IF(smaller_PEs>0) CALL MPI_WAITALL(smaller_PEs,req2,sta2,ierr)
     
    
    !C order of PE_list is assumed.
    DO i=1,smaller_PEs
       Neib_index  =PE_list(i,2)
       Pe_num_index=PE_list(i,3)
       k=in_send_buf(i-1)
       IF(Neib_index>0) THEN
          DO j=STACK_EXPORT(Neib_index-1)+1,STACK_EXPORT(Neib_index)
             k=k+1
             node_index(NOD_EXPORT(j)+ZERO_ORIGIN)=recv_buf2(k)
          END DO
          DO j=STACK_IMPORT(Neib_index-1)+1,STACK_IMPORT(Neib_index)
             k=k+1
             node_index(NOD_IMPORT(j)+ZERO_ORIGIN)=recv_buf2(k)
          END DO
       END IF
       IF(Pe_num_index>0) THEN
          DO j=in_nodes_for_each_PE(Pe_num_index-1)+1,in_nodes_for_each_PE(Pe_num_index)
             k=k+1
             node_index(nodes_for_each_PE(j)+ZERO_ORIGIN)=recv_buf2(k)
          END DO
       END IF
    END DO
    DEALLOCATE(recv_buf2)
    !C-- end:receive the message from smaller PEs and assign node_index

    
    !C-- aggregation on EXPORT_NOD.
    DO i=1,NEIBPETOT
       DO j=STACK_EXPORT(i-1)+1,STACK_EXPORT(i)
          node=NOD_EXPORT(j)+ZERO_ORIGIN
          IF(node_index(node)==0) THEN
             l=0
             DO k=PNI(2*node-2)+1,PNI(2*node-1)
                index=IAL(NI(k))+ZERO_ORIGIN
                l=l+node_index(index)
             END DO
             DO k=PNI(2*node-1)+1,PNI(2*node)
                index=IAU(NI(k))+ZERO_ORIGIN
                l=l+node_index(index)
             END DO
             IF(l==0) THEN
                node_index(node)=-2
                DO k=PNI(2*node-2)+1,PNI(2*node-1)
                   node_index(IAL(NI(k))+ZERO_ORIGIN)=-1
                END DO
                DO k=PNI(2*node-1)+1,PNI(2*node)
                   node_index(IAU(NI(k))+ZERO_ORIGIN)=-1
                END DO
             END IF
          END IF
       END DO
    END DO
    !C-- end:aggregation.
    


    !C-- send the message of node_index to bigger PEs
    send_buf_size=in_send_buf(PE_list_size)-in_send_buf(smaller_PEs)+1
    ALLOCATE(send_buf(send_buf_size))
    k=0
    DO i=smaller_PEs+1,PE_list_size
       Neib_index  =PE_list(i,2)
       PE_num_index=PE_list(i,3)
       des_PE=PE_list(i,1)
       k=in_send_buf(i-1)-in_send_buf(smaller_PEs)
       IF(Neib_index>0) THEN
          DO j=STACK_IMPORT(Neib_index-1)+1,STACK_IMPORT(Neib_index)
             k=k+1
             send_buf(k) = node_index(NOD_IMPORT(j)+ZERO_ORIGIN)
          END DO
          DO j=STACK_EXPORT(Neib_index-1)+1,STACK_EXPORT(Neib_index)          
             k=k+1
             send_buf(k) = node_index(NOD_EXPORT(j)+ZERO_ORIGIN)
          END DO
       END IF
       IF(Pe_num_index>0) THEN
          DO j=in_nodes_for_each_PE(Pe_num_index-1)+1,in_nodes_for_each_PE(Pe_num_index)
             k=k+1
             send_buf(k) = node_index(nodes_for_each_PE(j)+ZERO_ORIGIN)
          END DO
       END IF

       istart = in_send_buf(i-1)-in_send_buf(smaller_PEs)
       inum   = in_send_buf(i)  -in_send_buf(smaller_PEs)-istart

       CALL MPI_ISEND(send_buf(istart+1),inum,LIS_MPI_INTEGER, &
            & des_PE,0,SOLVER_COMM,req1(i-smaller_PEs),ierr)
       
    END DO

    


    IF(PE_list_size-smaller_PEs>0) THEN
       CALL MPI_WAITALL(PE_list_size-smaller_PEs,req1,sta1,ierr)
    END IF
    DEALLOCATE(send_buf)
    !C-- end:send the message of node_index to bigger PEs
    !C-- end:overlooked nodes in export table by dependency of PEs are checked.


    !C-- bigger PE has correct node_Index. So make node_index correct by communication 
    !C   the result of aggregation is communicated
    k = 0
    l = 0
    send_buf_size = STACK_IMPORT(NEIBPETOT)
    ALLOCATE(send_buf(send_buf_size))
    recv_buf_size = STACK_EXPORT(NEIBPETOT) + 1
    ALLOCATE(recv_buf2(recv_buf_size))

    DO i = 1, NEIBPETOT
       des_PE = NEIBPE(i)
       IF(des_pe < my_rank) THEN
          l = l + 1
          istart = STACK_IMPORT(i - 1)
          inum   = STACK_IMPORT(i) - istart
          DO j = istart + 1, istart + inum
             send_buf(j) = node_index(NOD_IMPORT(j)+ZERO_ORIGIN)
          END DO
          !C MPI_ISEND
          CALL MPI_ISEND(send_buf(istart+1), inum, LIS_MPI_INTEGER, &
               &  des_PE,0,SOLVER_COMM,req1(l),ierr)
       ELSE
          k = k + 1
          istart = STACK_EXPORT(i - 1)
          inum   = STACK_EXPORT(i) - istart
          !C MPI_IRECV
          CALL MPI_IRECV(recv_buf2(istart+1),inum,LIS_MPI_INTEGER,  &
               &  des_PE, 0, SOLVER_COMM, req2(k), ierr)
       END IF
    END DO
    IF(k > 0) CALL MPI_WAITALL(k, req2, sta2, ierr)
    IF(l > 0) CALL MPI_WAITALL(l, req1, sta1, ierr)
    DEALLOCATE(send_buf)    




    !C-- assign by small to large number order of PE 
    DO i = 1, NEIBPETOT
       sni = sorted_NEIBPE(i)
       IF(NEIBPE(sni) > my_rank) THEN
          istart = STACK_EXPORT(sni - 1)
          inum = STACK_EXPORT(sni) - istart
          DO j = istart + 1, istart + inum
             k = NOD_EXPORT(j)+ZERO_ORIGIN
             node_index(k) = recv_buf2(j)
          END DO
       END IF
    END DO
    DEALLOCATE(recv_buf2)
    !C-- end:bigger PE has correct node_index. So make node_index correct by communication 


    
    !C node_index is correct in the own domain, not in the domain of NOD_IMPORT.

    DEALLOCATE(PE_in_charge)
    DEALLOCATE(req1)
    DEALLOCATE(req2)
    DEALLOCATE(sta1)
    DEALLOCATE(sta2)

  END SUBROUTINE border_aggregate

  SUBROUTINE aggregate(PNI, NI, LEVEL_NO, node_index, in_aggregates_result,           &
       &               aggregates_result, in_aggregates_result_size, my_rank,         &
       &               aggregate_table_size, aggregate_table_array,                   &
       &               aggregate_number_in_table, SOLVER_COMM, NPROCS, GIN_aggregate, &
       &               PE_num_size, PE_nums, nodes_for_each_PE, in_nodes_for_each_PE, &
       &               PE_list_size, PE_list, kinds)
    
    USE data_structure_for_AMG
    USE queue_mod
    IMPLICIT NONE
    INCLUDE  'mpif.h'

    INTEGER(kind=kint), INTENT(in) :: NI(:),PNI(0:),LEVEL_NO
    INTEGER(kind=kint), INTENT(inout) :: node_index(:)
    INTEGER(kind=kint), pointer :: in_aggregates_result(:)
    INTEGER(kind=kint), pointer :: aggregates_result(:)
    INTEGER(kind=kint), INTENT(out) :: in_aggregates_result_size
    INTEGER(kind=kint), INTENT(in) :: my_rank, PE_list_size
    INTEGER(kind=kint), INTENT(in) :: PE_list(:, :)
    !C kinds = 1,or 2,or 3
    INTEGER(kind=kint), INTENT(in) :: kinds

    REAL   (kind=kreal), POINTER :: AL(:), AU(:), D(:)
    INTEGER(kind=kint), POINTER :: INL(:), IAU(:)
    INTEGER(kind=kint), POINTER :: IAL(:), INU(:)
    INTEGER(kind=kint), POINTER :: NOD_EXPORT(:), NOD_IMPORT(:), NEIBPE(:)
    INTEGER(kind=kint), POINTER :: STACK_EXPORT(:), STACK_IMPORT(:)
    INTEGER(kind=kint) :: NEIBPETOT, N, NP

    INTEGER(kind=kint) :: i, j, k, l, aggr, size_of_aggregates
    INTEGER(kind=kint), allocatable :: Temp(:,:)

    INTEGER(kind=kint),parameter :: aggregate_number = 1
    INTEGER(kind=kint),parameter :: aggregate_index  = 2
    
    INTEGER(kind=kint), ALLOCATABLE :: P_queue(:,:)
    INTEGER(kind=kint), ALLOCATABLE :: aggregate_flags(:)
    INTEGER(kind=kint), ALLOCATABLE :: node_record(:)
    INTEGER(kind=kint) :: queue_size
    
    INTEGER(kind=kint),intent(in)  :: aggregate_table_size, NPROCS
!!$    INTEGER(kind=kint),intent(out) :: aggregate_table_array(aggregate_table_size * 2, NPROCS)
!!$    INTEGER(kind=kint),intent(out) :: aggregate_number_in_table(NPROCS)
    INTEGER(kind=kint),intent(out) :: aggregate_table_array(:, :)
    INTEGER(kind=kint),intent(out) :: aggregate_number_in_table(:)

    INTEGER(kind=kint),intent(in)  :: SOLVER_COMM
    INTEGER(kind=kint)             :: ierr
    INTEGER(kind=kint),intent(out) :: GIN_aggregate(0:NPROCS)

    INTEGER(kind=kint), intent(in) :: PE_num_size
    INTEGER(kind=kint), pointer    :: PE_nums(:)
    INTEGER(kind=kint), pointer    :: nodes_for_each_PE(:)
    INTEGER(kind=kint), pointer    :: in_nodes_for_each_PE(:)
    
    INTEGER(kind=kint), allocatable :: send_buf(:)
    INTEGER(kind=kint), allocatable :: in_send_buf_for_PEs(:)
    INTEGER(kind=kint), allocatable :: recv_buf(:)
    INTEGER(kind=kint), allocatable :: req1(:),req2(:)
    INTEGER(kind=kint), allocatable :: sta1(:,:), sta2(:,:)

    
    INTEGER(kind=kint) :: inum,istart,maxbufsize,position
    INTEGER(kind=kint) :: Neib_index,PE_num_index,dest,m,lc
    INTEGER(kind=kint), allocatable :: inner_nd(:)

    N  = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP 
    
    D   => HIERARCHICAL_DATA(LEVEL_NO-1) % D
    INL => HIERARCHICAL_DATA(LEVEL_NO-1) % INL
    INU => HIERARCHICAL_DATA(LEVEL_NO-1) % INU
    IAL => HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU => HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL  => HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU  => HIERARCHICAL_DATA(LEVEL_NO-1) % AU

    NEIBPETOT   = HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPETOT
    if(NEIBPETOT > 0) then
       NEIBPE      =>HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NEIBPE
       NOD_EXPORT  =>HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_EXPORT
       STACK_EXPORT=>HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_EXPORT
       NOD_IMPORT  =>HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % NOD_IMPORT
       STACK_IMPORT=>HIERARCHICAL_DATA(LEVEL_NO-1) % COMM_TABLE % STACK_IMPORT
    END if


    if(kinds < 1 .or. kinds > 3) stop "error in aggregate()."    

    ALLOCATE(Temp(0:NP, 2))
    Temp = 0    

    if(kinds == 1) then
       !C greedy case
       !C== creating independent aggregate
       size_of_aggregates = 0
       DO i = 1, N
          IF(node_index(i) == -2) THEN
             size_of_aggregates=size_of_aggregates+1
             DO k=PNI(2*i-2)+1,PNI(2*i-1)
                j=IAL(NI(k))+ZERO_ORIGIN
                Temp(j, aggregate_number) = size_of_aggregates
                node_index(j)=-1
             END DO
             DO k=PNI(2*i-1)+1,PNI(2*i)
                j=IAU(NI(k))+ZERO_ORIGIN
                Temp(j, aggregate_number) = size_of_aggregates
                node_index(j)=-1
             END DO
             Temp(i, aggregate_number) = size_of_aggregates
          END IF
       END DO
       ALLOCATE(inner_nd(N))
       inner_nd = 1
       DO i = 1, STACK_EXPORT(NEIBPETOT)
          inner_nd(NOD_EXPORT(i)+ZERO_ORIGIN) = 0
       END DO
       DO i = 1, N
          IF(node_index(i) == 0 .AND. inner_nd(i) == 1) THEN
             l = 0
             DO k = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
                j = IAL(NI(k))+ZERO_ORIGIN
                l = l + node_index(j)
             END DO
             DO k = PNI(2 * i - 1) + 1, PNI(2 * i)
                j = IAU(NI(k))+ZERO_ORIGIN
                l = l + node_index(j)
             END DO
             IF(l == 0) THEN
                
                size_of_aggregates = size_of_aggregates + 1

                DO k = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
                   j = IAL(NI(k))+ZERO_ORIGIN
                   Temp(j, aggregate_number) = size_of_aggregates
                   node_index(j) = -1
                END DO
                
                DO k = PNI(2 * i - 1) + 1, PNI(2 * i)
                   j = IAU(NI(k))+ZERO_ORIGIN
                   Temp(j, aggregate_number) = size_of_aggregates
                   node_index(j) = -1
                END DO

                Temp(i, aggregate_number) = size_of_aggregates
                node_index(i) = -1

             END IF
          END IF
       END DO

       deallocate(inner_nd)
       !C== end: creating independent aggregates
    else if(kinds == 2 .or. kinds == 3 ) then
       !C These cases are priority queue is utilized

       ALLOCATE(P_queue(2,N))
       ALLOCATE(aggregate_flags(N))
       ALLOCATE(node_record(N))

       P_queue(1:2, 1:N)=0

       !C== enqueue internal nodes using aggregate_flags
       !C== node_index -2 -1 0: 
       !C== root, determined, undtermined
       !C== node_index 0 -> part of the queue data >0
       aggregate_flags(1:N) = 1
       DO i = 1, STACK_EXPORT(NEIBPETOT)
          aggregate_flags(NOD_EXPORT(i)+ZERO_ORIGIN) = 0
       END DO


       
       queue_size = 0
       DO i = 1, N 
          if(aggregate_flags(i) == 1 .and. node_index(i) /= -1) then
             queue_size = queue_size + 1
             node_index(i) = queue_size
             P_queue(1, queue_size) = 1
             P_queue(2, queue_size) = i

             aggregate_flags(i) = 0
          end if
       END DO


       !C== dequeue on nodes around the nodes on EXPORT TABLE with score -1 
       !C== node_index -2 -1 0 >0
       !C== root, determined, undetermined but not root, undetermined and part of queue
       DO l = 1, STACK_EXPORT(NEIBPETOT)
          i = NOD_EXPORT(l)+ZERO_ORIGIN
          if(node_index(i) == -1) then
             DO k = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
                j = IAL( NI(k) ) + ZERO_ORIGIN
                if(node_index(j) > 0) then

                   call dequeue(j, P_queue, node_index, queue_size)
                   node_index(j) = 0
                end if
             END DO
             DO k = PNI(2 * i - 1) + 1, PNI(2 * i)
                j=IAU( NI(k) )+ZERO_ORIGIN
                if(node_index(j) > 0) then
                   call dequeue(j, P_queue, node_index,queue_size)
                   node_index(j) = 0
                end if
             END DO
          end if

       END DO
       

       !C== creating aggregates which are predetermined on borders
       !C== dequeue and scoreup on nodes around the nodes on EXPORT TABLE with score -2 
       size_of_aggregates = 0
       DO i = 1, N
          if(node_index(i) == -2) then
             size_of_aggregates = size_of_aggregates + 1
             DO k = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
                j=IAL(NI(k)) + ZERO_ORIGIN
                Temp(j, aggregate_number) = size_of_aggregates
                if(node_index(j)>0) call dequeue(j, P_queue, node_index,queue_size)
                node_index(j) = -1
             END DO
             DO k = PNI(2 * i - 1) + 1,PNI(2 * i)
                j=IAU(NI(k))+ZERO_ORIGIN
                Temp(j, aggregate_number) = size_of_aggregates
                if(node_index(j) > 0) call dequeue(j, P_queue, node_index,queue_size)
                node_index(j) = -1
             END DO
             Temp(i, aggregate_number) = size_of_aggregates

             if(kinds==3) then
                !C-- this function weight the nodes around the aggregates on borders
                call dequeue_scoreup_Layer23(i, N, queue_size, P_queue,  & 
                     & node_index, aggregate_flags, node_record, NI, PNI, IAL, IAU)
             end if
          END if
       END DO

       
       !C== create aggregates internal part
       if(N > 0) then
          DO while(P_queue(1, 1) > 0)


             i = P_queue(2, 1)
             size_of_aggregates = size_of_aggregates + 1


             !C----layer1 start
             DO k = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
                j = IAL(NI(k)) + ZERO_ORIGIN

                Temp(j, aggregate_number) = size_of_aggregates
                

                if(node_index(j)>0) then
                   !C dequeue node j
                   call dequeue(j, P_queue, node_index, queue_size)
                end if
                node_index(j)=-1
                aggregate_flags(j)=i

             END DO

             DO k=PNI(2*i-1)+1,PNI(2*i)
                j=IAU( NI(k) )+ZERO_ORIGIN
                Temp(j, aggregate_number) = size_of_aggregates


                if(node_index(j)>0) then
                   !C dequeue node j
                   call dequeue(j, P_queue, node_index,queue_size)
                end if
                node_index(j)=-1
                aggregate_flags(j)=i
             END DO
             Temp(i, aggregate_number) = size_of_aggregates



             !C dequeue node i
             call dequeue(i, P_queue, node_index,queue_size)

             node_index(i) = -1
             aggregate_flags(i) = i
             !C----layer1 end


             !C----layer2 and 3 
             call dequeue_scoreup_Layer23(i, N, queue_size, P_queue,  & 
                  & node_index, aggregate_flags, node_record, NI, PNI, IAL, IAU)

          end DO
       end if
       deallocate(node_record,P_queue,aggregate_flags)
    end if


    !C== sendrecv independent aggregates
    CALL MPI_ALLGATHER(size_of_aggregates,1,LIS_MPI_INTEGER, &
         & GIN_aggregate(1),1,LIS_MPI_INTEGER,SOLVER_COMM,ierr)
    GIN_aggregate(0)=0
    DO i=2,NPROCS
       GIN_aggregate(i)=GIN_aggregate(i-1)+GIN_aggregate(i)
    END DO
    

    aggregate_number_in_table=0
    maxbufsize = STACK_IMPORT(NEIBPETOT)+STACK_EXPORT(NEIBPETOT)
    if(PE_num_size /= 0) maxbufsize=maxbufsize+in_nodes_for_each_PE(PE_num_size)+1
    
    allocate(req1(PE_list_size))
    allocate(sta1(MPI_STATUS_SIZE,PE_list_size))
    allocate(req2(PE_list_size))
    allocate(sta2(MPI_STATUS_SIZE,PE_list_size))

    allocate(send_buf(maxbufsize))
    allocate(recv_buf(maxbufsize))
    allocate(in_send_buf_for_PEs(0:PE_list_size))
    in_send_buf_for_PEs(0)=0
    position=0
    DO i=1,PE_list_size
       Neib_index=PE_list(i,2)
       PE_num_index=PE_list(i,3)
       if(Neib_index>0) then
          istart= STACK_EXPORT(Neib_index-1)
          inum  = STACK_EXPORT(Neib_index)-istart
          DO j=istart+1,istart+inum
             k=Temp(NOD_EXPORT(j)+ZERO_ORIGIN,aggregate_number)
             if(k>0) k=k+GIN_aggregate(my_rank)
             position=position+1
             send_buf(position)=k
          END DO
          istart= STACK_IMPORT(Neib_index-1)
          inum  = STACK_IMPORT(Neib_index)-istart
          DO j=istart+1,istart+inum
             k=Temp(NOD_IMPORT(j)+ZERO_ORIGIN,aggregate_number)
             if(k>0) k=k+GIN_aggregate(my_rank)
             position=position+1
             send_buf(position)=k
          END DO
       end if
       if(PE_num_index>0) then
          istart= in_nodes_for_each_PE(PE_num_index-1)
          inum = in_nodes_for_each_PE(PE_num_index)-istart
          DO j=istart+1,istart+inum
             k=Temp(nodes_for_each_PE(j)+ZERO_ORIGIN,aggregate_number)
             if(k>0) k=k+GIN_aggregate(my_rank)
             position=position+1
             send_buf(position)=k
          END DO
       end if
       in_send_buf_for_PEs(i)=position
    end DO
    DO i=1,PE_list_size
       dest=PE_list(i,1)
       istart=in_send_buf_for_PEs(i-1)
       inum=in_send_buf_for_PEs(i)-in_send_buf_for_PEs(i-1)
       CALL MPI_ISEND(send_buf(istart+1),inum,LIS_MPI_INTEGER, &
            & dest,0,SOLVER_COMM,req1(i),ierr)

       CALL MPI_IRECV(recv_buf(istart+1),inum,LIS_MPI_INTEGER, &
            & dest,0,SOLVER_COMM,req2(i),ierr)
    END DO

    CALL MPI_WAITALL(PE_list_size, req2, sta2, ierr)
    CALL MPI_WAITALL(PE_list_size, req1, sta1, ierr)
    deallocate(send_buf)
    deallocate(in_send_buf_for_PEs)

    position = 0
    DO i = 1, PE_list_size
       dest = PE_list(i, 1)
       Neib_index = PE_list(i, 2)
       PE_num_index = PE_list(i, 3)
       if(Neib_index > 0) then
          istart = STACK_IMPORT(Neib_index - 1)
          inum   = STACK_IMPORT(Neib_index) - STACK_IMPORT(Neib_index - 1)
          DO j = istart + 1, istart + inum
             position = position + 1

             !C== assigning part
             k = recv_buf(position)
             if(k > 0) then
                CALL assign_aggregate_table(aggregate_number_in_table,     &
                     & aggregate_table_array, aggregate_table_size, k, lc, &
                     & i, NPROCS, size_of_aggregates)

                Temp(NOD_IMPORT(j)+ZERO_ORIGIN, aggregate_number) = lc
                node_index(NOD_IMPORT(j)+ZERO_ORIGIN) = -1
             end if
             !C== end:assigning part             
          END DO
          istart = STACK_EXPORT(Neib_index - 1)
          inum   = STACK_EXPORT(Neib_index) - STACK_EXPORT(Neib_index - 1)
          DO j = istart + 1, istart + inum
             position =position + 1
             
             !C== assigning part
             k = recv_buf(position)
             if(k > 0) then
                CALL assign_aggregate_table(aggregate_number_in_table,     &
                     & aggregate_table_array, aggregate_table_size, k, lc, &
                     & i, NPROCS, size_of_aggregates)

                Temp(NOD_EXPORT(j)+ZERO_ORIGIN, aggregate_number) = lc
                node_index(NOD_EXPORT(j)+ZERO_ORIGIN) = -1
             end if
             !C== end:assigning part             
          END DO
       end if
       if(PE_num_index > 0) then
          istart = in_nodes_for_each_PE(PE_num_index - 1)
          inum = in_nodes_for_each_PE(PE_num_index) - istart
          DO j = istart + 1, istart + inum
             position = position + 1

             !C== assigning part
             k = recv_buf(position)
             if(k > 0) then
                CALL assign_aggregate_table(aggregate_number_in_table,     &
                     & aggregate_table_array, aggregate_table_size, k, lc, &
                     & i, NPROCS, size_of_aggregates)

                Temp(nodes_for_each_PE(j)+ZERO_ORIGIN, aggregate_number) = lc
                node_index(nodes_for_each_PE(j)+ZERO_ORIGIN) = -1
             end if
             !C== end:assigning part             

          END DO
       end if
    END DO

    deallocate(recv_buf)
    deallocate(req1)
    deallocate(req2)
    deallocate(sta1)
    deallocate(sta2)

    !C== end: sendrecv independent aggregates

    !C== left nodes are associated to independent aggregates
    DO i = 1, N
       if(node_index(i) == 0) then
          aggr = -1
          k = PNI(2 * i - 2) + 1
          DO while(aggr <= 0 .and. k <= PNI(2 * i - 1))
             aggr = Temp(IAL( NI(k) )+ZERO_ORIGIN, aggregate_number)
             k = k + 1
          END DO
          if(aggr <= 0) then
             k = PNI(2 * i - 1) + 1
             DO while( aggr <= 0 .and. k <= PNI(2 * i))
                aggr = Temp(IAU( NI(k) )+ZERO_ORIGIN, aggregate_number)
                k = k + 1
             END DO
          end if
          !C in case of aggr<0, neiboring node_index is -1 for some reason.
          if(aggr > 0) Temp(i, aggregate_number) = aggr
       end if
    END DO
    !C==end: left nodes are associated to independent aggregates
#ifdef PRINT_DATA_CREATE    
    write(*,*) GIN_aggregate(my_rank + 1) - GIN_aggregate(my_rank),                  &
         & size_of_aggregates - GIN_aggregate(my_rank + 1) + GIN_aggregate(my_rank), &
         & N, my_rank
#endif    
    !C== Temp is written to the aggregates_results.
    !C   size_of_aggregates is the number of aggregates in 1:NP
    in_aggregates_result_size = size_of_aggregates
    ALLOCATE(in_aggregates_result(0:size_of_aggregates))
    ALLOCATE(aggregates_result(N))

    DO i=1, N
       aggr=Temp(i,aggregate_number)
       if(aggr>0) Temp(aggr,aggregate_index) = Temp(aggr,aggregate_index)+1
    END DO
    in_aggregates_result(0)=0    
    in_aggregates_result_size=size_of_aggregates
    DO i=1,size_of_aggregates
       Temp(i,aggregate_index)=Temp(i,aggregate_index)+Temp(i-1,aggregate_index)
       in_aggregates_result(i)=Temp(i-1,aggregate_index)
    END DO

    DO i=1,N
       aggr=Temp(i,aggregate_number)
       if(aggr>0) then
          k=in_aggregates_result(aggr)+1
          in_aggregates_result(aggr)=k
          aggregates_result(k)=i
       end if
    END DO
    
    DEALLOCATE(Temp)
    !C== end: Temp is written to the aggregates_results.

  END SUBROUTINE aggregate

  !C-- insert node into aggregate table
  !C-- global_nd -> aggregate_table(:,neib) and local aggregate number is in local_nd
  SUBROUTINE assign_aggregate_table(aggregate_number_in_table, aggregate_table_array, &
       & aggregate_table_size, global_nd, local_nd, neib, NPROCS, size_of_aggregates)
    implicit NONE
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    INTEGER(kind=kint), intent(in)    :: global_nd, NPROCS, aggregate_table_size, neib
    INTEGER(kind=kint), intent(inout) :: local_nd
!!$    INTEGER(kind=kint), intent(inout) :: aggregate_number_in_table(NPROCS), size_of_aggregates
!!$    INTEGER(kind=kint), intent(inout) :: aggregate_table_array(aggregate_table_size * 2, NPROCS)
    INTEGER(kind=kint), intent(inout) :: aggregate_number_in_table(:), size_of_aggregates
    INTEGER(kind=kint), intent(inout) :: aggregate_table_array(:, :)
    
    INTEGER(kind=kint) :: m, l

    m = aggregate_number_in_table(neib)
    DO l = 1, m
       if(aggregate_table_array(2 * l - 1, neib) == global_nd) exit
    END DO
    if (l > m) then
       size_of_aggregates = size_of_aggregates + 1

#ifdef ALLOCATION_CHECK       
       if(aggregate_table_size < l) then
          write(*,*) "enlarge aggregate_table_size"
          stop 'error in assign_aggregate_table():'
       end if
#endif
       local_nd = size_of_aggregates
       
       aggregate_table_array(2 * l - 1, neib) = global_nd
       aggregate_table_array(2 * l, neib) = local_nd
       aggregate_number_in_table(neib) = l

    else

       local_nd = aggregate_table_array(2 * l, neib)
       
    end if
  END SUBROUTINE assign_aggregate_table
END MODULE aggregate_mod
