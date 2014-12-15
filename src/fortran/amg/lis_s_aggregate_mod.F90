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
!C   * MODULE aggregate_mod
!C     CONTAINS
!C   * SOUBROUTINE indpdt_agrgt
!C   ************************************************

MODULE aggregate_mod
CONTAINS
SUBROUTINE indpdt_agrgt (PNI,NI,LEVEL_NO,node_index,in_aggregates_result, &
     & aggregates_result,in_aggregates_result_size)
  USE data_structure_for_AMG
    USE queue_mod
    IMPLICIT NONE

    INTEGER(kind=kint), DIMENSION(: ), INTENT(in) :: NI
    INTEGER(kind=kint), DIMENSION(0:), INTENT(in) :: PNI
    INTEGER(kind=kint), INTENT(in)                :: LEVEL_NO
    INTEGER(kind=kint), DIMENSION(:), INTENT(inout) :: node_index
    INTEGER(kind=kint), DIMENSION(:),pointer :: in_aggregates_result
    INTEGER(kind=kint), DIMENSION(:),pointer :: aggregates_result
    INTEGER(kind=kint), INTENT(out) :: in_aggregates_result_size


    !C    TYPE(INTER_LEVEL_OPERATOR),POINTER :: P
    INTEGER(kind=kint), DIMENSION(:,:), ALLOCATABLE :: Temp
    
    INTEGER(kind=kint), DIMENSION(:,:), ALLOCATABLE :: P_queue
    INTEGER(kind=kint), DIMENSION(:),   ALLOCATABLE :: aggregate_flags

    INTEGER(kind=kint) :: head,tail,queue_size

    INTEGER(kind=kint ) :: N,NP
    REAL   (kind=kreal),POINTER, DIMENSION(:)  ::  D
    REAL   (kind=kreal),POINTER, DIMENSION(:)  ::  AU
    REAL   (kind=kreal),POINTER, DIMENSION(:)  ::  AL

    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  INU
    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  IAU
    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  INL
    INTEGER(kind=kint ),POINTER, DIMENSION(:) ::  IAL
    
    INTEGER(kind=kint ),allocatable, DIMENSION(:) :: node_record

    INTEGER(kind=kint) :: i,s,s_u,k,f,aggr,i_CN,size_of_aggregates,count,l,j,m,front,back
    INTEGER(kind=kint) :: node_record_counter,node,score,insert_point,max_score

    REAL   (kind=kreal) :: r
    INTEGER(kind=kint),parameter :: aggregate_number=1
    INTEGER(kind=kint),parameter :: aggregate_index=2
    logical :: flag

    s  = HIERARCHICAL_DATA(LEVEL_NO-1) % N 
    NP = HIERARCHICAL_DATA(LEVEL_NO-1) % NP 

    ALLOCATE(Temp(0:s,2))
    ALLOCATE(P_queue(2,s))
    ALLOCATE(aggregate_flags(s))
    ALLOCATE(node_record(s))
    
    D  =>HIERARCHICAL_DATA(LEVEL_NO-1) % D 
    INL=>HIERARCHICAL_DATA(LEVEL_NO-1) % INL 
    INU=>HIERARCHICAL_DATA(LEVEL_NO-1) % INU 
    IAL=>HIERARCHICAL_DATA(LEVEL_NO-1) % IAL
    IAU=>HIERARCHICAL_DATA(LEVEL_NO-1) % IAU
    AL =>HIERARCHICAL_DATA(LEVEL_NO-1) % AL
    AU =>HIERARCHICAL_DATA(LEVEL_NO-1) % AU

    s_u=INU(s)

    !C aggregates is stored in P temporarily 
    !C    P=>HIERARCHICAL_DATA(LEVEL_NO) % P    
    
    !C Temp =0
    Temp(1:s,aggregate_number) = 0
    Temp(0:s,aggregate_index) = 0
    

    DO i=1,s
       if(node_index(i)==0) then 
          node_index(i)=1
       end if
    END DO
    
    !C initialize Priority queue
    !C P_queue has score and node index. node_index has the index of P_queue
    !C -1 in node_index : selected as part of aggregate
    !C  0 in node_index : not selected as part of aggregate, but it's neighbor is aggregate.

    P_queue(1:2,1:s)=0
    aggregate_flags(1:s)=0
    
    
    !C interface part is scored, and enqueued
    !C Temp(,aggregate_number) is used as buffer
    queue_size=0
    
    !C enqueue nodes whose score is 1
    DO count=1,s
       if(Temp(count,aggregate_number)==1)then
          Temp(count,aggregate_number)=0
          !C this case is enqueued
       else if(node_index(count)==1)then
          queue_size=queue_size+1
          node_index(count)=queue_size
          P_queue(1,queue_size)=1
          P_queue(2,queue_size)=count
       end if
    END DO
    size_of_aggregates=0

    DO while(P_queue(1,1) > 0)
       
       i=P_queue(2,1)
       size_of_aggregates = size_of_aggregates+1

       !C----layer1 start
       DO k=PNI(i-1)+1,PNI(i)
          j=IAL( NI(k) )+ZERO_ORIGIN
          
          Temp(j, aggregate_number) = size_of_aggregates
          
          if(node_index(j)>0) then
             !C dequeue node j
             call dequeue(j, P_queue, node_index,queue_size)
          end if
          node_index(j)=-1
          aggregate_flags(j)=i
          
       END DO
       DO k=PNI(i+NP-1)+1,PNI(i+NP)
          if(IAU(NI(k))+ZERO_ORIGIN>s) cycle
          j=IAU( NI(k) )+ZERO_ORIGIN
          Temp(j, aggregate_number) = size_of_aggregates
          
          if(node_index(j)>0) then
             !C dequeue node j
             call dequeue(j, P_queue, node_index,queue_size)
          end if
          node_index(j)=-1
          aggregate_flags(j)=i
       END DO
       
       Temp(i , aggregate_number) = size_of_aggregates
       
       !C dequeue node i
       call dequeue(i, P_queue, node_index,queue_size)
       
       node_index(i) = -1
       aggregate_flags(i) = i
       !C----layer1 end

       !C----layer2 and 3 
       call dequeue_scoreup_Layer23(i, s, NP, queue_size, P_queue,  & 
            & node_index, aggregate_flags, node_record, NI, PNI, IAL, IAU)
    END DO


    deallocate( node_record, aggregate_flags )

    aggregate_check: DO count=1,s
       if(node_index(count)==0) then
          lower_part: DO k=PNI(count-1)+1,PNI(count)
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
          
          upper_part: DO k=PNI(count+NP-1)+1,PNI(count+NP)
             !C aggr is the aggregate number of next to i
             if(IAU(NI(k))+ZERO_ORIGIN>s) cycle
             aggr = Temp(IAU(NI(k))+ZERO_ORIGIN, aggregate_number)
             
             IF(aggr > 0) Then
                Temp(count,aggregate_number) = aggr
                exit upper_part
             END IF
          END DO upper_part
          IF(Temp(count, aggregate_number) == 0) Then
             
             write(*,*) "independent node in aggregate()",node_index(count),count
          END IF
       end if
    END DO aggregate_check

    DEALLOCATE(P_queue)
  
    in_aggregates_result_size=size_of_aggregates
    ALLOCATE( in_aggregates_result(0:size_of_aggregates),aggregates_result(s) )

    in_aggregates_result(0)=0

    DO i=1,s
       aggr=Temp(i,aggregate_number)
       if(aggr>0) Temp(aggr,aggregate_index) = Temp(aggr,aggregate_index)+1
    END DO

    !C size_of_aggregates has the size of aggregates
    !C P % IN(aggr)+1 points the place for new data of aggr
    DO i=1,size_of_aggregates
       Temp(i,aggregate_index) = Temp(i,aggregate_index)+Temp(i-1,aggregate_index)
       in_aggregates_result(i) = Temp(i-1,aggregate_index)
    END DO

    !C unneeded nodes belong to aggregate NO 0.
    DO i=1,s
       aggr = Temp(i, aggregate_number)
       if(aggr>0) then
          k=in_aggregates_result(aggr)+1
          in_aggregates_result(aggr)=k
          aggregates_result(k)=i
       end if
    END DO

!C    write(*,*) size_of_aggregates, in_aggregates_result(size_of_aggregates)
    
    DEALLOCATE(Temp)
    
  end subroutine indpdt_agrgt

  SUBROUTINE agrgt_dealloc ( in_aggregates_result,aggregates_result )
    implicit none
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif

    INTEGER(kind=kint), pointer :: in_aggregates_result(:)
    INTEGER(kind=kint), pointer :: aggregates_result(:)
    
    deallocate(in_aggregates_result, aggregates_result)
  END SUBROUTINE agrgt_dealloc
    

END MODULE aggregate_mod
