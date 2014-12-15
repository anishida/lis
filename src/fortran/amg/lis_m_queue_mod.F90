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
!C   * MODULE queue_mod
!C     CONTAINS
!C   * SUBROUTINE score_up
!C   * SUBROUTINE dequeue
!C   * SUBROUTINE dequeue_scoreup_Layer23
!C   ************************************************

MODULE queue_mod
contains
  SUBROUTINE score_up(node, P_queue, node_index, n)
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    INTEGER(kind=kint), intent(in)    :: node, n
    INTEGER(kind=kint), intent(inout) :: P_queue(:,:)
    !C-DIMENSION(2,n)
    INTEGER(kind=kint), intent(inout) :: node_index(:)
    
    INTEGER(kind=kint) :: i, w1, w2, p

    i = node_index(node)

    P_queue(1, i) = P_queue(1, i) + 1

    P = i / 2
    DO while(i > 1)
       if(P_queue(1, p) >= P_queue(1, i)) exit
       w1 = P_queue(1, p)
       P_queue(1, p) = P_queue(1, i)
       P_queue(1, i) = w1

       w2 = P_queue(2, p)
       w1 = P_queue(2, i)
       P_queue(2, p) = w1
       P_queue(2, i) = w2

       !C node_index(P_queue(2, i)) = i
       node_index(w2) = i
       !C node_index(P_queue(2, p)) = p
       node_index(w1) = p        

       i = p
       p = i / 2
    END DO
  END SUBROUTINE score_up

  SUBROUTINE dequeue(node, P_queue, node_index, n)
    IMPLICIT NONE
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif

    INTEGER(kind=kint), intent(in) :: node, n
    INTEGER(kind=kint), intent(inout) :: P_queue(:,:)
    !C- DIMENSION(2,n)
    INTEGER(kind=kint), intent(inout) :: node_index(:)

    INTEGER(kind=kint) :: i, w, c
    !C dequeue i
    i = node_index(node)



    P_queue(1, i) = -1
    P_queue(2, i) = -1
    c = i * 2
    DO while(c <= n)
       if(c + 1 <= n) then
          if(P_queue(1, c + 1) > P_queue(1, c)) c = c + 1
       end if
       if(P_queue(1, c) == -1) exit



       !C swap
       w = P_queue(2, c)
#ifdef DEBUG
     if(w == 0) write(*,*) i, P_queue(1, i), P_queue(2, i), c, P_queue(1, c), P_queue(2, c)
#endif       
       P_queue(1, i) = P_queue(1, c)
       P_queue(2, i) = w
       node_index(w) = i

       P_queue(1, c) = -1
       P_queue(2, c) = -1

       
       i = c
       c = 2 * i
    end DO
  END SUBROUTINE dequeue

  SUBROUTINE dequeue_scoreup_Layer23(i, N, queue_size, P_queue, node_index, & 
       &  aggregate_flags, node_record, NI, PNI, IAL, IAU)
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    INTEGER(kind=kint), intent(in)    :: i, N, queue_size
    INTEGER(kind=kint), intent(inout) :: P_queue(:,:)
    !C- DIMENSION(2,N)
    INTEGER(kind=kint), intent(inout) :: node_index(:)
    INTEGER(kind=kint), intent(inout) :: aggregate_flags(:)
    !C- DIMENSION(1:N)
    INTEGER(kind=kint), intent(inout) :: node_record(:)
    !C- DIMENSION(1:N)
    INTEGER(kind=kint), intent(in)    :: NI(:)
    INTEGER(kind=kint), intent(in)    :: PNI(0:)
    INTEGER(kind=kint), intent(in)    :: IAL(:),IAU(:)

    INTEGER(kind=kint) :: j,k,l,node_record_counter,node

    node_record_counter = 0



    !C----layer2 start
    DO k = PNI(2 * i - 2) + 1, PNI(2 * i - 1)
       !C layer1: j
       j=IAL( NI(k) ) + ZERO_ORIGIN
       DO l = PNI(2 * j - 2) + 1, PNI(2 * j - 1)
          !C layer2: node
          node = IAL(NI(l)) + ZERO_ORIGIN
          if(aggregate_flags(node) /= i) then
             !C record the node
             aggregate_flags(node) = i
             node_record_counter = node_record_counter + 1
             node_record(node_record_counter) = node
             if(node_index(node) > 0) then
                !C dequeue node
                call dequeue(node, P_queue, node_index, queue_size)
                node_index(node) = 0
             end if
          end if
       END DO
       DO l = PNI(2 * j - 1) + 1, PNI(2 * j)
          !C layer2: node
          if(IAU(NI(l))+ZERO_ORIGIN > N) cycle
          node = IAU(NI(l))+ZERO_ORIGIN
          if(aggregate_flags(node) /= i)then
             !C record the node
             aggregate_flags(node) = i
             node_record_counter = node_record_counter + 1
             node_record(node_record_counter) = node
             if(node_index(node) > 0) then
                !C dequeue node
                call dequeue(node, P_queue, node_index, queue_size)
                node_index(node) = 0
             end if
          end if
       END DO
    ENDDO


    DO k=PNI(2*i-1)+1,PNI(2*i)
       !C layer1: j
       if(IAU(NI(k))+ZERO_ORIGIN > N) cycle
       j = IAU( NI(k) )+ZERO_ORIGIN
       DO l = PNI(2*j-2)+1, PNI(2*j-1)
          !C layer2: node
          node = IAL(NI(l)) + ZERO_ORIGIN
          if(aggregate_flags(node) /= i)then
             !C record the node
             aggregate_flags(node) = i
             node_record_counter = node_record_counter+1
             node_record(node_record_counter) = node
             if(node_index(node) > 0) then
                !C dequeue node

                call dequeue(node, P_queue, node_index, queue_size)
                node_index(node)=0

                
             end if
          end if
       END DO
       DO l=PNI(2*j-1)+1,PNI(2*j)
          !C layer2: node
          if(IAU(NI(l))+ZERO_ORIGIN>N) cycle
          node=IAU(NI(l))+ZERO_ORIGIN
          if(aggregate_flags(node)/= i)then
             !C record the node
             aggregate_flags(node)=i
             node_record_counter=node_record_counter+1
             node_record(node_record_counter)=node
             if(node_index(node)>0) then
                !C dequeue node

                call dequeue(node, P_queue, node_index,queue_size)
                node_index(node)=0



             end if
          end if
       END DO
    ENDDO
    !C----layer2 end



    !C----layer3 start
    DO k=1,node_record_counter
       !C layer2: j
       j=node_record(k)
       DO l=PNI(2*j-2)+1,PNI(2*j-1)
          !C layer3: node
          node=IAL(NI(l)) + ZERO_ORIGIN
          if(aggregate_flags(node)/= i)then
             !C record the node
             aggregate_flags(node)=i
             if(node_index(node)>0) then
                call score_up(node, P_queue, node_index, queue_size)
             end if
          end if
       END DO
       DO l=PNI(2*j-1)+1,PNI(2*j)
          !C layer3: node
          if(IAU(NI(l))+ZERO_ORIGIN>N) cycle
          node=IAU(NI(l))+ZERO_ORIGIN
          if(aggregate_flags(node)/= i)then
             !C record the node
             aggregate_flags(node)=i
             if(node_index(node)>0) then
                call score_up(node, P_queue, node_index, queue_size)
             end if
          end if
       END DO
    END DO
    !C----layer3 end


  END SUBROUTINE dequeue_scoreup_Layer23
end MODULE queue_mod
