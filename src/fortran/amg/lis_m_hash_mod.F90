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
!C   * MODULE hash_mod
!C     CONTAINS
!C   * SUBROUTINE hash_zero
!C   * SUBROUTINE hash
!C   ************************************************

MODULE hash_mod
#ifdef DEBUG
  INTEGER(kind=4) :: hash_conflict_count = 0 
#endif
    
  contains
    SUBROUTINE hash_zero(hash_size, hash_array, val, position, id)
    implicit none
    
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    INTEGER(kind=kint), intent(in)  :: hash_size, val
    INTEGER(kind=kint), intent(in)  :: hash_array(1:hash_size)
    CHARACTER(len=*), intent(in) :: id
    INTEGER(kind=kint), intent(out) :: position

    INTEGER(kind=kint) :: i

    position = MOD(val, hash_size) + 1
    if(hash_array(position) /= 0) then
       position = MOD(val * MOD(val, 1001), hash_size) + 1
       if(hash_array(position) /= 0) then
#ifdef DEBUG          
          hash_conflict_count = hash_conflict_count + 1
#endif
          do i = position + 1, hash_size
             if(hash_array(i) == 0) exit
          end do
          if(i > hash_size)then
             do i = 1, position - 1
                if(hash_array(i) == 0) exit
             end do

#ifdef ALLOCATION_CHECK
             if(i >= position) then
                write(*,*) "enlarge ", id
                stop "error in hash_zero():"
             end if
#endif
          end if
          position = i
       end if
    end if
  END SUBROUTINE hash_zero

  SUBROUTINE hash(hash_size, hash_array, val, position, id)
    implicit none
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    
    INTEGER(kind=kint), intent(in)  :: hash_size, val, hash_array(:)
    CHARACTER(len=*), intent(in) :: id
    INTEGER(kind=kint), intent(out) :: position

    INTEGER(kind=4) :: i
    
    position = MOD(val, hash_size) + 1
    if(hash_array(position) /= 0 .and. hash_array(position) /= val) then
       position = MOD(val*MOD(val, 1001), hash_size) + 1

       if(hash_array(position) /= 0 .and. hash_array(position) /= val) then
#ifdef DEBUG          
          hash_conflict_count = hash_conflict_count + 1
#endif
          do i = position + 1, hash_size
             if(hash_array(i) == 0 .or. hash_array(i) == val) exit
          end do

          if(i > hash_size) then
             do i = 1, position - 1
                if(hash_array(i) == 0 .or. hash_array(i) == val) exit
             end do
#ifdef ALLOCATION_CHECK
             if(i >= position) then
                write(*,*) "enlarge ", id
                stop "error in hash_zero():"
             end if
#endif
          end if
          position = i
       end if
    end if
  END SUBROUTINE hash
  
END MODULE hash_mod
