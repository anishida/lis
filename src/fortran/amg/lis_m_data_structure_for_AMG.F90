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
!C   * MODULE data_structure_for_AMG
!C   ************************************************

module data_structure_for_AMG
  implicit none
#ifdef LONGLONG
  include 'precision_longlong.inc'
#else
  include 'precision.inc'
#endif

  
  integer, parameter :: SWITCH_SIZE    = 1700
  integer, parameter :: MAX_LEVEL_SIZE = 10
  integer, parameter :: MIN_NODE_SIZE  = 90
  
  type row_node
     INTEGER(kind=kint) :: column
     REAL(kind=kreal)   :: value
  end type row_node
  
  type INTER_LEVEL_OPERATOR
     integer(kind=kint) :: ROW_SIZE
     integer(kind=kint),  pointer :: IN(:)
     integer(kind=kint),  pointer :: CN(:)
     real   (kind=kreal), pointer :: V(:)
  end type INTER_LEVEL_OPERATOR
  
  type COMMUNICATION_TABLE
     integer(kind=kint) :: NEIBPETOT
     integer(kind=kint), pointer :: NEIBPE(:)
     integer(kind=kint), pointer :: STACK_IMPORT(:)
     integer(kind=kint), pointer :: STACK_EXPORT(:)
     integer(kind=kint), pointer :: NOD_IMPORT(:)
     integer(kind=kint), pointer :: NOD_EXPORT(:)
  end type COMMUNICATION_TABLE
  
  type DATA_FOR_EACH_LEV
     integer(kind=kint ) :: N,NP,NPL,NPU
     integer(kind=kint ), pointer :: INU(:),INL(:)
     integer(kind=kint ), pointer :: IAU(:),IAL(:)
     real   (kind=kreal), pointer :: AU(:),AL(:),D(:),X(:),B(:)
     
     type(COMMUNICATION_TABLE) :: COMM_TABLE
     type(COMMUNICATION_TABLE) :: INT_LVL_TABLE
     type(INTER_LEVEL_OPERATOR), pointer :: R, P

  end type DATA_FOR_EACH_LEV
  
  type (DATA_FOR_EACH_LEV) :: HIERARCHICAL_DATA(MAX_LEVEL_SIZE)
  
!C Rate_Of_Space * mesh size is the width of the next coarser 
!C level's restriction

  real (kind=kreal),  parameter :: EPS = 1.E-10
  real (kind=kreal),  parameter :: EPS_COARSEST = 1.E-5
  INTEGER(kind=kint), parameter :: NODE_RECORD_SIZE=300
  LOGICAL :: SYMMETRIC_FLAG
  
  INTEGER*8 :: pspcommf
  INTEGER(kind=kint),pointer :: rowdistbx(:)
  INTEGER(kind=kint) :: ldb,ldx,nrhs
  INTEGER(kind=kint) :: ioptions(16)
  REAL(kind=kreal)   :: doptions(16)
  REAL(kind=kreal),pointer        :: coarsest_X(:,:),coarsest_B(:,:)

end module data_structure_for_AMG






