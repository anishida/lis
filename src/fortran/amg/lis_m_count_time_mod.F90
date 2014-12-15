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
!C   * MODULE count_time_mod
!C     CONTAINS
!C   * SUBROUTINE count_time
!C   ************************************************

MODULE count_time_mod
contains
  SUBROUTINE count_time(op, SOLVER_COMM, my_rank, time_kind)
    IMPLICIT NONE
#ifdef LONGLONG
    include 'precision_longlong.inc'
#else
    include 'precision.inc'
#endif
    INCLUDE 'mpif.h'

    INTEGER(kind=kint), INTENT(in) :: op, solver_comm, my_rank, time_kind
    INTEGER(kind=kint), parameter :: TIME_KINDS = 20
    INTEGER(kind=kint) :: ierr, NPROCS, i, num_of_kinds
    REAL   (kind=kreal), SAVE :: time_table(TIME_KINDS), start_table(TIME_KINDS)
    REAL   (kind=kreal) :: recv_time, send_time, out_time
    LOGICAL :: prt_ech_PE(TIME_KINDS)
    DATA time_table/TIME_KINDS*0/
    DATA start_table/TIME_KINDS*0/

    !C num_of_kinds must be less than or equal to TIME_KINDS
    num_of_kinds = 15

    !C prt_ech_PE records time_kind at which each PE's time is needed.
    prt_ech_PE = .FALSE.

    
    !C op:1 start 2 stop 3 print time 0 reset time_table
    SELECT CASE(op)
    CASE(0)
       time_table(time_kind) = 0
       
    CASE(1)
       start_table(time_kind) = MPI_WTIME()

    CASE(2)
       IF(start_table(time_kind) /= 0) THEN
          time_table(time_kind) = time_table(time_kind) + &
               &                  MPI_WTIME() - start_table(time_kind)
          start_table(time_kind) = 0
       ELSE
          WRITE(*,*) "timer error...:", time_kind, start_table(time_kind)
          STOP "timer error!!"
       END IF
    CASE(3)
       send_time = time_table(time_kind)
       recv_time = 0
       CALL MPI_REDUCE(send_time, recv_time, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, SOLVER_COMM, ierr)
       CALL MPI_COMM_SIZE(SOLVER_COMM, NPROCS, ierr)
       IF(my_rank == 0) THEN
          out_time = recv_time / NPROCS
          WRITE(*,'("***", i3, ":", 1pe16.6, " sec")') time_kind, out_time
       END IF

    CASE(4)
       CALL MPI_COMM_SIZE(SOLVER_COMM, NPROCS, ierr)
       IF(my_rank == 0) THEN
          WRITE(*,*) "# 1:sr_rep 2:sr_aggre 3:sr_elem 4:RAP 5:neib 6:agr 7:sm_agr #"
          WRITE(*,*) "# 8:sgs    9:v_cycle 10:calculation part of sgs             #"
          WRITE(*,*) "# 11:matrix transform 12:ParMetis 13:matrix redistribution  #"
          WRITE(*,*) "# 14:Total time in matrix dist 15:calculation in Redist()   #"
       END IF
       DO i = 1, num_of_kinds
          send_time = time_table(i)
          recv_time = 0

          CALL MPI_barrier(solver_comm, ierr)
          IF(prt_ech_PE(i)) WRITE(*,'( i3,"--#", i3, ":", 1pe16.6, " sec")') my_rank, i, send_time
          CALL MPI_REDUCE(send_time, recv_time, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, SOLVER_COMM, ierr)

          IF(my_rank == 0) THEN
             out_time = recv_time / NPROCS
             WRITE(*,'("***", i3, ":", 1pe16.6, " sec")') i, out_time
          END IF
       END DO
    CASE DEFAULT
       STOP 'error start_time stop time'
    END SELECT
  END SUBROUTINE count_time
END MODULE count_time_mod
