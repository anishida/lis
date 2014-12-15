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

!C   ************************************************
!C   * SUBROUTINE finit_data_creation
!C   * SUBROUTINE finit_data_creation_unsym
!C   * SUBROUTINE finit_v_cycle
!C   * SUBROUTINE finit_clear_matrix
!C   ************************************************

SUBROUTINE finit_data_creation(cinit_data_creation)
  use data_creation_AMGCG
!  extern cinit
  call cinit_data_creation(data_creation)
END SUBROUTINE finit_data_creation

SUBROUTINE finit_data_creation_unsym(cinit_data_creation_unsym)
  use data_creation_AMGCG
!  extern cinit
  call cinit_data_creation_unsym(data_creation_unsym)
END SUBROUTINE finit_data_creation_unsym

SUBROUTINE finit_v_cycle(cinit_v_cycle)
  use solver_AMGCG
!  extern cinit
  call cinit_v_cycle(v_cycle)
END SUBROUTINE finit_v_cycle

SUBROUTINE finit_clear_matrix(cinit_clear_matrix)
  use solver_AMGCG
!  extern cinit
  call cinit_clear_matrix(clear_matrix)
END SUBROUTINE finit_clear_matrix
