C   Copyright (C) 2005 The Scalable Software Infrastructure Project. 
C   All rights reserved.
C
C   Redistribution and use in source and binary forms, with or without
C   modification, are permitted provided that the following conditions
C   are met:
C   1. Redistributions of source code must retain the above copyright
C      notice, this list of conditions and the following disclaimer.
C   2. Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C   3. Neither the name of the project nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C
C   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE
C   PROJECT ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
C   BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
C   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
C   THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT BE LIABLE FOR ANY
C   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
C   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
C   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
C   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
C   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
C   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

      implicit none
      
#include "lisf.h"

      LIS_COMPLEX z
      LIS_VECTOR v
      LIS_INTEGER n,ln,comm_world,ierr
      LIS_REAL nrm2

      n = 10
      ln = 0
      comm_world = LIS_COMM_WORLD      

      z = (1.0, 2.0);

      write(*,*) 'complex number z = ', real(z), '+', imag(z), ' * I'

#ifdef COMPLEX      
      call lis_initialize(ierr)
      call lis_vector_create(comm_world,v,ierr)
      call lis_vector_set_size(v,ln,n,ierr)
      call lis_vector_set_all(z,v,ierr)
      call lis_vector_print(v,ierr)
      call lis_vector_nrm2(v,nrm2,ierr)
      write(*,*) '2-norm of v = ', nrm2
      write(*,*) 'abs(z) = ', abs(z)      
      call lis_vector_destroy(v,ierr)
      call lis_finalize(ierr)
#endif      

      stop
      end
      
