
/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
        #include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
        #include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lis.h"


#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(LIS_INT argc, char* argv[])
{
	LIS_COMPLEX z,dot;
	LIS_VECTOR v;
	LIS_INT n,ln;
	LIS_REAL nrm2;

#ifdef _Complex_I
	z = 1.0 + 2.0 * _Complex_I;
#else
	z[0] = 1.0;
	z[1] = 2.0;
#endif
	
#ifdef HAVE_COMPLEX_H
#ifdef _LONG__DOUBLE
	printf("complex number z = (%Lf, %Lf)\n", creall(z), cimagl(z));
#else
	printf("complex number z = (%f, %f)\n", creal(z), cimag(z));
#endif
#else	
#ifdef _Complex_I
#ifdef _LONG__DOUBLE
	printf("complex number z = (%Lf, %Lf)\n", __real__ z, __imag__ z);
#else
	printf("complex number z = (%f, %f)\n", __real__ z, __imag__ z);
#endif
#else
#ifdef _LONG__DOUBLE
	printf("complex number z = (%Lf, %Lf)\n", z[0], z[1]);
#else
	printf("complex number z = (%f, %f)\n", z[0], z[1]);
#endif
#endif
#endif

#ifdef _COMPLEX	
	lis_initialize(&argc, &argv);
	n = 10;
	ln = 0;
	lis_vector_create(LIS_COMM_WORLD,&v);
	lis_vector_set_size(v,ln,n);
	lis_vector_set_all(z,v);
	printf("complex vector v =\n");	
	lis_vector_print(v);
	lis_vector_conjugate(v);
	printf("conj(v) =\n");		
	lis_vector_print(v);	
	lis_vector_dot(v,v,&dot);
	lis_vector_nrm2(v,&nrm2);
#ifdef _LONG__DOUBLE
	printf("inner product (v,v) = (%Lf, %Lf)\n", creall(dot), cimagl(dot));
	printf("2-norm of v = %Lf\n", nrm2);
	printf("abs(z) = %Lf\n", fabs(z));
#else
	printf("inner product (v,v) = (%f, %f)\n", creal(dot), cimag(dot));
	printf("2-norm of v = %f\n", nrm2);
	printf("abs(z) = %f\n", fabs(z));	
#endif
	lis_vector_destroy(v);
	lis_finalize();
#endif

	return 0;
}


