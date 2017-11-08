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
#include <math.h>
#include "lis.h"

#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(int argc, char* argv[])
{
	LIS_Comm comm;
	LIS_SCALAR *a,*x,*b,*u,*w;
	LIS_INT m,n,nn;
	LIS_INT	i,j,ii,jj,nnz;
	double time,time0;
	LIS_REAL resid_r,resid_b;	

	LIS_DEBUG_FUNC_IN;

	lis_initialize(&argc, &argv);

	comm = LIS_COMM_WORLD;

	if( argc < 3 )
	{
	  lis_printf(comm,"Usage: %s m n\n", argv[0]);
	  CHKERR(1);	  	    	    	  	  
	}

	m  = atoi(argv[1]);
	n  = atoi(argv[2]);
	if( m<=0 || n<=0 )
	{
	  lis_printf(comm,"m=%D <=0 or n=%D <=0\n", m,n);
	  CHKERR(1);
	}
	
	lis_printf(comm,"\n");

	/* create arrays */

	nn = m*n;

	a = (LIS_SCALAR *)malloc(nn*nn*sizeof(LIS_SCALAR));
	b = (LIS_SCALAR *)malloc(nn*sizeof(LIS_SCALAR));
	x = (LIS_SCALAR *)malloc(nn*sizeof(LIS_SCALAR));
	u = (LIS_SCALAR *)malloc(nn*sizeof(LIS_SCALAR));
	w = (LIS_SCALAR *)malloc(nn*nn*sizeof(LIS_SCALAR));

	/* define two-dimensional Laplacian */

	lis_array_set_all(nn*nn,(LIS_SCALAR)0.0,a);

	nnz = 0;
	for(ii=0;ii<nn;ii++)
	  {
	    i = ii/m;
	    j = ii - i*m;
	    if( i>0 )   { jj = ii - m; a[ii + jj * nn] = -1.0; nnz++;} 
	    if( i<n-1 ) { jj = ii + m; a[ii + jj * nn] = -1.0; nnz++;}
	    if( j>0 )   { jj = ii - 1; a[ii + jj * nn] = -1.0; nnz++;}
	    if( j<m-1 ) { jj = ii + 1; a[ii + jj * nn] = -1.0; nnz++;}
	    jj = ii; a[ii + jj * nn] = 4.0; nnz++;
	  }

	lis_printf(comm,"matrix size = %D x %D (%D nonzero entries)\n\n", nn,nn,nnz);

	lis_array_set_all(nn,(LIS_SCALAR)1.0,u);
 	lis_array_matvec(nn,a,u,b,LIS_INS_VALUE);

	/* solve linear system */	

	time0 = lis_wtime();
	lis_array_solve(nn,a,b,x,w);
	time = lis_wtime() - time0;

	lis_array_xpay(nn,x,(LIS_SCALAR)-1.0,u);	
	lis_array_nrm2(nn,u,&resid_r);
	lis_array_nrm2(nn,b,&resid_b);

	/*
	for (i=0;i<nn;i++)
	  {
#ifdef _COMPLEX
	    lis_printf(comm,"x(%D) = (%e, %e)\n", i, (double)creal(x[i]), (double)cimag(x[i]));
#else
	    lis_printf(comm,"x(%D) = %e\n", i, (double)x[i]);	    
#endif	    
	  }
	lis_printf(comm,"\n");
	*/

	lis_printf(comm,"Direct: elapsed time         = %e sec.\n", time);
	lis_printf(comm,"Direct:   linear solver      = %e sec.\n", time);
	lis_printf(comm,"Direct: relative residual    = %e\n\n", (double)(resid_r/resid_b));

	free(a);
	free(b);
	free(x);
	free(u);
	free(w);

	lis_finalize();

	LIS_DEBUG_FUNC_OUT;

	return 0;
}


