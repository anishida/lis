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

#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include "lislib.h"

/************************************************
 * lis_sort_i
 * lis_sort_id
 * lis_sort_id_block
 * lis_sort_ii
 * lis_sort_iid
 * lis_sort_iiid
 * lis_sortr_ii
 * lis_bswap_int
 * lis_bswap_double
 * lis_bswap_scalar
 * lis_bswap_size_t
 * lis_sort_jad
 * lis_sort_di
 * lis_sort_d
 * lis_sort_dd
 ************************************************/

#define lis_swap(a,b,t) {t = a; a = b; b = t;}

void lis_sort_i(LIS_INT is, LIS_INT ie, LIS_INT *i1)
{
	LIS_INT i,j;
	LIS_INT p,t1;
	LIS_INT v;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = i1[p];
	lis_swap(i1[p],i1[ie],t1);

	i = is; j = ie;
	while(i<=j)
	{
		while(i1[i] < v) { i++; }
		while(i1[j] > v) { j--; }
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			i++; j--;
		}
	}
	lis_sort_i(is,j ,i1);
	lis_sort_i(i ,ie,i1);
}

void lis_sort_id(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_SCALAR *d1)
{
	LIS_INT i,j;
	LIS_INT p,t1;
	LIS_INT v;
	LIS_SCALAR s1;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = i1[p];
	lis_swap(i1[p],i1[ie],t1);
	lis_swap(d1[p],d1[ie],s1);

	i = is; j = ie;
	while(i<=j)
	{
		while(i1[i] < v) { i++; }
		while(i1[j] > v) { j--; }
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			lis_swap(d1[i],d1[j],s1);
			i++; j--;
		}
	}
	lis_sort_id(is,j ,i1,d1);
	lis_sort_id(i ,ie,i1,d1);
}

void lis_sort_id_block(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_SCALAR *d1, LIS_INT bs)
{
	LIS_INT i,j;
	LIS_INT p,t1;
	LIS_INT v;
	LIS_SCALAR s1[9];

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = i1[p];
	lis_swap(i1[p],i1[ie],t1);
	memcpy(s1,&d1[bs*p],bs*sizeof(LIS_SCALAR));
	memcpy(&d1[bs*p],&d1[bs*ie],bs*sizeof(LIS_SCALAR));
	memcpy(&d1[bs*ie],s1,bs*sizeof(LIS_SCALAR));

	i = is; j = ie;
	while(i<=j)
	{
		while(i1[i] < v) { i++; }
		while(i1[j] > v) { j--; }
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			memcpy(s1,&d1[bs*i],bs*sizeof(LIS_SCALAR));
			memcpy(&d1[bs*i],&d1[bs*j],bs*sizeof(LIS_SCALAR));
			memcpy(&d1[bs*j],s1,bs*sizeof(LIS_SCALAR));
			i++; j--;
		}
	}
	lis_sort_id_block(is,j ,i1,d1,bs);
	lis_sort_id_block(i ,ie,i1,d1,bs);
}

void lis_sort_ii(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2)
{
	LIS_INT i,j;
	LIS_INT p,t1,t2;
	LIS_INT v;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = i1[p];
	lis_swap(i1[p],i1[ie],t1);
	lis_swap(i2[p],i2[ie],t2);

	i = is; j = ie;
	while(i<=j)
	{
		while(i1[i] < v) { i++; }
		while(i1[j] > v) { j--; }
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			lis_swap(i2[i],i2[j],t2);
			i++; j--;
		}
	}
	lis_sort_ii(is,j ,i1,i2);
	lis_sort_ii(i ,ie,i1,i2);
}

void lis_sort_iid(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2, LIS_SCALAR *d1)
{
	LIS_INT i,j;
	LIS_INT p,t1,t2;
	LIS_INT v;
	LIS_SCALAR s1;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = i1[p];
	lis_swap(i1[p],i1[ie],t1);
	lis_swap(i2[p],i2[ie],t2);
	lis_swap(d1[p],d1[ie],s1);

	i = is; j = ie;
	while(i<=j)
	{
		while(i1[i] < v) { i++; }
		while(i1[j] > v) { j--; }
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			lis_swap(i2[i],i2[j],t2);
			lis_swap(d1[i],d1[j],s1);
			i++; j--;
		}
	}
	lis_sort_iid(is,j ,i1,i2,d1);
	lis_sort_iid(i ,ie,i1,i2,d1);
}

void lis_sort_iiid(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2, LIS_INT *i3, LIS_SCALAR *d1)
{
	LIS_INT i,j;
	LIS_INT p,t1,t2,t3;
	LIS_INT v;
	LIS_SCALAR s1;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = i1[p];
	lis_swap(i1[p],i1[ie],t1);
	lis_swap(i2[p],i2[ie],t2);
	lis_swap(i3[p],i3[ie],t3);
	lis_swap(d1[p],d1[ie],s1);

	i = is; j = ie;
	while(i<=j)
	{
		while(i1[i] < v) { i++; }
		while(i1[j] > v) { j--; }
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			lis_swap(i2[i],i2[j],t2);
			lis_swap(i3[i],i3[j],t3);
			lis_swap(d1[i],d1[j],s1);
			i++; j--;
		}
	}
	lis_sort_iiid(is,j ,i1,i2,i3,d1);
	lis_sort_iiid(i ,ie,i1,i2,i3,d1);
}

void lis_sortr_ii(LIS_INT is, LIS_INT ie, LIS_INT *i1, LIS_INT *i2)
{
	LIS_INT i,j;
	LIS_INT p,t1,t2;
	LIS_INT v;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = i1[p];
	lis_swap(i1[p],i1[ie],t1);
	lis_swap(i2[p],i2[ie],t2);

	i = is; j = ie;
	while(i<=j)
	{
		while(i1[i] > v) { i++; }
		while(i1[j] < v) { j--; }
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			lis_swap(i2[i],i2[j],t2);
			i++; j--;
		}
	}
	lis_sortr_ii(is,j ,i1,i2);
	lis_sortr_ii(i ,ie,i1,i2);
}

LIS_INT lis_bswap_int(LIS_INT n, LIS_INT *buf)
{
	LIS_INT i;
	LIS_INT t;
	LIS_INT *p;
	char *pp, *pt;

	p = buf;
	for(i=0;i<n;i++)
	{
		t = *p;
		pp = (char *)p;
		pt = (char *)&t;
		pp[0] = pt[3];
		pp[1] = pt[2];
		pp[2] = pt[1];
		pp[3] = pt[0];
		p++;
	}
	return LIS_SUCCESS;
}

LIS_INT lis_bswap_double(LIS_INT n, double *buf)
{
	LIS_INT i;
	double t;
	double *p;
	char *pp, *pt;

	p = buf;
	for(i=0;i<n;i++)
	{
		t = *p;
		pp = (char *)p;
		pt = (char *)&t;
		pp[0] = pt[7];
		pp[1] = pt[6];
		pp[2] = pt[5];
		pp[3] = pt[4];
		pp[4] = pt[3];
		pp[5] = pt[2];
		pp[6] = pt[1];
		pp[7] = pt[0];
		p++;
	}
	return LIS_SUCCESS;
}

LIS_INT lis_bswap_scalar(LIS_INT n, LIS_SCALAR *buf)
{
	LIS_INT i;
	LIS_SCALAR t;
	LIS_SCALAR *p;
	char *pp, *pt;

	p = buf;
	for(i=0;i<n;i++)
	{
		t = *p;
		pp = (char *)p;
		pt = (char *)&t;
		pp[0] = pt[7];
		pp[1] = pt[6];
		pp[2] = pt[5];
		pp[3] = pt[4];
		pp[4] = pt[3];
		pp[5] = pt[2];
		pp[6] = pt[1];
		pp[7] = pt[0];
		p++;
	}
	return LIS_SUCCESS;
}

LIS_INT lis_bswap_size_t(LIS_INT n, size_t *buf)
{
	LIS_INT i;
	size_t t;
	size_t *p;
	char *pp, *pt;

	p = buf;
	for(i=0;i<n;i++)
	{
		t = *p;
		pp = (char *)p;
		pt = (char *)&t;
		pp[0] = pt[7];
		pp[1] = pt[6];
		pp[2] = pt[5];
		pp[3] = pt[4];
		pp[4] = pt[3];
		pp[5] = pt[2];
		pp[6] = pt[1];
		pp[7] = pt[0];
		p++;
	}
	return LIS_SUCCESS;
}

void lis_sort_jad(LIS_INT is, LIS_INT ie, LIS_INT maxnzr, LIS_INT *i1, LIS_INT *i2)
{
	LIS_INT i,j;
	LIS_INT *iw,*iw2;

	iw  = (LIS_INT *)lis_malloc((maxnzr+2)*sizeof(LIS_INT),"lis_sort_jad::iw");
	iw2 = (LIS_INT *)lis_malloc((maxnzr+2)*sizeof(LIS_INT),"lis_sort_jad::iw2");

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0;i<maxnzr+2;i++)
	{
		iw[i] = 0;
	}
	for(i=is;i<ie;i++)
	{
		iw[(maxnzr - i1[i])+1]++;
	}
	iw[0] = is;
	for(i=0;i<maxnzr+1;i++)
	{
		iw[i+1] += iw[i];
	}
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	for(i=0;i<maxnzr+2;i++)
	{
		iw2[i] = iw[i];
	}

	for(i=is;i<ie;i++)
	{
		i2[iw[maxnzr - i1[i]]] = i;
		iw[maxnzr - i1[i]]++;
	}
	for(i=0;i<maxnzr+1;i++)
	{
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		for(j=iw2[i];j<iw2[i+1];j++)
		{
			i1[j] = maxnzr - i;
		}
	}
	lis_free2(2,iw,iw2);
}

void lis_sort_di(LIS_INT is, LIS_INT ie, LIS_SCALAR *d1, LIS_INT *i1)
{
	LIS_INT i,j;
	LIS_INT p,t1;
	LIS_SCALAR v;
	LIS_SCALAR s1;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = d1[p];
	lis_swap(i1[p],i1[ie],t1);
	lis_swap(d1[p],d1[ie],s1);

	i = is; j = ie;
	while(i<=j)
	{
#ifdef _COMPLEX
#ifdef _LONG__DOUBLE
	  while(creall(d1[i]) < creall(v)) { i++; }
	  while(creall(d1[j]) > creall(v)) { j--; }
#else
	  while(creal(d1[i]) < creal(v)) { i++; }
	  while(creal(d1[j]) > creal(v)) { j--; }
#endif
#else
	  while(d1[i] < v) { i++; }
	  while(d1[j] > v) { j--; }
#endif	  
		if( i<=j )
		{
			lis_swap(i1[i],i1[j],t1);
			lis_swap(d1[i],d1[j],s1);
			i++; j--;
		}
	}
	lis_sort_di(is,j ,d1,i1);
	lis_sort_di(i ,ie,d1,i1);
}

void lis_sort_d(LIS_INT is, LIS_INT ie, LIS_SCALAR *d1)
{
	LIS_INT i,j;
	LIS_INT p;
	LIS_SCALAR v;
	LIS_SCALAR s1;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = d1[p];
	lis_swap(d1[p],d1[ie],s1);

	i = is; j = ie;
	while(i<=j)
	{
#ifdef _COMPLEX
#ifdef _LONG__DOUBLE	  
	  while(creall(d1[i]) < creall(v)) { i++; }
	  while(creall(d1[j]) > creall(v)) { j--; }
#else
	  while(creal(d1[i]) < creal(v)) { i++; }
	  while(creal(d1[j]) > creal(v)) { j--; }
#endif
#else
	  while(d1[i] < v) { i++; }
	  while(d1[j] > v) { j--; }
#endif
		if( i<=j )
		{
			lis_swap(d1[i],d1[j],s1);
			i++; j--;
		}
	}
	lis_sort_d(is,j ,d1);
	lis_sort_d(i ,ie,d1);
}

void lis_sort_dd(LIS_INT is, LIS_INT ie, LIS_SCALAR *d1, LIS_VECTOR *d2)
{
	LIS_INT i,j;
	LIS_INT p;
	LIS_SCALAR v;
	LIS_SCALAR s1;
	LIS_VECTOR t1;

	if( ie <= is ) return;

	p = (is+ie)/2;
	v = d1[p];
	lis_swap(d1[p],d1[ie],s1);
	lis_swap(d2[p],d2[ie],t1);

	i = is; j = ie;
	while(i<=j)
	{
#ifdef _COMPLEX
#ifdef _LONG__DOUBLE	  	  
	  while(creall(d1[i]) < creall(v)) { i++; }
	  while(creall(d1[j]) > creall(v)) { j--; }
#else
	  while(creal(d1[i]) < creal(v)) { i++; }
	  while(creal(d1[j]) > creal(v)) { j--; }
#endif
#else
	  while(d1[i] < v) { i++; }
	  while(d1[j] > v) { j--; }
#endif
		if( i<=j )
		{
			lis_swap(d1[i],d1[j],s1);
			lis_swap(d2[i],d2[j],t1);
			i++; j--;
		}
	}
	lis_sort_dd(is,j ,d1,d2);
	lis_sort_dd(i ,ie,d1,d2);
}

