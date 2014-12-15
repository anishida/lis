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
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_hashtable_create
 * lis_hashtable_destroy
 * lis_hashtable_clear
 * lis_hashtable_search
 * lis_hashtable_set_value
 * lis_hashtable_get_value
 ************************************************/

#define LIS_HASHTABLE_SIZE 1021

#undef __FUNC__
#define __FUNC__ "lis_hashtable_create"
LIS_INT lis_hashtable_create(LIS_HASHTABLE *hashtable)
{
	LIS_HASHTABLE table;

	LIS_DEBUG_FUNC_IN;

	*hashtable = NULL;

	table = (LIS_HASHTABLE)malloc(LIS_HASHTABLE_SIZE*sizeof(struct LIS_HASH_STRUCT *));
	if( table==NULL )
	{
		LIS_SETERR_MEM(LIS_HASHTABLE_SIZE*sizeof(struct LIS_HASH_STRUCT *));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	memset(table,0,LIS_HASHTABLE_SIZE*sizeof(struct LIS_HASH_STRUCT *));

	*hashtable = table;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_hashtable_destroy"
LIS_INT lis_hashtable_destroy(LIS_HASHTABLE hashtable)
{
	LIS_INT	i;
	LIS_HASH p,t;

	LIS_DEBUG_FUNC_IN;

	for(i=0;i<LIS_HASHTABLE_SIZE;i++)
	{
		p = hashtable[i];
		while( p!=NULL )
		{
			t = p;
			p = p->next;
			free(t);
		}
	}
	free(hashtable);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_hashtable_clear"
LIS_INT lis_hashtable_clear(LIS_HASHTABLE hashtable)
{
	LIS_INT i;
	LIS_HASH p,t;

	LIS_DEBUG_FUNC_IN;

	for(i=0;i<LIS_HASHTABLE_SIZE;i++)
	{
		p = hashtable[i];
		while( p!=NULL )
		{
			t = p;
			p = p->next;
			free(t);
		}
		hashtable[i] = NULL;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_hashtable_search"
LIS_HASH lis_hashtable_search(LIS_HASHTABLE hashtable, LIS_INT index)
{
	LIS_INT	hash;
	LIS_HASH p;


	hash = index%LIS_HASHTABLE_SIZE;
	for(p=hashtable[hash];p!=NULL;p=p->next)
	{
		if( p->index==index ) return p;
	}

	return NULL;
}

#undef __FUNC__
#define __FUNC__ "lis_hashtable_set_value"
LIS_INT lis_hashtable_set_value(LIS_HASHTABLE hashtable, LIS_INT index, LIS_INT value)
{
	LIS_INT hashval;
	LIS_HASH p;

	LIS_DEBUG_FUNC_IN;

	p = (LIS_HASH)malloc(sizeof(struct LIS_HASH_STRUCT));
	if( p==NULL )
	{
		LIS_SETERR_MEM(LIS_HASHTABLE_SIZE*sizeof(struct LIS_HASH_STRUCT *));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	hashval            = index%LIS_HASHTABLE_SIZE;
	p->next            = hashtable[hashval];
	p->index           = index;
	p->value           = value;
	hashtable[hashval] = p;


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_hashtable_get_value"
LIS_INT lis_hashtable_get_value(LIS_HASHTABLE hashtable, LIS_INT index)
{
	LIS_INT hashval;
	LIS_HASH p;


	hashval = index%LIS_HASHTABLE_SIZE;
	for(p=hashtable[hashval];p!=NULL;p=p->next)
	{
		if( p->index==index ) return p->value;
	}

	return 0;
}
