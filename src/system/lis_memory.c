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
/*#undef malloc*/

#include <stdlib.h>
#include <stdint.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include "lislib.h"

/************************************************
 * *rpl_malloc
 * *lis_malloc
 * *lis_calloc
 * *lis_realloc
 * lis_free
 * lis_free_mat
 * lis_free2
 * lis_free_all
 * lis_is_malloc
 * lis_free2
 * lis_is_malloc
 ************************************************/

#ifdef _DEBUG
	#define USE_MALLOC_TAG	1
#endif

typedef struct _malloc_address {
	struct _malloc_address *next, *prev;
	struct _malloc_address *hash_next;
	void *address, *real_address;
	size_t size;
	#ifdef USE_MALLOC_TAG
		char *tag;
	#endif
} malloc_address;

malloc_address malloc_address_top =
  {&malloc_address_top, &malloc_address_top,
   NULL, NULL, NULL, 0};

#define LIS_MALLOC_ADDRESS_HASH_INITIAL_SIZE 1024

static malloc_address **malloc_address_table = NULL;
static size_t malloc_address_table_size = 0;
static size_t malloc_address_count = 0;


static size_t
lis_malloc_address_hash_value(
	const void *p,
	size_t table_size)
{
	uintptr_t x;

	x = (uintptr_t)p;
	x >>= 4;
	x ^= x >> 7;
	x ^= x >> 17;

	return (size_t)(x & (table_size-1));
}


static void
lis_malloc_address_hash_rebuild(size_t new_size)
{
	malloc_address **new_table;
	malloc_address *ma;
	size_t h;

	new_table = (malloc_address **)calloc(
		new_size,
		sizeof(malloc_address *));

	if( new_table==NULL )
	{
		return;
	}

	for(
		ma=malloc_address_top.next;
		ma->address;
		ma=ma->next)
	{
		h = lis_malloc_address_hash_value(
			ma->address,
			new_size);

		ma->hash_next = new_table[h];
		new_table[h] = ma;
	}

	free(malloc_address_table);

	malloc_address_table = new_table;
	malloc_address_table_size = new_size;
}


static void
lis_malloc_address_hash_prepare(void)
{
	if( malloc_address_table==NULL )
	{
		lis_malloc_address_hash_rebuild(
			LIS_MALLOC_ADDRESS_HASH_INITIAL_SIZE);

		return;
	}

	if(
		malloc_address_table_size<=((size_t)-1)/2
		&&
		malloc_address_count>=
			2*malloc_address_table_size)
	{
		lis_malloc_address_hash_rebuild(
			2*malloc_address_table_size);
	}
}


static void
lis_malloc_address_hash_insert(
	malloc_address *ma)
{
	size_t h;

	if( malloc_address_table==NULL )
	{
		return;
	}

	h = lis_malloc_address_hash_value(
		ma->address,
		malloc_address_table_size);

	ma->hash_next = malloc_address_table[h];
	malloc_address_table[h] = ma;
}


static void
lis_malloc_address_hash_remove(
	malloc_address *ma)
{
	malloc_address **link;
	size_t h;

	if(
		malloc_address_table==NULL
		||
		ma==NULL)
	{
		return;
	}

	h = lis_malloc_address_hash_value(
		ma->address,
		malloc_address_table_size);

	link = &malloc_address_table[h];

	while(
		*link!=NULL
		&&
		*link!=ma)
	{
		link = &(*link)->hash_next;
	}

	if( *link==ma )
	{
		*link = ma->hash_next;
	}

	ma->hash_next = NULL;
}


static malloc_address *
lis_malloc_address_find(void *p)
{
	malloc_address *ma;
	size_t h;

	if( p==NULL )
	{
		return NULL;
	}

	if( malloc_address_table!=NULL )
	{
		h = lis_malloc_address_hash_value(
			p,
			malloc_address_table_size);

		for(
			ma=malloc_address_table[h];
			ma!=NULL;
			ma=ma->hash_next)
		{
			if( ma->address==p )
			{
				return ma;
			}
		}

		return NULL;
	}

	for(
		ma=malloc_address_top.next;
		ma->address;
		ma=ma->next)
	{
		if( ma->address==p )
		{
			return ma;
		}
	}

	return NULL;
}


static void
lis_malloc_address_hash_reset(void)
{
	free(malloc_address_table);

	malloc_address_table = NULL;
	malloc_address_table_size = 0;
	malloc_address_count = 0;
}

/*
void *malloc ();
 
void *rpl_malloc( size_t n )
{
	if( n == 0 )
	n = 1;
	return malloc( n );
}
*/

void *lis_malloc( size_t size, char *tag )
{
	union {
		void *ptr;
		unsigned long l;
	} addr;
	malloc_address *ma;
	LIS_INT sz;

	addr.ptr = malloc(size + 16);
	ma = (malloc_address *)malloc(sizeof(malloc_address));
	lis_malloc_address_hash_prepare();
	ma->next = &malloc_address_top;
	ma->prev = malloc_address_top.prev;
	ma->prev->next = ma;
	ma->next->prev = ma;
	ma->hash_next = NULL;
	ma->real_address = addr.ptr;
	ma->size = size;
	addr.l += 0xf;
	addr.l &= ~0xf;
	ma->address = addr.ptr;

	#ifdef USE_MALLOC_TAG
		sz = (LIS_INT)strlen(tag);
		ma->tag = (char *)malloc(sz+1);
		strcpy(ma->tag,tag);
	#endif
 
  lis_malloc_address_hash_insert(ma);
  malloc_address_count++;

  return addr.ptr;
}

void *lis_calloc( size_t size, char *tag )
{
  union {
    void *ptr;
    unsigned long l;
  } addr;
  malloc_address *ma;
  LIS_INT sz;

	addr.ptr = malloc(size + 16);
	ma = (malloc_address *)malloc(sizeof(malloc_address));
	lis_malloc_address_hash_prepare();
	ma->next = &malloc_address_top;
	ma->prev = malloc_address_top.prev;
	ma->prev->next = ma;
	ma->next->prev = ma;
	ma->hash_next = NULL;
	ma->real_address = addr.ptr;
	ma->size = size;
	addr.l += 0xf;
	addr.l &= ~0xf;
	ma->address = addr.ptr;
	memset(addr.ptr,0,size);

	#ifdef USE_MALLOC_TAG
		sz = (LIS_INT)strlen(tag);
		ma->tag = (char *)malloc(sz+1);
		strcpy(ma->tag,tag);
	#endif

  lis_malloc_address_hash_insert(ma);
  malloc_address_count++;

  return addr.ptr;
}

void *lis_realloc( void *p, size_t size )
{
	union {
		void *ptr;
		unsigned long l;
	} addr;
	malloc_address *ma;
	size_t old_size;
	void *real_address,*address;

	ma = lis_malloc_address_find(p);
	if( ma!=NULL )
	{
		old_size       = ma->size;
		addr.ptr       = malloc(size + 16);
		real_address   = addr.ptr;
		addr.l        += 0xf;
		addr.l        &= ~0xf;
		address        = addr.ptr;

		memcpy(address,p,old_size<size ? old_size : size);
		lis_malloc_address_hash_remove(ma);
		free(ma->real_address);

		ma->address      = address;
		ma->real_address = real_address;
		ma->size         = size;
		lis_malloc_address_hash_insert(ma);
		return address;
	}
	address = realloc(p,size);
	return address;
}

void lis_free(void *p)
{
  malloc_address *ma;

  ma = lis_malloc_address_find(p);
	if( ma!=NULL ) {
    lis_malloc_address_hash_remove(ma);
		if( malloc_address_count>0 )
		{
			malloc_address_count--;
		}
		ma->next->prev = ma->prev;
    ma->prev->next = ma->next;
    free(ma->real_address);
#ifdef USE_MALLOC_TAG
    free(ma->tag);
#endif
    free(ma); 
    return;
  }
  free(p);
  return;
}

void lis_free_mat(LIS_MATRIX A)
{
	malloc_address *ma;
	malloc_address *ma_bak;
	LIS_INT i;

	ma = malloc_address_top.next;
	for(i=0;i<A->n;i++)
	{
		if( A->w_index[i]==NULL ) continue;
		while (ma->address) {
			if (ma->address == A->w_index[i]) break;
			ma = ma->next;
		}
		if (ma->address) {
			lis_malloc_address_hash_remove(ma);
			if( malloc_address_count>0 )
			{
				malloc_address_count--;
			}
			ma->next->prev = ma->prev;
			ma->prev->next = ma->next;
			ma_bak         = ma->next;
			free(ma->real_address);
			#ifdef USE_MALLOC_TAG
				free(ma->tag);
			#endif
			free(ma);
			ma = ma_bak;
		}
		while (ma->address) {
			if (ma->address == A->w_value[i]) break;
			ma = ma->next;
		}
		if (ma->address) {
			lis_malloc_address_hash_remove(ma);
			if( malloc_address_count>0 )
			{
				malloc_address_count--;
			}
			ma->next->prev = ma->prev;
			ma->prev->next = ma->next;
			ma_bak         = ma->next;
			free(ma->real_address);
			#ifdef USE_MALLOC_TAG
				free(ma->tag);
			#endif
			free(ma);
			ma = ma_bak;
		}
	}
	return;
}

void lis_free2(LIS_INT n, ...)
{
	va_list vvlist;
	LIS_INT	i;
	char *p;

	va_start( vvlist, n );
	for(i=0;i<n;i++)
	{
		p = va_arg(vvlist, char *);
		if( p )
		{
			lis_free(p);
		}
	}
	va_end(vvlist);
}

void lis_free_all()
{
	malloc_address *ma;
	malloc_address *ma_bak;

	ma = malloc_address_top.next;
	while (ma->address) {
#ifdef USE_MALLOC_TAG
		lis_printf(LIS_COMM_WORLD,"memory leak: address = %x size=%d (%s)\n",ma->address,ma->size,ma->tag);
#endif
		ma_bak = ma->next;
		/*  if (ma->real_address) free(ma->real_address); */
		ma->real_address = NULL;
#ifdef USE_MALLOC_TAG
		if (ma->tag) lis_free(ma->tag);
		ma->tag = NULL;
#endif
		ma = ma_bak;
	}
	malloc_address_top.next         = &malloc_address_top;
	malloc_address_top.prev         = &malloc_address_top;
	malloc_address_top.address      = NULL;
	malloc_address_top.real_address = NULL;
	malloc_address_top.size         = 0;
	lis_malloc_address_hash_reset();
	return;
}

LIS_INT lis_is_malloc( void *p )
{
	if( lis_malloc_address_find(p)!=NULL )
	{
		return LIS_TRUE;
	}

	return LIS_FALSE;
}
