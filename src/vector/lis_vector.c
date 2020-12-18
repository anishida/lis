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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_vector_init
 * lis_vector_create
 * lis_vector_duplicate
 * lis_vector_destroy
 * lis_vector_set_value
 * lis_vector_set_values
 * lis_vector_set_values2
 * lis_vector_get_value
 * lis_vector_get_values
 * lis_vector_get_range
 * lis_vector_get_size
 * lis_vector_scatter
 * lis_vector_gather
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_vector_init"
LIS_INT lis_vector_init(LIS_VECTOR *vec)
{
	LIS_DEBUG_FUNC_IN;

	memset(*vec,0,sizeof(struct LIS_VECTOR_STRUCT));
	(*vec)->status = LIS_VECTOR_NULL;
	(*vec)->is_destroy = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_check"
LIS_INT lis_vector_check(LIS_VECTOR v, LIS_INT level)
{
	LIS_DEBUG_FUNC_IN;

	switch( level )
	{
	case LIS_VECTOR_CHECK_NULL:
		if( !lis_is_malloc(v) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"vector v is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	default:
		if( !lis_is_malloc(v) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"vector v is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( v->status<=LIS_VECTOR_ASSEMBLING )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"vector v is assembling\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_vector_create"
LIS_INT lis_vector_create(LIS_Comm comm, LIS_VECTOR *vec)
{
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_createex(LIS_PRECISION_DEFAULT,comm,vec);
	if( err ) return err;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_createex"
LIS_INT lis_vector_createex(LIS_INT precision, LIS_Comm comm, LIS_VECTOR *vec)
{

	LIS_DEBUG_FUNC_IN;

	*vec = NULL;

	*vec = (LIS_VECTOR)lis_malloc( sizeof(struct LIS_VECTOR_STRUCT),"lis_vector_createex::vec" );
	if( NULL==*vec )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_VECTOR_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_vector_init(vec);
	
	(*vec)->status      = LIS_VECTOR_NULL;
	(*vec)->precision   = precision;
	(*vec)->comm        = comm;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_size"
LIS_INT lis_vector_set_size(LIS_VECTOR vec, LIS_INT local_n, LIS_INT global_n)
{
	LIS_INT nprocs,my_rank;
	LIS_INT is,ie;
	LIS_INT i,err,precision;
	LIS_INT *ranges;

	LIS_DEBUG_FUNC_IN;

	if( global_n>0 && local_n>global_n )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) is larger than global n(=%D)\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}
	if( local_n<0 || global_n<0 )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) or global n(=%D) are less than 0\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}
//NEH
/*    if( local_n==0 && global_n==0 )*/
/*    {*/
/*        LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) and global n(=%D) are 0\n",local_n,global_n);*/
/*        return LIS_ERR_ILL_ARG;*/
/*    }*/


	err = lis_ranges_create(vec->comm,&local_n,&global_n,&ranges,&is,&ie,&nprocs,&my_rank);
	if( err )
	{
		return err;
	}
	vec->ranges      = ranges;

	precision = vec->precision;
	if( !precision )
	{
		vec->value = (LIS_SCALAR *)lis_malloc( local_n*sizeof(LIS_SCALAR),"lis_vector_set_size::vec->value" );
		if( NULL==vec->value )
		{
			LIS_SETERR_MEM(local_n*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<local_n;i++)
		{
			vec->value[i] = 0.0;
		}
	}
	else
	{
		vec->value = (LIS_SCALAR *)lis_malloc( (2*local_n + local_n%2)*sizeof(LIS_SCALAR),"lis_vector_set_size::vec->value" );
		if( NULL==vec->value )
		{
			LIS_SETERR_MEM((2*local_n+local_n%2)*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		vec->value_lo = vec->value + local_n + local_n%2;
		vec->work = (LIS_SCALAR *)lis_malloc( 32*sizeof(LIS_SCALAR),"lis_vector_set_size::vec->work" );
		if( NULL==vec->work )
		{
			LIS_SETERR_MEM(32*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<local_n;i++)
		{
			vec->value[i]    = 0.0;
			vec->value_lo[i] = 0.0;
		}
	}
	
	vec->is_copy     = LIS_TRUE;
	vec->status      = LIS_VECTOR_ASSEMBLED;
	vec->n           = local_n;
	vec->gn          = global_n;
	vec->np          = local_n;
	vec->my_rank     = my_rank;
	vec->nprocs      = nprocs;
	vec->is          = is;
	vec->ie          = ie;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*NEH support for extended "solve_kernel" workflow*/
#undef __FUNC__
#define __FUNC__ "lis_vector_psd_reset_scale"
LIS_INT lis_vector_psd_reset_scale(LIS_VECTOR vec)
{
	vec->is_scaled=LIS_FALSE;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_reuse"
LIS_INT lis_vector_reuse(LIS_VECTOR *vec)
{
	LIS_INT	err,np,precision;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(*vec,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	np = (*vec)->np;
	if( (*vec)->status==LIS_VECTOR_NULL )
	{
		precision = ((LIS_VECTOR)*vec)->precision;
		if( !precision )
		{
			(*vec)->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_reuse::vec->value" );
			if( NULL==(*vec)->value )
			{
				LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
			(*vec)->is_copy = LIS_TRUE;
		}
		else
		{
			(*vec)->value = (LIS_SCALAR *)lis_malloc( (2*np+np%2)*sizeof(LIS_SCALAR),"lis_vector_reuse::vec->value" );
			if( NULL==(*vec)->value )
			{
				LIS_SETERR_MEM((2*np+np%2)*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
			(*vec)->is_copy = LIS_TRUE;
			(*vec)->value_lo = (*vec)->value + np + np%2;
			(*vec)->work = (LIS_SCALAR *)lis_malloc( 32*sizeof(LIS_SCALAR),"lis_vector_reuse::vec->work" );
			if( NULL==(*vec)->work )
			{
				LIS_SETERR_MEM(32*sizeof(LIS_SCALAR));
				lis_vector_destroy(*vec);
				*vec = NULL;
				return LIS_OUT_OF_MEMORY;
			}
		}
	}

	(*vec)->status = LIS_VECTOR_ASSEMBLED;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_unset"
LIS_INT lis_vector_unset(LIS_VECTOR vec)
{
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(vec,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	if( vec->is_copy ) lis_free(vec->value);
	vec->value  = NULL;
	vec->status = LIS_VECTOR_NULL;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set"
LIS_INT lis_vector_set(LIS_VECTOR vec, LIS_SCALAR *value)
{
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(vec,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	if( vec->is_destroy ) lis_free(vec->value);
	vec->value   = value;
	vec->is_copy = LIS_FALSE;

	vec->status = LIS_VECTOR_ASSEMBLING;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_destroy"
LIS_INT lis_vector_destroy(LIS_VECTOR vec)
{
	LIS_DEBUG_FUNC_IN;
	if( lis_is_malloc(vec) )
	{
		if( vec->value && vec->is_destroy ) lis_free( vec->value );
		if( vec->work ) lis_free( vec->work );
		if( vec->ranges ) lis_free( vec->ranges );
		if( vec ) lis_free(vec);
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_vector_duplicate"
LIS_INT lis_vector_duplicate(void *vin, LIS_VECTOR *vout)
{
	LIS_INT precision,err;

	LIS_DEBUG_FUNC_IN;

	precision = LIS_PRECISION_DEFAULT;
	if( ((LIS_VECTOR)vin)->label==LIS_LABEL_VECTOR)
	{
		precision = ((LIS_VECTOR)vin)->precision;
	}
	else if( ((LIS_VECTOR)vin)->label!=LIS_LABEL_MATRIX)
	{
		LIS_SETERR(LIS_ERR_ILL_ARG, "First argument is not LIS_VECTOR or LIS_MATRIX\n");
		return LIS_ERR_ILL_ARG;
	}
	err = lis_vector_duplicateex(precision,vin,vout);

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_duplicateex"
LIS_INT lis_vector_duplicateex(LIS_INT precision, void *A, LIS_VECTOR *vout)
{
	LIS_INT np,pad;
	LIS_INT nprocs;
	LIS_INT i;
	#ifdef USE_MPI
		LIS_INT *ranges;
	#endif
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	if( ((LIS_VECTOR)A)->label!=LIS_LABEL_VECTOR && ((LIS_VECTOR)A)->label!=LIS_LABEL_MATRIX)
	{
		LIS_SETERR(LIS_ERR_ILL_ARG, "Second argument is not LIS_VECTOR or LIS_MATRIX\n");
		return LIS_ERR_ILL_ARG;
	}
	nprocs = ((LIS_VECTOR)A)->nprocs;
	np     = ((LIS_VECTOR)A)->np;
	pad    = ((LIS_VECTOR)A)->pad;
	*vout  = NULL;
	*vout  = (LIS_VECTOR)lis_malloc( sizeof(struct LIS_VECTOR_STRUCT),"lis_vector_duplicateex::vout" );
	if( NULL==*vout )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_VECTOR_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_vector_init(vout);


	if( !precision )
	{
		value = (LIS_SCALAR *)lis_malloc( (np+pad)*sizeof(LIS_SCALAR),"lis_vector_duplicateex::value" );
		if( NULL==value )
		{
			LIS_SETERR_MEM((np+pad)*sizeof(LIS_SCALAR));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*vout)->value = value;
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<np+pad;i++)
		{
			(*vout)->value[i] = 0.0;
		}
	}
	else
	{
		value = (LIS_SCALAR *)lis_malloc( (2*(np+pad) + (np+pad)%2)*sizeof(LIS_SCALAR),"lis_vector_duplicateex::value" );
		if( NULL==value )
		{
			LIS_SETERR_MEM((2*(np+pad) + (np+pad)%2)*sizeof(LIS_SCALAR));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*vout)->value = value;
		(*vout)->value_lo = value + np+pad + (np+pad)%2;
		(*vout)->work = (LIS_SCALAR *)lis_malloc( 32*sizeof(LIS_SCALAR),"lis_vector_duplicateex::vout->work" );
		if( NULL==(*vout)->work )
		{
			LIS_SETERR_MEM(32*sizeof(LIS_SCALAR));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#pragma _NEC ivdep
		#endif
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<np+pad;i++)
		{
			(*vout)->value[i]    = 0.0;
			(*vout)->value_lo[i] = 0.0;
		}
	}

	#ifdef USE_MPI
		ranges = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_vector_duplicateex::ranges" );
		if( ranges==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*sizeof(LIS_INT));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		for(i=0;i<nprocs+1;i++) ranges[i] = ((LIS_VECTOR)A)->ranges[i];
		(*vout)->ranges      = ranges;
	#else
		(*vout)->ranges      = NULL;
	#endif


	(*vout)->is_copy     = LIS_TRUE;
	(*vout)->status      = LIS_VECTOR_ASSEMBLED;
	(*vout)->precision   = precision;
	(*vout)->n           = ((LIS_VECTOR)A)->n;
	(*vout)->gn          = ((LIS_VECTOR)A)->gn;
	(*vout)->np          = ((LIS_VECTOR)A)->np;
	(*vout)->pad         = ((LIS_VECTOR)A)->pad;
	(*vout)->comm        = ((LIS_VECTOR)A)->comm;
	(*vout)->my_rank     = ((LIS_VECTOR)A)->my_rank;
	(*vout)->nprocs      = ((LIS_VECTOR)A)->nprocs;
	(*vout)->is          = ((LIS_VECTOR)A)->is;
	(*vout)->ie          = ((LIS_VECTOR)A)->ie;
	(*vout)->origin      = ((LIS_VECTOR)A)->origin;
	(*vout)->is_destroy  = ((LIS_VECTOR)A)->is_destroy;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_value0"
LIS_INT lis_vector_set_value0(LIS_INT flag, LIS_INT i, LIS_SCALAR value, LIS_VECTOR v)
{
	LIS_INT np,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	is  = v->is;
	ie  = v->ie;
	if( v->origin ) i--;
	if( i<is || i>=ie )
	{
		if( v->origin )
		{
			is++;
			ie++;
			i++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "i(=%D) is less than %D or not less than %D\n",i,is,ie);
		return LIS_ERR_ILL_ARG;
	}

	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_value::v->value" );
		if( NULL==v->value )
		{
			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		v->value[i-is] = value;
	}
	else
	{
		v->value[i-is] += value;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_value"
LIS_INT lis_vector_set_value(LIS_INT flag, LIS_INT i, LIS_SCALAR value, LIS_VECTOR v)
{
	LIS_INT np,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	is  = v->is;
	ie  = v->ie;
	if( v->origin ) i--;
	if( i<is || i>=ie )
	{
		if( v->origin )
		{
			is++;
			ie++;
			i++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "i(=%D) is less than %D or larger than %D\n",i,is,ie);
		return LIS_ERR_ILL_ARG;
	}

	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_value::v->value" );
		if( NULL==v->value )
		{

			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		v->value[i-is] = value;
	}
	else
	{
		v->value[i-is] += value;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values"
LIS_INT lis_vector_set_values(LIS_INT flag, LIS_INT count, LIS_INT index[], LIS_SCALAR value[], LIS_VECTOR v)
{
	LIS_INT np,i,ii,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	is  = v->is;
	ie  = v->ie;
	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_values::v->value" );
		if( NULL==v->value )
		{
			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		for(i=0;i<count;i++)
		{
			ii = index[i];
			if( v->origin ) ii--;
			if( ii<is || ii>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					ii++;
					i++;
				}
				LIS_SETERR4(LIS_ERR_ILL_ARG, "index[%D](=%D) is less than %D or not less than %D\n",i,ii,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[ii-is] = value[i];
		}
	}
	else
	{
		for(i=0;i<count;i++)
		{
			ii = index[i];
			if( v->origin ) ii++;
			if( ii<is || ii>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					ii++;
					i++;
				}
				LIS_SETERR4(LIS_ERR_ILL_ARG, "index[%D](=%D) is less than %D or not less than %D\n",i,ii,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[ii-is] += value[i];
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values2"
LIS_INT lis_vector_set_values2(LIS_INT flag, LIS_INT start, LIS_INT count, LIS_SCALAR value[], LIS_VECTOR v)
{
	LIS_INT np,i,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	is  = v->is;
	ie  = v->ie;
	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_values::v->value" );
		if( NULL==v->value )
		{
			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		for(i=0;i<count;i++)
		{
			start = i;
			if( v->origin ) start--;
			if( start<is || start>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					start++;
					i++;
				}
				LIS_SETERR3(LIS_ERR_ILL_ARG, "%D is less than %D or not less than %D\n",start,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[start-is] = value[i];
		}
	}
	else
	{
		for(i=0;i<count;i++)
		{
			start = i;
			if( v->origin ) start++;
			if( start<is || start>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					start++;
					i++;
				}
				LIS_SETERR3(LIS_ERR_ILL_ARG, "%D is less than %D or not less than %D\n",start,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[start-is] += value[i];
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_size"
LIS_INT lis_vector_get_size(LIS_VECTOR v, LIS_INT *local_n, LIS_INT *global_n)
{
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	*local_n  = v->n;
	*global_n = v->gn;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_range"
LIS_INT lis_vector_get_range(LIS_VECTOR v, LIS_INT *is, LIS_INT *ie)
{
	LIS_INT	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	*is = v->is;
	*ie = v->ie;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_value"
LIS_INT lis_vector_get_value(LIS_VECTOR v, LIS_INT i, LIS_SCALAR *value)
{
	LIS_INT err,is,ie;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	is  = v->is;
	ie  = v->ie;
	if( v->origin ) i--;
	if( i<is || i>=ie )
	{
		if( v->origin )
		{
			i++;
			is++;
			ie++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "i(=%D) is less than %D or not less than %D\n",i,is,ie);
		return LIS_ERR_ILL_ARG;
	}

	*value = v->value[i-is];
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_values"
LIS_INT lis_vector_get_values(LIS_VECTOR v, LIS_INT start, LIS_INT count, LIS_SCALAR value[])
{
	LIS_INT err,n,i,is,ie;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n   = v->n;
	is  = v->is;
	ie  = v->ie;
	if( v->origin ) start--;
	if( start<is || start>=ie )
	{
		if( v->origin )
		{
			start++;
			is++;
			ie++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "start(=%D) is less than %D or not less than %D\n",start,is,ie);
		return LIS_ERR_ILL_ARG;
	}
	if( (start-is+count)>n )
	{
		LIS_SETERR3(LIS_ERR_ILL_ARG, "start(=%D) + count(=%D) exceeds the range of vector v(=%D).\n",start,count,ie);
		return LIS_ERR_ILL_ARG;
	}
	for(i=0;i<count;i++)
	{
		value[i] = v->value[start-is + i];
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

LIS_INT lis_vector_print(LIS_VECTOR x)
{
#ifdef USE_MPI
	LIS_INT err,i,ii,is,n,k,nprocs,my_rank;

	err = lis_vector_check(x,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	nprocs  = x->nprocs;
	my_rank = x->my_rank;
	n       = x->n;
	is      = x->is;

	for(k=0;k<nprocs;k++)
	{
		if( k==my_rank )
		{
			for(i=0;i<n;i++)
			{
				ii = i+is;
				if( x->origin ) ii++;
#ifdef _COMPLEX
#ifdef _LONG__LONG
				printf("%6lld  (%e, %e)\n",ii,(double)creal(x->value[i]),(double)cimag(x->value[i]));
#else
				printf("%6d  (%e, %e)\n",ii,(double)creal(x->value[i]),(double)cimag(x->value[i]));
#endif
#else
#ifdef _LONG__LONG
				printf("%6lld  %e\n",ii,(double)x->value[i]);
#else
				printf("%6d  %e\n",ii,(double)x->value[i]);
#endif
#endif				
			}
		}
		MPI_Barrier(x->comm);
	}
	return LIS_SUCCESS;
#else
	LIS_INT err,i,ii,n;

	err = lis_vector_check(x,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n = x->n;

	for(i=0; i<n; i++)
	{
		ii = i;
		if( x->origin ) ii++;
		if( x->precision==LIS_PRECISION_DEFAULT )
		{
#ifdef _COMPLEX
#ifdef _LONG__LONG
			printf("%6lld  (%e, %e)\n",ii,(double)creal(x->value[i]),(double)cimag(x->value[i]));
#else
			printf("%6d  (%e, %e)\n",ii,(double)creal(x->value[i]),(double)cimag(x->value[i]));
#endif
#else		  
#ifdef _LONG__LONG
			printf("%6lld  %e\n",ii,(double)x->value[i]);
#else
			printf("%6d  %e\n",ii,(double)x->value[i]);
#endif
#endif			
		}
		else
		{
#ifdef _COMPLEX		  
#ifdef _LONG__LONG
			printf("%6lld  %e,%e\n",ii,(double)creal(x->value[i]),(double)x->value_lo[i]);
#else
			printf("%6d  %e,%e\n",ii,(double)creal(x->value[i]),(double)x->value_lo[i]);
#endif
#else
#ifdef _LONG__LONG
			printf("%6lld  %e,%e\n",ii,(double)x->value[i],(double)x->value_lo[i]);
#else
			printf("%6d  %e,%e\n",ii,(double)x->value[i],(double)x->value_lo[i]);
#endif
#endif			
		}
	}

	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scatter"
LIS_INT lis_vector_scatter(LIS_SCALAR value[], LIS_VECTOR v)
{
#ifdef USE_MPI
	LIS_INT err,i,is,n,nprocs,my_rank,*sendcounts;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	nprocs  = v->nprocs;
	my_rank = v->my_rank;
	n       = v->n;
	is      = v->is;

	sendcounts = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_vector_scatter::sendcounts" );

	for(i=0; i<nprocs; i++)
	{
	  sendcounts[i] = v->ranges[i+1] - v->ranges[i];
	}

	if(my_rank == 0) 
	{
	  MPI_Scatterv(&value[0],sendcounts,v->ranges,LIS_MPI_SCALAR,MPI_IN_PLACE,n,LIS_MPI_SCALAR,0,v->comm);
	}
	else
	{
	  MPI_Scatterv(&value[0],sendcounts,v->ranges,LIS_MPI_SCALAR,&value[is],n,LIS_MPI_SCALAR,0,v->comm);
	}

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  v->value[i] = value[i+is];
	}

	return LIS_SUCCESS;
#else
	LIS_INT err,i,n;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n = v->n;

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  v->value[i] = value[i];
	}

	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_vector_gather"
LIS_INT lis_vector_gather(LIS_VECTOR v, LIS_SCALAR value[])
{
#ifdef USE_MPI
	LIS_INT err,i,is,n,nprocs,*recvcounts;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	nprocs  = v->nprocs;
	n       = v->n;
	is      = v->is;

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  value[i+is] = v->value[i];
	}
	recvcounts = (LIS_INT *)lis_malloc( (nprocs+1)*sizeof(LIS_INT),"lis_vector_gather::recvcounts" );
	for(i=0; i<nprocs; i++)
	{
	  recvcounts[i] = v->ranges[i+1] - v->ranges[i];
	}
	MPI_Allgatherv(MPI_IN_PLACE,n,LIS_MPI_SCALAR,&value[0],recvcounts,v->ranges,LIS_MPI_SCALAR,v->comm);

	return LIS_SUCCESS;
#else
	LIS_INT err,i,n;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n = v->n;

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#pragma _NEC ivdep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  value[i] = v->value[i];
	}

	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_vector_is_null"
LIS_INT lis_vector_is_null(LIS_VECTOR v)
{
	LIS_DEBUG_FUNC_IN;

	if( v->status==LIS_VECTOR_NULL ) return LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_FALSE;
}
