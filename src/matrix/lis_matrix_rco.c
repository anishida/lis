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
 * lis_matrix_malloc
 * lis_matrix_realloc
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_rco"
LIS_INT lis_matrix_create_rco(LIS_INT local_n, LIS_INT global_n, LIS_Comm comm, LIS_INT annz, LIS_INT *nnz, LIS_MATRIX *Amat)
{
	LIS_INT nprocs,my_rank;
	LIS_INT is,ie,err;
	LIS_INT i,k;
	LIS_INT *ranges;

	LIS_DEBUG_FUNC_IN;

	*Amat = NULL;

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
	if( local_n==0 && global_n==0 )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%D) and global n(=%D) are 0\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}

	*Amat = (LIS_MATRIX)lis_malloc( sizeof(struct LIS_MATRIX_STRUCT),"lis_matrix_create_rco::Amat" );
	if( NULL==*Amat )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_MATRIX_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_matrix_init(Amat);

	err = lis_ranges_create(comm,&local_n,&global_n,&ranges,&is,&ie,&nprocs,&my_rank);
	if( err )
	{
		lis_matrix_destroy(*Amat);
		*Amat = NULL;
		return err;
	}
	(*Amat)->ranges      = ranges;

	(*Amat)->w_nnz  = (LIS_INT *)lis_malloc(local_n*sizeof(LIS_INT),"lis_matrix_create_rco::Amat->w_nnz");
	if( (*Amat)->w_nnz==NULL )
	{
		LIS_SETERR_MEM(local_n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	if( nnz==NULL )
	{
		(*Amat)->w_annz = annz;
		for(k=0;k<local_n;k++) (*Amat)->w_nnz[k] = (*Amat)->w_annz;
	}
	else
	{
		i = 0;
		for(k=0;k<local_n;k++)
		{
			(*Amat)->w_nnz[k]  = nnz[k];
			i                 += nnz[k];
		}
		(*Amat)->w_annz = i/local_n;
	}
	err = lis_matrix_malloc_rco(local_n,(*Amat)->w_nnz,&(*Amat)->w_row,&(*Amat)->w_index,&(*Amat)->w_value);
	if( err )
	{
		lis_free((*Amat)->w_nnz);
		return err;
	}
	(*Amat)->status = LIS_MATRIX_ASSEMBLING;

	(*Amat)->n           = local_n;
	(*Amat)->gn          = global_n;
	(*Amat)->np          = local_n;
	(*Amat)->comm        = comm;
	(*Amat)->my_rank     = my_rank;
	(*Amat)->nprocs      = nprocs;
	(*Amat)->is          = is;
	(*Amat)->ie          = ie;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_rco"
LIS_INT lis_matrix_malloc_rco(LIS_INT n, LIS_INT nnz[], LIS_INT **row, LIS_INT ***index, LIS_SCALAR ***value)
{
	LIS_INT	i,j;
	LIS_INT *w_row,**w_index;
	LIS_SCALAR **w_value;

	LIS_DEBUG_FUNC_IN;

	w_row     = NULL;
	w_index   = NULL;
	w_value   = NULL;

	w_row = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_malloc_rco::w_row" );
	if( w_row==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	w_index = (LIS_INT **)lis_malloc( n*sizeof(LIS_INT *),"lis_matrix_malloc_rco::w_index" );
	if( w_index==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT *));
		lis_free2(3,w_row,w_index,w_value);
		return LIS_OUT_OF_MEMORY;
	}
	w_value = (LIS_SCALAR **)lis_malloc( n*sizeof(LIS_SCALAR *),"lis_matrix_malloc_rco::w_value" );
	if( w_value==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR *));
		lis_free2(3,w_row,w_index,w_value);
		return LIS_OUT_OF_MEMORY;
	}
	if( nnz!=NULL )
	{
		for(i=0;i<n;i++)
		{
			w_index[i] = NULL;
			w_value[i] = NULL;
			if( nnz[i]==0 ) continue;
			w_index[i] = (LIS_INT *)lis_malloc( nnz[i]*sizeof(LIS_INT),"lis_matrix_malloc_rco::w_index[i]" );
			if( w_index[i]==NULL )
			{
				LIS_SETERR_MEM(nnz[i]*sizeof(LIS_INT));
				break;
			}
			w_value[i] = (LIS_SCALAR *)lis_malloc( nnz[i]*sizeof(LIS_SCALAR),"lis_matrix_malloc_rco::w_value[i]" );
			if( w_value[i]==NULL )
			{
				LIS_SETERR_MEM(nnz[i]*sizeof(LIS_SCALAR));
				break;
			}
		}
		if(i<n)
		{
			for(j=0;j<i;j++)
			{
				if( w_index[i] ) lis_free(w_index[i]);
				if( w_value[i] ) lis_free(w_value[i]);
			}
			lis_free2(3,w_row,w_index,w_value);
			return LIS_OUT_OF_MEMORY;
		}
	}
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n;i++) w_row[i] = 0;
	*row   = w_row;
	*index = w_index;
	*value = w_value;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_realloc_rco"
LIS_INT lis_matrix_realloc_rco(LIS_INT row, LIS_INT nnz, LIS_INT ***index, LIS_SCALAR ***value)
{
	LIS_INT **w_index;
	LIS_SCALAR **w_value;

	LIS_DEBUG_FUNC_IN;

	w_index = *index;
	w_value = *value;

	w_index[row] = (LIS_INT *)lis_realloc(w_index[row],nnz*sizeof(LIS_INT));
	if( w_index[row]==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	w_value[row] = (LIS_SCALAR *)lis_realloc(w_value[row],nnz*sizeof(LIS_SCALAR));
	if( w_value[row]==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	*index = w_index;
	*value = w_value;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_rco2csr"
LIS_INT lis_matrix_convert_rco2csr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k,n,nnz,err;
	LIS_INT *ptr,*index;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;

	n       = Ain->n;
	nnz     = 0;
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+:nnz) private(i)
	#endif
	for(i=0;i<n;i++)
	{
		nnz += Ain->w_row[i];
	}

	err = lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
	if( err )
	{
		return err;
	}

	#ifdef _NUMA
		#pragma omp parallel for private(i)
		for(i=0;i<n+1;i++) ptr[i] = 0;
	#else
		ptr[0] = 0;
	#endif
	for(i=0;i<n;i++)
	{
		ptr[i+1] = ptr[i] + Ain->w_row[i];
	}
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k)
	#endif
	for(i=0;i<n;i++)
	{
		k = ptr[i];
		for(j=0;j<Ain->w_row[i];j++)
		{
			index[k] = Ain->w_index[i][j];
			value[k] = Ain->w_value[i][j];
			k++;
		}
	}

	err = lis_matrix_set_csr(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(3,ptr,index,value);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_rco2bsr"
LIS_INT lis_matrix_convert_rco2bsr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k,n,gn,nnz,bnnz,nr,nc,bnr,bnc,err;
	LIS_INT ii,jj,kk,bj,jpos,ij,kv,bi;
	LIS_INT *iw,*iw2;
	LIS_INT *bptr,*bindex;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	bnr     = Ain->conv_bnr;
	bnc     = Ain->conv_bnc;
	n       = Ain->n;
	gn      = Ain->gn;
	nr      = 1 + (n-1)/bnr;
	nc      = 1 + (gn-1)/bnc;
	bptr    = NULL;
	bindex  = NULL;
	value   = NULL;
	iw      = NULL;
	iw2     = NULL;


	bptr = (LIS_INT *)lis_malloc( (nr+1)*sizeof(LIS_INT),"lis_matrix_convert_rco2bsr::bptr" );
	if( bptr==NULL )
	{
		LIS_SETERR_MEM((nr+1)*sizeof(LIS_INT));
		lis_free2(5,bptr,bindex,value,iw,iw2);
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel private(i,k,ii,j,bj,kk,ij,jj,iw,iw2,kv,jpos)
	#endif
	{
		iw    = (LIS_INT *)lis_malloc( nc*sizeof(LIS_INT),"lis_matrix_convert_rco2bsr::iw" );
		iw2   = (LIS_INT *)lis_malloc( nc*sizeof(LIS_INT),"lis_matrix_convert_rco2bsr::iw2" );
		memset(iw,0,nc*sizeof(LIS_INT));

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr;i++)
		{
			k = 0;
			kk   = bnr*i;
			jj   = 0;
			for(ii=0;ii+kk<n&&ii<bnr;ii++)
			{
				for(j=0;j<Ain->w_row[kk+ii];j++)
				{
					bj   = Ain->w_index[kk+ii][j]/bnc;
					jpos = iw[bj];
					if( jpos==0 )
					{
						iw[bj] = 1;
						iw2[jj] = bj;
						jj++;
					}
				}
			}
			for(bj=0;bj<jj;bj++)
			{
				k++;
				ii = iw2[bj];
				iw[ii]=0;
			}
			bptr[i+1] = k;
		}
		lis_free(iw);
		lis_free(iw2);
	}

	bptr[0] = 0;
	for(i=0;i<nr;i++)
	{
		bptr[i+1] += bptr[i];
	}
	bnnz = bptr[nr];
	nnz  = bnnz*bnr*bnc;
	
	bindex = (LIS_INT *)lis_malloc( bnnz*sizeof(LIS_INT),"lis_matrix_convert_rco2bsr::bindex" );
	if( bindex==NULL )
	{
		LIS_SETERR_MEM((nr+1)*sizeof(LIS_INT));
		lis_free2(3,bptr,bindex,value);
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_rco2bsr::value" );
	if( value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		lis_free2(3,bptr,bindex,value);
		return LIS_OUT_OF_MEMORY;
	}

	/* convert bsr */
	#ifdef _OPENMP
	#pragma omp parallel private(bi,i,ii,k,j,bj,jpos,kv,kk,ij,jj,iw)
	#endif
	{
		iw = (LIS_INT *)lis_malloc( nc*sizeof(LIS_INT),"lis_matrix_convert_rco2bsr::iw" );
		memset(iw,0,nc*sizeof(LIS_INT));

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(bi=0;bi<nr;bi++)
		{
			i  = bi*bnr;
			ii = 0;
			kk = bptr[bi];
			while( i+ii<n && ii<=bnr-1 )
			{
				for( k=0;k<Ain->w_row[i+ii];k++)
				{
					j    = Ain->w_index[i+ii][k];
					bj   = j/bnc;
					j    = j%bnc;
					jpos = iw[bj];
					if( jpos==0 )
					{
						kv     = kk * bnr * bnc;
						iw[bj] = kv+1;
						bindex[kk]  = bj;
						for(jj=0;jj<bnr*bnc;jj++) value[kv+jj] = 0.0;
						ij = j*bnr + ii;
						value[kv+ij]   = Ain->w_value[i+ii][k];
						kk = kk+1;
					}
					else
					{
						ij = j*bnr + ii;
						value[jpos+ij-1]   = Ain->w_value[i+ii][k];
					}
				}
				ii = ii+1;
			}
			for(j=bptr[bi];j<bptr[bi+1];j++)
			{
				iw[bindex[j]] = 0;
			}
		}
		lis_free(iw);
	}

	err = lis_matrix_set_bsr(bnr,bnc,bnnz,bptr,bindex,value,Aout);
	if( err )
	{
		lis_free2(3,bptr,bindex,value);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_rco2csc"
LIS_INT lis_matrix_convert_rco2csc(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	LIS_INT i,j,k,l,n,nnz,err;
	LIS_INT *ptr,*index,*iw;
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;
	iw      = NULL;
	n       = Ain->n;


	iw = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_convert_rco2csc::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		lis_free2(4,ptr,index,value,iw);
		return LIS_OUT_OF_MEMORY;
	}
	ptr = (LIS_INT *)lis_malloc((n+1)*sizeof(LIS_INT),"lis_matrix_convert_rco2csc::ptr");
	if( ptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		lis_free2(4,ptr,index,value,iw);
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<n;i++) iw[i] = 0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<Ain->w_row[i];j++)
		{
			iw[Ain->w_index[i][j]]++;
		}
	}
	ptr[0] = 0;
	for(i=0;i<n;i++)
	{
		ptr[i+1] = ptr[i] + iw[i];
		iw[i]    = ptr[i];
	}
	nnz = ptr[n];

	index = (LIS_INT *)lis_malloc( nnz*sizeof(LIS_INT),"lis_matrix_convert_rco2csc::index" );
	if( index==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_INT));
		lis_free2(4,ptr,index,value,iw);
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_rco2csc::value" );
	if( value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		lis_free2(4,ptr,index,value,iw);
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<Ain->w_row[i];j++)
		{
			k        = Ain->w_index[i][j];
			l        = iw[k];
			value[l] = Ain->w_value[i][j];
			index[l] = i;
			iw[k]++;
		}
	}

	err = lis_matrix_set_csc(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(4,ptr,index,value,iw);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}

	lis_free(iw);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*
#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_brco"
LIS_INT lis_matrix_malloc_brco(LIS_INT n, LIS_INT **row, LIS_INT ***index, LIS_SCALAR ***value)
{
	LIS_INT	i,j;
	LIS_INT *w_row,**w_index;
	LIS_SCALAR **w_value;

	LIS_DEBUG_FUNC_IN;

	w_row     = NULL;
	w_index   = NULL;
	w_value   = NULL;

	w_row = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_malloc_brco::w_row" );
	if( w_row==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	w_index = (LIS_INT **)lis_malloc( n*sizeof(LIS_INT *),"lis_matrix_malloc_brco::w_index" );
	if( w_index==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT *));
		lis_free2(3,w_row,w_index,w_value);
		return LIS_OUT_OF_MEMORY;
	}
	w_value = (LIS_SCALAR **)lis_malloc( n*sizeof(LIS_SCALAR *),"lis_matrix_malloc_brco::w_value" );
	if( w_value==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR *));
		lis_free2(3,w_row,w_index,w_value);
		return LIS_OUT_OF_MEMORY;
	}
	#pragma omp parallel for private(i)
	for(i=0;i<n;i++)
	{
		w_row[i] = 0;
		w_index[i] = NULL;
	}
	*row   = w_row;
	*index = w_index;
	*value = w_value;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_vrco"
LIS_INT lis_matrix_malloc_vrco(LIS_INT n, LIS_INT **nnz, LIS_INT **row, LIS_INT ***index, LIS_SCALAR ****value)
{
	LIS_INT	i,j;
	LIS_INT *w_nnz, *w_row, **w_index;
	LIS_SCALAR ***w_value;

	LIS_DEBUG_FUNC_IN;

	w_nnz     = NULL;
	w_row     = NULL;
	w_index   = NULL;
	w_value   = NULL;

	w_nnz = (LIS_INT *)lis_malloc( n*sizeof(LIS_INT),"lis_matrix_malloc_vrco::w_ptr" );
	if( w_nnz==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	w_row = (LIS_INT *)lis_malloc( (n+1)*sizeof(LIS_INT),"lis_matrix_malloc_vrco::w_row" );
	if( w_row==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(LIS_INT));
		return LIS_OUT_OF_MEMORY;
	}
	w_index = (LIS_INT **)lis_malloc( n*sizeof(LIS_INT *),"lis_matrix_malloc_vrco::w_index" );
	if( w_index==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_INT *));
		lis_free2(3,w_row,w_index,w_value);
		return LIS_OUT_OF_MEMORY;
	}
	w_value = (LIS_SCALAR ***)lis_malloc( n*sizeof(LIS_SCALAR **),"lis_matrix_malloc_vrco::w_value" );
	if( w_value==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR **));
		lis_free2(3,w_row,w_index,w_value);
		return LIS_OUT_OF_MEMORY;
	}
	#pragma omp parallel for private(i)
	for(i=0;i<n;i++)
	{
		w_nnz[i] = 0;
		w_index[i] = NULL;
	}
	*nnz   = w_nnz;
	*row   = w_row;
	*index = w_index;
	*value = w_value;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

*/
