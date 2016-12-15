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
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_output_mm_vec
 * lis_output_mm_header
 * lis_output_mm_csr
 * lis_output_mm_csc
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_output_mm_vec"
#ifdef USE_MPI
LIS_INT lis_output_mm_vec(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT format, char *path)
{
	LIS_INT	n,nprocs,my_rank;
	LIS_INT	i,is;
	LIS_INT	pe;
	LIS_INT	err,ret;
	FILE *file;
	LIS_MM_VECFMT fmt;

	LIS_DEBUG_FUNC_IN;


	nprocs    = A->nprocs;
	my_rank   = A->my_rank;
	n         = A->n;
	is        = A->is;

	if( !lis_vector_is_null(b) )
	{
		err = 0;
		for(pe=0;pe<nprocs;pe++)
		{
			if( my_rank==pe )
			{
				if( format==LIS_FMT_MM )
				{
					file = fopen(path, "a");
				}
				else
				{
					file = fopen(path, "ab");
				}
				if( file==NULL )
				{
					LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
					err = 1;
				}
				else
				{
					if( format==LIS_FMT_MM )
					{
						for(i=0;i<n;i++)
						{
						  if (b->intvalue)
						    {
#ifdef _LONG__LONG						     
						      fprintf(file, "%lld %28lld\n", is+i+1,(long long int)b->value[i]);
#else
						      fprintf(file, "%d %28d\n", is+i+1,(int)b->value[i]);
#endif
						    }
						  else
						    {
#ifdef _COMPLEX
#ifdef _LONG__LONG
							fprintf(file, "%lld %28.20e %28.20e\n", is+i+1, (double)creal(b->value[i]), (double)cimag(b->value[i]));
#else
							fprintf(file, "%d %28.20e %28.20e\n", is+i+1, (double)creal(b->value[i]), (double)cimag(b->value[i]));
#endif
#else
#ifdef _LONG__LONG
							fprintf(file, "%lld %28.20e\n", is+i+1,(double)b->value[i]);
#else
							fprintf(file, "%d %28.20e\n", is+i+1,(double)b->value[i]);
#endif
#endif							
						    }
						}
					}
					else
					{
						for(i=0;i<n;i++)
						{
							fmt.i     = is+i+1;
							fmt.value = b->value[i];
							fwrite(&fmt,sizeof(fmt),1,file);
						}
					}
					fclose(file);
				}
			}
			MPI_Allreduce(&err,&ret,1,LIS_MPI_INT,MPI_SUM,A->comm);
			if( ret )
			{
				return LIS_ERR_FILE_IO;
			}
		}
	}
	if( !lis_vector_is_null(x) )
	{
		err = 0;
		for(pe=0;pe<nprocs;pe++)
		{
			if( my_rank==pe )
			{
				if( format==LIS_FMT_MM )
				{
					file = fopen(path, "a");
				}
				else
				{
					file = fopen(path, "ab");
				}
				if( file==NULL )
				{
					LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
					err = 1;
				}
				else
				{
					if( format==LIS_FMT_MM )
					{
					  
						for(i=0;i<n;i++)
						{
						  if (x->intvalue)
						    {
#ifdef _LONG__LONG						     
						      fprintf(file, "%lld %28lld\n", is+i+1,(long long int)x->value[i]);
#else
						      fprintf(file, "%d %28d\n", is+i+1,(int)x->value[i]);
#endif
						    }
						  else
						    {
#ifdef _COMPLEX
#ifdef _LONG__LONG
							fprintf(file, "%lld %28.20e %28.20e\n", is+i+1, (double)creal(x->value[i]), (double)cimag(x->value[i]));
#else
							fprintf(file, "%d %28.20e %28.20e\n", is+i+1, (double)creal(x->value[i]), (double)cimag(x->value[i]));
#endif
#else
#ifdef _LONG__LONG
							fprintf(file, "%lld %28.20e\n", is+i+1,(double)x->value[i]);
#else
							fprintf(file, "%d %28.20e\n", is+i+1,(double)x->value[i]);
#endif
#endif							
						    }
						}
					}
					else
					{
						for(i=0;i<n;i++)
						{
							fmt.i     = is+i+1;
							fmt.value = x->value[i];
							fwrite(&fmt,sizeof(fmt),1,file);
						}
					}
					fclose(file);
				}
			}
			MPI_Allreduce(&err,&ret,1,LIS_MPI_INT,MPI_SUM,A->comm);
			if( ret )
			{
				return LIS_ERR_FILE_IO;
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#else
LIS_INT lis_output_mm_vec(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT format, FILE *file)
{
	LIS_INT	n;
	LIS_INT	i;
	LIS_MM_VECFMT fmt;

	LIS_DEBUG_FUNC_IN;

	n   = A->n;

	if( !lis_vector_is_null(b) )
	{
		if( format==LIS_FMT_MM )
		{
			for(i=0;i<n;i++)
			{
			  if (b->intvalue)
			    {
#ifdef _LONG__LONG						     
			      fprintf(file, "%lld %28lld\n", i+1,(long long int)b->value[i]);
#else
			      fprintf(file, "%d %28d\n", i+1,(int)b->value[i]);
#endif
			    }
			  else
			    {
#ifdef _COMPLEX
#ifdef _LONG__LONG
				fprintf(file, "%lld %28.20e %28.20e\n", i+1, (double)creal(b->value[i]), (double)cimag(b->value[i]));
#else
				fprintf(file, "%d %28.20e %28.20e\n", i+1, (double)creal(b->value[i]), (double)cimag(b->value[i]));
#endif
#else			      
#ifdef _LONG__LONG
				fprintf(file, "%lld %28.20e\n", i+1,(double)b->value[i]);
#else
				fprintf(file, "%d %28.20e\n", i+1,(double)b->value[i]);
#endif
#endif				
			    }
			}
		}
		else
		{
			for(i=0;i<n;i++)
			{
				fmt.i     = i+1;
				fmt.value = b->value[i];
				fwrite(&fmt,sizeof(fmt),1,file);
			}
		}
	}
	if( !lis_vector_is_null(x) )
	{
		if( format==LIS_FMT_MM )
		{			      
			for(i=0;i<n;i++)
			{
			  if (x->intvalue)
			    {
#ifdef _LONG__LONG						     
			      fprintf(file, "%lld %28lld\n", i+1,(long long int)x->value[i]);
#else
			      fprintf(file, "%d %28d\n", i+1,(int)x->value[i]);
#endif
			    }
			  else
			    {
#ifdef _COMPLEX
#ifdef _LONG__LONG
				fprintf(file, "%lld %28.20e %28.20e\n", i+1, (double)creal(x->value[i]), (double)cimag(x->value[i]));
#else
				fprintf(file, "%d %28.20e %28.20e\n", i+1, (double)creal(b->value[i]), (double)cimag(b->value[i]));
#endif
#else			      
#ifdef _LONG__LONG
				fprintf(file, "%lld %28.20e\n", i+1,(double)x->value[i]);
#else
				fprintf(file, "%d %28.20e\n", i+1,(double)b->value[i]);
#endif
#endif				
			    }
			}
		}
		else
		{
			for(i=0;i<n;i++)
			{
				fmt.i     = i+1;
				fmt.value = x->value[i];
				fwrite(&fmt,sizeof(fmt),1,file);
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_output_mm_header"
#ifdef USE_MPI
LIS_INT lis_output_mm_header(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT format, char *path, FILE **file)
{
	LIS_INT	local_nnz,nnz,my_rank;
	LIS_INT	err,ret;
	LIS_INT	isb,isx,endian;

	LIS_DEBUG_FUNC_IN;

	my_rank   = A->my_rank;
	local_nnz = A->nnz;
	isb = 0;
	isx = 0;
	if( !lis_vector_is_null(b) ) isb = 1;
	if( !lis_vector_is_null(x) ) isx = 1;
	endian = 1;
	endian = *(char *)&endian;

	err       = 0;
	MPI_Allreduce(&local_nnz,&nnz,1,LIS_MPI_INT,MPI_SUM,A->comm);
	if( my_rank==0 )
	{
		if( format==LIS_FMT_MM )
		{
			*file = fopen(path, "w");
		}
		else
		{
			*file = fopen(path, "wb");
		}
		if( file==NULL )
		{
			LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n",path);
			err = 1;
		}
		else
		{
			fprintf(*file, "%%%%MatrixMarket matrix coordinate real general\n");

			if( format==LIS_FMT_MM )
			{
				if( isb==0 && isx==0 )
				{
#ifdef _LONG__LONG
					fprintf(*file, "%lld %lld %lld\n", A->gn, A->gn, nnz);
#else
					fprintf(*file, "%d %d %d\n", A->gn, A->gn, nnz);
#endif
				}
				else
				{
#ifdef _LONG__LONG
					fprintf(*file, "%lld %lld %lld %lld %lld\n", A->gn, A->gn, nnz, isb, isx);
#else
					fprintf(*file, "%d %d %d %d %d\n", A->gn, A->gn, nnz, isb, isx);
#endif
				}
			}
			else
			{
#ifdef _LONG__LONG
				fprintf(*file, "%lld %lld %lld %lld %lld %lld\n", A->gn, A->gn, nnz, isb, isx, endian+1);
#else
				fprintf(*file, "%d %d %d %d %d %d\n", A->gn, A->gn, nnz, isb, isx, endian+1);
#endif
			}
			fclose(*file);
		}
	}
	MPI_Allreduce(&err,&ret,1,LIS_MPI_INT,MPI_SUM,A->comm);
	if( ret )
	{
		return LIS_ERR_FILE_IO;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#else
LIS_INT lis_output_mm_header(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT format, char *path, FILE **file)
{
	LIS_INT	nnz;
	LIS_INT	isb,isx,endian;

	LIS_DEBUG_FUNC_IN;

	nnz = A->nnz;
	isb = 0;
	isx = 0;
	if( !lis_vector_is_null(b) ) isb = 1;
	if( !lis_vector_is_null(x) ) isx = 1;
	endian = 1;
	endian = *(char *)&endian;

	if( format==LIS_FMT_MM )
	{
		*file = fopen(path, "w");
	}
	else
	{
		*file = fopen(path, "wb");
	}
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n",path);
		return LIS_ERR_FILE_IO;
	}
#ifdef _COMPLEX	
	fprintf(*file, "%%%%MatrixMarket matrix coordinate complex general\n");
#else
	fprintf(*file, "%%%%MatrixMarket matrix coordinate real general\n");
#endif	
	if( format==LIS_FMT_MM )
	{
		if( isb==0 && isx==0 )
		{
#ifdef _LONG__LONG
			fprintf(*file, "%lld %lld %lld\n", A->gn, A->gn, nnz);
#else
			fprintf(*file, "%d %d %d\n", A->gn, A->gn, nnz);
#endif
		}
		else
		{
#ifdef _LONG__LONG
			fprintf(*file, "%lld %lld %lld %lld %lld\n", A->gn, A->gn, nnz, isb, isx);
#else
			fprintf(*file, "%d %d %d %d %d\n", A->gn, A->gn, nnz, isb, isx);
#endif
		}
	}
	else
	{
#ifdef _LONG__LONG
		fprintf(*file, "%lld %lld %lld %lld %lld %lld\n", A->gn, A->gn, nnz, isb, isx, endian+1);
#else
		fprintf(*file, "%d %d %d %d %d %d\n", A->gn, A->gn, nnz, isb, isx, endian+1);
#endif
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_output_mm_csr"
#ifdef USE_MPI
LIS_INT lis_output_mm_csr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT format, char *path)
{
	LIS_INT	n,nprocs,my_rank;
	LIS_INT	i,j,jj,is;
	LIS_INT	pe;
	LIS_INT	err;
	FILE *file;
	LIS_MM_MATFMT fmt;

	LIS_DEBUG_FUNC_IN;

	nprocs    = A->nprocs;
	my_rank   = A->my_rank;
	n         = A->n;
	is        = A->is;

	err = lis_output_mm_header(A,b,x,format,path,&file);
	if( err ) return err;


	for(pe=0;pe<nprocs;pe++)
	{
		if( my_rank==pe )
		{
			file = fopen(path, "a");
			if( file==NULL )
			{
				LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
				return LIS_FAILS;
			}
			if( A->matrix_type==LIS_MATRIX_CSR )
			{
				if( format==LIS_FMT_MM )
				{
					for(i=0;i<n;i++)
					{
						for(j=A->ptr[i];j<A->ptr[i+1];j++)
						{
							if( A->index[j]>=n )
							{
								jj = A->l2g_map[A->index[j]-n]+1;
							}
							else
							{
								jj = is + A->index[j]+1;
							}
#ifdef _COMPLEX
#ifdef _LONG__LONG
							fprintf(file, "%lld %lld %28.20e %28.20e\n", is+i+1, jj, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#else
							fprintf(file, "%d %d %28.20e %28.20e\n", is+i+1, jj, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#endif
#else							
#ifdef _LONG__LONG
							fprintf(file, "%lld %lld %28.20e\n", is+i+1,jj,(double)A->value[j]);
#else
							fprintf(file, "%d %d %28.20e\n", is+i+1,jj,(double)A->value[j]);
#endif
#endif							
						}
					}
				}
				else
				{
					for(i=0;i<n;i++)
					{
						for(j=A->ptr[i];j<A->ptr[i+1];j++)
						{
							if( A->index[j]>=n )
							{
								jj = A->l2g_map[A->index[j]-n]+1;
							}
							else
							{
								jj = is + A->index[j]+1;
							}
							fmt.i     = is+i+1;
							fmt.j     = jj;
							fmt.value = A->value[j];
							fwrite(&fmt,sizeof(fmt),1,file);
						}
					}
				}
			}
			else
			{
				if( format==LIS_FMT_MM )
				{
					for(i=0;i<n;i++)
					{
						for(j=A->ptr[i];j<A->ptr[i+1];j++)
						{
							if( A->index[j]>=n )
							{
								jj = A->l2g_map[A->index[j]-n]+1;
							}
							else
							{
								jj = is + A->index[j]+1;
							}
#ifdef _COMPLEX
#ifdef _LONG__LONG
							fprintf(file, "%lld %lld %28.20e %28.20e\n", jj, is+i+1, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#else
							fprintf(file, "%d %d %28.20e %28.20e\n", jj, is+i+1, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#endif
#else							
#ifdef _LONG__LONG
							fprintf(file, "%lld %lld %28.20e\n", jj,is+i+1,(double)A->value[j]);
#else
							fprintf(file, "%d %d %28.20e\n", jj,is+i+1,(double)A->value[j]);
#endif
#endif							
						}
					}
				}
				else
				{
					for(i=0;i<n;i++)
					{
						for(j=A->ptr[i];j<A->ptr[i+1];j++)
						{
							if( A->index[j]>=n )
							{
								jj = A->l2g_map[A->index[j]-n]+1;
							}
							else
							{
								jj = is + A->index[j]+1;
							}
							fmt.j     = is+i+1;
							fmt.i     = jj;
							fmt.value = A->value[j];
							fwrite(&fmt,sizeof(fmt),1,file);
						}
					}
				}
			}
			fclose(file);
		}
		MPI_Barrier(A->comm);
	}
	err = lis_output_mm_vec(A,b,x,format,path);
	if( err ) return err;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#else
LIS_INT lis_output_mm_csr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT format, char *path)
{
	LIS_INT	n,nnz;
	LIS_INT	i,j,jj;
	LIS_INT	err;
	FILE *file;
	LIS_MM_MATFMT fmt;

	LIS_DEBUG_FUNC_IN;

	n   = A->n;
	nnz = A->nnz;

	err = lis_output_mm_header(A,b,x,format,path,&file);
	if( err ) return err;

	if( A->matrix_type==LIS_MATRIX_CSR )
	{
		if( format==LIS_FMT_MM )
		{
			for(i=0;i<n;i++)
			{
				for(j=A->ptr[i];j<A->ptr[i+1];j++)
				{
					jj = A->index[j]+1;
#ifdef _COMPLEX
#ifdef _LONG__LONG
					fprintf(file, "%lld %lld %28.20e %28.20e\n", i+1, jj, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#else
					fprintf(file, "%d %d %28.20e %28.20e\n", i+1, jj, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#endif
#else					
#ifdef _LONG__LONG
					fprintf(file, "%lld %lld %28.20e\n", i+1,jj,(double)A->value[j]);
#else
					fprintf(file, "%d %d %28.20e\n", i+1,jj,(double)A->value[j]);
#endif
#endif					
				}
			}
		}
		else
		{
			for(i=0;i<n;i++)
			{
				for(j=A->ptr[i];j<A->ptr[i+1];j++)
				{
					jj = A->index[j]+1;
					fmt.i     = i+1;
					fmt.j     = jj;
					fmt.value = A->value[j];
					fwrite(&fmt,sizeof(fmt),1,file);
				}
			}
		}
	}
	else
	{
		if( format==LIS_FMT_MM )
		{
			for(i=0;i<n;i++)
			{
				for(j=A->ptr[i];j<A->ptr[i+1];j++)
				{
					jj = A->index[j]+1;
#ifdef _COMPLEX
#ifdef _LONG__LONG
					fprintf(file, "%lld %lld %28.20e %28.20e\n", jj, i+1, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#else
					fprintf(file, "%d %d %28.20e %28.20e\n", jj, i+1, (double)creal(A->value[j]), (double)cimag(A->value[j]));
#endif
#else
#ifdef _LONG__LONG
					fprintf(file, "%lld %lld %28.20e\n", jj,i+1,(double)A->value[j]);
#else
					fprintf(file, "%d %d %28.20e\n", jj,i+1,(double)A->value[j]);
#endif
#endif					
				}
			}
		}
		else
		{
			for(i=0;i<n;i++)
			{
				for(j=A->ptr[i];j<A->ptr[i+1];j++)
				{
					jj = A->index[j]+1;
					fmt.j     = i+1;
					fmt.i     = jj;
					fmt.value = A->value[j];
					fwrite(&fmt,sizeof(fmt),1,file);
				}
			}
		}
	}
	lis_output_mm_vec(A,b,x,format,file);
	fclose(file);
	
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_output_mm_csc"
LIS_INT lis_output_mm_csc(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT format, char *path)
{
	LIS_INT err;

	LIS_DEBUG_FUNC_IN;

	err = lis_output_mm_csr(A,b,x,format,path);

	LIS_DEBUG_FUNC_OUT;
	return err;

}

