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

#ifndef __LIS_H__
#define __LIS_H__

#ifdef __cplusplus
#undef HAVE_COMPLEX_H
extern "C" {
#endif

/**************************************/
#define LIS_VERSION	"2.1.6"
/**************************************/
#include <stdio.h>
#ifdef USE_COMPLEX
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#else
#ifdef HAVE_COMPLEX_H
#undef HAVE_COMPLEX_H
#endif
#endif
#ifdef HAVE_QUADMATH_H
#include <quadmath.h>
#endif

#define _max(a,b) ((a) >= (b) ? (a) : (b))
#define _min(a,b) ((a) <= (b) ? (a) : (b))

#define LIS_FMT_AUTO 0
#define LIS_FMT_PLAIN 1
#define LIS_FMT_MM 2
#define LIS_FMT_LIS 3
#define LIS_FMT_LIS_ASCII 3
#define LIS_FMT_LIS_BINARY 4
#define LIS_FMT_FREE 5
#define LIS_FMT_ITBL 6
#define LIS_FMT_HB 7
#define LIS_FMT_MMB 8

#define LIS_BINARY_BIG 0
#define LIS_BINARY_LITTLE 1


#define LIS_OPTIONS_LEN 27
#define LIS_OPTIONS_SOLVER 0
#define LIS_OPTIONS_PRECON 1
#define LIS_OPTIONS_MAXITER 2
#define LIS_OPTIONS_OUTPUT 3
#define LIS_OPTIONS_RESTART 4
#define LIS_OPTIONS_ELL 5
#define LIS_OPTIONS_SCALE 6
#define LIS_OPTIONS_FILL 7
#define LIS_OPTIONS_M 8
#define LIS_OPTIONS_PSOLVER 9
#define LIS_OPTIONS_PMAXITER 10
#define LIS_OPTIONS_PRESTART 11
#define LIS_OPTIONS_PELL 12
#define LIS_OPTIONS_PPRECON 13
#define LIS_OPTIONS_ISLEVEL 14
#define LIS_OPTIONS_INITGUESS_ZEROS 15
#define LIS_OPTIONS_ADDS 16
#define LIS_OPTIONS_ADDS_ITER 17
#define LIS_OPTIONS_PRECISION 18
#define LIS_OPTIONS_USE_AT 19
#define LIS_OPTIONS_SWITCH_MAXITER 20
#define LIS_OPTIONS_SAAMG_UNSYM 21
#define LIS_OPTIONS_STORAGE 22
#define LIS_OPTIONS_STORAGE_BLOCK 23
#define LIS_OPTIONS_CONV_COND 24
#define LIS_OPTIONS_INIT_SHADOW_RESID 25
#define LIS_OPTIONS_IDRS_RESTART 26

#define LIS_EOPTIONS_LEN 13
#define LIS_EOPTIONS_ESOLVER 0
#define LIS_EOPTIONS_MAXITER 1
#define LIS_EOPTIONS_SUBSPACE 2
#define LIS_EOPTIONS_MODE 3
#define LIS_EOPTIONS_OUTPUT 4
#define LIS_EOPTIONS_INITGUESS_ONES 5
#define LIS_EOPTIONS_INNER_ESOLVER 6
#define LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER 7
#define LIS_EOPTIONS_STORAGE 8
#define LIS_EOPTIONS_STORAGE_BLOCK 9
#define LIS_EOPTIONS_PRECISION 10
#define LIS_EOPTIONS_SWITCH_MAXITER 11
#define LIS_EOPTIONS_RVAL 12

#define LIS_PARAMS_LEN 15
#define LIS_PARAMS_RESID LIS_OPTIONS_LEN+0
#define LIS_PARAMS_OMEGA LIS_OPTIONS_LEN+1
#define LIS_PARAMS_RELAX LIS_OPTIONS_LEN+2
#define LIS_PARAMS_DROP LIS_OPTIONS_LEN+3
#define LIS_PARAMS_ALPHA LIS_OPTIONS_LEN+4
#define LIS_PARAMS_TAU LIS_OPTIONS_LEN+5
#define LIS_PARAMS_SIGMA LIS_OPTIONS_LEN+6
#define LIS_PARAMS_GAMMA LIS_OPTIONS_LEN+7
#define LIS_PARAMS_SSOR_OMEGA LIS_OPTIONS_LEN+8
#define LIS_PARAMS_PRESID LIS_OPTIONS_LEN+9
#define LIS_PARAMS_POMEGA LIS_OPTIONS_LEN+10
#define LIS_PARAMS_SWITCH_RESID LIS_OPTIONS_LEN+11
#define LIS_PARAMS_RATE LIS_OPTIONS_LEN+12
#define LIS_PARAMS_RESID_WEIGHT LIS_OPTIONS_LEN+13
#define LIS_PARAMS_SAAMG_THETA LIS_OPTIONS_LEN+14

#define LIS_EPARAMS_LEN 3
#define LIS_EPARAMS_RESID LIS_EOPTIONS_LEN+0
#define LIS_EPARAMS_SHIFT LIS_EOPTIONS_LEN+1
#define LIS_EPARAMS_SHIFT_IM LIS_EOPTIONS_LEN+2

#define LIS_OPTIONS_FILE -1
#define LIS_OPTIONS_HELP -2
#define LIS_OPTIONS_VER -3


#define LIS_PRINT_NONE 0
#define LIS_PRINT_MEM 1
#define LIS_PRINT_OUT 2
#define LIS_PRINT_ALL 3

#define LIS_EPRINT_NONE 0
#define LIS_EPRINT_MEM 1
#define LIS_EPRINT_OUT 2
#define LIS_EPRINT_ALL 3

#define LIS_SCALE_NONE 0
#define LIS_SCALE_JACOBI 1
#define LIS_SCALE_SYMM_DIAG 2

#define LIS_CONV_COND_DEFAULT 0
#define LIS_CONV_COND_NRM2_R 0
#define LIS_CONV_COND_NRM2_B 1
#define LIS_CONV_COND_NRM1_B 2

#define LIS_SOLVER_LEN 25
#define LIS_SOLVER_CG 1
#define LIS_SOLVER_BICG 2
#define LIS_SOLVER_CGS 3
#define LIS_SOLVER_BICGSTAB 4
#define LIS_SOLVER_BICGSTABL 5
#define LIS_SOLVER_GPBICG 6
#define LIS_SOLVER_QMR 7
#define LIS_SOLVER_TFQMR 7
#define LIS_SOLVER_ORTHOMIN 8
#define LIS_SOLVER_GMRES 9
#define LIS_SOLVER_JACOBI 10
#define LIS_SOLVER_GS 11
#define LIS_SOLVER_SOR 12
#define LIS_SOLVER_BICGSAFE 13
#define LIS_SOLVER_CR 14
#define LIS_SOLVER_BICR 15
#define LIS_SOLVER_CRS 16
#define LIS_SOLVER_BICRSTAB 17
#define LIS_SOLVER_GPBICR 18
#define LIS_SOLVER_BICRSAFE 19
#define LIS_SOLVER_FGMRES 20
#define LIS_SOLVER_IDRS 21
#define LIS_SOLVER_IDR1 22
#define LIS_SOLVER_MINRES 23
#define LIS_SOLVER_COCG 24
#define LIS_SOLVER_COCR 25

#define LIS_ESOLVER_LEN 16
#define LIS_ESOLVER_PI 1
#define LIS_ESOLVER_II 2
#define LIS_ESOLVER_RQI 3
#define LIS_ESOLVER_CG 4
#define LIS_ESOLVER_CR 5
#define LIS_ESOLVER_SI 6
#define LIS_ESOLVER_LI 7
#define LIS_ESOLVER_AI 8
#define LIS_ESOLVER_GPI 9
#define LIS_ESOLVER_GII 10
#define LIS_ESOLVER_GRQI 11
#define LIS_ESOLVER_GCG 12
#define LIS_ESOLVER_GCR 13
#define LIS_ESOLVER_GSI 14
#define LIS_ESOLVER_GLI 15
#define LIS_ESOLVER_GAI 16

#define LIS_INS_VALUE 0
#define LIS_ADD_VALUE 1
#define LIS_SUB_VALUE 2

#define LIS_MATRIX_LOWER 0
#define LIS_MATRIX_UPPER 1
#define LIS_MATRIX_SSOR 2

#define LIS_ORIGIN_0 0
#define LIS_ORIGIN_1 1

#define LIS_RESID 0
#define LIS_RANDOM 1

#define LIS_PRECISION_DEFAULT 0
#define LIS_PRECISION_DOUBLE 0
#define LIS_PRECISION_QUAD 1
#define LIS_PRECISION_SWITCH 2

#define LIS_LABEL_VECTOR 0
#define LIS_LABEL_MATRIX 1

#define LIS_VEC_TMP_PADD 128

#define LIS_VECTOR_NULL -1
#define LIS_VECTOR_ASSEMBLING 0
#define LIS_VECTOR_ASSEMBLED 1

#define LIS_PRECONNAME_MAX 10
#define LIS_PRECON_REGISTER_MAX 10

#define LIS_PRECON_TYPE_LEN 12
#define LIS_PRECON_TYPE_NONE 0
#define LIS_PRECON_TYPE_JACOBI 1
#define LIS_PRECON_TYPE_ILU 2
#define LIS_PRECON_TYPE_SSOR 3
#define LIS_PRECON_TYPE_HYBRID 4
#define LIS_PRECON_TYPE_IS 5
#define LIS_PRECON_TYPE_SAI 6
#define LIS_PRECON_TYPE_SAAMG 7
#define LIS_PRECON_TYPE_ILUC 8
#define LIS_PRECON_TYPE_ILUT 9
#define LIS_PRECON_TYPE_BJACOBI 10
#define LIS_PRECON_TYPE_ADDS 11
#define LIS_PRECON_TYPE_USERDEF LIS_PRECON_TYPE_LEN

#define LIS_MATRIX_ASSEMBLING 0
#define LIS_MATRIX_CSR 1
#define LIS_MATRIX_CSC 2
#define LIS_MATRIX_MSR 3
#define LIS_MATRIX_DIA 4
#define LIS_MATRIX_CDS 4
#define LIS_MATRIX_ELL 5
#define LIS_MATRIX_JAD 6
#define LIS_MATRIX_BSR 7
#define LIS_MATRIX_BSC 8
#define LIS_MATRIX_VBR 9
#define LIS_MATRIX_COO 10
#define LIS_MATRIX_DENSE 11
#define LIS_MATRIX_DNS 11
#define LIS_MATRIX_RCO 255

#define LIS_MATRIX_TJAD 12
#define LIS_MATRIX_BJAD 13
#define LIS_MATRIX_BCR 14
#define LIS_MATRIX_CJAD 15
#define LIS_MATRIX_PCSR 16
#define LIS_MATRIX_LCSR 17
#define LIS_MATRIX_LJAD 18
#define LIS_MATRIX_LBSR 19
#define LIS_MATRIX_CDIA 20
#define LIS_MATRIX_MSC 21
#define LIS_MATRIX_DECIDING_SIZE -(LIS_MATRIX_RCO+1)
#define LIS_MATRIX_NULL -(LIS_MATRIX_RCO+2)

#define LIS_MATRIX_DEFAULT LIS_MATRIX_CSR
#define LIS_MATRIX_POINT LIS_MATRIX_CSR
#define LIS_MATRIX_BLOCK LIS_MATRIX_BSR


#if defined(_DEBUG)
#define LIS_DEBUG_FUNC_IN lis_debug_trace_func(1,__FUNC__)
#define LIS_DEBUG_FUNC_OUT lis_debug_trace_func(0,__FUNC__)
#else
#define LIS_DEBUG_FUNC_IN
#define LIS_DEBUG_FUNC_OUT
#endif
/****************************************/
/****************************************/
typedef struct
{
	double	hi;
	double	lo;
} LIS_DOUBLE_DOUBLE;

typedef struct
{
	double	hi[2];
	double	lo[2];
} LIS_DOUBLE_DOUBLE_PD;

typedef struct
{
	double	*hi;
	double	*lo;
} LIS_DOUBLE_DOUBLE_PTR;

#if defined(_LONG__DOUBLE)
#ifdef HAVE_COMPLEX_H
typedef long double complex LIS_COMPLEX;
#else
typedef long double LIS_COMPLEX[2];
#endif
#ifdef _COMPLEX
#ifdef HAVE_QUADMATH_H
typedef __complex128 LIS_SCALAR;
typedef __float128 LIS_REAL;
#define creal(x) crealq(x)
#define cimag(x) cimagq(x)
#define acos(x) cacosq(x)
#define acosh(x) cacoshq(x)
#define asin(x) casinq(x)
#define asinh(x) casinhq(x)
#define atan(x) catanq(x)
#define atanh(x) catanhq(x)
#define cos(x) ccosq(x)
#define cosh(x) ccoshq(x)
#define exp(x) cexpq(x)
#define fabs(x) cabsq(x)
#define log(x) clogq(x)
#define conj(x) conjq(x)
#define pow(x,y) cpowq(x,y)
#define proj(x) cprojq(x)
#define sin(x) csinq(x)
#define sinh(x) csinhq(x)
#define sqrt(x) csqrtq(x)
#define tan(x) ctanq(x)
#define tanh(x) ctanhq(x)
#define LIS_MPI_SCALAR MPI_C_LONG_DOUBLE_COMPLEX
#define LIS_MPI_REAL MPI_LONG_DOUBLE
#else
typedef LIS_COMPLEX LIS_SCALAR;
typedef long double LIS_REAL;
#define creal(x) creall(x)
#define cimag(x) cimagl(x)
#define acos(x) cacosl(x)
#define acosh(x) cacoshl(x)
#define asin(x) casinl(x)
#define asinh(x) casinhl(x)
#define atan(x) catanl(x)
#define atanh(x) catanhl(x)
#define cos(x) ccosl(x)
#define cosh(x) ccoshl(x)
#define exp(x) cexpl(x)
#define fabs(x) cabsl(x)
#define log(x) clogl(x)
#define conj(x) conjl(x)
#define pow(x,y) cpowl(x,y)
#define proj(x) cprojl(x)
#define sin(x) csinl(x)
#define sinh(x) csinhl(x)
#define sqrt(x) csqrtl(x)
#define tan(x) ctanl(x)
#define tanh(x) ctanhl(x)
#define LIS_MPI_SCALAR MPI_C_LONG_DOUBLE_COMPLEX
#define LIS_MPI_REAL MPI_LONG_DOUBLE
#endif
#else
#ifdef HAVE_QUADMATH_H
typedef __float128 LIS_SCALAR;
typedef __float128 LIS_REAL;  
#define sin(x) sinq(x)
#define cos(x) cosq(x)
#define tan(x) tanq(x)
#define asin(x) asinq(x)
#define acos(x) acosq(x)
#define atan(x) atanq(x)
#define atan2(x,y) atan2q((x),(y))
#define fmod(x,y) fmodq((x),(y))
#define fabs(x) fabsq(x)
#define sqrt(x) sqrtq(x)
#define log(x) logq(x)
#define log10(x) log10q(x)
#define conj(x) x
#define exp(x) expq(x)
#define pow(x,y) powq((x),(y))
#define LIS_MPI_SCALAR MPI_LONG_DOUBLE
#define LIS_MPI_REAL MPI_LONG_DOUBLE
#else
typedef long double LIS_SCALAR;
typedef long double LIS_REAL;
#define sin(x) sinl(x)
#define cos(x) cosl(x)
#define tan(x) tanl(x)
#define asin(x) asinl(x)
#define acos(x) acosl(x)
#define atan(x) atanl(x)
#define atan2(x,y) atan2l((x),(y))
#define fmod(x,y) fmodl((x),(y))
#define fabs(x) fabsl(x)
#define sqrt(x) sqrtl(x)
#define log(x) logl(x)
#define log10(x) log10l(x)
#define conj(x) x
#define exp(x) expl(x)
#define pow(x,y) powl((x),(y))
#define LIS_MPI_SCALAR MPI_LONG_DOUBLE
#define LIS_MPI_REAL MPI_LONG_DOUBLE
#endif
#endif
#else
#ifdef HAVE_COMPLEX_H
typedef double complex LIS_COMPLEX;
#else
typedef double LIS_COMPLEX[2];
#endif
#ifdef _COMPLEX
typedef LIS_COMPLEX LIS_SCALAR;
typedef double LIS_REAL;
#define acos(x) cacos(x)
#define acosh(x) cacosh(x)
#define asin(x) casin(x)
#define asinh(x) casinh(x)
#define atan(x) catan(x)
#define atanh(x) catanh(x)
#define cos(x) ccos(x)
#define cosh(x) ccosh(x)
#define exp(x) cexp(x)
#define fabs(x) cabs(x)
#define log(x) clog(x)
#define pow(x,y) cpow(x,y)
#define proj(x) cproj(x)
#define sin(x) csin(x)
#define sinh(x) csinh(x)
#define sqrt(x) csqrt(x)
#define tan(x) ctan(x)
#define tanh(x) ctanh(x)
#define LIS_MPI_SCALAR MPI_C_DOUBLE_COMPLEX
#define LIS_MPI_REAL MPI_DOUBLE
#else
typedef double LIS_SCALAR;
typedef double LIS_REAL;
#define conj(x) x
#define LIS_MPI_SCALAR MPI_DOUBLE
#define LIS_MPI_REAL MPI_DOUBLE
#endif
#endif
typedef LIS_DOUBLE_DOUBLE LIS_QUAD;
typedef LIS_DOUBLE_DOUBLE_PD LIS_QUAD_PD;
typedef LIS_DOUBLE_DOUBLE_PTR LIS_QUAD_PTR;

#if defined(_LONG__LONG)
typedef long long int LIS_INT;
typedef unsigned long long int LIS_UNSIGNED_INT;
#else
typedef int LIS_INT;
typedef unsigned int LIS_UNSIGNED_INT;
#endif

#ifdef USE_MAIN__
	#define main MAIN__
#endif

#ifdef _OPENMP
	#include <omp.h>
#endif

#ifdef USE_MPI
	#include <mpi.h>
	typedef MPI_Comm LIS_Comm;
	#define LIS_COMM_WORLD	((LIS_Comm)MPI_COMM_WORLD)
#if defined(_LONG__LONG) 
	#define LIS_MPI_INT     MPI_LONG_LONG_INT
	#define LIS_MPI_Fint    MPI_Fint
#else
	#define LIS_MPI_INT     MPI_INT
	#define LIS_MPI_Fint    MPI_Fint
#endif
#else
	typedef LIS_INT LIS_Comm;
	#define LIS_COMM_WORLD	((LIS_Comm)0x1)
#endif

struct LIS_COMMTABLE_STRUCT
{
	LIS_Comm comm;
	LIS_INT pad;
	LIS_INT neibpetot;
	LIS_INT imnnz;
	LIS_INT exnnz;
	LIS_INT wssize;
	LIS_INT wrsize;
	LIS_INT *neibpe;
	LIS_INT *import_ptr;
	LIS_INT *import_index;
	LIS_INT *export_ptr;
	LIS_INT *export_index;
	LIS_SCALAR *ws;
	LIS_SCALAR *wr;
#ifdef USE_MPI
	MPI_Request  *req1,*req2;
	MPI_Status *sta1,*sta2;
#endif
};

typedef struct LIS_COMMTABLE_STRUCT *LIS_COMMTABLE;

struct LIS_VECTOR_STRUCT
{
	LIS_INT label;
	LIS_INT status;
	LIS_INT precision;
	LIS_INT gn;
	LIS_INT n;
	LIS_INT np;
	LIS_INT pad;
	LIS_INT origin;
	LIS_INT is_copy;
	LIS_INT is_destroy;
	LIS_INT is_scaled;
	LIS_INT my_rank;
	LIS_INT nprocs;
	LIS_Comm comm;
	LIS_INT is;
	LIS_INT ie;
	LIS_INT *ranges;
	LIS_SCALAR *value;
	LIS_SCALAR *value_lo;
	LIS_SCALAR *work;
	LIS_INT intvalue;
};
typedef struct LIS_VECTOR_STRUCT *LIS_VECTOR;

struct LIS_VECTOR_S_STRUCT
{
	LIS_INT label;
	LIS_INT status;
	LIS_INT precision;
	LIS_INT gn;
	LIS_INT n;
	LIS_INT np;
	LIS_INT pad;
	LIS_INT origin;
	LIS_INT is_copy;
	LIS_INT is_destroy;
	LIS_INT is_scaled;
	LIS_INT my_rank;
	LIS_INT nprocs;
	LIS_Comm comm;
	LIS_INT is;
	LIS_INT ie;
	LIS_INT *ranges;
	LIS_SCALAR *value;
	LIS_SCALAR *value_lo;
	LIS_QUAD_PTR value_hl;
	LIS_SCALAR *work;
	LIS_INT *index;
	LIS_INT nnz;
};
typedef struct LIS_VECTOR_S_STRUCT *LIS_VECTOR_S;

#define LIS_MATRIX_OPTION_LEN 10

struct LIS_MATRIX_CORE_STRUCT
{
	LIS_INT nnz;
	LIS_INT ndz;
	LIS_INT bnr;
	LIS_INT bnc;
	LIS_INT nr;
	LIS_INT nc;
	LIS_INT bnnz;
	LIS_INT nnd;
	LIS_INT maxnzr;
	LIS_INT *ptr;
	LIS_INT *row;
	LIS_INT *col;
	LIS_INT *index;
	LIS_INT *bptr;
	LIS_INT *bindex;
	LIS_SCALAR *value;
	LIS_SCALAR *work;
};
typedef struct LIS_MATRIX_CORE_STRUCT *LIS_MATRIX_CORE;

struct LIS_MATRIX_DIAG_STRUCT
{
	LIS_INT label;
	LIS_INT status;
	LIS_INT precision;
	LIS_INT gn;
	LIS_INT n;
	LIS_INT np;
	LIS_INT pad;
	LIS_INT origin;
	LIS_INT is_copy;
	LIS_INT is_destroy;
	LIS_INT is_scaled;
	LIS_INT my_rank;
	LIS_INT nprocs;
	LIS_Comm comm;
	LIS_INT is;
	LIS_INT ie;
	LIS_INT *ranges;
	LIS_SCALAR *value;
	LIS_SCALAR *work;

	LIS_INT bn;
	LIS_INT nr;
	LIS_INT *bns;
	LIS_INT *ptr;
	LIS_SCALAR **v_value;
};
typedef struct LIS_MATRIX_DIAG_STRUCT *LIS_MATRIX_DIAG;

struct LIS_MATRIX_STRUCT
{
	LIS_INT label;
	LIS_INT status;
	LIS_INT precision;
	LIS_INT gn;
	LIS_INT n;
	LIS_INT np;
	LIS_INT pad;
	LIS_INT origin;
	LIS_INT is_copy;
	LIS_INT is_destroy;
	LIS_INT is_scaled;
	LIS_INT my_rank;
	LIS_INT nprocs;
	LIS_Comm comm;
	LIS_INT is;
	LIS_INT ie;
	LIS_INT *ranges;

	LIS_INT matrix_type;
	LIS_INT nnz; /* CSR,CSC,MSR,JAD,VBR,COO */
	LIS_INT ndz; /* MSR */
	LIS_INT bnr; /* BSR,BSC */
	LIS_INT bnc; /* BSR,BSC */
	LIS_INT nr; /* BSR,BSC,VBR */
	LIS_INT nc; /* BSR,BSC,VBR */
	LIS_INT bnnz; /* BSR,BSC,VBR */
	LIS_INT nnd; /* DIA */
	LIS_INT maxnzr; /* ELL,JAD */
	LIS_INT *ptr; /* CSR,CSC,JAD */
	LIS_INT *row; /* JAD,VBR,COO */
	LIS_INT *col; /* JAD,VBR,COO */
	LIS_INT *index; /* CSR,CSC,MSR,DIA,ELL,JAD */
	LIS_INT *bptr; /* BSR,BSC,VBR */
	LIS_INT *bindex; /* BSR,BSC,VBR */
	LIS_SCALAR *value; /* CSR,CSC,MSR,DIA,ELL,JAD,BSR,BSC,VBR,DNS,COO */
	LIS_SCALAR *work;

	LIS_MATRIX_CORE L;
	LIS_MATRIX_CORE U;
	LIS_MATRIX_DIAG D;
	LIS_MATRIX_DIAG WD;

	LIS_INT is_block;
	LIS_INT pad_comm;
	LIS_INT is_pmat;
	LIS_INT is_sorted;
	LIS_INT is_splited;
	LIS_INT is_save;
	LIS_INT is_comm;
	LIS_INT is_fallocated;
	LIS_INT use_wd;
	LIS_INT conv_bnr;
	LIS_INT conv_bnc;
	LIS_INT *conv_row;
	LIS_INT *conv_col;
	LIS_INT options[LIS_MATRIX_OPTION_LEN];

	LIS_INT w_annz;
	LIS_INT *w_nnz;
	LIS_INT *w_row;
	LIS_INT **w_index;
	LIS_SCALAR **w_value;
	LIS_SCALAR ***v_value;

	LIS_INT *l2g_map;
	LIS_COMMTABLE commtable;
};
typedef struct LIS_MATRIX_STRUCT *LIS_MATRIX;


struct LIS_MATRIX_ILU_STRUCT
{
	LIS_INT n;
	LIS_INT bs;
	LIS_INT *nnz_ma;
	LIS_INT *nnz;
	LIS_INT *bsz;
	LIS_INT **index;
	LIS_SCALAR **value;
	LIS_SCALAR ***values;
};
typedef struct LIS_MATRIX_ILU_STRUCT *LIS_MATRIX_ILU;

struct LIS_PRECON_STRUCT
{
	LIS_INT precon_type;
	LIS_MATRIX A; /* SSOR */
	LIS_MATRIX Ah;
	LIS_MATRIX_ILU L; /* ilu(k),ilut,iluc,sainv */
	LIS_MATRIX_ILU U; /* ilu(k),ilut,iluc,sainv */
	LIS_MATRIX_DIAG WD; /* bilu(k),bilut,biluc,bjacobi */
	LIS_VECTOR D; /* ilu(k),ilut,iluc,jacobi,sainv */
	LIS_VECTOR Pb; /* i+s */
	LIS_VECTOR temp; /* saamg */
	LIS_REAL theta; /* saamg */
	LIS_VECTOR *work; /* adds */
	struct LIS_SOLVER_STRUCT *solver; /* hybrid */
	LIS_INT worklen; /* adds */
	LIS_INT level_num; /* saamg */
	LIS_INT wsize; /* saamg */
	LIS_INT solver_comm; /* saamg */
	LIS_INT my_rank; /* saamg */
	LIS_INT nprocs; /* saamg */
	LIS_INT is_copy;
	LIS_COMMTABLE commtable; /* saamg */
};
typedef struct LIS_PRECON_STRUCT *LIS_PRECON;


struct LIS_SOLVER_STRUCT
{
	LIS_MATRIX A,Ah;
	LIS_VECTOR b,x,xx,d;
	LIS_MATRIX_DIAG WD;
	LIS_PRECON precon;
	LIS_VECTOR *work;
	LIS_REAL *rhistory;
	LIS_INT worklen;
	LIS_INT options[LIS_OPTIONS_LEN];
	LIS_SCALAR params[LIS_PARAMS_LEN];
	LIS_INT retcode;
	LIS_INT iter;
	LIS_INT iter2;
	LIS_REAL resid;
	double time;
	double itime;
	double ptime;
	double p_c_time;
	double p_i_time;
	LIS_INT precision;
	LIS_REAL bnrm;
	LIS_REAL tol;
	LIS_REAL tol_switch;
	LIS_INT setup;
};
typedef struct LIS_SOLVER_STRUCT *LIS_SOLVER;

struct LIS_ESOLVER_STRUCT
{
        LIS_MATRIX A,B;
        LIS_VECTOR x,xx,d; 
        LIS_SCALAR *evalue;
        LIS_VECTOR *evector;
	LIS_REAL *resid;
	LIS_VECTOR *work;
	LIS_REAL *rhistory;
	LIS_INT worklen;
	LIS_INT options[LIS_EOPTIONS_LEN];
	LIS_SCALAR params[LIS_EPARAMS_LEN];
	LIS_INT retcode;
        LIS_INT *iter;       
	LIS_INT *iter2;
        double time;
        LIS_INT *nesol;
        double itime;
        double ptime;
        double p_c_time;
        double p_i_time;
	LIS_INT eprecision;
        LIS_SCALAR ishift;
	LIS_REAL nrm2;
        LIS_REAL tol;
};
typedef struct LIS_ESOLVER_STRUCT *LIS_ESOLVER;

struct LIS_CONV_OPTIONS_STRUCT
{
	LIS_INT bnr;
	LIS_INT bnc;
	LIS_INT *row;
	LIS_INT *col;
};
typedef struct LIS_CONV_OPTIONS_STRUCT LIS_CONV_OPTIONS;

typedef LIS_INT (*LIS_PRECON_CREATE_XXX)(LIS_SOLVER solver, LIS_PRECON precon);
/*NEH support for extended "solve_kernel" workflow*/
typedef LIS_INT (*LIS_PRECON_PSD_CREATE_XXX)(LIS_SOLVER solver, LIS_PRECON precon);
/*NEH support for extended "solve_kernel" workflow*/
typedef LIS_INT (*LIS_PRECON_PSD_UPDATE_XXX)(LIS_SOLVER solver, LIS_PRECON precon);
typedef LIS_INT (*LIS_PSOLVE_XXX)(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
typedef LIS_INT (*LIS_PSOLVEH_XXX)(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);

typedef struct LIS_PRECON_REGISTER_STRUCT
{
	LIS_INT precon_type;
	char name[LIS_PRECONNAME_MAX+1];
	LIS_PRECON_CREATE_XXX pcreate;
	LIS_PSOLVE_XXX psolve;
	LIS_PSOLVEH_XXX psolveh;
} LIS_PRECON_REGISTER;


#ifdef __cplusplus
extern "C"
{
#endif

/****************************/
/* Vector Operations        */
/****************************/

	extern LIS_INT lis_vector_create(LIS_Comm comm, LIS_VECTOR *vec); 
	extern LIS_INT lis_vector_set_size(LIS_VECTOR vec, LIS_INT local_n, LIS_INT global_n); 
/*NEH support for extended "solve_kernel" workflow*/
	extern LIS_INT lis_vector_psd_reset_scale(LIS_VECTOR vec); 
	extern LIS_INT lis_vector_destroy(LIS_VECTOR vec);
	extern LIS_INT lis_vector_duplicate(void *vin, LIS_VECTOR *vout);
	extern LIS_INT lis_vector_get_size(LIS_VECTOR v, LIS_INT *local_n, LIS_INT *global_n);
	extern LIS_INT lis_vector_get_range(LIS_VECTOR v, LIS_INT *is, LIS_INT *ie);
	extern LIS_INT lis_vector_get_value(LIS_VECTOR v, LIS_INT i, LIS_SCALAR *value);
	extern LIS_INT lis_vector_get_values(LIS_VECTOR v, LIS_INT start, LIS_INT count, LIS_SCALAR value[]);
	extern LIS_INT lis_vector_set_value(LIS_INT flag, LIS_INT i, LIS_SCALAR value, LIS_VECTOR v);
	extern LIS_INT lis_vector_set_values(LIS_INT flag, LIS_INT count, LIS_INT index[], LIS_SCALAR value[], LIS_VECTOR v);
	extern LIS_INT lis_vector_set_values2(LIS_INT flag, LIS_INT start, LIS_INT count, LIS_SCALAR value[], LIS_VECTOR v);
	extern LIS_INT lis_vector_print(LIS_VECTOR x);
	extern LIS_INT lis_vector_scatter(LIS_SCALAR value[], LIS_VECTOR v);
        extern LIS_INT lis_vector_gather(LIS_VECTOR v, LIS_SCALAR value[]);
	extern LIS_INT lis_vector_is_null(LIS_VECTOR v);
	extern LIS_INT lis_vector_swap(LIS_VECTOR vsrc, LIS_VECTOR vdst);
	extern LIS_INT lis_vector_copy(LIS_VECTOR vsrc, LIS_VECTOR vdst);
	extern LIS_INT lis_vector_axpy(LIS_SCALAR alpha, LIS_VECTOR vx, LIS_VECTOR vy);
	extern LIS_INT lis_vector_xpay(LIS_VECTOR vx, LIS_SCALAR alpha, LIS_VECTOR vy);
	extern LIS_INT lis_vector_axpyz(LIS_SCALAR alpha, LIS_VECTOR vx, LIS_VECTOR vy, LIS_VECTOR vz);
	extern LIS_INT lis_vector_scale(LIS_SCALAR alpha, LIS_VECTOR vx);
	extern LIS_INT lis_vector_pmul(LIS_VECTOR vx,LIS_VECTOR vy,LIS_VECTOR vz);
	extern LIS_INT lis_vector_pdiv(LIS_VECTOR vx,LIS_VECTOR vy,LIS_VECTOR vz);
	extern LIS_INT lis_vector_set_all(LIS_SCALAR alpha, LIS_VECTOR vx);
	extern LIS_INT lis_vector_abs(LIS_VECTOR vx);
	extern LIS_INT lis_vector_reciprocal(LIS_VECTOR vx);
	extern LIS_INT lis_vector_conjugate(LIS_VECTOR vx);
	extern LIS_INT lis_vector_shift(LIS_SCALAR sigma, LIS_VECTOR vx);
	extern LIS_INT lis_vector_dot(LIS_VECTOR vx, LIS_VECTOR vy, LIS_SCALAR *value);
	extern LIS_INT lis_vector_nhdot(LIS_VECTOR vx, LIS_VECTOR vy, LIS_SCALAR *value);  
	extern LIS_INT lis_vector_nrm1(LIS_VECTOR vx, LIS_REAL *value);
	extern LIS_INT lis_vector_nrm2(LIS_VECTOR vx, LIS_REAL *value);
	extern LIS_INT lis_vector_nrmi(LIS_VECTOR vx, LIS_REAL *value);
	extern LIS_INT lis_vector_sum(LIS_VECTOR vx, LIS_SCALAR *value);

/****************************/
/* Matrix Operations        */
/****************************/

	extern LIS_INT lis_matrix_create(LIS_Comm comm, LIS_MATRIX *Amat);
	extern LIS_INT lis_matrix_destroy(LIS_MATRIX Amat);
	extern LIS_INT lis_matrix_assemble(LIS_MATRIX A);
	extern LIS_INT lis_matrix_is_assembled(LIS_MATRIX A);
	extern LIS_INT lis_matrix_duplicate(LIS_MATRIX Ain, LIS_MATRIX *Aout);
	extern LIS_INT lis_matrix_set_size(LIS_MATRIX A, LIS_INT local_n, LIS_INT global_n);
	extern LIS_INT lis_matrix_get_size(LIS_MATRIX A, LIS_INT *local_n, LIS_INT *global_n);
	extern LIS_INT lis_matrix_get_range(LIS_MATRIX A, LIS_INT *is, LIS_INT *ie);
	extern LIS_INT lis_matrix_get_nnz(LIS_MATRIX A, LIS_INT *nnz);
	extern LIS_INT lis_matrix_set_type(LIS_MATRIX A, LIS_INT matrix_type);
	extern LIS_INT lis_matrix_get_type(LIS_MATRIX A, LIS_INT *matrix_type);
	extern LIS_INT lis_matrix_set_value(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_set_value_new(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_set_values(LIS_INT flag, LIS_INT n, LIS_SCALAR value[], LIS_MATRIX A);
	extern LIS_INT lis_matrix_psd_set_value(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A);  
	extern LIS_INT lis_matrix_malloc(LIS_MATRIX A, LIS_INT nnz_row, LIS_INT nnz[]);
	extern LIS_INT lis_matrix_get_diagonal(LIS_MATRIX A, LIS_VECTOR d);
	extern LIS_INT lis_matrix_scale(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR D, LIS_INT action);
/*NEH support for extended "solve_kernel" workflow*/
	extern LIS_INT lis_matrix_psd_reset_scale(LIS_MATRIX A);
	extern LIS_INT lis_matrix_convert(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_copy(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern LIS_INT lis_matrix_set_blocksize(LIS_MATRIX A, LIS_INT bnr, LIS_INT bnc, LIS_INT row[], LIS_INT col[]);
	extern LIS_INT lis_matrix_unset(LIS_MATRIX A);  

	extern LIS_INT lis_matrix_malloc_csr(LIS_INT n, LIS_INT nnz, LIS_INT **ptr, LIS_INT **index, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_value_csr(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_psd_set_value_csr(LIS_INT flag, LIS_INT i, LIS_INT j, LIS_SCALAR value, LIS_MATRIX A);  
	extern LIS_INT lis_matrix_set_csr(LIS_INT nnz, LIS_INT *row, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_csc(LIS_INT n, LIS_INT nnz, LIS_INT **ptr, LIS_INT **index, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_csc(LIS_INT nnz, LIS_INT *row, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_bsr(LIS_INT n, LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT **bptr, LIS_INT **bindex, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_bsr(LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_set_msr(LIS_INT nnz, LIS_INT ndz, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_msr(LIS_INT n, LIS_INT nnz, LIS_INT ndz, LIS_INT **index, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_ell(LIS_INT maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_ell(LIS_INT n, LIS_INT maxnzr, LIS_INT **index, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_jad(LIS_INT nnz, LIS_INT maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_jad(LIS_INT n, LIS_INT nnz, LIS_INT maxnzr, LIS_INT **perm, LIS_INT **ptr, LIS_INT **index, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_dia(LIS_INT nnd, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_dia(LIS_INT n, LIS_INT nnd, LIS_INT **index, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_malloc_bsc(LIS_INT n, LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT **bptr, LIS_INT **bindex, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_bsc(LIS_INT bnr, LIS_INT bnc, LIS_INT bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_get_vbr_rowcol(LIS_MATRIX Ain, LIS_INT *nr, LIS_INT *nc, LIS_INT **row, LIS_INT **col);
	extern LIS_INT lis_matrix_malloc_vbr(LIS_INT n, LIS_INT nnz, LIS_INT nr, LIS_INT nc, LIS_INT bnnz, LIS_INT **row, LIS_INT **col, LIS_INT **ptr, LIS_INT **bptr, LIS_INT **bindex, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_vbr(LIS_INT nnz, LIS_INT nr, LIS_INT nc, LIS_INT bnnz, LIS_INT *row, LIS_INT *col, LIS_INT *ptr, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_coo(LIS_INT nnz, LIS_INT **row, LIS_INT **col, LIS_SCALAR **value);
	extern LIS_INT lis_matrix_set_coo(LIS_INT nnz, LIS_INT *row, LIS_INT *col, LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_set_dns(LIS_SCALAR *value, LIS_MATRIX A);
	extern LIS_INT lis_matrix_malloc_dns(LIS_INT n, LIS_INT gn, LIS_SCALAR **value);

/****************************/
/* Matrix-Vector Operations */
/****************************/

	extern LIS_INT lis_matvec(LIS_MATRIX A, LIS_VECTOR x, LIS_VECTOR y);
	extern LIS_INT lis_matvech(LIS_MATRIX A, LIS_VECTOR x, LIS_VECTOR y);  

/****************************/
/* Array Operations         */
/****************************/

        extern LIS_INT lis_array_swap(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y);
        extern LIS_INT lis_array_copy(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y);
        extern LIS_INT lis_array_axpy(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x, LIS_SCALAR *y);
        extern LIS_INT lis_array_xpay(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR alpha, LIS_SCALAR *y);
        extern LIS_INT lis_array_axpyz(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z);
        extern LIS_INT lis_array_scale(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x);
        extern LIS_INT lis_array_pmul(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z);
        extern LIS_INT lis_array_pdiv(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z);
        extern LIS_INT lis_array_set_all(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x);
        extern LIS_INT lis_array_abs(LIS_INT n, LIS_SCALAR *x);
        extern LIS_INT lis_array_reciprocal(LIS_INT n, LIS_SCALAR *x);
        extern LIS_INT lis_array_conjugate(LIS_INT n, LIS_SCALAR *x);  
        extern LIS_INT lis_array_shift(LIS_INT n, LIS_SCALAR sigma, LIS_SCALAR *x);
        extern LIS_INT lis_array_dot(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *value);
        extern LIS_INT lis_array_nhdot(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *value);  
	extern LIS_INT lis_array_nrm1(LIS_INT n, LIS_SCALAR *x, LIS_REAL *value);
	extern LIS_INT lis_array_nrm2(LIS_INT n, LIS_SCALAR *x, LIS_REAL *value);
	extern LIS_INT lis_array_nrmi(LIS_INT n, LIS_SCALAR *x, LIS_REAL *value);
	extern LIS_INT lis_array_sum(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *value);
	extern LIS_INT lis_array_matvec(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT op);
	extern LIS_INT lis_array_matvech(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT op);
	extern LIS_INT lis_array_matvec_ns(LIS_INT m, LIS_INT n, LIS_SCALAR *a, LIS_INT lda, LIS_SCALAR *x, LIS_SCALAR *y, LIS_INT op);
	extern LIS_INT lis_array_matmat(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, LIS_INT op);
	extern LIS_INT lis_array_matmat_ns(LIS_INT m, LIS_INT n, LIS_INT k, LIS_SCALAR *a, LIS_INT lda, LIS_SCALAR *b, LIS_INT ldb, LIS_SCALAR *c, LIS_INT ldc, LIS_INT op);
	extern LIS_INT lis_array_ge(LIS_INT n, LIS_SCALAR *a);
	extern LIS_INT lis_array_solve(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *x, LIS_SCALAR *w);
        extern LIS_INT lis_array_cgs(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r);
        extern LIS_INT lis_array_mgs(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r);
        extern LIS_INT lis_array_qr(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *q, LIS_SCALAR *r, LIS_INT *qriter, LIS_REAL *qrerr);

/****************************/
/* Linear Solvers           */
/****************************/

    extern LIS_INT lis_solver_create(LIS_SOLVER *solver);
	extern LIS_INT lis_solver_destroy(LIS_SOLVER solver);
	extern LIS_INT lis_solver_get_iter(LIS_SOLVER solver, LIS_INT *iter);
	extern LIS_INT lis_solver_get_iterex(LIS_SOLVER solver, LIS_INT *iter, LIS_INT *iter_double, LIS_INT *iter_quad);
	extern LIS_INT lis_solver_get_time(LIS_SOLVER solver, double *time);
	extern LIS_INT lis_solver_get_timeex(LIS_SOLVER solver, double *time, double *itime, double *ptime, double *p_c_time, double *p_i_time);
	extern LIS_INT lis_solver_get_residualnorm(LIS_SOLVER solver, LIS_REAL *residual);
	extern LIS_INT lis_solver_get_solver(LIS_SOLVER solver, LIS_INT *nsol);
	extern LIS_INT lis_solver_get_precon(LIS_SOLVER solver, LIS_INT *precon_type);
	extern LIS_INT lis_solver_get_status(LIS_SOLVER solver, LIS_INT *status);
        extern LIS_INT lis_solver_get_rhistory(LIS_SOLVER solver, LIS_VECTOR v);
	extern LIS_INT lis_solver_set_option(char *text, LIS_SOLVER solver);
	extern LIS_INT lis_solver_set_optionC(LIS_SOLVER solver);
/*NEH support for extended "solve_kernel" workflow*/
    extern LIS_INT lis_solver_set_matrix(LIS_MATRIX A, LIS_SOLVER solver);
	extern LIS_INT lis_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_SOLVER solver);
	extern LIS_INT lis_solve_setup(LIS_MATRIX A, LIS_SOLVER solver);
	extern LIS_INT lis_solve_kernel(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_SOLVER solver, LIS_PRECON precon);
	extern LIS_PRECON_REGISTER *precon_register_top;
	extern LIS_INT precon_register_type;
	extern LIS_INT lis_precon_register(char *name, LIS_PRECON_CREATE_XXX pcreate, LIS_PSOLVE_XXX psolve, LIS_PSOLVEH_XXX psolveh);
	extern LIS_INT lis_precon_register_free(void);
	extern LIS_INT lis_solver_get_solvername(LIS_INT solver, char *solvername);
	extern LIS_INT lis_solver_get_preconname(LIS_INT precon_type, char *preconname);

/****************************/
/* Eigensolvers             */
/****************************/

	extern LIS_INT lis_esolver_create(LIS_ESOLVER *esolver);
	extern LIS_INT lis_esolver_destroy(LIS_ESOLVER esolver);
	extern LIS_INT lis_iesolver_destroy(LIS_ESOLVER esolver);
	extern LIS_INT lis_esolver_set_option(char *text, LIS_ESOLVER esolver);
	extern LIS_INT lis_esolver_set_optionC(LIS_ESOLVER esolver);
	extern LIS_INT lis_esolve(LIS_MATRIX A, LIS_VECTOR x, LIS_SCALAR *evalue0, LIS_ESOLVER esolver);
	extern LIS_INT lis_gesolve(LIS_MATRIX A, LIS_MATRIX B, LIS_VECTOR x, LIS_SCALAR *evalue0, LIS_ESOLVER esolver);  
	extern LIS_INT lis_esolver_get_iter(LIS_ESOLVER esolver, LIS_INT *iter);
	extern LIS_INT lis_esolver_get_iterex(LIS_ESOLVER esolver, LIS_INT *iter, LIS_INT *iter_double, LIS_INT *iter_quad);
	extern LIS_INT lis_esolver_get_time(LIS_ESOLVER esolver, double *time);
	extern LIS_INT lis_esolver_get_timeex(LIS_ESOLVER esolver, double *time, double *itime, double *ptime, double *p_c_time, double *p_i_time);
	extern LIS_INT lis_esolver_get_residualnorm(LIS_ESOLVER esolver, LIS_REAL *residual);
	extern LIS_INT lis_esolver_get_status(LIS_ESOLVER esolver, LIS_INT *status);
        extern LIS_INT lis_esolver_get_rhistory(LIS_ESOLVER esolver, LIS_VECTOR v);
        extern LIS_INT lis_esolver_get_evalues(LIS_ESOLVER esolver, LIS_VECTOR v);
        extern LIS_INT lis_esolver_get_specific_evalue(LIS_ESOLVER esolver, LIS_INT mode, LIS_SCALAR *evalue);  
        extern LIS_INT lis_esolver_get_evectors(LIS_ESOLVER esolver, LIS_MATRIX M);
        extern LIS_INT lis_esolver_get_specific_evector(LIS_ESOLVER esolver, LIS_INT mode, LIS_VECTOR x);  
        extern LIS_INT lis_esolver_get_residualnorms(LIS_ESOLVER esolver, LIS_VECTOR v);
        extern LIS_INT lis_esolver_get_specific_residualnorm(LIS_ESOLVER esolver, LIS_INT mode, LIS_REAL *residual);  
        extern LIS_INT lis_esolver_get_iters(LIS_ESOLVER esolver, LIS_VECTOR v);
        extern LIS_INT lis_esolver_get_specific_iter(LIS_ESOLVER esolver, LIS_INT mode, LIS_INT *iter);  
	extern LIS_INT lis_esolver_get_esolver(LIS_ESOLVER esolver, LIS_INT *nesol);
	extern LIS_INT lis_esolver_get_esolvername(LIS_INT esolver, char *esolvername);

/****************************/
/* I/O Functions            */
/****************************/

	extern LIS_INT lis_input(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, char *filename);
	extern LIS_INT lis_input_matrix(LIS_MATRIX A, char *filename);
	extern LIS_INT lis_input_vector(LIS_VECTOR v, char *filename);
	extern LIS_INT lis_output(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT mode, char *path);
	extern LIS_INT lis_output_matrix(LIS_MATRIX A, LIS_INT mode, char *path);
	extern LIS_INT lis_output_vector(LIS_VECTOR v, LIS_INT format, char *filename);
	extern LIS_INT lis_solver_output_rhistory(LIS_SOLVER solver, char *filename);
	extern LIS_INT lis_esolver_output_rhistory(LIS_ESOLVER esolver, char *filename);

/****************************/
/* Utilities                */
/****************************/

	extern LIS_INT lis_initialize(int *argc, char** argv[]);
	extern LIS_INT lis_finalize(void);
	extern double lis_wtime(void);
	extern void CHKERR(LIS_INT err);
	extern LIS_INT lis_debug_trace_func(LIS_INT flag, char *func);
	extern void *lis_malloc( size_t size, char *tag );
	extern void *lis_calloc( size_t size, char *tag );
	extern void *lis_realloc( void *p, size_t size );
	extern void lis_free(void *p);
	extern void lis_free2(LIS_INT n, ...);
	extern LIS_INT lis_is_malloc( void *p );
	extern void lis_date(char *date);
	extern LIS_INT lis_printf(LIS_Comm comm, const char *mess, ...);


#ifdef __cplusplus
}
#endif

#define LIS_TRUE 1
#define LIS_FALSE 0
#define LIS_FAILS -1
#define LIS_SUCCESS 0
#define LIS_ILL_OPTION 1
#define LIS_ERR_ILL_ARG 1
#define LIS_BREAKDOWN 2
#define LIS_OUT_OF_MEMORY 3
#define LIS_ERR_OUT_OF_MEMORY 3
#define LIS_MAXITER 4
#define LIS_ERR_NOT_IMPLEMENTED	5
#define LIS_ERR_FILE_IO	6


#if 1
#define LIS_GET_ISIE(id,nprocs,n,is,ie) \
		if( (id) < (n)%(nprocs) ) \
		{ \
			(ie) = (n)/(nprocs)+1; \
			(is) = (ie)*(id); \
		} \
		else \
		{ \
			(ie) = (n)/(nprocs); \
			(is) = (ie)*(id) + (n)%(nprocs); \
		} \
		(ie) = (ie)+(is);
#else
#define LIS_GET_ISIE(id,nprocs,n,is,ie) \
		(ie) = (n)/(nprocs)+1; \
		(is) = (id)*(ie)<((n)+1)?(id)*(ie):(n)+1; \
		(ie) = (is)+(ie)-1<(n)?(is)+(ie)-1:(n);
#endif

#ifdef __cplusplus
}
#endif

#endif
