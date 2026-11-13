/*
 * Optional half precision preconditioner storage support.
 *
 * Binary16 conversion uses IEEE 754 round-to-nearest, ties-to-even.
 * Solver arithmetic remains in the configured LIS scalar precision.
 */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__GNUC__) && \
	(defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define LIS_HALF_X86_F16C_DISPATCH 1
#else
#define LIS_HALF_X86_F16C_DISPATCH 0
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lislib.h"


static uint64_t
lis_half_double_bits(double value)
{
	uint64_t bits;

	memcpy(&bits,&value,sizeof(bits));
	return bits;
}


static double
lis_half_double_from_bits(uint64_t bits)
{
	double value;

	memcpy(&value,&bits,sizeof(value));
	return value;
}


uint16_t
lis_half_from_double(double value)
{
	uint64_t bits;
	uint64_t frac;
	uint64_t mant;
	uint64_t kept;
	uint64_t rem;
	uint64_t halfway;
	uint16_t sign;
	uint16_t payload;
	unsigned int exp;
	int e;
	int shift;

	bits = lis_half_double_bits(value);

	sign = (uint16_t)((bits >> 48) & 0x8000u);
	exp = (unsigned int)((bits >> 52) & 0x7ffu);
	frac = bits & UINT64_C(0x000fffffffffffff);

	/* Infinity and NaN. */
	if( exp==0x7ffu )
	{
		if( frac==0 )
		{
			return (uint16_t)(sign | 0x7c00u);
		}

		payload = (uint16_t)(frac >> 42);

		/*
		 * Quiet NaNs on conversion while retaining the upper payload.
		 */
		payload = (uint16_t)((payload | 0x0200u) & 0x03ffu);

		return (uint16_t)(sign | 0x7c00u | payload);
	}

	/*
	 * Zero and every binary64 subnormal are far below the binary16
	 * representable range.
	 */
	if( exp==0 )
	{
		return sign;
	}

	e = (int)exp - 1023;
	mant = UINT64_C(0x0010000000000000) | frac;

	/* Definite overflow. */
	if( e>15 )
	{
		return (uint16_t)(sign | 0x7c00u);
	}

	/*
	 * Binary16 normal result.
	 *
	 * Reduce the binary64 53-bit significand to the binary16
	 * 11-bit significand and apply round-to-nearest, ties-to-even.
	 */
	if( e>=-14 )
	{
		kept = mant >> 42;
		rem = mant & ((UINT64_C(1) << 42) - 1);
		halfway = UINT64_C(1) << 41;

		if( rem>halfway ||
			(rem==halfway && (kept & 1u)) )
		{
			kept++;
		}

		if( kept==2048u )
		{
			kept = 1024u;
			e++;

			if( e>15 )
			{
				return (uint16_t)(sign | 0x7c00u);
			}
		}

		return (uint16_t)(
			sign |
			((uint16_t)(e + 15) << 10) |
			(uint16_t)(kept - 1024u));
	}

	/*
	 * Values below 2^-25 cannot round to a nonzero binary16 value.
	 * Exactly 2^-25 is a tie between zero and the minimum subnormal;
	 * ties-to-even selects zero.
	 */
	if( e < -25 )
	{
		return sign;
	}

	/*
	 * Binary16 subnormal:
	 *
	 *     value = integer * 2^-24
	 *
	 * so round the binary64 significand directly to that integer.
	 */
	shift = 28 - e;

	kept = mant >> shift;
	rem = mant & ((UINT64_C(1) << shift) - 1);
	halfway = UINT64_C(1) << (shift - 1);

	if( rem>halfway ||
		(rem==halfway && (kept & 1u)) )
	{
		kept++;
	}

	/* Rounding may promote the largest subnormal to min normal. */
	if( kept>=1024u )
	{
		return (uint16_t)(sign | 0x0400u);
	}

	return (uint16_t)(sign | (uint16_t)kept);
}


double
lis_half_to_double(uint16_t value)
{
	uint64_t bits;
	uint64_t sign;
	uint64_t frac;
	unsigned int exp;
	double result;

	sign = ((uint64_t)(value & 0x8000u)) << 48;
	exp = (unsigned int)((value >> 10) & 0x1fu);
	frac = (uint64_t)(value & 0x03ffu);

	if( exp==0 )
	{
		if( frac==0 )
		{
			return lis_half_double_from_bits(sign);
		}

		result = ldexp((double)frac,-24);

		return sign ? -result : result;
	}

	if( exp==0x1fu )
	{
		bits =
			sign |
			UINT64_C(0x7ff0000000000000) |
			(frac << 42);

		return lis_half_double_from_bits(bits);
	}

	bits =
		sign |
		((uint64_t)(exp + 1008u) << 52) |
		(frac << 42);

	return lis_half_double_from_bits(bits);
}



struct LIS_HALF_ILU_MATRIX_STRUCT
{
	LIS_INT n;

	/*
	 * Non-owning symbolic ILU pattern.
	 * Only numerical factor values are stored in binary16.
	 */
	LIS_MATRIX_ILU pattern;

	uint16_t **value;
	int16_t *scale;
};

struct LIS_HALF_ILU_STRUCT
{
	struct LIS_HALF_ILU_MATRIX_STRUCT L;
	struct LIS_HALF_ILU_MATRIX_STRUCT U;
	LIS_REAL underflow_ratio;
};


static void
lis_half_matrix_destroy(
	struct LIS_HALF_ILU_MATRIX_STRUCT *A)
{
	LIS_INT i;

	if( A==NULL )
	{
		return;
	}

	if( A->value!=NULL )
	{
		for(i=0;i<A->n;i++)
		{
			if( A->value[i]!=NULL )
			{
				lis_free(A->value[i]);
			}
		}

		lis_free(A->value);
	}

	if( A->scale!=NULL )
	{
		lis_free(A->scale);
	}

	memset(A,0,sizeof(*A));
}


#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

#undef __FUNC__
#define __FUNC__ "lis_half_matrix_compress"
static LIS_INT
lis_half_matrix_compress(
	LIS_MATRIX_ILU A,
	struct LIS_HALF_ILU_MATRIX_STRUCT *H,
	LIS_REAL *worst_underflow)
{
	LIS_INT i,j,n,nnz;
	uint16_t *value;

	double v;
	double a;
	double maxabs;
	double scaled;
	double decoded;
	double row_total;
	double row_lost;
	double row_ratio;

	int exponent;

	if( A==NULL || H==NULL ||
		worst_underflow==NULL ||
		A->bs!=1 ||
		A->n<0 ||
		A->nnz==NULL ||
		A->index==NULL ||
		A->value==NULL )
	{
		LIS_SETERR(
			LIS_ERR_ILL_ARG,
			"invalid ILU factors for half storage\n");

		return LIS_ERR_ILL_ARG;
	}

	memset(H,0,sizeof(*H));

	n = A->n;

	H->n = n;
	H->pattern = A;

	if( n==0 )
	{
		return LIS_SUCCESS;
	}

	H->value =
		(uint16_t **)lis_malloc(
			n*sizeof(uint16_t *),
			"lis_half_matrix_compress::value");

	H->scale =
		(int16_t *)lis_malloc(
			n*sizeof(int16_t),
			"lis_half_matrix_compress::scale");

	if( H->value==NULL ||
		H->scale==NULL )
	{
		lis_half_matrix_destroy(H);

		LIS_SETERR_MEM(
			n*sizeof(uint16_t *));

		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<n;i++)
	{
		H->value[i] = NULL;
		H->scale[i] = 0;
	}

	for(i=0;i<n;i++)
	{
		nnz = A->nnz[i];

		if( nnz<0 ||
			(nnz>0 &&
			 (A->index[i]==NULL ||
			  A->value[i]==NULL)) )
		{
			lis_half_matrix_destroy(H);

			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"invalid ILU row for half storage\n");

			return LIS_ERR_ILL_ARG;
		}

		if( nnz==0 )
		{
			continue;
		}

		maxabs = 0.0;

		for(j=0;j<nnz;j++)
		{
			v = (double)A->value[i][j];

			if( !isfinite(v) )
			{
				lis_half_matrix_destroy(H);

				LIS_SETERR(
					LIS_ERR_ILL_ARG,
					"non-finite ILU factor cannot use half storage\n");

				return LIS_ERR_ILL_ARG;
			}

			a = fabs(v);

			if( a>maxabs )
			{
				maxabs = a;
			}
		}

		value =
			(uint16_t *)lis_malloc(
				nnz*sizeof(uint16_t),
				"lis_half_matrix_compress::row_value");

		if( value==NULL )
		{
			lis_half_matrix_destroy(H);

			LIS_SETERR_MEM(
				nnz*sizeof(uint16_t));

			return LIS_OUT_OF_MEMORY;
		}

		exponent = 0;

		if( maxabs>0.0 )
		{
			(void)frexp(
				maxabs,
				&exponent);

			/*
			 * The largest row value is placed around
			 * 16384..32768.  Scaling remains an exact
			 * power-of-two operation.
			 */
			exponent -= 15;
		}

		H->scale[i] =
			(int16_t)exponent;

		row_total = 0.0;
		row_lost = 0.0;

		for(j=0;j<nnz;j++)
		{
			v = (double)A->value[i][j];
			a = fabs(v);

			scaled =
				ldexp(
					v,
					-exponent);

			if( !isfinite(scaled) )
			{
				lis_free(value);
				lis_half_matrix_destroy(H);

				LIS_SETERR(
					LIS_ERR_ILL_ARG,
					"half factor scaling overflow\n");

				return LIS_ERR_ILL_ARG;
			}

			value[j] =
				lis_half_from_double(
					scaled);

			/*
			 * binary16 exponent=31 means Inf or NaN.
			 */
			if( (value[j] & 0x7c00u)==0x7c00u )
			{
				lis_free(value);
				lis_half_matrix_destroy(H);

				LIS_SETERR(
					LIS_ERR_ILL_ARG,
					"half factor storage overflow\n");

				return LIS_ERR_ILL_ARG;
			}

			row_total += a;

			if( v!=0.0 &&
				(value[j] & 0x7fffu)==0 )
			{
				row_lost += a;
			}

			decoded =
				ldexp(
					lis_half_to_double(
						value[j]),
					exponent);

			if( !isfinite(decoded) )
			{
				lis_free(value);
				lis_half_matrix_destroy(H);

				LIS_SETERR(
					LIS_ERR_ILL_ARG,
					"half factor restoration overflow\n");

				return LIS_ERR_ILL_ARG;
			}
		}

		H->value[i] = value;

		if( row_total>0.0 )
		{
			row_ratio =
				row_lost /
				row_total;

			if( row_ratio>
				*worst_underflow )
			{
				*worst_underflow =
					(LIS_REAL)row_ratio;
			}
		}
	}

	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_half_matrix_expand"
static LIS_INT
lis_half_matrix_expand(
	const struct LIS_HALF_ILU_MATRIX_STRUCT *H,
	LIS_MATRIX_ILU *A)
{
	LIS_INT i,j,nnz,err;
	LIS_INT *index;
	LIS_SCALAR *value;

	if( H==NULL ||
		A==NULL ||
		H->pattern==NULL )
	{
		return LIS_ERR_ILL_ARG;
	}

	*A = NULL;

	err = lis_matrix_ilu_create(
		H->n,
		1,
		A);

	if( err )
	{
		return err;
	}

	err = lis_matrix_ilu_setCR(*A);

	if( err )
	{
		lis_matrix_ilu_destroy(*A);
		*A = NULL;
		return err;
	}

	for(i=0;i<H->n;i++)
	{
		(*A)->value[i] = NULL;
	}

	for(i=0;i<H->n;i++)
	{
		nnz = H->pattern->nnz[i];

		if( nnz==0 )
		{
			continue;
		}

		index = (LIS_INT *)malloc(
			nnz*sizeof(LIS_INT));

		value = (LIS_SCALAR *)malloc(
			nnz*sizeof(LIS_SCALAR));

		if( index==NULL || value==NULL )
		{
			free(index);
			free(value);

			lis_matrix_ilu_destroy(*A);
			*A = NULL;

			LIS_SETERR_MEM(
				nnz*sizeof(LIS_SCALAR));

			return LIS_OUT_OF_MEMORY;
		}

		for(j=0;j<nnz;j++)
		{
			index[j] = H->pattern->index[i][j];

			value[j] = (LIS_SCALAR)ldexp(
				lis_half_to_double(
					H->value[i][j]),
				(int)H->scale[i]);
		}

		(*A)->index[i] = index;
		(*A)->value[i] = value;
		(*A)->nnz[i] = nnz;
	}

	return LIS_SUCCESS;
}

#endif


#undef __FUNC__
#define __FUNC__ "lis_half_ilu_destroy"
LIS_INT
lis_half_ilu_destroy(
	LIS_HALF_ILU half)
{
	if( half==NULL )
	{
		return LIS_SUCCESS;
	}

	lis_half_matrix_destroy(&half->L);
	lis_half_matrix_destroy(&half->U);

	lis_free(half);

	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_half_ilu_compress"
LIS_INT
lis_half_ilu_compress(
	LIS_MATRIX_ILU L,
	LIS_MATRIX_ILU U,
	LIS_HALF_ILU *half,
	LIS_REAL *underflow_ratio)
{
#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

	LIS_HALF_ILU object;
	LIS_REAL worst;
	LIS_INT err;

	if( half==NULL ||
		L==NULL ||
		U==NULL ||
		L->n!=U->n ||
		L->bs!=1 ||
		U->bs!=1 )
	{
		return LIS_ERR_ILL_ARG;
	}

	*half = NULL;
	worst = 0.0;

	object = (LIS_HALF_ILU)lis_malloc(
		sizeof(struct LIS_HALF_ILU_STRUCT),
		"lis_half_ilu_compress::half");

	if( object==NULL )
	{
		LIS_SETERR_MEM(
			sizeof(struct LIS_HALF_ILU_STRUCT));

		return LIS_OUT_OF_MEMORY;
	}

	memset(
		object,
		0,
		sizeof(struct LIS_HALF_ILU_STRUCT));

	err = lis_half_matrix_compress(
		L,
		&object->L,
		&worst);

	if( err )
	{
		lis_half_ilu_destroy(object);
		return err;
	}

	err = lis_half_matrix_compress(
		U,
		&object->U,
		&worst);

	if( err )
	{
		lis_half_ilu_destroy(object);
		return err;
	}

	object->underflow_ratio = worst;

	if( underflow_ratio!=NULL )
	{
		*underflow_ratio = worst;
	}

	*half = object;

	return LIS_SUCCESS;

#else

	if( half!=NULL )
	{
		*half = NULL;
	}

	if( underflow_ratio!=NULL )
	{
		*underflow_ratio = 0.0;
	}

	return LIS_ERR_NOT_IMPLEMENTED;

#endif
}


#undef __FUNC__
#define __FUNC__ "lis_half_ilu_expand"
LIS_INT
lis_half_ilu_expand(
	LIS_HALF_ILU half,
	LIS_MATRIX_ILU *L,
	LIS_MATRIX_ILU *U)
{
#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

	LIS_INT err;

	if( half==NULL || L==NULL || U==NULL )
	{
		return LIS_ERR_ILL_ARG;
	}

	*L = NULL;
	*U = NULL;

	err = lis_half_matrix_expand(
		&half->L,
		L);

	if( err )
	{
		return err;
	}

	err = lis_half_matrix_expand(
		&half->U,
		U);

	if( err )
	{
		lis_matrix_ilu_destroy(*L);
		*L = NULL;
		return err;
	}

	return LIS_SUCCESS;

#else

	if( L!=NULL )
	{
		*L = NULL;
	}

	if( U!=NULL )
	{
		*U = NULL;
	}

	return LIS_ERR_NOT_IMPLEMENTED;

#endif
}



#define LIS_HALF_UNDERFLOW_LIMIT 1.0e-3
#define LIS_HALF_PROBE_TOL_DEFAULT 5.0e-2


#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

static double
lis_half_probe_tolerance(void)
{
	const char *text;
	char *end;
	double value;

	text = getenv("LIS_HALF_PROBE_TOL");

	if( text==NULL || *text=='\0' )
	{
		return LIS_HALF_PROBE_TOL_DEFAULT;
	}

	errno = 0;
	end = NULL;

	value = strtod(text,&end);

	if( errno!=0 ||
		end==text ||
		*end!='\0' ||
		!isfinite(value) ||
		value<=0.0 ||
		value>1.0 )
	{
		return LIS_HALF_PROBE_TOL_DEFAULT;
	}

	return value;
}


static void
lis_half_fill_probe(
	LIS_INT probe,
	LIS_INT n,
	double *z)
{
	LIS_INT i;
	uint32_t state;

	if( probe==0 )
	{
		state = UINT32_C(0x6d2b79f5);

		for(i=0;i<n;i++)
		{
			state =
				UINT32_C(1664525)*state +
				UINT32_C(1013904223);

			z[i] =
				((double)((state >> 8) & 0xffffu)
				 / 32768.0) - 1.0;
		}

		return;
	}

	if( probe==1 )
	{
		for(i=0;i<n;i++)
		{
			z[i] = (i & 1) ? -1.0 : 1.0;
		}

		return;
	}

	for(i=0;i<n;i++)
	{
		z[i] = 0.0;
	}

	if( n>0 )
	{
		z[0] = 1.0;
	}

	if( n>2 )
	{
		z[n/2] = -0.5;
	}

	if( n>1 )
	{
		z[n-1] = 0.25;
	}
}


static int
lis_half_apply_reference(
	LIS_MATRIX_ILU L,
	LIS_MATRIX_ILU U,
	LIS_VECTOR D,
	const double *z,
	double *x)
{
	LIS_INT i,j,jj,row,n;

	n = L->n;

	for(i=0;i<n;i++)
	{
		x[i] = z[i];
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			jj = L->index[i][j];

			if( jj<0 || jj>=n )
			{
				return -1;
			}

			x[i] -=
				(double)L->value[i][j] *
				x[jj];
		}

		if( !isfinite(x[i]) )
		{
			return 0;
		}
	}

	for(i=n;i>0;i--)
	{
		row = i - 1;

		for(j=0;j<U->nnz[row];j++)
		{
			jj = U->index[row][j];

			if( jj<0 || jj>=n )
			{
				return -1;
			}

			x[row] -=
				(double)U->value[row][j] *
				x[jj];
		}

		x[row] =
			(double)D->value[row] *
			x[row];

		if( !isfinite(x[row]) )
		{
			return 0;
		}
	}

	return 1;
}


static int
lis_half_apply_compressed(
	LIS_HALF_ILU half,
	LIS_VECTOR D,
	const double *z,
	double *x)
{
	LIS_INT i,j,jj,row,n;
	double value;

	n = half->L.n;

	for(i=0;i<n;i++)
	{
		x[i] = z[i];
	}

	for(i=0;i<n;i++)
	{
		for(j=0;j<half->L.pattern->nnz[i];j++)
		{
			jj = half->L.pattern->index[i][j];

			if( jj<0 || jj>=n )
			{
				return -1;
			}

			value = ldexp(
				lis_half_to_double(
					half->L.value[i][j]),
				(int)half->L.scale[i]);

			x[i] -= value*x[jj];
		}

		if( !isfinite(x[i]) )
		{
			return 0;
		}
	}

	for(i=n;i>0;i--)
	{
		row = i - 1;

		for(j=0;j<half->U.pattern->nnz[row];j++)
		{
			jj = half->U.pattern->index[row][j];

			if( jj<0 || jj>=n )
			{
				return -1;
			}

			value = ldexp(
				lis_half_to_double(
					half->U.value[row][j]),
				(int)half->U.scale[row]);

			x[row] -= value*x[jj];
		}

		x[row] =
			(double)D->value[row] *
			x[row];

		if( !isfinite(x[row]) )
		{
			return 0;
		}
	}

	return 1;
}


static double
lis_half_relative_probe_error(
	const double *reference,
	const double *candidate,
	LIS_INT n)
{
	LIS_INT i;
	double ref_norm;
	double diff_norm;
	double diff;

	ref_norm = 0.0;
	diff_norm = 0.0;

	for(i=0;i<n;i++)
	{
		if( !isfinite(reference[i]) ||
			!isfinite(candidate[i]) )
		{
			return HUGE_VAL;
		}

		diff = candidate[i] - reference[i];

		if( !isfinite(diff) )
		{
			return HUGE_VAL;
		}

		ref_norm = hypot(
			ref_norm,
			reference[i]);

		diff_norm = hypot(
			diff_norm,
			diff);
	}

	if( !isfinite(ref_norm) ||
		!isfinite(diff_norm) )
	{
		return HUGE_VAL;
	}

	if( ref_norm==0.0 )
	{
		return diff_norm==0.0 ?
			0.0 : HUGE_VAL;
	}

	return diff_norm/ref_norm;
}

#endif


#undef __FUNC__
#define __FUNC__ "lis_half_ilu_assess"
LIS_INT
lis_half_ilu_assess(
	LIS_MATRIX_ILU L,
	LIS_MATRIX_ILU U,
	LIS_VECTOR D,
	LIS_HALF_ILU half,
	LIS_INT mode,
	LIS_INT *accepted,
	LIS_REAL *probe_error,
	LIS_REAL *probe_tolerance)
{
	if( accepted==NULL )
	{
		return LIS_ERR_ILL_ARG;
	}

	*accepted = LIS_FALSE;

	if( probe_error!=NULL )
	{
		*probe_error = 0.0;
	}

	if( probe_tolerance!=NULL )
	{
		*probe_tolerance =
			LIS_HALF_PROBE_TOL_DEFAULT;
	}

	if( mode<LIS_PRECON_HALF_OFF ||
		mode>LIS_PRECON_HALF_NOPROBE )
	{
		return LIS_ERR_ILL_ARG;
	}

	if( mode==LIS_PRECON_HALF_OFF )
	{
		return LIS_SUCCESS;
	}

#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

	{
		LIS_INT n;
		LIS_INT probe;
		LIS_INT result;
		double *z;
		double *reference;
		double *candidate;
		double tolerance;
		double error;
		double worst;

		if( L==NULL ||
			U==NULL ||
			D==NULL ||
			half==NULL ||
			L->n!=U->n ||
			L->n!=half->L.n ||
			L->n!=half->U.n ||
			D->n<L->n )
		{
			return LIS_ERR_ILL_ARG;
		}

		/*
		 * Safety checks are never bypassed, including NOPROBE.
		 */
		if( !isfinite(
				(double)half->underflow_ratio) ||
			half->underflow_ratio >
				(LIS_REAL)LIS_HALF_UNDERFLOW_LIMIT )
		{
			return LIS_SUCCESS;
		}

		tolerance =
			lis_half_probe_tolerance();

		if( probe_tolerance!=NULL )
		{
			*probe_tolerance =
				(LIS_REAL)tolerance;
		}

		if( mode==
			LIS_PRECON_HALF_NOPROBE )
		{
			*accepted = LIS_TRUE;
			return LIS_SUCCESS;
		}

		n = L->n;

		if( n==0 )
		{
			*accepted = LIS_TRUE;
			return LIS_SUCCESS;
		}

		z = (double *)lis_malloc(
			n*sizeof(double),
			"lis_half_ilu_assess::z");

		reference = (double *)lis_malloc(
			n*sizeof(double),
			"lis_half_ilu_assess::reference");

		candidate = (double *)lis_malloc(
			n*sizeof(double),
			"lis_half_ilu_assess::candidate");

		if( z==NULL ||
			reference==NULL ||
			candidate==NULL )
		{
			if( z!=NULL )
			{
				lis_free(z);
			}

			if( reference!=NULL )
			{
				lis_free(reference);
			}

			if( candidate!=NULL )
			{
				lis_free(candidate);
			}

			return LIS_OUT_OF_MEMORY;
		}

		worst = 0.0;

		for(probe=0;probe<3;probe++)
		{
			lis_half_fill_probe(
				probe,
				n,
				z);

			result =
				lis_half_apply_reference(
					L,U,D,z,reference);

			if( result<0 )
			{
				lis_free(z);
				lis_free(reference);
				lis_free(candidate);

				return LIS_ERR_ILL_ARG;
			}

			if( result==0 )
			{
				worst = HUGE_VAL;
				break;
			}

			result =
				lis_half_apply_compressed(
					half,D,z,candidate);

			if( result<0 )
			{
				lis_free(z);
				lis_free(reference);
				lis_free(candidate);

				return LIS_ERR_ILL_ARG;
			}

			if( result==0 )
			{
				worst = HUGE_VAL;
				break;
			}

			error =
				lis_half_relative_probe_error(
					reference,
					candidate,
					n);

			if( error>worst )
			{
				worst = error;
			}
		}

		lis_free(z);
				lis_free(reference);
				lis_free(candidate);

		if( probe_error!=NULL )
		{
			*probe_error =
				(LIS_REAL)worst;
		}

		if( isfinite(worst) &&
			worst<=tolerance )
		{
			*accepted = LIS_TRUE;
		}

		return LIS_SUCCESS;
	}

#else

	/*
	 * Unsupported scalar modes fall back to normal factor storage.
	 */
	return LIS_SUCCESS;

#endif
}



#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

#if LIS_HALF_X86_F16C_DISPATCH

__attribute__((target("avx,f16c")))
static double
lis_half_to_double_f16c(
	uint16_t value)
{
	/*
	 * Every finite binary16 value is exactly representable
	 * in float, so half -> float -> double is exact.
	 */
	return (double)_cvtsh_ss(value);
}

#endif


static LIS_INT
lis_half_runtime_use_f16c(void)
{
	const char *disable;

	disable =
		getenv("LIS_HALF_NO_SIMD");

	if( disable!=NULL &&
		strcmp(disable,"1")==0 )
	{
		return LIS_FALSE;
	}

#if LIS_HALF_X86_F16C_DISPATCH

	__builtin_cpu_init();

	if( __builtin_cpu_supports("avx") &&
		__builtin_cpu_supports("f16c") )
	{
		return LIS_TRUE;
	}

#endif

	return LIS_FALSE;
}


static double
lis_half_ilu_scale_factor(
	const struct LIS_HALF_ILU_MATRIX_STRUCT *A,
	LIS_INT row,
	LIS_INT *direct_scale)
{
	int exponent;
	uint64_t bits;

	exponent =
		(int)A->scale[row];

	/*
	 * Build an exact binary64 power-of-two scale for the common
	 * apply path without a per-entry ldexp call.
	 */
	if( exponent>=-1022 &&
		exponent<=1023 )
	{
		bits =
			(uint64_t)
			(exponent + 1023)
			<< 52;

		*direct_scale = LIS_TRUE;

		return
			lis_half_double_from_bits(
				bits);
	}

	if( exponent>=-1074 &&
		exponent<-1022 )
	{
		bits =
			((uint64_t)1)
			<< (exponent + 1074);

		*direct_scale = LIS_TRUE;

		return
			lis_half_double_from_bits(
				bits);
	}

	/*
	 * Extremely small row exponents retain the original ldexp
	 * path so the existing numerical range contract is unchanged.
	 */
	*direct_scale = LIS_FALSE;

	return 0.0;
}


static LIS_SCALAR
lis_half_ilu_value(
	const struct LIS_HALF_ILU_MATRIX_STRUCT *A,
	LIS_INT row,
	LIS_INT entry,
	LIS_INT use_f16c,
	double scale_factor,
	LIS_INT direct_scale)
{
	double decoded;

#if LIS_HALF_X86_F16C_DISPATCH

	if( use_f16c )
	{
		decoded =
			lis_half_to_double_f16c(
				A->value[row][entry]);
	}
	else

#endif

	{
		decoded =
			lis_half_to_double(
				A->value[row][entry]);
	}

	if( direct_scale )
	{
		return
			(LIS_SCALAR)
			(decoded * scale_factor);
	}

	return
		(LIS_SCALAR)
		ldexp(
			decoded,
			(int)A->scale[row]);
}


#if LIS_HALF_X86_F16C_DISPATCH

__attribute__((target("avx,f16c")))
static void
lis_half_ilu_row_subtract_f16c(
	const struct LIS_HALF_ILU_MATRIX_STRUCT *A,
	LIS_INT row,
	LIS_SCALAR *x,
	double scale_factor)
{
	LIS_INT j,jj;
	double decoded;

	for(j=0;j<A->pattern->nnz[row];j++)
	{
		jj =
			A->pattern->index[row][j];

		decoded =
			(double)_cvtsh_ss(
				A->value[row][j]);

		x[row] -=
			(LIS_SCALAR)
			(decoded * scale_factor) *
			x[jj];
	}
}


__attribute__((target("avx,f16c")))
static void
lis_half_ilu_row_scatter_f16c(
	const struct LIS_HALF_ILU_MATRIX_STRUCT *A,
	LIS_INT row,
	LIS_SCALAR *x,
	double scale_factor)
{
	LIS_INT j,jj;
	double decoded;

	for(j=0;j<A->pattern->nnz[row];j++)
	{
		jj =
			A->pattern->index[row][j];

		decoded =
			(double)_cvtsh_ss(
				A->value[row][j]);

		x[jj] -=
			(LIS_SCALAR)
			(decoded * scale_factor) *
			x[row];
	}
}

#endif

#endif


#undef __FUNC__
#define __FUNC__ "lis_half_ilu_psolve"
LIS_INT
lis_half_ilu_psolve(
	LIS_SOLVER solver,
	LIS_VECTOR B,
	LIS_VECTOR X)
{
#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

	LIS_PRECON precon;
	LIS_HALF_ILU half;
	LIS_VECTOR D;
	LIS_SCALAR *x;
	LIS_INT n;
	LIS_INT i,j,jj;
	LIS_INT use_f16c;
	LIS_INT direct_scale;
	double scale_factor;

#ifdef _OPENMP
	LIS_INT is,ie,my_rank,nprocs;
#endif

	if( solver==NULL ||
		solver->precon==NULL ||
		solver->A==NULL ||
		B==NULL ||
		X==NULL )
	{
		return LIS_ERR_ILL_ARG;
	}

	precon = solver->precon;
	half = precon->half;
	D = precon->D;

	if( half==NULL ||
		D==NULL ||
		half->L.n!=solver->A->n ||
		half->U.n!=solver->A->n )
	{
		return LIS_ERR_ILL_ARG;
	}

	n = solver->A->n;

	use_f16c = lis_half_runtime_use_f16c();

	lis_vector_copy(B,X);
	x = X->value;

#ifdef _OPENMP

	nprocs = omp_get_max_threads();

	#pragma omp parallel private(i,j,jj,is,ie,my_rank,direct_scale,scale_factor)
	{
		my_rank = omp_get_thread_num();

		LIS_GET_ISIE(
			my_rank,
			nprocs,
			n,
			is,
			ie);

		for(i=is;i<ie;i++)
		{
			scale_factor =
				lis_half_ilu_scale_factor(
					&half->L,
					i,
					&direct_scale);

			#if LIS_HALF_X86_F16C_DISPATCH
			if( use_f16c && direct_scale )
			{
				lis_half_ilu_row_subtract_f16c(
					&half->L,
					i,
					x,
					scale_factor);
			}
			else
			#endif
			{
				for(j=0;j<half->L.pattern->nnz[i];j++)
				{
					jj = half->L.pattern->index[i][j];

					x[i] -=
						lis_half_ilu_value(
							&half->L,
							i,
							j,
							use_f16c,
							scale_factor,
							direct_scale) *
						x[jj];
				}
			}
		}

		for(i=ie-1;i>=is;i--)
		{
			scale_factor =
				lis_half_ilu_scale_factor(
					&half->U,
					i,
					&direct_scale);

			#if LIS_HALF_X86_F16C_DISPATCH
			if( use_f16c && direct_scale )
			{
				lis_half_ilu_row_subtract_f16c(
					&half->U,
					i,
					x,
					scale_factor);
			}
			else
			#endif
			{
				for(j=0;j<half->U.pattern->nnz[i];j++)
				{
					jj = half->U.pattern->index[i][j];

					x[i] -=
						lis_half_ilu_value(
							&half->U,
							i,
							j,
							use_f16c,
							scale_factor,
							direct_scale) *
						x[jj];
				}
			}

			x[i] =
				D->value[i] *
				x[i];
		}
	}

#else

	for(i=0;i<n;i++)
	{
		scale_factor =
			lis_half_ilu_scale_factor(
				&half->L,
				i,
				&direct_scale);

		#if LIS_HALF_X86_F16C_DISPATCH
		if( use_f16c && direct_scale )
		{
			lis_half_ilu_row_subtract_f16c(
				&half->L,
				i,
				x,
				scale_factor);
		}
		else
		#endif
		{
			for(j=0;j<half->L.pattern->nnz[i];j++)
			{
				jj = half->L.pattern->index[i][j];

				x[i] -=
					lis_half_ilu_value(
						&half->L,
						i,
						j,
						use_f16c,
						scale_factor,
						direct_scale) *
					x[jj];
			}
		}
	}

	for(i=n-1;i>=0;i--)
	{
		scale_factor =
			lis_half_ilu_scale_factor(
				&half->U,
				i,
				&direct_scale);

		#if LIS_HALF_X86_F16C_DISPATCH
		if( use_f16c && direct_scale )
		{
			lis_half_ilu_row_subtract_f16c(
				&half->U,
				i,
				x,
				scale_factor);
		}
		else
		#endif
		{
			for(j=0;j<half->U.pattern->nnz[i];j++)
			{
				jj = half->U.pattern->index[i][j];

				x[i] -=
					lis_half_ilu_value(
						&half->U,
						i,
						j,
						use_f16c,
						scale_factor,
						direct_scale) *
					x[jj];
			}
		}

		x[i] =
			D->value[i] *
			x[i];
	}

#endif

	return LIS_SUCCESS;

#else

	(void)solver;
	(void)B;
	(void)X;

	return LIS_ERR_NOT_IMPLEMENTED;

#endif
}


#undef __FUNC__
#define __FUNC__ "lis_half_ilu_psolveh"
LIS_INT
lis_half_ilu_psolveh(
	LIS_SOLVER solver,
	LIS_VECTOR B,
	LIS_VECTOR X)
{
#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

	LIS_PRECON precon;
	LIS_HALF_ILU half;
	LIS_VECTOR D;
	LIS_SCALAR *x;
	LIS_INT n;
	LIS_INT i,j,jj;
	LIS_INT use_f16c;
	LIS_INT direct_scale;
	double scale_factor;

#ifdef _OPENMP
	LIS_INT is,ie,my_rank,nprocs;
#endif

	if( solver==NULL ||
		solver->precon==NULL ||
		solver->A==NULL ||
		B==NULL ||
		X==NULL )
	{
		return LIS_ERR_ILL_ARG;
	}

	precon = solver->precon;
	half = precon->half;
	D = precon->D;

	if( half==NULL ||
		D==NULL ||
		half->L.n!=solver->A->n ||
		half->U.n!=solver->A->n )
	{
		return LIS_ERR_ILL_ARG;
	}

	n = solver->A->n;

	use_f16c = lis_half_runtime_use_f16c();

	lis_vector_copy(B,X);
	x = X->value;

#ifdef _OPENMP

	nprocs = omp_get_max_threads();

	#pragma omp parallel private(i,j,jj,is,ie,my_rank,direct_scale,scale_factor)
	{
		my_rank = omp_get_thread_num();

		LIS_GET_ISIE(
			my_rank,
			nprocs,
			n,
			is,
			ie);

		for(i=is;i<ie;i++)
		{
			x[i] =
				D->value[i] *
				x[i];

			scale_factor =
				lis_half_ilu_scale_factor(
					&half->U,
					i,
					&direct_scale);

			#if LIS_HALF_X86_F16C_DISPATCH
			if( use_f16c && direct_scale )
			{
				lis_half_ilu_row_scatter_f16c(
					&half->U,
					i,
					x,
					scale_factor);
			}
			else
			#endif
			{
				for(j=0;j<half->U.pattern->nnz[i];j++)
				{
					jj = half->U.pattern->index[i][j];

					x[jj] -=
						lis_half_ilu_value(
							&half->U,
							i,
							j,
							use_f16c,
							scale_factor,
							direct_scale) *
						x[i];
				}
			}
		}

		for(i=ie-1;i>=is;i--)
		{
			scale_factor =
				lis_half_ilu_scale_factor(
					&half->L,
					i,
					&direct_scale);

			#if LIS_HALF_X86_F16C_DISPATCH
			if( use_f16c && direct_scale )
			{
				lis_half_ilu_row_scatter_f16c(
					&half->L,
					i,
					x,
					scale_factor);
			}
			else
			#endif
			{
				for(j=0;j<half->L.pattern->nnz[i];j++)
				{
					jj = half->L.pattern->index[i][j];

					x[jj] -=
						lis_half_ilu_value(
							&half->L,
							i,
							j,
							use_f16c,
							scale_factor,
							direct_scale) *
						x[i];
				}
			}
		}
	}

#else

	for(i=0;i<n;i++)
	{
		x[i] =
			D->value[i] *
			x[i];

		scale_factor =
			lis_half_ilu_scale_factor(
				&half->U,
				i,
				&direct_scale);

		#if LIS_HALF_X86_F16C_DISPATCH
		if( use_f16c && direct_scale )
		{
			lis_half_ilu_row_scatter_f16c(
				&half->U,
				i,
				x,
				scale_factor);
		}
		else
		#endif
		{
			for(j=0;j<half->U.pattern->nnz[i];j++)
			{
				jj = half->U.pattern->index[i][j];

				x[jj] -=
					lis_half_ilu_value(
						&half->U,
						i,
						j,
						use_f16c,
						scale_factor,
						direct_scale) *
					x[i];
			}
		}
	}

	for(i=n-1;i>=0;i--)
	{
		scale_factor =
			lis_half_ilu_scale_factor(
				&half->L,
				i,
				&direct_scale);

		#if LIS_HALF_X86_F16C_DISPATCH
		if( use_f16c && direct_scale )
		{
			lis_half_ilu_row_scatter_f16c(
				&half->L,
				i,
				x,
				scale_factor);
		}
		else
		#endif
		{
			for(j=0;j<half->L.pattern->nnz[i];j++)
			{
				jj = half->L.pattern->index[i][j];

				x[jj] -=
					lis_half_ilu_value(
						&half->L,
						i,
						j,
						use_f16c,
						scale_factor,
						direct_scale) *
					x[i];
			}
		}
	}

#endif

	return LIS_SUCCESS;

#else

	(void)solver;
	(void)B;
	(void)X;

	return LIS_ERR_NOT_IMPLEMENTED;

#endif
}


#undef __FUNC__
#define __FUNC__ "lis_half_get_capability"
LIS_INT lis_half_get_capability(
	LIS_INT *available)
{
	LIS_DEBUG_FUNC_IN;

	if( available==NULL )
	{
		LIS_SETERR(
			LIS_ERR_ILL_ARG,
			"half precision capability output pointer is NULL\n");

		return LIS_ERR_ILL_ARG;
	}

#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

	*available = LIS_TRUE;

#else

	*available = LIS_FALSE;

#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
