#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "lis.h"

/* Internal C1 conversion interface. */
extern uint16_t lis_half_from_double(double value);
extern double lis_half_to_double(uint16_t value);


static int
expect_encode(
	const char *name,
	double value,
	uint16_t expected)
{
	uint16_t got;

	got = lis_half_from_double(value);

	if( got!=expected )
	{
		fprintf(
			stderr,
			"%s: got=0x%04x expected=0x%04x\n",
			name,
			(unsigned int)got,
			(unsigned int)expected);

		return 1;
	}

	return 0;
}


int
main(void)
{
	unsigned int i;
	unsigned int exp;
	unsigned int frac;
	uint16_t h;
	uint16_t back;
	double d;
	double smallest;
	int failed;

	failed = 0;

	failed |= expect_encode("+zero", 0.0, 0x0000u);
	failed |= expect_encode("-zero", -0.0, 0x8000u);
	failed |= expect_encode("+one", 1.0, 0x3c00u);
	failed |= expect_encode("-two", -2.0, 0xc000u);

	failed |= expect_encode(
		"min-normal",
		ldexp(1.0,-14),
		0x0400u);

	failed |= expect_encode(
		"max-subnormal",
		ldexp(1023.0,-24),
		0x03ffu);

	failed |= expect_encode(
		"min-subnormal",
		ldexp(1.0,-24),
		0x0001u);

	failed |= expect_encode(
		"zero-tie-even",
		ldexp(1.0,-25),
		0x0000u);

	failed |= expect_encode(
		"subnormal-tie-1.5",
		ldexp(3.0,-25),
		0x0002u);

	failed |= expect_encode(
		"subnormal-tie-2.5",
		ldexp(5.0,-25),
		0x0002u);

	failed |= expect_encode(
		"normal-tie-even-low",
		1.0 + ldexp(1.0,-11),
		0x3c00u);

	failed |= expect_encode(
		"normal-tie-even-high",
		1.0 + 3.0*ldexp(1.0,-11),
		0x3c02u);

	failed |= expect_encode(
		"max-finite",
		65504.0,
		0x7bffu);

	failed |= expect_encode(
		"overflow-tie",
		65520.0,
		0x7c00u);

	failed |= expect_encode(
		"below-overflow-tie",
		nextafter(65520.0,0.0),
		0x7bffu);

	failed |= expect_encode(
		"+inf",
		HUGE_VAL,
		0x7c00u);

	failed |= expect_encode(
		"-inf",
		-HUGE_VAL,
		0xfc00u);

	smallest = nextafter(0.0,1.0);

	failed |= expect_encode(
		"binary64-min-subnormal",
		smallest,
		0x0000u);

	failed |= expect_encode(
		"binary64-max",
		DBL_MAX,
		0x7c00u);

	if( !isnan(
			lis_half_to_double(
				lis_half_from_double(NAN))) )
	{
		fprintf(stderr,"NaN conversion failed\n");
		failed = 1;
	}

	/*
	 * Exhaust every binary16 bit pattern.
	 *
	 * Finite values and infinities must round-trip exactly.
	 * NaNs are allowed to be quieted but must remain NaNs.
	 */
	for(i=0;i<65536u;i++)
	{
		h = (uint16_t)i;
		exp = (i >> 10) & 0x1fu;
		frac = i & 0x03ffu;

		d = lis_half_to_double(h);

		if( exp==0x1fu && frac!=0 )
		{
			if( !isnan(d) )
			{
				fprintf(
					stderr,
					"NaN decode failed: 0x%04x\n",
					i);

				failed = 1;
				break;
			}

			continue;
		}

		back = lis_half_from_double(d);

		if( back!=h )
		{
			fprintf(
				stderr,
				"roundtrip failed: "
				"0x%04x -> %.17g -> 0x%04x\n",
				i,
				d,
				(unsigned int)back);

			failed = 1;
			break;
		}
	}

	if( signbit(lis_half_to_double(0x8000u))==0 )
	{
		fprintf(stderr,"negative zero sign lost\n");
		failed = 1;
	}

	if( failed )
	{
		return 1;
	}

	printf(
		"HALF_BINARY16_C1 PASSED patterns=65536\n");

	return 0;
}
