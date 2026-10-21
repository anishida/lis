#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <stdio.h>

#include "lis.h"


static int
check_mode(
	LIS_SOLVER solver,
	LIS_INT expected,
	const char *label)
{
	LIS_INT err;
	LIS_INT mode;

	mode = -1;

	err =
		lis_solver_get_precon_half(
			solver,
			&mode);

	if( err!=LIS_SUCCESS ||
		mode!=expected )
	{
		fprintf(
			stderr,
			"%s: rc=%d mode=%d expected=%d\n",
			label,
			(int)err,
			(int)mode,
			(int)expected);

		return 1;
	}

	return 0;
}


int
main(
	int argc,
	char **argv)
{
	LIS_SOLVER solver;

	LIS_INT err;
	LIS_INT capability;
	LIS_INT expected_capability;

	int failed;

	solver = NULL;
	capability = -1;
	failed = 0;

	err =
		lis_initialize(
			&argc,
			&argv);

	if( err )
	{
		return 1;
	}


#if defined(USE_HALF_PRECISION) && \
	!defined(_LONG__DOUBLE) && \
	!defined(USE_COMPLEX) && \
	!defined(_COMPLEX)

	expected_capability = LIS_TRUE;

#else

	expected_capability = LIS_FALSE;

#endif


	err =
		lis_half_get_capability(
			&capability);

	if( err!=LIS_SUCCESS ||
		capability!=expected_capability )
	{
		fprintf(
			stderr,
			"CAPABILITY: rc=%d got=%d expected=%d\n",
			(int)err,
			(int)capability,
			(int)expected_capability);

		failed = 1;
		goto cleanup;
	}

	if( lis_half_get_capability(NULL)!=
		LIS_ERR_ILL_ARG )
	{
		fprintf(
			stderr,
			"NULL capability query was not rejected\n");

		failed = 1;
		goto cleanup;
	}


	err =
		lis_solver_create(
			&solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}


	if( check_mode(
			solver,
			LIS_PRECON_HALF_OFF,
			"DEFAULT") )
	{
		failed = 1;
		goto cleanup;
	}


	err =
		lis_solver_set_precon_half(
			solver,
			-1);

	if( err!=LIS_ERR_ILL_ARG )
	{
		fprintf(
			stderr,
			"NEGATIVE mode rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}


	err =
		lis_solver_set_precon_half(
			solver,
			3);

	if( err!=LIS_ERR_ILL_ARG )
	{
		fprintf(
			stderr,
			"OUT-OF-RANGE mode rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}


	err =
		lis_solver_set_option(
			"-precon_half invalid",
			solver);

	if( err!=LIS_ERR_ILL_ARG )
	{
		fprintf(
			stderr,
			"MALFORMED parser rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}


	err =
		lis_solver_set_precon_half(
			solver,
			LIS_PRECON_HALF_ON);

	if( err )
	{
		fprintf(
			stderr,
			"DIRECT ON failed rc=%d capability=%d\n",
			(int)err,
			(int)capability);

		failed = 1;
		goto cleanup;
	}

	if( check_mode(
			solver,
			LIS_PRECON_HALF_ON,
			"DIRECT-ON") )
	{
		failed = 1;
		goto cleanup;
	}


	err =
		lis_solver_set_option(
			"-precon_half 2",
			solver);

	if( err )
	{
		fprintf(
			stderr,
			"PARSER NOPROBE failed rc=%d capability=%d\n",
			(int)err,
			(int)capability);

		failed = 1;
		goto cleanup;
	}

	if( check_mode(
			solver,
			LIS_PRECON_HALF_NOPROBE,
			"PARSER-NOPROBE") )
	{
		failed = 1;
		goto cleanup;
	}


	err =
		lis_solver_set_option(
			"-precon_half 0",
			solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	if( check_mode(
			solver,
			LIS_PRECON_HALF_OFF,
			"PARSER-OFF") )
	{
		failed = 1;
		goto cleanup;
	}


	printf(
		"HALF_CONFIG PASSED "
		"capability=%d policy_independent=1\n",
		(int)capability);


cleanup:

	if( solver!=NULL )
	{
		lis_solver_destroy(
			solver);
	}

	lis_finalize();

	return failed ? 1 : 0;
}
