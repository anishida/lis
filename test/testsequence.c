#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"
#include "lis_solver.h"

#include <stdio.h>

static int create_count = 0;
static int destroy_count = 0;

static LIS_INT
sequence_create(LIS_SOLVER solver, LIS_PRECON precon)
{
	(void)solver;
	(void)precon;

	create_count++;
	return LIS_SUCCESS;
}

static LIS_INT
sequence_psolve(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	(void)solver;
	return lis_vector_copy(b,x);
}

static LIS_INT
sequence_destroy(LIS_PRECON precon)
{
	(void)precon;

	destroy_count++;
	return LIS_SUCCESS;
}

static LIS_INT
create_matrix(LIS_SCALAR diagonal, LIS_MATRIX *A)
{
	LIS_INT err;
	LIS_INT i,is,ie;

	err = lis_matrix_create(LIS_COMM_WORLD,A);
	if( err ) return err;

	err = lis_matrix_set_size(*A,0,2);
	if( err ) return err;

	err = lis_matrix_get_range(*A,&is,&ie);
	if( err ) return err;

	for(i=is;i<ie;i++)
	{
		err = lis_matrix_set_value(
			LIS_INS_VALUE,
			i,i,
			diagonal,
			*A);
		if( err ) return err;
	}

	err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
	if( err ) return err;

	return lis_matrix_assemble(*A);
}

static int
check_counts(
	int creates,
	int destroys,
	const char *label)
{
	if( create_count!=creates ||
	    destroy_count!=destroys )
	{
		fprintf(
			stderr,
			"%s: create/destroy=%d/%d expected=%d/%d\n",
			label,
			create_count,
			destroy_count,
			creates,
			destroys);
		return 1;
	}

	return 0;
}

int main(int argc, char **argv)
{
	LIS_MATRIX A0 = NULL;
	LIS_MATRIX A1 = NULL;
	LIS_SOLVER solver = NULL;
	LIS_SOLVER legacy_solver = NULL;
	LIS_SOLVER sequence_solver = NULL;
	LIS_PRECON sequence_precon0 = NULL;
	void *sequence_work0 = NULL;
	LIS_VECTOR b = NULL;
	LIS_VECTOR x = NULL;
	LIS_PRECON p0 = NULL;
	LIS_PRECON p1 = NULL;
	LIS_INT err;
	int failed = 0;

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

	err = lis_precon_register_ex(
		"seqs4",
		sequence_create,
		sequence_psolve,
		sequence_psolve,
		sequence_destroy);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = create_matrix((LIS_SCALAR)2.0,&A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = create_matrix((LIS_SCALAR)3.0,&A1);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_create(&solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_set_option(
		"-i gmres -p seqs4 "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	/* BUILD */
	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p0);
	if( err || p0==NULL )
	{
		fprintf(stderr,"BUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(1,0,"BUILD");
	if( failed ) goto cleanup;

	/* REUSE */
	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1!=p0 )
	{
		fprintf(stderr,"REUSE failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(1,0,"REUSE");
	if( failed ) goto cleanup;

	/*
	 * Configuration change must invalidate the persistent
	 * preconditioner even when matrix generations are unchanged.
	 */
	err = lis_solver_set_option(
		"-ilu_fill 1",
		solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"CONFIG_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(2,1,"CONFIG_REBUILD");
	if( failed ) goto cleanup;

	/* Numerical generation change -> REBUILD. */
	err = lis_matrix_mark_values_changed(A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"VALUE_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(3,2,"VALUE_REBUILD");
	if( failed ) goto cleanup;

	/* Structural generation change -> REBUILD. */
	err = lis_matrix_mark_structure_changed(A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"STRUCTURE_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(
		4,3,"STRUCTURE_REBUILD");
	if( failed ) goto cleanup;

	/*
	 * Different matrix object with compatible generations
	 * must also force a rebuild.
	 */
	err = lis_solver_sequence_prepare_precon(
		A1,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"MATRIX_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(5,4,"MATRIX_REBUILD");
	if( failed ) goto cleanup;

	/* Explicit reset. */
	err = lis_solver_sequence_precon_destroy(solver);
	if( err )
	{
		fprintf(stderr,"RESET failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(5,5,"RESET");
	if( failed ) goto cleanup;

	if( solver->sequence_precon!=NULL ||
	    solver->precon!=NULL )
	{
		fprintf(
			stderr,
			"RESET left live preconditioner state\n");
		failed = 1;
		goto cleanup;
	}

	/*
	 * Final BUILD is intentionally left alive.
	 * lis_solver_destroy() must release it.
	 */
	err = lis_solver_sequence_prepare_precon(
		A1,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"FINAL_BUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(6,5,"FINAL_BUILD");
	if( failed ) goto cleanup;

	lis_solver_destroy(solver);
	solver = NULL;

	failed |= check_counts(
		6,6,"SOLVER_DESTROY");

	/*
	 * Legacy lis_solve() still owns a temporary preconditioner.
	 * It must destroy it and leave solver->precon clear.
	 */
	err = lis_solver_create(&legacy_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_set_option(
		"-i gmres -p seqs4 "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		legacy_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A0,&b);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A0,&x);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_set_all((LIS_SCALAR)1.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve(A0,b,x,legacy_solver);
	if( err )
	{
		fprintf(stderr,"LEGACY_SOLVE failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(7,7,"LEGACY_SOLVE");
	if( failed ) goto cleanup;

	if( legacy_solver->precon!=NULL )
	{
		fprintf(stderr,
			"legacy lis_solve left dangling preconditioner\n");
		failed = 1;
		goto cleanup;
	}

	printf("LEGACY_SOLVE PASSED\n");

	/* ---------------------------------------------------- */
	/* S5 public sequence-solve API                         */
	/* ---------------------------------------------------- */

	err = lis_solver_create(&sequence_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_set_option(
		"-i gmres -p seqs4 "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		sequence_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	/* First solve -> BUILD. */
	lis_vector_set_all((LIS_SCALAR)1.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(
		A0,b,x,sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_0 failed\n");
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->retcode!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_0 iterative status=%d\n",
			(int)sequence_solver->retcode);
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(8,7,"SEQUENCE_SOLVE_0");
	if( failed ) goto cleanup;

	if( sequence_solver->sequence_precon==NULL ||
	    sequence_solver->work==NULL )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_0 did not retain reusable state\n");
		failed = 1;
		goto cleanup;
	}

	sequence_precon0 = sequence_solver->sequence_precon;
	sequence_work0 = (void *)sequence_solver->work;

	/* Second solve -> REUSE preconditioner and workspace. */
	lis_vector_set_all((LIS_SCALAR)2.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(
		A0,b,x,sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_1 failed\n");
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->retcode!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_1 iterative status=%d\n",
			(int)sequence_solver->retcode);
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->sequence_precon!=sequence_precon0 )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_1 rebuilt preconditioner\n");
		failed = 1;
		goto cleanup;
	}

	if( (void *)sequence_solver->work!=sequence_work0 )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_1 rebuilt workspace\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(8,7,"SEQUENCE_SOLVE_1");
	if( failed ) goto cleanup;

	/* Matrix value generation change -> REBUILD. */
	err = lis_matrix_mark_values_changed(A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_set_all((LIS_SCALAR)3.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(
		A0,b,x,sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->retcode!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_REBUILD iterative status=%d\n",
			(int)sequence_solver->retcode);
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(
		9,8,"SEQUENCE_SOLVE_REBUILD");
	if( failed ) goto cleanup;

	if( sequence_solver->sequence_precon==NULL )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_REBUILD lost preconditioner\n");
		failed = 1;
		goto cleanup;
	}

	/* Explicit sequence reset. */
	err = lis_solver_reset_sequence(sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_RESET failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(9,9,"SEQUENCE_RESET");
	if( failed ) goto cleanup;

	if( sequence_solver->sequence_precon!=NULL ||
	    sequence_solver->precon!=NULL ||
	    sequence_solver->work!=NULL ||
	    sequence_solver->rhistory!=NULL ||
	    sequence_solver->A!=NULL ||
	    sequence_solver->b!=NULL ||
	    sequence_solver->x!=NULL ||
	    sequence_solver->xx!=NULL )
	{
		fprintf(
			stderr,
			"SEQUENCE_RESET left persistent/transient state\n");
		failed = 1;
		goto cleanup;
	}

	/* Reset must be idempotent. */
	err = lis_solver_reset_sequence(sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_RESET_2 failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(9,9,"SEQUENCE_RESET_2");
	if( failed ) goto cleanup;


    /*
     * A value-generation change with a PSD-capable built-in
     * preconditioner must perform a numeric update without replacing
     * the persistent preconditioner object.
     */
    {
        LIS_PRECON update_precon0;
        LIS_UNSIGNED_INT structure0;
        LIS_UNSIGNED_INT value0;
        LIS_UNSIGNED_INT structure1;
        LIS_UNSIGNED_INT value1;
        LIS_SEQUENCE_DIAGNOSTICS diag;
        void *user_data;
        int user_sentinel;

        user_sentinel = 12345;
        user_data = NULL;

        err = lis_solver_set_option(
                "-i gmres -p ilu -ilu_fill 0 "
                "-scale none -print none "
                "-restart 4 -maxiter 50 "
                "-maxiter_noimp 0 -tol 1e-12 "
                "-initx_zeros true",
                sequence_solver);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE option setup failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }

        lis_vector_set_all((LIS_SCALAR)1.0,b);
        lis_vector_set_all((LIS_SCALAR)0.0,x);

        err = lis_solve_sequence(
                A0,b,x,sequence_solver);
        if( err!=LIS_SUCCESS ||
            sequence_solver->retcode!=LIS_SUCCESS )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE first solve failed "
                    "api=%d status=%d\n",
                    (int)err,
                    (int)sequence_solver->retcode);
            failed = 1;
            goto cleanup;
        }

        update_precon0 =
                sequence_solver->sequence_precon;

        if( update_precon0==NULL )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE no persistent "
                    "preconditioner after first solve\n");
            failed = 1;
            goto cleanup;
        }

        err = lis_matrix_get_generation(
                A0,&structure0,&value0);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE generation query 0 "
                    "failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }

        /*
         * Mark the outer preconditioner object so an accidental
         * destroy/create cycle cannot masquerade as pointer reuse.
         */
        err = lis_precon_set_user_data(
                update_precon0,
                (void *)&user_sentinel);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE user-data setup "
                    "failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }

        err = lis_matrix_mark_values_changed(A0);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE mark-values "
                    "failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }

        err = lis_matrix_get_generation(
                A0,&structure1,&value1);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE generation query 1 "
                    "failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }

        if( structure1!=structure0 ||
            value1!=value0+1 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE bad generation "
                    "before=%llu/%llu after=%llu/%llu\n",
                    (unsigned long long)structure0,
                    (unsigned long long)value0,
                    (unsigned long long)structure1,
                    (unsigned long long)value1);
            failed = 1;
            goto cleanup;
        }

        lis_vector_set_all((LIS_SCALAR)2.0,b);
        lis_vector_set_all((LIS_SCALAR)0.0,x);

        err = lis_solve_sequence(
                A0,b,x,sequence_solver);
        if( err!=LIS_SUCCESS ||
            sequence_solver->retcode!=LIS_SUCCESS )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE second solve failed "
                    "api=%d status=%d\n",
                    (int)err,
                    (int)sequence_solver->retcode);
            failed = 1;
            goto cleanup;
        }

        if( sequence_solver->sequence_precon!=
            update_precon0 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE replaced "
                    "persistent preconditioner\n");
            failed = 1;
            goto cleanup;
        }

        err = lis_precon_get_user_data(
                sequence_solver->sequence_precon,
                &user_data);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE user-data query "
                    "failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }

        if( user_data!=(void *)&user_sentinel )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE preconditioner "
                    "identity marker was lost\n");
            failed = 1;
            goto cleanup;
        }

        if( sequence_solver->
                sequence_precon_structure_generation!=
                    structure1 ||
            sequence_solver->
                sequence_precon_value_generation!=
                    value1 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE cached generation "
                    "is stale\n");
            failed = 1;
            goto cleanup;
        }

        err = lis_solver_get_sequence_diagnostics(
                sequence_solver,&diag);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE diagnostics "
                    "failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }

        if( diag.precon_rebuilds!=0 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE unexpected recovery "
                    "rebuilds=%d\n",
                    (int)diag.precon_rebuilds);
            failed = 1;
            goto cleanup;
        }

        printf(
                "SEQUENCE_VALUE_UPDATE PASSED "
                "generation=%llu/%llu->%llu/%llu "
                "precon_preserved=1 rebuilds=%d\n",
                (unsigned long long)structure0,
                (unsigned long long)value0,
                (unsigned long long)structure1,
                (unsigned long long)value1,
                (int)diag.precon_rebuilds);

        err = lis_solver_reset_sequence(
                sequence_solver);
        if( err )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_VALUE_UPDATE reset failed: %d\n",
                    (int)err);
            failed = 1;
            goto cleanup;
        }
    }

    printf("LIS_SEQUENCE_VALUE_UPDATE_S10 PASSED\n");

    /*
     * A nonzero sequence preconditioner max age permits bounded
     * numerical reuse across value-generation changes.
     */
    {
        LIS_SOLVER age_solver;
        LIS_PRECON age_precon;
        LIS_UNSIGNED_INT sg0,vg0;
        LIS_UNSIGNED_INT sg1,vg1;
        LIS_INT configured_age;

        age_solver = NULL;
        age_precon = NULL;

        err = lis_solver_create(&age_solver);
        if( err )
        {
            failed = 1;
            goto cleanup;
        }

        err = lis_solver_set_option(
                "-i gmres -p ilu -ilu_fill 0 "
                "-scale none -print none "
                "-restart 4 -maxiter 50 "
                "-maxiter_noimp 0 -tol 1e-12 "
                "-initx_zeros true",
                age_solver);
        if( err )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        /*
         * Newly created solvers must preserve legacy eager
         * preconditioner updates until max-age reuse is enabled.
         */
        configured_age = -1;

        err = lis_solver_get_sequence_precon_max_age(
                age_solver,&configured_age);

        if( err ||
            configured_age!=0 ||
            age_solver->sequence_precon_age!=0 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_PRECON_MAX_AGE "
                    "default initialization failed "
                    "max_age=%d age=%d err=%d\n",
                    (int)configured_age,
                    (int)age_solver->sequence_precon_age,
                    (int)err);

            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        err = lis_solver_set_sequence_precon_max_age(
                age_solver,2);
        if( err )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        configured_age = -1;

        err = lis_solver_get_sequence_precon_max_age(
                age_solver,&configured_age);
        if( err || configured_age!=2 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_PRECON_MAX_AGE config failed "
                    "value=%d err=%d\n",
                    (int)configured_age,
                    (int)err);
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        lis_vector_set_all((LIS_SCALAR)1.0,b);
        lis_vector_set_all((LIS_SCALAR)0.0,x);

        err = lis_solve_sequence(
                A0,b,x,age_solver);

        if( err!=LIS_SUCCESS ||
            age_solver->retcode!=LIS_SUCCESS )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        age_precon = age_solver->sequence_precon;

        err = lis_matrix_get_generation(
                A0,&sg0,&vg0);

        if( err )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        if( age_solver->sequence_precon_age!=0 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_PRECON_MAX_AGE initial age=%d\n",
                    (int)age_solver->sequence_precon_age);
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        /* First changed generation: reuse stale preconditioner. */
        err = lis_matrix_mark_values_changed(A0);
        if( err )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        lis_vector_set_all((LIS_SCALAR)2.0,b);
        lis_vector_set_all((LIS_SCALAR)0.0,x);

        err = lis_solve_sequence(
                A0,b,x,age_solver);

        if( err!=LIS_SUCCESS ||
            age_solver->retcode!=LIS_SUCCESS ||
            age_solver->sequence_precon!=age_precon ||
            age_solver->sequence_precon_age!=1 ||
            age_solver->sequence_precon_value_generation!=vg0 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_PRECON_MAX_AGE reuse1 failed "
                    "age=%d factor_gen=%llu expected=%llu\n",
                    (int)age_solver->sequence_precon_age,
                    (unsigned long long)
                        age_solver->sequence_precon_value_generation,
                    (unsigned long long)vg0);
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        /* Second changed generation: one more bounded stale reuse. */
        err = lis_matrix_mark_values_changed(A0);
        if( err )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        lis_vector_set_all((LIS_SCALAR)3.0,b);
        lis_vector_set_all((LIS_SCALAR)0.0,x);

        err = lis_solve_sequence(
                A0,b,x,age_solver);

        if( err!=LIS_SUCCESS ||
            age_solver->retcode!=LIS_SUCCESS ||
            age_solver->sequence_precon!=age_precon ||
            age_solver->sequence_precon_age!=2 ||
            age_solver->sequence_precon_value_generation!=vg0 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_PRECON_MAX_AGE reuse2 failed "
                    "age=%d factor_gen=%llu expected=%llu\n",
                    (int)age_solver->sequence_precon_age,
                    (unsigned long long)
                        age_solver->sequence_precon_value_generation,
                    (unsigned long long)vg0);
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        /*
         * Third changed generation: max age reached, so perform
         * the PSD numeric update and refresh the generation snapshot.
         */
        err = lis_matrix_mark_values_changed(A0);
        if( err )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        err = lis_matrix_get_generation(
                A0,&sg1,&vg1);

        if( err )
        {
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        lis_vector_set_all((LIS_SCALAR)4.0,b);
        lis_vector_set_all((LIS_SCALAR)0.0,x);

        err = lis_solve_sequence(
                A0,b,x,age_solver);

        if( err!=LIS_SUCCESS ||
            age_solver->retcode!=LIS_SUCCESS ||
            age_solver->sequence_precon!=age_precon ||
            age_solver->sequence_precon_age!=0 ||
            age_solver->sequence_precon_value_generation!=vg1 )
        {
            fprintf(
                    stderr,
                    "SEQUENCE_PRECON_MAX_AGE refresh failed "
                    "age=%d factor_gen=%llu current=%llu\n",
                    (int)age_solver->sequence_precon_age,
                    (unsigned long long)
                        age_solver->sequence_precon_value_generation,
                    (unsigned long long)vg1);
            failed = 1;
            lis_solver_destroy(age_solver);
            goto cleanup;
        }

        printf(
                "SEQUENCE_PRECON_MAX_AGE PASSED "
                "max_age=2 generation=%llu->%llu "
                "age_sequence=1,2,0 precon_preserved=1\n",
                (unsigned long long)vg0,
                (unsigned long long)vg1);

        lis_solver_destroy(age_solver);
    }

    printf("LIS_SEQUENCE_PRECON_MAX_AGE_S11 PASSED\n");




	/*
	 * The sequence solve operator and the matrix used to build the
	 * persistent preconditioner may be distinct.
	 *
	 * Preconditioner generation tracking must follow precon_A, while
	 * Krylov, residual certification, and solver state continue to use A.
	 */
	{
		LIS_SOLVER split_solver;
		LIS_PRECON split_precon;
		LIS_UNSIGNED_INT psg0,pvg0;
		LIS_UNSIGNED_INT psg1,pvg1;

		split_solver = NULL;
		split_precon = NULL;

		err = lis_solver_create(&split_solver);
		if( err )
		{
			failed = 1;
			goto cleanup;
		}

		err = lis_solver_set_option(
			"-i gmres -p ilu -ilu_fill 0 "
			"-scale none -print none "
			"-restart 4 -maxiter 50 "
			"-maxiter_noimp 0 -tol 1e-12 "
			"-initx_zeros true",
			split_solver);
		if( err )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"option setup failed: %d\n",
				(int)err);
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		err = lis_matrix_get_generation(
			A1,&psg0,&pvg0);
		if( err )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"initial generation query failed: %d\n",
				(int)err);
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		lis_vector_set_all((LIS_SCALAR)1.0,b);
		lis_vector_set_all((LIS_SCALAR)0.0,x);

		err = lis_solve_sequence_with_precon_matrix(
			A0,A1,b,x,split_solver);

		if( err!=LIS_SUCCESS ||
		    split_solver->retcode!=LIS_SUCCESS )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"first solve failed api=%d status=%d\n",
				(int)err,
				(int)split_solver->retcode);
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		if( split_solver->sequence_precon==NULL ||
		    split_solver->sequence_precon_matrix!=A1 ||
		    split_solver->A!=A0 )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"wrong solve/preconditioner ownership\n");
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		if( split_solver->
				sequence_precon_structure_generation!=psg0 ||
		    split_solver->
				sequence_precon_value_generation!=pvg0 )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"initial preconditioner generation mismatch\n");
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		split_precon =
			split_solver->sequence_precon;

		/*
		 * Changing only the solve-operator value generation must
		 * not invalidate the preconditioner matrix snapshot.
		 */
		err = lis_matrix_mark_values_changed(A0);
		if( err )
		{
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		lis_vector_set_all((LIS_SCALAR)2.0,b);
		lis_vector_set_all((LIS_SCALAR)0.0,x);

		err = lis_solve_sequence_with_precon_matrix(
			A0,A1,b,x,split_solver);

		if( err!=LIS_SUCCESS ||
		    split_solver->retcode!=LIS_SUCCESS ||
		    split_solver->sequence_precon!=split_precon ||
		    split_solver->sequence_precon_matrix!=A1 ||
		    split_solver->
				sequence_precon_structure_generation!=psg0 ||
		    split_solver->
				sequence_precon_value_generation!=pvg0 )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"solve-A generation incorrectly "
				"invalidated preconditioner state\n");
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		/*
		 * Changing precon_A must be observed by the persistent
		 * preconditioner lifecycle.  With the default max age of
		 * zero, a PSD-capable ILU performs an eager numeric update.
		 */
		err = lis_matrix_mark_values_changed(A1);
		if( err )
		{
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		err = lis_matrix_get_generation(
			A1,&psg1,&pvg1);
		if( err )
		{
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		if( psg1!=psg0 || pvg1!=pvg0+1 )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"unexpected preconditioner generation "
				"%llu/%llu -> %llu/%llu\n",
				(unsigned long long)psg0,
				(unsigned long long)pvg0,
				(unsigned long long)psg1,
				(unsigned long long)pvg1);
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		lis_vector_set_all((LIS_SCALAR)3.0,b);
		lis_vector_set_all((LIS_SCALAR)0.0,x);

		err = lis_solve_sequence_with_precon_matrix(
			A0,A1,b,x,split_solver);

		if( err!=LIS_SUCCESS ||
		    split_solver->retcode!=LIS_SUCCESS )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"preconditioner update solve failed "
				"api=%d status=%d\n",
				(int)err,
				(int)split_solver->retcode);
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		if( split_solver->sequence_precon!=split_precon ||
		    split_solver->sequence_precon_matrix!=A1 ||
		    split_solver->A!=A0 ||
		    split_solver->
				sequence_precon_structure_generation!=psg1 ||
		    split_solver->
				sequence_precon_value_generation!=pvg1 )
		{
			fprintf(
				stderr,
				"SEQUENCE_SEPARATE_PRECON "
				"preconditioner update ownership/generation "
				"failed\n");
			failed = 1;
			lis_solver_destroy(split_solver);
			goto cleanup;
		}

		printf(
			"SEQUENCE_SEPARATE_PRECON PASSED "
			"solve_matrix_distinct=1 "
			"precon_generation=%llu->%llu "
			"precon_preserved=1\n",
			(unsigned long long)pvg0,
			(unsigned long long)pvg1);

		lis_solver_destroy(split_solver);
	}

	printf("LIS_SEQUENCE_SEPARATE_PRECON_MATRIX PASSED\n");

	printf("LIS_SEQUENCE_SOLVE_S5 PASSED\n");

cleanup:

	if( sequence_solver )
	{
		lis_solver_destroy(sequence_solver);
		sequence_solver = NULL;
	}

	if( x )
	{
		lis_vector_destroy(x);
		x = NULL;
	}

	if( b )
	{
		lis_vector_destroy(b);
		b = NULL;
	}

	if( legacy_solver )
	{
		lis_solver_destroy(legacy_solver);
		legacy_solver = NULL;
	}

	if( solver )
	{
		lis_solver_destroy(solver);
		solver = NULL;
	}

	if( A1 ) lis_matrix_destroy(A1);
	if( A0 ) lis_matrix_destroy(A0);

	printf(
		"LIS_SEQUENCE_PRECON_S4_S5 %s "
		"create=%d destroy=%d\n",
		failed ? "FAILED" : "PASSED",
		create_count,
		destroy_count);

	lis_finalize();

	return failed ? 1 : 0;
}
