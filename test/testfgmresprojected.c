/*
 * Regression for the relative projected FGMRES residual.
 *
 * Verify first-step projected residual reporting against an independent
 * two-dimensional Arnoldi/Givens calculation. Do not assert an exact
 * iteration count: that depends on the orthogonalization implementation.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#include <stdio.h>
#include "lislib.h"

#ifdef _COMPLEX
int main(void)
{
    puts("SKIP: real-valued FGMRES projected-residual regression");
    return 77;
}
#else

struct projected_sample
{
    LIS_INT seen;
    LIS_REAL first;
};

static void projected_monitor(
    LIS_SOLVER solver, LIS_INT iter, LIS_REAL residual,
    LIS_REAL tolerance, LIS_INT kind, void *user_data)
{
    struct projected_sample *sample=(struct projected_sample *)user_data;
    (void)solver;
    (void)iter;
    (void)tolerance;
    if(kind==LIS_SOLVER_RESIDUAL_PROJECTED && !sample->seen)
    {
        sample->seen=LIS_TRUE;
        sample->first=residual;
    }
}

static LIS_REAL expected_projected(int warm)
{
    /* b=(1,1), A=diag(2,3); x0=(0,0.01) if warm.
     * For supported NRM2_R, normalize by ||r0||. */
    LIS_REAL r0=1.0, r1=warm ? 0.97 : 1.0;
    LIS_REAL denom=sqrt(r0*r0+r1*r1);
    LIS_REAL v0=r0/denom, v1=r1/denom;
    LIS_REAL av0=2.0*v0, av1=3.0*v1;
    LIS_REAL h00=v0*av0+v1*av1;
    LIS_REAL h10=hypot(av0-h00*v0,av1-h00*v1);
    return h10/hypot(h00,h10);
}

static int run_case(
    LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_VECTOR work,
    LIS_REAL magnitude, int warm, int quad)
{
    LIS_SOLVER solver=NULL;
    struct projected_sample sample;
    LIS_REAL initialnorm, actual, expected, expected_bnrm;
    LIS_REAL first_initial1=warm ? 0.97*magnitude : magnitude;
    LIS_INT rc;
    int failed=1;
    const char *opts=quad
        ? "-i fgmres -p none -scale none -print none -f quad "
          "-tol 0.3 -maxiter 20 -restart 2 -conv_cond nrm2_r -initx_zeros false"
        : "-i fgmres -p none -scale none -print none "
          "-tol 0.3 -maxiter 20 -restart 2 -conv_cond nrm2_r -initx_zeros false";
    sample.seen=LIS_FALSE;
    sample.first=-1.0;
    initialnorm=hypot(magnitude,first_initial1);

#define CHECK(call) do { rc=(call); if(rc!=LIS_SUCCESS) { \
    fprintf(stderr,"%s failed (%d)\n",#call,(int)rc); goto cleanup; \
} } while(0)

    CHECK(lis_vector_set_all(0.0,b));
    CHECK(lis_vector_set_value(LIS_INS_VALUE,0,magnitude,b));
    CHECK(lis_vector_set_value(LIS_INS_VALUE,1,magnitude,b));
    CHECK(lis_vector_set_all(0.0,x));
    if(warm)
        CHECK(lis_vector_set_value(LIS_INS_VALUE,1,0.01*magnitude,x));
    CHECK(lis_solver_create(&solver));
    CHECK(lis_solver_set_option((char *)opts,solver));
    CHECK(lis_solver_set_monitor(solver,projected_monitor,&sample));
    CHECK(lis_solve(A,b,x,solver));

    expected_bnrm=1.0/initialnorm;
    expected=expected_projected(warm);
    if(!sample.seen || !isfinite(sample.first) ||
       fabs(sample.first-expected)>1.0e-8 ||
       fabs(solver->bnrm/expected_bnrm-1.0)>1.0e-10)
    {
        fprintf(stderr,"Projected mismatch warm=%d scale=%.3e "
            "observed=%.12e expected=%.12e bnrm=%.12e expected_bnrm=%.12e quad=%d\n",
            warm,(double)magnitude,(double)sample.first,
            (double)expected,(double)solver->bnrm,(double)expected_bnrm,quad);
        goto cleanup;
    }

    CHECK(lis_matvec(A,x,work));
    CHECK(lis_vector_xpay(b,-1.0,work));
    CHECK(lis_vector_nrm2(work,&actual));
    actual/=initialnorm;
    if(!isfinite(actual) || actual>=0.3 || solver->retcode!=LIS_SUCCESS)
    {
        fprintf(stderr,"Physical residual failed warm=%d scale=%.3e "
            "relative=%.12e\n",warm,(double)magnitude,(double)actual);
        goto cleanup;
    }
    printf("FGMRES_PROJECTED_PASSED warm=%d quad=%d scale=%.3e "
        "projected=%.12e physical=%.12e\n",
        warm,quad,(double)magnitude,(double)sample.first,(double)actual);
    failed=0;

cleanup:
    if(solver) lis_solver_destroy(solver);
#undef CHECK
    return failed;
}

int main(int argc, char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,work=NULL;
    LIS_INT rc;
    int failed=1,warm,i,quad;
    const LIS_REAL scales[]={1.0e4,1.0e-2};
    rc=lis_initialize(&argc,&argv);
    if(rc!=LIS_SUCCESS) return 1;
#ifdef USE_MPI
    {
        int ranks=1;
        MPI_Comm_size(LIS_COMM_WORLD,&ranks);
        if(ranks!=1) { lis_finalize(); return 77; }
    }
#endif
#define CHECK(call) do { rc=(call); if(rc!=LIS_SUCCESS) { \
    fprintf(stderr,"%s failed (%d)\n",#call,(int)rc); goto cleanup; \
} } while(0)
    CHECK(lis_matrix_create(LIS_COMM_WORLD,&A));
    CHECK(lis_matrix_set_size(A,0,2));
    CHECK(lis_matrix_set_value(LIS_INS_VALUE,0,0,2.0,A));
    CHECK(lis_matrix_set_value(LIS_INS_VALUE,1,1,3.0,A));
    CHECK(lis_matrix_assemble(A));
    CHECK(lis_vector_duplicate(A,&b));
    CHECK(lis_vector_duplicate(A,&x));
    CHECK(lis_vector_duplicate(A,&work));
    for(quad=0;quad<=1;quad++)
    {
#ifndef USE_QUAD_PRECISION
        if(quad) break;
#endif
        for(warm=0;warm<=1;warm++)
            for(i=0;i<2;i++)
                if(run_case(A,b,x,work,scales[i],warm,quad))
                    goto cleanup;
    }
    failed=0;
cleanup:
    if(work) lis_vector_destroy(work);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
#undef CHECK
}
#endif
