/* Regression: opt-in first-level multivector aggregation.
 * Compare the default and larger fine-level targets with 2 and 3 modes.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#include <stdio.h>
#include "lislib.h"
#include "lis_solver_options.h"

#define TEST_N 240
#define TEST_LEVEL_CAP 32

extern LIS_INT lis_precon_mfamg_get_hierarchy_info(
    LIS_PRECON precon, LIS_INT max_levels, LIS_INT *nlevels,
    LIS_INT *coarse_n, LIS_INT *level_n, LIS_INT *owns_A);

static int check_profile(void)
{
    LIS_SOLVER solver=NULL;
    LIS_INT err,changed=0;
    int failed=1;
    const char *two=
        "-i gmres -p mfamg -scale none -print none -mfamg_fine_agg_factor 2";
    const char *three=
        "-i gmres -p mfamg -scale none -print none -mfamg_fine_agg_factor 3";

    err=lis_solver_create(&solver);
    if(err) goto done;
    if(solver->mfamg_fine_agg_factor!=0) goto done;
    err=lis_solver_set_option_replace(two,solver);
    if(err || solver->mfamg_fine_agg_factor!=2) goto done;
    err=lis_solver_set_option_replace_if_changed(two,solver,&changed);
    if(err || changed!=LIS_FALSE) goto done;
    err=lis_solver_set_option_replace_if_changed(three,solver,&changed);
    if(err || changed!=LIS_TRUE ||
       solver->mfamg_fine_agg_factor!=3) goto done;
    err=lis_solver_set_option_replace_if_changed(two,solver,&changed);
    if(err || changed!=LIS_TRUE ||
       solver->mfamg_fine_agg_factor!=2) goto done;
    if(lis_solver_set_option("-mfamg_fine_agg_factor 5",solver)==LIS_SUCCESS ||
       solver->mfamg_fine_agg_factor!=2) goto done;
    failed=0;
done:
    if(failed) fprintf(stderr,"MFAMG FINE AGG: option-profile regression failed\n");
    if(solver) lis_solver_destroy(solver);
    return failed;
}

static int run_case(LIS_MATRIX A, LIS_VECTOR basis[3],
                    int mode_count, int factor, LIS_INT *fine_coarse)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_VECTOR b=NULL,x=NULL,ax=NULL;
    LIS_INT level_n[TEST_LEVEL_CAP],owned[TEST_LEVEL_CAP];
    LIS_INT nlevels=0,coarse_n=0,err=0,i;
    LIS_REAL b2=0.0,r2=0.0,bi,ri,defect;
    char option[180];
    int failed=1;

    err=lis_solver_create(&solver);
    if(err) goto done;
    if(solver->mfamg_fine_agg_factor!=0) goto done;
    if(factor==0)
        snprintf(option,sizeof(option),
                 "-i gmres -p mfamg -scale none -print none");
    else
        snprintf(option,sizeof(option),
                 "-i gmres -p mfamg -scale none -print none -mfamg_fine_agg_factor %d",
                 factor);
    err=lis_solver_set_option(option,solver);
    if(err || solver->mfamg_fine_agg_factor!=factor) goto done;
    err=lis_solver_set_near_nullspace(solver,mode_count,basis);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_precon_create(solver,&precon);
    if(err) goto done;
    err=lis_precon_mfamg_get_hierarchy_info(
        precon,TEST_LEVEL_CAP,&nlevels,&coarse_n,level_n,owned);
    if(err || nlevels<2 || level_n[0]!=TEST_N ||
       level_n[1]<=0 || level_n[1]>=TEST_N || coarse_n<=0)
    {
        fprintf(stderr,
                "MFAMG FINE AGG: hierarchy modes=%d factor=%d levels=%d fine_coarse=%d\n",
                mode_count,factor,(int)nlevels,(int)level_n[1]);
        goto done;
    }
    *fine_coarse=level_n[1];

    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&ax);
    if(err) goto done;
    for(i=0;i<TEST_N;i++)
        b->value[i]=1.0+((LIS_REAL)(i%7))/7.0;

    solver->precon=precon;
    err=lis_psolve(solver,b,x);
    if(err) goto done;
    err=lis_matvec(A,x,ax);
    if(err) goto done;
    for(i=0;i<TEST_N;i++)
    {
#ifdef _COMPLEX
        bi=(LIS_REAL)cabs(b->value[i]);
        ri=(LIS_REAL)cabs(b->value[i]-ax->value[i]);
#else
        bi=(LIS_REAL)fabs(b->value[i]);
        ri=(LIS_REAL)fabs(b->value[i]-ax->value[i]);
#endif
        b2+=bi*bi;
        r2+=ri*ri;
    }
    if(!isfinite((double)b2) || !isfinite((double)r2) ||
       b2<=0.0 || r2<0.0) goto done;
    defect=sqrt(r2/b2);
    if(!isfinite((double)defect) || defect>=0.99)
    {
        fprintf(stderr,
            "MFAMG FINE AGG: no physical defect reduction modes=%d factor=%d defect=%.9e\n",
            mode_count,factor,(double)defect);
        goto done;
    }
    printf("MFAMG FINE AGG modes=%d factor=%d level1=%d levels=%d defect=%.9e PASS\n",
           mode_count,factor,(int)level_n[1],(int)nlevels,(double)defect);
    failed=0;
done:
    if(err) fprintf(stderr,"MFAMG FINE AGG: LIS error=%d modes=%d factor=%d\n",
                    (int)err,mode_count,factor);
    if(solver) solver->precon=NULL;
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    return failed;
}

static int run_test(void)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR basis[3]={NULL,NULL,NULL};
    LIS_INT i,err=0;
    int mode_count,failed=1;
    LIS_INT coarse_default,coarse_two,coarse_three;

    if(check_profile()!=0) return 1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,TEST_N);
    if(err) goto done;
    for(i=0;i<TEST_N;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
            if(err) goto done;
        }
        if(i+1<TEST_N)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    for(i=0;i<3;i++)
    {
        err=lis_vector_duplicate(A,&basis[i]);
        if(err) goto done;
    }
    /*
     * Each block of eight unknowns has a constant, locally rank-one
     * near-nullspace. Across blocks, the two- and three-mode families
     * are globally independent. This permits an actual default
     * hierarchy, making its dimension comparable to opt-in targets.
     */
    for(i=0;i<TEST_N;i++)
    {
        LIS_REAL coordinate=((LIS_REAL)(i/8))/30.0;
        basis[0]->value[i]=1.0;
        basis[1]->value[i]=coordinate;
        basis[2]->value[i]=coordinate*coordinate;
    }
    for(mode_count=2;mode_count<=3;mode_count++)
    {
        if(run_case(A,basis,mode_count,0,&coarse_default)!=0 ||
           run_case(A,basis,mode_count,2,&coarse_two)!=0 ||
           run_case(A,basis,mode_count,3,&coarse_three)!=0) goto done;
        if(!(coarse_three<=coarse_two &&
             coarse_two<coarse_default &&
             coarse_two<=TEST_N/2))
        {
            fprintf(stderr,
                    "MFAMG FINE AGG: no coarsening gain modes=%d default=%d factor2=%d factor3=%d\n",
                    mode_count,(int)coarse_default,
                    (int)coarse_two,(int)coarse_three);
            goto done;
        }
    }
    failed=0;
done:
    if(err) fprintf(stderr,"MFAMG FINE AGG setup LIS error=%d\n",(int)err);
    for(i=0;i<3;i++) if(basis[i]) lis_vector_destroy(basis[i]);
    if(A) lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
    int failed;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#ifdef USE_MPI
    {
        int nprocs=1;
        MPI_Comm_size(LIS_COMM_WORLD,&nprocs);
        if(nprocs!=1)
        {
            lis_finalize();
            return 0; /* One-rank numerical hierarchy comparison. */
        }
    }
#endif
    failed=run_test();
    lis_finalize();
    return failed;
}
