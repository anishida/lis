#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_solver_options.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Full-fill reference: nonsymmetric, strongly row/column imbalanced.
 * The test exercises A, A^T, RCM, row-only and two-sided Equil.
 */
#define N 5
static const LIS_REAL core[N][N]={
    {7.,2.,0.,0.,1.},
    {1.,8.,1.,0.,0.},
    {0.,1.,9.,2.,0.},
    {0.,0.,1.,10.,1.},
    {2.,0.,0.,1.,6.}
};
/* Condition number ~2e10: ill-scaled without placing exact-solve tests
 * at the ~1e16 double-precision limit of the initial stress fixture. */
static const LIS_REAL row_factor[N]={1.e-4,1.e3,1.e-1,10.,1.e2};
static const LIS_REAL col_factor[N]={1.e3,1.e-2,10.,1.e-1,1.e3};
static const LIS_REAL truth[N]={0.5,-1.,2.,0.25,-0.75};

static LIS_REAL entry(LIS_INT i,LIS_INT j)
{
    return row_factor[i]*core[i][j]*col_factor[j];
}

static LIS_INT construct(LIS_MATRIX *A)
{
    LIS_INT err,i,j;
    *A=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,A);if(err)return err;
    err=lis_matrix_set_size(*A,0,N);if(err)return err;
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            if(core[i][j]!=0.)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,j,
                                          (LIS_SCALAR)entry(i,j),*A);
                if(err)return err;
            }
    err=lis_matrix_set_type(*A,LIS_MATRIX_CSR);if(err)return err;
    return lis_matrix_assemble(*A);
}

static int run_case(LIS_INT rcm,LIS_INT row_scale,LIS_INT equil)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,bt=NULL,x=NULL;
    LIS_INT err=0,i,j,k;
    LIS_REAL maxv,scaled,diff;
    int fail=1;
    char opt[256];

    err=construct(&A);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&bt);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    for(i=0;i<N;i++)
    {
        b->value[i]=bt->value[i]=0.;
        for(j=0;j<N;j++)
        {
            b->value[i]+=entry(i,j)*truth[j];
            bt->value[i]+=entry(j,i)*truth[j];
        }
    }
    err=lis_solver_create(&solver);if(err)goto cleanup;
    solver->A=A;
    snprintf(opt,sizeof(opt),
        "-p ilut -ilut_rcm %d -ilut_row_scale %d -ilut_equil %d "
        "-ilut_drop 0 -ilut_rate 100 -ilu_pivot_tol 0 "
        "-scale none -omp_num_threads 1",
        (int)rcm,(int)row_scale,(int)equil);
    err=lis_solver_set_option(opt,solver);if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);if(err){pc=NULL;goto cleanup;}
    if(!pc->D || !pc->L || !pc->U)goto cleanup;
    if(equil)
    {
        if(!pc->ilut_col_scale || !pc->ilut_perm || !pc->ilut_row_scale)
            goto cleanup;
        /* Check actual two-sided equilibration: every nonempty
           scaled column must have maximum absolute entry ~1. */
        for(j=0;j<N;j++)
        {
            maxv=0.;
            for(i=0;i<N;i++)
            {
                scaled=fabs(entry(pc->ilut_perm[i],
                                  pc->ilut_perm[j]) *
                            pc->ilut_row_scale[i] *
                            pc->ilut_col_scale[j]);
                if(scaled>maxv)maxv=scaled;
            }
            if(fabs(maxv-1.)>1.e-10)
            {
                fprintf(stderr,"Equil column %d max %.17g\n",
                        (int)j,(double)maxv);
                goto cleanup;
            }
        }
    }
    else if(pc->ilut_col_scale!=NULL)goto cleanup;

    solver->precon=pc;
    err=lis_psolve_ilut_csr(solver,b,x);if(err)goto cleanup;
    for(k=0;k<N;k++)
    {
        diff=fabs(x->value[k]-truth[k]);
        if(!isfinite((double)diff) || diff>1.e-7)
        {
            fprintf(stderr,"forward mismatch %d: %.12g\n",
                    (int)k,(double)diff);
            goto cleanup;
        }
    }

    err=lis_psolveh_ilut_csr(solver,bt,x);if(err)goto cleanup;
    for(k=0;k<N;k++)
    {
        diff=fabs(x->value[k]-truth[k]);
        if(!isfinite((double)diff) || diff>1.e-7)
        {
            fprintf(stderr,"transpose mismatch %d: %.12g\n",
                    (int)k,(double)diff);
            goto cleanup;
        }
    }
    fail=0;

cleanup:
    if(fail)
        fprintf(stderr,"ILUT Equil failure rcm=%d row=%d equil=%d err=%d\n",
                (int)rcm,(int)row_scale,(int)equil,(int)err);
    if(solver)solver->precon=NULL;
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(bt)lis_vector_destroy(bt);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return fail;
}

int main(int argc,char **argv)
{
    LIS_INT err,rcm,row,eq,changed;
    LIS_SOLVER solver=NULL;
    int failed=0;

    err=lis_initialize(&argc,&argv);if(err)return 1;
#ifdef _OPENMP
    /* Direct precon_create bypasses the solver's thread dispatcher. */
    omp_set_dynamic(0);
    omp_set_num_threads(1);
#endif
#if defined(_COMPLEX) || defined(USE_COMPLEX)
    puts("testilutequil: SKIP (real-double-only option)");
#else
    for(rcm=0;rcm<=1;rcm++)
        for(row=0;row<=1;row++)
            for(eq=0;eq<=1;eq++)
                failed|=run_case(rcm,row,eq);

    err=lis_solver_create(&solver);
    if(err)failed=1;
    else
    {
        if(lis_solver_set_option("-p ilut -ilut_equil 1",solver)!=LIS_SUCCESS ||
           !solver->ilut_equil ||
           lis_solver_set_option("-ilut_equil 2",solver)==LIS_SUCCESS)
            failed=1;
        err=lis_solver_set_option_replace_if_changed(
            "-p ilut -ilut_equil 1 -ilut_rcm 1",solver,&changed);
        if(err || changed!=LIS_TRUE ||
           solver->ilut_equil!=1 || solver->ilut_rcm!=1)failed=1;
        err=lis_solver_set_option_replace_if_changed(
            "-p ilut -ilut_equil 1 -ilut_rcm 1",solver,&changed);
        if(err || changed!=LIS_FALSE)failed=1;
        err=lis_solver_set_option_replace_if_changed(
            "-p ilut -ilut_equil 0 -ilut_row_scale 1",solver,&changed);
        if(err || changed!=LIS_TRUE ||
           solver->ilut_equil!=0 || solver->ilut_row_scale!=1)
            failed=1;
        lis_solver_destroy(solver);
    }
#endif
    lis_finalize();
    if(!failed)puts("testilutequil: PASS");
    return failed;
}
