/*
 * Reproducible nonsymmetric convection-diffusion Equil experiment.
 * A(i,j) = stencil(i,j) * variable_units(j): a manufactured change
 * of the units of adjacent unknowns. Compare ordinary ILUT, Equil, ILUTP, and their combination.
 * True residual always uses original, unscaled A.
 *
 * Usage: ./test/benchilutpequil [unit_amplitude] [fill_rate] [drop]
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define SIDE 12
#define N (SIDE*SIDE)
static LIS_REAL units(LIS_INT j,LIS_REAL amp)
{
    LIS_INT x=j%SIDE,y=j/SIDE;
    return ((x+y)&1) ? (LIS_REAL)(1.0/amp) : amp;
}
static LIS_REAL coeff(LIS_INT i,LIS_INT j)
{
    LIS_INT x=i%SIDE,y=i/SIDE;
    if(i==j)return 5.;
    if(x+1<SIDE && j==i+1)return -1.65;
    if(x>0 && j==i-1)return -0.35;
    if(y+1<SIDE && j==i+SIDE)return -1.;
    if(y>0 && j==i-SIDE)return -1.;
    return 0.;
}
static int run_case(LIS_INT mode,LIS_REAL amp,LIS_REAL rate,LIS_REAL drop)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,ax=NULL;
    LIS_SOLVER s=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err=0,i,j,iter=-1;
    char options[450];
    LIS_REAL nr=0.,nb=0.,rr,actual_rate,actual_drop;
    double elapsed=0.,solve_seconds=0.,setup_seconds=0.,create_seconds=0.,init_seconds=0.;
    LIS_INT swaps=0;
    long long nnz_factor=0;
    const char *name=mode==0?"ilut":mode==1?"equil":
                     mode==2?"ilutp":mode==3?"ilutp_equil":"row_only";
    int failed=1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,N);if(err)goto done;
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            LIS_REAL v=coeff(i,j);
            if(v!=0.)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,j,
                                         (LIS_SCALAR)(v*units(j,amp)),A);
                if(err)goto done;
            }
        }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto done;
    err=lis_matrix_assemble(A);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&ax);if(err)goto done;
    for(i=0;i<N;i++)
    {
        LIS_REAL sum=0.;
        for(j=0;j<N;j++)sum+=coeff(i,j)*units(j,amp)*
                                 (LIS_REAL)cos((double)j*0.24);
        b->value[i]=(LIS_SCALAR)sum;
    }
    err=lis_solver_create(&s);if(err)goto done;
    snprintf(options,sizeof(options),
        "-i fgmres -restart 50 -p ilut -ilut_rcm 0 "
        "-ilut_row_scale %d -ilut_equil %d -ilutp %d "
        "-ilutp_threshold 0.5 "
        "-ilut_drop %.9g -ilut_rate %.7g -ilu_pivot_tol 0.001 "
        "-scale none -tol 1e-6 -maxiter 300 -omp_num_threads 1",
        (int)(mode==4),(int)(mode==1 || mode==3),
        (int)(mode==2 || mode==3),(double)drop,(double)rate);
    err=lis_solver_set_option(options,s);if(err)goto done;
    actual_rate=s->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
    actual_drop=s->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
    /* Explicitly collect fill: lis_solve destroys its preconditioner. */
    s->A=A;
    err=lis_precon_create(s,&pc);if(err){pc=NULL;goto done;}
    if(!pc->L || !pc->U || !pc->D)goto done;
    swaps=pc->ilutp_swap_count;
    for(i=0;i<N;i++)
        nnz_factor+=(long long)pc->L->nnz[i]+
                     (long long)pc->U->nnz[i];
    lis_precon_destroy(pc);pc=NULL;
    err=lis_solve(A,b,x,s);
    if(lis_solver_get_iter(s,&iter))iter=-1;
    (void)lis_solver_get_timeex(s,&elapsed,&solve_seconds,
                                &setup_seconds,&create_seconds,&init_seconds);
    if(lis_matvec(A,x,ax))goto done;
    for(i=0;i<N;i++)
    {
        LIS_REAL d=(LIS_REAL)(ax->value[i]-b->value[i]);
        nr+=d*d;
        nb+=(LIS_REAL)(b->value[i]*b->value[i]);
    }
    rr=nb>0. ? (LIS_REAL)sqrt((double)(nr/nb)) : 0.;
    printf("ILUTP_EQUIL_BENCH mode=%s n=%d amplitude=%.9g rate=%.7g drop=%.9g actual_rate=%.7g actual_drop=%.7g nnz_A=%lld lfil_budget=%d "
           "status=%d iterations=%d true_residual=%.12e "
           "nnz_L_plus_U=%lld swaps=%d wall_seconds=%.6g "
           "precon_seconds=%.6g solve_seconds=%.6g\n",name,N,(double)amp,
           (double)rate,(double)drop,(double)actual_rate,(double)actual_drop,
           (long long)A->nnz,
           (int)((double)A->nnz/(2.0*N)*actual_rate),
           (int)err,(int)iter,(double)rr,nnz_factor,
           (int)swaps,elapsed,setup_seconds,solve_seconds);
    fflush(stdout);
    failed=0;
done:
    if(failed)fprintf(stderr,"ILUTP_EQUIL_BENCH error mode=%s err=%d\n",
                      name,(int)err);
    if(pc)lis_precon_destroy(pc);
    if(s)lis_solver_destroy(s);
    if(ax)lis_vector_destroy(ax);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}
int main(int argc,char **argv)
{
    LIS_REAL amp=100.,rate=1.,drop=0.005;
    LIS_INT err;
    int failed=0;
    if(argc>1)amp=(LIS_REAL)atof(argv[1]);
    if(argc>2)rate=(LIS_REAL)atof(argv[2]);
    if(argc>3)drop=(LIS_REAL)atof(argv[3]);
    if(!(amp>=1. && amp<=1000.) || !(rate>0. && rate<=10.) ||
       !(drop>=0. && drop<=1.))return 2;
    err=lis_initialize(&argc,&argv);if(err)return 1;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(1);
#endif
#if defined(_COMPLEX) || defined(USE_COMPLEX)
    puts("ILUTP_EQUIL_BENCH skip: real-double only");
#else
    failed|=run_case(0,amp,rate,drop); /* classic ILUT */
    failed|=run_case(1,amp,rate,drop); /* Equil alone */
    failed|=run_case(2,amp,rate,drop); /* ILUTP alone */
    failed|=run_case(3,amp,rate,drop); /* Equil + ILUTP */
    failed|=run_case(4,amp,rate,drop); /* row-only reference */
#endif
    lis_finalize();
    return failed;
}
