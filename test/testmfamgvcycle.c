/* Isolated MFAMG per-solver V-cycle smoothing regression. */
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
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string.h>
#include "lislib.h"

#define N 96
#define MAX_LEVELS 32

extern LIS_INT lis_precon_mfamg_apply_block_smoother_for_test(
    LIS_MATRIX A, LIS_INT block_size, LIS_VECTOR b, LIS_VECTOR x,
    LIS_INT initial_zero, LIS_INT *active_blocks, LIS_INT *fallback_blocks,
    LIS_REAL *omega);

extern LIS_INT lis_precon_mfamg_get_omp_min_rows_for_test(void);
extern LIS_INT lis_precon_mfamg_get_multivector_levels_for_test(
    LIS_PRECON precon, LIS_INT *count);

static int mfamg_force_omp_test_threshold(const char *threshold)
{
    /* Set before workers start. Explicitly cover OpenMP loops even when
     * production is compiled with a much higher default threshold. */
#ifdef _WIN32
    return _putenv_s("LIS_MFAMG_OMP_MIN_ROWS",threshold);
#else
    return setenv("LIS_MFAMG_OMP_MIN_ROWS",threshold,1);
#endif
}

extern LIS_INT lis_precon_mfamg_get_hierarchy_info(
    LIS_PRECON precon, LIS_INT max_levels, LIS_INT *nlevels,
    LIS_INT *coarse_n, LIS_INT *level_n, LIS_INT *owns_A);

static int run_case(LIS_MATRIX A,LIS_VECTOR basis,
                    LIS_VECTOR b,LIS_VECTOR x,int expected, int explicit_option,
                    LIS_SCALAR *solution, LIS_REAL *relative_defect)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_VECTOR near[1],ax=NULL;
    LIS_INT err=0,i,nlevels=0,coarse_n=0;
    LIS_INT level_n[MAX_LEVELS],owns_A[MAX_LEVELS];
    LIS_REAL r2=0.0,b2=0.0,bi,ri;
    int failed=1;
    char opts[180];
    near[0]=basis;
    err=lis_solver_create(&solver);
    if(err) goto done;
    if(solver->mfamg_vcycle_smooth_steps!=1) goto done;
    if(explicit_option)
        snprintf(opts,sizeof(opts),
            "-i gmres -p mfamg -scale none -print none -mfamg_vcycle_smooth_steps %d",
            expected);
    else
        snprintf(opts,sizeof(opts),
            "-i gmres -p mfamg -scale none -print none");
    err=lis_solver_set_option(opts,solver);
    if(err || solver->mfamg_vcycle_smooth_steps!=expected) goto done;
    if(explicit_option &&
       (lis_solver_set_option("-mfamg_vcycle_smooth_steps 0",solver)==LIS_SUCCESS ||
        solver->mfamg_vcycle_smooth_steps!=expected)) goto done;
    err=lis_solver_set_near_nullspace(solver,1,near);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_precon_create(solver,&precon);
    if(err) {precon=NULL;goto done;}
    solver->precon=precon;
    err=lis_precon_mfamg_get_hierarchy_info(
        precon,MAX_LEVELS,&nlevels,&coarse_n,level_n,owns_A);
    if(err || nlevels<4 || level_n[0]!=N || coarse_n<=0)
    {
        fprintf(stderr,"MFAMG VCYCLE unexpected hierarchy nlevels=%d coarse=%d\n",
                (int)nlevels,(int)coarse_n);
        goto done;
    }
    err=lis_psolve(solver,b,x);
    if(err) goto done;
    err=lis_psolve(solver,b,x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&ax);
    if(err) goto done;
    err=lis_matvec(A,x,ax);
    if(err) goto done;
    for(i=0;i<N;i++) {
#ifdef _COMPLEX
        bi=(LIS_REAL)cabs(b->value[i]);
        ri=(LIS_REAL)cabs(b->value[i]-ax->value[i]);
#else
        bi=(LIS_REAL)fabs(b->value[i]);
        ri=(LIS_REAL)fabs(b->value[i]-ax->value[i]);
#endif
        r2+=ri*ri;
        b2+=bi*bi;
    }
    if(!isfinite((double)r2) || !isfinite((double)b2) ||
       b2<=0.0 || r2<0.0) goto done;
    *relative_defect = sqrt(r2/b2);
    if(!isfinite((double)*relative_defect) ||
       *relative_defect>=0.90)
    {
        fprintf(stderr,"MFAMG VCYCLE no useful defect reduction steps=%d rel=%.9e\n",
                expected,(double)*relative_defect);
        goto done;
    }
    memcpy(solution,x->value,(size_t)N*sizeof(LIS_SCALAR));
    printf("MFAMG VCYCLE steps=%d explicit=%d levels=%d relative_defect=%.9e PASS\n",
        expected,explicit_option,(int)nlevels,(double)*relative_defect);
    failed=0;
done:
    if(err) fprintf(stderr,"MFAMG VCYCLE test rc=%d steps=%d\n",
                    (int)err,expected);
    if(ax) lis_vector_destroy(ax);
    if(solver) solver->precon=NULL;
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    return failed;
}


/*
 * Execute MFAMG with enough rows to activate the parallel smoother and
 * prolongation loops.  Compare every output bit across 1,2,4,8 threads;
 * both weighted and block Jacobi use the same original matrix and RHS.
 * This intentionally tests one V-cycle (the preconditioner) rather than
 * hiding numerical discrepancies in the outer Krylov solver.
 */
static int mfamg_omp_stress(void)
{
#ifndef _OPENMP
    puts("MFAMG_OMP_STRESS SKIP: built without OpenMP");
    return 0;
#else
    enum { STRESS_N=4096 };
    const int threads[4]={1,2,4,8};
    /* Golden residuals measured from the unchanged develop commit
     * 14516697c38286e002c51c66867a8b1de9ea565a using
     * OMP_DYNAMIC=FALSE OMP_NUM_THREADS=8
     * ./test/testmfamgvcycle --omp-stress
     * on the baseline-only audit branch. */
    const LIS_REAL baseline_relres[2]={
        (LIS_REAL)0.3499046495,(LIS_REAL)0.4079550598
    };
    LIS_MATRIX A=NULL;
    LIS_VECTOR basis=NULL,b=NULL,x=NULL,Ax=NULL,near[1];
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_SCALAR *reference=NULL;
    LIS_REAL rel=0.0,normb,normr,diff;
    double setup_started=0,setup_seconds=0,apply_started=0,apply_seconds=0;
    LIS_INT err=0,i;
    int failed=1,mode=-1,t=-1,repeat=0;
    const char *stage="init";
    int original_threads=omp_get_max_threads();
    int original_dynamic=omp_get_dynamic();
    const char *options[2]={
        "-i gmres -p mfamg -scale none -print none "
        "-mfamg_vcycle_smooth_steps 2",
        "-i gmres -p mfamg -scale none -print none "
        "-mfamg_block_size 4 -mfamg_vcycle_smooth_steps 2"
    };
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    if(mfamg_force_omp_test_threshold("4096")!=0 ||
       lis_precon_mfamg_get_omp_min_rows_for_test()>STRESS_N)
    {
        fprintf(stderr,"MFAMG_OMP_STRESS parallel threshold not enabled\n");
        goto done;
    }
    reference=(LIS_SCALAR *)malloc((size_t)STRESS_N*sizeof(LIS_SCALAR));
    if(!reference) goto done;
    stage="matrix-create";
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,STRESS_N);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<STRESS_N;i++)
    {
        /* Near-Poisson coupling stays strong across many AMG levels. */
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)2.001,A);
        if(err) goto done;
        if(i>0) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
        if(i+1<STRESS_N) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
    }
    stage="matrix-set-type";
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    stage="matrix-assemble";
    err=lis_matrix_assemble(A);
    if(err) goto done;
    stage="vector-duplicate";
    if(lis_vector_duplicate(A,&basis) || lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) || lis_vector_duplicate(A,&Ax))
        goto done;
    for(i=0;i<STRESS_N;i++)
    {
        basis->value[i]=(LIS_SCALAR)1;
        b->value[i]=(LIS_SCALAR)(1.0+(i%7)*0.0625);
    }
    near[0]=basis;
    for(mode=0;mode<2;mode++)
    {
        for(t=0;t<4;t++)
        {
            if(omp_get_thread_limit()<threads[t])
            {
                printf("MFAMG_OMP_STRESS mode=%d threads=%d "
                       "SKIP: OMP_THREAD_LIMIT=%d\n",
                       mode,threads[t],omp_get_thread_limit());
                continue;
            }
            omp_set_num_threads(threads[t]);
            stage="solver-create";
            setup_started=omp_get_wtime();
            err=lis_solver_create(&solver);
            if(err) goto done;
            stage="solver-options";
            err=lis_solver_set_option((char *)options[mode],solver);
            if(err) goto done;
            stage="near-nullspace";
            err=lis_solver_set_near_nullspace(solver,1,near);
            if(err) goto done;
            stage="solver-matrix";
            err=lis_solver_set_matrix(A,solver);
            if(err) goto done;
            stage="preconditioner-create";
            err=lis_precon_create(solver,&precon);
            if(err) {precon=NULL;goto done;}
            solver->precon=precon;
            setup_seconds=omp_get_wtime()-setup_started;
            apply_started=omp_get_wtime();

            for(repeat=0;repeat<3;repeat++)
            {
                stage="psolve";
                err=lis_psolve(solver,b,x);
                if(err) goto done;
                err=lis_matvec(A,x,Ax);
                if(err) goto done;
                normr=0;
                normb=0;
                for(i=0;i<STRESS_N;i++)
                {
#ifdef _COMPLEX
                    diff=(LIS_REAL)cabs(Ax->value[i]-b->value[i]);
                    normb+=(LIS_REAL)pow(cabs(b->value[i]),2.0);
#else
                    diff=(LIS_REAL)fabs(Ax->value[i]-b->value[i]);
                    normb+=(LIS_REAL)(b->value[i]*b->value[i]);
#endif
                    normr+=diff*diff;
                }
                rel=sqrt(normr/normb);
                if(!isfinite((double)rel) || !(rel<0.60) ||
                   fabs(rel-baseline_relres[mode])>(LIS_REAL)1.0e-7)
                {
                    fprintf(stderr,
                        "MFAMG_OMP_STRESS bad residual mode=%d threads=%d "
                        "repeat=%d residual=%.8e\n",
                        mode,threads[t],repeat,(double)rel);
                    goto done;
                }
                if(t==0 && repeat==0)
                    memcpy(reference,x->value,(size_t)STRESS_N*sizeof(LIS_SCALAR));
                else if(memcmp(reference,x->value,
                               (size_t)STRESS_N*sizeof(LIS_SCALAR))!=0)
                {
                    fprintf(stderr,
                        "MFAMG_OMP_STRESS nondeterministic mode=%d "
                        "threads=%d repeat=%d\n",mode,threads[t],repeat);
                    goto done;
                }
            }
            if(omp_get_max_threads()!=threads[t])
            {
                fprintf(stderr,"MFAMG_OMP_STRESS thread ICV mutated\n");
                goto done;
            }
            apply_seconds=omp_get_wtime()-apply_started;
            /* Timings are informational; shared CI runners are not
             * stable enough to enforce a speedup threshold. */
            printf("MFAMG_OMP_STRESS mode=%d threads=%d relres=%.9e "
                   "setup_s=%.6f apply3_s=%.6f PASS\n",
                   mode,threads[t],(double)rel,
                   setup_seconds,apply_seconds);
            solver->precon=NULL;
            lis_precon_destroy(precon);precon=NULL;
            lis_solver_destroy(solver);solver=NULL;
        }
    }
    failed=0;
done:
    if(failed || err) fprintf(stderr,"MFAMG_OMP_STRESS stage=%s rc=%ld mode=%d threads=%d\n",stage,(long)err,mode,(mode>=0 && t>=0 && t<4)?threads[t]:-1);
    if(solver) solver->precon=NULL;
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    if(Ax) lis_vector_destroy(Ax);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(basis) lis_vector_destroy(basis);
    if(A) lis_matrix_destroy(A);
    free(reference);
    omp_set_num_threads(original_threads);
    omp_set_dynamic(original_dynamic);
    return failed;
#endif
}


/* Exercise both the short terminal block and scalar fallback under real
 * OpenMP block factorization. The first 4x4 block has rank one but nonzero
 * diagonal, while the remaining rows have an irregular CSR stencil. */
static int mfamg_omp_block_edges(void)
{
#ifndef _OPENMP
    return 0;
#else
    enum { EDGE_N=4099 };
    const int threads[4]={1,2,4,8};
    const int block_sizes[2]={4,5};
    /* Edge signatures measured with source at upstream develop
     * 14516697c38286e002c51c66867a8b1de9ea565a,
     * OMP_DYNAMIC=FALSE OMP_NUM_THREADS=8
     * ./test/testmfamgvcycle --omp-stress
     * The reference is the weighted real-part sum of x[i]. */
    const double edge_signature[2][2]={
        {1.503476740167112e4,1.804837753128949e4},
        {1.592697674017924e4,1.875507539588175e4}
    };
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SCALAR *reference=NULL;
    LIS_INT err=0,i,j,active=0,fallback=0;
    LIS_REAL omega=0.0;
    int old_threads=omp_get_max_threads();
    int old_dynamic=omp_get_dynamic();
    int failed=1,block=-1,phase=-1,t=-1,repeat=0;
    double start=0.0,elapsed=0.0;
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    if(mfamg_force_omp_test_threshold("4096")!=0 ||
       lis_precon_mfamg_get_omp_min_rows_for_test()>EDGE_N)
    {
        fprintf(stderr,"MFAMG_OMP_EDGES parallel threshold not enabled\n");
        goto done;
    }

    reference=(LIS_SCALAR *)malloc((size_t)EDGE_N*sizeof(LIS_SCALAR));
    if(!reference) goto done;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,EDGE_N);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<EDGE_N;i++)
    {
        if(i<4)
        {
            for(j=0;j<4;j++)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,j,
                                         (LIS_SCALAR)1.0,A);
                if(err) goto done;
            }
        }
        else
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
                                     (LIS_SCALAR)2.5,A);
            if(err) goto done;
            if(i>4)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                         (LIS_SCALAR)-0.7,A);
                if(err) goto done;
            }
            if(i+1<EDGE_N)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                         (LIS_SCALAR)-0.75,A);
                if(err) goto done;
            }
            if((i%11)==0 && i+2<EDGE_N)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i+2,
                                         (LIS_SCALAR)-0.08,A);
                if(err) goto done;
            }
            if((i%17)==0 && i>=6)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i-2,
                                         (LIS_SCALAR)-0.05,A);
                if(err) goto done;
            }
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x))
        goto done;
    for(i=0;i<EDGE_N;i++)
        b->value[i]=(LIS_SCALAR)(0.25+0.125*(i%7));

    for(block=0;block<2;block++)
    for(phase=0;phase<2;phase++)
    for(t=0;t<4;t++)
    {
        if(omp_get_thread_limit()<threads[t])
        {
            printf("MFAMG_OMP_EDGES block=%d phase=%d threads=%d "
                   "SKIP: OMP_THREAD_LIMIT=%d\n",
                   block_sizes[block],phase,threads[t],
                   omp_get_thread_limit());
            continue;
        }
        omp_set_num_threads(threads[t]);
        start=omp_get_wtime();
        for(repeat=0;repeat<3;repeat++)
        {
            for(i=0;i<EDGE_N;i++)
                x->value[i]=phase ? (LIS_SCALAR)(0.125*(i%3)) : 0.0;
            err=lis_precon_mfamg_apply_block_smoother_for_test(
                A,block_sizes[block],b,x,phase ? LIS_FALSE : LIS_TRUE,
                &active,&fallback,&omega);
            if(err) goto done;
            if(fallback<1 || active<1 ||
               active+fallback!=(EDGE_N+block_sizes[block]-1)/
                                 block_sizes[block] ||
               !(omega>0.0 && omega<=1.0))
            {
                fprintf(stderr,"MFAMG_OMP_EDGES invalid block status "
                               "block=%d threads=%d active=%ld fallback=%ld\n",
                        block_sizes[block],threads[t],(long)active,
                        (long)fallback);
                goto done;
            }
            for(i=0;i<EDGE_N;i++)
            {
#ifdef _COMPLEX
                if(!isfinite((double)cabs(x->value[i])))
#else
                if(!isfinite((double)x->value[i]))
#endif
                {
                    fprintf(stderr,"MFAMG_OMP_EDGES nonfinite x[%ld] "
                                   "block=%d phase=%d threads=%d\n",
                            (long)i,block_sizes[block],phase,threads[t]);
                    goto done;
                }
            }
            if(t==0 && repeat==0)
            {
                double signature=0.0;
                for(i=0;i<EDGE_N;i++)
                {
#ifdef _COMPLEX
                    signature+=(double)(1+i%23)*creal(x->value[i]);
#else
                    signature+=(double)(1+i%23)*(double)x->value[i];
#endif
                }
                if(fabs(signature-edge_signature[block][phase])>1.0e-6)
                {
                    fprintf(stderr,"MFAMG_OMP_EDGES baseline mismatch "
                                   "block=%d phase=%d actual=%.15e expected=%.15e\n",
                            block_sizes[block],phase,signature,
                            edge_signature[block][phase]);
                    goto done;
                }
                memcpy(reference,x->value,
                       (size_t)EDGE_N*sizeof(LIS_SCALAR));
            }
            else if(memcmp(reference,x->value,
                           (size_t)EDGE_N*sizeof(LIS_SCALAR))!=0)
            {
                fprintf(stderr,"MFAMG_OMP_EDGES differs block=%d "
                               "phase=%d threads=%d repeat=%d\n",
                        block_sizes[block],phase,threads[t],repeat);
                goto done;
            }
        }
        elapsed=omp_get_wtime()-start;
        if(omp_get_max_threads()!=threads[t])
        {
            fprintf(stderr,"MFAMG_OMP_EDGES OpenMP ICV changed "
                           "block=%d phase=%d expected=%d actual=%d\n",
                    block_sizes[block],phase,threads[t],
                    omp_get_max_threads());
            goto done;
        }
        printf("MFAMG_OMP_EDGES block=%d phase=%d threads=%d "
               "active=%ld fallback=%ld seconds=%.6f PASS\n",
               block_sizes[block],phase,threads[t],
               (long)active,(long)fallback,elapsed);
    }
    failed=0;
done:
    if(failed || err)
        fprintf(stderr,"MFAMG_OMP_EDGES failed=%d error=%ld "
                       "block=%d phase=%d threads=%d\n",
                failed,(long)err,
                (block>=0 && block<2)?block_sizes[block]:-1,
                phase,(t>=0 && t<4)?threads[t]:-1);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    free(reference);
    omp_set_num_threads(old_threads);
    omp_set_dynamic(old_dynamic);
    return failed;
#endif
}



/* Exercise the rank-adaptive multivector prolong_row_ptr branch with
 * two near-nullspace modes and genuine OpenMP loops at 24 DOF. The
 * two-mode setup mirrors the pre-OpenMP serial MFAMG regression.
 * A runtime row threshold of one activates the small hierarchy without
 * changing the production OpenMP threshold or the library build flags. */
static int mfamg_omp_multivector(void)
{
#ifndef _OPENMP
    return 0;
#else
    enum { MV_N=24 };
    const int thread_counts[4]={1,2,4,8};
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,ax=NULL,z0=NULL,z1=NULL,near[2];
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_SCALAR reference[MV_N];
    LIS_INT i,err=0,mvlevels=0;
    LIS_REAL bn=0.0,rn=0.0;
    int failed=1,t=-1,repeat;
    int old_threads=omp_get_max_threads();
    int old_dynamic=omp_get_dynamic();

    omp_set_dynamic(0);
    omp_set_num_threads(1);
    if(mfamg_force_omp_test_threshold("1")!=0 ||
       lis_precon_mfamg_get_omp_min_rows_for_test()!=1)
    {
        fprintf(stderr,"MFAMG_OMP_MULTIVECTOR cannot force low OMP threshold\n");
        goto done;
    }
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,MV_N);
    if(err) goto done;
    if(A->nprocs!=1) {failed=0;goto done;}
    for(i=0;i<MV_N;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)2.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
        if(i+1<MV_N)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(lis_vector_duplicate(A,&z0) || lis_vector_duplicate(A,&z1) ||
       lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&ax))
        goto done;
    for(i=0;i<MV_N;i++)
    {
        z0->value[i]=(LIS_SCALAR)1.0;
        z1->value[i]=(LIS_SCALAR)(i<10 ? (i%2) : 0);
        b->value[i]=(LIS_SCALAR)(1.0+0.2*(i%5));
    }
    near[0]=z0;near[1]=z1;
    for(t=0;t<4;t++)
    {
        if(omp_get_thread_limit()<thread_counts[t])
        {
            printf("MFAMG_OMP_MULTIVECTOR threads=%d SKIP: OMP_THREAD_LIMIT\n",
                   thread_counts[t]);
            continue;
        }
        omp_set_num_threads(thread_counts[t]);
        err=lis_solver_create(&solver);
        if(err) goto done;
        err=lis_solver_set_option(
            "-i gmres -p mfamg -scale none -print none",solver);
        if(err) goto done;
        err=lis_solver_set_near_nullspace(solver,2,near);
        if(err) goto done;
        err=lis_solver_set_matrix(A,solver);
        if(err) goto done;
        err=lis_precon_create(solver,&pc);
        if(err) {pc=NULL;goto done;}
        solver->precon=pc;
        err=lis_precon_mfamg_get_multivector_levels_for_test(pc,&mvlevels);
        if(err || mvlevels<1)
        {
            fprintf(stderr,"MFAMG_OMP_MULTIVECTOR no rank-adaptive "
                           "prolongation threads=%d levels=%ld\n",
                    thread_counts[t],(long)mvlevels);
            goto done;
        }
        for(repeat=0;repeat<3;repeat++)
        {
            err=lis_psolve(solver,b,x);
            if(err) goto done;
            err=lis_matvec(A,x,ax);
            if(err) goto done;
            bn=0.0;rn=0.0;
            for(i=0;i<MV_N;i++)
            {
#ifdef _COMPLEX
                LIS_REAL bb=(LIS_REAL)cabs(b->value[i]);
                LIS_REAL rr=(LIS_REAL)cabs(ax->value[i]-b->value[i]);
#else
                LIS_REAL bb=(LIS_REAL)fabs(b->value[i]);
                LIS_REAL rr=(LIS_REAL)fabs(ax->value[i]-b->value[i]);
#endif
                bn+=bb*bb; rn+=rr*rr;
            }
            /* A preconditioner correction is not required to contract
             * the original residual for every right-hand side. This
             * check rejects non-finite/catastrophic results; deterministic
             * repeated and cross-thread equivalence is checked below. */
            if(!(bn>0.0) || !isfinite((double)rn) || !(rn/bn<25.0))
            {
                fprintf(stderr,"MFAMG_OMP_MULTIVECTOR residual invalid "
                               "threads=%d rel=%.8e\n",
                        thread_counts[t],(double)sqrt(rn/bn));
                goto done;
            }
            if(t==0 && repeat==0)
                memcpy(reference,x->value,sizeof(reference));
            else if(memcmp(reference,x->value,sizeof(reference))!=0)
            {
                fprintf(stderr,"MFAMG_OMP_MULTIVECTOR correction mismatch "
                               "threads=%d repeat=%d\n",
                        thread_counts[t],repeat);
                goto done;
            }
        }
        printf("MFAMG_OMP_MULTIVECTOR threads=%d levels=%ld "
               "relres=%.9e PASS\n",
               thread_counts[t],(long)mvlevels,(double)sqrt(rn/bn));
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
    failed=0;
done:
    if(failed || err)
        fprintf(stderr,"MFAMG_OMP_MULTIVECTOR failed=%d rc=%ld "
                       "thread_index=%d\n",failed,(long)err,t);
    if(solver) solver->precon=NULL;
    if(pc) lis_precon_destroy(pc);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(z1) lis_vector_destroy(z1);
    if(z0) lis_vector_destroy(z0);
    if(A) lis_matrix_destroy(A);
    omp_set_num_threads(old_threads);
    omp_set_dynamic(old_dynamic);
    return failed;
#endif
}

/*
 * Optional, non-gating 3D Poisson benchmark: one million equations and
 * approximately seven million CSR entries. This isolates the *combined*
 * block-Jacobi construction plus one apply via the existing production
 * regression hook, not the time of the entire multilevel V-cycle.
 *
 * Run explicitly:
 * OMP_NUM_THREADS=8 ./test/testmfamgvcycle --omp-benchmark
 *
 * Shared CI runners can be oversubscribed; interpret timings as samples.
 */
static int mfamg_omp_large_benchmark(void)
{
#ifndef _OPENMP
    puts("MFAMG_OMP_BENCH SKIP: built without OpenMP");
    return 0;
#else
    enum { GRID=100, BENCH_N=GRID*GRID*GRID };
    const int counts[4]={1,2,4,8};
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SCALAR *reference=NULL;
    LIS_INT i,err=0,active=0,fallback=0;
    LIS_REAL omega=0.0;
    int t,repeat,failed=1;
    int original_threads=omp_get_max_threads();
    int original_dynamic=omp_get_dynamic();
    double start,elapsed,baseline=0.0,speedup;
    omp_set_dynamic(0);
    omp_set_num_threads(1);

    reference=(LIS_SCALAR *)malloc((size_t)BENCH_N*sizeof(LIS_SCALAR));
    if(!reference) goto done;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,BENCH_N);
    if(err) goto done;
    if(A->nprocs!=1)
    {
        puts("MFAMG_OMP_BENCH SKIP: only single-rank CSR supported");
        failed=0;
        goto done;
    }
    for(i=0;i<BENCH_N;i++)
    {
        LIS_INT ix=i%GRID,iy=(i/GRID)%GRID,iz=i/(GRID*GRID);
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)6.001,A);
        if(err) goto done;
        if(ix>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
        if(ix+1<GRID)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
        if(iy>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-GRID,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
        if(iy+1<GRID)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+GRID,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
        if(iz>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-GRID*GRID,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
        if(iz+1<GRID)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+GRID*GRID,(LIS_SCALAR)-1,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x))
        goto done;
    for(i=0;i<BENCH_N;i++)
        b->value[i]=(LIS_SCALAR)(1.0+(i%11)*0.03125);

    for(t=0;t<4;t++)
    {
        double mean=0.0;
        omp_set_num_threads(counts[t]);
        for(repeat=0;repeat<2;repeat++)
        {
            start=omp_get_wtime();
            err=lis_precon_mfamg_apply_block_smoother_for_test(
                A,4,b,x,LIS_TRUE,&active,&fallback,&omega);
            elapsed=omp_get_wtime()-start;
            if(err) goto done;
            if(active+fallback!=(BENCH_N+3)/4)
            {
                fprintf(stderr,"MFAMG_OMP_BENCH bad block count\n");
                goto done;
            }
            mean+=elapsed;
            if(t==0 && repeat==0)
                memcpy(reference,x->value,(size_t)BENCH_N*sizeof(LIS_SCALAR));
            else if(memcmp(reference,x->value,
                           (size_t)BENCH_N*sizeof(LIS_SCALAR))!=0)
            {
                fprintf(stderr,"MFAMG_OMP_BENCH non-deterministic correction "
                               "threads=%d repeat=%d\n",counts[t],repeat);
                goto done;
            }
        }
        mean/=2;
        if(t==0) baseline=mean;
        speedup=mean>0 ? baseline/mean : 0.0;
        printf("MFAMG_OMP_BENCH grid=100x100x100 n=%ld threads=%d "
               "avg_setup_apply_s=%.6f relative_speedup=%.3f "
               "active=%ld fallback=%ld PASS\n",
               (long)BENCH_N,counts[t],mean,speedup,
               (long)active,(long)fallback);
        fflush(stdout);
    }
    failed=0;
done:
    if(failed || err)
        fprintf(stderr,"MFAMG_OMP_BENCH failed=%d error=%ld\n",
                failed,(long)err);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    free(reference);
    omp_set_num_threads(original_threads);
    omp_set_dynamic(original_dynamic);
    return failed;
#endif
}


/*
 * A collection of disconnected, positive-definite 16-row chains with 12
 * locally independent near-nullspace vectors stalls the multivector
 * hierarchy: the first transfer is 576 -> 432, while subsequent local
 * QR transfers cannot reduce the 12-dimensional block spaces further.
 * Exercise the opt-in sparse terminal without an O(n^2) coarse factor,
 * and verify deterministic finite corrections under OpenMP.
 */
static int mfamg_sparse_terminal_regression(LIS_REAL diagonal)
{
    enum { CHAIN=16, BLOCKS=36, MODES=12, NTEST=CHAIN*BLOCKS };
    LIS_MATRIX A=NULL;
    LIS_VECTOR near[MODES]={NULL}, b=NULL, x=NULL, Ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_SCALAR *reference=NULL, *baseline8=NULL;
    LIS_REAL normb,normr,diff,relative8=0.0;
    LIS_INT i,m,err=0,nlevels=0,coarse_n=0;
    LIS_INT sizes[MAX_LEVELS],owned[MAX_LEVELS];
    int run,repeat,failed=1,pass_case;
    const char *stage="setup";
#ifdef _OPENMP
    const int threads[2]={1,8};
    int original_threads=omp_get_max_threads();
    int original_dynamic=omp_get_dynamic();
    omp_set_dynamic(0);
#else
    const int threads[2]={1,1};
#endif
    if(mfamg_force_omp_test_threshold("1")!=0) goto done;
#ifdef _WIN32
    if(_putenv_s("LIS_MFAMG_SPARSE_TERMINAL","1")!=0) goto done;
#else
    if(setenv("LIS_MFAMG_SPARSE_TERMINAL","1",1)!=0) goto done;
#endif
    reference=(LIS_SCALAR *)malloc((size_t)NTEST*sizeof(LIS_SCALAR));
    baseline8=(LIS_SCALAR *)malloc((size_t)NTEST*sizeof(LIS_SCALAR));
    if(reference==NULL || baseline8==NULL) goto done;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,NTEST);
    if(err) goto done;
    if(A->nprocs!=1) {failed=0;goto done;}
    for(i=0;i<NTEST;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)diagonal,A);
        if(err) goto done;
        if(i%CHAIN>0)
        {
            err=lis_matrix_set_value(
                LIS_INS_VALUE,i,i-1,(LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
        if(i%CHAIN+1<CHAIN)
        {
            err=lis_matrix_set_value(
                LIS_INS_VALUE,i,i+1,(LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&Ax)) goto done;
    for(m=0;m<MODES;m++)
    {
        err=lis_vector_duplicate(A,&near[m]);
        if(err) goto done;
    }
    for(i=0;i<NTEST;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+(i%7)*0.0625);
        for(m=0;m<MODES;m++)
            near[m]->value[i]=
                (LIS_SCALAR)((i%CHAIN)==m ? 1.0 : 0.0);
    }

    /* The unset environment must be bitwise compatible with 8 passes.
     * 16 passes are opt-in and should improve this manufactured SPD case.
     * Compare repeated applications across serial and OpenMP workers. */
    for(pass_case=0;pass_case<2;pass_case++)
    {
#ifdef _WIN32
        if(_putenv_s("LIS_MFAMG_SPARSE_TERMINAL_PASSES",
                     pass_case ? "16" : "")!=0) goto done;
#else
        if(pass_case ?
           setenv("LIS_MFAMG_SPARSE_TERMINAL_PASSES","16",1) :
           unsetenv("LIS_MFAMG_SPARSE_TERMINAL_PASSES")) goto done;
#endif
    for(run=0;run<2;run++)
    {
#ifdef _OPENMP
        if(omp_get_thread_limit()<threads[run]) continue;
        omp_set_num_threads(threads[run]);
#endif
        stage="solver-create";
        err=lis_solver_create(&solver);
        if(err) goto done;
        stage="solver-options";
        err=lis_solver_set_option(
            "-i fgmres -p mfamg -scale none -print none "
            "-mfamg_block_size 4 -mfamg_fine_agg_factor 2",
            solver);
        if(err) goto done;
        stage="near-nullspace";
        err=lis_solver_set_near_nullspace(solver,MODES,near);
        if(err) goto done;
        err=lis_solver_set_matrix(A,solver);
        if(err) goto done;
        stage="preconditioner-create";
        err=lis_precon_create(solver,&pc);
        if(err) {pc=NULL;goto done;}
        solver->precon=pc;
        err=lis_precon_mfamg_get_hierarchy_info(
            pc,MAX_LEVELS,&nlevels,&coarse_n,sizes,owned);
        if(err || nlevels<2 || coarse_n<=8*MODES || coarse_n>=NTEST)
        {
            fprintf(stderr,
                "MFAMG_SPARSE_TERMINAL unexpected terminal n=%d levels=%d\n",
                (int)coarse_n,(int)nlevels);
            goto done;
        }
        for(repeat=0;repeat<2;repeat++)
        {
            stage="psolve";
            err=lis_psolve(solver,b,x);
            if(err) goto done;
            err=lis_matvec(A,x,Ax);
            if(err) goto done;
            normb=0.0;
            normr=0.0;
            for(i=0;i<NTEST;i++)
            {
#ifdef _COMPLEX
                diff=(LIS_REAL)cabs(Ax->value[i]-b->value[i]);
                normb+=(LIS_REAL)pow(cabs(b->value[i]),2.0);
#else
                diff=(LIS_REAL)fabs(Ax->value[i]-b->value[i]);
                normb+=(LIS_REAL)(b->value[i]*b->value[i]);
#endif
                normr+=diff*diff;
            }
            if(!isfinite((double)normr) || !isfinite((double)normb) ||
               normb<=0.0 || normr>=normb)
            {
                fprintf(stderr,
                    "MFAMG_SPARSE_TERMINAL defect not reduced: %.9e\n",
                    (double)(normr/normb));
                goto done;
            }
            if(run==0 && repeat==0)
            {
                memcpy(reference,x->value,
                       (size_t)NTEST*sizeof(LIS_SCALAR));
                if(pass_case==0)
                {
                    memcpy(baseline8,x->value,
                           (size_t)NTEST*sizeof(LIS_SCALAR));
                    relative8=normr/normb;
                }
                else if(memcmp(baseline8,x->value,
                               (size_t)NTEST*sizeof(LIS_SCALAR))==0 ||
                        normr/normb >= relative8 *
                           (diagonal<3.0 ? (LIS_REAL)0.99 : (LIS_REAL)1.0))
                {
                    fprintf(stderr,
                        "MFAMG_SPARSE_TERMINAL_PASSES insufficient improvement "
                        "diagonal=%.3f 8pass=%.9e 16pass=%.9e\n",
                        (double)diagonal,(double)relative8,
                        (double)(normr/normb));
                    goto done;
                }
            }
            else if(memcmp(reference,x->value,
                           (size_t)NTEST*sizeof(LIS_SCALAR))!=0)
            {
                fprintf(stderr,
                    "MFAMG_SPARSE_TERMINAL correction differs "
                    "passes=%d run=%d repeat=%d\n",
                    pass_case ? 16 : 8,run,repeat);
                goto done;
            }
        }
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
        printf("MFAMG_SPARSE_TERMINAL_PASSES passes=%d diagonal=%.3f "
               "relative_defect=%.9e PASS\n",
               pass_case ? 16 : 8,(double)diagonal,(double)(normr/normb));
    }
    if(diagonal<3.0)
        printf("MFAMG_SPARSE_TERMINAL_SLOW diagonal=%.3f "
               "relative_improvement=%.6f PASS\n",
               (double)diagonal,
               (double)(1.0-(normr/normb)/relative8));
    /* Reject invalid work budgets rather than silently ignoring typos. */
    for(pass_case=0;pass_case<5;pass_case++)
    {
        const char *invalid[5]={"0","33","invalid","16abc","-1"};
#ifdef _WIN32
        if(_putenv_s("LIS_MFAMG_SPARSE_TERMINAL_PASSES",
                     invalid[pass_case])!=0) goto done;
#else
        if(setenv("LIS_MFAMG_SPARSE_TERMINAL_PASSES",
                  invalid[pass_case],1)!=0) goto done;
#endif
        stage="invalid-pass-solver-create";
        err=lis_solver_create(&solver);
        if(err) goto done;
        err=lis_solver_set_option(
            "-i fgmres -p mfamg -scale none -print none "
            "-mfamg_block_size 4 -mfamg_fine_agg_factor 2",solver);
        if(err) goto done;
        err=lis_solver_set_near_nullspace(solver,MODES,near);
        if(err) goto done;
        err=lis_solver_set_matrix(A,solver);
        if(err) goto done;
        stage="invalid-pass-preconditioner";
        err=lis_precon_create(solver,&pc);
        if(err!=LIS_ERR_ILL_ARG)
        {
            fprintf(stderr,
                "MFAMG_SPARSE_TERMINAL invalid passes=%s got rc=%d\n",
                invalid[pass_case],(int)err);
            goto done;
        }
        pc=NULL; /* Failed factory already destroys the partial preconditioner. */
        lis_solver_destroy(solver);solver=NULL;
        err=LIS_SUCCESS;
    }
    puts("MFAMG_SPARSE_TERMINAL_PASSES invalid range PASS");
    puts("MFAMG_SPARSE_TERMINAL stalled multivector hierarchy PASS");
    failed=0;
done:
    if(failed)
        fprintf(stderr,"MFAMG_SPARSE_TERMINAL stage=%s rc=%d\n",
                stage,(int)err);
    if(solver) solver->precon=NULL;
    if(pc) lis_precon_destroy(pc);
    if(solver) lis_solver_destroy(solver);
    if(Ax) lis_vector_destroy(Ax);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    for(m=0;m<MODES;m++)
        if(near[m]) lis_vector_destroy(near[m]);
    if(A) lis_matrix_destroy(A);
    free(reference);
    free(baseline8);
#ifdef _WIN32
    (void)_putenv_s("LIS_MFAMG_SPARSE_TERMINAL_PASSES","");
    (void)_putenv_s("LIS_MFAMG_SPARSE_TERMINAL","");
#else
    (void)unsetenv("LIS_MFAMG_SPARSE_TERMINAL_PASSES");
    (void)unsetenv("LIS_MFAMG_SPARSE_TERMINAL");
#endif
    (void)mfamg_force_omp_test_threshold("16384");
#ifdef _OPENMP
    omp_set_num_threads(original_threads);
    omp_set_dynamic(original_dynamic);
#endif
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR z=NULL,b=NULL,x=NULL;
    LIS_INT i,err=0;
    LIS_SCALAR x_default[N],x_one[N],x_two[N],x_three[N];
    LIS_REAL d_default=0.0,d_one=0.0,d_two=0.0,d_three=0.0;
    LIS_REAL difference_12=0.0,difference_23=0.0;
    int failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    if(argc==2 && strcmp(argv[1],"--omp-benchmark")==0)
    {
        failed=mfamg_omp_large_benchmark();
        lis_finalize();
        return failed ? 1 : 0;
    }
    if(argc==2 && strcmp(argv[1],"--omp-stress")==0)
    {
        failed=mfamg_omp_stress();
        if(!failed) failed=mfamg_omp_block_edges();
        if(!failed) failed=mfamg_omp_multivector();
        if(!failed) failed=mfamg_sparse_terminal_regression((LIS_REAL)4.0);
        if(!failed) failed=mfamg_sparse_terminal_regression((LIS_REAL)2.25);
        lis_finalize();
        return failed ? 1 : 0;
    }
#ifdef USE_MPI
    {
        int ranks=1;
        MPI_Comm_size(LIS_COMM_WORLD,&ranks);
        if(ranks!=1) {lis_finalize();return 0;}
    }
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,N);
    if(err) goto done;
    for(i=0;i<N;i++) {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A);
        if(err) goto done;
        if(i>0) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
            if(err) goto done;
        }
        if(i+1<N) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_vector_duplicate(A,&z);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    for(i=0;i<N;i++) {
        z->value[i]=1.0;
        b->value[i]=1.0+((LIS_REAL)(i%5))/5.0;
    }
    if(run_case(A,z,b,x,1,0,x_default,&d_default)!=0 ||
       run_case(A,z,b,x,1,1,x_one,&d_one)!=0 ||
       run_case(A,z,b,x,2,1,x_two,&d_two)!=0 ||
       run_case(A,z,b,x,3,1,x_three,&d_three)!=0) goto done;

    /* The default and the explicit one-pass setting must agree bitwise. */
    if(memcmp(x_default,x_one,sizeof(x_default))!=0 ||
       d_default!=d_one)
    {
        fprintf(stderr,"MFAMG VCYCLE default differs from explicit steps=1\n");
        goto done;
    }
    for(i=0;i<N;i++)
    {
#ifdef _COMPLEX
        difference_12 += (LIS_REAL)pow(cabs(x_two[i]-x_one[i]),2.0);
        difference_23 += (LIS_REAL)pow(cabs(x_three[i]-x_two[i]),2.0);
#else
        difference_12 += (LIS_REAL)pow(fabs(x_two[i]-x_one[i]),2.0);
        difference_23 += (LIS_REAL)pow(fabs(x_three[i]-x_two[i]),2.0);
#endif
    }
    /* Verify the option actually changes the correction on this matrix. */
    if(!isfinite((double)difference_12) ||
       !isfinite((double)difference_23) ||
       difference_12<=1.0e-18 || difference_23<=1.0e-18)
    {
        fprintf(stderr,"MFAMG VCYCLE additional passes had no measurable effect\n");
        goto done;
    }
    printf("MFAMG VCYCLE default=%.9e steps2=%.9e steps3=%.9e PASS\n",
           (double)d_default,(double)d_two,(double)d_three);
    failed=0;
#ifdef _OPENMP
    /* Included in make check, not only the custom CI smoke command. */
    if(mfamg_omp_stress()!=0 || mfamg_omp_block_edges()!=0 ||
       mfamg_omp_multivector()!=0)
        failed=1;
#endif
    if(!failed && (mfamg_sparse_terminal_regression((LIS_REAL)4.0)!=0 ||
                   mfamg_sparse_terminal_regression((LIS_REAL)2.25)!=0))
        failed=1;
done:
    if(err) fprintf(stderr,"MFAMG VCYCLE setup rc=%d\n",(int)err);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(z) lis_vector_destroy(z);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
