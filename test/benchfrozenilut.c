#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include "lis_solver.h"
#include "lis_matrix_operator.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

/*
 * Reproducible 2D diffusion-stencil proxy for changing sparse FE systems.
 * Fixed topology, diagonal coefficient changes, zero initial guess per step.
 * This is a standalone performance experiment, not a CCX benchmark.
 * Intentionally excluded from TESTS (make check).
 */
static double seconds_now(void)
{
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)t.tv_sec + 1.0e-6*(double)t.tv_usec;
}

static LIS_INT create_diffusion(int side,LIS_MATRIX *out)
{
    LIS_MATRIX A=NULL;
    LIS_INT i,n=(LIS_INT)(side*side),err;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) return err;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto failed;
    for(i=0;i<n;i++)
    {
        LIS_INT row=i/side,col=i%side;
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.8,A);
        if(err) goto failed;
        if(row>0){err=lis_matrix_set_value(LIS_INS_VALUE,i,i-side,(LIS_SCALAR)-1.0,A);if(err)goto failed;}
        if(row<side-1){err=lis_matrix_set_value(LIS_INS_VALUE,i,i+side,(LIS_SCALAR)-1.0,A);if(err)goto failed;}
        if(col>0){err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,(LIS_SCALAR)-1.0,A);if(err)goto failed;}
        if(col<side-1){err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,(LIS_SCALAR)-1.0,A);if(err)goto failed;}
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto failed;
    err=lis_matrix_assemble(A);
    if(err) goto failed;
    *out=A;
    return LIS_SUCCESS;
failed:
    lis_matrix_destroy(A);
    return err;
}

static LIS_INT create_identity(int n,LIS_MATRIX *out)
{
    LIS_MATRIX A=NULL;
    LIS_INT i,err;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)return err;
    err=lis_matrix_set_size(A,0,n);if(err)goto fail;
    for(i=0;i<n;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,A);
        if(err)goto fail;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto fail;
    err=lis_matrix_assemble(A);if(err)goto fail;
    *out=A;return LIS_SUCCESS;
fail:
    lis_matrix_destroy(A);return err;
}

static int run_sequence(LIS_MATRIX C,int side,int steps,LIS_INT policy,
                        double *setup,double *solve,double *wall,
                        long *iterations,LIS_REAL *worst)
{
    LIS_SOLVER solver=NULL;
    LIS_VECTOR truth=NULL,b=NULL,x=NULL;
    LIS_INT i,step,err,status,iter,n=(LIS_INT)(side*side);
    LIS_PRECON precon=NULL;
    double t0,t1,tbegin;
    LIS_REAL norm;
    int bad=1;
    *setup=0;*solve=0;*wall=0;*iterations=0;*worst=0;
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(
        "-i gmres -restart 25 -p ilut -ilut_drop 1e-3 "
        "-ilut_rate 2.0 -tol 1e-10 -maxiter 400 "
        "-scale none -initx_zeros true",solver);
    if(err)goto done;
    err=lis_solver_set_sequence_policy(solver,LIS_SEQUENCE_POLICY_AUTO);
    if(err)goto done;
    err=lis_solver_set_sequence_precon_refresh_policy(solver,policy);
    if(err)goto done;
    err=lis_vector_duplicate(C,&truth);if(err)goto done;
    err=lis_vector_duplicate(C,&b);if(err)goto done;
    err=lis_vector_duplicate(C,&x);if(err)goto done;
    for(i=0;i<n;i++)
        truth->value[i]=(LIS_SCALAR)(0.9+0.003*(double)(i%side)+0.001*(double)(i/side));
    tbegin=seconds_now();
    for(step=0;step<steps;step++)
    {
        /* Vary only the matrix coefficients, never the graph. */
        LIS_REAL shift=(LIS_REAL)(0.18 + 0.12*sin(0.32*(double)step));
        err=lis_matrix_operator_set_coefficients(
            C,(LIS_SCALAR)1.0,(LIS_SCALAR)(-shift));
        if(err)goto done;
        err=lis_matvec(C,truth,b);if(err)goto done;
        err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto done;
        t0=seconds_now();
        err=lis_solver_sequence_prepare_precon(C,solver,&precon);
        t1=seconds_now();
        if(err||precon==NULL)goto done;
        *setup+=t1-t0;
        t0=seconds_now();
        err=lis_solve_sequence(C,b,x,solver);
        t1=seconds_now();
        if(err)goto done;
        *solve+=t1-t0;
        err=lis_solver_get_status(solver,&status);
        if(err||status!=LIS_SUCCESS)goto done;
        err=lis_solver_get_iter(solver,&iter);
        if(err)goto done;
        *iterations+=(long)iter;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,truth,x);
        if(err)goto done;
        err=lis_vector_nrm2(x,&norm);
        if(err || !isfinite((double)norm) || norm>(LIS_REAL)1e-6)goto done;
        if(norm>*worst)*worst=norm;
    }
    *wall=seconds_now()-tbegin;
    bad=0;
done:
    if(bad)fprintf(stderr,
        "benchmark failed: mode=%d step=%d err=%d (n=%d)\n",
        (int)policy,step,(int)err,(int)n);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(truth)lis_vector_destroy(truth);
    return bad;
}

int main(int argc,char **argv)
{
    int side=24,steps=24,failed=1;
    LIS_MATRIX A=NULL,B=NULL,C=NULL;
    LIS_INT err;
    double setup[2],solve[2],wall[2];
    long iterations[2];
    LIS_REAL worst[2];
    if(argc>1)side=atoi(argv[1]);
    if(argc>2)steps=atoi(argv[2]);
    if(side<10||side>96||steps<8||steps>200)
    {
        fprintf(stderr,"usage: %s [side 10..96] [steps 8..200]\n",argv[0]);
        return 2;
    }
    err=lis_initialize(&argc,&argv);if(err)return 1;
    err=create_diffusion(side,&A);if(err)goto done;
    err=create_identity(side*side,&B);if(err)goto done;
    err=lis_matrix_create_operator(
        (LIS_SCALAR)1.0,A,(LIS_SCALAR)-0.18,B,&C);
    if(err)goto done;
    if(run_sequence(C,side,steps,LIS_SEQUENCE_PRECON_REFRESH_SAFE,
                    &setup[0],&solve[0],&wall[0],&iterations[0],&worst[0]))
        goto done;
    if(run_sequence(C,side,steps,LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN,
                    &setup[1],&solve[1],&wall[1],&iterations[1],&worst[1]))
        goto done;
    printf("BENCH_FROZEN_ILUT side=%d n=%d steps=%d\n",side,side*side,steps);
    printf("policy,setup_seconds,solve_seconds,total_seconds,iterations,max_solution_error\n");
    printf("SAFE,%.6f,%.6f,%.6f,%ld,%.12g\n",
           setup[0],solve[0],wall[0],iterations[0],(double)worst[0]);
    printf("FROZEN,%.6f,%.6f,%.6f,%ld,%.12g\n",
           setup[1],solve[1],wall[1],iterations[1],(double)worst[1]);
    if(setup[1]>0)
        printf("SETUP_SPEEDUP=%.6f\n",setup[0]/setup[1]);
    if(wall[1]>0)
        printf("TOTAL_SPEEDUP=%.6f\n",wall[0]/wall[1]);
    printf("BENCH_FROZEN_ILUT PASS\n");
    failed=0;
done:
    if(C)lis_matrix_destroy(C);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}