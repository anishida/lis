#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

/* Regression: BiCGSTAB(l) once added the nonzero initial x twice:
 * for A=4I, b[0]=1, x0[0]=0.15 it returned x[0]=0.40 instead of 0.25.
 *
 * Require LIS_SUCCESS and the exact physical iterate. The solver's
 * true-residual certification must reject an incorrectly doubled x0.
 */
typedef struct { LIS_INT n,calls; int varying; } OPERATOR;
static LIS_INT ax(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    OPERATOR *o=(OPERATOR *)ptr;
    LIS_INT i;
    ++o->calls;
    for(i=0;i<o->n;i++) y[i]=(LIS_SCALAR)(4.0+(o->varying?0.125*(double)i:0.0))*x[i];
    return LIS_SUCCESS;
}

static int run(LIS_INT ell,int init_zeros,int varying,int stop)
{
    OPERATOR op;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=LIS_SUCCESS,status=LIS_SUCCESS,i;
    LIS_REAL residual=0.0;
    double initial2=0.0,physical2=0.0;
    int failed=1;
    op.n=8;op.calls=0;op.varying=varying;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,op.n);if(err)goto done;
    err=lis_matrix_set_user(A,&op,ax,ax);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<op.n;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
        x->value[i]=(LIS_SCALAR)0.15;
        if(stop==1)
        {
            x->value[i]=(LIS_SCALAR)0.25;
            b->value[i]=(LIS_SCALAR)(4.0+(varying?0.125*(double)i:0.0))*x->value[i];
        }
        {
            double r=(double)b->value[i]-(init_zeros?0.0:
                     (4.0+(varying?0.125*(double)i:0.0))*(double)x->value[i]);
            initial2+=r*r;
        }
    }
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(
        "-i bicgstabl -p none -scale none -maxiter 120 -tol 1e-9 -print none",
        solver);
    if(err)goto done;
    if(stop==2)
    {
        err=lis_solver_set_option("-maxiter 1 -tol 1e-14",solver);
        if(err)goto done;
    }
    solver->options[LIS_OPTIONS_ELL]=ell;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=init_zeros?LIS_TRUE:LIS_FALSE;
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status) || lis_solver_get_residualnorm(solver,&residual))goto done;
    if(stop==2)
    {
        if(err!=LIS_SUCCESS || status!=LIS_MAXITER)goto done;
        for(i=0;i<op.n;i++)
        {
            double r=(double)b->value[i]-
              (4.0+(varying?0.125*(double)i:0.0))*(double)x->value[i];
            physical2+=r*r;
        }
        if(!(fabs(sqrt(physical2/initial2)-(double)residual)<=1e-9))
        {
            fprintf(stderr,"MAXITER physical residual mismatch: %.17g vs %.17g\n",
                    sqrt(physical2/initial2),(double)residual);
            goto done;
        }
        failed=0;
        goto done;
    }
    /* The solver's true-residual post-certification must agree with
     * the numerical solution once x0 is added exactly once. */
    if(err!=LIS_SUCCESS || status!=LIS_SUCCESS)
    {
        fprintf(stderr,"BiCGSTAB(l) api=%d status=%d ell=%d init_zeros=%d\n",
                (int)err,(int)status,(int)ell,init_zeros);
        goto done;
    }
    for(i=0;i<op.n;i++)
    {
        double exact=(double)b->value[i]/(4.0+(varying?0.125*(double)i:0.0));
        if(!(fabs((double)x->value[i]-exact)<=1e-8))
        {
            fprintf(stderr,"BiCGSTAB(l) x[%d]=%.17g exact=%.17g ell=%d init_zeros=%d status=%d\n",
                    (int)i,(double)x->value[i],exact,
                    (int)ell,init_zeros,(int)status);
            goto done;
        }
    }
    if(stop==1 && op.calls!=2)goto done;
    if(stop==0 && varying && op.calls<=ell+2)
    {
        fprintf(stderr,"BiCGSTAB(l) nontrivial diagonal did not exercise multiple applications: calls=%d\n",(int)op.calls);
        goto done;
    }
    failed=0;
done:
    if(failed)fprintf(stderr,"BICGSTABL_NONZERO_INIT FAILED ell=%d zero=%d\n",
                      (int)ell,init_zeros);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}
int main(int argc,char **argv)
{
#ifdef USE_MPI
    (void)argc;(void)argv;
    puts("BICGSTABL_NONZERO_INIT SKIP: serial callback test");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(1,0,0,0);
    failed|=run(1,0,1,0);
    failed|=run(2,0,0,0);
    failed|=run(2,0,1,0);
    failed|=run(1,1,0,0);
    failed|=run(1,1,1,0);
    failed|=run(2,1,0,0);
    failed|=run(2,1,1,0);
    failed|=run(1,0,1,1);
    failed|=run(2,0,1,1);
    failed|=run(1,0,1,2);
    failed|=run(2,0,1,2);
    lis_finalize();
    printf("BICGSTABL_NONZERO_INIT %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
