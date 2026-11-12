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
#include <string.h>
#include "lislib.h"

#define N 96
#define MAX_LEVELS 32

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
done:
    if(err) fprintf(stderr,"MFAMG VCYCLE setup rc=%d\n",(int)err);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(z) lis_vector_destroy(z);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
