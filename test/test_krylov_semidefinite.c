/* Regression: positive-semidefinite generalized metric B.
 * Independent analytic reference: diagonal generalized eigenpairs,
 * and 2-by-2 blocks with a massless constrained DOF.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 64
#define NMODE 4
#define NACTIVE 32

typedef struct {
    double a[NACTIVE];
    double c[NACTIVE];
    double d[NACTIVE];
    double b[NACTIVE];
    double sigma;
    int solve_count;
} blocks;

static int apply_a(void *ctx,int n,const double *x,double *y)
{
    blocks *s=(blocks*)ctx;
    int i;
    if(n!=N) return -1;
    for(i=0;i<NACTIVE;i++) {
        y[i]=s->a[i]*x[i]+s->c[i]*x[NACTIVE+i];
        y[NACTIVE+i]=s->c[i]*x[i]+s->d[i]*x[NACTIVE+i];
    }
    return 0;
}
static int apply_b(void *ctx,int n,const double *x,double *y)
{
    blocks *s=(blocks*)ctx;
    int i;
    if(n!=N) return -1;
    for(i=0;i<NACTIVE;i++) {
        y[i]=s->b[i]*x[i];
        y[NACTIVE+i]=0.;
    }
    return 0;
}
static int solve_shifted(void *ctx,int n,const double *x,double *y)
{
    blocks *s=(blocks*)ctx;
    int i;
    if(n!=N) return -1;
    ++s->solve_count;
    for(i=0;i<NACTIVE;i++) {
        double aa=s->a[i]-s->sigma*s->b[i];
        double dd=s->d[i],c=s->c[i];
        double det=aa*dd-c*c;
        if(det==0.) return -1;
        y[i]=(dd*x[i]-c*x[NACTIVE+i])/det;
        y[NACTIVE+i]=(aa*x[NACTIVE+i]-c*x[i])/det;
    }
    return 0;
}

static int run_psd(blocks *s,const char *label,unsigned seed,
                   const double *start,double expected[NMODE])
{
    lis_krylov_problem p;
    lis_krylov_options opt;
    lis_krylov_stats st;
    double vals[NMODE],vec[N*NMODE],res[NMODE];
    int i,j,status;
    memset(&p,0,sizeof(p));
    p.n=N;
    p.transform=LIS_KRYLOV_SHIFT_INVERT_B;
    p.sigma=s->sigma;
    p.ctx=s;
    p.apply_a=apply_a;
    p.apply_b=apply_b;
    p.solve=solve_shifted;
    lis_krylov_options_default(&opt);
    opt.nev=NMODE;
    opt.ncv=24;
    opt.maxiter=450;
    opt.tol=1e-9;
    opt.allow_semidefinite_metric=1;
    opt.initial_vector=start;
    opt.random_seed=seed;
    opt.inner_tol=1e-11;
    status=lis_krylov_symmetric(&p,&opt,vals,vec,res,&st);
    if(status!=0 || st.nconv!=NMODE || st.nreturned!=NMODE) {
        fprintf(stderr,"FAIL %s status=%d conv=%d returned=%d solve=%d restarts=%d\n",
                label,status,st.nconv,st.nreturned,st.solve_count,st.restarts);
        return 1;
    }
    for(i=0;i<NMODE;i++){
        double nrmb=0.;
        if(fabs(vals[i]-expected[i]) > 2e-8*fmax(1.,fabs(expected[i]))) {
            fprintf(stderr,"FAIL %s eigen[%d]=%.16g expected %.16g\n",
                    label,i,vals[i],expected[i]);
            return 1;
        }
        for(j=0;j<N;j++){
            double z=vec[(size_t)i*N+j];
            if(!isfinite(z)) return 1;
            if(j<NACTIVE)nrmb+=s->b[j]*z*z;
        }
        if(fabs(nrmb-1.)>2e-7 || res[i]>opt.tol) {
            fprintf(stderr,"FAIL %s norm %.6g res %.6g\n",label,nrmb,res[i]);
            return 1;
        }
        for(j=0;j<NMODE;j++){
            double overlap=0.;int k;
            for(k=0;k<NACTIVE;k++)
                overlap+=s->b[k]*vec[(size_t)i*N+k]*vec[(size_t)j*N+k];
            if(fabs(overlap-(i==j?1.:0.))>1e-6){
                fprintf(stderr,"FAIL %s B-orthogonality i=%d j=%d %.4g\n",
                        label,i,j,overlap);
                return 1;
            }
        }
    }
    printf("PASS %-19s dim=%d rank=%d iter=%d solve=%d maxres=%.3g\n",
           label,N,NACTIVE,st.iterations,st.solve_count,st.max_residual);
    return 0;
}
int main(void)
{
    blocks s;
    double expected[NMODE];
    double start[N]={0.};
    int i,err;
    lis_krylov_problem p;
    lis_krylov_options opt;
    lis_krylov_stats st;
    double vals[NMODE],vec[N*NMODE],res[NMODE];
    memset(&s,0,sizeof(s));
    s.sigma=-0.25;
    for(i=0;i<NACTIVE;i++){
        s.a[i]=1.0+0.05*i;
        s.d[i]=2.0;
        s.c[i]=0.25;
        s.b[i]=1.0;
    }
    for(i=0;i<NMODE;i++) expected[i]=s.a[i]-s.c[i]*s.c[i]/s.d[i];
    if(run_psd(&s,"coupled PSD B",1,NULL,expected))return 1;
    s.sigma = expected[0] + 1e-6;
    if(run_psd(&s,"near-pole PSD B",11,NULL,expected))return 1;
    s.sigma = -0.25;
    /* Initial guess lives entirely in ker B; a range-projected random
     * replacement is necessary for nonzero metric norm. */
    for(i=0;i<NACTIVE;i++)start[NACTIVE+i]=i+1.;
    if(run_psd(&s,"start in kernel B",3,start,expected))return 1;
    /* A singular stiffness matrix gives a rigid-body (lambda=0) mode.
     * A nonzero shift makes the shifted operator nonsingular. */
    s.a[0] = s.c[0]*s.c[0]/s.d[0];
    expected[0] = 0.0;
    if(run_psd(&s,"zero mode and PSD B",5,NULL,expected))return 1;
    s.a[0] = 1.0;
    expected[0] = s.a[0]-s.c[0]*s.c[0]/s.d[0];
    /* Arbitrary scaling of A, B and sigma by the same constant preserves
     * eigenvalues and eigenvectors; do not compare absolute thresholds. */
    for(i=0;i<NACTIVE;i++){
        s.a[i]*=1e9;s.d[i]*=1e9;s.c[i]*=1e9;s.b[i]*=1e9;
    }
    s.sigma=-0.25;
    if(run_psd(&s,"scaled PSD metric",7,NULL,expected))return 1;
    for(i=0;i<NACTIVE;i++){
        s.a[i]*=1e-18;s.d[i]*=1e-18;s.c[i]*=1e-18;s.b[i]*=1e-18;
    }
    if(run_psd(&s,"scaled PSD metric 1e-9",9,NULL,expected))return 1;

    /* Rank(B)<nev cannot produce enough independent finite eigenpairs.
     * We must not manufacture bogus eigendirections. */
    for(i=0;i<NACTIVE;i++)s.b[i]=(i==0?1.0:0.0);
    memset(&p,0,sizeof(p));p.n=N;p.transform=LIS_KRYLOV_SHIFT_INVERT_B;
    p.sigma=s.sigma;p.ctx=&s;p.apply_a=apply_a;
    p.apply_b=apply_b;p.solve=solve_shifted;
    lis_krylov_options_default(&opt);opt.nev=NMODE;opt.ncv=20;
    opt.maxiter=80;opt.allow_semidefinite_metric=1;
    err=lis_krylov_symmetric(&p,&opt,vals,vec,res,&st);
    if(err!=LIS_KRYLOV_NO_CONVERGENCE || st.nreturned>NMODE ||
       st.nconv>NMODE) {
        fprintf(stderr,"FAIL low rank: rc=%d returned=%d conv=%d\n",
                err,st.nreturned,st.nconv);
        return 1;
    }
    for(i=0;i<st.nreturned;i++){
        if(!isfinite(vals[i]) || !isfinite(res[i]))return 1;
    }
    puts("PASS low-rank PSD rejects missing eigenpairs");
    opt.allow_semidefinite_metric=2;
    if(lis_krylov_symmetric(&p,&opt,vals,vec,res,&st)!=
       LIS_KRYLOV_BAD_ARGUMENT) return 1;
    puts("PASS invalid PSD option validation");
    return 0;
}
