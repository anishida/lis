/* Spectrum-aware restart regression and operator-count benchmark.
 * Compare the same fixture against develop and this change; verify
 * independently in the original eigenproblem, not only solver stats. */
#include "lis_krylov.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 96
#define NEV 4

typedef struct { double a[N], b[N], sigma; long ac,bcall,sc; } context;

static int action_a(void *ctx,int n,const double *x,double *y)
{
    context *c=(context*)ctx; int j;
    if(n!=N)return 1;
    ++c->ac;
    for(j=0;j<n;j++)y[j]=c->a[j]*x[j];
    return 0;
}
static int action_b(void *ctx,int n,const double *x,double *y)
{
    context *c=(context*)ctx; int j;
    if(n!=N)return 1;
    ++c->bcall;
    for(j=0;j<n;j++)y[j]=c->b[j]*x[j];
    return 0;
}
static int action_solve(void *ctx,int n,const double *x,double *y)
{
    context *c=(context*)ctx; int j;
    if(n!=N)return 1;
    ++c->sc;
    for(j=0;j<n;j++){
        double denom=c->a[j]-c->sigma*c->b[j];
        if(fabs(denom)<1e-13)return 1;
        y[j]=x[j]/denom;
    }
    return 0;
}
static int run(int spectrum,int ncv,int mode,int api)
{
    context c={0};
    lis_krylov_problem p={0};
    lis_krylov_options o;
    lis_krylov_stats st;
    lis_krylov_ext_options ex;
    lis_krylov_ext_stats ext_stats;
    lis_krylov_workspace *ws=NULL;
    double values[NEV]={0},residuals[NEV]={0},vectors[N*NEV]={0};
    int i,j,rc;
    for(i=0;i<N;i++){
        switch(spectrum){
        case 0: c.a[i]=1.0+0.037*i;break;
        case 1: c.a[i]=1.0+0.002*i+0.02*(i/8);break;
        case 2: c.a[i]=1.0+0.00002*i+0.1*(i/12);break;
        case 3: c.a[i]=1.0+0.000001*i+0.08*(i/16);break;
        default:c.a[i]=1.0+0.00015*i+0.05*(i/6);break;
        }
        c.b[i]=1.0+0.005*(i%7);
    }
    c.sigma=1.06;
    p.n=N;p.transform=mode;p.sigma=mode==LIS_KRYLOV_STANDARD?0.0:c.sigma;
    p.ctx=&c;p.apply_a=action_a;p.apply_b=action_b;p.solve=action_solve;
    lis_krylov_options_default(&o);
    o.nev=NEV;o.ncv=ncv;o.tol=1e-9;o.maxiter=160;
    o.random_seed=5;o.which="LM";o.hard_locking=0;
    lis_krylov_ext_options_default(&ex);
    lis_krylov_ext_stats_default(&ext_stats);
    /* Equal active and reserved capacities isolate Ritz retention
     * from growth. The legacy API remains an unchanged control. */
    if (api != 0) {
        ex.adaptive_ncv = 1;
        ex.ncv_max = ncv;
    }
    if(api==2){
        rc=lis_krylov_workspace_create(&p,&o,&ex,&ws);
        if(rc!=LIS_KRYLOV_SUCCESS)return 1;
    }
    rc=api==0?
       lis_krylov_symmetric(&p,&o,values,vectors,residuals,&st):
       lis_krylov_symmetric_ex(&p,&o,&ex,ws,values,vectors,residuals,&st,&ext_stats);
    lis_krylov_workspace_destroy(ws);
    if(rc!=LIS_KRYLOV_SUCCESS&&rc!=LIS_KRYLOV_NO_CONVERGENCE){
        fprintf(stderr,"FAIL solver rc=%d\n",rc);return 1;
    }
    if(st.nconv<0||st.nreturned>NEV||st.nconv>st.nreturned)return 1;
    for(i=0;i<st.nreturned;i++){
        double nr2=0,na2=0,nb2=0,den;
        for(j=0;j<N;j++){
            const double ax=c.a[j]*vectors[i*N+j];
            const double bx=(mode==LIS_KRYLOV_STANDARD?1.0:c.b[j])*vectors[i*N+j];
            nr2+=(ax-values[i]*bx)*(ax-values[i]*bx);
            na2+=ax*ax;nb2+=bx*bx;
        }
        den=sqrt(na2)+fabs(values[i])*sqrt(nb2);
        if(!isfinite(residuals[i]) || !isfinite(values[i]) ||
           (den>0 && sqrt(nr2)/den > 12*residuals[i]+1e-11)){
            fprintf(stderr,"FAIL physical residual %d %.17g\n",i,residuals[i]);return 1;
        }
    }
    printf("CASE %d %d %d %d %d %d %d %d %d %ld %ld %ld",
           spectrum,ncv,mode,api,rc,st.nconv,st.nreturned,
           st.iterations,st.restarts,c.ac,c.bcall,c.sc);
    for(i=0;i<st.nreturned;i++)printf(" %.17g %.17g",values[i],residuals[i]);
    puts("");
    return 0;
}
int main(void)
{
    int spectrum,ncv,mode,api;
    for(spectrum=0;spectrum<5;spectrum++)
        for(ncv=9;ncv<=17;ncv+=4)
            for(mode=0;mode<=2;mode++)
                for(api=0;api<3;api++)
                    if(run(spectrum,ncv,mode,api))return 1;
    puts("PASS spectrum-aware restart eigenpairs");
    return 0;
}
