/* Analytic non-diagonal spring-chain reference and factorization-reuse test. */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define N 160
#define NEV 5
#define PI 3.14159265358979323846
typedef struct { double pivot[N], c[N], tmp[N]; int factor,solves; } chain;
static int apply_a(void *context,int n,const double *x,double *y) {
    int i;(void)context;
    for(i=0;i<n;i++)y[i]=2*x[i]-(i?x[i-1]:0)-(i+1<n?x[i+1]:0);
    return 0;
}
static int apply_b(void *context,int n,const double *x,double *y) {
    int i;(void)context;
    for(i=0;i<n;i++)y[i]=x[i];
    return 0;
}
static void factor(chain *c) {
    int i;c->pivot[0]=2.;
    for(i=1;i<N;i++) {c->c[i-1]=-1./c->pivot[i-1];c->pivot[i]=2.+c->c[i-1];}
    c->factor++;
}
static int solve(void *context,int n,const double *x,double *y) {
    chain *c=(chain*)context;int i;
    c->tmp[0]=x[0];
    for(i=1;i<n;i++)c->tmp[i]=x[i]+c->tmp[i-1]/c->pivot[i-1];
    y[n-1]=c->tmp[n-1]/c->pivot[n-1];
    for(i=n-2;i>=0;i--)y[i]=c->tmp[i]/c->pivot[i]-c->c[i]*y[i+1];
    c->solves++;return 0;
}
int main(void) {
    chain c={0};
    lis_krylov_problem p={N,LIS_KRYLOV_SHIFT_INVERT_B,0.,&c,apply_a,apply_b,solve};
    lis_krylov_options o={0};
    lis_krylov_stats s={0};
    double eigen[NEV],evec[N*NEV],res[NEV];int i,j,rc;
    lis_krylov_options_default(&o);
    o.nev=NEV; o.ncv=25; o.maxiter=600; o.tol=1e-9;
    factor(&c);
    rc=lis_krylov_symmetric(&p,&o,eigen,evec,res,&s);
    if(rc||s.nconv!=NEV||c.factor!=1||c.solves!=s.solve_count){
        fprintf(stderr,"spring failed rc=%d conv=%d iter=%d restart=%d solve=%d\n",rc,s.nconv,s.iterations,s.restarts,s.solve_count);return 1;
    }
    for(i=0;i<NEV;i++) {
        double exact=2.-2.*cos((i+1)*PI/(N+1));double best=1e100;
        for(j=0;j<NEV;j++)if(fabs(eigen[j]-exact)<best)best=fabs(eigen[j]-exact);
        if(best>1e-6*fabs(exact) || res[i]>1e-9) {
            fprintf(stderr,"spring mode %d best err=%g res=%g value=%.15g\n",i,best,res[i],eigen[i]);return 1;
        }
    }
    printf("PASS spring N=%d nev=%d iter=%d restarts=%d solve=%d factorizations=%d maxres=%.3g\n",N,NEV,s.iterations,s.restarts,s.solve_count,c.factor,s.max_residual);
    return 0;
}
