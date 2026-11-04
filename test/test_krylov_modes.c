/* Native Krylov tests: independently known generalized eigenpairs, mixed
 * spectral signs, standard mode, validation and callback failure handling.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 64
#define NEV 4
typedef struct { double a[N], b[N], shift[N]; int fail, solves, factor; } diag_ctx;
static int apply_a(void *ctx,int n,const double *x,double *y) {
    diag_ctx *d=(diag_ctx*)ctx;int i;
    if(d->fail)return -1;
    for(i=0;i<n;i++)y[i]=d->a[i]*x[i];
    return 0;
}
static int apply_b(void *ctx,int n,const double *x,double *y) {
    diag_ctx *d=(diag_ctx*)ctx;int i;
    for(i=0;i<n;i++)y[i]=d->b[i]*x[i];
    return 0;
}
static int solve(void *ctx,int n,const double *x,double *y) {
    diag_ctx *d=(diag_ctx*)ctx;int i;
    d->solves++;
    for(i=0;i<n;i++)y[i]=x[i]/d->shift[i];
    return 0;
}
static int check(int transform,double sigma,const double *expect,int nev) {
    diag_ctx d={0};
    lis_krylov_problem p;
    lis_krylov_options o={0};
    lis_krylov_stats s={0};
    double val[NEV],vec[N*NEV],res[NEV];
    int i,j,rc;
    lis_krylov_options_default(&o);
    o.nev=nev; o.ncv=24; o.maxiter=400; o.tol=1e-9;
    for(i=0;i<N;i++) {
        if(transform==LIS_KRYLOV_SHIFT_INVERT_B) {
            d.b[i]=1.+0.025*i;
            d.a[i]=d.b[i]*(4.+3.*i);
        } else if(transform==LIS_KRYLOV_SHIFT_INVERT_A) {
            double lambda=i%2==0?2.5+0.4*i:-(3.+0.5*i);
            d.a[i]=10.+i;
            d.b[i]=d.a[i]/lambda;
        } else {
            d.a[i]=(i%2? -2.: 1.)*(i+1);
            d.b[i]=1.;
        }
        d.shift[i]=d.a[i]-sigma*d.b[i];
    }
    p.n=N;p.transform=transform;p.sigma=sigma;p.ctx=&d;
    p.apply_a=apply_a;p.apply_b=apply_b;p.solve=solve;
    rc=lis_krylov_symmetric(&p,&o,val,vec,res,&s);
    if(rc||s.nconv!=nev||s.solve_count!=d.solves) {
        fprintf(stderr,"FAIL transform=%d sigma=%g rc=%d nconv=%d it=%d restart=%d solve=%d\n",transform,sigma,rc,s.nconv,s.iterations,s.restarts,s.solve_count);
        return 1;
    }
    for(j=0;j<nev;j++) {
        double best=1e100;
        for(i=0;i<nev;i++)if(fabs(expect[j]-val[i])<best)best=fabs(expect[j]-val[i]);
        if(best>1e-7||res[j]>1e-9) {
            fprintf(stderr,"FAIL: transform=%d eigenpair %d value=%.15g best=%g resid=%g\n",transform,j,val[j],best,res[j]);
            return 1;
        }
    }
    for(i=0;i<nev;i++)for(j=0;j<nev;j++) {
        double dot=0.;int k;
        for(k=0;k<N;k++) {
            double metric=(transform==LIS_KRYLOV_SHIFT_INVERT_A?d.a[k]:transform==LIS_KRYLOV_SHIFT_INVERT_B?d.b[k]:1.);
            dot+=vec[i*N+k]*vec[j*N+k]*metric;
        }
        if(fabs(dot-(i==j?1.:0.))>1e-7){fprintf(stderr,"B-orthog FAIL %d %d %g\n",i,j,dot);return 1;}
    }
    printf("PASS mode=%d sigma=%g nev=%d iter=%d restarts=%d solve=%d residual=%.3g\n",transform,sigma,nev,s.iterations,s.restarts,s.solve_count,s.max_residual);
    o.nev=N;
    if(lis_krylov_symmetric(&p,&o,val,vec,res,&s)!=LIS_KRYLOV_BAD_ARGUMENT)return 1;
    o.nev=nev;o.maxiter=1;
    if(lis_krylov_symmetric(&p,&o,val,vec,res,&s)!=LIS_KRYLOV_NO_CONVERGENCE)return 1;
    o.maxiter=400;o.which="XX";
    if(lis_krylov_symmetric(&p,&o,val,vec,res,&s)!=LIS_KRYLOV_BAD_ARGUMENT)return 1;
    o.which="LM";d.fail=1;
    if(lis_krylov_symmetric(&p,&o,val,vec,res,&s)!=LIS_KRYLOV_CALLBACK_ERROR)return 1;
    return 0;
}
int main(void){
    const double freq[]={4.,7.,10.,13.};
    const double pos[]={2.5,3.3,4.1};
    const double neg[]={-3.5,-4.5,-5.5};
    const double standard[]={-128.,-124.,-120.};
    return check(LIS_KRYLOV_SHIFT_INVERT_B,0.,freq,4)||
           check(LIS_KRYLOV_SHIFT_INVERT_A,2.,pos,3)||
           check(LIS_KRYLOV_SHIFT_INVERT_A,-3.,neg,3)||
           check(LIS_KRYLOV_STANDARD,0.,standard,3);
}
