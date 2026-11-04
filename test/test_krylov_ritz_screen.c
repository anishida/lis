/* Projected Ritz prefilter: verify full physical results and require
 * fewer user matrix callbacks in deterministic hard-locking fixtures.
 * This regression is independent of any removed CI comparison script.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 64
#define NEV 4

typedef struct {
    double diagonal[N];
    int calls_a;
} context;

static int apply_a(void *data, int n, const double *x, double *y)
{
    context *c = (context *)data;
    int i;
    if (n != N) return 1;
    ++c->calls_a;
    for (i = 0; i < n; ++i) y[i] = c->diagonal[i] * x[i];
    return 0;
}

static int verify(const context *c, const double *v,
                  double lambda, double reported, double tol, double scale)
{
    double nr=0.0, na=0.0, nx=0.0, den;
    int i;
    for (i=0; i<N; ++i) {
        double av=c->diagonal[i]*v[i], r=av-lambda*v[i];
        nr += r*r;
        na += av*av;
        nx += v[i]*v[i];
    }
    den=sqrt(na)+fabs(lambda)*sqrt(nx);
    if (scale*sqrt(nx)>den) den=scale*sqrt(nx);
    if (!isfinite(reported) || !isfinite(lambda) ||
        reported<0.0 || reported>1.0 ||
        sqrt(nr)>den*(8.0*reported+1e-11) ||
        (reported<=tol && sqrt(nr)>1e-8*fmax(den,1.0))) {
        fprintf(stderr,"FAIL true residual lambda=%g reported=%g\n",
                lambda,reported);
        return 1;
    }
    return 0;
}

static int test(int long_run, int api, double scale)
{
    context c;
    lis_krylov_problem p={0};
    lis_krylov_options opt;
    lis_krylov_stats st;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *workspace=NULL;
    double values[NEV], vectors[N*NEV], residuals[NEV];
    int rc, i;
    memset(&c,0,sizeof(c));
    for(i=0;i<N;++i) c.diagonal[i]=1.0+0.027*i;
    p.n=N; p.transform=LIS_KRYLOV_STANDARD;
    p.ctx=&c; p.apply_a=apply_a;
    lis_krylov_options_default(&opt);
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    opt.nev=NEV;
    opt.ncv=long_run?24:10;
    opt.maxiter=long_run?250:24;
    opt.tol=1e-9;
    opt.which="LA";
    opt.random_seed=3;
    opt.hard_locking=1;
    opt.residual_scale=scale;

    if(api==2) {
        rc=lis_krylov_workspace_create(&p,&opt,&ext,&workspace);
        if(rc!=LIS_KRYLOV_SUCCESS || !workspace) {
            fprintf(stderr,"FAIL workspace allocation %d\n",rc);
            return 1;
        }
    }
    if(api==0)
        rc=lis_krylov_symmetric(&p,&opt,values,vectors,residuals,&st);
    else
        rc=lis_krylov_symmetric_ex(&p,&opt,&ext,workspace,
                                   values,vectors,residuals,&st,&extra);
    lis_krylov_workspace_destroy(workspace);

    if(rc!=LIS_KRYLOV_SUCCESS && rc!=LIS_KRYLOV_NO_CONVERGENCE) {
        fprintf(stderr,"FAIL solver rc=%d api=%d long=%d\n",rc,api,long_run);
        return 1;
    }
    if(st.nreturned<0 || st.nreturned>NEV ||
       st.nconv<0 || st.nconv>st.nreturned ||
       (long_run && (rc!=LIS_KRYLOV_SUCCESS || st.nconv!=NEV))) {
        fprintf(stderr,"FAIL returned/count rc=%d nconv=%d nreturned=%d\n",
                rc,st.nconv,st.nreturned);
        return 1;
    }
    for(i=0;i<st.nreturned;++i) {
        if(verify(&c,vectors+(size_t)i*N,
                  values[i],residuals[i],opt.tol,scale)) return 1;
    }
    /* Unoptimized develop: 221/1228 A callbacks for short/long.
     * Guarded prefilter: 131/804 calls. These budgets leave room for
     * compiler variation but still fail on a lost optimization. */
    if(scale==4.0 && c.calls_a>(long_run?1000:180)) {
        fprintf(stderr,"FAIL callback budget api=%d long=%d A=%d\n",
                api,long_run,c.calls_a);
        return 1;
    }
    printf("PASS Ritz prefilter api=%d long=%d scale=%g A=%d rc=%d nconv=%d\n",
           api,long_run,scale,c.calls_a,rc,st.nconv);
    return 0;
}

int main(void)
{
    int api,long_run,k;
    const double scales[] = {4.0, 0.25, 25.0};
    /* Verify both ordinary and widely different fixed positive scales.
     * The prefilter must use exactly the physical residual denominator
     * instead of assuming that residual_scale is always 4.0.
     * Keep callback budgets for the original calibrated fixture only. */
    for(api=0;api<3;++api)
        for(long_run=0;long_run<2;++long_run)
            for(k=0;k<3;++k)
                if(test(long_run,api,scales[k])) return 1;
    return 0;
}
