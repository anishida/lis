/* BothEnds transformed Ritz selection and independent original-value
 * output ordering. This fixture never assumes source-level Ritz ordering
 * is the same as ordering of original generalized eigenvalues.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 16
#define NEV 4

typedef struct {
    double a[N], b[N], shift[N], lambda[N];
    int calls_a, calls_b, solves;
} diagonal_ctx;

static const double spectrum[N] = {
    -18.0, -13.0, -9.0, -4.0, 0.5, 1.0, 3.0, 4.0,
    5.0, 7.0, 11.0, 16.0, 23.0, 30.0, -2.0, -1.0
};

static int act_a(void *v, int n, const double *x, double *y)
{
    diagonal_ctx *c = (diagonal_ctx*)v;int i;
    if(n!=N)return 1;
    c->calls_a++;
    for(i=0;i<n;i++)y[i]=c->a[i]*x[i];
    return 0;
}
static int act_b(void *v, int n, const double *x, double *y)
{
    diagonal_ctx *c = (diagonal_ctx*)v;int i;
    if(n!=N)return 1;
    c->calls_b++;
    for(i=0;i<n;i++)y[i]=c->b[i]*x[i];
    return 0;
}
static int shifted_solve(void *v, int n, const double *x, double *y)
{
    diagonal_ctx *c = (diagonal_ctx*)v;int i;
    if(n!=N)return 1;
    c->solves++;
    for(i=0;i<n;i++){
        if(fabs(c->shift[i])<1e-12)return 1;
        y[i]=x[i]/c->shift[i];
    }
    return 0;
}
static void setup(diagonal_ctx *c, lis_krylov_problem *p, int mode)
{
    int i;
    memset(c,0,sizeof(*c));memset(p,0,sizeof(*p));
    p->n=N;p->ctx=c;p->transform=mode;
    p->apply_a=act_a;p->apply_b=act_b;p->solve=shifted_solve;
    p->sigma=mode==LIS_KRYLOV_SHIFT_INVERT_A ? 2.0 : 0.0;
    for(i=0;i<N;i++){
        c->lambda[i]=spectrum[i];
        if(mode==LIS_KRYLOV_STANDARD){
            c->a[i]=c->lambda[i];
            c->b[i]=1.0;
        }else if(mode==LIS_KRYLOV_SHIFT_INVERT_B){
            c->b[i]=1.0+0.015*i;
            c->a[i]=c->lambda[i]*c->b[i];
        }else{
            c->a[i]=7.0+i;
            c->b[i]=c->a[i]/c->lambda[i];
        }
        c->shift[i]=c->a[i]-p->sigma*c->b[i];
    }
}
static double mu_of(double lambda,int mode,double sigma)
{
    if(mode==LIS_KRYLOV_STANDARD)return lambda;
    if(mode==LIS_KRYLOV_SHIFT_INVERT_B)return 1.0/(lambda-sigma);
    return lambda/(lambda-sigma);
}
static void reference_expected(const diagonal_ctx *c,int mode,
                               double sigma,const char *spec,
                               double out[NEV])
{
    int order[N],selected[NEV],i,j,tmp;
    int both=spec[0]=='B',kind=spec[0],magnitude=spec[1]=='M';
    int asc=strstr(spec,":ASC")!=NULL,desc=strstr(spec,":DESC")!=NULL;
    for(i=0;i<N;i++)order[i]=i;
    for(i=1;i<N;i++){
        int key=order[i];
        double x=mu_of(c->lambda[key],mode,sigma);
        j=i-1;
        while(j>=0){
            double y=mu_of(c->lambda[order[j]],mode,sigma);
            double xa=magnitude?fabs(x):x,ya=magnitude?fabs(y):y;
            int before=both||kind=='S' ? xa<ya : xa>ya;
            if(!before)break;
            order[j+1]=order[j];j--;
        }
        order[j+1]=key;
    }
    if(both){
        for(i=0;i<NEV;i++)
            selected[i]=(i&1)?order[i/2]:order[N-1-i/2];
    }else for(i=0;i<NEV;i++)selected[i]=order[i];
    for(i=0;i<NEV;i++)out[i]=c->lambda[selected[i]];
    if(asc||desc)for(i=1;i<NEV;i++){
        double key=out[i];
        j=i-1;
        while(j>=0&&(asc?out[j]>key:out[j]<key)){
            out[j+1]=out[j];j--;
        }
        out[j+1]=key;
    }
    (void)tmp;
}
static int verify_physical(const diagonal_ctx *c,int mode,
                           double lambda,const double *x,double reported,
                           double tol)
{
    int i;double ra=0,nax=0,nbx=0,nx=0,r,den,scale=0.0;
    for(i=0;i<N;i++){
        double ax=c->a[i]*x[i];
        double bx=mode==LIS_KRYLOV_STANDARD?x[i]:c->b[i]*x[i];
        double diff=ax-lambda*bx;
        ra+=diff*diff;nax+=ax*ax;nbx+=bx*bx;nx+=x[i]*x[i];
    }
    for(i=0;i<N;i++)if(fabs(c->a[i])>scale)scale=fabs(c->a[i]);
    den=sqrt(nax)+fabs(lambda)*sqrt(nbx);
    if(scale*sqrt(nx)>den)den=scale*sqrt(nx);
    r=den>0?sqrt(ra)/den:sqrt(ra);
    if(!isfinite(lambda)||!isfinite(reported)||!isfinite(r)||
       r>12.0*reported+1e-10||
       (reported<=tol&&r>5.0e-8)){
        fprintf(stderr,"FAIL physical lambda %.17g reported %.17g actual %.17g\n",
                lambda,reported,r);
        return 1;
    }
    return 0;
}

static int run(int mode,int locking,int api,const char *spec)
{
    diagonal_ctx c;
    lis_krylov_problem p;
    lis_krylov_options o;
    lis_krylov_stats stats;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *ws=NULL;
    double expected[NEV],values[NEV],residuals[NEV],vectors[NEV*N];
    int rc,i,j;
    setup(&c,&p,mode);reference_expected(&c,mode,p.sigma,spec,expected);
    lis_krylov_options_default(&o);
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    o.nev=NEV;o.ncv=15;o.maxiter=380;o.tol=5e-9;
    o.which=spec;o.hard_locking=locking;o.random_seed=11;
    for(i=0;i<NEV;i++){values[i]=-999;residuals[i]=-999;}
    memset(vectors,0,sizeof(vectors));
    if(api==2){
        rc=lis_krylov_workspace_create(&p,&o,&ext,&ws);
        if(rc!=LIS_KRYLOV_SUCCESS||!ws)return 1;
    }
    rc=api==0 ?
       lis_krylov_symmetric(&p,&o,values,vectors,residuals,&stats):
       lis_krylov_symmetric_ex(&p,&o,&ext,ws,values,vectors,
                               residuals,&stats,&extra);
    lis_krylov_workspace_destroy(ws);
    if(rc!=LIS_KRYLOV_SUCCESS||stats.nconv!=NEV||
       stats.nreturned!=NEV){
        fprintf(stderr,"FAIL solver mode=%d locking=%d api=%d spec=%s rc=%d returned=%d conv=%d iter=%d\n",
                mode,locking,api,spec,rc,stats.nreturned,stats.nconv,stats.iterations);
        return 1;
    }
    for(i=0;i<NEV;i++){
        if(fabs(values[i]-expected[i])>1e-5){
            fprintf(stderr,"FAIL ordering mode=%d locking=%d api=%d spec=%s index=%d got=%.15g want=%.15g\n",
                    mode,locking,api,spec,i,values[i],expected[i]);
            return 1;
        }
        if(verify_physical(&c,mode,values[i],vectors+(size_t)i*N,
                           residuals[i],o.tol))return 1;
        for(j=0;j<i;j++){
            double dot=0;int k;
            for(k=0;k<N;k++){
                double metric=(mode==LIS_KRYLOV_SHIFT_INVERT_A?
                               c.a[k]:mode==LIS_KRYLOV_SHIFT_INVERT_B?
                               c.b[k]:1.0);
                dot+=metric*vectors[(size_t)i*N+k]*vectors[(size_t)j*N+k];
            }
            if(fabs(dot)>1e-6){
                fprintf(stderr,"FAIL metric orthogonality %s %g\n",spec,dot);
                return 1;
            }
        }
    }
    printf("PASS CASE mode=%d hard=%d api=%d which=%s A=%d B=%d solves=%d\n",
           mode,locking,api,spec,c.calls_a,c.calls_b,c.solves);
    return 0;
}
static int check_invalid(void)
{
    static const char *bad[]={"","B","BEx","BE:BAD","BE:ASCjunk",
                               "LM:ASC:BAD","BE:","SM:"};
    diagonal_ctx c;
    lis_krylov_problem p;
    lis_krylov_options o;
    lis_krylov_stats st;
    double values[NEV],vec[N*NEV],res[NEV];
    int i;
    setup(&c,&p,LIS_KRYLOV_STANDARD);
    lis_krylov_options_default(&o);
    o.nev=NEV;o.ncv=15;o.maxiter=40;
    for(i=0;i<(int)(sizeof bad/sizeof bad[0]);i++){
        o.which=bad[i];
        if(lis_krylov_symmetric(&p,&o,values,vec,res,&st)!=
           LIS_KRYLOV_BAD_ARGUMENT){
            fprintf(stderr,"FAIL accepted invalid which '%s'\n",bad[i]);
            return 1;
        }
    }
    return 0;
}
/* A reusable workspace must not retain any selection or output
 * permutation from a preceding solve. In particular, BE's temporary
 * interleaving buffer shares storage with convergence flags, and
 * :ASC/:DESC reuses output permutation scratch at publication. */
static int check_repeated_workspace(void)
{
    diagonal_ctx c;
    lis_krylov_problem p;
    lis_krylov_options o;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *ws=NULL;
    lis_krylov_stats st;
    double values[NEV],vectors[N*NEV],residuals[NEV],expected[NEV];
    const char *specs[]={"BE","BE:ASC","BE:DESC","LA:ASC"};
    int i,j,rc;
    setup(&c,&p,LIS_KRYLOV_STANDARD);
    lis_krylov_options_default(&o);
    lis_krylov_ext_options_default(&ext);
    o.nev=NEV;o.ncv=15;o.maxiter=380;o.tol=5e-9;
    o.hard_locking=1;o.random_seed=11;
    rc=lis_krylov_workspace_create(&p,&o,&ext,&ws);
    if(rc||!ws)return 1;
    for(i=0;i<4;i++){
        o.which=specs[i];
        lis_krylov_ext_stats_default(&extra);
        rc=lis_krylov_symmetric_ex(&p,&o,&ext,ws,
                                   values,vectors,residuals,&st,&extra);
        if(rc!=LIS_KRYLOV_SUCCESS||st.nreturned!=NEV){
            fprintf(stderr,"FAIL workspace reuse with %s rc=%d returned=%d\n",
                    o.which,rc,st.nreturned);
            lis_krylov_workspace_destroy(ws);
            return 1;
        }
        reference_expected(&c,p.transform,p.sigma,o.which,expected);
        for(j=0;j<NEV;j++)if(fabs(values[j]-expected[j])>1e-5||
            verify_physical(&c,p.transform,values[j],
                            vectors+(size_t)j*N,residuals[j],o.tol)){
            fprintf(stderr,"FAIL reused workspace changed results for %s\n",
                    o.which);
            lis_krylov_workspace_destroy(ws);
            return 1;
        }
    }
    lis_krylov_workspace_destroy(ws);
    puts("PASS repeated workspace with alternating selection and output sorts");
    return 0;
}
/* Partial results must be reordered as well, but only on a nonfatal
 * NO_CONVERGENCE return. Early termination leaves the requested
 * source modes untouched; no uninitialized output is inspected. */
static int check_partial_output(void)
{
    diagonal_ctx c;
    lis_krylov_problem p;
    lis_krylov_options o;
    lis_krylov_stats st;
    double values[NEV],vectors[N*NEV],residuals[NEV];
    int rc,i;
    setup(&c,&p,LIS_KRYLOV_STANDARD);
    lis_krylov_options_default(&o);
    o.which="BE:ASC";o.nev=NEV;o.ncv=15;
    o.maxiter=3;o.tol=5e-9;
    rc=lis_krylov_symmetric(&p,&o,values,vectors,residuals,&st);
    if(rc!=LIS_KRYLOV_NO_CONVERGENCE||
       st.nreturned<0||st.nreturned>NEV)return 1;
    for(i=0;i<st.nreturned;i++){
        if(i&&values[i]<values[i-1]){
            fprintf(stderr,"FAIL unsorted partial result\n");
            return 1;
        }
        if(verify_physical(&c,p.transform,values[i],
                           vectors+(size_t)i*N,residuals[i],o.tol))
            return 1;
    }
    puts("PASS original-value sorting on partial NO_CONVERGENCE");
    return 0;
}

int main(void)
{
    static const char *specs[]={"BE","BE:ASC","BE:DESC","LA:ASC",
                                 "SA:DESC","LM:ASC","SM:DESC"};
    int mode,locking,api,which;
    if(check_invalid())return 1;
    if(check_repeated_workspace()||check_partial_output())return 1;
    for(mode=0;mode<=2;mode++)
        for(locking=0;locking<=1;locking++)
            for(api=0;api<=2;api++)
                for(which=0;which<(int)(sizeof specs/sizeof specs[0]);which++)
                    if(run(mode,locking,api,specs[which]))return 1;
    puts("PASS BothEnds and independent original-eigenvalue sorting");
    return 0;
}
