/* Standalone benchmark: native restarted Krylov vs diagonal Davidson.
 * These methods use different starting subspaces and are NOT a controlled
 * performance comparison. Verify both original-problem residuals.
 */
#define main davidson_prototype_main
#include "test_krylov_davidson_prototype.c"
#undef main
#include "lis_krylov.h"

typedef struct {
    double a[DV_NMAX][DV_NMAX];
    int n, calls;
} NativeCtx;
static int native_apply(void *v,int n,const double *x,double *y)
{
    NativeCtx *c=(NativeCtx*)v;
    int i,j;
    if(n!=c->n)return 1;
    c->calls++;
    for(i=0;i<n;i++){
        double sum=0.0;
        for(j=0;j<n;j++)sum+=c->a[i][j]*x[j];
        y[i]=sum;
    }
    return 0;
}
static int check_eigenpair(const double a[DV_NMAX][DV_NMAX],
                           int n,const double *x,double lambda,
                           double claimed,double tol)
{
    int i,j;double r2=0,nax2=0,nx2=0;
    for(i=0;i<n;i++){
        double ax=0.0;
        for(j=0;j<n;j++)ax+=a[i][j]*x[j];
        r2+=(ax-lambda*x[i])*(ax-lambda*x[i]);
        nax2+=ax*ax;nx2+=x[i]*x[i];
    }
    {
        double den=sqrt(nax2)+fabs(lambda)*sqrt(nx2);
        double direct=sqrt(r2)/fmax(den,DBL_MIN);
        if(!isfinite(direct)||!isfinite(claimed)||
           direct>tol||direct>12*claimed+1e-11){
            fprintf(stderr,"FAIL independent physical residual %.17g reported %.17g\n",
                    direct,claimed);
            return 1;
        }
    }
    return 0;
}
int main(void)
{
    static NativeCtx c;
    static const double expected[3][3]={
        {6.04069946274,5.95961541770,5.87976745035},
        {1.61549914705,1.60937686515,1.60522013312},
        {1.55986684425,1.55315504490,1.54939876293}
    };
    double init[DV_NMAX]={0};
    int kind,i,j;
    for(kind=0;kind<3;kind++){
        Davidson *d=(Davidson*)malloc(sizeof(Davidson));
        lis_krylov_options o;
        lis_krylov_problem problem;
        lis_krylov_stats st;
        double vals[3]={0},res[3]={0},vecs[3*64]={0};
        int status,rc;
        if(!d)return 2;
        memset(&c,0,sizeof c);c.n=64;
        dv_matrix(c.a,64,kind);
        status=dv_solve(c.a,64,24,4,3,2e-8,220,d);
        memset(&problem,0,sizeof problem);
        lis_krylov_options_default(&o);
        memset(init,0,sizeof init);init[63]=1.0;
        problem.n=64;problem.transform=LIS_KRYLOV_STANDARD;
        problem.ctx=&c;problem.apply_a=native_apply;
        o.nev=3;o.ncv=24;o.maxiter=220;o.tol=2e-8;
        o.which="LA";o.initial_vector=init;
        rc=lis_krylov_symmetric(&problem,&o,vals,vecs,res,&st);
        printf("COMPARE kind=%d native_rc=%d native_calls=%d native_restarts=%d native_nconv=%d davidson_rc=%d davidson_calls=%d davidson_restarts=%d davidson_fallbacks=%d\n",
               kind,rc,c.calls,st.restarts,st.nconv,
               status,d->calls,d->restarts,d->fallbacks);
        if(status!=1)return 1;
        for(i=0;i<3;i++){
            if(fabs(d->result[i]-expected[kind][i])>3e-6||
               check_eigenpair(c.a,64,d->vectors[i],
                               d->result[i],d->residual[i],2e-8))return 1;
            for(j=0;j<i;j++)
                if(fabs(dv_dot(64,d->vectors[i],d->vectors[j]))>2e-7)return 1;
        }
        if(rc==LIS_KRYLOV_SUCCESS){
            if(st.nconv!=3)return 1;
            for(i=0;i<3;i++)
                if(fabs(vals[i]-expected[kind][i])>3e-5||
                   check_eigenpair(c.a,64,vecs+(size_t)i*64,
                                   vals[i],res[i],2e-8))return 1;
        }else if(rc!=LIS_KRYLOV_NO_CONVERGENCE)return 1;
        free(d);
    }
    puts("PASS independent native-versus-Davidson physical residual checks");
    return 0;
}
