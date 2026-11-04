/* Compare a standalone shifted-QR Lanczos prototype with the existing
 * symmetric solver on precisely the same three constructed problems.
 * This comparison does not change either public interface.
 */
#define main irl_prototype_main
#include "test_krylov_irl_prototype.c"
#undef main
#include "lis_krylov.h"

typedef struct {
    double A[MAXN][MAXN];
    int n;
    int calls;
} matrix_context;

static int baseline_op(void *ctx,int n,const double *x,double *y)
{
    matrix_context *c=(matrix_context *)ctx;
    int i,j;
    if(n!=c->n) return 1;
    c->calls++;
    for(i=0;i<n;i++){
        double sum=0.0;
        for(j=0;j<n;j++)sum+=c->A[i][j]*x[j];
        y[i]=sum;
    }
    return 0;
}
int main(void)
{
    static matrix_context c;
    static const double expected[3][3] = {
        {6.04069946274,5.95961541770,5.87976745035},
        {1.61549914705,1.60937686515,1.60522013312},
        {1.55986684425,1.55315504490,1.54939876293}
    };
    lis_krylov_problem p;
    lis_krylov_options o;
    lis_krylov_stats st;
    double values[3],resids[3],vectors[MAXN*3],irlvals[3],irlres[3];
    int kind,ops,restarts,rc,proto_rc,i;
    for(kind=0;kind<3;kind++){
        memset(&c,0,sizeof c);
        c.n=64;
        make_matrix(c.A,c.n,kind);
        memset(&p,0,sizeof p);
        lis_krylov_options_default(&o);
        p.n=c.n;
        p.ctx=&c;
        p.apply_a=baseline_op;
        p.transform=LIS_KRYLOV_STANDARD;
        o.nev=3;
        o.ncv=24;
        o.maxiter=150;
        o.tol=2e-8;
        o.which="LA";
        o.random_seed=17;
        rc=lis_krylov_symmetric(&p,&o,values,vectors,resids,&st);
        proto_rc=solve(c.A,c.n,24,9,3,2e-8,150,irlvals,irlres,&ops,&restarts);
        printf("COMPARE kind=%d baseline_rc=%d baseline_A=%d baseline_restarts=%d baseline_conv=%d irl_rc=%d irl_A=%d irl_restarts=%d\n",
               kind,rc,c.calls,st.restarts,st.nconv,proto_rc,ops,restarts);
        if(proto_rc!=1)return 1;
        for(i=0;i<3;i++)
            if(fabs(irlvals[i]-expected[kind][i])>2e-6 ||
               !isfinite(irlres[i]) || irlres[i]>2e-8) return 1;
        /* The existing solver has independent stopping rules. Do not
         * treat its callback counts as an apples-to-apples speedup. */
        if(rc==LIS_KRYLOV_SUCCESS){
            for(i=0;i<3;i++)
                if(fabs(values[i]-expected[kind][i])>2e-5 ||
                   !isfinite(resids[i]) || resids[i]>2e-8) return 1;
        }else if(rc!=LIS_KRYLOV_NO_CONVERGENCE) return 1;
    }
    puts("PASS baseline-versus-IRL physical eigenpair checks");
    return 0;
}
