/* Independent C99 diagonal-preconditioned Davidson experiment.
 * Not a production LIS solver. No Spectra code is copied.
 * Standard real symmetric problems only; maximum dimensions are fixed
 * by this test harness, not by the production Krylov workspace.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DV_NMAX 72
#define DV_MMAX 28
#define DV_NEV 3

typedef struct {
    int n, m, cap, nev, calls, restarts, iterations, fallbacks;
    int oblique_steps, oblique_fallbacks;
    double v[DV_MMAX][DV_NMAX];
    double av[DV_MMAX][DV_NMAX];
    double d[DV_NMAX];
    double result[DV_NEV], residual[DV_NEV], vectors[DV_NEV][DV_NMAX];
} Davidson;

static double dv_dot(int n, const double *x, const double *y)
{
    int j; double sum=0.0;
    for(j=0;j<n;j++) sum+=x[j]*y[j];
    return sum;
}
static double dv_norm(int n, const double *x)
{
    return sqrt(fmax(0.0,dv_dot(n,x,x)));
}
static void dv_matvec(Davidson *p, const double a[DV_NMAX][DV_NMAX],
                      const double *x, double *y)
{
    int i,j;
    for(i=0;i<p->n;i++){
        double sum=0.0;
        for(j=0;j<p->n;j++) sum+=a[i][j]*x[j];
        y[i]=sum;
    }
    p->calls++;
}
/* Two-pass MGS prevents repeated and tightly clustered eigenvectors
 * from silently entering an ill-conditioned search subspace. */
static int dv_append(Davidson *p,const double a[DV_NMAX][DV_NMAX],
                     double *v)
{
    int pass,i,j; double norm0=dv_norm(p->n,v),norm;
    if(!isfinite(norm0)||norm0<DBL_MIN||p->m>=p->cap) return 1;
    for(pass=0;pass<2;pass++)for(j=0;j<p->m;j++){
        double proj=dv_dot(p->n,p->v[j],v);
        for(i=0;i<p->n;i++)v[i]-=proj*p->v[j][i];
    }
    norm=dv_norm(p->n,v);
    if(!isfinite(norm)||norm<1e-11*norm0)return 1;
    for(i=0;i<p->n;i++)p->v[p->m][i]=v[i]/norm;
    dv_matvec(p,a,p->v[p->m],p->av[p->m]);
    p->m++;
    return 0;
}
/* Symmetric Jacobi diagonalization of the small projected matrix.
 * The eigenvectors are columns of u. This routine is deliberately
 * self-contained and uses no external numerical dependencies. */
static int dv_eigen(int m,double h[DV_MMAX][DV_MMAX],
                    double values[DV_MMAX],double u[DV_MMAX][DV_MMAX])
{
    int i,j,p,q,step,t;
    memset(u,0,(size_t)DV_MMAX*DV_MMAX*sizeof(double));
    for(i=0;i<m;i++)u[i][i]=1.0;
    for(step=0;step<50*m*m;step++){
        double largest=0.0;
        p=0;q=0;
        for(i=0;i<m;i++)for(j=i+1;j<m;j++)
            if(fabs(h[i][j])>largest){
                largest=fabs(h[i][j]);p=i;q=j;
            }
        if(largest<1e-13)break;
        {
            double app=h[p][p],aqq=h[q][q],apq=h[p][q];
            double tau=(aqq-app)/(2.0*apq);
            double z=copysign(1.0,tau)/
                     (fabs(tau)+hypot(1.0,tau));
            double c=1.0/sqrt(1.0+z*z),s=z*c;
            for(t=0;t<m;t++)if(t!=p&&t!=q){
                double x=h[t][p],y=h[t][q];
                h[p][t]=h[t][p]=c*x-s*y;
                h[q][t]=h[t][q]=s*x+c*y;
            }
            h[p][p]=app-z*apq;
            h[q][q]=aqq+z*apq;
            h[p][q]=h[q][p]=0.0;
            for(t=0;t<m;t++){
                double x=u[t][p],y=u[t][q];
                u[t][p]=c*x-s*y;
                u[t][q]=s*x+c*y;
            }
        }
    }
    if(step==50*m*m)return -1;
    for(i=0;i<m;i++)values[i]=h[i][i];
    return 0;
}
/* One-step diagonal approximation to a Jacobi-Davidson correction:
 *   M = regularized(diag(A) - theta)
 *   t = -M^(-1)r + epsilon M^(-1)x
 *   epsilon = (x^T M^(-1)r) / (x^T M^(-1)x).
 * This changes the expansion direction even after orthogonalization
 * against V; it is NOT a full JD inner solve. If the denominator is
 * cancellation-prone for an indefinite M, use ordinary Davidson.
 */
static int dv_correction(int n,const double *diag,double theta,
                         const double *x,const double *r,double *out,
                         int *oblique_used)
{
    double inv_x[DV_NMAX],plain[DV_NMAX];
    double scale=1.0,cutoff,norm_plain,norm,projection;
    long double numerator=0.0L,denominator=0.0L,abs_denominator=0.0L;
    int i;
    if(oblique_used)*oblique_used=0;
    if(n<1||n>DV_NMAX||!isfinite(theta))return -1;
    for(i=0;i<n;i++){
        if(!isfinite(diag[i])||!isfinite(x[i])||!isfinite(r[i]))return -1;
        if(fabs(diag[i])>scale)scale=fabs(diag[i]);
    }
    cutoff=1e-6*(scale+fabs(theta));
    if(!isfinite(cutoff)||cutoff<=0.0)return -1;
    for(i=0;i<n;i++){
        double delta=diag[i]-theta;
        if(fabs(delta)<cutoff)
            delta=copysign(cutoff,delta==0.0?1.0:delta);
        plain[i]=-r[i]/delta;
        inv_x[i]=x[i]/delta;
        if(!isfinite(plain[i])||!isfinite(inv_x[i]))return -1;
        numerator+=(long double)x[i]*(long double)(-plain[i]);
        denominator+=(long double)x[i]*(long double)inv_x[i];
        abs_denominator+=fabsl((long double)x[i]*(long double)inv_x[i]);
    }
    norm_plain=dv_norm(n,plain);
    if(!isfinite(norm_plain))return -1;
    if(abs_denominator>0.0L &&
       fabsl(denominator)>128.0L*DBL_EPSILON*abs_denominator){
        double epsilon=(double)(numerator/denominator);
        int valid=isfinite(epsilon);
        if(valid){
            for(i=0;i<n;i++){
                out[i]=plain[i]+epsilon*inv_x[i];
                if(!isfinite(out[i]))valid=0;
            }
        }
        if(valid){
            norm=dv_norm(n,out);
            if(isfinite(norm)&&norm>1e-11*fmax(norm_plain,DBL_MIN)){
                if(oblique_used)*oblique_used=1;
                return 0;
            }
        }
    }
    /* Indefinite or nearly singular projected diagonal: do not divide
     * by a cancelled x^T M^(-1)x. Standard Davidson fallback. */
    projection=dv_dot(n,x,plain);
    for(i=0;i<n;i++)out[i]=plain[i]-projection*x[i];
    norm=dv_norm(n,out);
    return isfinite(norm)&&norm>DBL_MIN?0:1;
}
static void dv_restart(Davidson *p,int keep,int m,
                       const int *order,const double u[DV_MMAX][DV_MMAX])
{
    double newv[DV_MMAX][DV_NMAX],newav[DV_MMAX][DV_NMAX];
    int i,j,k;
    memset(newv,0,sizeof newv);
    memset(newav,0,sizeof newav);
    for(i=0;i<keep;i++)for(j=0;j<m;j++){
        double c=u[j][order[i]];
        for(k=0;k<p->n;k++){
            newv[i][k]+=c*p->v[j][k];
            newav[i][k]+=c*p->av[j][k];
        }
    }
    for(i=0;i<keep;i++){
        memcpy(p->v[i],newv[i],(size_t)p->n*sizeof(double));
        memcpy(p->av[i],newav[i],(size_t)p->n*sizeof(double));
    }
    p->m=keep;
    p->restarts++;
}
/* For the prototype, the target is the largest algebraic eigenvalues.
 * Return 1 on verified convergence, 0 on budget exhaustion and negative
 * on numerical failure. The authoritative convergence test uses A*x.
 */
static int dv_solve(const double a[DV_NMAX][DV_NMAX],int n,int cap,
                    int initial,int nev,double tol,int maxit,
                    Davidson *p)
{
    double h[DV_MMAX][DV_MMAX],eval[DV_MMAX],u[DV_MMAX][DV_MMAX];
    double x[DV_NEV][DV_NMAX],ax[DV_NEV][DV_NMAX];
    double residual[DV_NEV][DV_NMAX],correction[DV_NMAX];
    int order[DV_MMAX],sorteddiag[DV_NMAX],i,j,k,t,converged,rc;
    if(!p||!a||n<4||n>DV_NMAX||cap<=nev+2||
       cap>DV_MMAX||cap>n||initial<1||initial>=cap||
       initial>n||nev<1||nev>DV_NEV||nev>=n||
       !isfinite(tol)||tol<=0.0||maxit<=0)return -10;
    memset(p,0,sizeof(*p));
    p->n=n;p->cap=cap;p->nev=nev;
    for(i=0;i<n;i++){
        p->d[i]=a[i][i];sorteddiag[i]=i;
        if(!isfinite(p->d[i]))return -10;
    }
    /* Start with a few coordinate vectors selected by the diagonal.
     * This is an intentional advantage of the diagonal method and must
     * be accounted for in comparisons with generic Krylov solvers. */
    for(i=0;i<n;i++)for(j=i+1;j<n;j++)
        if(p->d[sorteddiag[j]]>p->d[sorteddiag[i]]){
            t=sorteddiag[i];sorteddiag[i]=sorteddiag[j];sorteddiag[j]=t;
        }
    for(i=0;i<initial;i++){
        memset(correction,0,(size_t)n*sizeof(double));
        correction[sorteddiag[i]]=1.0;
        if(dv_append(p,a,correction))return -11;
    }
    for(t=0;t<maxit;t++){
        int candidate=-1;
        p->iterations=t+1;
        memset(h,0,sizeof h);
        for(i=0;i<p->m;i++)for(j=i;j<p->m;j++)
            h[i][j]=h[j][i]=0.5*
              (dv_dot(n,p->v[i],p->av[j])+
               dv_dot(n,p->v[j],p->av[i]));
        rc=dv_eigen(p->m,h,eval,u);
        if(rc)return -12;
        for(i=0;i<p->m;i++)order[i]=i;
        for(i=0;i<p->m;i++)for(j=i+1;j<p->m;j++)
            if(eval[order[j]]>eval[order[i]]){
                k=order[i];order[i]=order[j];order[j]=k;
            }
        converged=0;
        memset(x,0,sizeof x);
        memset(ax,0,sizeof ax);
        for(i=0;i<nev;i++){
            int idx=order[i];
            double nr,scale;
            for(j=0;j<p->m;j++)for(k=0;k<n;k++){
                x[i][k]+=u[j][idx]*p->v[j][k];
                ax[i][k]+=u[j][idx]*p->av[j][k];
            }
            for(k=0;k<n;k++)
                residual[i][k]=ax[i][k]-eval[idx]*x[i][k];
            nr=dv_norm(n,residual[i]);
            scale=dv_norm(n,ax[i])+fabs(eval[idx])*dv_norm(n,x[i]);
            p->result[i]=eval[idx];
            p->residual[i]=nr/fmax(scale,DBL_MIN);
            if(!isfinite(p->residual[i]))return -13;
            if(p->residual[i]<=tol)converged++;
            else if(candidate<0)candidate=i;
            memcpy(p->vectors[i],x[i],(size_t)n*sizeof(double));
        }
        if(converged==nev){
            for(i=0;i<nev;i++)for(j=0;j<i;j++)
                if(fabs(dv_dot(n,x[i],x[j]))>1e-7)return -14;
            return 1;
        }
        if(t==maxit-1)break;
        if(candidate<0)return -15;
        {
            int idx=order[candidate],ok,used_oblique=0;
            ok=dv_correction(n,p->d,eval[idx],x[candidate],
                             residual[candidate],correction,&used_oblique);
            if(ok<0)return -16;
            if(used_oblique)++p->oblique_steps;
            else ++p->oblique_fallbacks;
            if(ok){
                memcpy(correction,residual[candidate],
                       (size_t)n*sizeof(double));
                p->fallbacks++;
            }
        }
        if(p->m>=p->cap){
            int keep=p->cap/2;
            if(keep<nev+1)keep=nev+1;
            if(keep>=p->cap)keep=p->cap-1;
            dv_restart(p,keep,p->m,order,u);
        }
        if(dv_append(p,a,correction)){
            /* A diagonal correction can become parallel to the search
             * basis near a repeated Ritz value. Fall back to its physical
             * residual and finally a deterministic independent seed. */
            p->fallbacks++;
            memcpy(correction,residual[candidate],
                   (size_t)n*sizeof(double));
            if(dv_append(p,a,correction)){
                int accepted=0;
                for(i=0;i<n;i++){
                    memset(correction,0,(size_t)n*sizeof(double));
                    correction[sorteddiag[i]]=1.0;
                    if(!dv_append(p,a,correction)){
                        accepted=1;break;
                    }
                }
                if(!accepted)return -17;
            }
        }
    }
    return 0;
}

static void dv_matrix(double a[DV_NMAX][DV_NMAX],int n,int kind)
{
    int i,j;
    for(i=0;i<DV_NMAX;i++)for(j=0;j<DV_NMAX;j++)a[i][j]=0;
    for(i=0;i<n;i++){
        if(kind==0)a[i][i]=1.0+0.08*i;
        else if(kind==1)a[i][i]=1.0+0.0001*i+0.06*(i/6);
        else a[i][i]=1.0+0.00001*i+0.11*(i/12);
        if(i+1<n)a[i][i+1]=a[i+1][i]=0.0025*(1+(i%5));
        if(i+3<n)a[i][i+3]=a[i+3][i]=0.0005;
    }
}
int main(void)
{
    static double a[DV_NMAX][DV_NMAX];
    static const double expected[3][3]={
        {6.04069946274,5.95961541770,5.87976745035},
        {1.61549914705,1.60937686515,1.60522013312},
        {1.55986684425,1.55315504490,1.54939876293}
    };
    int kind,cap,init,i,successful=0,failures=0;
    for(kind=0;kind<3;kind++)for(cap=12;cap<=24;cap+=4)
        for(init=3;init<=5;init++){
            Davidson *p=(Davidson*)malloc(sizeof(Davidson));
            int status;
            if(!p)return 2;
            dv_matrix(a,64,kind);
            status=dv_solve(a,64,cap,init,3,2e-8,220,p);
            printf("CASE kind=%d cap=%d init=%d rc=%d ops=%d restarts=%d iter=%d fallbacks=%d oblique=%d standard=%d",
                   kind,cap,init,status,p->calls,p->restarts,
                   p->iterations,p->fallbacks,p->oblique_steps,p->oblique_fallbacks);
            for(i=0;i<3;i++)printf(" lambda%d=%.13g residual%d=%.3g",
                   i,p->result[i],i,p->residual[i]);
            puts("");
            if(status==1){
                for(i=0;i<3;i++)
                    if(fabs(p->result[i]-expected[kind][i])>3e-6 ||
                       !isfinite(p->residual[i])||p->residual[i]>2e-8){
                        status=-18;break;
                    }
            }
            if(status==1)successful++;
            else {
                fprintf(stderr,"FAIL Davidson case kind=%d cap=%d init=%d status=%d\n",
                        kind,cap,init,status);
                failures++;
            }
            free(p);
        }
    /* Equal diagonal entries make the diagonal correction a scaled
     * physical residual. This checks the difficult, weak-preconditioner
     * regime rather than only near-diagonal matrices. A symmetric
     * tridiagonal Toeplitz matrix has exact eigenvalues
     * a + 2*b*cos(j*pi/(n+1)), for j=1,...,n. */
    {
        const double pi=3.14159265358979323846264338327950288;
        int trial,challenged=0;
        for(trial=0;trial<2;trial++){
            int n=trial?40:32,cap=trial?20:16;
            double offdiag=trial?0.30:0.20;
            Davidson *p=(Davidson*)malloc(sizeof(Davidson));
            int j,k,rc;
            if(!p)return 2;
            for(j=0;j<DV_NMAX;j++)for(k=0;k<DV_NMAX;k++)
                a[j][k]=0.0;
            for(j=0;j<n;j++){
                a[j][j]=1.0;
                if(j+1<n)a[j][j+1]=a[j+1][j]=offdiag;
            }
            rc=dv_solve(a,n,cap,3,3,1e-7,450,p);
            printf("TOEPLITZ n=%d cap=%d rc=%d calls=%d restarts=%d fallbacks=%d\n",
                   n,cap,rc,p->calls,p->restarts,p->fallbacks);
            if(rc!=1){fprintf(stderr,"FAIL Toeplitz convergence\n");free(p);return 4;}
            for(j=0;j<3;j++){
                double lambda=1.0+2.0*offdiag*cos(pi*(j+1)/(n+1));
                double norm2=0.0,res2=0.0,ax2=0.0;
                if(fabs(p->result[j]-lambda)>4e-6) {
                    fprintf(stderr,"FAIL Toeplitz eigenvalue %d %.15g vs %.15g\n",
                            j,p->result[j],lambda);
                    free(p);return 5;
                }
                for(k=0;k<n;k++){
                    double x=p->vectors[j][k];
                    double ax=x+offdiag*((k>0?p->vectors[j][k-1]:0.0)+
                                         (k+1<n?p->vectors[j][k+1]:0.0));
                    double res=ax-p->result[j]*x;
                    norm2+=x*x;res2+=res*res;ax2+=ax*ax;
                }
                if(sqrt(res2)/fmax(sqrt(ax2)+fabs(p->result[j])*sqrt(norm2),
                                   DBL_MIN)>1e-7){
                    fprintf(stderr,"FAIL independent Toeplitz residual\n");
                    free(p);return 6;
                }
            }
            if(p->restarts>0)challenged++;
            free(p);
        }
        if(challenged==0){
            fprintf(stderr,"FAIL no difficult case exercised a restart\n");
            return 7;
        }
        puts("PASS Toeplitz equal-diagonal restart and analytic spectrum");
    }

    /* Reject impossible workspace dimensions and nonfinite tolerance
     * before initializing any Ritz basis or writing output arrays. */
    {
        Davidson *p=(Davidson*)malloc(sizeof(Davidson));
        if(!p)return 2;
        if(dv_solve(a,4,12,3,3,1e-8,100,p)!=-10 ||
           dv_solve(a,64,24,64,3,1e-8,100,p)!=-10 ||
           dv_solve(a,64,24,3,3,NAN,100,p)!=-10){
            fprintf(stderr,"FAIL Davidson invalid input validation\n");
            free(p);return 8;
        }
        free(p);
        puts("PASS invalid workspace and tolerance rejection");
    }

    /* Exact zero denominators must remain finite and perpendicular
     * to the current Ritz vector after diagonal regularization. */
    {
        double diag[4]={1,1,2,3},x[4]={1,0,0,0};
        double r[4]={0,1,0.1,-0.2},z[4]={0};
        if(dv_correction(4,diag,1.0,x,r,z,NULL)!=0||
           !isfinite(dv_norm(4,z))||fabs(dv_dot(4,x,z))>1e-12)
            return 3;
        puts("PASS regularized zero-diagonal denominator");
    }

    /* A true preconditioned correction differs from subtracting x(x^T t).
     * Check the formula and a nontrivial difference after x-projection. */
    {
        const double q=1.0/sqrt(3.0),h=1.0/sqrt(2.0);
        double diag[4]={2.0,3.0,5.0,7.0};
        double x[4]={q,q,q,0.0},r[4]={h,-h,0.0,0.0};
        double z[4]={0},simple[4]={0},difference[4]={0};
        double numerator=0.0,denominator=0.0,epsilon,proj;
        int used=0,j;
        for(j=0;j<4;j++){
            numerator+=x[j]*r[j]/(diag[j]-1.0);
            denominator+=x[j]*x[j]/(diag[j]-1.0);
            simple[j]=-r[j]/(diag[j]-1.0);
        }
        epsilon=numerator/denominator;
        proj=dv_dot(4,x,simple);
        if(dv_correction(4,diag,1.0,x,r,z,&used)!=0||!used)
            return 9;
        for(j=0;j<4;j++){
            double expected=(-r[j]+epsilon*x[j])/(diag[j]-1.0);
            if(fabs(z[j]-expected)>2e-13)return 10;
            simple[j]-=proj*x[j];
            difference[j]=z[j]-simple[j];
        }
        if(fabs(dv_dot(4,x,z))>2e-13||
           dv_norm(4,difference)<0.02*dv_norm(4,z))
            return 11;
        puts("PASS oblique diagonal JD correction differs from projected Davidson");
    }
    /* For an indefinite M with x^T M^(-1)x exactly zero, use a
     * bounded standard Davidson direction rather than divide by zero. */
    {
        const double q=1.0/sqrt(2.0);
        double diag[4]={2.0,0.0,3.0,4.0};
        double x[4]={q,q,0.0,0.0},r[4]={0.0,0.0,1.0,0.0};
        double z[4]={0};
        int used=1;
        if(dv_correction(4,diag,1.0,x,r,z,&used)!=0||used||
           !isfinite(dv_norm(4,z))||fabs(dv_dot(4,x,z))>2e-13)
            return 12;
        puts("PASS zero-denominator fallback for indefinite diagonal M");
    }
    printf("DAVIDSON SUMMARY %d successes %d failures\n",successful,failures);
    return failures==0&&successful==36?0:1;
}
