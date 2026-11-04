/* Independent C99 proof of concept for implicitly restarted Lanczos.
 * A dense QR similarity explicitly realizes the implicit spectral shifts.
 * This is intentionally NOT a new public LIS solver or production kernel.
 * No Spectra source code is copied; this is an independent implementation.
 */
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXN 80
#define MAXM 32

typedef struct {
    int n, m, k, nev, op_calls, qr_steps, restarts, checks;
    double V[MAXM+1][MAXN], T[MAXM][MAXM], extra;
} Proto;

static double scalar(int n,const double*x,const double*y) {
    int i; double a=0; for(i=0;i<n;i++)a+=x[i]*y[i];return a;
}
static double norm(int n,const double*x){return sqrt(fmax(0,scalar(n,x,x)));}
static void matvec(Proto*p,const double A[MAXN][MAXN],const double*x,double*y){
    int i,j; for(i=0;i<p->n;i++){double a=0;for(j=0;j<p->n;j++)a+=A[i][j]*x[j];y[i]=a;}p->op_calls++;
}
static int eigen_symmetric(int m,const double X[MAXM][MAXM],double eval[MAXM],double U[MAXM][MAXM]){
    double W[MAXM][MAXM], max=0;int i,j,t,round;
    memset(W,0,sizeof W); memset(U,0,MAXM*MAXM*sizeof(double));
    for(i=0;i<m;i++){U[i][i]=1;for(j=0;j<m;j++)W[i][j]=X[i][j];}
    for(round=0;round<60*m*m;round++){
        int p=0,q=0;max=0;
        for(i=0;i<m;i++)for(j=i+1;j<m;j++)if(fabs(W[i][j])>max){max=fabs(W[i][j]);p=i;q=j;}
        if(max<3e-14)break;
        double theta=0.5*atan2(2*W[p][q],W[q][q]-W[p][p]);
        double c=cos(theta),s=sin(theta), app=W[p][p],aqq=W[q][q],apq=W[p][q];
        W[p][p]=c*c*app-2*c*s*apq+s*s*aqq;
        W[q][q]=s*s*app+2*c*s*apq+c*c*aqq;
        W[p][q]=W[q][p]=0;
        for(t=0;t<m;t++)if(t!=p&&t!=q){double x=W[t][p],y=W[t][q];W[t][p]=W[p][t]=c*x-s*y;W[t][q]=W[q][t]=s*x+c*y;}
        for(t=0;t<m;t++){double x=U[t][p],y=U[t][q];U[t][p]=c*x-s*y;U[t][q]=s*x+c*y;}
    }
    if(round==60*m*m) return -1;
    for(i=0;i<m;i++) eval[i]=W[i][i];
    return 0;
}
static int shifted_qr(Proto*p,int m,double T[MAXM][MAXM],double shift,double Qtotal[MAXM][MAXM]){
    double R[MAXM][MAXM],Q[MAXM][MAXM],temp[MAXM][MAXM],next[MAXM][MAXM],newQ[MAXM][MAXM];
    int i,j,l;
    memset(R,0,sizeof R);memset(Q,0,sizeof Q);
    for(i=0;i<m;i++){Q[i][i]=1;for(j=0;j<m;j++)R[i][j]=T[i][j]-(i==j?shift:0);}
    /* QR by Givens rotations, accumulating the orthogonal right factor. */
    for(j=0;j<m-1;j++)for(i=m-1;i>j;i--){
        double a=R[i-1][j],b=R[i][j],h=hypot(a,b);
        if(h<DBL_MIN)continue;
        double c=a/h,s=b/h;
        for(l=j;l<m;l++){double x=R[i-1][l],y=R[i][l];R[i-1][l]=c*x+s*y;R[i][l]=-s*x+c*y;}
        for(l=0;l<m;l++){double x=Q[l][i-1],y=Q[l][i];Q[l][i-1]=c*x+s*y;Q[l][i]=-s*x+c*y;}
    }
    memset(temp,0,sizeof temp);memset(next,0,sizeof next);memset(newQ,0,sizeof newQ);
    for(i=0;i<m;i++)for(j=0;j<m;j++)for(l=0;l<m;l++)temp[i][j]+=T[i][l]*Q[l][j];
    for(i=0;i<m;i++)for(j=0;j<m;j++)for(l=0;l<m;l++)next[i][j]+=Q[l][i]*temp[l][j];
    for(i=0;i<m;i++)for(j=0;j<m;j++)for(l=0;l<m;l++)newQ[i][j]+=Qtotal[i][l]*Q[l][j];
    /* A symmetric tridiagonal QR similarity remains tridiagonal in
     * exact arithmetic. Roundoff-only fill must be removed AFTER
     * EACH shift, otherwise the next Givens QR can rotate nearly-zero
     * non-Hessenberg entries into large offband elements. Refuse to
     * hide significant fill: this is a structural invariant, not
     * an unconditional band truncation. */
    {
        double scale=1.0;
        for(i=0;i<m;i++)for(j=0;j<m;j++)
            if(fabs(next[i][j])>scale)scale=fabs(next[i][j]);
        for(i=0;i<m;i++)for(j=0;j<m;j++){
            double v=0.5*(next[i][j]+next[j][i]);
            if(abs(i-j)>1){
                if(fabs(v)>1e-10*scale) return -1;
                v=0.0;
            }
            T[i][j]=v;
            Qtotal[i][j]=newQ[i][j];
        }
    }
    p->qr_steps++;
    return 0;
}
static int extend(Proto*p,const double A[MAXN][MAXN],int start){
    int j,i,t;double w[MAXN];
    for(j=start;j<p->m;j++){
        matvec(p,A,p->V[j],w);
        if(j>0)for(i=0;i<p->n;i++)w[i]-=p->T[j-1][j]*p->V[j-1][i];
        p->T[j][j]=scalar(p->n,p->V[j],w);
        for(i=0;i<p->n;i++)w[i]-=p->T[j][j]*p->V[j][i];
        /* Two-pass full reorthogonalization on clustered spectra. */
        for(t=0;t<2;t++)for(i=0;i<=j;i++){
            double coeff=scalar(p->n,p->V[i],w);
            int row;for(row=0;row<p->n;row++)w[row]-=coeff*p->V[i][row];
        }
        double beta=norm(p->n,w);
        if(!isfinite(beta)||beta < 1e-14){fprintf(stderr,"Lanczos breakdown step=%d beta=%g\n",j,beta);return -1;}
        if(j+1<p->m)p->T[j][j+1]=p->T[j+1][j]=beta;
        else p->extra=beta;
        for(i=0;i<p->n;i++)p->V[j+1][i]=w[i]/beta;
    }
    return 0;
}
static int restart(Proto*p,const double eval[MAXM]){
    double mat[MAXM][MAXM],Qtotal[MAXM][MAXM],tmp[MAXM+1][MAXN],f[MAXN], normf;
    int i,j,t,n=p->n,m=p->m,k=p->k,shifts=m-k;
    memset(mat,0,sizeof mat);memset(Qtotal,0,sizeof Qtotal);
    for(i=0;i<m;i++){Qtotal[i][i]=1;for(j=0;j<m;j++)mat[i][j]=p->T[i][j];}
    /* Shift by unwanted Ritz values (the smallest m-k eigenvalues). */
    for(i=0;i<shifts;i++) if(shifted_qr(p,m,mat,eval[i],Qtotal))
        return -8;
    for(i=0;i<k+1;i++)for(t=0;t<n;t++){
        double z=0;for(j=0;j<m;j++)z+=Qtotal[j][i]*p->V[j][t];
        tmp[i][t]=z;
    }
    for(t=0;t<n;t++)
        f[t]=mat[k][k-1]*tmp[k][t]+p->extra*Qtotal[m-1][k-1]*p->V[m][t];
    normf=norm(n,f);
    if(!isfinite(normf)||normf<1e-14) return -1;
    /* The transformed residual must affect only the last retained column. */
    for(i=0;i<k-1;i++)if(fabs(Qtotal[m-1][i])>1e-8){
        fprintf(stderr,"restart residual leaks into column %d: %g\n",i,Qtotal[m-1][i]);return -2;
    }
    /* QR similarity should preserve tridiagonality (up to rounding).
     * The rotated retained Lanczos basis must remain orthonormal. */
    for(i=0;i<m;i++)for(j=0;j<m;j++)
        if(abs(i-j)>1 && fabs(mat[i][j]) > 3e-9)
            return -3;
    for(i=0;i<k;i++)for(j=0;j<=i;j++){
        double dot=scalar(n,tmp[i],tmp[j]);
        if(fabs(dot-(i==j ? 1.0 : 0.0))>2e-8)
            return -4;
    }
    for(i=0;i<k;i++)if(fabs(scalar(n,tmp[i],f))/normf>2e-8)
        return -5;
    for(i=0;i<k;i++)for(t=0;t<n;t++)p->V[i][t]=tmp[i][t];
    for(t=0;t<n;t++)p->V[k][t]=f[t]/normf;
    memset(p->T,0,sizeof p->T);
    for(i=0;i<k;i++){
        p->T[i][i]=mat[i][i];
        if(i+1<k)p->T[i+1][i]=p->T[i][i+1]=mat[i+1][i];
    }
    p->T[k][k-1]=p->T[k-1][k]=normf;
    return 0;
}
static int solve(const double A[MAXN][MAXN],int n,int m,int k,int nev,double tol,int maxrestart,double *lambdas,double *resids,int *ops,int *restarts){
    Proto *p=calloc(1,sizeof *p);
    double eig[MAXM],U[MAXM][MAXM],sorted[MAXM],x[MAXN],y[MAXN];
    double accepted[MAXM][MAXN];
    int ix[MAXM],i,j,l,r,done=0,verified;
    if(!p)return -10;
    if(n>MAXN||m>MAXM||m>=n||k>=m||nev>=k||k<1){free(p);return -10;}
    p->n=n;p->m=m;p->k=k;p->nev=nev;
    for(i=0;i<n;i++)p->V[0][i]=sin((double)(i+1)*0.771)+cos((double)(i+1)*1.109);
    double len=norm(n,p->V[0]);for(i=0;i<n;i++)p->V[0][i]/=len;
    if(extend(p,A,0)){free(p);return -6;}
    for(r=0;r<=maxrestart;r++){
        if(eigen_symmetric(m,p->T,eig,U)){done=-2;break;}
        for(i=0;i<m;i++)ix[i]=i;
        for(i=0;i<m;i++)for(j=i+1;j<m;j++)if(eig[ix[j]]<eig[ix[i]]){l=ix[j];ix[j]=ix[i];ix[i]=l;}
        verified=0;
        for(i=0;i<nev;i++){
            int id=ix[m-1-i];
            double estimate=fabs(p->extra*U[m-1][id]),residue,norma;
            /* This is a valid Lanczos estimate; the real physical residual
             * must still be computed before a Ritz value is accepted. */
            if(estimate>tol*(1+fabs(eig[id]))*0.5) continue;
            for(j=0;j<n;j++){double a=0;for(l=0;l<m;l++)a+=U[l][id]*p->V[l][j];x[j]=a;}
            matvec(p,A,x,y);p->checks++;
            norma=norm(n,y)+fabs(eig[id])*norm(n,x);
            for(j=0;j<n;j++)y[j]-=eig[id]*x[j];
            residue=norm(n,y)/fmax(norma,DBL_MIN);
            lambdas[i]=eig[id];resids[i]=residue;
            if(residue<tol){
                memcpy(accepted[i],x,(size_t)n*sizeof(double));
                verified++;
            }
        }
        if(verified==nev){
            for(i=0;i<nev;i++)for(j=0;j<i;j++)
                if(fabs(scalar(n,accepted[i],accepted[j]))>2e-8)
                    verified=-1;
            if(verified==nev){done=1;break;}
            done=-7;break;
        }
        if(r==maxrestart)break;
        for(i=0;i<m;i++)sorted[i]=eig[ix[i]];
        int ret=restart(p,sorted);
        if(ret){done=ret;break;}
        if(extend(p,A,k)){done=-6;break;}
        p->restarts++;
    }
    *ops=p->op_calls;*restarts=p->restarts;
    free(p);return done;
}
static void make_matrix(double A[MAXN][MAXN],int n,int kind){
    int i;memset(A,0,sizeof(double)*MAXN*MAXN);
    for(i=0;i<n;i++){
        if(kind==0)A[i][i]=1+0.08*i;
        else if(kind==1)A[i][i]=1+0.0001*i+0.06*(i/6);
        else A[i][i]=1+0.00001*i+0.11*(i/12);
        if(i+1<n)A[i][i+1]=A[i+1][i]=0.0025*(1+(i%5));
        if(i+3<n)A[i][i+3]=A[i+3][i]=0.0005;
    }
}
int main(void){
    static double A[MAXN][MAXN];
    const double expected[3][3] = {
        {6.04069946274,5.95961541770,5.87976745035},
        {1.61549914705,1.60937686515,1.60522013312},
        {1.55986684425,1.55315504490,1.54939876293}
    };
    int kind,n=64,m,k,nev=3,ok=0,failures=0;
    for(m=16;m<=28;m+=4)for(k=5;k<=11;k+=3)for(kind=0;kind<3;kind++){
        double evals[3]={0},res[3]={0};int ops=0,restarts=0,i;
        make_matrix(A,n,kind);
        int code=solve(A,n,m,k,nev,2e-8,150,evals,res,&ops,&restarts);
        printf("CASE %d m=%d k=%d rc=%d ops=%d restarts=%d",
               kind,m,k,code,ops,restarts);
        for(i=0;i<nev;i++)printf(" val%d=%.12g res%d=%.3g",i,evals[i],i,res[i]);
        puts("");
        if(code!=1){fprintf(stderr,"FAIL IRL case %d m=%d k=%d\n",kind,m,k);failures++;continue;}
        for(i=0;i<nev;i++)if(!isfinite(evals[i])||!isfinite(res[i])||
           fabs(evals[i]-expected[kind][i])>2e-6||res[i]>2e-8){
            fprintf(stderr,"FAIL eigenpair kind=%d m=%d k=%d eigen=%d\n",kind,m,k,i);
            failures++;
        }
        ok++;
    }
    printf("IRL SUMMARY %d successful cases; %d failures\n",ok,failures);
    return failures==0&&ok==36?0:1;
}
