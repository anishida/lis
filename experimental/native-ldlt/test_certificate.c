#include "mode_certificate.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define N 12
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c);exit(1);}}while(0)
static unsigned tests=0;static long double max_res=0;static double max_lambda=0;
/* Construct K = Q diag(g_i*lambda_i) Q^T, G = Q diag(g_i) Q^T.
   Q is composed of exact +/- 1/sqrt(2) pair rotations, which gives
   known generalized eigenvalues 1..12 and non-diagonal G and K. */
typedef struct {p69_csr a; size_t p[N+1];uint32_t j[N*N];double v[N*N];} mat;
static void make(mat *m,const double *a,p69_storage st,unsigned base){size_t z=0;for(unsigned i=0;i<N;i++){m->p[i]=z+base;for(unsigned j=0;j<N;j++)if(a[i*N+j]!=0&&(st==P69_FULL||(st==P69_LOWER&&j<=i)||(st==P69_UPPER&&j>=i))){m->j[z]=j+base;m->v[z++]=a[i*N+j];}}m->p[N]=z+base;m->a=(p69_csr){N,z,m->p,m->j,m->v,base,st};}
static void mv(const double *a,const double *x,double *y){for(unsigned i=0;i<N;i++){long double t=0;for(unsigned j=0;j<N;j++)t+=(long double)a[i*N+j]*x[j];y[i]=(double)t;}}
static long double dot(const double*x,const double*y){long double t=0;for(unsigned i=0;i<N;i++)t+=(long double)x[i]*y[i];return t;}
static void pencil(double *k,double*g,double scale){for(unsigned i=0;i<N;i+=2){double a=(i+1)*scale,b=(i+2)*scale,c=a*(i+1),d=b*(i+2);g[i*N+i]=g[(i+1)*N+i+1]=(a+b)/2;g[i*N+i+1]=g[(i+1)*N+i]=(a-b)/2;k[i*N+i]=k[(i+1)*N+i+1]=(c+d)/2;k[i*N+i+1]=k[(i+1)*N+i]=(c-d)/2;}}
static void analytical(void){for(unsigned t=0;t<36;t++){double k[N*N]={0},g[N*N]={0};pencil(k,g,t%3==0?1e-100:t%3==1?1e100:1);mat a,b;make(&a,k,(p69_storage)(t%3),t%2);make(&b,g,(p69_storage)((t/3)%3),(t/2)%2);for(unsigned j=0;j<N;j++){double x[N]={0};double s=t%2?1e-200:1e200;x[j&~1u]=s;x[(j&~1u)+1]=j%2?-s:s;p70_certificate r;CHECK(!p70_certify(&a.a,&b.a,j+1,x,N,1e-12L,4096,&r)&&r.accepted);if(r.relative_residual>max_res)max_res=r.relative_residual;tests++;CHECK(!p70_certify(&a.a,&b.a,j+1.5,x,N,1e-12L,4096,&r)&&!r.accepted);tests++;}}}
/* Shift-invert iteration with G-inner-product deflation, only for this SPD-G
   test pencil. Not a production eigensolver. */
static void modes(double sigma){double k[N*N]={0},g[N*N]={0};pencil(k,g,1);mat a,b;make(&a,k,P69_LOWER,1);make(&b,g,P69_UPPER,0);p69_options o;p69_defaults(&o);p56_context *f=NULL;CHECK(!p69_factorize(&f,&a.a,&b.a,sigma,&o,NULL));double basis[6][N]={{0}},gb[6][N]={{0}};unsigned used[N]={0};for(unsigned m=0;m<6;m++){double x[N],rhs[N],y[N],w[N];for(unsigned i=0;i<N;i++)x[i]=1.+((i*7+m*3)%13)/13.;for(unsigned it=0;it<1200;it++){mv(g,x,rhs);CHECK(!p56_apply(f,rhs,y,w,N));for(unsigned repeat=0;repeat<2;repeat++)for(unsigned j=0;j<m;j++){long double p=dot(gb[j],y);for(unsigned i=0;i<N;i++)y[i]-=(double)(p*basis[j][i]);}mv(g,y,rhs);long double norm=sqrtl(dot(y,rhs));CHECK(norm>0&&isfinite(norm));for(unsigned i=0;i<N;i++)x[i]=(double)(y[i]/norm);}
  mv(g,x,gb[m]);mv(k,x,rhs);double lam=(double)(dot(x,rhs)/dot(x,gb[m]));unsigned expected=0;double best=INFINITY;for(unsigned j=0;j<N;j++)if(!used[j]&&fabs(j+1-sigma)<best){best=fabs(j+1-sigma);expected=j;}used[expected]=1;double err=fabs(lam-(expected+1));CHECK(err<1e-10);if(err>max_lambda)max_lambda=err;for(unsigned i=0;i<N;i++)basis[m][i]=x[i];p70_certificate r;CHECK(!p70_certify(&a.a,&b.a,lam,x,N,1e-10L,4096,&r)&&r.accepted);if(r.relative_residual>max_res)max_res=r.relative_residual;printf("sigma=%.2f mode=%u lambda=%.15g original_residual=%.3Le PASS\n",sigma,m+1,lam,r.relative_residual);tests++;}p56_destroy(f);}
static void invalid(void){double k[N*N]={0},g[N*N]={0},x[N]={0};pencil(k,g,1);mat a,b;make(&a,k,P69_FULL,0);make(&b,g,P69_FULL,0);p70_certificate r={42,42};
#define BAD(e) do{CHECK((e)!=P55_OK&&r.relative_residual==42&&r.accepted==42);tests++;}while(0)
 BAD(p70_certify(&a.a,&b.a,1,x,N,1e-10L,4096,&r));x[0]=1;BAD(p70_certify(&a.a,&b.a,NAN,x,N,1e-10L,4096,&r));BAD(p70_certify(&a.a,&b.a,1,x,N-1,1e-10L,4096,&r));BAD(p70_certify(&a.a,&b.a,1,x,N,0,4096,&r));BAD(p70_certify(&a.a,&b.a,1,x,N,1e-10L,1,&r));x[0]=INFINITY;BAD(p70_certify(&a.a,&b.a,1,x,N,1e-10L,4096,&r));x[0]=1;a.v[1]+=1;BAD(p70_certify(&a.a,&b.a,1,x,N,1e-10L,4096,&r));a.v[1]-=1;double zero[N*N]={0};make(&a,zero,P69_LOWER,0);make(&b,zero,P69_LOWER,0);BAD(p70_certify(&a.a,&b.a,1,x,N,1e-10L,4096,&r));
#undef BAD
}
int main(void){analytical();modes(.25);modes(5.3);invalid();printf("certificate %u/%u PASS max_original_residual=%.3Le max_eigenvalue_error=%.3e\n",tests,tests,max_res,max_lambda);return 0;}
