#include "packed_ldlt.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 24
static unsigned long long rng=0x53ab98d;
static double randv(void){rng=rng*6364136223846793005ULL+1;return ((double)(rng>>32)/4294967296.0)*2-1;}
static unsigned tests,passed;static double worst_err,worst_res;
static void check(int ok,const char *name){tests++;if(ok)passed++;else fprintf(stderr,"FAIL %s\n",name);}
static void run_case(int mode,double scale,int serial){
    double L[N][N]={{0}},d[N],off[N]={0},x[N],rhs[N],values[N][N*N];
    long double A[N][N]={{0}}; uint8_t kind[N];uint32_t rows[N][N];p55_panel_input p[N];
    size_t np=0;unsigned expected=0,ib=0;uint32_t first=0;
    for(unsigned i=0;i<N;i++){L[i][i]=1;x[i]=randv();kind[i]=1;d[i]=(1+fabs(randv()))*scale*(i%2?1:-1);}
    while(first<N){uint32_t w=1+(unsigned)(fabs(randv())*5);if(w>N-first)w=N-first;uint32_t end=first+w,nr=0;
        for(uint32_t r=end;r<N;r++)if(fabs(randv())>.35)rows[np][nr++]=r;
        for(uint32_t j=first;j+1<end;j++)if(mode && (mode==2 || fabs(randv())>.45)){
            kind[j]=2;kind[j+1]=0;d[j]=(serial%2?0:.1)*scale;d[j+1]=-.2*scale;off[j]=(1+fabs(randv()))*scale;expected++;j++;}
        size_t z=0;for(uint32_t j=first;j<end;j++){
            for(uint32_t r=j+1;r<end;r++){double v=(kind[j]==2&&r==j+1)?0:randv()*.15;L[r][j]=v;values[np][z++]=v;}
            for(uint32_t r=0;r<nr;r++){double v=randv()*.15;L[rows[np][r]][j]=v;values[np][z++]=v;}}
        p[np]=(p55_panel_input){first,w,nr,rows[np],values[np],z};np++;ib+=2*nr;first=end;
    }
    /* Independent dense reconstruction of A; no factor solve used to form RHS. */
    for(unsigned i=0;i<N;i++)for(unsigned j=0;j<N;j++)for(unsigned k=0;k<N;k++){
        A[i][j]+=(long double)L[i][k]*d[k]*L[j][k];
        if(kind[k]==2)A[i][j]+=(long double)off[k]*( (long double)L[i][k]*L[j][k+1]+(long double)L[i][k+1]*L[j][k]);}
    for(unsigned i=0;i<N;i++){long double sum=0;for(unsigned j=0;j<N;j++)sum+=A[i][j]*x[j];rhs[i]=(double)sum;}
    double sol[N];memcpy(sol,rhs,sizeof sol);p55_factor *f=NULL;
    int rc=p55_pack(&f,N,np,p,kind,d,off);int ok=rc==P55_OK&&p55_block_count(f)==expected&&p55_index_bytes(f)==ib;
    if(ok)ok=p55_solve(f,sol)==P55_OK;
    if(ok){long double rn=0,an=0,xn=0,bn=0;double err=0;
        for(unsigned i=0;i<N;i++){long double ax=0,ar=0;for(unsigned j=0;j<N;j++){ax+=A[i][j]*sol[j];ar+=fabsl(A[i][j]);}rn=fmaxl(rn,fabsl(ax-rhs[i]));an=fmaxl(an,ar);xn=fmaxl(xn,fabsl(sol[i]));bn=fmaxl(bn,fabsl(rhs[i]));err=fmax(err,fabs(sol[i]-x[i]));}
        double res=(double)(rn/(an*xn+bn));worst_err=fmax(worst_err,err);worst_res=fmax(worst_res,res);ok=err<1e-12&&res<1e-13;
    }
    check(ok,"dense reference");p55_free(f);
}
static void edge(void){
    p55_factor *f=NULL;uint8_t k[]={2,0};double d[]={0,0},o[]={2,0},v[]={0},x[]={4,6};
    p55_panel_input p={0,2,0,NULL,v,1};
    check(p55_pack(&f,2,1,&p,k,d,o)==0&&p55_solve(f,x)==0&&x[0]==3&&x[1]==2,"zero diagonal");p55_free(f);
    d[0]=d[1]=o[0]=1;check(p55_pack(&f,2,1,&p,k,d,o)==P55_SINGULAR&&f==NULL,"singular block");
    d[0]=d[1]=0;o[0]=2;v[0]=1;check(p55_pack(&f,2,1,&p,k,d,o)==P55_INVALID&&f==NULL,"nonzero L inside block");v[0]=0;
    p55_panel_input split[]={{0,1,0,NULL,NULL,0},{1,1,0,NULL,NULL,0}};
    check(p55_pack(&f,2,2,split,k,d,o)==P55_INVALID&&f==NULL,"block crosses panel");
    k[0]=1;k[1]=1;check(p55_pack(&f,2,1,&p,k,d,o)==P55_SINGULAR,"zero scalar");
    d[0]=d[1]=1;p.nvalues=0;check(p55_pack(&f,2,1,&p,k,d,o)==P55_INVALID,"bad value count");p.nvalues=1;
    v[0]=NAN;check(p55_pack(&f,2,1,&p,k,d,o)==P55_INVALID,"nan factor");v[0]=0;
    check(p55_pack(&f,2,1,&p,k,d,o)==0,"pack scalar");x[0]=NAN;check(p55_solve(f,x)==P55_RANGE,"nan rhs");p55_free(f);
    k[0]=2;k[1]=0;d[0]=d[1]=0;o[0]=DBL_TRUE_MIN;x[0]=1;x[1]=1;
    check(p55_pack(&f,2,1,&p,k,d,o)==0&&p55_solve(f,x)==P55_RANGE,"overflow block solve");p55_free(f);
    o[0]=DBL_TRUE_MIN;x[0]=DBL_TRUE_MIN;x[1]=-DBL_TRUE_MIN;
    check(p55_pack(&f,2,1,&p,k,d,o)==0&&p55_solve(f,x)==0&&x[0]==-1&&x[1]==1,"subnormal block");p55_free(f);
}
static void compression(void){
    const uint32_t n=65538; p55_panel_input *p=calloc(n,sizeof *p);uint8_t *k=malloc(n);double *d=malloc(n*sizeof *d),*x=malloc(n*sizeof *x);uint32_t row=n-1;double val=.25;p55_factor *f=NULL;
    if(!p||!k||!d||!x){check(0,"compression alloc");goto done;}
    for(uint32_t i=0;i<n;i++){p[i].first=i;p[i].width=1;k[i]=1;d[i]=1;x[i]=1;}
    p[0].nrows=1;p[0].rows=&row;p[0].values=&val;p[0].nvalues=1;
    x[0]=1.25;x[n-1]=1.3125;
    check(p55_pack(&f,n,n,p,k,d,NULL)==0&&p55_index_bytes(f)==4&&p55_solve(f,x)==0&&x[0]==1&&x[n-1]==1,"32 bit fallback");p55_free(f);f=NULL;
    row=n-2;for(uint32_t i=0;i<n;i++)x[i]=1;x[0]=1.25;x[row]=1.3125;
    check(p55_pack(&f,n,n,p,k,d,NULL)==0&&p55_index_bytes(f)==2&&p55_solve(f,x)==0&&x[0]==1&&x[row]==1,"16 bit boundary");p55_free(f);
    row=0;check(p55_pack(&f,n,n,p,k,d,NULL)==P55_INVALID,"invalid external row");
    done:free(p);free(k);free(d);free(x);
}
int main(void){double scales[]={1e-100,1e-50,1,1e50,1e100};for(int m=0;m<3;m++)for(unsigned s=0;s<5;s++)for(int t=0;t<20;t++)run_case(m,scales[s],t);edge();compression();printf("PASS %u/%u\nmax_abs_solution_error=%.17g\nmax_scaled_backward_error=%.17g\n",passed,tests,worst_err,worst_res);return passed!=tests;}
