#include "lis_adapter.h"
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
    double sol[N], orig_rhs[N], orig_sol[N], scratch[N]; uint32_t perm[N];
    for(unsigned i=0;i<N;i++)perm[i]=i;
    for(unsigned i=N-1;i>0;i--){unsigned j=(unsigned)(fabs(randv())*(i+1));uint32_t t=perm[i];perm[i]=perm[j];perm[j]=t;}
    if(serial%4==0)for(unsigned i=0;i<N;i++)perm[i]=i;
    for(unsigned i=0;i<N;i++)orig_rhs[perm[i]]=rhs[i];
    p56_context *f=NULL;
    int rc=p56_rebuild(&f,N,np,p,kind,d,off,serial%4==0?NULL:perm);
    int ok=rc==P55_OK&&p56_block_count(f)==expected&&p56_size(f)==N;
    if(ok)ok=p56_apply(f,orig_rhs,orig_sol,scratch,N)==P55_OK;
    if(ok)for(unsigned i=0;i<N;i++)sol[i]=orig_sol[perm[i]];
    if(ok) {
        memcpy(orig_sol,orig_rhs,sizeof orig_sol);
        ok=p56_apply(f,orig_sol,orig_sol,scratch,N)==P55_OK;
        for(unsigned i=0;i<N;i++)if(orig_sol[perm[i]]!=sol[i])ok=0;
    }
    (void)ib;
    if(ok){long double rn=0,an=0,xn=0,bn=0;double err=0;
        for(unsigned i=0;i<N;i++){long double ax=0,ar=0;for(unsigned j=0;j<N;j++){ax+=A[i][j]*sol[j];ar+=fabsl(A[i][j]);}rn=fmaxl(rn,fabsl(ax-rhs[i]));an=fmaxl(an,ar);xn=fmaxl(xn,fabsl(sol[i]));bn=fmaxl(bn,fabsl(rhs[i]));err=fmax(err,fabs(sol[i]-x[i]));}
        double res=(double)(rn/(an*xn+bn));worst_err=fmax(worst_err,err);worst_res=fmax(worst_res,res);ok=err<1e-12&&res<1e-13;
    }
    check(ok,"permuted dense reference, out-of-place and in-place");p56_destroy(f);
}
static void edge(void){
    p56_context *f=NULL; uint8_t k[]={2,0}; double d[]={0,0},o[]={2,0},v[]={0};
    uint32_t perm[]={1,0}; p55_panel_input p={0,2,0,NULL,v,1};
    double rhs[]={4,6},out[]={9,9},scratch[2],overlap[]={4,6,99};
    check(p56_rebuild(&f,2,1,&p,k,d,o,perm)==P55_OK,"initial build");
    check(p56_apply(f,rhs,out,scratch,2)==P55_OK&&out[0]==3&&out[1]==2,"block solve");
    check(p56_apply(f,overlap,overlap+1,scratch,2)==P55_OK&&overlap[1]==3&&overlap[2]==2,"partially overlapping rhs/out");
    p56_context *old=f;
    perm[1]=1;
    check(p56_rebuild(&f,2,1,&p,k,d,o,perm)==P55_INVALID&&f==old,"duplicate permutation rollback");
    perm[1]=2;
    check(p56_rebuild(&f,2,1,&p,k,d,o,perm)==P55_INVALID&&f==old,"out-of-range permutation rollback");
    perm[1]=0;d[0]=d[1]=o[0]=1;
    check(p56_rebuild(&f,2,1,&p,k,d,o,perm)==P55_SINGULAR&&f==old,"singular rollback");
    check(p56_apply(f,rhs,out,scratch,2)==P55_OK&&out[0]==3&&out[1]==2,"old factor retained");
    out[0]=out[1]=9;
    check(p56_apply(f,rhs,out,scratch,1)==P55_INVALID&&out[0]==9&&out[1]==9,"short workspace");
    check(p56_apply(f,rhs,out,rhs,2)==P55_INVALID&&rhs[0]==4&&out[0]==9,"scratch aliases rhs");
    check(p56_apply(f,rhs,out,out,2)==P55_INVALID&&out[0]==9,"scratch aliases output");
    double buf[]={9,9,9};
    check(p56_apply(f,rhs,buf,buf+1,2)==P55_INVALID&&buf[0]==9&&buf[1]==9,"partial scratch overlap");
    rhs[0]=NAN;
    check(p56_apply(f,rhs,out,scratch,2)==P55_RANGE&&out[0]==9&&out[1]==9,"nonfinite rhs preserves output");
    rhs[0]=4;k[0]=k[1]=1;d[0]=2;d[1]=3;perm[0]=0;perm[1]=1;
    check(p56_rebuild(&f,2,1,&p,k,d,NULL,perm)==P55_OK&&p56_block_count(f)==0,"replace by scalar factor");
    d[0]=0;v[0]=99;perm[0]=1;
    check(p56_apply(f,rhs,out,scratch,2)==P55_OK&&out[0]==2&&out[1]==2,"inputs copied");
    k[0]=k[1]=1;d[0]=DBL_TRUE_MIN;d[1]=1;v[0]=0;
    check(p56_rebuild(&f,2,1,&p,k,d,NULL,NULL)==P55_OK,"subnormal scalar build");
    out[0]=out[1]=9;
    check(p56_apply(f,rhs,out,scratch,2)==P55_RANGE&&out[0]==9&&out[1]==9,"overflow preserves output");
    check(p56_apply(f,NULL,out,scratch,2)==P55_INVALID,"null rhs");
    check(p56_rebuild(NULL,2,1,&p,k,d,NULL,NULL)==P55_INVALID,"null slot");
    p56_destroy(f);p56_destroy(NULL);
}
int main(void){double scales[]={1e-100,1e-50,1,1e50,1e100};for(int m=0;m<3;m++)for(unsigned s=0;s<5;s++)for(int t=0;t<20;t++)run_case(m,scales[s],t);edge();printf("PASS %u/%u\nmax_abs_solution_error=%.17g\nmax_scaled_backward_error=%.17g\n",passed,tests,worst_err,worst_res);return passed!=tests;}
