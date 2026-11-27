#define _POSIX_C_SOURCE 200809L
#include "factor_cache.h"
#include <unistd.h>
#include <sys/stat.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 24
static unsigned long long rng=0x53ab98d;
static double randv(void){rng=rng*6364136223846793005ULL+1;return ((double)(rng>>32)/4294967296.0)*2-1;}
static char path[256]; static uint8_t key[32]={0x42,0x57};
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
    if(ok){
        double before[N];memcpy(before,orig_sol,sizeof before);
        p56_context *loaded=NULL;
        ok=p57_save(f,path,key)==P55_OK && p57_load(&loaded,path,key,1000000,N)==P55_OK;
        if(ok)ok=p56_size(loaded)==N&&p56_block_count(loaded)==expected&&p56_apply(loaded,orig_rhs,orig_sol,scratch,N)==P55_OK&&memcmp(before,orig_sol,sizeof before)==0;
        p56_destroy(loaded);
    }
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

static int write_bytes(const unsigned char *data,size_t n){FILE *fp=fopen(path,"wb");if(!fp)return 0;int ok=fwrite(data,1,n,fp)==n;return fclose(fp)==0&&ok;}
static void set32(unsigned char *p,uint32_t v){for(unsigned i=0;i<4;i++)p[i]=(unsigned char)(v>>(8*i));}
static void checksum(unsigned char *p,size_t n){uint32_t c=UINT32_MAX;for(size_t i=0;i<n-4;i++){c^=p[i];for(int j=0;j<8;j++)c=(c>>1)^((c&1)?UINT32_C(0xedb88320):0);}set32(p+n-4,c^UINT32_MAX);}
static void edge(void){
    p56_context *f=NULL; uint8_t k[]={2,0}; double d[]={0,0},o[]={2,0},v[]={0};uint32_t perm[]={1,0};p55_panel_input p={0,2,0,NULL,v,1};
    check(p56_rebuild(&f,2,1,&p,k,d,o,perm)==P55_OK,"initial");
    check(p57_save(f,path,key)==P55_OK,"save");
    unsigned char data[256],mut[257];FILE *fp=fopen(path,"rb");if(!fp){check(0,"open");p56_destroy(f);return;}size_t n=fread(data,1,sizeof data,fp);fclose(fp);
    check(n==116,"canonical bytes"); if(n!=116){p56_destroy(f);return;}
    p56_context *old=f;uint8_t other[32]={0};
    check(p57_load(&f,path,other,10000,2)==P57_MISS&&f==old,"key mismatch");
    check(p57_load(&f,path,key,n-1,2)==P57_LIMIT&&f==old,"byte limit");
    check(p57_load(&f,path,key,n,1)==P57_LIMIT&&f==old,"dimension limit");
    for(size_t i=0;i<n;i++){
        memcpy(mut,data,n);mut[i]^=1;
        check(write_bytes(mut,n)&&p57_load(&f,path,key,10000,2)!=P55_OK&&f==old,"single-bit corruption rollback");
    }
    for(size_t i=0;i<n;i++)check(write_bytes(data,i)&&p57_load(&f,path,key,10000,2)!=P55_OK&&f==old,"truncated rollback");
    memcpy(mut,data,n);mut[n]=0;check(write_bytes(mut,n+1)&&p57_load(&f,path,key,10000,2)!=P55_OK&&f==old,"trailing byte");
    /* Valid checksum with malformed fields: independent semantic validation. */
    const struct {size_t pos;uint32_t val;} cases[]={
      {40,0},{44,0},{48,2},{52,2},{72,1},{84,2},{88,1},{92,0},{92,1},{96,1},{100,3},{100,4},{108,0x7ff80000},{80,0x7ff80000},{80,0}};
    for(size_t i=0;i<sizeof cases/sizeof *cases;i++){
        memcpy(mut,data,n);set32(mut+cases[i].pos,cases[i].val);checksum(mut,n);
        check(write_bytes(mut,n)&&p57_load(&f,path,key,10000,2)!=P55_OK&&f==old,"valid CRC invalid semantic field");
    }
    check(write_bytes(data,n)&&p57_load(&f,path,key,n,2)==P55_OK,"successful replacement");
    double rhs[]={4,6},out[2],scratch[2];check(p56_apply(f,rhs,out,scratch,2)==P55_OK&&out[0]==3&&out[1]==2,"restored block");
    check(p57_save(f,"/nonexistent-p57-dir/test",key)==P57_IO,"save IO failure");
    unlink(path);check(p57_load(&f,path,key,10000,2)==P57_MISS,"missing file");
    check(mkfifo(path,0600)==0&&p57_load(&f,path,key,10000,2)==P57_IO,"FIFO rejected without blocking");unlink(path);
    p56_destroy(f);
}
static void wide(void){
    const uint32_t n=65538; p55_panel_input *p=calloc(n,sizeof *p);uint8_t *k=malloc(n);double *d=malloc((size_t)n*8),*rhs=malloc((size_t)n*8),*out=malloc((size_t)n*8),*tmp=malloc((size_t)n*8),*before=malloc((size_t)n*8);
    if(!p||!k||!d||!rhs||!out||!tmp||!before){check(0,"wide alloc");goto done;}
    uint32_t rows[]={65536,65537};double v[]={.25,-.125};
    for(uint32_t i=0;i<n;i++){p[i]=(p55_panel_input){i,1,0,NULL,NULL,0};k[i]=1;d[i]=2;rhs[i]=1;}
    p[0]=(p55_panel_input){0,1,2,rows,v,2};p[1]=(p55_panel_input){1,1,1,rows+1,v,1}; /* offset 65535 */
    p56_context *f=NULL,*g=NULL;
    int ok=p56_rebuild(&f,n,n,p,k,d,NULL,NULL)==P55_OK&&p56_apply(f,rhs,before,tmp,n)==P55_OK&&p57_save(f,path,key)==P55_OK&&p57_load(&g,path,key,4000000,n)==P55_OK&&p56_apply(g,rhs,out,tmp,n)==P55_OK&&memcmp(before,out,(size_t)n*8)==0;
    check(ok,"wide 16/32bit indexes roundtrip");p56_destroy(f);p56_destroy(g);
done:free(p);free(k);free(d);free(rhs);free(out);free(tmp);free(before);
}
int main(void){char dir[]="/tmp/lis57.XXXXXX";if(!mkdtemp(dir))return 2;snprintf(path,sizeof path,"%s/factor",dir);
    double scales[]={1e-100,1e-50,1,1e50,1e100};for(int m=0;m<3;m++)for(unsigned s=0;s<5;s++)for(int t=0;t<20;t++)run_case(m,scales[s],t);edge();wide();unlink(path);rmdir(dir);
    printf("PASS %u/%u\nmax_abs_solution_error=%.17g\nmax_scaled_backward_error=%.17g\n",passed,tests,worst_err,worst_res);return passed!=tests;}
