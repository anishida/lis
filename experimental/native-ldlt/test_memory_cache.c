#define _POSIX_C_SOURCE 200809L
#include "parallel_ldlt.h"
#include "factor_cache.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#define CHECK(c) do {if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);exit(1);}}while(0)
int scalar_save(const p56_context *,const char *,const uint8_t[32]);
int scalar_load(p56_context **,const char *,const uint8_t[32],uint64_t,uint32_t);
static unsigned tests;
static void same(p56_context *a,p56_context *b,uint32_t n){
 double *r=malloc((size_t)n*8),*x=malloc((size_t)n*8),*y=malloc((size_t)n*8),*w=malloc((size_t)n*8);CHECK(r&&x&&y&&w);
 for(uint32_t i=0;i<n;i++)r[i]=((int)(i%23)-11)/13.;
 CHECK(!p56_apply(a,r,x,w,n)&&!p56_apply(b,r,y,w,n)&&!memcmp(x,y,(size_t)n*8));free(r);free(x);free(y);free(w);
}
static void chain(uint32_t n,uint32_t width){
 size_t *cp=malloc(((size_t)n+1)*sizeof *cp);uint32_t *ri=malloc((size_t)n*2*sizeof *ri);double *a=malloc((size_t)n*2*sizeof *a);CHECK(cp&&ri&&a);size_t z=0;
 for(uint32_t i=0;i<n;i++){cp[i]=z;ri[z]=i;a[z++]=5;if(i+1<n){ri[z]=i+1;a[z++]=-1;}}cp[n]=z;
 p56_context *x=NULL,*y=NULL;p67_stats s,t;p62_limits sl={128u*1024*1024,100000000};p61_limits nl={64,128u*1024*1024,10000000,1000000,0};
 CHECK(!p67_factorize(&x,n,z,cp,ri,a,&sl,&nl,width,1,&s));nl.max_work_bytes=s.numeric.peak_work_bytes;
 CHECK(!p67_factorize(&y,n,z,cp,ri,a,&sl,&nl,width,4,&t));CHECK(t.threads_used==1&&t.max_ready==1&&t.numeric.peak_work_bytes==s.numeric.peak_work_bytes);same(x,y,n);tests++;
 printf("chain n=%u width=%u one=%zu requested_four=%zu levels=%zu\n",n,width,s.numeric.peak_work_bytes,t.numeric.peak_work_bytes,t.levels);
 p56_context *old=y;nl.max_work_bytes--;CHECK(p67_factorize(&y,n,z,cp,ri,a,&sl,&nl,width,4,&t)==P60_LIMIT&&y==old);tests++;
 p56_destroy(x);p56_destroy(y);free(cp);free(ri);free(a);
}
static int files_same(const char *a,const char *b){FILE *x=fopen(a,"rb"),*y=fopen(b,"rb");CHECK(x&&y);unsigned char p[4096],q[4096];int ok=1;size_t n,m;do{n=fread(p,1,sizeof p,x);m=fread(q,1,sizeof q,y);if(n!=m||memcmp(p,q,n)){ok=0;break;}}while(n);CHECK(!ferror(x)&&!ferror(y));fclose(x);fclose(y);return ok;}
static void cache(uint32_t n,uint32_t width,int mixed){
 size_t np=n/width,nv=(size_t)width*(width-1)/2;
 double *d=malloc((size_t)n*8),*off=calloc(n,8),*v=calloc(nv?nv:1,8);uint8_t *kind=malloc(n);uint32_t *perm=malloc((size_t)n*4);p55_panel_input *p=calloc(np,sizeof *p);CHECK(d&&off&&v&&kind&&perm&&p);
 for(uint32_t i=0;i<n;i++){d[i]=2;kind[i]=1;perm[i]=n-1-i;}
 for(size_t i=0;i<np;i++){uint32_t first=(uint32_t)i*width;p[i]=(p55_panel_input){first,width,0,NULL,v,nv};if(mixed){kind[first]=2;kind[first+1]=0;d[first]=d[first+1]=0;off[first]=2;}}
 p56_context *x=NULL,*y=NULL;CHECK(!p56_rebuild(&x,n,np,p,kind,d,off,perm));
 char tmp[]="/tmp/lis68-XXXXXX";CHECK(mkdtemp(tmp));char a[128],b[128];snprintf(a,sizeof a,"%s/a",tmp);snprintf(b,sizeof b,"%s/b",tmp);uint8_t key[32]={68};double ts[3],tb[3],ls[3],lb[3];
 for(int rep=0;rep<3;rep++){
 double t=omp_get_wtime();CHECK(!scalar_save(x,a,key));ts[rep]=omp_get_wtime()-t;t=omp_get_wtime();CHECK(!p57_save(x,b,key));tb[rep]=omp_get_wtime()-t;CHECK(files_same(a,b));
 t=omp_get_wtime();CHECK(!scalar_load(&y,b,key,UINT64_MAX,n));ls[rep]=omp_get_wtime()-t;same(x,y,n);p56_destroy(y);y=NULL;
 t=omp_get_wtime();CHECK(!p57_load(&y,a,key,UINT64_MAX,n));lb[rep]=omp_get_wtime()-t;same(x,y,n);p56_destroy(y);y=NULL;tests++;
 }
 struct stat st;CHECK(!stat(a,&st));
 double *arrays[]={ts,tb,ls,lb};for(int k=0;k<4;k++)for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)if(arrays[k][j]<arrays[k][i]){double t=arrays[k][i];arrays[k][i]=arrays[k][j];arrays[k][j]=t;}
 printf("cache n=%u panel=%u pairs=%d bytes=%lld save_scalar=%.6f save_bulk=%.6f load_scalar=%.6f load_bulk=%.6f\n",n,width,mixed,(long long)st.st_size,ts[1],tb[1],ls[1],lb[1]);
 CHECK(!unlink(a)&&!unlink(b)&&!rmdir(tmp));p56_destroy(x);free(d);free(off);free(v);free(kind);free(perm);free(p);
}
static void concurrent(void){int ok=1;omp_set_max_active_levels(2);
 #pragma omp parallel for num_threads(2) reduction(&:ok)
 for(int i=0;i<60;i++){
  size_t cp[]={0,2,4,5,6};uint32_t ri[]={0,2,1,3,2,3};double a[]={0,2,0,3,0,0};p56_context *p=NULL;p67_stats st={0};p62_limits sl={1000000,1000000};p61_limits nl={16,1000000,1000,1000,0};double r[]={2,3,2,3},x[4],w[4];
  int rc=p67_factorize(&p,4,6,cp,ri,a,&sl,&nl,1,2,&st);ok&=!rc&&st.numeric.pairs==2&&st.threads_used==2;if(!rc)ok&=!p56_apply(p,r,x,w,4)&&x[0]==1&&x[1]==1&&x[2]==1&&x[3]==1;p56_destroy(p);
 }
 CHECK(ok);tests+=60;omp_set_max_active_levels(1);
}
int main(void){omp_set_dynamic(0);chain(152592,16);chain(512,1);chain(1,16);cache(65536,64,0);cache(8192,32,1);concurrent();printf("stage68 %u/%u PASS\n",tests,tests);return 0;}
