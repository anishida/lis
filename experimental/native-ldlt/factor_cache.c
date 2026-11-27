#define _POSIX_C_SOURCE 200809L
#include "factor_cache.h"
#include "factor_private.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
_Static_assert(CHAR_BIT==8 && sizeof(double)==8 && DBL_MANT_DIG==53 && DBL_MAX_EXP==1024 && DBL_MIN_EXP==-1021, "binary64 required");
static const unsigned char magic[8]={'L','D','L','T','C','0','0','1'};
typedef struct { FILE *fp; uint32_t crc,tab[256]; uint64_t left; int error; } stream;
static void init(stream *s,FILE *fp,uint64_t left){
    s->fp=fp;s->crc=UINT32_MAX;s->left=left;s->error=0;
    for(unsigned i=0;i<256;i++){uint32_t v=i;for(int j=0;j<8;j++)v=(v>>1)^((v&1)?UINT32_C(0xedb88320):0);s->tab[i]=v;}
}
static void crc(stream *s,const unsigned char *p,size_t n){for(size_t i=0;i<n;i++)s->crc=(s->crc>>8)^s->tab[(s->crc^p[i])&255];}
static void put(stream *s,const void *p,size_t n){if(s->error)return;if(fwrite(p,1,n,s->fp)!=n){s->error=1;return;}crc(s,p,n);}
static void get(stream *s,void *p,size_t n){if(s->error||n>s->left){s->error=1;return;}if(fread(p,1,n,s->fp)!=n){s->error=1;return;}s->left-=n;crc(s,p,n);}
static void putint(stream *s,uint64_t x,unsigned n){unsigned char b[8];for(unsigned i=0;i<n;i++)b[i]=(unsigned char)(x>>(8*i));put(s,b,n);}
static uint64_t getint(stream *s,unsigned n){unsigned char b[8]={0};get(s,b,n);uint64_t v=0;for(unsigned i=0;i<n;i++)v|=(uint64_t)b[i]<<(8*i);return v;}
static void putd(stream *s,double d){uint64_t v;memcpy(&v,&d,8);putint(s,v,8);}
static double getd(stream *s){uint64_t v=getint(s,8);double d;memcpy(&d,&v,8);return d;}
static int layout_ok(void){double d=1.;uint64_t v;memcpy(&v,&d,8);return v==UINT64_C(0x3ff0000000000000);}
/* Preserve LDLTC001 bytes; bulk I/O only on native little endian binary64. */
static int native_le(void){
#ifdef P68_CACHE_SCALAR
    return 0;
#else
    uint16_t x=1;unsigned char b[2];memcpy(b,&x,2);return b[0]==1;
#endif
}
static void putda(stream *s,const double *v,size_t n){if(!n)return;if(native_le())put(s,v,n*sizeof *v);else for(size_t j=0;j<n;j++)putd(s,v[j]);}
static void getda(stream *s,double *v,size_t n){if(!n)return;if(native_le())get(s,v,n*sizeof *v);else for(size_t j=0;j<n;j++)v[j]=getd(s);}
static void putua(stream *s,const void *v,size_t n,unsigned sz){
    if(!n)return;
    if(native_le()){put(s,v,n*sz);return;}
    for(size_t j=0;j<n;j++)putint(s,sz==2?((const uint16_t*)v)[j]:((const uint32_t*)v)[j],sz);
}
int p57_save(const p56_context *ctx,const char *path,const uint8_t key[32]){
    if(!ctx||!path||!*path||!key||!layout_ok())return P55_INVALID;
    size_t len=strlen(path);if(len>SIZE_MAX-12)return P55_INVALID;
    char *tmp=malloc(len+12);if(!tmp)return P55_NOMEM;
    memcpy(tmp,path,len);memcpy(tmp+len,".tmp.XXXXXX",12);
    int fd=mkstemp(tmp);if(fd<0){free(tmp);return P57_IO;}
    FILE *fp=fdopen(fd,"wb");if(!fp){close(fd);unlink(tmp);free(tmp);return P57_IO;}
    stream s;init(&s,fp,0);const p55_factor *f=ctx->factor;
    put(&s,magic,8);put(&s,key,32);putint(&s,f->n,4);putint(&s,f->np,4);putint(&s,f->nb,4);putint(&s,ctx->perm!=NULL,4);
    putda(&s,f->diag,f->n);
    for(size_t j=0;j<f->nb;j++){putint(&s,f->blocks[j].k,4);putd(&s,f->blocks[j].original_b);}
    if(ctx->perm)putua(&s,ctx->perm,f->n,4);
    for(size_t j=0;j<f->np;j++){
        const panel *p=&f->p[j];putint(&s,p->width,4);putint(&s,p->nrows,4);putint(&s,p->short_index?2:4,4);
        putua(&s,p->rows,p->nrows,p->short_index?2:4);
        putda(&s,p->values,p->nvalues);
    }
    uint32_t sum=s.crc^UINT32_MAX;putint(&s,sum,4);
    int ok=!s.error && fflush(fp)==0 && fsync(fd)==0;
    if(fclose(fp)!=0)ok=0;
    if(ok && rename(tmp,path)!=0)ok=0;
    if(!ok)unlink(tmp);
    free(tmp);return ok?P55_OK:P57_IO;
}
static void release(p55_panel_input *p,size_t np){if(p)for(size_t i=0;i<np;i++){free((void*)p[i].rows);free((void*)p[i].values);}free(p);}
int p57_load(p56_context **slot,const char *path,const uint8_t key[32],uint64_t max_bytes,uint32_t max_n){
    if(!slot||!path||!*path||!key||!max_bytes||!max_n||!layout_ok())return P55_INVALID;
    int fd=open(path,O_RDONLY|O_NONBLOCK);if(fd<0)return errno==ENOENT?P57_MISS:P57_IO;
    FILE *fp=fdopen(fd,"rb");if(!fp){close(fd);return P57_IO;}
    struct stat st;int rc=P57_CORRUPT;double *d=NULL,*off=NULL;uint8_t *kind=NULL;uint32_t *perm=NULL; p55_panel_input *p=NULL;uint32_t np=0; p56_context *fresh=NULL;
    if(fstat(fileno(fp),&st)!=0||!S_ISREG(st.st_mode)){rc=P57_IO;goto done;}
    if(st.st_size<60)goto done;
    if((uint64_t)st.st_size>max_bytes){rc=P57_LIMIT;goto done;}
    stream s;init(&s,fp,(uint64_t)st.st_size-4);unsigned char hdr[8]={0},gotkey[32]={0};get(&s,hdr,8);get(&s,gotkey,32);
    if(s.error)goto done;
    if(memcmp(hdr,magic,8)||memcmp(gotkey,key,32)){rc=P57_MISS;goto done;}
    uint32_t n=(uint32_t)getint(&s,4);np=(uint32_t)getint(&s,4);uint32_t nb=(uint32_t)getint(&s,4),flag=(uint32_t)getint(&s,4);
    if(s.error||!n||!np||np>n||nb>n/2||flag>1)goto done;
    if(n>max_n){rc=P57_LIMIT;goto done;}
    if((uint64_t)n*8+(uint64_t)nb*12+(uint64_t)flag*n*4+(uint64_t)np*12>s.left)goto done;
    if(SIZE_MAX/n<sizeof(double)||SIZE_MAX/np<sizeof *p){rc=P57_LIMIT;goto done;}
    d=malloc((size_t)n*sizeof *d);kind=malloc(n);if(nb)off=calloc(n,sizeof *off);p=calloc(np,sizeof *p);if(flag)perm=malloc((size_t)n*sizeof *perm);
    if(!d||!kind||!p||(nb&&!off)||(flag&&!perm)){rc=P55_NOMEM;goto done;}
    memset(kind,1,n);getda(&s,d,n);
    uint32_t prev=0;
    for(uint32_t j=0;j<nb;j++){uint32_t k=(uint32_t)getint(&s,4);double b=getd(&s);if(s.error||k>=n-1||(j&&k<prev+2))goto done;kind[k]=2;kind[k+1]=0;off[k]=b;prev=k;}
    if(flag){if(native_le())get(&s,perm,(size_t)n*4);else for(uint32_t k=0;k<n;k++)perm[k]=(uint32_t)getint(&s,4); }
    uint32_t first=0;
    for(uint32_t j=0;j<np;j++){
        uint32_t w=(uint32_t)getint(&s,4),nr=(uint32_t)getint(&s,4),sz=(uint32_t)getint(&s,4);
        if(s.error||!w||w>n-first||nr>n-first-w||(sz!=2&&sz!=4))goto done;
        uint64_t nv=(uint64_t)w*(w-1)/2+(uint64_t)w*nr;
        if(nv>SIZE_MAX/sizeof(double)||(uint64_t)nr*sz>s.left||nv>(s.left-(uint64_t)nr*sz)/8|| (nr&&SIZE_MAX/nr<sizeof(uint32_t)))goto done;
        uint32_t *rows=nr?malloc((size_t)nr*sizeof *rows):NULL;double *values=nv?malloc((size_t)nv*sizeof *values):NULL;
        p[j]=(p55_panel_input){first,w,nr,rows,values,(size_t)nv};
        if((nr&&!rows)||(nv&&!values)){rc=P55_NOMEM;goto done;}
        uint32_t end=first+w;
        for(uint32_t r=0;r<nr;r++){uint32_t v=(uint32_t)getint(&s,sz);if(sz==2){if(v>=n-end)goto done;v+=end;}if(v<end||v>=n||(r&&v<=rows[r-1]))goto done;rows[r]=v;}
        if((!nr||rows[nr-1]-end<=UINT16_MAX)!=(sz==2))goto done;
        getda(&s,values,(size_t)nv);
        if(s.error)goto done;
        first=end;
    }
    if(first!=n||s.left||s.error)goto done;
    unsigned char tail[4];if(fread(tail,1,4,fp)!=4||fgetc(fp)!=EOF||ferror(fp))goto done;
    uint32_t expected=0;for(unsigned i=0;i<4;i++)expected|=(uint32_t)tail[i]<<(8*i);
    if(expected!=(s.crc^UINT32_MAX))goto done;
    rc=p56_rebuild(&fresh,n,np,p,kind,d,off,perm);
    if(rc!=P55_OK && rc!=P55_NOMEM)rc=P57_CORRUPT;
done:
    if(fclose(fp)!=0&&rc==P55_OK)rc=P57_IO;
    release(p,np);free(d);free(off);free(kind);free(perm);
    if(rc==P55_OK){p56_context *old=*slot;*slot=fresh;p56_destroy(old);}else p56_destroy(fresh);
    return rc;
}
