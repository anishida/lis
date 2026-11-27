#include "parallel_ldlt.h"
#include <omp.h>
#include <stdatomic.h>
#include <fenv.h>
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif
#include "schur_simd.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stddef.h>
#define NONE UINT32_MAX
typedef union {max_align_t align;size_t size;} hdr;
typedef struct {size_t live,peak,limit;atomic_flag lock;} budget;
static void lock_budget(budget *b){while(atomic_flag_test_and_set_explicit(&b->lock,memory_order_acquire)){} }
static void unlock_budget(budget *b){atomic_flag_clear_explicit(&b->lock,memory_order_release);}
static _Thread_local int alloc_error;
static void *alloc(budget *b,size_t n,size_t s){
 if(!n)return NULL;
 if(s>SIZE_MAX/n||n*s>SIZE_MAX-sizeof(hdr)){alloc_error=P60_LIMIT;return NULL;}
 size_t z=n*s;
 lock_budget(b);
 if(z>b->limit-b->live){unlock_budget(b);alloc_error=P60_LIMIT;return NULL;}
 b->live+=z;if(b->live>b->peak)b->peak=b->live;
 unlock_budget(b);
 hdr *h=calloc(1,sizeof *h+z);
 if(!h){lock_budget(b);b->live-=z;unlock_budget(b);alloc_error=P55_NOMEM;return NULL;}
 h->size=z;return h+1;
}
static void release(budget *b,void *p){if(p){hdr *h=(hdr*)p-1;size_t z=h->size;free(h);lock_budget(b);b->live-=z;unlock_budget(b);}}
typedef struct {uint32_t m,pending;uint32_t *id;double *a;} update;
typedef struct {uint32_t row;double x,y;} entry;
typedef struct {uint32_t first,w,nr;entry *e;} pivot;
static void clear(budget *b,update *u){release(b,u->id);release(b,u->a);*u=(update){0};}
static int cv(long double x,double *v){if(!isfinite(x)||fabsl(x)>DBL_MAX)return 0;*v=(double)x;return 1;}
static void sw(double *a,uint32_t *id,size_t m,size_t k,size_t i,size_t j){
 if(i==j)return;
 for(size_t r=k;r<m;r++){double v=a[i*m+r];a[i*m+r]=a[j*m+r];a[j*m+r]=v;}
 for(size_t r=k;r<m;r++){double v=a[r*m+i];a[r*m+i]=a[r*m+j];a[r*m+j]=v;}
 uint32_t t=id[i];id[i]=id[j];id[j]=t;
}
static int cmp(const void *a,const void *b){uint32_t x=((const entry*)a)->row,y=((const entry*)b)->row;return (x>y)-(x<y);}
static int add(uint32_t id,uint32_t *ids,uint32_t *map,size_t *m,size_t cap){
 if(map[id]!=NONE)return 1;
 if(*m>=cap)return 0;
 map[id]=(uint32_t)*m;ids[(*m)++]=id;return 1;
}
typedef struct {uint32_t first,last;int rc;p61_stats st;} front_result;
static int reserve(budget *b,p61_stats *total,const p61_limits *lim,size_t d,size_t slots){
 int ok=0;lock_budget(b);
 if(d<=lim->max_delays-total->delays&&slots<=lim->max_factor_slots-total->factor_slots){total->delays+=d;total->factor_slots+=slots;ok=1;}
 unlock_budget(b);return ok;
}
static int front(uint32_t lo,uint32_t f,const size_t *col,const uint32_t *row,const double *val,const uint32_t *parent,const uint32_t *head,const uint32_t *next,update *u,pivot *pv,uint8_t *kind,double *diag,double *off,uint32_t *link,front_result *fr,uint32_t *map,uint32_t *ids,size_t cap,budget *b,p61_stats *total,const p61_limits *lim,long double initial,p66_kernel schur){
 int rc=0;double *a=NULL;long double peak=initial;const long double alpha=(1+sqrtl(17))/8;
 fr->first=fr->last=NONE;
 #define NEW(p,count) do{p=alloc(b,(count),sizeof *(p));if(!(p)){rc=alloc_error;goto exit;}}while(0)
  size_t m=0,k=0,eligible;
  /* Candidate prefix: original column and each child's delayed variables. */
  for(uint32_t j=lo;j<=f;j++)if(!add(j,ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  for(uint32_t c=head[f];c!=NONE;c=next[c])for(uint32_t i=0;i<u[c].pending;i++)if(!add(u[c].id[i],ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  eligible=m;
  for(size_t z=col[lo];z<col[f+1];z++)if(!add(row[z],ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  for(uint32_t c=head[f];c!=NONE;c=next[c])for(uint32_t i=u[c].pending;i<u[c].m;i++)if(!add(u[c].id[i],ids,map,&m,cap)){rc=P60_LIMIT;goto exit;}
  fr->st.peak_front=m;
  if(m>SIZE_MAX/m){rc=P60_LIMIT;goto exit;}NEW(a,m*m);
  for(uint32_t j=lo;j<=f;j++)for(size_t z=col[j];z<col[j+1];z++){size_t i=map[row[z]],h=map[j];a[i*m+h]=a[h*m+i]=val[z];}
  for(uint32_t c=head[f];c!=NONE;c=next[c]){
   update *q=&u[c];for(size_t i=0;i<q->m;i++)for(size_t j=0;j<=i;j++){
    size_t x=map[q->id[i]],y=map[q->id[j]];
    if(!cv((long double)a[x*m+y]+q->a[i*q->m+j],&a[x*m+y])){rc=P55_RANGE;goto exit;}a[y*m+x]=a[x*m+y];
    peak=fmaxl(peak,fabsl(a[x*m+y]));
   }
   clear(b,q);
  }
  for(size_t i=0;i<m;i++)map[ids[i]]=NONE;
  if(lim->max_growth&&peak/initial>lim->max_growth){rc=P59_GROWTH;goto exit;}
  size_t candidates=eligible;
  while(k<eligible){
   long double akk=fabsl(a[k*m+k]),colmax=0;size_t r=k,w=1,chosen=k;
   for(size_t i=k+1;i<m;i++)if(fabsl(a[i*m+k])>colmax){colmax=fabsl(a[i*m+k]);r=i;}
   if(akk==0&&colmax==0){rc=P58_NO_PIVOT;goto exit;}
   if(akk<alpha*colmax){long double rm=0;for(size_t i=k;i<m;i++)if(i!=r)rm=fmaxl(rm,fabsl(a[r*m+i]));
    if(akk<alpha*colmax*(colmax/rm)){chosen=r;if(fabsl(a[r*m+r])<alpha*rm)w=2;}}
   size_t target=k+w-1;
   if(chosen>=eligible||target>=eligible){
    if(parent[f]==NONE){rc=P58_NO_PIVOT;goto exit;}
    if(!reserve(b,total,lim,1,0)){rc=P60_LIMIT;goto exit;}--eligible;
    if(k!=eligible){sw(a,ids,m,k,k,eligible);fr->st.swaps++;}fr->st.delays++;continue;
   }
   if(chosen!=target){sw(a,ids,m,k,target,chosen);fr->st.swaps++;}
   long double aa=a[k*m+k],bb=0,cc=0,s=1,det=1;
   if(w==2){bb=a[(k+1)*m+k];cc=a[(k+1)*m+k+1];s=fmaxl(fabsl(aa),fmaxl(fabsl(bb),fabsl(cc)));aa/=s;bb/=s;cc/=s;det=fmal(aa,cc,-bb*bb);if(!isfinite(det)||fabsl(det)<=8*LDBL_EPSILON*(fabsl(aa*cc)+bb*bb)){rc=P58_NO_PIVOT;goto exit;}}
   uint32_t first=ids[k];kind[first]=(uint8_t)w;diag[first]=a[k*m+k];
   if(fr->first==NONE)fr->first=first;else link[fr->last]=first;
   fr->last=first;link[first]=NONE;
   if(w==2){uint32_t second=ids[k+1];diag[second]=a[(k+1)*m+k+1];off[first]=a[(k+1)*m+k];link[first]=second;link[second]=NONE;fr->last=second;fr->st.pairs++;}else fr->st.scalars++;
   pivot *p=&pv[first];p->first=first;p->w=(uint32_t)w;
   /* Store multipliers in unused upper part of the pivot rows; lower
      columns retain original coupling for the Schur update. */
   for(size_t i=k+w;i<m;i++){
    long double x=a[i*m+k],y=w==2?a[i*m+k+1]:0;
    if(w==1)x/=aa;else{long double z=x/s,v=y/s;x=fmal(cc,z,-bb*v)/det;y=fmal(aa,v,-bb*z)/det;}
    if(!cv(x,&a[k*m+i])||(w==2&&!cv(y,&a[(k+1)*m+i]))){rc=P55_RANGE;goto exit;}
    fr->st.max_multiplier=fmaxl(fr->st.max_multiplier,fmaxl(fabsl(x),fabsl(y)));
    if(a[k*m+i]!=0||(w==2&&a[(k+1)*m+i]!=0))p->nr++;
   }
   size_t slots=(size_t)p->nr*w+(w-1);
   if(!reserve(b,total,lim,0,slots)){rc=P60_LIMIT;goto exit;}fr->st.factor_slots+=slots;
   if(p->nr){NEW(p->e,p->nr);size_t at=0;for(size_t i=k+w;i<m;i++)if(a[k*m+i]!=0||(w==2&&a[(k+1)*m+i]!=0))p->e[at++]=(entry){ids[i],a[k*m+i],w==2?a[(k+1)*m+i]:0};}
   if(!schur(a,m,k,w,&peak)){rc=P55_RANGE;goto exit;}
   if(lim->max_growth&&peak/initial>lim->max_growth){rc=P59_GROWTH;goto exit;}
   k+=w;
  }
  if(k<m){
   if(parent[f]==NONE){rc=P58_NO_PIVOT;goto exit;}
   update *q=&u[f];q->m=(uint32_t)(m-k);q->pending=(uint32_t)(candidates-k);NEW(q->id,q->m);NEW(q->a,(size_t)q->m*q->m);
   memcpy(q->id,ids+k,(size_t)q->m*sizeof *ids);for(size_t i=0;i<q->m;i++)memcpy(q->a+i*q->m,a+(i+k)*m+k,(size_t)q->m*sizeof *a);
  }

 fr->st.growth=peak/initial;fr->st.fronts=1;
 exit:
 for(size_t i=0;i<m;i++)map[ids[i]]=NONE;
 release(b,a);return rc;
 #undef NEW
}
int p67_numeric(p56_context **slot,uint32_t n,const size_t *col,const uint32_t *row,const double *val,const uint32_t *parent,const uint32_t *bounds,size_t ng,const p61_limits *lim,unsigned threads,p67_stats *stats){
 if(!threads||threads>64)return P55_INVALID;
 if(threads==1){p61_stats st={0};int rc=p65_numeric(slot,n,col,row,val,parent,bounds,ng,lim,&st);if(!rc&&stats)*stats=(p67_stats){st,1,0,1};return rc;}
 if(!slot||!n||n==NONE||!col||!parent||!lim||!lim->max_front||!lim->max_work_bytes||col[0]||(!row&&col[n])||(!val&&col[n])||!isfinite(lim->max_growth)||(lim->max_growth&&lim->max_growth<1))return P55_INVALID;
 /* A contracted chain has no concurrent front. Avoid allocating any
    parallel planning arrays, even when a large team was requested. */
 size_t nf0=bounds?ng:n;int chain=nf0>0&&nf0<=n;
 if(bounds&&chain)chain=bounds[0]==0&&bounds[ng]==n;
 for(size_t g=0;chain&&g<nf0;g++){
  uint32_t lo=bounds?bounds[g]:(uint32_t)g,hi=bounds?bounds[g+1]:(uint32_t)g+1;
  if(lo>=hi||hi>n){chain=0;break;}
  for(uint32_t j=lo;j+1<hi;j++)if(parent[j]!=j+1){chain=0;break;}
  if(g+1==nf0)chain=chain&&parent[hi-1]==NONE;
  else{uint32_t end=bounds?bounds[g+2]:hi+1;chain=chain&&end>hi&&end<=n&&parent[hi-1]>=hi&&parent[hi-1]<end;}
 }
 if(chain){p61_stats st={0};int r=p65_numeric(slot,n,col,row,val,parent,bounds,ng,lim,&st);if(!r&&stats)*stats=(p67_stats){st,1,nf0,1};return r;}
 budget b={0,0,lim->max_work_bytes,ATOMIC_FLAG_INIT};p61_stats st={0};p67_stats result={0};int rc=P55_INVALID;
 uint32_t *head=NULL,*next=NULL,*tin=NULL,*tout=NULL,*stack=NULL,*perm=NULL,*inv=NULL,*group=NULL;uint8_t *kind=NULL;
 double *diag=NULL,*off=NULL;update *u=NULL;pivot *pv=NULL;p55_panel_input *out=NULL;
 uint32_t *link=NULL,*levels=NULL,*queue=NULL,*maps=NULL,*allids=NULL,*lo_of=NULL;size_t *offsets=NULL,*cursor=NULL;front_result *fr=NULL;uint8_t *ordered_kind=NULL;double *ordered_diag=NULL,*ordered_off=NULL;size_t used=0,done=0,maxlevel=0,width=0;long double initial=0;size_t cap=lim->max_front<n?lim->max_front:n;
 #define NEW(p,count) do{p=alloc(&b,(count),sizeof *(p));if(!(p)){rc=alloc_error;goto exit;}}while(0)
 NEW(head,n);NEW(next,n);NEW(tin,n);NEW(tout,n);NEW(stack,n);NEW(inv,n);
 for(uint32_t i=0;i<n;i++){head[i]=NONE; if(parent[i]!=NONE&&(parent[i]<=i||parent[i]>=n))goto exit; if(col[i]>col[i+1])goto exit;}
 for(uint32_t i=0;i<n;i++){next[i]=NONE;if(parent[i]!=NONE){next[i]=head[parent[i]];head[parent[i]]=i;}}
 /* DFS entry/exit intervals, without recursion. Reuse inv for child cursor. */
 uint32_t tick=0;
 for(uint32_t root=0;root<n;root++)if(parent[root]==NONE){size_t depth=0;stack[depth++]=root;tin[root]=tick++;inv[root]=head[root];while(depth){uint32_t v=stack[depth-1],c=inv[v];if(c==NONE){tout[v]=tick;depth--;continue;}inv[v]=next[c];tin[c]=tick++;inv[c]=head[c];stack[depth++]=c;}}
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++){
  uint32_t i=row[z];if(i<j||i>=n||(z>col[j]&&i<=row[z-1])||!isfinite(val[z])||!(tin[i]<=tin[j]&&tin[j]<tout[i]))goto exit;
  initial=fmaxl(initial,fabsl(val[z]));
 }
 if(bounds){
  if(!ng||ng>n||bounds[0]!=0||bounds[ng]!=n)goto exit;
  NEW(group,n);
  for(size_t g=0;g<ng;g++){
   uint32_t lo=bounds[g],hi=bounds[g+1];if(lo>=hi||hi>n)goto exit;
   for(uint32_t j=lo;j<hi;j++){if(j+1<hi&&parent[j]!=j+1)goto exit;group[j]=hi-1;}
  }
  /* Contract chains into representative (last column) front nodes. */
  for(uint32_t i=0;i<n;i++)head[i]=next[i]=NONE;
  for(size_t g=0;g<ng;g++){uint32_t f=bounds[g+1]-1,p=parent[f];if(p!=NONE){p=group[p];next[f]=head[p];head[p]=f;}}
 }
 if(initial==0){rc=P58_NO_PIVOT;goto exit;}
 release(&b,tin);tin=NULL;release(&b,tout);tout=NULL;release(&b,stack);stack=NULL;
 /* Contracted front level = max(child level + 1). Stable counting sort. */
 NEW(levels,n);NEW(offsets,(size_t)n+1);
 size_t nf=bounds?ng:n;
 for(size_t g=0;g<nf;g++){uint32_t lo=bounds?bounds[g]:(uint32_t)g,f=bounds?bounds[g+1]-1:(uint32_t)g;(void)lo;size_t lv=0;for(uint32_t c=head[f];c!=NONE;c=next[c])if((size_t)levels[c]+1>lv)lv=(size_t)levels[c]+1;levels[f]=(uint32_t)lv;offsets[lv+1]++;if(lv>maxlevel)maxlevel=lv;}
 for(size_t i=0;i<=maxlevel;i++)offsets[i+1]+=offsets[i];
 for(size_t l=0;l<=maxlevel;l++){size_t w=offsets[l+1]-offsets[l];if(w>width)width=w;}
 if(threads>width)threads=(unsigned)width;
 /* Avoid barriers and per-thread workspace on purely serial elimination trees. */
 if(threads==1){rc=0;goto serial_fallback;}
 NEW(queue,nf);NEW(cursor,maxlevel+1);NEW(lo_of,n);NEW(fr,n);NEW(link,n);
 memcpy(cursor,offsets,(maxlevel+1)*sizeof *cursor);
 for(size_t g=0;g<nf;g++){uint32_t f=bounds?bounds[g+1]-1:(uint32_t)g;queue[cursor[levels[f]]++]=f;lo_of[f]=bounds?bounds[g]:(uint32_t)g;}
 release(&b,cursor);cursor=NULL;release(&b,levels);levels=NULL;release(&b,group);group=NULL;
 NEW(perm,n);NEW(kind,n);NEW(diag,n);NEW(off,n);NEW(u,n);NEW(pv,n);
 if((size_t)threads>SIZE_MAX/n||(size_t)threads>SIZE_MAX/cap){rc=P60_LIMIT;goto exit;}
 NEW(maps,(size_t)threads*n);NEW(allids,(size_t)threads*cap);
 memset(maps,0xff,(size_t)threads*n*sizeof *maps);
 p66_kernel schur=p66_select_kernel();fenv_t env;if(fegetenv(&env)){rc=P55_INVALID;goto exit;}
 #if defined(__x86_64__) || defined(__i386__)
 unsigned mxcsr=_mm_getcsr();
 #endif
 result.levels=maxlevel+1;result.max_ready=width;result.threads_used=1;
 atomic_int failed=0;
 #pragma omp parallel num_threads(threads) shared(failed)
 {
  unsigned tid=(unsigned)omp_get_thread_num();fenv_t oldenv;int have_env=fegetenv(&oldenv)==0;
  #if defined(__x86_64__) || defined(__i386__)
  unsigned oldmx=_mm_getcsr();
  #endif
  int env_ok=have_env&&fesetenv(&env)==0;
  #if defined(__x86_64__) || defined(__i386__)
  _mm_setcsr(mxcsr);
  #endif
  if(!env_ok)atomic_store(&failed,P55_INVALID);
  #pragma omp barrier
  #pragma omp single
  {result.threads_used=(unsigned)omp_get_num_threads();}
  for(size_t l=0;l<=maxlevel;l++){
   if(failed)break;
   #pragma omp for schedule(dynamic,1)
   for(size_t j=offsets[l];j<offsets[l+1];j++){
    uint32_t f=queue[j];fr[f].rc=front(lo_of[f],f,col,row,val,parent,head,next,u,pv,kind,diag,off,link,&fr[f],maps+(size_t)tid*n,allids+(size_t)tid*cap,cap,&b,&st,lim,initial,schur);
   }
   #pragma omp single
   {for(size_t j=offsets[l];j<offsets[l+1];j++)if(fr[queue[j]].rc){failed=fr[queue[j]].rc;break;}}
  }
  if(have_env)fesetenv(&oldenv);
  #if defined(__x86_64__) || defined(__i386__)
  _mm_setcsr(oldmx);
  #endif
 }
 if(failed){rc=failed;goto exit;}
 release(&b,maps);maps=NULL;release(&b,allids);allids=NULL;
 release(&b,head);head=NULL;release(&b,next);next=NULL;release(&b,queue);queue=NULL;release(&b,offsets);offsets=NULL;release(&b,lo_of);lo_of=NULL;
 NEW(ordered_kind,n);NEW(ordered_diag,n);NEW(ordered_off,n);
 /* Assemble in original topological front order, independent of scheduling. */
 for(size_t g=0;g<nf;g++){
  uint32_t f=bounds?bounds[g+1]-1:(uint32_t)g;
  st.fronts+=fr[f].st.fronts;st.pairs+=fr[f].st.pairs;st.scalars+=fr[f].st.scalars;st.swaps+=fr[f].st.swaps;
  if(fr[f].st.peak_front>st.peak_front)st.peak_front=fr[f].st.peak_front;
  st.growth=fmaxl(st.growth,fr[f].st.growth);st.max_multiplier=fmaxl(st.max_multiplier,fr[f].st.max_multiplier);
  for(uint32_t id=fr[f].first;id!=NONE;id=link[id]){if(done>=n){rc=P55_INVALID;goto exit;}perm[done]=id;ordered_kind[done]=kind[id];ordered_diag[done]=diag[id];ordered_off[done]=off[id];done++;}
 }
 if(done!=n){rc=P58_NO_PIVOT;goto exit;}
 release(&b,fr);fr=NULL;release(&b,link);link=NULL;release(&b,kind);kind=NULL;release(&b,diag);diag=NULL;release(&b,off);off=NULL;
 release(&b,u);u=NULL; /* all child contributions consumed, roots have no update */
 NEW(out,n);
 for(uint32_t i=0;i<n;i++)inv[perm[i]]=i;
 for(size_t j=0;j<n;){
  pivot *p=&pv[perm[j]];if(p->w!=1&&p->w!=2){rc=P55_INVALID;goto exit;}
  for(uint32_t i=0;i<p->nr;i++)p->e[i].row=inv[p->e[i].row];
  if(p->nr)qsort(p->e,p->nr,sizeof *p->e,cmp);
  uint32_t *rr=NULL;double *vv=NULL;size_t nv=(size_t)p->nr*p->w+p->w-1;
  out[used]=(p55_panel_input){(uint32_t)j,p->w,p->nr,NULL,NULL,nv};used++;
  if(p->nr){NEW(rr,p->nr);out[used-1].rows=rr;}
  if(nv){NEW(vv,nv);out[used-1].values=vv;}
  for(uint32_t i=0;i<p->nr;i++){rr[i]=p->e[i].row;vv[i+p->w-1]=p->e[i].x;if(p->w==2)vv[1+p->nr+i]=p->e[i].y;}
  release(&b,p->e);p->e=NULL;j+=p->w;
 }
 st.peak_work_bytes=b.peak;result.numeric=st;
 rc=p56_rebuild(slot,n,used,out,ordered_kind,ordered_diag,ordered_off,perm);if(!rc&&stats)*stats=result;
 goto exit;
 serial_fallback: rc=-999;
 exit:
 if(u)for(uint32_t i=0;i<n;i++)clear(&b,&u[i]);
 if(pv)for(uint32_t i=0;i<n;i++)release(&b,pv[i].e);
 if(out)for(size_t i=0;i<used;i++){release(&b,(void*)out[i].rows);release(&b,(void*)out[i].values);}

 release(&b,group);release(&b,head);release(&b,next);release(&b,tin);release(&b,tout);release(&b,stack);release(&b,perm);release(&b,inv);release(&b,kind);release(&b,diag);release(&b,off);release(&b,u);release(&b,pv);release(&b,out);release(&b,link);release(&b,levels);release(&b,queue);release(&b,maps);release(&b,allids);release(&b,lo_of);release(&b,offsets);release(&b,cursor);release(&b,fr);release(&b,ordered_kind);release(&b,ordered_diag);release(&b,ordered_off);
 if(rc==-999){rc=p65_numeric(slot,n,col,row,val,parent,bounds,ng,lim,&st);if(!rc){if(b.peak>st.peak_work_bytes)st.peak_work_bytes=b.peak;if(stats)*stats=(p67_stats){st,1,maxlevel+1,width};}}return rc;
 #undef NEW
}
int p67_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *col,const uint32_t *row,const double *val,const p62_limits *sl,const p61_limits *nl,uint32_t max_width,unsigned threads,p67_stats *stats){
 if(!slot||!nl||!max_width||(nnz&&!val))return P55_INVALID;

 p62_symbolic *s=NULL;int rc=p62_analyze(&s,n,nnz,col,row,sl);if(rc)return rc;
 /* Reuse the already allocated n+1 boundary array. Fill backwards so
    expansion cannot overwrite an unread original boundary. */
 size_t ng=0;
 for(size_t g=0;g<s->nsuper;g++){uint32_t w=s->super[g+1]-s->super[g];ng+=1+(w-1)/max_width;}
 if(ng==n){rc=p67_numeric(slot,n,col,row,val,s->parent,NULL,0,nl,threads,stats);p62_free(s);return rc;}
 size_t k=ng;uint32_t hi=n;s->super[k]=n;
 for(size_t g=s->nsuper;g>0;g--){uint32_t lo=s->super[g-1];uint32_t w=hi-lo;size_t pieces=1+(w-1)/max_width;
  for(size_t j=pieces;j>0;j--)s->super[--k]=lo+(uint32_t)(j-1)*max_width;
  hi=lo;
 }
 rc=p67_numeric(slot,n,col,row,val,s->parent,s->super,ng,nl,threads,stats);p62_free(s);return rc;
}

#include "factor_private.h"
typedef struct {uint32_t c,r;double v;} order_entry;
static int order_cmp(const void *a,const void *b){const order_entry *x=a,*y=b;if(x->c!=y->c)return x->c<y->c?-1:1;return (x->r>y->r)-(x->r<y->r);}
int p67_ordered_factorize(p56_context **slot,uint32_t n,size_t nnz,const size_t *cp,const uint32_t *ri,const double *v,const p64_limits *ol,const p62_limits *sl,const p61_limits *nl,p64_method method,uint32_t max_width,unsigned threads,p64_stats *stats,p67_stats *parallel){
 if(!slot||!nl||!max_width||!threads||threads>64||(nnz&&!v))return P55_INVALID;
 p64_stats st={0};p67_stats ps={0};p63_order *o=NULL;int rc=p64_order(n,nnz,cp,ri,ol,sl,method,&o,&st);if(rc)return rc;
 p56_context *ctx=NULL;size_t *c=NULL;uint32_t *r=NULL;double *a=NULL;order_entry *e=NULL;
 if(st.selected==P64_NATURAL)rc=p67_factorize(&ctx,n,nnz,cp,ri,v,sl,nl,max_width,threads,&ps);
 else{
  budget b={0,0,ol->max_work_bytes,ATOMIC_FLAG_INIT};
  /* Charge the retained permutation too; arrays released on every exit. */
  if((size_t)n>b.limit/(2*sizeof(uint32_t))){rc=P60_LIMIT;goto end;}
  b.live=(size_t)n*2*sizeof(uint32_t);
  c=alloc(&b,(size_t)n+1,sizeof *c);if(!c){rc=alloc_error;goto ordered_end;}
  if(nnz){r=alloc(&b,nnz,sizeof *r);if(!r){rc=alloc_error;goto ordered_end;}a=alloc(&b,nnz,sizeof *a);if(!a){rc=alloc_error;goto ordered_end;}e=alloc(&b,nnz,sizeof *e);if(!e){rc=alloc_error;goto ordered_end;}}
  for(uint32_t j=0;j<n;j++)for(size_t z=cp[j];z<cp[j+1];z++){uint32_t x=o->inv[j],y=o->inv[ri[z]];e[z]=(order_entry){x<y?x:y,x<y?y:x,v[z]};}
  if(nnz)qsort(e,nnz,sizeof *e,order_cmp);
  for(size_t z=0;z<nnz;z++){c[e[z].c+1]++;r[z]=e[z].r;a[z]=e[z].v;}
  for(uint32_t j=0;j<n;j++)c[j+1]+=c[j];
  release(&b,e);e=NULL;
  rc=p67_factorize(&ctx,n,nnz,c,r,a,sl,nl,max_width,threads,&ps);
  if(!rc){int identity=1;for(uint32_t j=0;j<n;j++){o->inv[j]=o->perm[ctx->perm?ctx->perm[j]:j];if(o->inv[j]!=j)identity=0;}free(ctx->perm);ctx->perm=identity?NULL:o->inv;if(!identity)o->inv=NULL;}
 ordered_end:release(&b,e);release(&b,c);release(&b,r);release(&b,a);
 }
 end:if(!rc){st.numeric=ps.numeric;p56_destroy(*slot);*slot=ctx;ctx=NULL;if(stats)*stats=st;if(parallel)*parallel=ps;}p56_destroy(ctx);p63_free(o);return rc;
}
