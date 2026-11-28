/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_restore.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#define C(e) do{if(!(e)){fprintf(stderr,"SPLIT FAIL %d: %s\n",__LINE__,#e);return 1;}}while(0)
typedef struct {uint64_t a,b;int calls,fail,generalized;} state;
static const double A[3][3]={{4,1,0},{-1,5,.5},{0,.5,6}};
static int generations(void *u,uint64_t *a,uint64_t *b){state *s=u;*a=s->a;*b=s->b;return 0;}
static int apply(void *u,uint32_t op,size_t n,const double *in,double *out)
{
 state *s=u;s->calls++;if(s->fail==s->calls)return 77;
 if(n!=3)return 99;
 for(size_t i=0;i<n;i++){
  out[i]=0;
  if(op==LIS_RCI_OP_B)out[i]=(i+1)*in[i];
  else if(op==LIS_RCI_OP_A)for(size_t j=0;j<n;j++)out[i]+=A[i][j]*in[j];
  else return 98;
 }
 return 0;
}
static double complex coeff(lis_rci_complex_coeff v){return (double)v.real+(double)v.imag*I;}
static void reference(int general,lis_rci_complex_coeff shift,const double complex b[3],double complex x[3])
{
 double complex a[3][4];
 for(int i=0;i<3;i++){
  for(int j=0;j<3;j++)a[i][j]=-A[i][j]+(i==j?coeff(shift)*(general?i+1:1):0);
  a[i][3]=b[i];
 }
 for(int k=0;k<3;k++){
  int p=k;for(int i=k+1;i<3;i++)if(cabs(a[i][k])>cabs(a[p][k]))p=i;
  for(int j=k;j<4;j++){double complex t=a[k][j];a[k][j]=a[p][j];a[p][j]=t;}
  for(int i=k+1;i<3;i++){double complex f=a[i][k]/a[k][k];for(int j=k;j<4;j++)a[i][j]-=f*a[k][j];}
 }
 for(int i=2;i>=0;i--){x[i]=a[i][3];for(int j=i+1;j<3;j++)x[i]-=a[i][j]*x[j];x[i]/=a[i][i];}
}
static int run(int general,lis_rci_method method)
{
 state s={.a=23,.b=40,.generalized=general};
 lis_rci_options opt={6,100,6,1e-12,0,method};
 lis_rci_split_restore *ctx=lis_rci_split_restore_create(3,general,&opt,&s,apply,generations);C(ctx);
 lis_rci_opaque_backend backend;C(!lis_rci_split_restore_backend(ctx,&backend));
 double real[10][3]={{0}},imag[10][3]={{0}};
 lis_rci_split_vector v[10];lis_rci_handle h[10];
 for(int i=0;i<10;i++){v[i]=(lis_rci_split_vector){3,real[i],imag[i]};h[i]=&v[i];}
 for(int j=0;j<2;j++)for(int i=0;i<3;i++){real[j][i]=1+i+j;imag[j][i]=.5*j-i;}
 lis_rci_complex_coeff shifts[]={{2,1.5},{7,-2},{3,3}},weights[]={{.1,.2},{.3,-.4},{-.2,.1}};
 lis_rci_projector_restore_options po={2,3,general,shifts,weights,h,h+2,h+4,h+6,h+8};
 lis_rci_projector_restore *p=lis_rci_projector_restore_create(&po);C(p);
 lis_rci_opaque_request q;int status=lis_rci_projector_restore_start(p,s.a,s.b,&q),code=0;
 while(status>0){C(!lis_rci_opaque_dispatch(&backend,&q,sizeof q,&code));status=lis_rci_projector_restore_advance(p,q.token,s.a,s.b,0,&q);}
 C(status==0);const lis_rci_handle *res=lis_rci_projector_restore_result(p);C(res);
 for(int j=0;j<2;j++){
  double complex rhs[3],want[3]={0},x[3];
  for(int i=0;i<3;i++)rhs[i]=(real[j][i]+I*imag[j][i])*(general?i+1:1);
  for(int k=0;k<3;k++){reference(general,shifts[k],rhs,x);for(int i=0;i<3;i++)want[i]+=coeff(weights[k])*x[i];}
  lis_rci_split_vector *o=res[j];
  for(int i=0;i<3;i++)C(cabs(o->real[i]+I*o->imag[i]-want[i])<2e-10);
 }
 /* Test failed shifted solve never publishes a partially computed column. */
 lis_rci_handle in=h[0],out=h[4];
 q=(lis_rci_opaque_request){.version=1,.size=sizeof q,.operation=LIS_RCI_OP_SHIFTED_SOLVE,.token=99,.generation_a=s.a,.generation_b=s.b,.columns=1,.in=&in,.out=&out,.shift=shifts[0]};
 for(int i=0;i<3;i++)real[4][i]=imag[4][i]=-9;
 s.fail=s.calls+1;
 C(lis_rci_opaque_dispatch(&backend,&q,sizeof q,&code)==LIS_RCI_CALLBACK&&code==77);
 C(lis_rci_split_restore_error(ctx,&code)==LIS_RCI_CALLBACK&&code==77);
 for(int i=0;i<3;i++)C(real[4][i]==-9&&imag[4][i]==-9);
 s.fail=0;C(!lis_rci_opaque_dispatch(&backend,&q,sizeof q,&code));
 lis_rci_projector_restore_destroy(p);lis_rci_split_restore_destroy(ctx);
 return 0;
}
int main(void)
{
 for(int general=0;general<2;general++)for(int m=LIS_RCI_BICGSTAB;m<=LIS_RCI_FGMRES;m++)C(!run(general,(lis_rci_method)m));
 puts("SPLIT PASS: six non-diagonal two-RHS contour projections, 36 shifted solves, callback failure/reuse");
 return 0;
}
