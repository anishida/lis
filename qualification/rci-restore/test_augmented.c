/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_augmented_restore.h"
#include "lis_rci_native_restore.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define N 20
#define C(e) do{if(!(e)){fprintf(stderr,"AUG FAIL %d: %s\n",__LINE__,#e);rc=1;goto cleanup;}}while(0)
static int run(int flex)
{
 LIS_MATRIX A=NULL;
 lis_rci_native_restore *native=NULL;
 lis_rci_augmented_restore *s=NULL;
 int rc=0;
 double exact[N],b[N],ax[N];
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 for(LIS_INT i=0;i<N;i++){
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,1+.4*i,A));
  if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.12,A));
  if(i+1<N)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,.2,A));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 native=lis_rci_native_restore_create(A,NULL,NULL);C(native);
 lis_rci_augmented_restore_options o={N,7,3,40,1e-11,0,1e-11,flex,lis_rci_restore_harmonic_native,NULL};
 s=lis_rci_augmented_restore_create(&o);C(s);
 lis_rci_recycle_request q;uint64_t last=0;
 for(int rep=0;rep<3;rep++){
  for(size_t i=0;i<N;i++)exact[i]=sin(.19*i+.3*rep)+.3*cos(.5*i);
  C(!lis_rci_native_restore_apply(native,LIS_RCI_OP_A,N,exact,b));
  int status=lis_rci_augmented_restore_start(s,b,NULL,1,1,&q);
  size_t count=0,pre=0,mat=0;
  while(status>0){
   C(++count<4000&&q.token>last);last=q.token;
   C(lis_rci_augmented_restore_advance(s,q.token-1,1,1,0,&q)==LIS_RCI_PROTOCOL);
   if(status==LIS_RCI_OP_A){C(!lis_rci_native_restore_apply(native,LIS_RCI_OP_A,N,q.in,q.out));mat++;}
   else {
    C(status==LIS_RCI_OP_PRECON);
    pre++;
    for(size_t i=0;i<N;i++)q.out[i]=q.in[i]/(1+.13*i+(flex?.09*((pre+i)%4):0));
   }
   status=lis_rci_augmented_restore_advance(s,q.token,1,1,0,&q);
  }
  lis_rci_augmented_restore_info info;C(!lis_rci_augmented_restore_get_info(s,&info));
  if(status)fprintf(stderr,"status=%d flex=%d cycles=%zu iterations=%zu residual=%Lg builds=%zu updates=%zu rejects=%zu\n",status,flex,info.cycles,info.iterations,info.residual,info.builds,info.updates,info.rejects);
  C(status==0&&info.status==0&&info.matvecs==mat&&info.precons==pre);
  C(info.builds+info.updates>0&&info.recycle_rank>0);
  if(rep)C(info.recycle_uses>0);
  const double *x=lis_rci_augmented_restore_solution(s);C(x);
  C(!lis_rci_native_restore_apply(native,LIS_RCI_OP_A,N,x,ax));
  for(size_t i=0;i<N;i++)C(fabs(x[i]-exact[i])<1e-8);
  long double res=0,bn=0;
  for(size_t i=0;i<N;i++){res=hypotl(res,(long double)b[i]-ax[i]);bn=hypotl(bn,b[i]);}
  C(res<=1.01e-11L*bn);
  lis_rci_recycle_basis basis;C(!lis_rci_augmented_restore_basis(s,&basis));
  for(size_t j=0;j<basis.count;j++){
   C(!lis_rci_native_restore_apply(native,LIS_RCI_OP_A,N,basis.physical+j*N,ax));
   for(size_t i=0;i<N;i++){
    C(fabs(ax[i]-basis.image[j*N+i])<1e-9);
    if(!flex)C(fabs(basis.physical[j*N+i]-basis.domain[j*N+i]/(1+.13*i))<1e-9);
   }
   for(size_t k=0;k<basis.count;k++){
    double d=0;for(size_t i=0;i<N;i++)d+=basis.image[j*N+i]*basis.image[k*N+i];
    C(fabs(d-(j==k))<1e-10);
   }
  }
  printf("AUG flex=%d rhs=%d cycles=%zu it=%zu A=%zu P=%zu builds=%zu updates=%zu rejects=%zu reuse=%zu residual=%Lg\n",flex,rep,info.cycles,info.iterations,mat,pre,info.builds,info.updates,info.rejects,info.recycle_uses,res/bn);
 }
 /* Failure before a cycle finishes must preserve the previously committed basis. */
 lis_rci_recycle_basis old,after;C(!lis_rci_augmented_restore_basis(s,&old));
 double saved[3*3*N];size_t len=old.count*N;
 memcpy(saved,old.domain,len*sizeof(double));memcpy(saved+len,old.physical,len*sizeof(double));memcpy(saved+2*len,old.image,len*sizeof(double));
 C(lis_rci_augmented_restore_start(s,b,NULL,1,1,&q)>0);
 C(lis_rci_augmented_restore_advance(s,q.token,1,1,73,&q)==LIS_RCI_CALLBACK);
 C(!lis_rci_augmented_restore_solution(s));
 C(!lis_rci_augmented_restore_basis(s,&after)&&after.count==old.count);
 C(!memcmp(saved,after.domain,len*sizeof(double))&&!memcmp(saved+len,after.physical,len*sizeof(double))&&!memcmp(saved+2*len,after.image,len*sizeof(double)));
 C(lis_rci_augmented_restore_start(s,b,NULL,1,1,&q)>0);
 C(lis_rci_augmented_restore_cancel(s,q.token)==LIS_RCI_CANCELLED);
 C(lis_rci_augmented_restore_start(s,b,NULL,1,1,&q)>0);
 C(lis_rci_augmented_restore_advance(s,q.token,2,1,0,&q)==LIS_RCI_PROTOCOL);
cleanup:
 lis_rci_augmented_restore_destroy(s);lis_rci_native_restore_destroy(native);if(A)lis_matrix_destroy(A);
 return rc;
}
static int harmonic(void)
{
 int rc=0;double g[9]={1,0,0,0,1,0,0,0,1},f[9]={.1,0,0,0,1,0,0,0,.5},coef[6];
 size_t count=0;C(!lis_rci_restore_harmonic_native(NULL,3,2,g,f,coef,&count));C(count==2);
 C(fabs(coef[1])>.999999&&fabs(coef[3+2])>.999999);
cleanup:return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 int rc=harmonic()|run(0)|run(1);lis_finalize();
 if(!rc)puts("AUGMENTED RESTORE PASS: six native sequence solves, fixed/flexible harmonic updates, preserved Arnoldi relations, failure rollback");
 return rc;
}
