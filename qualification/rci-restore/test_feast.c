/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_native_restore.h"
#include "lis_rci_feast_restore.h"
#include "lis_esolver.h"
#include <math.h>
#include <stdio.h>
#include <limits.h>
#define N 4
#define M 3
#define P 16
#define C(e) do{if(!(e)){fprintf(stderr,"FEAST FAIL %d: %s\n",__LINE__,#e);rc=1;goto cleanup;}}while(0)
static int dense(void *user,size_t n,double *h,double *ev,double *coef)
{
 (void)user;LIS_INT iter=0;LIS_REAL error=0;
 return (int)lis_eigen_projected_solve((LIS_INT)n,h,ev,coef,&iter,&error);
}
static int run(int general)
{
 LIS_MATRIX A=NULL,B=NULL;lis_rci_native_restore *native=NULL;
 lis_rci_split_restore *split=NULL;lis_rci_feast_restore *f=NULL;int rc=0;
 const double T[N][N]={{1.5,.5,0,0},{.5,1.5,0,0},{0,0,6.5,1.5},{0,0,1.5,6.5}};
 double diag[N]={1,1,1,1},initial[N*M]={1,.3,.7,.2,.4,1,-.2,.6,.3,.2,1,1};
 if(general)for(int i=0;i<N;i++)diag[i]=i+1;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 if(general){C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,N));}
 for(LIS_INT i=0;i<N;i++){
  if(B)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,diag[i],B));
  for(LIS_INT j=0;j<N;j++)if(T[i][j]!=0)
   C(!lis_matrix_set_value(LIS_INS_VALUE,i,j,T[i][j]*sqrt(diag[i]*diag[j]),A));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 if(B){C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));}
 native=lis_rci_native_restore_create(A,B,NULL);C(native);
 lis_rci_options child={2*N,100,2*N,1e-13,0,LIS_RCI_FGMRES};
 split=lis_rci_split_restore_create(N,general,&child,native,lis_rci_native_restore_apply,lis_rci_native_restore_generations);C(split);
 lis_rci_opaque_backend back;C(!lis_rci_split_restore_backend(split,&back));
 lis_rci_complex_coeff z[P],w[P];long double pi=acosl(-1.L);
 for(int k=0;k<P;k++){
  long double theta=2*pi*(k+.5L)/P,dr=1.1L*cosl(theta),di=1.1L*sinl(theta);
  z[k]=(lis_rci_complex_coeff){1.5L+dr,di};w[k]=(lis_rci_complex_coeff){dr/P,di/P};
 }
 lis_rci_feast_restore_options o={N,M,P,12,general,.5,2.5,1e-10,1e-10,1e-12,1e-9,z,w,dense,NULL};
 f=lis_rci_feast_restore_create(&o);C(f);
 uint64_t a=0,b=0,last=0;C(!lis_rci_native_restore_generations(native,&a,&b));
 lis_rci_opaque_request q;lis_rci_feast_restore_result r;
 for(int rep=0;rep<2;rep++){
  int status=lis_rci_feast_restore_start(f,initial,a,b,&q);
  C(status>0);C(lis_rci_feast_restore_get_result(f,&r)==LIS_RCI_STATE);
  int count=0,code=0;
  while(status>0){
   C(++count<2000&&q.token>last);last=q.token;
   C(lis_rci_feast_restore_advance(f,q.token-1,a,b,0,&q)==LIS_RCI_PROTOCOL);
   int d=lis_rci_opaque_dispatch(&back,&q,sizeof q,&code);
   if(d)fprintf(stderr,"dispatch=%d code=%d op=%u\n",d,code,q.operation);
   C(!d);
   status=lis_rci_feast_restore_advance(f,q.token,a,b,0,&q);
  }
  if(status)fprintf(stderr,"outer status=%d general=%d\n",status,general);
  C(status==0);C(!lis_rci_feast_restore_get_result(f,&r));C(r.count==2&&r.cycles>=2);
  C(fabs(r.values[0]-1)<1e-9&&fabs(r.values[1]-2)<1e-9);
  for(size_t j=0;j<r.count;j++){
   C(r.residuals[j]<1e-10);
   const double *x=r.vectors+j*N;double norm=0,err=0;
   for(int i=0;i<N;i++){double ax=0;for(int k=0;k<N;k++)ax+=T[i][k]*sqrt(diag[i]*diag[k])*x[k];
    err=hypot(err,ax-r.values[j]*diag[i]*x[i]);norm+=diag[i]*x[i]*x[i];}
   C(err<1e-8&&fabs(norm-1)<1e-8);
  }
 }
 C(lis_rci_feast_restore_start(f,initial,a,b,&q)>0);
 C(lis_rci_feast_restore_cancel(f,q.token-1)==LIS_RCI_PROTOCOL);
 C(lis_rci_feast_restore_cancel(f,q.token)==LIS_RCI_CANCELLED);
 C(lis_rci_feast_restore_get_result(f,&r)==LIS_RCI_STATE);
 C(lis_rci_feast_restore_start(f,initial,a,b,&q)>0);
 C(lis_rci_feast_restore_advance(f,q.token,a,b,73,&q)==LIS_RCI_CALLBACK);
 {int code=0;C(lis_rci_feast_restore_error(f,&code)==LIS_RCI_CALLBACK&&code==73);}
 C(lis_rci_feast_restore_start(f,initial,a,b,&q)>0);
 C(!lis_matrix_mark_values_changed(A));C(!lis_rci_native_restore_generations(native,&a,&b));
 C(lis_rci_feast_restore_advance(f,q.token,a,b,0,&q)==LIS_RCI_PROTOCOL);
cleanup:
 lis_rci_feast_restore_destroy(f);lis_rci_split_restore_destroy(split);lis_rci_native_restore_destroy(native);
 if(B)lis_matrix_destroy(B);if(A)lis_matrix_destroy(A);return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 int rc=run(0)|run(1);lis_finalize();
 if(!rc)puts("FEAST RESTORE PASS: 4 native standard/generalized interval solves, physical residuals, repeat, stale token/cancel/error/generation checks");
 return rc;
}
