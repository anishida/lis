/* Installed public API: independent residual and classic parity on nontrivial systems. */
#define LIS_USE_INSTALLED_CONFIG
#include <lis_rci_native.h>
#include <math.h>
#include <stdio.h>
#define C(x) do{if(!(x)){fprintf(stderr,"RELEASE FAIL %d: %s\n",__LINE__,#x);return 1;}}while(0)
static lis_rci_native_context *context(int kind)
{
 LIS_ABI_INFO a;lis_abi_get_client_info(&a);
 return lis_rci_native_create(kind,LIS_RCI_NATIVE_API_VERSION,sizeof(lis_rci_native_session_request),&a,sizeof a);
}
static int run(lis_rci_native_context *s,int st,lis_rci_native_session_request q)
{
 uint64_t last=0;
 while(st>0){
  C(q.token>last);last=q.token;
  LIS_INT rc=lis_rci_native_execute_default(s,last);
  C(rc==LIS_SUCCESS);st=lis_rci_native_advance(s,last,rc,&q);
 }
 LIS_INT native;C(st==LIS_RCI_DONE&&lis_rci_native_status(s,&native)==st&&native==LIS_SUCCESS);
 return 0;
}
static int linear(void)
{
 enum{N=24};
 LIS_MATRIX A;LIS_VECTOR b,x,y,z;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 for(int i=0;i<N;i++){
  LIS_SCALAR lower=-.7,upper=-.7;
#ifdef _COMPLEX
  lower-=.15*I;upper+=.15*I;
#endif
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,3+.01*i,A));
  if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,lower,A));
  if(i<N-1)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,upper,A));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));C(!lis_vector_duplicate(A,&y));C(!lis_vector_duplicate(A,&z));
 const char *m[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 for(int i=0;i<6;i++){
  LIS_SOLVER s,r;C(!lis_solver_create(&s));C(!lis_solver_create(&r));char opts[240];
  snprintf(opts,sizeof opts,"-i %s -p jacobi -restart 8 -gcrodr_keep 2 -fgcrodr_k 2 -tol 1e-10 -maxiter 300 -initx_zeros false",m[i]);
  C(!lis_solver_set_option(opts,s));C(!lis_solver_set_option(opts,r));
  lis_rci_native_context *K=context(LIS_NR_LINEAR);C(K);
  for(int rep=0;rep<3;rep++){
   LIS_SCALAR exact[N];
   for(int j=0;j<N;j++){
    exact[j]=1+.02*j+.05*rep;
#ifdef _COMPLEX
    exact[j]+=(.13+.007*j)*I;
#endif
    C(!lis_vector_set_value(LIS_INS_VALUE,j,exact[j],z));
   }
   C(!lis_matvec(A,z,b));C(!lis_vector_set_all(.125,x));C(!lis_vector_set_all(.125,y));
   C(!lis_solve(A,b,y,r));C(r->retcode==LIS_SUCCESS);
   lis_rci_native_session_request q;int st=lis_rci_native_start_linear(K,A,b,x,s,&q);
   C(st>0);C(!run(K,st,q));C(s->iter==r->iter&&s->resid==r->resid);
   LIS_SCALAR xx[N],yy[N];C(!lis_vector_get_values(x,0,N,xx));C(!lis_vector_get_values(y,0,N,yy));
   for(int j=0;j<N;j++){C(xx[j]==yy[j]);C(fabs(xx[j]-exact[j])<1e-8);}
   C(!lis_matvec(A,x,z));C(!lis_vector_axpy(-1,b,z));LIS_REAL residual,rhs;
   C(!lis_vector_nrm2(z,&residual));C(!lis_vector_nrm2(b,&rhs));C(residual/rhs<1e-9);
  }
  lis_rci_native_destroy(K);lis_solver_destroy(s);lis_solver_destroy(r);printf("RELEASE LINEAR PASS %s\n",m[i]);
 }
 lis_vector_destroy(b);lis_vector_destroy(x);lis_vector_destroy(y);lis_vector_destroy(z);lis_matrix_destroy(A);return 0;
}
static int feast(void)
{
 LIS_MATRIX A,B;LIS_VECTOR x,y;LIS_ESOLVER e,r;LIS_SCALAR v,w;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_create(LIS_COMM_WORLD,&B));
 C(!lis_matrix_set_size(A,0,4));C(!lis_matrix_set_size(B,0,4));double ev[]={1,2,4,7};
 for(int i=0;i<4;i++){double b=1+.2*i;C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b*ev[i],A));C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b,B));}
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));C(!lis_matrix_assemble(B));
 C(!lis_vector_duplicate(A,&x));C(!lis_vector_duplicate(A,&y));C(!lis_esolver_create(&e));C(!lis_esolver_create(&r));
 char opts[]="-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16";
 C(!lis_esolver_set_option(opts,e));C(!lis_esolver_set_option(opts,r));
 lis_rci_native_context *K=context(LIS_NR_FEAST);C(K);
 for(int rep=0;rep<2;rep++){
  C(!lis_vector_set_all(1,x));C(!lis_vector_set_all(1,y));v=-123;
  C(!lis_gesolve(A,B,y,&w,r));C(r->retcode==LIS_SUCCESS);
  lis_rci_native_session_request q;int st=lis_rci_native_start_feast(K,A,B,x,&v,e,&q);C(st>0);
  if(!rep){++B->value_generation;C(lis_rci_native_advance(K,q.token,0,&q)==LIS_NR_STALE&&v==-123&&e->nev==0);
   st=lis_rci_native_start_feast(K,A,B,x,&v,e,&q);C(st>0);}
  C(!run(K,st,q));C(e->nev==2&&e->nev==r->nev&&v==w);
  int found=0;
  for(int i=0;i<2;i++){
   C(e->evalue[i]==r->evalue[i]&&e->resid[i]==r->resid[i]);
   if(fabs(e->evalue[i]-2)<1e-6)found|=1;else if(fabs(e->evalue[i]-4)<1e-6)found|=2;else C(0);
   for(int j=0;j<4;j++)C(e->evector[i]->value[j]==r->evector[i]->value[j]);
  }
  C(found==3);printf("RELEASE FEAST GENERALIZED PASS rep=%d\n",rep);
 }
 lis_rci_native_destroy(K);lis_esolver_destroy(e);lis_esolver_destroy(r);lis_vector_destroy(x);lis_vector_destroy(y);lis_matrix_destroy(A);lis_matrix_destroy(B);return 0;
}
int main(int argc,char **argv){C(!lis_initialize(&argc,&argv));C(!linear());C(!feast());lis_finalize();return 0;}
