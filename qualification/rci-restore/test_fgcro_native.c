/* AI-assisted continuation reconstruction; proposed under LIS terms in
 * COPYING, subject to contributor confirmation. Qualification only. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_native_fgcro_restore.h"
#include "lislib.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define N 48
extern LIS_INT lis_fgcrodr_validate_recycle(LIS_SOLVER,LIS_REAL *,LIS_REAL *);
#define C(e) do { if (!(e)) { fprintf(stderr,"FGCRO FAIL %d: %s\n",__LINE__,#e); rc=1; goto cleanup; } } while(0)
typedef struct {
 LIS_INT native,iter,dim,action;
 LIS_REAL residual,quality,minquality,orth,relation;
 LIS_SCALAR x[N],basis[3*3*N];
} answer;
static int equal_answer(const answer *a,const answer *b)
{
 if(a->native!=b->native||a->iter!=b->iter||a->dim!=b->dim||a->action!=b->action||a->residual!=b->residual||a->quality!=b->quality||a->minquality!=b->minquality||a->orth!=b->orth||a->relation!=b->relation)return 0;
 for(size_t i=0;i<N;i++)if(a->x[i]!=b->x[i])return 0;
 for(size_t i=0;i<3*3*N;i++)if(a->basis[i]!=b->basis[i])return 0;
 return 1;
}

static uint64_t stop, observed;
static int fault, terminal, stale;
static LIS_INT external(LIS_SOLVER s)
{
 lis_rci_native_fgcro_restore *k=lis_rci_native_fgcro_restore_create();
 if(!k)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_linear_request q;
 int st=lis_rci_native_fgcro_restore_start(k,s,1,&q);
 uint64_t last=0;
 observed=0;stale=0;
 while(st>0){
  if(q.token<=last || lis_rci_native_fgcro_restore_advance(k,q.token-1,0,&q)!=LIS_RCI_PROTOCOL){stale=1;break;}
  last=q.token;
  if(++observed==stop){
   if(fault==1)st=lis_rci_native_fgcro_restore_cancel(k,q.token);
   else st=lis_rci_native_fgcro_restore_advance(k,q.token,0,&q);
  } else {
   LIS_INT e=q.operation==LIS_RCI_OP_A?lis_matvec(q.A,q.in,q.out):lis_psolve(q.solver,q.in,q.out);
   st=lis_rci_native_fgcro_restore_advance(k,q.token,e,&q);
  }
 }
 terminal=st;
 LIS_INT native=0;uint64_t a,p;
 if(lis_rci_native_fgcro_restore_status(k,&native,&a,&p)!=st || stale)native=LIS_BREAKDOWN;
 lis_rci_native_fgcro_restore_destroy(k);
 return native;
}
static int sequence(int mode,answer out[3],int fail_rep,uint64_t fail_at,int kind,uint64_t counts[3])
{
 LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL;LIS_SOLVER s=NULL;
 int rc=0;
 memset(out,0,3*sizeof(*out));
 C(!lis_solver_create(&s));
 C(!lis_solver_set_option("-i fgcrodr -p jacobi -restart 8 -fgcrodr_k 3 -maxiter 300 -tol 1e-11 -initx_zeros false",s));
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));
 C(!lis_matrix_set_size(A,0,N));
 for(LIS_INT i=0;i<N;i++){
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,2.6+.008*i,A));
  if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.7,A));
  if(i+1<N)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.1,A));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));
 lis_rci_native_fgcro_test_select(mode==0,mode==2?external:NULL);
 for(int rep=0;rep<3;rep++){
  LIS_SCALAR exact[N];
  for(int i=0;i<N;i++)exact[i]=sin(.19*i+.1*rep)+.4;
#ifdef _COMPLEX
  for(int j=0;j<N;j++)exact[j]+=I*(.15+.007*j);
#endif
  for(LIS_INT i=0;i<N;i++){
   LIS_SCALAR rhs=(2.6+.008*i)*exact[i]-(i?.7*exact[i-1]:0)-(i+1<N?1.1*exact[i+1]:0);
   C(!lis_vector_set_value(LIS_INS_VALUE,i,rhs,b));C(!lis_vector_set_value(LIS_INS_VALUE,i,.1,x));
  }
  stop=rep==fail_rep?fail_at:0;fault=kind;terminal=0;stale=0;observed=0;
  (void)lis_solve(A,b,x,s);
  answer *a=out+rep;
  C(!lis_solver_get_status(s,&a->native));
  C(!lis_solver_get_fgcrodr_recycle(s,&a->dim,&a->action,&a->quality,&a->minquality));
  if(stop){
   C(terminal==(kind==1?LIS_RCI_CANCELLED:LIS_RCI_NONFINITE));
   C(a->native==LIS_BREAKDOWN&&a->dim==0);
   goto cleanup;
  }
  C(!stale&&a->native==LIS_SUCCESS);
  C(!lis_solver_get_iter(s,&a->iter));C(!lis_solver_get_residualnorm(s,&a->residual));
  C(a->dim>0&&a->dim<=3);
  C(!lis_fgcrodr_validate_recycle(s,&a->orth,&a->relation));
  if(!(a->orth<1e-8&&a->relation<1e-8))fprintf(stderr,"FGCRO RECYCLE mode=%d rep=%d orth=%.17g relation=%.17g\n",mode,rep,(double)a->orth,(double)a->relation);
  C(a->orth<1e-8&&a->relation<1e-8);
  C(!lis_rci_native_fgcro_test_basis(s,a->basis,3*3*N));
  C(!lis_vector_get_values(x,0,N,a->x));
  for(int i=0;i<N;i++)C(fabs(a->x[i]-exact[i])<1e-8);
  if(counts)counts[rep]=observed;
 }
cleanup:
 lis_rci_native_fgcro_test_select(0,NULL);
 if(s)lis_solver_destroy(s);
 if(x)lis_vector_destroy(x);
 if(b)lis_vector_destroy(b);
 if(A)lis_matrix_destroy(A);
 return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 int rc=0;answer a[3],b[3],c[3],scratch[3];uint64_t counts[3],faults=0;
 C(!sequence(0,a,-1,0,0,NULL));
 C(!sequence(1,b,-1,0,0,NULL));
 C(!sequence(2,c,-1,0,0,counts));
 for(int i=0;i<3;i++){
  C(counts[i]>0);
  C(equal_answer(a+i,b+i)&&equal_answer(a+i,c+i));
 }
 for(int rep=0;rep<2;rep++)for(uint64_t f=1;f<=counts[rep];f++)
  for(int kind=1;kind<=2;kind++){C(!sequence(2,scratch,rep,f,kind,NULL));faults++;}
 printf("NATIVE FGCRO PASS: three-way sequence; %llu boundary faults\n",(unsigned long long)faults);
cleanup:
 lis_finalize();return rc;
}
