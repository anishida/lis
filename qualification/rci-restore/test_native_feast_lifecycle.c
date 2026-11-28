#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_esolver.h"
#include "lis_rci_native_feast_lifecycle.h"
#include <stdio.h>
#define C(x) do{if(!(x)){fprintf(stderr,"FEAST LIFECYCLE FAIL %d: %s\n",__LINE__,#x);return 1;}}while(0)
static int test(int gen)
{
 LIS_MATRIX A,B=NULL;LIS_VECTOR x,y;LIS_ESOLVER e,r;LIS_SCALAR value,ref;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,4));
 if(gen){C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,4));}
 double ev[]={1,2,4,7};
 for(int i=0;i<4;i++){double b=gen?1+.2*i:1;C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b*ev[i],A));if(B)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b,B));}
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 if(B){C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));}
 C(!lis_vector_duplicate(A,&x));C(!lis_vector_duplicate(A,&y));
 C(!lis_esolver_create(&e));C(!lis_esolver_create(&r));
 char opts[]="-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16";
 C(!lis_esolver_set_option(opts,e));C(!lis_esolver_set_option(opts,r));
 lis_rci_native_feast_lifecycle *K=lis_rci_native_feast_lifecycle_create();C(K);
 uint64_t last=0;
 for(int rep=0;rep<3;rep++){
  C(!lis_vector_set_all(.125,x));C(!lis_vector_set_all(.125,y));value=-123;
  lis_rci_native_session_request q;
  int st=lis_rci_native_feast_lifecycle_start(K,A,B,x,&value,e,&q);
  if(st<=0){LIS_INT n;lis_rci_native_feast_lifecycle_status(K,&n);fprintf(stderr,"FEAST START gen=%d rep=%d st=%d native=%lld\n",gen,rep,st,(long long)n);}
  C(st>0&&q.token>last);
  if(rep==0){last=q.token;C(lis_rci_native_feast_lifecycle_cancel(K,last-1)==LIS_RCI_PROTOCOL);C(lis_rci_native_feast_lifecycle_cancel(K,last)==LIS_RCI_CANCELLED);C(e->nev==0&&value==-123);for(int i=0;i<4;i++)C(x->value[i]==.125);continue;}
  C(!lis_gesolve(A,B,y,&ref,r));C(r->retcode==LIS_SUCCESS);
  int calls=0,shifted=0;
  while(st>0){
   C(q.token>last);last=q.token;C(lis_rci_native_feast_lifecycle_advance(K,last-1,0,&q)==LIS_RCI_PROTOCOL);
   C(value==-123);for(int i=0;i<4;i++)C(x->value[i]==.125);
   LIS_INT rc;
   if(q.representation==LIS_NR_SPLIT){C(q.shift_sign==1);shifted++;lis_vector_set_all(0,q.out_real);lis_vector_set_all(0,q.out_imag);rc=lis_efeast_shifted_solve_split(q.A,q.B,q.zr,q.zi,q.in_real,q.in_imag,q.out_real,q.out_imag,q.solver);}
   else{C(q.representation==LIS_NR_VECTOR);rc=q.block?lis_matvec_block(q.A,q.count,q.in,q.out):lis_matvec(q.A,q.in[0],q.out[0]);}
   calls++;st=lis_rci_native_feast_lifecycle_advance(K,last,rc,&q);
  }
  LIS_INT n;C(lis_rci_native_feast_lifecycle_status(K,&n)==LIS_RCI_DONE&&n==LIS_SUCCESS&&st==LIS_RCI_DONE);
  C(shifted>0&&e->nev==2&&e->nev==r->nev&&value==ref);
  for(int i=0;i<4;i++)C(x->value[i]==y->value[i]);
  for(int i=0;i<e->nev;i++){C(e->evalue[i]==r->evalue[i]&&e->resid[i]==r->resid[i]);for(int j=0;j<4;j++)C(e->evector[i]->value[j]==r->evector[i]->value[j]);}
  printf("FEAST FULL LIFECYCLE PASS generalized=%d rep=%d requests=%d\n",gen,rep,calls);
 }
 lis_rci_native_feast_lifecycle_destroy(K);lis_esolver_destroy(e);lis_esolver_destroy(r);lis_vector_destroy(x);lis_vector_destroy(y);if(B)lis_matrix_destroy(B);lis_matrix_destroy(A);return 0;
}
int main(int argc,char **argv){C(!lis_initialize(&argc,&argv));C(!test(0));C(!test(1));lis_finalize();return 0;}
