/* Differential native full FEAST qualification. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_esolver.h"
#include "lis_rci_native_feast.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define N 4
#define C(x) do{if(!(x)){fprintf(stderr,"FEAST FULL FAIL %d: %s\n",__LINE__,#x);rc=1;goto cleanup;}}while(0)
typedef struct {LIS_INT nev,iter[4],iter2[4];LIS_REAL value[4],resid[4];LIS_SCALAR vectors[16];} answer;
static int equal_answer(const answer *a,const answer *b)
{
 if(a->nev!=b->nev)return 0;
 for(size_t i=0;i<4;i++)if(a->iter[i]!=b->iter[i])return 0;
 for(size_t i=0;i<4;i++)if(a->iter2[i]!=b->iter2[i])return 0;
 for(size_t i=0;i<4;i++)if(a->value[i]!=b->value[i])return 0;
 for(size_t i=0;i<4;i++)if(a->resid[i]!=b->resid[i])return 0;
 for(size_t i=0;i<16;i++)if(a->vectors[i]!=b->vectors[i])return 0;
 return 1;
}

static int stop,kind,observed,terminal;
static int first[5],last_request[5],seen[5];
static LIS_INT external(LIS_ESOLVER e)
{
 lis_rci_native_feast *K=lis_rci_native_feast_create();if(!K)return LIS_OUT_OF_MEMORY;
 lis_rci_native_feast_request q;
 int st=lis_rci_native_feast_start(K,e,1,&q);uint64_t last=0;observed=0;
 while(st>0){
  if(q.token<=last||lis_rci_native_feast_advance(K,q.token-1,0,&q)!=NF_PROTOCOL){st=NF_ERROR;break;}
  last=q.token;
  ++observed;
  int bucket=q.operation==NF_SHIFT?4:q.operation==NF_B?(q.block?3:1):(q.block?2:0);
  if(!stop){if(!first[bucket])first[bucket]=observed;last_request[bucket]=observed;seen[bucket]++;}
  if(observed==stop){
   if(kind==1)st=lis_rci_native_feast_cancel(K,q.token);
   else if(kind==2)st=lis_rci_native_feast_advance(K,q.token,0,&q);
   else if(kind==3)st=lis_rci_native_feast_advance(K,q.token,LIS_ERR_ILL_ARG,&q);
   else {
    /* Test-only counter injection, without changing numerical values. */
    if(kind==4)e->A->value_generation++;
    else if(kind==5)e->A->structure_generation++;
    else e->B->value_generation++;
    st=lis_rci_native_feast_advance(K,q.token,0,&q);
   }
   break;
  }
  LIS_INT err;
  if(q.operation==NF_SHIFT){
   lis_vector_set_all(0,q.out_real);lis_vector_set_all(0,q.out_imag);
   err=lis_efeast_shifted_solve_split(q.A,q.B,q.zr,q.zi,q.in_real,q.in_imag,q.out_real,q.out_imag,q.solver);
  }else err=!q.block?lis_matvec(q.A,q.in[0],q.out[0]):lis_matvec_block(q.A,q.count,q.in,q.out);
  st=lis_rci_native_feast_advance(K,q.token,err,&q);
 }
 terminal=st;
 LIS_INT rc=lis_rci_native_feast_native(K);lis_rci_native_feast_destroy(K);return rc;
}
static int solve(int mode,int generalized,int fail,int fault,answer *a)
{
 LIS_MATRIX A=NULL,B=NULL;LIS_VECTOR x=NULL;LIS_ESOLVER e=NULL;
 int rc=0;LIS_SCALAR value=0;memset(a,0,sizeof(*a));
 const double lambda[N]={1,2,4,7};
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 if(generalized){C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,N));}
 for(LIS_INT i=0;i<N;i++){
  double b=generalized?1+.2*i:1;
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b*lambda[i],A));
  if(B)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b,B));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 if(B){C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));}
 C(!lis_vector_duplicate(A,&x));C(!lis_vector_set_all(1,x));
 C(!lis_esolver_create(&e));
 C(!lis_esolver_set_option("-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16",e));
 stop=fail;kind=fault;terminal=0;observed=0;
 lis_rci_native_feast_test_select(mode==0,mode==2?external:NULL);
 LIS_INT err=B?lis_gesolve(A,B,x,&value,e):lis_esolve(A,x,&value,e);
 if(fail){
  /* lis_gesolve reports dispatch/setup errors directly; numerical status
   * is stored in esolver->retcode even when its public return is SUCCESS. */
  if(e->retcode==LIS_SUCCESS||e->nev!=0)fprintf(stderr,"FEAST fault gen=%d site=%d kind=%d api=%lld native=%lld nev=%lld terminal=%d\n",generalized,fail,fault,(long long)err,(long long)e->retcode,(long long)e->nev,terminal);
  C(e->retcode==(fault==3?LIS_ERR_ILL_ARG:LIS_BREAKDOWN)&&e->nev==0);
  C(terminal==(fault==1?NF_CANCELLED:fault==2?NF_NONFINITE:fault==3?NF_ERROR:NF_STALE));goto cleanup;
 }
 C(!err&&e->retcode==LIS_SUCCESS);C(e->nev==2);
 a->nev=e->nev;
 for(LIS_INT i=0;i<e->nev;i++){
  #ifdef _COMPLEX
  C(fabs(cimag(e->evalue[i]))<1e-7);
  a->value[i]=(LIS_REAL)creal(e->evalue[i]);
#else
  a->value[i]=e->evalue[i];
#endif
  a->resid[i]=e->resid[i];a->iter[i]=e->iter[i];a->iter2[i]=e->iter2[i];
  C(e->evector[i]);memcpy(a->vectors+N*i,e->evector[i]->value,N*sizeof(LIS_SCALAR));
  printf("FEAST pair gen=%d mode=%d index=%lld value=%.17g residual=%.17g\n",generalized,mode,(long long)i,(double)a->value[i],(double)a->resid[i]);
  C((fabs(a->value[i]-2)<1e-6||fabs(a->value[i]-4)<1e-6)&&a->resid[i]<=1e-7);
  /* Independent physical residual, not an echo of LIS certification. */
  double ax2=0,bx2=0,r2=0;
  for(int j=0;j<N;j++){
   double bj=generalized?1+.2*j:1;
   LIS_SCALAR bx=bj*a->vectors[N*i+j],ax=lambda[j]*bx;
   double abx=(double)fabs(bx),aax=(double)fabs(ax);
   ax2+=aax*aax;bx2+=abx*abx;
   double delta=(double)fabs(ax-a->value[i]*bx);r2+=delta*delta;
  }
  C((double)fabs(sqrt(r2))/( (double)fabs(sqrt(ax2))+(double)fabs(a->value[i])*(double)fabs(sqrt(bx2)))<1e-7);
 }
 C((fabs(a->value[0]-2)<1e-6&&fabs(a->value[1]-4)<1e-6)||(fabs(a->value[0]-4)<1e-6&&fabs(a->value[1]-2)<1e-6));
cleanup:
 lis_rci_native_feast_test_select(0,NULL);
 if(e)lis_esolver_destroy(e);if(x)lis_vector_destroy(x);if(B)lis_matrix_destroy(B);if(A)lis_matrix_destroy(A);return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 int rc=0;answer a,b,c,scratch;
 for(int gen=0;gen<=1;gen++){
  memset(first,0,sizeof(first));memset(last_request,0,sizeof(last_request));memset(seen,0,sizeof(seen));
  C(!solve(0,gen,0,0,&a));C(!solve(1,gen,0,0,&b));C(!solve(2,gen,0,0,&c));
  C(equal_answer(&a,&b)&&equal_answer(&a,&c));
  int n=observed;C(n>0);
  C(seen[0]>0&&seen[2]>0&&seen[4]>0);
  C(gen?(seen[1]>0&&seen[3]>0):(!seen[1]&&!seen[3]));
  int faults=0;
  /* First and last of every request type, plus first 32 boundaries. */
  for(int f=1;f<=n;f++){
   int selected=f<=32||f==n;
   for(int j=0;j<5;j++)selected|=f==first[j]||f==last_request[j];
   if(selected)for(int k=1;k<=3;k++){C(!solve(2,gen,f,k,&scratch));faults++;}
  }
  for(int k=4;k<=(gen?6:5);k++){C(!solve(2,gen,1,k,&scratch));C(!solve(2,gen,n,k,&scratch));faults+=2;}
  printf("NATIVE FULL FEAST PASS generalized=%d requests=%d faults=%d A=%d B=%d ABLOCK=%d BBLOCK=%d SHIFT=%d\n",gen,n,faults,seen[0],seen[1],seen[2],seen[3],seen[4]);
 }
cleanup:
 lis_finalize();return rc;
}
