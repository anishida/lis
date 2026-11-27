#include "lis.h"
#include "lis_solver.h"
#include "lis_precon.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static LIS_MATRIX load(char *path){FILE *f=fopen(path,"rb"); if(!f)exit(11); int n,z; if(fread(&n,4,1,f)!=1||fread(&z,4,1,f)!=1)exit(12); LIS_INT *p=lis_malloc((n+1)*sizeof(*p),"row"),*j=lis_malloc(z*sizeof(*j),"col"); LIS_SCALAR *v=lis_malloc(z*sizeof(*v),"val"); if(sizeof(LIS_INT)!=4||sizeof(LIS_SCALAR)!=8)exit(13); if(fread(p,4,n+1,f)!=n+1||fread(j,4,z,f)!=z||fread(v,8,z,f)!=z)exit(14);fclose(f);LIS_MATRIX M;lis_matrix_create(0,&M);lis_matrix_set_size(M,0,n); if(lis_matrix_set_csr(z,p,j,v,M)||lis_matrix_assemble(M))exit(15);return M;}
typedef struct LIS_ILUTP_REUSE_STRUCT *LIS_ILUTP_REUSE;
extern LIS_INT lis_ilutp_reuse_create(LIS_SOLVER,LIS_MATRIX,LIS_ILUTP_REUSE*);
extern LIS_INT lis_ilutp_reuse_refresh(LIS_ILUTP_REUSE,LIS_SOLVER,LIS_MATRIX);
extern LIS_PRECON lis_ilutp_reuse_factors(LIS_ILUTP_REUSE);
extern void lis_ilutp_reuse_destroy(LIS_ILUTP_REUSE);
int main(int argc,char **argv){
 if(argc<4)return 10;
 lis_initialize(&argc,&argv);LIS_SOLVER s;lis_solver_create(&s);if(lis_solver_set_optionC(s))return 2;
 LIS_ILUTP_REUSE c=NULL;LIS_PRECON p=NULL;int fresh=atoi(argv[1]);int negative=getenv("NEGATIVE")?atoi(getenv("NEGATIVE")):0;
 for(int k=2;k<argc && argv[k][0]!='-';k++){
  LIS_MATRIX A=load(argv[k]);LIS_VECTOR b,x,r,t;LIS_INT rc;LIS_REAL rr,br;double start;
  s->A=A;lis_vector_duplicate(A,&b);lis_vector_duplicate(A,&x);lis_vector_duplicate(A,&r);lis_vector_duplicate(A,&t);
  for(int i=0;i<A->n;i++)t->value[i]=sin(.23*i)+cos(.11*i);
  if(k==3 && negative==1){ A->index[0]=(A->index[0]+1)%A->n; }
  if(k==3 && negative==2){for(int j=0;j<A->nnz;j++)A->value[j]=0;}
  if(k==3 && negative==3){A->value[0]=NAN;}
  if(k==3 && negative==4){s->ilutp_full_limit=1;}
  start=lis_wtime();
  if(fresh){rc=lis_precon_create(s,&p);if(rc)p=NULL;}
  else {rc=c?lis_ilutp_reuse_refresh(c,s,A):lis_ilutp_reuse_create(s,A,&c);p=lis_ilutp_reuse_factors(c);}
  printf("SETUP case=%d fresh=%d rc=%d factors=%d wall=%.9g\n",k-2,fresh,rc,p!=NULL,lis_wtime()-start);fflush(stdout);
  if(k==3 && negative){if(!rc||p)return 5;lis_vector_destroy(b);lis_vector_destroy(x);lis_vector_destroy(r);lis_vector_destroy(t);lis_matrix_destroy(A);s->A=NULL;break;}
  if(rc||!p)return 3;
  s->precon=p;
  for(int adj=0;adj<2;adj++){
   if(getenv("LIS_ACTUAL_RHS") && !adj){FILE *f=fopen(getenv("LIS_ACTUAL_RHS"),"rb");if(!f||fseek(f,(long)A->n*8,SEEK_SET)||fread(b->value,8,A->n,f)!=A->n)return 7;fclose(f);}else if(adj)lis_matvech(A,t,b);else lis_matvec(A,t,b);
   rc=adj?lis_psolveh_ilut_csr(s,b,x):lis_psolve_ilut_csr(s,b,x);if(rc)return 4;
   if(adj)lis_matvech(A,x,r);else lis_matvec(A,x,r);lis_vector_axpy(-1,b,r);lis_vector_nrm2(r,&rr);lis_vector_nrm2(b,&br);
   printf("RES case=%d adj=%d true_res=%.17g\n",k-2,adj,(double)(rr/br));fflush(stdout);if(!isfinite(rr)||rr/br>1e-8)return 6;
  }
  s->precon=NULL;s->A=NULL;if(fresh){lis_precon_destroy(p);p=NULL;}
  lis_vector_destroy(b);lis_vector_destroy(x);lis_vector_destroy(r);lis_vector_destroy(t);lis_matrix_destroy(A);
 }
 lis_ilutp_reuse_destroy(c);lis_solver_destroy(s);lis_finalize();return 0;
}
