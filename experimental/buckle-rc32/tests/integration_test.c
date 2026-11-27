#include "lis.h"
#include "lis_solver.h"
#include "lis_precon.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static LIS_MATRIX load(char *path){FILE *f=fopen(path,"rb"); if(!f)exit(11); int n,z; if(fread(&n,4,1,f)!=1||fread(&z,4,1,f)!=1)exit(12); LIS_INT *p=lis_malloc((n+1)*sizeof(*p),"row"),*j=lis_malloc(z*sizeof(*j),"col"); LIS_SCALAR *v=lis_malloc(z*sizeof(*v),"val"); if(sizeof(LIS_INT)!=4||sizeof(LIS_SCALAR)!=8)exit(13); if(fread(p,4,n+1,f)!=n+1||fread(j,4,z,f)!=z||fread(v,8,z,f)!=z)exit(14);fclose(f);LIS_MATRIX M;lis_matrix_create(0,&M);lis_matrix_set_size(M,0,n); if(lis_matrix_set_csr(z,p,j,v,M)||lis_matrix_assemble(M))exit(15);return M;}
extern LIS_INT lis_esolver_rebuild_precon(LIS_MATRIX,LIS_SOLVER,LIS_PRECON*);
int main(int argc,char **argv){
 if(argc<4)return 10;
 lis_initialize(&argc,&argv);LIS_SOLVER s;lis_solver_create(&s);if(lis_solver_set_optionC(s))return 2;
 int sequence=atoi(argv[1]);LIS_MATRIX A=NULL;LIS_PRECON p=NULL;
 if(sequence && lis_solver_set_sequence_precon_refresh_policy(s,LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN))return 8;
 for(int k=2;k<argc && argv[k][0]!='-';k++){
  LIS_MATRIX B=load(argv[k]),old=NULL;LIS_VECTOR b,x,r,t;LIS_INT rc;LIS_REAL rr,br;double start;
  if(sequence && A){if(A->n!=B->n||A->nnz!=B->nnz||memcmp(A->ptr,B->ptr,(A->n+1)*sizeof(LIS_INT))||memcmp(A->index,B->index,A->nnz*sizeof(LIS_INT)))return 9;memcpy(A->value,B->value,A->nnz*sizeof(LIS_SCALAR));lis_matrix_mark_values_changed(A);lis_matrix_destroy(B);}else{old=A;A=B;}
  lis_vector_duplicate(A,&b);lis_vector_duplicate(A,&x);lis_vector_duplicate(A,&r);lis_vector_duplicate(A,&t);
  for(int i=0;i<A->n;i++)t->value[i]=sin(.23*i)+cos(.11*i);
  if(getenv("LIS_ACTUAL_RHS")){FILE *f=fopen(getenv("LIS_ACTUAL_RHS"),"rb");if(!f||fseek(f,(long)A->n*8,SEEK_SET)||fread(b->value,8,A->n,f)!=A->n)return 7;fclose(f);}else lis_matvec(A,t,b);
  lis_vector_set_all(0,x);start=lis_wtime();
  if(sequence){rc=lis_solve_sequence(A,b,x,s);p=s->precon;}
  else {rc=lis_esolver_rebuild_precon(A,s,&p);if(!rc){s->precon=p;rc=lis_psolve_ilut_csr(s,b,x);}}
  if(old)lis_matrix_destroy(old);
  printf("RESULT case=%d sequence=%d rc=%d wall=%.9g\n",k-2,sequence,rc,lis_wtime()-start);fflush(stdout);if(rc||!p)return 3;
  void *user=NULL;int get=lis_precon_get_user_data(p,&user);int set=lis_precon_set_user_data(p,NULL);printf("ACCESS get=%d set=%d\n",get,set);if(getenv("LIS_EXPERIMENT_ILUTP_REUSE") && (!get||!set))return 11;
  lis_matvec(A,x,r);lis_vector_axpy(-1,b,r);lis_vector_nrm2(r,&rr);lis_vector_nrm2(b,&br);printf("RES case=%d true_res=%.17g\n",k-2,(double)(rr/br));fflush(stdout);if(!isfinite(rr)||rr/br>1e-8)return 6;
  lis_vector_destroy(b);lis_vector_destroy(x);lis_vector_destroy(r);lis_vector_destroy(t);
 }
 if(!sequence){s->precon=NULL;lis_precon_destroy(p);}lis_solver_destroy(s);lis_matrix_destroy(A);lis_finalize();return 0;
}
