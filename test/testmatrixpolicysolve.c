#include "lis.h"
#include "lis_matrix_materialize.h"
#include <math.h>
#include <stdio.h>
typedef struct { int dense; LIS_REAL scale; } CTX;
static LIS_INT row_cb(void *ptr,LIS_INT r,LIS_INT *nnz,LIS_INT *cols,LIS_SCALAR *vals) {
  const double m[9]={4,-1,0.25,-1,4,-1,0.25,-1,3};
  CTX *c=(CTX*)ptr; LIS_INT j,n=0,cap;
  if(!c||!nnz||r<0||r>=3)return LIS_ERR_ILL_ARG;
  cap=*nnz;
  if((cols==NULL)!=(vals==NULL))return LIS_ERR_ILL_ARG;
  for(j=0;j<3;j++) {
    if(!c->dense&&((r==0&&j==2)||(r==2&&j==0)))continue;
    if(m[r*3+j]==0)continue;
    if(cols){if(n>=cap)return LIS_ERR_ILL_ARG;cols[n]=j;vals[n]=(LIS_SCALAR)(m[r*3+j]*c->scale);}
    n++;
  }
  *nnz=n;return LIS_SUCCESS;
}
static LIS_INT mv(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y) {
  LIS_INT i,j,n,idx[3],err;LIS_SCALAR val[3];
  if(!ptr||!x||!y)return LIS_ERR_ILL_ARG;
  for(i=0;i<3;i++){n=3;err=row_cb(ptr,i,&n,idx,val);if(err)return err;y[i]=0;for(j=0;j<n;j++)y[i]+=val[j]*x[idx[j]];}
  return LIS_SUCCESS;
}
static int solve(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,const char *pc) {
  LIS_SOLVER s=NULL;LIS_INT e,status;int fail=1;char opts[128];
  snprintf(opts,sizeof(opts),"-i gmres -p %s -restart 10 -tol 1e-11 -maxiter 80 -scale none",pc);
  e=lis_solver_create(&s);if(e)goto done;
  e=lis_solver_set_option(opts,s);if(e)goto done;
  e=lis_vector_set_all((LIS_SCALAR)0,x);if(e)goto done;
  e=lis_solve(A,b,x,s);if(e)goto done;
  e=lis_solver_get_status(s,&status);if(e||status!=LIS_SUCCESS)goto done;
  fail=0;
done:
  if(fail)fprintf(stderr,"solve %s failed: %lld\n",pc,(long long)e);
  if(s)lis_solver_destroy(s);return fail;
}
static int check(LIS_VECTOR x,LIS_VECTOR truth,const char *label) {
  LIS_VECTOR diff=NULL;LIS_REAL norm;LIS_INT e;int fail=1;
  e=lis_vector_duplicate(x,&diff);if(e)goto done;
  e=lis_vector_copy(x,diff);if(e)goto done;
  e=lis_vector_axpy((LIS_SCALAR)-1,truth,diff);if(e)goto done;
  e=lis_vector_nrm2(diff,&norm);if(e)goto done;
  printf("%s error=%.15g\n",label,(double)norm);
  fail= !(isfinite((double)norm)&&norm<(LIS_REAL)1e-8);
done:
  if(diff)lis_vector_destroy(diff);return fail;
}
int main(int argc,char **argv) {
 CTX ctx={0,(LIS_REAL)1};LIS_MATRIX A=NULL,sel=NULL;
 LIS_VECTOR truth=NULL,b=NULL,xd=NULL,xc=NULL,xs=NULL;
 LIS_SOLVER ilut=NULL,ssor=NULL;LIS_MATRIX_CSR_CACHE cache=NULL;
 LIS_INT e,policy,action,expect[3]={2,1,2};int step,fail=1;
 e=lis_initialize(&argc,&argv);if(e)return 1;
 e=lis_matrix_create(LIS_COMM_WORLD,&A);if(e)goto done;
 e=lis_matrix_set_size(A,0,3);if(e)goto done;
 e=lis_matrix_set_user(A,&ctx,mv,mv);if(e)goto done;
 e=lis_matrix_set_user_row(A,row_cb);if(e)goto done;
 e=lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);if(e)goto done;
 e=lis_vector_duplicate(A,&truth);if(e)goto done;
 e=lis_vector_duplicate(A,&b);if(e)goto done;
 e=lis_vector_duplicate(A,&xd);if(e)goto done;
 e=lis_vector_duplicate(A,&xc);if(e)goto done;
 e=lis_vector_duplicate(A,&xs);if(e)goto done;
 truth->value[0]=1;truth->value[1]=-2;truth->value[2]=3;
 e=lis_solver_create(&ilut);if(e)goto done;
 e=lis_solver_create(&ssor);if(e)goto done;
 {char d[]="-i gmres -p ilut",s[]="-i gmres -p ssor";
  e=lis_solver_set_option(d,ilut);if(e)goto done;
  e=lis_solver_set_option(s,ssor);if(e)goto done;}
 for(step=0;step<3;step++){
   if(step==1){ctx.scale=(LIS_REAL)1.5;e=lis_matrix_mark_values_changed(A);if(e)goto done;}
   if(step==2){ctx.dense=1;e=lis_matrix_mark_structure_changed(A);if(e)goto done;}
   e=lis_matvec(A,truth,b);if(e)goto done;
   e=lis_matrix_prepare_solver_matrix(A,ilut,&cache,&sel,&policy,&action);if(e)goto done;
   if(sel!=A||policy!=LIS_MATRIX_POLICY_DIRECT||action!=LIS_MATRIX_POLICY_ACTION_DIRECT)goto done;
   if(solve(sel,b,xd,"ilut"))goto done;
   e=lis_matrix_prepare_solver_matrix(A,ssor,&cache,&sel,&policy,&action);if(e)goto done;
   if(sel==A||policy!=LIS_MATRIX_POLICY_CACHED_CSR||action!=expect[step]){fprintf(stderr,"bad cache action %lld at %d\n",(long long)action,step);goto done;}
   {
     LIS_VECTOR tmp=NULL;
     LIS_REAL mismatch;
     e=lis_vector_duplicate(A,&tmp);if(e)goto done;
     e=lis_matvec(sel,truth,tmp);if(e){lis_vector_destroy(tmp);goto done;}
     e=lis_vector_axpy((LIS_SCALAR)-1,b,tmp);
     if(e){lis_vector_destroy(tmp);goto done;}
     e=lis_vector_nrm2(tmp,&mismatch);
     lis_vector_destroy(tmp);
     if(e)goto done;
     printf("STEP %d CSR_MATVEC_RHS_DIFF=%.15g\n",step,(double)mismatch);
   }
   if(solve(sel,b,xc,"ilut")||solve(sel,b,xs,"ssor"))goto done;
   if(check(xd,truth,"direct USER/ILUT")||check(xc,truth,"cached CSR/ILUT")||check(xs,truth,"cached CSR/SSOR"))goto done;
   printf("MATRIX_POLICY STEP %d PASS\n",step);
 }
 fail=0;
done:
 printf("MATRIX_POLICY_NUMERICAL=%s\n",fail?"FAIL":"PASS");
 if(ilut)lis_solver_destroy(ilut);if(ssor)lis_solver_destroy(ssor);
 if(truth)lis_vector_destroy(truth);if(b)lis_vector_destroy(b);
 if(xd)lis_vector_destroy(xd);if(xc)lis_vector_destroy(xc);if(xs)lis_vector_destroy(xs);
 lis_matrix_csr_cache_destroy(cache);if(A)lis_matrix_destroy(A);lis_finalize();return fail;
}
