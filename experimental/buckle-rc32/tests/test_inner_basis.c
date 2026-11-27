#include "lis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define CHECK(c) do { if(!(c)){fprintf(stderr,"FAIL line=%d: %s\n",__LINE__,#c);return 1;} }while(0)
#define OK(c) CHECK((c)==LIS_SUCCESS)
static int calls=0,errors=0;
extern LIS_INT __real_lis_precon_create(LIS_SOLVER,LIS_PRECON*);
/* Observation only: no injection, no modification of solver or preconditioner. */
LIS_INT __wrap_lis_precon_create(LIS_SOLVER solver,LIS_PRECON *precon){
 int p=solver->options[LIS_OPTIONS_PRECON];
 if(p==LIS_PRECON_TYPE_MFAMG || p==LIS_PRECON_TYPE_MFAMG_ILUTP){
  calls++; if(solver->A->matrix_type!=LIS_MATRIX_CSR || solver->near_nullspace_dim!=1 || !solver->near_nullspace || fabs(solver->near_nullspace[0]->value[0]-1)>1e-14)errors++;
  printf("OBS call=%d precon=%d matrix=%d dim=%d\n",calls,p,solver->A->matrix_type,solver->near_nullspace_dim);
 }
 return __real_lis_precon_create(solver,precon);
}
static LIS_MATRIX diag(int n,int which){LIS_MATRIX A=NULL;int i;lis_matrix_create(0,&A);lis_matrix_set_size(A,0,n);for(i=0;i<n;i++){double a=which?(i==0?1.0:(i==1?1.0/1.01:1.0/(3.0+i*.15))):1.0;lis_matrix_set_value(LIS_INS_VALUE,i,i,a,A);}lis_matrix_set_type(A,LIS_MATRIX_CSR);lis_matrix_assemble(A);return A;}
int main(int argc,char **argv){
 int n=12,b=atoi(argv[1]),status=-1,rc,i,req,ss; double lam=0; char opts[256]; LIS_MATRIX A,B; LIS_ESOLVER e; LIS_VECTOR x,v,bad;
 const char *pre=argv[2],*method=argv[3],*transform=argc>4?argv[4]:"buckling";
 lis_initialize(&argc,&argv); A=diag(n,0);B=diag(n,1); OK(lis_vector_duplicate(A,&x)); OK(lis_vector_duplicate(A,&v)); OK(lis_esolver_create(&e));
 OK(lis_vector_set_all(1,v)); OK(lis_esolver_set_inner_near_nullspace(e,1,&v));
 /* The caller's vectors may be overwritten and destroyed after the setter. */
 OK(lis_vector_set_all(7,v)); OK(lis_vector_destroy(v));
 /* Failed replacement must preserve the previous valid copy. */
 CHECK(lis_esolver_set_inner_near_nullspace(e,1,NULL)==LIS_ERR_ILL_ARG);
 OK(lis_vector_duplicate(A,&bad)); OK(lis_vector_set_all(NAN,bad)); CHECK(lis_esolver_set_inner_near_nullspace(e,1,&bad)==LIS_ERR_ILL_ARG);OK(lis_vector_destroy(bad));
 req=!strcmp(method,"gii")?1:2;ss=req==1?1:12;if(getenv("TEST_SS"))ss=atoi(getenv("TEST_SS"));
 snprintf(opts,sizeof opts,"-e %s -etransform %s -shift 0.98 -ss %d -etol 1e-8 -emaxiter 500",method,transform,ss);OK(lis_esolver_set_option(opts,e));
 snprintf(opts,sizeof opts,"-nev %d",req);OK(lis_esolver_set_option(opts,e));e->block_size=b;
 snprintf(opts,sizeof opts,"-i %s -p %s -tol 1e-13 -maxiter 500",!strcmp(pre,"mfamg")?"bicgstab":"fgmres",pre);OK(lis_esolver_set_linear_options(opts,e));
 for(i=0;i<2;i++){
  rc=lis_gesolve(A,B,x,&lam,e);OK(rc); OK(lis_esolver_get_status(e,&status));
  printf("SOLVE repeat=%d status=%d lambda=%.17g count=%d residual=%.17g\n",i,status,lam,e->nev,e->resid[0]);
  CHECK(status==LIS_SUCCESS); CHECK(fabs(lam-1)<1e-7); CHECK(e->nev>=req); CHECK(e->resid[0]<1e-8);if(req==2){CHECK(fabs(e->evalue[1]-1.01)<1e-7);CHECK(e->resid[1]<1e-8);}
 }
 CHECK(calls>=2 && errors==0);
 /* A different matrix/vector dimension must fail before inner setup. */
 { LIS_MATRIX C=diag(13,0),D=diag(13,1); LIS_VECTOR y;int before=calls;OK(lis_vector_duplicate(C,&y));CHECK(lis_gesolve(C,D,y,&lam,e)==LIS_ERR_ILL_ARG);CHECK(calls==before);lis_vector_destroy(y);lis_matrix_destroy(C);lis_matrix_destroy(D); }
 OK(lis_esolver_set_option("-e gai",e));CHECK(lis_gesolve(A,B,x,&lam,e)==LIS_ERR_NOT_IMPLEMENTED);
 OK(lis_esolver_set_inner_near_nullspace(e,0,NULL)); CHECK(e->inner_near_nullspace_owner==NULL); OK(lis_esolver_set_inner_near_nullspace(e,0,NULL));
 lis_esolver_destroy(e);lis_vector_destroy(x);lis_matrix_destroy(A);lis_matrix_destroy(B);lis_finalize();puts("PASS");return 0;
}
