/* AI-assisted reconstruction qualification; subject to LIS COPYING confirmation. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_native_linear_restore.h"
#include "lislib.h"
#include <stdio.h>
#include <float.h>
#define N 9
#define H 201
#define C(e) do{if(!(e)){fprintf(stderr,"PRECISION FAIL %d: %s\n",__LINE__,#e);rc=1;goto cleanup;}}while(0)
typedef struct { LIS_INT status,iter; LIS_REAL residual,history[H]; LIS_SCALAR x[N]; } answer;
static LIS_REAL mag(LIS_SCALAR x){return (LIS_REAL)fabs(x);}
static LIS_REAL tiny(void)
{
 LIS_REAL v=1;
#if defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
 for(int i=0;i<90;i++)v/=2;
#elif defined(_LONG__DOUBLE)
 for(int i=0;i<58;i++)v/=2;
#else
 for(int i=0;i<40;i++)v/=2;
#endif
 return v;
}
static const char *tol_text(void)
{
#if defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
 return "1e-25";
#elif defined(_LONG__DOUBLE)
 return "1e-15";
#else
 return "1e-11";
#endif
}
static LIS_REAL solution_tol(void)
{
#if defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
 return (LIS_REAL)1e-23;
#elif defined(_LONG__DOUBLE)
 return (LIS_REAL)1e-13;
#else
 return (LIS_REAL)1e-9;
#endif
}
static LIS_SCALAR exact_at(int i,int variant)
{
 LIS_SCALAR v=(LIS_REAL)1+(LIS_REAL)i/16+tiny()*(i+1);
#ifdef _COMPLEX
 v+=I*((LIS_REAL)(i+1)/32+(variant?tiny():0));
#else
 (void)variant;
#endif
 return v;
}
static LIS_SCALAR upper(int variant)
{
 LIS_SCALAR v=-(LIS_REAL).25;
#ifdef _COMPLEX
 if(variant)v+=I*(LIS_REAL).125;
#else
 (void)variant;
#endif
 return v;
}
static LIS_SCALAR rhs(int i,int variant)
{
 LIS_SCALAR u=upper(variant);
 LIS_SCALAR v=((LIS_REAL)3+(LIS_REAL)i/8)*exact_at(i,variant);
 if(i)v+=conj(u)*exact_at(i-1,variant);
 if(i+1<N)v+=u*exact_at(i+1,variant);
 return v;
}
static LIS_INT external(LIS_SOLVER s,LIS_INT method)
{
 lis_rci_native_linear_restore *k=lis_rci_native_linear_restore_create();if(!k)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_linear_request q;int st=lis_rci_native_linear_restore_start_method(k,s,method,1,&q);
 uint64_t last=0;
 while(st>0){
  if(q.token<=last){lis_rci_native_linear_restore_destroy(k);return 99;}
  last=q.token;
  LIS_INT err=q.operation==LIS_RCI_OP_A?lis_matvec(q.A,q.in,q.out):lis_psolve(q.solver,q.in,q.out);
  st=lis_rci_native_linear_restore_advance(k,q.token,err,&q);
 }
 LIS_INT native=0;uint64_t a=0,p=0;lis_rci_native_linear_restore_status(k,&native,&a,&p);
 if(!a||!p)native=98;
 lis_rci_native_linear_restore_destroy(k);return native;
}
static int solve(const char *method,int mode,int variant,answer *out)
{
 LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL;LIS_SOLVER s=NULL;
 int rc=0;memset(out,0,sizeof(*out));
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 for(LIS_INT i=0;i<N;i++){
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_REAL)3+(LIS_REAL)i/8,A));
  if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,conj(upper(variant)),A));
  if(i+1<N)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,upper(variant),A));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));C(!lis_solver_create(&s));
 for(LIS_INT i=0;i<N;i++){C(!lis_vector_set_value(LIS_INS_VALUE,i,rhs(i,variant),b));C(!lis_vector_set_value(LIS_INS_VALUE,i,(LIS_REAL).125,x));}
 char opts[160];snprintf(opts,sizeof(opts),"-i %s -p jacobi -tol %s -maxiter 200 -restart 9 -initx_zeros false -print mem",method,tol_text());
 C(!lis_solver_set_option(opts,s));lis_rci_native_linear_test_select(mode==0,mode==2?external:NULL);
 (void)lis_solve(A,b,x,s);
 C(!lis_solver_get_status(s,&out->status));C(!lis_solver_get_iter(s,&out->iter));C(!lis_solver_get_residualnorm(s,&out->residual));C(!lis_vector_get_values(x,0,N,out->x));
 LIS_INT count=out->iter;if(count>200)count=200;
 if(s->rhistory&&count>0)for(LIS_INT i=0;i<count;i++)out->history[i]=s->rhistory[i];
 if(out->status==LIS_SUCCESS){
  LIS_REAL error=0;
  for(int i=0;i<N;i++){
   LIS_REAL d=mag(out->x[i]-exact_at(i,variant));if(d>error)error=d;
  }
  C(error<=solution_tol());
 }
cleanup:
 lis_rci_native_linear_test_select(0,NULL);
 if(s)lis_solver_destroy(s);if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(A)lis_matrix_destroy(A);return rc;
}
static int compare(const char *method,int variant,int *converged)
{
 int rc=0;answer a,b,c;
 C(!solve(method,0,variant,&a));C(!solve(method,1,variant,&b));C(!solve(method,2,variant,&c));
 printf("PRECISION method=%s variant=%d status=%lld/%lld/%lld iterations=%lld/%lld/%lld residual=%Lg\n",method,variant,(long long)a.status,(long long)b.status,(long long)c.status,(long long)a.iter,(long long)b.iter,(long long)c.iter,(long double)a.residual);
 C(a.status==b.status&&a.status==c.status&&a.iter==b.iter&&a.iter==c.iter);
 C((a.residual==b.residual&&a.residual==c.residual)||(isnan(a.residual)&&isnan(b.residual)&&isnan(c.residual)));
 for(int i=0;i<N;i++)C(a.x[i]==b.x[i]&&a.x[i]==c.x[i]);
 for(int i=0;i<H;i++)C(a.history[i]==b.history[i]&&a.history[i]==c.history[i]);
 if(a.status==LIS_SUCCESS)(*converged)++;
cleanup:return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
#if defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
 _Static_assert(_Generic((LIS_REAL)0,__float128:1,default:0),"binary128 required");
 puts("PRECISION binary128");
#elif defined(_LONG__DOUBLE)
 _Static_assert(_Generic((LIS_REAL)0,long double:1,default:0),"long double required");
 if(LDBL_MANT_DIG<=DBL_MANT_DIG)return 1;
 printf("PRECISION extended mantissa=%d\n",LDBL_MANT_DIG);
#else
 puts("PRECISION double");
#endif
#ifdef _COMPLEX
 puts("PRECISION complex");
#endif
#if defined(_LONG__DOUBLE)
 if((LIS_REAL)(1+tiny())==(LIS_REAL)1)return 1;
#endif
 const char *methods[]={"cg","bicgstab","gmres","fgmres"};int rc=0,n=0,converged=0;
 for(int m=0;m<4;m++){
  rc|=compare(methods[m],0,&converged);n++;
#ifdef _COMPLEX
  rc|=compare(methods[m],1,&converged);n++;
#endif
 }
 lis_finalize();
 printf("PRECISION DIFFERENTIAL cases=%d converged=%d\n",n,converged);
 if(!rc&&converged==n)puts("PRECISION FULL PASS");
 return rc||converged!=n;
}
