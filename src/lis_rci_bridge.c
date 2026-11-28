/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include "lis_rci_bridge.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>
/* These optimizations invalidate the nonfinite-output protocol. */
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
#error "experimental RCI requires IEEE nonfinite checks (disable fast/finite-only math)"
#endif
#if defined(USE_MPI)
#error Bridge is qualified only for non-MPI LIS
#endif
_Static_assert(_Generic((LIS_SCALAR)0,double:1,default:0),"real-double LIS only");
_Static_assert(sizeof(LIS_INT)<=sizeof(int64_t),"unsupported LIS_INT");
int lis_rci_bridge_check_abi(unsigned ib,unsigned sb,unsigned rb,unsigned kind){
 LIS_ABI_INFO a={0};
 if(lis_abi_get_info(&a,sizeof(a))!=LIS_ABI_SUCCESS)return LIS_RCI_BADARG;
 if(a.flags&LIS_ABI_FLAG_MPI || kind!=LIS_ABI_SCALAR_REAL || sb!=sizeof(double) || rb!=sizeof(double) ||
    !lis_abi_is_compatible(&a,ib,sb,rb,kind))return LIS_RCI_BADARG;
 return LIS_RCI_DONE;
}
struct lis_rci_bridge { LIS_MATRIX A; LIS_VECTOR u,v; lis_rci_context *ctx; LIS_INT n; double *buf,*diag; int jacobi,cg; };
static int lis_error(LIS_INT e){int64_t k=(int64_t)e;return (k>=INT_MIN&&k<=INT_MAX&&k!=0)?(int)k:1;}
void lis_rci_bridge_destroy(lis_rci_bridge*p){if(!p)return;if(p->v)lis_vector_destroy(p->v);if(p->u)lis_vector_destroy(p->u);lis_rci_destroy(p->ctx);free(p->buf);free(p->diag);free(p);}
int lis_rci_bridge_create(LIS_MATRIX A,const lis_rci_options *o,int jacobi,lis_rci_bridge **out){
 if(!out)return LIS_RCI_BADARG;
 *out=NULL;
 if(LIS_RCI_CHECK_CLIENT_ABI()!=LIS_RCI_DONE)return LIS_RCI_BADARG;
 if(!A||!o||!o->n||o->n>SIZE_MAX/sizeof(double)||(jacobi!=0&&jacobi!=1))return LIS_RCI_BADARG;
 LIS_INT n=0,g=0;
 if(lis_matrix_get_size(A,&n,&g)!=0||n<=0||n!=g||(uintmax_t)n!=(uintmax_t)o->n)return LIS_RCI_BADARG;
 lis_rci_bridge*p=calloc(1,sizeof(*p));if(!p)return LIS_RCI_BADARG;
 p->A=A;p->n=n;p->jacobi=jacobi;p->cg=(o->method==LIS_RCI_CG);
 p->ctx=lis_rci_create(o);p->buf=malloc(o->n*sizeof(double));if(jacobi)p->diag=malloc(o->n*sizeof(double));
 if(!p->ctx||!p->buf||(jacobi&&!p->diag)){lis_rci_bridge_destroy(p);return LIS_RCI_BADARG;}
 if(lis_vector_duplicate(A,&p->u)!=0||lis_vector_duplicate(A,&p->v)!=0){lis_rci_bridge_destroy(p);return LIS_RCI_CALLBACK;}
 *out=p;return 0;
}
int lis_rci_bridge_solve(lis_rci_bridge*p,const double*b,double*x,lis_rci_info*info){
 if(!info)return LIS_RCI_BADARG;
 *info=(lis_rci_info){.residual=NAN,.status=LIS_RCI_BADARG};
 if(!p||!b||!x)return LIS_RCI_BADARG;
 LIS_INT n=p->n,e=0,local=0,global=0;
 if(lis_matrix_get_size(p->A,&local,&global)!=0||local!=n||global!=n)return LIS_RCI_BADARG;
 if(p->jacobi){
  e=lis_matrix_get_diagonal(p->A,p->v);if(!e)e=lis_vector_get_values(p->v,0,n,p->diag);
  if(e){info->status=LIS_RCI_CALLBACK;info->callback_code=lis_error(e);return LIS_RCI_CALLBACK;}
  for(LIS_INT i=0;i<n;i++)if(!isfinite(p->diag[i])||p->diag[i]==0||(p->cg&&p->diag[i]<0))return LIS_RCI_BADARG;
 }
 lis_rci_request q;int s=lis_rci_start(p->ctx,b,x,&q);
 while(s==LIS_RCI_MATVEC||s==LIS_RCI_PRECON){
  e=0;
  if(s==LIS_RCI_MATVEC){memcpy(p->buf,q.in,q.n*sizeof(double));
   e=lis_vector_set_values2(LIS_INS_VALUE,0,n,p->buf,p->u);
   if(!e)e=lis_matvec(p->A,p->u,p->v);
   if(!e)e=lis_vector_get_values(p->v,0,n,q.out);
  }else for(size_t i=0;i<q.n;i++)q.out[i]=p->jacobi?q.in[i]/p->diag[i]:q.in[i];
  s=lis_rci_advance(p->ctx,e?lis_error(e):0,&q);
 }
 lis_rci_get_info(p->ctx,info);if(s==LIS_RCI_DONE)memcpy(x,lis_rci_solution(p->ctx),(size_t)n*sizeof(double));return s;
}
