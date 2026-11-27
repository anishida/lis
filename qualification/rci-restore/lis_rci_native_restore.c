/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_native_restore.h"
#include "lis_precon.h"
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#if defined(USE_MPI)
#error "Native restoration adapter is non-MPI"
#endif
_Static_assert(_Generic((LIS_SCALAR)0,double:1,default:0),"real-double only");
struct lis_rci_native_restore {
 LIS_MATRIX A,B;
 LIS_SOLVER solver;
 LIS_PRECON precon;
 LIS_VECTOR in,out;
 LIS_INT n;
 LIS_UNSIGNED_INT as,av,bs,bv;
 int bad_a,bad_b;
 double *buf;
};
static int code(LIS_INT e)
{
 return e==0?0:((int64_t)e>=INT_MIN&&(int64_t)e<=INT_MAX?(int)e:1);
}
static int dimensions(LIS_MATRIX a,LIS_INT *n)
{
 LIS_INT g=0;
 return a&&lis_matrix_is_assembled(a)&&!lis_matrix_get_size(a,n,&g)&&*n>0&&*n==g;
}
lis_rci_native_restore *lis_rci_native_restore_create(LIS_MATRIX A,LIS_MATRIX B,LIS_SOLVER solver)
{
 LIS_INT n=0,m=0;
 if(LIS_RCI_CHECK_CLIENT_ABI()!=0||!dimensions(A,&n)||(B&&(!dimensions(B,&m)||m!=n))||
    (uintmax_t)n>SIZE_MAX/sizeof(double)||(solver&&(solver->A!=A||!solver->precon)))
  return NULL;
 lis_rci_native_restore *s=calloc(1,sizeof(*s));if(!s)return NULL;
 s->A=A;s->B=B;s->n=n;s->solver=solver;s->precon=solver?solver->precon:NULL;
 if(lis_matrix_get_generation(A,&s->as,&s->av)||
    (B&&lis_matrix_get_generation(B,&s->bs,&s->bv))||
    lis_vector_duplicate(A,&s->in)||lis_vector_duplicate(A,&s->out)){
  lis_rci_native_restore_destroy(s);return NULL;
 }
 s->buf=malloc((size_t)n*sizeof(double));if(!s->buf){lis_rci_native_restore_destroy(s);return NULL;}
 return s;
}
void lis_rci_native_restore_destroy(lis_rci_native_restore *s)
{if(s){if(s->in)lis_vector_destroy(s->in);if(s->out)lis_vector_destroy(s->out);free(s->buf);free(s);}}
int lis_rci_native_restore_generations(void *user,uint64_t *a,uint64_t *b)
{
 lis_rci_native_restore *s=user;
 if(!s||!a||!b)return LIS_RCI_BADARG;
 LIS_UNSIGNED_INT sg=0,vg=0;LIS_INT e=lis_matrix_get_generation(s->A,&sg,&vg);
 if(e)return code(e);
 if(sg!=s->as||vg!=s->av||(s->solver&&(s->solver->A!=s->A||s->solver->precon!=s->precon)))s->bad_a=1;
 if(s->B){
  e=lis_matrix_get_generation(s->B,&sg,&vg);if(e)return code(e);
  if(sg!=s->bs||vg!=s->bv)s->bad_b=1;
 }
 *a=s->bad_a?2:1;*b=s->bad_b?2:1;return 0;
}
static int unchanged(lis_rci_native_restore *s)
{
 uint64_t a=0,b=0;int e=lis_rci_native_restore_generations(s,&a,&b);
 return e?e:a==1&&b==1?0:LIS_RCI_PROTOCOL;
}
int lis_rci_native_restore_apply(void *user,uint32_t op,size_t n,const double *in,double *out)
{
 lis_rci_native_restore *s=user;
 if(!s||!in||!out||in==out||n!=(size_t)s->n||
    (op!=LIS_RCI_OP_A&&op!=LIS_RCI_OP_B&&op!=LIS_RCI_OP_PRECON))return LIS_RCI_BADARG;
 int g=unchanged(s);if(g)return g;
 for(size_t i=0;i<n;i++)if(!isfinite(in[i]))return LIS_RCI_NONFINITE;
 if((op==LIS_RCI_OP_B&&!s->B)||(op==LIS_RCI_OP_PRECON&&!s->solver)){
  memcpy(out,in,n*sizeof(double));return 0;
 }
 memcpy(s->buf,in,n*sizeof(double));
 LIS_INT e=lis_vector_set_values2(LIS_INS_VALUE,0,s->n,s->buf,s->in);
 if(!e)e=lis_vector_set_all(NAN,s->out);
 if(!e)e=op==LIS_RCI_OP_PRECON?lis_psolve(s->solver,s->in,s->out):
  lis_matvec(op==LIS_RCI_OP_A?s->A:s->B,s->in,s->out);
 if(!e)e=lis_vector_get_values(s->out,0,s->n,s->buf);
 if(e)return code(e);
 for(size_t i=0;i<n;i++)if(!isfinite(s->buf[i]))return LIS_RCI_NONFINITE;
 g=unchanged(s);if(g)return g;
 memcpy(out,s->buf,n*sizeof(double));return 0;
}
int lis_rci_native_restore_solve(lis_rci_native_restore *s,lis_rci_context *c,
 const double *b,double *x,lis_rci_info *info)
{
 if(!info)return LIS_RCI_BADARG;
 *info=(lis_rci_info){.residual=NAN,.status=LIS_RCI_BADARG};
 extern size_t lis_rci_restore_dimension(const lis_rci_context *);
 if(!s||!c||!b||!x||lis_rci_restore_dimension(c)!=(size_t)s->n)return LIS_RCI_BADARG;
 int g=unchanged(s);if(g){info->status=g==LIS_RCI_PROTOCOL?g:LIS_RCI_CALLBACK;info->callback_code=g==LIS_RCI_PROTOCOL?0:g;return info->status;}
 lis_rci_request q;int status=lis_rci_start(c,b,x,&q);
 while(status==LIS_RCI_MATVEC||status==LIS_RCI_PRECON){
  int e=lis_rci_native_restore_apply(s,status==LIS_RCI_MATVEC?LIS_RCI_OP_A:LIS_RCI_OP_PRECON,q.n,q.in,q.out);
  g=unchanged(s);
  if(g==LIS_RCI_PROTOCOL){
   lis_rci_request_full pending;
   if(lis_rci_get_pending(c,&pending)>0)(void)lis_rci_cancel(c,pending.request.token);
   lis_rci_get_info(c,info);info->status=LIS_RCI_PROTOCOL;return LIS_RCI_PROTOCOL;
  }
  status=lis_rci_advance(c,e?e:g,&q);
 }
 lis_rci_get_info(c,info);
 if(status==LIS_RCI_DONE)memcpy(x,lis_rci_solution(c),(size_t)s->n*sizeof(double));
 return status;
}
