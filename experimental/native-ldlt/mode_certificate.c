#include "mode_certificate.h"
#include <math.h>
#include <stdlib.h>
static void product(const p69_csr *a,const double *x,long double scale,long double *y){
 for(uint32_t i=0;i<a->n;i++)for(size_t z=a->ptr[i]-a->base;z<a->ptr[i+1]-a->base;z++){
  uint32_t j=a->index[z]-a->base;long double v=a->value[z];
  y[i]+=v*((long double)x[j]/scale);
  if(a->storage!=P69_FULL&&i!=j)y[j]+=v*((long double)x[i]/scale);
 }
}
int p70_certify(const p69_csr *k,const p69_csr *g,double lambda,const double *x,size_t nx,long double tol,size_t max_bytes,p70_certificate *out){
 if(!x||!out||!isfinite(lambda)||!isfinite(tol)||tol<=0||tol>=1||p69_validate_csr(k)||p69_validate_csr(g)||k->n!=g->n||nx!=k->n)return P55_INVALID;
 if(nx>max_bytes/(2*sizeof(long double)))return P60_LIMIT;
 long double scale=0;for(size_t i=0;i<nx;i++){if(!isfinite(x[i]))return P55_INVALID;scale=fmaxl(scale,fabsl(x[i]));}
 if(scale==0)return P55_INVALID;
 long double *y=calloc(2*nx,sizeof *y);if(!y)return P55_NOMEM;
 long double *z=y+nx;product(k,x,scale,y);product(g,x,scale,z);
 long double nk=0,ng=0,nr=0;int rc=P55_OK;
 for(size_t i=0;i<nx;i++){long double b=(long double)lambda*z[i],r=y[i]-b;if(!isfinite(y[i])||!isfinite(z[i])||!isfinite(b)||!isfinite(r)){rc=P55_RANGE;goto done;}nk=hypotl(nk,y[i]);ng=hypotl(ng,b);nr=hypotl(nr,r);}
 long double den=nk+ng;if(!isfinite(den)||!isfinite(nr)){rc=P55_RANGE;goto done;}if(den==0){rc=P55_SINGULAR;goto done;}
 long double rel=nr/den;*out=(p70_certificate){rel,rel<=tol};
 done:free(y);return rc;
}
