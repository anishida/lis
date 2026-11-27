#ifndef MODE_CERTIFICATE_H
#define MODE_CERTIFICATE_H
#include "buckle_bridge.h"
/* Independent original-pencil residual, not the residual of the shifted solve.
   Relative residual = ||Kx-lambda Gx||_2 / (||Kx||_2 + |lambda| ||Gx||_2).
   No assertion about completeness, multiplicity or eigenvalue ordering.
   On error output unchanged. P55_OK means computed, not necessarily accepted.
   max_bytes limits requested scratch only. */
typedef struct {long double relative_residual; int accepted;} p70_certificate;
int p70_certify(const p69_csr *k,const p69_csr *g,double lambda,const double *x,size_t nx,long double tol,size_t max_bytes,p70_certificate *out);
#endif
