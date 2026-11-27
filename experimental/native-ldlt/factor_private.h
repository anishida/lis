#ifndef FACTOR_PRIVATE_H
#define FACTOR_PRIVATE_H
#include "lis_adapter.h"
typedef struct {
    uint32_t first, width, nrows;
    int short_index;
    size_t nvalues;
    void *rows;
    double *values;
} panel;
typedef struct { uint32_t k; long double a,b,c,det,scale; double original_b; } block;
struct p55_factor {
    uint32_t n; size_t np, nb, index_bytes;
    panel *p; double *diag; block *blocks;
};
struct p56_context {
    p55_factor *factor;
    uint32_t n;
    uint32_t *perm; /* NULL for identity */
};
#endif
