/* Numerical regression for the FGCRO-DR flexible overlap F = Q^H Z.
 * Tests partial 8x4 tiles, 512-element strip tails and column-major
 * offsets against the original scalar lis_vector_dot() implementation.
 * Equality is tolerance-based, not bitwise (OpenMP reductions differ).
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <math.h>
#include <stdio.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#include "lis.h"
#include "lis_solver.h"
#include "lislib.h"

#define NQ 11
#define NZ 7
#define LD 15
#define COLUMNS 11
#define ROW_OFFSET 1
#define COL_OFFSET 2

static LIS_REAL scalar_abs(LIS_SCALAR x)
{
#ifdef _COMPLEX
    return (LIS_REAL)cabs(x);
#else
    return (LIS_REAL)fabs(x);
#endif
}

static int test_case(LIS_INT n, LIS_INT nq, LIS_INT nz)
{
    LIS_VECTOR q[NQ],z[NZ];
    LIS_SCALAR f[LD*COLUMNS],ref;
    LIS_INT i,j,r,c,index,err=0;
    int failed=1;
    for(r=0;r<NQ;r++) q[r]=NULL;
    for(c=0;c<NZ;c++) z[c]=NULL;

    for(r=0;r<nq;r++)
    {
        if((err=lis_vector_create(LIS_COMM_WORLD,&q[r]))) goto done;
        if((err=lis_vector_set_size(q[r],0,n))) goto done;
        for(i=0;i<q[r]->n;i++)
        {
            LIS_REAL t=(LIS_REAL)(q[r]->is+i+1);
            q[r]->value[i]=(LIS_SCALAR)(sin(0.011*t*(r+1))+
                                        0.003*(r+1));
#ifdef _COMPLEX
            q[r]->value[i]+=(LIS_SCALAR)(I*cos(0.007*t*(r+2)));
#endif
        }
    }
    for(c=0;c<nz;c++)
    {
        if((err=lis_vector_create(LIS_COMM_WORLD,&z[c]))) goto done;
        if((err=lis_vector_set_size(z[c],0,n))) goto done;
        for(i=0;i<z[c]->n;i++)
        {
            LIS_REAL t=(LIS_REAL)(z[c]->is+i+1);
            z[c]->value[i]=(LIS_SCALAR)(cos(0.013*t*(c+1))-
                                        0.002*(c+2));
#ifdef _COMPLEX
            z[c]->value[i]+=(LIS_SCALAR)(I*sin(0.009*t*(c+3)));
#endif
        }
    }
    for(j=0;j<LD*COLUMNS;j++) f[j]=(LIS_SCALAR)-9.5;
    err=lis_fgcrodr_flexible_overlap(
        q,nq,z,nz,f,LD,ROW_OFFSET,COL_OFFSET);
    if(err) goto done;
    for(c=0;c<nz;c++)
    {
        for(r=0;r<nq;r++)
        {
            err=lis_vector_dot(q[r],z[c],&ref);
            if(err) goto done;
            index=ROW_OFFSET+r+(COL_OFFSET+c)*LD;
            if(scalar_abs(f[index]-ref)>
               (LIS_REAL)1e-10*(1.0+scalar_abs(ref)))
            {
                fprintf(stderr,
                    "FGCRODR_BLOCKED_MISMATCH n=%lld r=%lld c=%lld error=%g ref=%g\n",
                    (long long)n,(long long)r,(long long)c,
                    (double)scalar_abs(f[index]-ref),
                    (double)scalar_abs(ref));
                goto done;
            }
        }
    }
    for(j=0;j<LD*COLUMNS;j++)
    {
        c=j/LD;
        r=j%LD;
        if(c>=COL_OFFSET && c<COL_OFFSET+nz &&
           r>=ROW_OFFSET && r<ROW_OFFSET+nq)
            continue;
        if(f[j]!=(LIS_SCALAR)-9.5)
        {
            fprintf(stderr,"FGCRODR_BLOCKED_OVERWRITE index=%lld\n",
                    (long long)j);
            goto done;
        }
    }
    if(lis_fgcrodr_flexible_overlap(
        q,nq,z,nz,f,ROW_OFFSET+nq-1,ROW_OFFSET,COL_OFFSET)
       !=LIS_ERR_ILL_ARG) goto done;
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_BLOCKED_FAIL n=%lld rc=%lld\n",
                (long long)n,(long long)err);
    for(c=0;c<NZ;c++) if(z[c]) lis_vector_destroy(z[c]);
    for(r=0;r<NQ;r++) if(q[r]) lis_vector_destroy(q[r]);
    return failed;
}

int main(int argc,char **argv)
{
    int failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    if(test_case(1103,11,7)) goto done; /* partial row/col tiles, 3 strips */
    if(test_case(1000,9,5)) goto done; /* partial 512-element strip */
    failed=0;
    puts("FGCRODR_BLOCKED_OVERLAP_PASS");
done:
    lis_finalize();
    return failed;
}
