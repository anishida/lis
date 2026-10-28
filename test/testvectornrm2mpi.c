#include "lis.h"

#include <math.h>
#include <stdio.h>

#ifdef _COMPLEX
#include <complex.h>
#endif


static LIS_SCALAR
test_value(
    LIS_INT index)
{
#ifdef _COMPLEX
    switch(index)
    {
    case 0:
        return (LIS_SCALAR)(3.0 + 4.0*I);

    case 1:
        return (LIS_SCALAR)0.0;

    case 2:
        return (LIS_SCALAR)(5.0 + 12.0*I);

    case 3:
        return (LIS_SCALAR)(8.0 + 15.0*I);

    default:
        return (LIS_SCALAR)0.0;
    }
#else
    switch(index)
    {
    case 0:
        return (LIS_SCALAR)3.0;

    case 1:
        return (LIS_SCALAR)0.0;

    case 2:
        return (LIS_SCALAR)5.0;

    case 3:
        return (LIS_SCALAR)8.0;

    default:
        return (LIS_SCALAR)0.0;
    }
#endif
}


int
main(
    int argc,
    char **argv)
{
    LIS_VECTOR v = NULL;
    LIS_REAL norm;
    LIS_REAL expected;
    LIS_REAL tolerance;
    LIS_INT is,ie,i,err;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if(err)
        return 1;

    err = lis_vector_create(LIS_COMM_WORLD,&v);
    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_set_size(v,0,4);
    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_set_all((LIS_SCALAR)0.0,v);
    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_get_range(v,&is,&ie);
    if(err)
    {
        failed = 1;
        goto cleanup;
    }

    for(i=is;i<ie;i++)
    {
        err = lis_vector_set_value(
            LIS_INS_VALUE,
            i,
            test_value(i),
            v);

        if(err)
        {
            failed = 1;
            goto cleanup;
        }
    }

    err = lis_vector_nrm2(v,&norm);
    if(err)
    {
        fprintf(stderr,
            "lis_vector_nrm2 returned error %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

#ifdef _COMPLEX
    expected = (LIS_REAL)sqrt((LIS_REAL)483.0);
#else
    expected = (LIS_REAL)sqrt((LIS_REAL)98.0);
#endif

    tolerance =
        (LIS_REAL)1.0e-12 *
        ((LIS_REAL)1.0 + expected);

    if(!isfinite((double)norm) ||
       fabs(norm-expected)>tolerance)
    {
        fprintf(stderr,
            "vector nrm2 mismatch: actual=%.17e expected=%.17e\n",
            (double)norm,
            (double)expected);

        failed = 1;
        goto cleanup;
    }

cleanup:

    if(v)
        lis_vector_destroy(v);

    lis_finalize();

    if(failed)
    {
        fprintf(stderr,
            "VECTOR_NRM2_MPI FAILED\n");

        return 1;
    }

    printf("VECTOR_NRM2_MPI PASSED\n");

    return 0;
}
