#include "lis.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#if defined(_COMPLEX) && defined(HAVE_COMPLEX_H)
#include <complex.h>
#endif

int main(int argc,char **argv)
{
    LIS_VECTOR v=NULL;
    LIS_INT err;
    const double src[9]={1,99,99,-2,99,99,3.5,99,99};
    double dst[9]={0};
    double bad[3]={1,NAN,3};
    double invalid[3]={7,8,9};
    int fail=1;
    err=lis_initialize(&argc,&argv);
    if(err) return 1;
    err=lis_vector_create(LIS_COMM_WORLD,&v);
    if(err)goto done;
    err=lis_vector_set_size(v,0,3);
    if(err)goto done;
    if(v->n!=3)goto done; /* single-process regression */
    err=lis_vector_import_real64(v,src,3,3);
    if(err)goto done;
    err=lis_vector_export_real64(
        v,dst,3,3,LIS_VECTOR_REAL64_REJECT_LOSS);
    if(err || dst[0]!=1 || dst[3]!=-2 || dst[6]!=3.5)
        goto done;

    err=lis_vector_import_real64(v,bad,3,1);
    if(err==LIS_SUCCESS)goto done;
    err=lis_vector_export_real64(
        v,invalid,3,1,LIS_VECTOR_REAL64_REJECT_LOSS);
    if(err || invalid[0]!=1 || invalid[1]!=-2 || invalid[2]!=3.5)
        goto done;
    /* Zero-length transfers may use NULL buffers. */
    {
        struct LIS_VECTOR_STRUCT empty;
        memset(&empty,0,sizeof(empty));
        empty.precision=LIS_PRECISION_DEFAULT;
        if(lis_vector_import_real64(&empty,NULL,0,1)!=LIS_SUCCESS ||
           lis_vector_export_real64(&empty,NULL,0,1,
               LIS_VECTOR_REAL64_REJECT_LOSS)!=LIS_SUCCESS)
            goto done;
    }
    if(lis_vector_import_real64(v,src,2,1)==LIS_SUCCESS)
        goto done;
    if(lis_vector_export_real64(
        v,dst,3,0,LIS_VECTOR_REAL64_REJECT_LOSS)==LIS_SUCCESS)
        goto done;
#ifdef _COMPLEX
#if defined(HAVE_COMPLEX_H)
    v->value[1] += (LIS_SCALAR)1.0*I;
#else
    v->value[1][1]=(LIS_REAL)1.0;
#endif
    invalid[0]=7;invalid[1]=8;invalid[2]=9;
    if(lis_vector_export_real64(
        v,invalid,3,1,LIS_VECTOR_REAL64_REJECT_LOSS)==LIS_SUCCESS ||
       invalid[0]!=7 || invalid[1]!=8 || invalid[2]!=9)
        goto done;
    err=lis_vector_export_real64(
        v,invalid,3,1,LIS_VECTOR_REAL64_ALLOW_LOSS);
    if(err || invalid[1]!=-2)goto done;
#endif
    fail=0;
done:
    if(v)lis_vector_destroy(v);
    lis_finalize();
    printf("REAL64_VECTOR_TRANSFER=%s\n",fail?"FAIL":"PASS");
    return fail;
}
