#include "lis.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#if defined(_COMPLEX) && defined(HAVE_COMPLEX_H)
#include <complex.h>
#endif

/*
 * Regression for #38: contiguous setters must honor the requested
 * global start index, including LIS_ORIGIN_1 and ADD semantics.
 */
static int test_vector_set_values_indices(void)
{
    LIS_VECTOR v=NULL;
    LIS_SCALAR in[2]={(LIS_SCALAR)11.0,(LIS_SCALAR)22.0};
    LIS_SCALAR add[2]={(LIS_SCALAR)1.0,(LIS_SCALAR)2.0};
    LIS_INT indices[2]={2,4};
    LIS_INT err;
    int failed=1;

    err=lis_vector_create(LIS_COMM_WORLD,&v);
    if(err) goto done;
    err=lis_vector_set_size(v,0,5);
    if(err || v->n!=5) goto done;
    err=lis_vector_set_all((LIS_SCALAR)0.0,v);
    if(err) goto done;

    /* Zero-based insertion must write at start, not at index zero. */
    err=lis_vector_set_values2(LIS_INS_VALUE,2,2,in,v);
    if(err || v->value[0]!=(LIS_SCALAR)0.0 ||
       v->value[1]!=(LIS_SCALAR)0.0 ||
       v->value[2]!=in[0] || v->value[3]!=in[1] ||
       v->value[4]!=(LIS_SCALAR)0.0) goto done;

    err=lis_vector_set_values2(LIS_ADD_VALUE,2,2,add,v);
    if(err || v->value[2]!=(LIS_SCALAR)12.0 ||
       v->value[3]!=(LIS_SCALAR)24.0) goto done;

    /* Reject an overlong segment without partially modifying v. */
    err=lis_vector_set_values2(LIS_INS_VALUE,4,2,in,v);
    if(err!=LIS_ERR_ILL_ARG ||
       v->value[4]!=(LIS_SCALAR)0.0) goto done;
    err=lis_vector_set_values2(LIS_INS_VALUE,1,-1,in,v);
    if(err!=LIS_ERR_ILL_ARG) goto done;

    v->origin=1;
    /* One-based start=2 means zero-based vector position 1. */
    err=lis_vector_set_values2(LIS_INS_VALUE,2,2,in,v);
    if(err || v->value[0]!=(LIS_SCALAR)0.0 ||
       v->value[1]!=in[0] || v->value[2]!=in[1] ||
       v->value[3]!=(LIS_SCALAR)24.0) goto done;

    /* Indexed ADD uses the same origin conversion as indexed INS. */
    err=lis_vector_set_values(LIS_ADD_VALUE,2,indices,add,v);
    if(err || v->value[1]!=(LIS_SCALAR)12.0 ||
       v->value[2]!=(LIS_SCALAR)22.0 ||
       v->value[3]!=(LIS_SCALAR)26.0) goto done;

    /* Out-of-range one-based access must leave all values unchanged. */
    err=lis_vector_set_values2(LIS_INS_VALUE,5,2,in,v);
    if(err!=LIS_ERR_ILL_ARG ||
       v->value[4]!=(LIS_SCALAR)0.0) goto done;
    err=lis_vector_set_values2(LIS_ADD_VALUE,6,0,NULL,v);
    if(err!=LIS_SUCCESS) goto done;

    failed=0;
done:
    if(v) lis_vector_destroy(v);
    return failed;
}

/*
 * Validate overflow rejection and preservation of tracked allocations when
 * lis_realloc() fails.  The old size+16 calculation could wrap on SIZE_MAX.
 */
static int test_lis_allocator_failures(void)
{
    unsigned char *p=NULL, *q=NULL, *zero=NULL;
    int failed=1;
    int i;

    p=(unsigned char *)lis_malloc(SIZE_MAX,"allocator-overflow");
    if(p!=NULL) goto done;
    zero=(unsigned char *)lis_calloc(SIZE_MAX,"calloc-overflow");
    if(zero!=NULL) goto done;

    p=(unsigned char *)lis_malloc(32,"allocator-realloc");
    if(p==NULL || !lis_is_malloc(p)) goto done;
    for(i=0;i<32;i++) p[i]=(unsigned char)(i+1);

    q=(unsigned char *)lis_realloc(p,SIZE_MAX);
    if(q!=NULL)
    {
        p=q;
        goto done;
    }
    if(!lis_is_malloc(p)) goto done;
    for(i=0;i<32;i++) if(p[i]!=(unsigned char)(i+1)) goto done;

    q=(unsigned char *)lis_realloc(p,80);
    if(q==NULL) goto done;
    p=q;
    if(!lis_is_malloc(p)) goto done;
    for(i=0;i<32;i++) if(p[i]!=(unsigned char)(i+1)) goto done;

    q=(unsigned char *)lis_realloc(p,8);
    if(q==NULL) goto done;
    p=q;
    if(!lis_is_malloc(p)) goto done;
    for(i=0;i<8;i++) if(p[i]!=(unsigned char)(i+1)) goto done;

    zero=(unsigned char *)lis_calloc(48,"calloc-zero");
    if(zero==NULL || !lis_is_malloc(zero)) goto done;
    for(i=0;i<48;i++) if(zero[i]!=0) goto done;

    failed=0;
done:
    if(p) lis_free(p);
    if(zero) lis_free(zero);
    return failed;
}

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
    if(test_lis_allocator_failures()) {
        fprintf(stderr,"ALLOCATOR_FAILURE_PATHS=FAIL\n");
        lis_finalize();
        return 1;
    }
    if(test_vector_set_values_indices()) {
        fprintf(stderr,"VECTOR_SET_VALUES_INDICES=FAIL\n");
        lis_finalize();
        return 1;
    }
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
