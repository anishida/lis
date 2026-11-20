/*
 * Diagnostic-only, allocation-overhead A/B for non-split CSR transpose
 * with one OpenMP worker. Run unchanged against both upstream base and
 * the optimized branch on the same host. Timings never gate CI.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#else
#error "CSR one-worker benchmark must be compiled with OpenMP"
#endif

static double monotonic_seconds(void)
{
    struct timespec ts;
    if(clock_gettime(CLOCK_MONOTONIC,&ts)!=0)return -1.;
    return (double)ts.tv_sec+1.e-9*(double)ts.tv_nsec;
}

static int benchmark(int n,int repeats,int batches,const char *variant)
{
    struct LIS_MATRIX_STRUCT a;
    LIS_INT *ptr=NULL,*idx=NULL;
    LIS_SCALAR *val=NULL,*x=NULL,*y=NULL,*expected=NULL;
    LIS_INT r,p=0,c,k;
    int b,iter,failed=1;
    double start,elapsed,checksum;
    size_t maximum=(size_t)n*7u;
    ptr=(LIS_INT *)malloc((size_t)(n+1)*sizeof(*ptr));
    idx=(LIS_INT *)malloc(maximum*sizeof(*idx));
    val=(LIS_SCALAR *)malloc(maximum*sizeof(*val));
    x=(LIS_SCALAR *)malloc((size_t)n*sizeof(*x));
    y=(LIS_SCALAR *)malloc((size_t)n*sizeof(*y));
    expected=(LIS_SCALAR *)calloc((size_t)n,sizeof(*expected));
    if(!ptr||!idx||!val||!x||!y||!expected)goto done;
    memset(&a,0,sizeof(a));
    a.n=n; a.np=n; a.nr=n; a.is_splited=LIS_FALSE;
    a.ptr=ptr; a.index=idx; a.value=val;
    for(r=0;r<n;r++)
    {
        ptr[r]=p;
        for(k=-3;k<=3;k++)
        {
            c=r+k*37;
            if(c<0||c>=n)continue;
            idx[p]=c;
            val[p]=(LIS_SCALAR)(k==0 ? 5.8 : (k%2 ? -0.17 : 0.11));
            p++;
        }
        x[r]=(LIS_SCALAR)(1.0+0.05*sin((double)r*0.031));
    }
    ptr[n]=p;
    for(r=0;r<n;r++)
        for(k=ptr[r];k<ptr[r+1];k++)
            expected[idx[k]]+=val[k]*x[r];
    if(omp_get_max_threads()!=1)
    {
        fprintf(stderr,"FAIL: omp_get_max_threads()=%d, expected 1\n",
                omp_get_max_threads());
        goto done;
    }
    for(iter=0;iter<128;iter++)lis_matvech_csr(&a,x,y);
    for(r=0;r<n;r++)
        if(fabs((double)(y[r]-expected[r]))>1.e-11)
        {
            fprintf(stderr,"FAIL: transpose mismatch i=%lld got=%.14g expected=%.14g\n",
                (long long)r,(double)y[r],(double)expected[r]);
            goto done;
        }
    checksum=0.;
    for(r=0;r<n;r++)checksum+=(double)y[r]*(r%17+1);
    for(b=0;b<batches;b++)
    {
        start=monotonic_seconds();
        if(start<0.)goto done;
        for(iter=0;iter<repeats;iter++)lis_matvech_csr(&a,x,y);
        elapsed=monotonic_seconds()-start;
        if(elapsed<=0. || !isfinite(elapsed))goto done;
        printf("%s,%d,%d,%d,%.12f,%.9f,%.12g\n",
               variant,n,repeats,b,elapsed,
               1.e9*elapsed/(double)repeats,checksum);
        fflush(stdout);
    }
    failed=0;
done:
    free(expected);free(y);free(x);free(val);free(idx);free(ptr);
    return failed;
}

int main(int argc,char **argv)
{
    int n,reps,batches,rc=1;
    if(argc!=5)
    {
        fprintf(stderr,"usage: %s VARIANT N REPEATS BATCHES\n",argv[0]);
        return 2;
    }
    n=atoi(argv[2]);reps=atoi(argv[3]);batches=atoi(argv[4]);
    if(n<128||n>8192||reps<10||reps>1000000||batches<2||batches>32)
        return 2;
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 3;
    rc=benchmark(n,reps,batches,argv[1]);
    lis_finalize();
    return rc;
}
