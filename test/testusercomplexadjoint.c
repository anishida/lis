#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

#if defined(_COMPLEX) && defined(HAVE_COMPLEX_H) && !defined(USE_MPI)
#include <complex.h>
#define USER_ADJOINT_N 3
typedef struct { LIS_SCALAR a[USER_ADJOINT_N][USER_ADJOINT_N]; LIS_INT bad_row; } CTX;
static LIS_INT row(void *raw,LIS_INT i,LIS_INT *nnz,LIS_INT *cols,LIS_SCALAR *vals)
{
    CTX *c=(CTX *)raw;
    LIS_INT k;
    if(!nnz || i<0 || i>=USER_ADJOINT_N)return LIS_ERR_ILL_ARG;
    if(i==c->bad_row)return LIS_ERR_ILL_ARG;
    /* Reverse the column order and split the diagonal into two entries. */
    if(!cols && !vals){*nnz=USER_ADJOINT_N+1;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<USER_ADJOINT_N+1)return LIS_ERR_ILL_ARG;
    for(k=0;k<USER_ADJOINT_N;k++){
        cols[k]=USER_ADJOINT_N-1-k;
        vals[k]=c->a[i][cols[k]]*(cols[k]==i?0.5:1.0);
    }
    cols[USER_ADJOINT_N]=i; vals[USER_ADJOINT_N]=0.5*c->a[i][i];
    *nnz=USER_ADJOINT_N+1;
    return LIS_SUCCESS;
}
static LIS_INT ax(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)raw;
    LIS_INT i,j;
    for(i=0;i<USER_ADJOINT_N;i++){
        y[i]=0;
        for(j=0;j<USER_ADJOINT_N;j++)y[i]+=c->a[i][j]*x[j];
    }
    return LIS_SUCCESS;
}
static LIS_INT ahx(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)raw;
    LIS_INT i,j;
    for(i=0;i<USER_ADJOINT_N;i++){
        y[i]=0;
        for(j=0;j<USER_ADJOINT_N;j++)y[i]+=conj(c->a[j][i])*x[j];
    }
    return LIS_SUCCESS;
}
static int close_scalar(LIS_SCALAR a,LIS_SCALAR b)
{
    return cabs(a-b)<=1e-11*(1+cabs(b));
}
static int check(LIS_MATRIX A,LIS_VECTOR x,LIS_VECTOR z,LIS_VECTOR y,LIS_VECTOR h)
{
    LIS_SCALAR ex[USER_ADJOINT_N]={0},eh[USER_ADJOINT_N]={0},lhs=0,rhs=0;
    LIS_INT i,j;
    /* Reference uses literal asymmetric coefficients, not USER callbacks. */
    const LIS_SCALAR m[USER_ADJOINT_N][USER_ADJOINT_N]={
        {4+I, 1+2*I, -I},
        {2-I, 5-2*I, 3+I},
        {1+I, -2*I, 6+3*I}};
    if(lis_matvec(A,x,y) || lis_matvech(A,z,h))return 1;
    for(i=0;i<USER_ADJOINT_N;i++){
        for(j=0;j<USER_ADJOINT_N;j++){
            ex[i]+=m[i][j]*x->value[j];
            eh[i]+=conj(m[j][i])*z->value[j];
        }
        if(!close_scalar(y->value[i],ex[i]) ||
           !close_scalar(h->value[i],eh[i]))return 1;
        lhs+=conj(z->value[i])*y->value[i];
        rhs+=conj(h->value[i])*x->value[i];
    }
    return !close_scalar(lhs,rhs);
}
static int exercise(void)
{
    CTX c={{{4+I,1+2*I,-I},{2-I,5-2*I,3+I},{1+I,-2*I,6+3*I}},-1};
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,z=NULL,y=NULL,h=NULL;
    LIS_INT i,threadsafe;
    int failed=1;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,USER_ADJOINT_N) ||
       lis_matrix_set_user(A,&c,ax,ahx) ||
       lis_vector_duplicate(A,&x) || lis_vector_duplicate(A,&z) ||
       lis_vector_duplicate(A,&y) || lis_vector_duplicate(A,&h))goto done;
    for(i=0;i<USER_ADJOINT_N;i++){
        x->value[i]=(1+i)+(2-i)*I;
        z->value[i]=(3-i)-(1+i)*I;
    }
    if(check(A,x,z,y,h) || lis_matrix_set_user_row(A,row))goto done;
    c.bad_row=1;
    /* Explicit forward and adjoint must win even if row access fails. */
    if(check(A,x,z,y,h))goto done;
    c.bad_row=-1;
    lis_matrix_destroy(A); A=NULL;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) ||
       lis_matrix_set_size(A,0,USER_ADJOINT_N) ||
       lis_matrix_set_user(A,&c,ax,NULL) ||
       lis_matrix_set_user_row(A,row) || check(A,x,z,y,h))goto done;
    c.bad_row=1;
    if(lis_matvec(A,x,y) || lis_matvech(A,z,h)!=LIS_ERR_ILL_ARG)goto done;
    c.bad_row=-1;
    /* USER descriptors are initialized once; create a separate row-only matrix. */
    lis_matrix_destroy(A); A=NULL;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,USER_ADJOINT_N) ||
       lis_matrix_set_user(A,&c,NULL,NULL) || lis_matrix_set_user_row(A,row))goto done;
    for(threadsafe=0;threadsafe<=1;threadsafe++){
        if(lis_matrix_set_user_row_thread_safe(A,threadsafe) ||
           check(A,x,z,y,h))goto done;
        if(lis_matrix_get_diagonal(A,h))goto done;
        for(i=0;i<USER_ADJOINT_N;i++)
            if(!close_scalar(h->value[i],c.a[i][i]))goto done;
        c.bad_row=1;
        if(lis_matvec(A,x,y)!=LIS_ERR_ILL_ARG ||
           lis_matvech(A,z,h)!=LIS_ERR_ILL_ARG)goto done;
        c.bad_row=-1;
        if(check(A,x,z,y,h))goto done;
    }
    failed=0;
done:
    if(h)lis_vector_destroy(h);
    if(y)lis_vector_destroy(y);
    if(z)lis_vector_destroy(z);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* A larger sparse operator exercises repeated per-thread row-buffer reuse.
 * The callback reads immutable data only; the expected result is evaluated
 * directly from the stencil, independently of row packing. */
#define USER_STRESS_N 257
static LIS_INT stress_row(void *raw,LIS_INT i,LIS_INT *nnz,LIS_INT *cols,LIS_SCALAR *vals)
{
    (void)raw;
    if(!nnz || i<0 || i>=USER_STRESS_N)return LIS_ERR_ILL_ARG;
    if(!cols && !vals){*nnz=4;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<4)return LIS_ERR_ILL_ARG;
    cols[0]=(i+1)%USER_STRESS_N; vals[0]=1+2*I;
    cols[1]=i; vals[1]=2+I;
    cols[2]=(i+USER_STRESS_N-1)%USER_STRESS_N; vals[2]=-1+I;
    cols[3]=i; vals[3]=2+I;
    *nnz=4;
    return LIS_SUCCESS;
}
static int stress(void)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL,h=NULL;
    LIS_INT i,repeat,threadsafe;
    int failed=1;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) ||
       lis_matrix_set_size(A,0,USER_STRESS_N) ||
       lis_matrix_set_user(A,NULL,NULL,NULL) ||
       lis_matrix_set_user_row(A,stress_row) ||
       lis_vector_duplicate(A,&x) || lis_vector_duplicate(A,&y) ||
       lis_vector_duplicate(A,&h))goto done;
    for(i=0;i<USER_STRESS_N;i++)
        x->value[i]=(1.0+i%13)+(double)(i%7-3)*I;
    for(threadsafe=0;threadsafe<=1;threadsafe++){
        if(lis_matrix_set_user_row_thread_safe(A,threadsafe))goto done;
        for(repeat=0;repeat<16;repeat++){
            if(lis_matvec(A,x,y) || lis_matvech(A,x,h))goto done;
            for(i=0;i<USER_STRESS_N;i++){
                LIS_INT prev=(i+USER_STRESS_N-1)%USER_STRESS_N;
                LIS_INT next=(i+1)%USER_STRESS_N;
                LIS_SCALAR ex=(4+2*I)*x->value[i]+(1+2*I)*x->value[next]+(-1+I)*x->value[prev];
                LIS_SCALAR eh=(4-2*I)*x->value[i]+(1-2*I)*x->value[prev]+(-1-I)*x->value[next];
                if(!close_scalar(y->value[i],ex) || !close_scalar(h->value[i],eh)){
                    fprintf(stderr,"USER stress mismatch row=%d repeat=%d opt-in=%d\n",
                            (int)i,(int)repeat,(int)threadsafe);
                    goto done;
                }
            }
        }
    }
    failed=0;
done:
    if(h)lis_vector_destroy(h);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}
#endif
int main(int argc,char **argv)
{
#if defined(_COMPLEX) && defined(HAVE_COMPLEX_H) && !defined(USE_MPI)
    int failed;
    if(lis_initialize(&argc,&argv))return 1;
    failed=exercise();
    failed|=stress();
    lis_finalize();
    printf("USER_COMPLEX_ADJOINT %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#else
    (void)argc;(void)argv;
    puts("USER_COMPLEX_ADJOINT SKIP: requires serial complex build");
    return 77;
#endif
}
