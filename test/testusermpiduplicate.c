#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <stdio.h>
#include <math.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

typedef struct { LIS_INT n,is; } USER_ROWS;
static LIS_INT get_row(void *data,LIS_INT row,LIS_INT *nnz,
                        LIS_INT *cols,LIS_SCALAR *vals)
{
    USER_ROWS *ctx=(USER_ROWS *)data;
    if(!nnz || row<0 || row>=ctx->n)return LIS_ERR_ILL_ARG;
    if(cols==NULL && vals==NULL){*nnz=1;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<1)return LIS_ERR_ILL_ARG;
    cols[0]=ctx->is+row;
    vals[0]=(LIS_SCALAR)(2.0+0.25*(double)(ctx->is+row));
    *nnz=1;
    return LIS_SUCCESS;
}
/* Exercise a real distributed CSR halo. Keep the copy usable after the
 * source is destroyed to detect shallowly copied communication metadata. */
static int test_csr_with_ghosts(int nranks)
{
    LIS_MATRIX src=NULL,copy=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT i,g,k,err;
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&src);
    if(err)goto done;
    err=lis_matrix_set_size(src,0,8);
    if(err)goto done;
    for(g=src->is;g<src->ie;g++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,g,g,(LIS_SCALAR)2.0,src);
        if(err)goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,g,(g+1)%8,
                                 (LIS_SCALAR)0.5,src);
        if(err)goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,g,(g+7)%8,
                                 (LIS_SCALAR)0.25,src);
        if(err)goto done;
    }
    err=lis_matrix_set_type(src,LIS_MATRIX_CSR);
    if(err)goto done;
    err=lis_matrix_assemble(src);
    if(err)goto done;
    if(nranks>1 && (src->np<=src->n || src->l2g_map==NULL ||
                    src->commtable==NULL))goto done;

    err=lis_matrix_duplicate(src,&copy);
    if(err || copy==NULL)goto done;
    /* A serial LIS build has NULL ranges on both matrices; that is
     * valid.  In MPI, ranges must be independently owned. */
    if((src->ranges==NULL && copy->ranges!=NULL) ||
       (src->ranges!=NULL &&
        (copy->ranges==NULL || copy->ranges==src->ranges)) ||
       copy->n!=src->n || copy->np!=src->np ||
       copy->gn!=src->gn)goto done;
    if(nranks>1)
    {
        if(copy->commtable==NULL || copy->commtable==src->commtable ||
           copy->l2g_map==NULL || copy->l2g_map==src->l2g_map)
            goto done;
        for(k=0;k<src->np-src->n;k++)
            if(copy->l2g_map[k]!=src->l2g_map[k])goto done;
    }
    err=lis_matrix_copy(src,copy);
    if(err)goto done;
    lis_matrix_destroy(src);
    src=NULL;
    err=lis_vector_duplicate(copy,&x);
    if(err)goto done;
    err=lis_vector_duplicate(copy,&y);
    if(err)goto done;
    for(i=0;i<copy->n;i++)
        x->value[i]=(LIS_SCALAR)(1.0+0.1*(double)(copy->is+i));
    err=lis_matvec(copy,x,y);
    if(err)goto done;
    for(i=0;i<copy->n;i++)
    {
        double row=(double)(copy->is+i);
        double next=(double)((copy->is+i+1)%8);
        double prev=(double)((copy->is+i+7)%8);
        double expected=2.0*(1.0+0.1*row)
                       +0.5*(1.0+0.1*next)
                       +0.25*(1.0+0.1*prev);
        if(fabs((double)y->value[i]-expected)>1e-11)goto done;
    }
    failed=0;
done:
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(copy)lis_matrix_destroy(copy);
    if(src)lis_matrix_destroy(src);
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    USER_ROWS state;
    LIS_INT err=LIS_SUCCESS,i;
    int rank=0,nranks=1,failed=0,all_failed=0;
    if(lis_initialize(&argc,&argv))return 1;
#ifdef USE_MPI
    MPI_Comm_rank(LIS_COMM_WORLD,&rank);
    MPI_Comm_size(LIS_COMM_WORLD,&nranks);
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err){failed=__LINE__;goto done;}
    err=lis_matrix_set_size(A,0,8);
    if(err){failed=__LINE__;goto done;}
    state.n=A->n;state.is=A->is;
    err=lis_matrix_set_user(A,&state,NULL,NULL);
    if(err){failed=__LINE__;goto done;}
    err=lis_matrix_set_user_row(A,get_row);
    if(err){failed=__LINE__;goto done;}
    if(A->commtable!=NULL || A->l2g_map!=NULL || A->np!=A->n)
    {failed=__LINE__;goto done;}
    err=lis_matrix_duplicate(A,&B);
    if(err || B==NULL){failed=__LINE__;goto done;}
    if(B->commtable!=NULL || B->l2g_map!=NULL ||
       B->n!=A->n || B->gn!=A->gn ||
       B->is!=A->is || B->ie!=A->ie ||
       (A->ranges!=NULL && B->ranges==A->ranges) ||
       (A->ranges==NULL && B->ranges!=NULL))
    {failed=__LINE__;goto done;}
    err=lis_matrix_set_type(B,LIS_MATRIX_CSR);
    if(err){failed=__LINE__;goto done;}
    err=lis_matrix_convert(A,B);
    if(err){failed=__LINE__;goto done;}
    err=lis_vector_duplicate(B,&x);
    if(err){failed=__LINE__;goto done;}
    err=lis_vector_duplicate(B,&y);
    if(err){failed=__LINE__;goto done;}
    lis_vector_set_all((LIS_SCALAR)1.0,x);
    err=lis_matvec(B,x,y);
    if(err){failed=__LINE__;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)(y->value[i]-(2.0+0.25*(double)(A->is+i))))>1e-12)
        {failed=__LINE__;goto done;}
    if(test_csr_with_ghosts(nranks)){failed=__LINE__;goto done;}
done:
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
#ifdef USE_MPI
    MPI_Allreduce(&failed,&all_failed,1,MPI_INT,MPI_MAX,LIS_COMM_WORLD);
#else
    all_failed=failed;
#endif
    lis_finalize();
    if(rank==0)printf("USER_MPI_DUPLICATE %s ranks=%d stage=%d err=%d\n",
      all_failed?"FAILED":"PASSED",nranks,all_failed,(int)err);
    return all_failed?1:0;
}
