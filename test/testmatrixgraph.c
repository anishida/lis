#include "lis.h"
#include "lis_matrix_graph.h"
#include <stdio.h>
#ifdef _COMPLEX
#include <complex.h>
#endif

typedef struct { int unused; } USER_CTX;

static LIS_INT user_row(void *user_data, LIS_INT row, LIS_INT *nnz,
                        LIS_INT *cols, LIS_SCALAR *vals)
{
    LIS_INT need;
    (void)user_data;

    if(nnz==NULL || row<0 || row>=4) return LIS_ERR_ILL_ARG;

    /*
     * row 0: four repeated entries cancel algebraically -> isolated row 0.
     * row 1: three repeated entries for column 2 remain nonzero.
     * row 2: reverse edge to row 1.
     * row 3: isolated.
     */
    need = row==0 ? 4 : (row==1 ? 3 : (row==2 ? 1 : 0));

    if(cols==NULL && vals==NULL)
    {
        *nnz = need;
        return LIS_SUCCESS;
    }
    if((cols==NULL)!=(vals==NULL) || *nnz<need) return LIS_ERR_ILL_ARG;

    if(row==0)
    {
        cols[0] = 1; vals[0] = (LIS_SCALAR) 2.0;
        cols[1] = 1; vals[1] = (LIS_SCALAR)-2.0;
        cols[2] = 1; vals[2] = (LIS_SCALAR) 1.0;
        cols[3] = 1; vals[3] = (LIS_SCALAR)-1.0;
    }
    else if(row==1)
    {
        cols[0] = 2; vals[0] = (LIS_SCALAR)1.0;
        cols[1] = 2; vals[1] = (LIS_SCALAR)1.0;
        cols[2] = 2; vals[2] = (LIS_SCALAR)1.0;
    }
    else if(row==2)
    {
        cols[0] = 1; vals[0] = (LIS_SCALAR)3.0;
    }

    *nnz = need;
    return LIS_SUCCESS;
}

static int labels_same_partition(const LIS_INT *a, const LIS_INT *b, LIS_INT n)
{
    LIS_INT i,j;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if((a[i]==a[j]) != (b[i]==b[j])) return 0;
    return 1;
}

static int check_distributed_rejection(void)
{
    LIS_MATRIX A = NULL;
    LIS_INT labels[4],count,first,last,i,err;
    int failed = 0;

    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) return 1;
    err = lis_matrix_set_size(A,0,4);
    if(err) { lis_matrix_destroy(A); return 1; }

    if(A->nprocs==1)
    {
        lis_matrix_destroy(A);
        return 0;
    }

    lis_matrix_get_range(A,&first,&last);
    for(i=first;i<last;i++)
        if(lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,A))
            failed = 1;

    if(!failed) err = lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) failed = 1;
    if(!failed) err = lis_matrix_assemble(A);
    if(err) failed = 1;

    if(!failed)
    {
        err = lis_matrix_connected_components(A,(LIS_REAL)0.0,labels,&count);
        if(err!=LIS_ERR_NOT_IMPLEMENTED)
        {
            fprintf(stderr,"distributed matrix was not rejected: err=%d\n",(int)err);
            failed = 1;
        }
    }

    lis_matrix_destroy(A);
    return failed;
}

int main(int argc, char **argv)
{
    LIS_MATRIX A=NULL,C=NULL,U=NULL,O=NULL;
#ifdef _COMPLEX
    LIS_MATRIX Z=NULL;
#endif
    LIS_INT err,count,i;
    LIS_INT labels[4],labels_csc[4],labels_user[4];
    LIS_INT expected_csr[4] = {0,0,1,1};
    LIS_INT expected_user[4] = {0,1,1,2};
    USER_CTX ctx;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if(err) return 1;

    failed |= check_distributed_rejection();

    /* The remainder is a serial regression. */
    err = lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) failed = 1;
    if(!failed && A->nprocs!=1)
    {
        lis_matrix_destroy(A);
        lis_finalize();
        if(failed) return 1;
        printf("MATRIX_CONNECTED_COMPONENTS MPI_REJECTION PASSED\n");
        return 0;
    }
    if(!failed) err = lis_matrix_set_size(A,0,4);
    if(err) failed = 1;

    if(!failed)
    {
        for(i=0;i<4;i++)
            if(lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)2.0,A)) failed=1;
        if(lis_matrix_set_value(LIS_INS_VALUE,0,1,(LIS_SCALAR)-1.0,A)) failed=1;
        if(lis_matrix_set_value(LIS_INS_VALUE,1,0,(LIS_SCALAR)-1.0,A)) failed=1;
        if(lis_matrix_set_value(LIS_INS_VALUE,2,3,(LIS_SCALAR)0.6,A)) failed=1;
        if(lis_matrix_set_value(LIS_INS_VALUE,3,2,(LIS_SCALAR)0.6,A)) failed=1;
    }
    if(!failed) err = lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) failed=1;
    if(!failed) err = lis_matrix_assemble(A);
    if(err) failed=1;

    /* Drop tolerance keeps 0.6 and 1.0 edges: two components. */
    if(!failed) err = lis_matrix_connected_components(A,(LIS_REAL)0.5,labels,&count);
    if(err || count!=2 || !labels_same_partition(labels,expected_csr,4))
    {
        fprintf(stderr,"CSR components mismatch: err=%d count=%d\n",(int)err,(int)count);
        failed=1;
    }

    /* A larger tolerance drops 0.6 but keeps 1.0: rows 2 and 3 isolate. */
    if(!failed) err = lis_matrix_connected_components(A,(LIS_REAL)0.75,labels,&count);
    if(err || count!=3)
    {
        fprintf(stderr,"drop-tolerance components mismatch: err=%d count=%d\n",(int)err,(int)count);
        failed=1;
    }

    if(!failed && lis_matrix_connected_components(A,(LIS_REAL)-1.0,labels,&count)!=LIS_ERR_ILL_ARG)
    {
        fprintf(stderr,"negative drop tolerance was not rejected\n");
        failed=1;
    }

    /* Another explicit storage format must use the normal conversion path. */
    if(!failed) err = lis_matrix_duplicate(A,&C);
    if(err) failed=1;
    if(!failed) err = lis_matrix_set_type(C,LIS_MATRIX_CSC);
    if(err) failed=1;
    if(!failed) err = lis_matrix_convert(A,C);
    if(err) failed=1;
    if(!failed) err = lis_matrix_connected_components(C,(LIS_REAL)0.5,labels_csc,&count);
    if(err || count!=2 || !labels_same_partition(labels_csc,expected_csr,4))
    {
        fprintf(stderr,"CSC conversion components mismatch: err=%d count=%d\n",(int)err,(int)count);
        failed=1;
    }

    /* USER path: repeated duplicates, cancellation and isolated rows. */
    if(!failed) err = lis_matrix_create(LIS_COMM_WORLD,&U);
    if(err) failed=1;
    if(!failed) err = lis_matrix_set_size(U,0,4);
    if(err) failed=1;
    if(!failed) err = lis_matrix_set_user(U,&ctx,NULL,NULL);
    if(err) failed=1;
    if(!failed) err = lis_matrix_set_user_row(U,user_row);
    if(err) failed=1;
    if(!failed) err = lis_matrix_connected_components(U,(LIS_REAL)0.0,labels_user,&count);
    if(err || count!=3 || !labels_same_partition(labels_user,expected_user,4))
    {
        fprintf(stderr,"USER duplicate components mismatch: err=%d count=%d\n",(int)err,(int)count);
        failed=1;
    }

#ifdef _COMPLEX
    if(!failed) err = lis_matrix_create(LIS_COMM_WORLD,&Z);
    if(err) failed=1;
    if(!failed) err = lis_matrix_set_size(Z,0,2);
    if(err) failed=1;
    if(!failed)
    {
        if(lis_matrix_set_value(LIS_INS_VALUE,0,0,(LIS_SCALAR)1.0,Z)) failed=1;
        if(lis_matrix_set_value(LIS_INS_VALUE,1,1,(LIS_SCALAR)1.0,Z)) failed=1;
        if(lis_matrix_set_value(LIS_INS_VALUE,0,1,(LIS_SCALAR)(0.0 + 0.75*I),Z)) failed=1;
    }
    if(!failed) err = lis_matrix_set_type(Z,LIS_MATRIX_CSR);
    if(err) failed=1;
    if(!failed) err = lis_matrix_assemble(Z);
    if(err) failed=1;
    if(!failed) err = lis_matrix_connected_components(Z,(LIS_REAL)0.5,labels,&count);
    if(err || count!=1)
    {
        fprintf(stderr,"complex magnitude components mismatch: err=%d count=%d\n",(int)err,(int)count);
        failed=1;
    }
#endif

    if(!failed)
    {
        err = lis_matrix_create_operator((LIS_SCALAR)1.0,A,(LIS_SCALAR)0.0,A,&O);
        if(err) failed=1;
    }
    if(!failed && lis_matrix_connected_components(O,(LIS_REAL)0.0,labels,&count)!=LIS_ERR_NOT_IMPLEMENTED)
    {
        fprintf(stderr,"operator components were not rejected\n");
        failed=1;
    }

    if(O) lis_matrix_destroy(O);
#ifdef _COMPLEX
    if(Z) lis_matrix_destroy(Z);
#endif
    if(U) lis_matrix_destroy(U);
    if(C) lis_matrix_destroy(C);
    if(A) lis_matrix_destroy(A);
    lis_finalize();

    if(failed)
    {
        fprintf(stderr,"MATRIX_CONNECTED_COMPONENTS FAILED\n");
        return 1;
    }

    printf("MATRIX_CONNECTED_COMPONENTS PASSED\n");
    return 0;
}
