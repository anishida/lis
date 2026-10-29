#include "lis.h"
#include "lis_precon.h"
#include "lis_esolver.h"
#include <math.h>
#include <stdio.h>
#ifdef _COMPLEX
#include <complex.h>
#endif

typedef struct
{
    LIS_INT which;
} USER_CTX;

static const double amat[9] = {
     4.0,-1.0, 0.0,
    -1.0, 4.0,-1.0,
     0.0,-1.0, 3.0
};
static const double bmat[9] = {
     1.0,0.0,0.0,
     0.0,2.0,0.0,
     0.0,0.0,1.0
};

static LIS_REAL scalar_abs(LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (LIS_REAL)cabs(value);
#else
    return (LIS_REAL)fabs((double)value);
#endif
}

static int near_scalar(LIS_SCALAR value,LIS_SCALAR expected)
{
    return scalar_abs(value-expected)<(LIS_REAL)1.0e-12;
}

static LIS_INT user_row(void *user_data,LIS_INT row,LIS_INT *nnz,
                        LIS_INT *cols,LIS_SCALAR *vals)
{
    USER_CTX *ctx=(USER_CTX *)user_data;
    const double *m=ctx->which ? bmat : amat;
    LIS_INT j,count=0;

    if(!nnz || row<0 || row>=3) return LIS_ERR_ILL_ARG;
    for(j=0;j<3;j++) if(m[3*row+j]!=0.0) count++;
    if(!cols && !vals){ *nnz=count; return LIS_SUCCESS; }
    if(!cols || !vals || *nnz<count) return LIS_ERR_ILL_ARG;
    count=0;
    for(j=0;j<3;j++) if(m[3*row+j]!=0.0){
        cols[count]=j;
        vals[count]=(LIS_SCALAR)m[3*row+j];
        count++;
    }
    *nnz=count;
    return LIS_SUCCESS;
}

static LIS_INT make_user(USER_CTX *ctx,LIS_MATRIX *A)
{
    LIS_INT err;
    err=lis_matrix_create(LIS_COMM_WORLD,A); if(err) return err;
    err=lis_matrix_set_size(*A,0,3); if(err) return err;
    err=lis_matrix_set_user(*A,ctx,NULL,NULL); if(err) return err;
    err=lis_matrix_set_user_row(*A,user_row); if(err) return err;
    return lis_matrix_set_user_row_thread_safe(*A,LIS_TRUE);
}

static int run_solver(LIS_MATRIX C,const char *precon)
{
    LIS_VECTOR truth=NULL,b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err,status;
    LIS_REAL nrm;
    char options[128];
    int rc=1;

    err=lis_vector_duplicate(C,&truth); if(err) goto done;
    err=lis_vector_duplicate(C,&b); if(err) goto done;
    err=lis_vector_duplicate(C,&x); if(err) goto done;
    truth->value[0]=(LIS_SCALAR)1.0;
    truth->value[1]=(LIS_SCALAR)2.0;
    truth->value[2]=(LIS_SCALAR)3.0;
    err=lis_matvec(C,truth,b); if(err) goto done;
    lis_vector_set_all((LIS_SCALAR)0.0,x);

    err=lis_solver_create(&solver); if(err) goto done;
    snprintf(options,sizeof(options),
             "-i gmres -p %s -restart 10 -tol 1e-11 -maxiter 100 -scale none",
             precon);
    err=lis_solver_set_option(options,solver); if(err) goto done;
    err=lis_solve(C,b,x,solver); if(err) goto done;
    err=lis_solver_get_status(solver,&status); if(err || status!=LIS_SUCCESS) goto done;
    lis_vector_axpy((LIS_SCALAR)-1.0,truth,x);
    err=lis_vector_nrm2(x,&nrm); if(err || nrm>(LIS_REAL)1.0e-8) goto done;
    rc=0;

done:
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(truth) lis_vector_destroy(truth);
    return rc;
}

static int check_cancellation_structure(LIS_MATRIX A,LIS_MATRIX B)
{
    LIS_MATRIX C=NULL;
    LIS_PRECON_ROW_SOURCE row;
    LIS_INT err,k,found=0;
    int rc=1;

    err=lis_matrix_create_operator(
        (LIS_SCALAR)1.0,A,(LIS_SCALAR)-2.0,B,&C);
    if(err) goto done;

    err=lis_precon_row_source_init(C,&row); if(err) goto done;
    err=lis_precon_row_source_get(&row,1);
    if(err || row.nnz!=3){ lis_precon_row_source_destroy(&row); goto done; }

    for(k=0;k<row.nnz;k++)
    {
        if(row.index[k]==1)
        {
            if(!near_scalar(row.value[k],(LIS_SCALAR)0.0))
            {
                lis_precon_row_source_destroy(&row);
                goto done;
            }
            found=1;
        }
    }
    lis_precon_row_source_destroy(&row);
    if(!found) goto done;

    rc=0;
done:
    if(C) lis_matrix_destroy(C);
    return rc;
}

#ifdef _COMPLEX
static int check_complex_operator(LIS_MATRIX A,LIS_MATRIX B)
{
    const LIS_SCALAR alpha=(LIS_SCALAR)(1.0+0.5*I);
    const LIS_SCALAR beta=(LIS_SCALAR)(-0.25+0.125*I);
    LIS_MATRIX C=NULL;
    LIS_PRECON_ROW_SOURCE row;
    LIS_INT err;
    int rc=1;

    err=lis_matrix_create_operator(alpha,A,beta,B,&C);
    if(err) goto done;

    err=lis_precon_row_source_init(C,&row); if(err) goto done;
    err=lis_precon_row_source_get(&row,1);
    if(err || row.nnz!=3){ lis_precon_row_source_destroy(&row); goto done; }

    if(row.index[0]!=0 ||
       !near_scalar(row.value[0],(LIS_SCALAR)(-1.0-0.5*I)) ||
       row.index[1]!=1 ||
       !near_scalar(row.value[1],(LIS_SCALAR)(3.5+2.25*I)) ||
       row.index[2]!=2 ||
       !near_scalar(row.value[2],(LIS_SCALAR)(-1.0-0.5*I)))
    {
        lis_precon_row_source_destroy(&row);
        goto done;
    }
    lis_precon_row_source_destroy(&row);

    if(run_solver(C,"ilu")) goto done;
    if(run_solver(C,"ilut")) goto done;
    if(run_solver(C,"iluc")) goto done;

    rc=0;
done:
    if(C) lis_matrix_destroy(C);
    return rc;
}
#endif


static LIS_INT build_csr_path_matrix(
    LIS_INT variant,
    LIS_MATRIX *A)
{
    LIS_INT *ptr=NULL,*index=NULL;
    LIS_SCALAR *value=NULL;
    LIS_INT err;

    err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return err;
    err=lis_matrix_set_size(*A,0,4);
    if(err) return err;

    if(variant==0)
    {
        err=lis_matrix_malloc_csr(4,10,&ptr,&index,&value);
        if(err) return err;
        ptr[0]=0; ptr[1]=2; ptr[2]=5; ptr[3]=8; ptr[4]=10;
        index[0]=0; value[0]=(LIS_SCALAR)4.0;
        index[1]=1; value[1]=(LIS_SCALAR)-1.0;
        index[2]=0; value[2]=(LIS_SCALAR)-1.0;
        index[3]=1; value[3]=(LIS_SCALAR)4.0;
        index[4]=2; value[4]=(LIS_SCALAR)-1.0;
        index[5]=1; value[5]=(LIS_SCALAR)-1.0;
        index[6]=2; value[6]=(LIS_SCALAR)4.0;
        index[7]=3; value[7]=(LIS_SCALAR)-1.0;
        index[8]=2; value[8]=(LIS_SCALAR)-1.0;
        index[9]=3; value[9]=(LIS_SCALAR)3.0;
        err=lis_matrix_set_csr(10,ptr,index,value,*A);
        if(err) return err;
        return lis_matrix_assemble(*A);
    }

    err=lis_matrix_malloc_csr(4,8,&ptr,&index,&value);
    if(err) return err;
    ptr[0]=0; ptr[1]=1; ptr[2]=3; ptr[3]=6; ptr[4]=8;
    index[0]=0; value[0]=(LIS_SCALAR)1.0;
    index[1]=0; value[1]=(LIS_SCALAR)0.25;
    index[2]=1; value[2]=(LIS_SCALAR)1.0;
    index[3]=1; value[3]=(LIS_SCALAR)0.25;
    index[4]=2; value[4]=(LIS_SCALAR)1.0;
    index[5]=3; value[5]=(LIS_SCALAR)0.25;
    index[6]=2; value[6]=(LIS_SCALAR)0.25;
    index[7]=3; value[7]=(LIS_SCALAR)1.0;
    err=lis_matrix_set_csr(8,ptr,index,value,*A);
    if(err) return err;
    return lis_matrix_assemble(*A);
}

static LIS_INT build_same_pattern_b(LIS_MATRIX source,LIS_MATRIX *B)
{
    LIS_INT err,i;
    err=lis_matrix_duplicate(source,B);
    if(err) return err;
    lis_matrix_set_type(*B,LIS_MATRIX_CSR);
    err=lis_matrix_convert(source,*B);
    if(err) return err;
    for(i=0;i<(*B)->nnz;i++)
        (*B)->value[i]=(LIS_SCALAR)0.125*source->value[i];
    return LIS_SUCCESS;
}

static int try_rebuild_precon(LIS_MATRIX op,const char *name)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_INT err;
    char options[128];
    int rc=1;

    err=lis_solver_create(&solver); if(err) goto done;
    snprintf(options,sizeof(options),"-i gmres -p %s -scale none",name);
    err=lis_solver_set_option(options,solver); if(err) goto done;
    err=lis_esolver_rebuild_precon(op,solver,&precon);
    if(err || precon==NULL) goto done;
    rc=0;

done:
    if(precon)
    {
        if(solver && solver->precon==precon) solver->precon=NULL;
        lis_precon_destroy(precon);
    }
    if(solver) lis_solver_destroy(solver);
    return rc;
}

static int check_csr_operator_paths(LIS_INT different_pattern)
{
    LIS_MATRIX A=NULL,B=NULL,C=NULL;
    LIS_INT err;
    int rc=1;

    err=build_csr_path_matrix(0,&A); if(err) goto done;
    if(different_pattern)
        err=build_csr_path_matrix(1,&B);
    else
        err=build_same_pattern_b(A,&B);
    if(err) goto done;

    err=lis_matrix_create_operator(
        (LIS_SCALAR)1.0,A,(LIS_SCALAR)-0.25,B,&C);
    if(err) goto done;

    if(try_rebuild_precon(C,"ilu")) goto done;
    if(try_rebuild_precon(C,"ilut")) goto done;
    if(try_rebuild_precon(C,"iluc")) goto done;

    rc=0;
done:
    if(C) lis_matrix_destroy(C);
    if(B) lis_matrix_destroy(B);
    if(A) lis_matrix_destroy(A);
    return rc;
}

int main(int argc,char **argv)
{
    USER_CTX actx={0},bctx={1};
    LIS_MATRIX A=NULL,B=NULL,C=NULL;
    LIS_PRECON_ROW_SOURCE row;
    LIS_INT err;
    int rc=1;

    err=lis_initialize(&argc,&argv); if(err) return 1;
    err=make_user(&actx,&A); if(err) goto done;
    err=make_user(&bctx,&B); if(err) goto done;
    err=lis_matrix_create_operator((LIS_SCALAR)1.0,A,(LIS_SCALAR)-0.25,B,&C);
    if(err) goto done;

    err=lis_precon_row_source_init(C,&row); if(err) goto done;
    err=lis_precon_row_source_get(&row,1);
    if(err || row.nnz!=3) { lis_precon_row_source_destroy(&row); goto done; }
    lis_precon_row_source_destroy(&row);

    if(run_solver(C,"ilu")) goto done;
    if(run_solver(C,"ilut")) goto done;
    if(run_solver(C,"iluc")) goto done;

    if(check_cancellation_structure(A,B)) goto done;

#ifdef _COMPLEX
    if(check_complex_operator(A,B)) goto done;
#endif

    if(check_csr_operator_paths(0)) goto done;
    if(check_csr_operator_paths(1)) goto done;

    printf("LIS_OPERATOR_ROW_PRECON PASSED\n");
    rc=0;

done:
    if(C) lis_matrix_destroy(C);
    if(B) lis_matrix_destroy(B);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return rc;
}
