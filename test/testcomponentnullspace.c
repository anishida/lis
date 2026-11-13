#include "lis.h"
#include "lis_nullspace_candidates.h"
#include <stdio.h>

static LIS_INT build_matrix(LIS_MATRIX *A)
{
    LIS_INT *ptr=NULL,*index=NULL;
    LIS_SCALAR *value=NULL;
    LIS_INT err;

    err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return err;
    err=lis_matrix_set_size(*A,0,4);
    if(err) return err;
    err=lis_matrix_malloc_csr(4,8,&ptr,&index,&value);
    if(err) return err;

    ptr[0]=0; ptr[1]=2; ptr[2]=4; ptr[3]=6; ptr[4]=8;
    index[0]=0; value[0]=(LIS_SCALAR) 1.0;
    index[1]=1; value[1]=(LIS_SCALAR)-1.0;
    index[2]=0; value[2]=(LIS_SCALAR)-1.0;
    index[3]=1; value[3]=(LIS_SCALAR) 1.0;
    index[4]=2; value[4]=(LIS_SCALAR) 1.0;
    index[5]=3; value[5]=(LIS_SCALAR)-1.0;
    index[6]=2; value[6]=(LIS_SCALAR)-1.0;
    index[7]=3; value[7]=(LIS_SCALAR) 1.0;

    err=lis_matrix_set_csr(8,ptr,index,value,*A);
    if(err) return err;
    return lis_matrix_assemble(*A);
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR candidate=NULL,q=NULL,aq=NULL;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR candidates[1];
    LIS_NULLSPACE_CANDIDATE_DIAGNOSTICS ordinary;
    LIS_COMPONENT_NULLSPACE_DIAGNOSTICS component;
    LIS_REAL norm;
    LIS_INT dim=0,i,err;
    int rc=1;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=build_matrix(&A);
    if(err) goto done;
    err=lis_vector_duplicate(A,&candidate);
    if(err) goto done;
    err=lis_vector_set_all((LIS_SCALAR)1.0,candidate);
    if(err) goto done;
    candidates[0]=candidate;

    err=lis_solver_create(&solver);
    if(err) goto done;

    err=lis_solver_classify_nullspace_candidates(
        A,solver,1,candidates,NULL,(LIS_REAL)2.0,
        (LIS_REAL)1.0e-12,(LIS_REAL)1.0e-8,0,&ordinary);
    if(err) goto done;
    if(ordinary.exact_count!=1)
    {
        fprintf(stderr,"ordinary classifier expected 1 exact mode, got %d\n",
            (int)ordinary.exact_count);
        goto done;
    }

    err=lis_solver_classify_component_nullspace_candidates(
        A,solver,1,candidates,(LIS_REAL)0.0,(LIS_REAL)2.0,
        (LIS_REAL)1.0e-12,(LIS_REAL)1.0e-8,0,&component);
    if(err) goto done;

    if(component.component_count!=2 ||
       component.localized_candidate_count!=2 ||
       component.classifier.exact_count!=2)
    {
        fprintf(stderr,
            "component classifier expected components=2 localized=2 exact=2; got %d %d %d\n",
            (int)component.component_count,
            (int)component.localized_candidate_count,
            (int)component.classifier.exact_count);
        goto done;
    }

    err=lis_solver_get_nullspace_dim(solver,&dim);
    if(err || dim!=2) goto done;

    err=lis_vector_duplicate(A,&q);
    if(err) goto done;
    err=lis_vector_duplicate(A,&aq);
    if(err) goto done;

    for(i=0;i<dim;i++)
    {
        err=lis_solver_copy_nullspace_vector(solver,i,q);
        if(err) goto done;
        err=lis_matvec(A,q,aq);
        if(err) goto done;
        err=lis_vector_nrm2(aq,&norm);
        if(err || norm>(LIS_REAL)1.0e-11) goto done;
    }

    printf("LIS_COMPONENT_NULLSPACE_CANDIDATES PASSED\n");
    rc=0;

done:
    if(aq) lis_vector_destroy(aq);
    if(q) lis_vector_destroy(q);
    if(solver) lis_solver_destroy(solver);
    if(candidate) lis_vector_destroy(candidate);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return rc;
}
