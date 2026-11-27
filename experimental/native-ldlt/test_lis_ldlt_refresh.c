#include "lis_ldlt_user_bridge.h"
#include "lis_precon.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"FAIL %d %s\n",__LINE__,#x);exit(1);}}while(0)
static void solve(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_VECTOR r,LIS_SOLVER s,LIS_PRECON p){
 CHECK(!lis_vector_set_all(0,x));CHECK(!lis_solve_kernel(A,b,x,s,p));LIS_INT status;CHECK(!lis_solver_get_status(s,&status)&&status==LIS_SUCCESS);CHECK(!lis_matvec(A,x,r));CHECK(!lis_vector_axpy(-1,b,r));LIS_REAL rn,bn;CHECK(!lis_vector_nrm2(r,&rn));CHECK(!lis_vector_nrm2(b,&bn));CHECK(isfinite(rn)&&rn/bn<1e-12);
}
int main(int argc,char **argv){
 LIS_MATRIX A;LIS_SOLVER s;LIS_PRECON p=NULL;LIS_VECTOR b,x,r;
 CHECK(!lis_initialize(&argc,&argv));CHECK(!lis_ldlt_experimental_register());CHECK(!lis_ldlt_experimental_set_threads(1));CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,4));
 for(LIS_INT i=0;i<4;i++)for(LIS_INT j=0;j<4;j++)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,j,i==j?0:i/2==j/2?4.:0.,A));
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));CHECK(!lis_vector_duplicate(A,&b));CHECK(!lis_vector_duplicate(A,&x));CHECK(!lis_vector_duplicate(A,&r));CHECK(!lis_vector_set_all(1,b));CHECK(!lis_solver_create(&s));CHECK(!lis_solver_set_option("-i fgmres -p ldltexpr -scale none -tol 1e-12 -maxiter 20",s));s->A=A;CHECK(!lis_precon_create(s,&p));solve(A,b,x,r,s,p);
 /* Refresh after values change. The new matrix may choose different pivots. */
 for(LIS_INT i=0;i<4;i++)for(LIS_INT z=A->ptr[i];z<A->ptr[i+1];z++)A->value[z]=A->index[z]==i?(i%2?5.:-5.):A->index[z]/2==i/2?1.:0.;
 CHECK(!lis_matrix_mark_values_changed(A));s->precon=p;CHECK(!lis_vector_set_all(7,x));CHECK(lis_psolve_dispatch(s,b,x)==LIS_ERR_ILL_ARG);for(LIS_INT i=0;i<4;i++)CHECK(x->value[i]==7.);
 CHECK(!lis_ldlt_experimental_refresh(s,p));solve(A,b,x,r,s,p);
 /* Failed numeric refresh must leave the existing object stale; recovery
    from a later value generation is possible using the same plan. */
 for(LIS_INT z=0;z<A->nnz;z++)A->value[z]=0.;
 CHECK(!lis_matrix_mark_values_changed(A));CHECK(lis_ldlt_experimental_refresh(s,p)!=LIS_SUCCESS);
 s->precon=p;CHECK(!lis_vector_set_all(7,x));CHECK(lis_psolve_dispatch(s,b,x)==LIS_ERR_ILL_ARG);for(LIS_INT i=0;i<4;i++)CHECK(x->value[i]==7.);
 for(LIS_INT i=0;i<4;i++)for(LIS_INT z=A->ptr[i];z<A->ptr[i+1];z++)A->value[z]=A->index[z]==i?(i%2?8.:-8.):0.;
 CHECK(!lis_matrix_mark_values_changed(A));CHECK(!lis_ldlt_experimental_refresh(s,p));solve(A,b,x,r,s,p);
 CHECK(!lis_matrix_mark_structure_changed(A));CHECK(lis_ldlt_experimental_refresh(s,p)==LIS_ERR_ILL_ARG);
 if(s->precon==p)s->precon=NULL;
 lis_precon_destroy(p);lis_solver_destroy(s);lis_vector_destroy(b);lis_vector_destroy(x);lis_vector_destroy(r);lis_matrix_destroy(A);lis_finalize();puts("LDLT_LIS_REFRESH=PASS");return 0;
}
