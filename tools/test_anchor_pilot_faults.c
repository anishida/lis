/* GNU ld --wrap: deterministic driver-boundary failures, not an exhaustive
 * native allocation sweep. Linked against real LIS in the dedicated CI. */
#include "lislib.h"
#include "lis_anchor_pilot.h"
#include <math.h>
#include <stdio.h>

static LIS_MATRIX candidate;
static LIS_VECTOR output;
static int fault, reached, fresh_done;
enum { FAIL_FACTOR=1, FAIL_SOLVE, FAIL_GENERATION, FAIL_COPY, FAIL_NUMERIC };
LIS_INT __real_lis_solver_sequence_prepare_precon(LIS_MATRIX,LIS_SOLVER,LIS_PRECON *);
LIS_INT __real_lis_solve_kernel(LIS_MATRIX,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER,LIS_PRECON);
LIS_INT __real_lis_matrix_get_generation(LIS_MATRIX,LIS_UNSIGNED_INT *,LIS_UNSIGNED_INT *);
LIS_INT __real_lis_vector_copy(LIS_VECTOR,LIS_VECTOR);

LIS_INT __wrap_lis_solver_sequence_prepare_precon(LIS_MATRIX A,LIS_SOLVER s,LIS_PRECON *pc)
{
    if(fault==FAIL_FACTOR && A==candidate){reached++;return LIS_ERR_OUT_OF_MEMORY;}
    return __real_lis_solver_sequence_prepare_precon(A,s,pc);
}
LIS_INT __wrap_lis_solve_kernel(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_SOLVER s,LIS_PRECON pc)
{
    LIS_INT rc;
    int fresh=A==candidate && s->sequence_precon_matrix==candidate;
    if(fresh && fault==FAIL_SOLVE){reached++;return LIS_ERR_OUT_OF_MEMORY;}
    rc=__real_lis_solve_kernel(A,b,x,s,pc);
    if(fresh && !rc){
        fresh_done=1;
        if(fault==FAIL_NUMERIC){
            /* Force a numerical rejection with a finite but wrong iterate. */
            reached++;s->retcode=LIS_MAXITER;
            lis_vector_set_all(0,x);
        }
    }
    return rc;
}
LIS_INT __wrap_lis_matrix_get_generation(LIS_MATRIX A,LIS_UNSIGNED_INT *s,LIS_UNSIGNED_INT *v)
{
    if(fault==FAIL_GENERATION && fresh_done && A==candidate){reached++;return LIS_ERR_ILL_ARG;}
    return __real_lis_matrix_get_generation(A,s,v);
}
LIS_INT __wrap_lis_vector_copy(LIS_VECTOR src,LIS_VECTOR dst)
{
    if(fault==FAIL_COPY && dst==output){reached++;return LIS_ERR_ILL_ARG;}
    return __real_lis_vector_copy(src,dst);
}
static LIS_INT diagonal(int hard,LIS_MATRIX *A)
{
    LIS_INT i,rc;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,A)) || (rc=lis_matrix_set_size(*A,0,80)))return rc;
    for(i=0;i<80;i++){
        double d=hard?pow(10.0,-8.0+8.0*i/79.0):1.0;
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,d,*A)))return rc;
    }
    if((rc=lis_matrix_set_type(*A,LIS_MATRIX_CSR)))return rc;
    return lis_matrix_assemble(*A);
}
static int run(int which)
{
    LIS_MATRIX P=NULL,A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    lis_anchor_pilot *ctx=NULL;
    lis_anchor_pilot_report r;
    LIS_INT rc,i,expected;
    int failed=1;
    fault=reached=fresh_done=0;candidate=NULL;output=NULL;
    if(diagonal(0,&P)||diagonal(1,&A)||lis_vector_duplicate(P,&b)||
       lis_vector_duplicate(P,&x)||lis_vector_set_all(1,b)||
       lis_vector_set_all(7,x)||lis_anchor_pilot_create_rolling(P,&ctx))goto done;
    candidate=A;output=x;fault=which;
    expected=which==FAIL_FACTOR||which==FAIL_SOLVE?LIS_ERR_OUT_OF_MEMORY:
             which==FAIL_NUMERIC?LIS_MAXITER:LIS_ERR_ILL_ARG;
    rc=lis_anchor_pilot_solve(ctx,A,b,x,&r);
    fault=0;
    if(rc!=expected||reached!=1||r.accepted||!r.fallback||
       lis_anchor_pilot_current_anchor(ctx)!=P)goto done;
    for(i=0;i<x->n;i++)if(x->value[i]!=7 || b->value[i]!=1)goto done;
    /* The failed candidate is not retained and can be released immediately. */
    lis_matrix_destroy(A);A=NULL;candidate=NULL;
    rc=lis_anchor_pilot_solve(ctx,P,b,x,&r);
    if(which==FAIL_NUMERIC){
        /* Numerical rejection leaves the old anchor eligible for rebuilding. */
        if(rc||!r.accepted||r.factorizations!=1||r.fallback)goto done;
    }else if(rc!=LIS_ERR_ILL_ARG||r.accepted)goto done;
    failed=0;
done:
    fault=0;candidate=NULL;output=NULL;
    lis_anchor_pilot_destroy(ctx);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    if(P)lis_matrix_destroy(P);
    printf("ROLLING_FAULT_%d=%s\n",which,failed?"FAIL":"PASS");
    return failed;
}
int main(int argc,char **argv)
{
    int i,fail=0;
    if(lis_initialize(&argc,&argv))return 1;
    for(i=FAIL_FACTOR;i<=FAIL_NUMERIC;i++)fail|=run(i);
    lis_finalize();
    return fail;
}
