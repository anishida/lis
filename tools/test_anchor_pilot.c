/* Deterministic native CSR regression for the experimental policy. */
#include "lislib.h"
#include "lis_anchor_pilot.h"
#include <math.h>
#include <stdio.h>

#if defined(USE_MPI) || defined(_COMPLEX) || defined(USE_COMPLEX) || defined(_LONG__DOUBLE) || defined(USE_QUAD_PRECISION)
int main(void){fputs("Unsupported configuration\n",stderr);return 77;}
#else
static LIS_INT diagonal(int kind,LIS_MATRIX *A)
{
    LIS_INT i,rc;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,A)) || (rc=lis_matrix_set_size(*A,0,80)))return rc;
    for(i=0;i<80;i++){
        double d=kind==0?1:kind==1?(double)(i+1):pow(10.0,-8.0+8.0*i/79.0);
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,d,*A)))return rc;
    }
    if((rc=lis_matrix_set_type(*A,LIS_MATRIX_CSR)))return rc;
    return lis_matrix_assemble(*A);
}
int main(int argc,char **argv)
{
    LIS_MATRIX P=NULL,A=NULL,B=NULL;LIS_VECTOR b=NULL,x=NULL,r=NULL;
    lis_anchor_pilot *ctx=NULL;lis_anchor_pilot_report report;
    LIS_REAL nr,nb;LIS_INT i;int fail=1;
    if(lis_initialize(&argc,&argv))return 1;
    if(diagonal(0,&P)||diagonal(1,&A)||diagonal(2,&B)||
       lis_vector_duplicate(P,&b)||lis_vector_duplicate(P,&x)||lis_vector_duplicate(P,&r)||
       lis_vector_set_all(1,b)||lis_anchor_pilot_create(P,&ctx))goto done;
    if(lis_anchor_pilot_solve(ctx,P,b,x,&report)||!report.accepted||report.factorizations||report.fallback)goto done;
    if(lis_anchor_pilot_solve(ctx,A,b,x,&report)||!report.accepted||!report.continued||report.fallback||report.factorizations)goto done;
    if(lis_matvec(A,x,r)||lis_vector_axpy(-1,b,r)||lis_vector_nrm2(r,&nr)||lis_vector_nrm2(b,&nb)||!(nr/nb<=1e-8))goto done;
    if(lis_anchor_pilot_solve(ctx,B,b,x,&report)||!report.accepted||!report.fallback||report.continued||report.factorizations!=1)goto done;
    if(lis_anchor_fresh_solve(B,b,x,&report)||!report.accepted||report.factorizations!=1)goto done;
    /* Fresh fallback did not replace the anchor. */
    if(lis_anchor_pilot_solve(ctx,P,b,x,&report)||!report.accepted||report.factorizations||report.fallback)goto done;
    lis_anchor_pilot_destroy(ctx);ctx=NULL;
    if(lis_anchor_pilot_create_low_memory(P,&ctx))goto done;
    if(lis_anchor_pilot_solve(ctx,A,b,x,&report)||!report.accepted||!report.continued||report.fallback||report.factorizations||report.pilot_iterations!=20)goto done;
    if(lis_anchor_pilot_solve(ctx,B,b,x,&report)||!report.accepted||!report.fallback||report.factorizations!=1)goto done;
    /* The next call reconstructs the original anchor exactly once. */
    if(lis_anchor_pilot_solve(ctx,P,b,x,&report)||!report.accepted||report.fallback||report.factorizations!=1)goto done;
    if(lis_anchor_pilot_solve(ctx,P,b,x,&report)||!report.accepted||report.factorizations)goto done;
    /* Also test invalidation while the anchor has no live factor. */
    if(lis_anchor_pilot_solve(ctx,B,b,x,&report)||!report.accepted||!report.fallback||report.factorizations!=1)goto done;
    if(lis_vector_set_all(7,x)||lis_matrix_mark_values_changed(P))goto done;
    if(lis_anchor_pilot_solve(ctx,A,b,x,&report)!=LIS_ERR_ILL_ARG||report.accepted)goto done;
    for(i=0;i<x->n;i++)if(x->value[i]!=7)goto done;
    if(lis_anchor_fresh_solve(A,b,b,&report)!=LIS_ERR_ILL_ARG||report.accepted)goto done;
    for(i=0;i<b->n;i++)if(b->value[i]!=1)goto done;
    lis_anchor_pilot_destroy(ctx);ctx=NULL;
    if(lis_anchor_pilot_create_rolling(P,&ctx))goto done;
    if(lis_anchor_pilot_solve(ctx,A,b,x,&report)||!report.accepted||!report.continued||report.fallback||report.factorizations||report.pilot_iterations!=20||lis_anchor_pilot_current_anchor(ctx)!=P)goto done;
    if(lis_anchor_pilot_solve(ctx,B,b,x,&report)||!report.accepted||!report.fallback||report.factorizations!=1||lis_anchor_pilot_current_anchor(ctx)!=B)goto done;
    lis_matrix_destroy(P);P=NULL; /* old anchor may now be released */
    if(lis_anchor_pilot_solve(ctx,B,b,x,&report)||!report.accepted||report.fallback||report.factorizations||lis_anchor_pilot_current_anchor(ctx)!=B)goto done;
    if(lis_matrix_mark_values_changed(B)||lis_vector_set_all(7,x))goto done;
    if(lis_anchor_pilot_solve(ctx,A,b,x,&report)!=LIS_ERR_ILL_ARG||report.accepted)goto done;
    for(i=0;i<x->n;i++)if(x->value[i]!=7)goto done;
    fail=0;
done:
    lis_anchor_pilot_destroy(ctx);
    if(r)lis_vector_destroy(r);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    if(P)lis_matrix_destroy(P);
    lis_finalize();
    printf("ANCHOR_PILOT_REGRESSION=%s\n",fail?"FAIL":"PASS");
    return fail;
}
#endif
