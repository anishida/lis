#include "lis_ldlt_user_bridge.h"
#include "symbolic_ldlt.h"
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#if !defined(USE_MPI) && !defined(_COMPLEX) && !defined(USE_COMPLEX) && !defined(_LONG__DOUBLE) && !defined(USE_QUAD_PRECISION)
typedef struct {
    p56_context *factor;
    LIS_MATRIX matrix; /* borrowed; caller keeps it alive */
    LIS_UNSIGNED_INT structure,value;
    double *scratch;
} bridge_state;

static LIS_INT drop(LIS_PRECON p)
{
    bridge_state *s=(bridge_state *)p->user_data;
    if(s){p56_destroy(s->factor);free(s->scratch);free(s);}
    p->user_data=NULL;
    return LIS_SUCCESS;
}

static LIS_INT create(LIS_SOLVER solver,LIS_PRECON p)
{
    LIS_MATRIX A=solver->A;
    LIS_INT i,k,j,lo,hi,mid,rc=LIS_ERR_ILL_ARG;
    size_t nnz,z=0,*col=NULL;
    uint32_t *row=NULL;
    double *val=NULL;
    bridge_state *s=NULL;
    const p62_limits sl={64u*1024*1024,4000000};
    const p61_limits nl={1024,128u*1024*1024,4000000,1000000,0};
    if(sizeof(LIS_SCALAR)!=sizeof(double)||!A||A->nprocs!=1||
       A->n!=A->gn||A->matrix_type!=LIS_MATRIX_CSR||A->is_scaled||
       solver->options[LIS_OPTIONS_SCALE]!=LIS_SCALE_NONE||
       solver->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE||
       solver->options[LIS_OPTIONS_ADDS]||
       (solver->options[LIS_OPTIONS_SOLVER]!=LIS_SOLVER_GMRES &&
        solver->options[LIS_OPTIONS_SOLVER]!=LIS_SOLVER_FGMRES))
        return LIS_ERR_NOT_IMPLEMENTED;
    /* Bounded first integration: do not silently reinterpret other storage,
       unsorted rows, duplicate entries or non-symmetric matrices. */
    if(A->n<=0||A->n>100000||!A->ptr||A->nnz<0||A->nnz>4000000||
       (A->nnz&&(!A->index||!A->value))||
       A->ptr[0]!=0||A->ptr[A->n]!=A->nnz) return LIS_ERR_ILL_ARG;
    nnz=(size_t)A->nnz;
    for(i=0;i<A->n;i++){
        if(A->ptr[i]<0||A->ptr[i+1]<A->ptr[i]||A->ptr[i+1]>A->nnz)
            return LIS_ERR_ILL_ARG;
        for(k=A->ptr[i];k<A->ptr[i+1];k++){
            j=A->index[k];
            if(j<0||j>=A->n||!isfinite(A->value[k])||
               (k>A->ptr[i]&&A->index[k-1]>=j))return LIS_ERR_ILL_ARG;
        }
    }
    for(i=0;i<A->n;i++)for(k=A->ptr[i];k<A->ptr[i+1];k++){
        j=A->index[k];lo=A->ptr[j];hi=A->ptr[j+1];
        while(lo<hi){mid=lo+(hi-lo)/2;if(A->index[mid]<i)lo=mid+1;else hi=mid;}
        if(lo==A->ptr[j+1]||A->index[lo]!=i||A->value[lo]!=A->value[k])
            return LIS_ERR_ILL_ARG;
    }
    s=calloc(1,sizeof *s);
    col=malloc(((size_t)A->n+1)*sizeof *col);
    row=malloc((nnz?nnz:1)*sizeof *row);
    val=malloc((nnz?nnz:1)*sizeof *val);
    if(!s||!col||!row||!val){rc=LIS_OUT_OF_MEMORY;goto done;}
    s->scratch=malloc((size_t)A->n*sizeof *s->scratch);
    if(!s->scratch){rc=LIS_OUT_OF_MEMORY;goto done;}
    /* For an exactly symmetric CSR, upper row i is lower CSC column i. */
    for(i=0;i<A->n;i++){
        col[i]=z;
        for(k=A->ptr[i];k<A->ptr[i+1];k++)if(A->index[k]>=i){
            row[z]=(uint32_t)A->index[k];val[z++]=A->value[k];
        }
    }
    col[A->n]=z;
    {
        int e=p62_factorize(&s->factor,(uint32_t)A->n,z,col,row,val,&sl,&nl,NULL);
        if(e){rc=e==P55_NOMEM?LIS_OUT_OF_MEMORY:LIS_BREAKDOWN;goto done;}
    }
    s->matrix=A;
    rc=lis_matrix_get_generation(A,&s->structure,&s->value);
    if(rc)goto done;
    p->user_data=s;s=NULL;rc=LIS_SUCCESS;
done:
    if(s){p56_destroy(s->factor);free(s->scratch);free(s);}
    free(col);free(row);free(val);
    return rc;
}

static LIS_INT apply(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    bridge_state *s;
    LIS_UNSIGNED_INT sg,vg;
    if(!solver||!solver->precon||!b||!x)return LIS_ERR_ILL_ARG;
    s=(bridge_state *)solver->precon->user_data;
    if(!s||solver->A!=s->matrix||b->n!=s->matrix->n||x->n!=b->n||
       !b->value||!x->value)return LIS_ERR_ILL_ARG;
    if(lis_matrix_get_generation(s->matrix,&sg,&vg)||
       sg!=s->structure||vg!=s->value)return LIS_ERR_ILL_ARG;
    return p56_apply(s->factor,b->value,x->value,s->scratch,(size_t)b->n)
        ?LIS_BREAKDOWN:LIS_SUCCESS;
}

LIS_INT lis_ldlt_experimental_register(void)
{
    return lis_precon_register_ex("ldltexp",create,apply,apply,drop);
}
#else
LIS_INT lis_ldlt_experimental_register(void){return LIS_ERR_NOT_IMPLEMENTED;}
#endif
