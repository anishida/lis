#include "lis_ldlt_user_bridge.h"
#include "qmd_plan.h"
#include "schur_simd.h"
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#if !defined(USE_MPI) && !defined(_COMPLEX) && !defined(USE_COMPLEX) && !defined(_LONG__DOUBLE) && !defined(USE_QUAD_PRECISION)
typedef struct {
    p56_context *factor;
    p71_plan *plan;
    LIS_MATRIX matrix; /* borrowed; caller keeps it alive */
    LIS_UNSIGNED_INT structure,value;
    double *scratch;
} bridge_state;

static LIS_INT drop(LIS_PRECON p)
{
    bridge_state *s=(bridge_state *)p->user_data;
    if(s){p56_destroy(s->factor);p71_plan_free(s->plan);free(s->scratch);free(s);}
    p->user_data=NULL;
    return LIS_SUCCESS;
}

static const p63_limits ol={512u*1024*1024,64000000,2000000000};
static const p62_limits sl={256u*1024*1024,64000000};
static const p61_limits nl={8192,1536u*1024*1024,64000000,1000000,0};
static LIS_INT validate(LIS_SOLVER solver)
{
    if(!solver)return LIS_ERR_ILL_ARG;
    LIS_MATRIX A=solver->A;
    LIS_INT i,k,j,lo,hi,mid;
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
    if(A->n<=0||A->n>200000||!A->ptr||A->nnz<0||A->nnz>10000000||
       (A->nnz&&(!A->index||!A->value))||
       A->ptr[0]!=0||A->ptr[A->n]!=A->nnz) return LIS_ERR_ILL_ARG;
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
    return LIS_SUCCESS;
}
static LIS_INT convert(LIS_MATRIX A,size_t **cp,uint32_t **ri,double **values,size_t *count)
{
    size_t nnz=(size_t)A->nnz,z=0;
    size_t *col=malloc(((size_t)A->n+1)*sizeof *col);
    uint32_t *row=malloc((nnz?nnz:1)*sizeof *row);
    double *val=malloc((nnz?nnz:1)*sizeof *val);
    if(!col||!row||!val){free(col);free(row);free(val);return LIS_OUT_OF_MEMORY;}
    for(LIS_INT i=0;i<A->n;i++){
        col[i]=z;
        for(LIS_INT k=A->ptr[i];k<A->ptr[i+1];k++)if(A->index[k]>=i){row[z]=(uint32_t)A->index[k];val[z++]=A->value[k];}
    }
    col[A->n]=z;*cp=col;*ri=row;*values=val;*count=z;return LIS_SUCCESS;
}
static LIS_INT factor_status(int e)
{
    if(!e)return LIS_SUCCESS;
    fprintf(stderr,"LDLT factorization status=%d\n",e);
    return e==P55_NOMEM?LIS_OUT_OF_MEMORY:e==P55_INVALID?LIS_ERR_ILL_ARG:LIS_BREAKDOWN;
}
static LIS_INT create_mode(LIS_SOLVER solver,LIS_PRECON p,int reuse)
{
    LIS_INT rc=validate(solver);if(rc)return rc;
    LIS_MATRIX A=solver->A;
    size_t z=0,*col=NULL;uint32_t *row=NULL;double *val=NULL;
    bridge_state *s=calloc(1,sizeof *s);if(!s)return LIS_OUT_OF_MEMORY;
    s->scratch=malloc((size_t)A->n*sizeof *s->scratch);
    if(!s->scratch){rc=LIS_OUT_OF_MEMORY;goto done;}
    rc=convert(A,&col,&row,&val,&z);if(rc)goto done;
    p63_stats st={0};int e;
    if(reuse){
        e=p71_plan_create(&s->plan,(uint32_t)A->n,z,col,row,&ol,&sl,32,256u*1024*1024);
        if(!e)e=p71_refactor(s->plan,&s->factor,(uint32_t)A->n,z,col,row,val,&nl,&st);
    }else e=p63q_factorize(&s->factor,(uint32_t)A->n,z,col,row,val,&ol,&sl,&nl,32,&st);
    rc=factor_status(e);if(rc)goto done;
    fprintf(stderr,"QMD slots=%zu front=%zu visits=%zu numeric_slots=%zu delays=%zu peak=%zu\n",st.summary.slots,st.summary.max_front,st.summary.pair_checks,st.numeric.factor_slots,st.numeric.delays,st.numeric.peak_work_bytes);
    s->matrix=A;
    rc=lis_matrix_get_generation(A,&s->structure,&s->value);
    if(rc)goto done;
    p->user_data=s;s=NULL;rc=LIS_SUCCESS;
done:
    if(s){p56_destroy(s->factor);p71_plan_free(s->plan);free(s->scratch);free(s);}
    free(col);free(row);free(val);return rc;
}
static LIS_INT create(LIS_SOLVER s,LIS_PRECON p){return create_mode(s,p,0);}
static LIS_INT create_reuse(LIS_SOLVER s,LIS_PRECON p){return create_mode(s,p,1);}
LIS_INT lis_ldlt_experimental_refresh(LIS_SOLVER solver,LIS_PRECON p)
{
    if(!p||p->user_destroy!=drop||!p->user_data)return LIS_ERR_ILL_ARG;
    bridge_state *s=p->user_data;
    if(!solver||!s->plan||solver->A!=s->matrix)return LIS_ERR_ILL_ARG;
    LIS_UNSIGNED_INT sg,vg;
    if(lis_matrix_get_generation(s->matrix,&sg,&vg)||sg!=s->structure)return LIS_ERR_ILL_ARG;
    LIS_INT rc=validate(solver);if(rc)return rc;
    if(vg==s->value)return LIS_SUCCESS;
    size_t z=0,*col=NULL;uint32_t *row=NULL;double *val=NULL;
    rc=convert(s->matrix,&col,&row,&val,&z);if(rc)return rc;
    p63_stats st={0};int e=p71_refactor(s->plan,&s->factor,(uint32_t)s->matrix->n,z,col,row,val,&nl,&st);
    free(col);free(row);free(val);rc=factor_status(e);if(rc)return rc;
    s->value=vg;
    fprintf(stderr,"LDLT refresh slots=%zu delays=%zu peak=%zu\n",st.numeric.factor_slots,st.numeric.delays,st.numeric.peak_work_bytes);
    return LIS_SUCCESS;
}
LIS_INT lis_ldlt_experimental_set_threads(unsigned n){return p66_set_threads(n)?LIS_SUCCESS:LIS_ERR_NOT_IMPLEMENTED;}
static LIS_INT apply(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    bridge_state *s;
    LIS_UNSIGNED_INT sg,vg;
    if(!solver||!solver->precon||!b||!x)return LIS_ERR_ILL_ARG;
    s=(bridge_state *)solver->precon->user_data;
    if(!s||solver->A!=s->matrix||b->n!=s->matrix->n||x->n!=b->n||
       !b->value||!x->value)return LIS_ERR_ILL_ARG;
    if(lis_matrix_get_generation(s->matrix,&sg,&vg)||
       sg!=s->structure||vg!=s->value){fprintf(stderr,"LDLT stale factor\n");return LIS_ERR_ILL_ARG;}
    return p56_apply(s->factor,b->value,x->value,s->scratch,(size_t)b->n)
        ?LIS_BREAKDOWN:LIS_SUCCESS;
}

LIS_INT lis_ldlt_experimental_register(void)
{
    LIS_INT rc=lis_precon_register_ex("ldltexp",create,apply,apply,drop);
    return rc?rc:lis_precon_register_ex("ldltexpr",create_reuse,apply,apply,drop);
}
#else
LIS_INT lis_ldlt_experimental_register(void){return LIS_ERR_NOT_IMPLEMENTED;}
LIS_INT lis_ldlt_experimental_refresh(LIS_SOLVER s,LIS_PRECON p){(void)s;(void)p;return LIS_ERR_NOT_IMPLEMENTED;}
LIS_INT lis_ldlt_experimental_set_threads(unsigned n){(void)n;return LIS_ERR_NOT_IMPLEMENTED;}
#endif
