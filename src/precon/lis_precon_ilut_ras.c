/*
 * Opt-in serial restricted additive Schwarz ILUT with diagonal-approximate
 * interface Schur correction. No MPI exchange: distributed matrices are
 * deliberately rejected, not silently treated as serial subdomains.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "lislib.h"

#define __FUNC__ "lis_precon_ilut_ras"
#if defined(_COMPLEX) || defined(USE_COMPLEX)
/* Explicitly unsupported: never drop imaginary parts in contact-Schur math. */
LIS_INT lis_precon_create_ilut_ras(LIS_SOLVER solver,LIS_PRECON precon)
{
    (void)solver; (void)precon;
    return LIS_ERR_NOT_IMPLEMENTED;
}
LIS_INT lis_psolve_ilut_ras(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    (void)solver; (void)b; (void)x;
    return LIS_ERR_NOT_IMPLEMENTED;
}
void lis_precon_ilut_ras_destroy(struct LIS_ILUT_RAS_STRUCT *s)
{
    (void)s;
}
#else
#define LIS_ILUT_RAS_SCHUR_LIMIT 256

struct LIS_ILUT_RAS_PART {
    LIS_INT start,end,nlocal;
    LIS_INT *global;         /* local index -> original row */
    LIS_MATRIX A;
    LIS_SOLVER solver;
    LIS_PRECON pc;
    LIS_VECTOR b,x;
};

struct LIS_ILUT_RAS_STRUCT {
    LIS_MATRIX A;            /* borrowed outer serial CSR */
    LIS_INT n,ndomains,overlap;
    struct LIS_ILUT_RAS_PART *parts;
    LIS_INT *owner;
    LIS_INT *gamma,*gamma_of,*pivots;
    LIS_INT ngamma;
    LIS_SCALAR *diag_inv,*schur_lu,*schur_rhs,*schur_delta;
    LIS_VECTOR residual;
};

void lis_precon_ilut_ras_destroy(struct LIS_ILUT_RAS_STRUCT *s)
{
    LIS_INT d;
    if(!s) return;
    if(s->parts)
    for(d=0;d<s->ndomains;d++)
    {
        struct LIS_ILUT_RAS_PART *p=&s->parts[d];
        if(p->pc)
        {
            if(p->solver) p->solver->precon=NULL;
            lis_precon_destroy(p->pc);
        }
        if(p->solver) lis_solver_destroy(p->solver);
        if(p->b) lis_vector_destroy(p->b);
        if(p->x) lis_vector_destroy(p->x);
        if(p->A) lis_matrix_destroy(p->A);
        lis_free(p->global);
    }
    lis_free(s->parts);
    lis_free(s->owner);
    lis_free(s->gamma);
    lis_free(s->gamma_of);
    lis_free(s->pivots);
    lis_free(s->diag_inv);
    lis_free(s->schur_lu);
    lis_free(s->schur_rhs);
    lis_free(s->schur_delta);
    if(s->residual) lis_vector_destroy(s->residual);
    lis_free(s);
}

static LIS_INT
lis_ilut_ras_make_subdomain(LIS_SOLVER outer, struct LIS_ILUT_RAS_STRUCT *s,
                            struct LIS_ILUT_RAS_PART *p)
{
    LIS_MATRIX A=s->A;
    LIS_INT n=s->n, depth,i,k,j,t,count=0,err;
    LIS_INT *selected=NULL,*local=NULL;
    LIS_INT new_start=p->start,new_end=p->end;

    selected=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"RAS selection");
    local=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"RAS global lookup");
    if(!selected || !local) { err=LIS_OUT_OF_MEMORY; goto done; }
    memset(selected,0,(size_t)n*sizeof(LIS_INT));
    for(i=new_start;i<new_end;i++) selected[i]=1;
    for(depth=1;depth<=s->overlap;depth++)
    for(i=0;i<n;i++)
    {
        if(selected[i]!=depth) continue;
        for(k=A->ptr[i];k<A->ptr[i+1];k++)
        {
            j=A->index[k];
            if(j<0 || j>=n) { err=LIS_ERR_ILL_ARG; goto done; }
            if(!selected[j]) selected[j]=depth+1;
        }
    }
    for(i=0;i<n;i++) if(selected[i]) count++;
    if(count<=0) { err=LIS_ERR_ILL_ARG; goto done; }
    p->global=(LIS_INT *)lis_malloc((size_t)count*sizeof(LIS_INT),"RAS local rows");
    if(!p->global) { err=LIS_OUT_OF_MEMORY; goto done; }
    p->nlocal=count;
    count=0;
    for(i=0;i<n;i++)
    {
        local[i]=-1;
        if(selected[i]) { local[i]=count; p->global[count++]=i; }
    }

    err=lis_matrix_create(LIS_COMM_WORLD,&p->A);
    if(err) goto done;
    err=lis_matrix_set_size(p->A,0,p->nlocal);
    if(err) goto done;
    err=lis_matrix_set_type(p->A,LIS_MATRIX_CSR);
    if(err) goto done;
    for(i=0;i<p->nlocal;i++)
    {
        t=p->global[i];
        for(k=A->ptr[t];k<A->ptr[t+1];k++)
        {
            j=local[A->index[k]];
            if(j<0) continue;
            err=lis_matrix_set_value(LIS_INS_VALUE,i,j,A->value[k],p->A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(p->A);
    if(err) goto done;

    err=lis_solver_create(&p->solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i gmres -p ilut -scale none -print none",p->solver);
    if(err) goto done;
    p->solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
    p->solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
    p->solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
    p->solver->params[LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN];
    /*
     * ILUTP is a local factorization policy, not a global RAS Schur
     * algorithm. Forward it explicitly, including its pivot threshold
     * and bounded full-fill controls; otherwise -ilutp 1 on the outer
     * solver would silently select plain local ILUT.
     */
    p->solver->ilutp_enabled=outer->ilutp_enabled;
    p->solver->ilutp_threshold=outer->ilutp_threshold;
    p->solver->ilutp_full=outer->ilutp_full;
    p->solver->ilutp_full_limit=outer->ilutp_full_limit;
    p->solver->ilut_rcm=outer->ilut_rcm;
    p->solver->ilut_row_scale=outer->ilut_row_scale;
    err=lis_solver_set_precon_half(p->solver,outer->precon_half_mode);
    if(err) goto done;
    err=lis_solver_set_matrix(p->A,p->solver);
    if(err) goto done;
    err=lis_precon_create(p->solver,&p->pc);
    if(err) { p->pc=NULL; goto done; }
    p->solver->precon=p->pc;
    err=lis_vector_duplicate(p->A,&p->b);
    if(err) goto done;
    err=lis_vector_duplicate(p->A,&p->x);

done:
    lis_free(selected);
    lis_free(local);
    return err;
}

/* Build S~=A_GG - A_GI diag(A_II)^-1 A_IG, a documented
 * approximate block Schur complement, never the exact Schur matrix.
 */
static LIS_INT
lis_ilut_ras_make_schur(struct LIS_ILUT_RAS_STRUCT *s)
{
    LIS_MATRIX A=s->A;
    LIS_INT i,j,k,g,q,t,p,best,ng=0;
    LIS_REAL val,peak,schur_norm=0.0;
    LIS_SCALAR tmp;
    size_t entries;

    s->gamma_of=(LIS_INT *)lis_malloc((size_t)s->n*sizeof(LIS_INT),
                                      "RAS gamma map");
    s->diag_inv=(LIS_SCALAR *)lis_malloc((size_t)s->n*sizeof(LIS_SCALAR),
                                       "RAS interior diag");
    if(!s->gamma_of || !s->diag_inv) return LIS_OUT_OF_MEMORY;
    for(i=0;i<s->n;i++) { s->gamma_of[i]=-1; s->diag_inv[i]=0.0; }
    for(i=0;i<s->n;i++)
    for(k=A->ptr[i];k<A->ptr[i+1];k++)
    {
        j=A->index[k];
        if(s->owner[i]!=s->owner[j])
        {
            s->gamma_of[i]=0;
            s->gamma_of[j]=0;
        }
    }
    for(i=0;i<s->n;i++) if(s->gamma_of[i]>=0) ++ng;
    if(ng==0) return LIS_SUCCESS;
    if(ng>LIS_ILUT_RAS_SCHUR_LIMIT)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "ILUT RAS Schur interface exceeds 256 rows");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
    s->ngamma=ng;
    s->gamma=(LIS_INT *)lis_malloc((size_t)ng*sizeof(LIS_INT),"RAS gamma rows");
    s->pivots=(LIS_INT *)lis_malloc((size_t)ng*sizeof(LIS_INT),"RAS pivot perm");
    s->schur_rhs=(LIS_SCALAR *)lis_malloc((size_t)ng*sizeof(LIS_SCALAR),
                                          "RAS schur rhs");
    s->schur_delta=(LIS_SCALAR *)lis_malloc((size_t)ng*sizeof(LIS_SCALAR),
                                            "RAS schur solution");
    entries=(size_t)ng*(size_t)ng;
    s->schur_lu=(LIS_SCALAR *)lis_malloc(entries*sizeof(LIS_SCALAR),
                                        "RAS Schur LU");
    if(!s->gamma || !s->pivots || !s->schur_rhs ||
       !s->schur_delta || !s->schur_lu) return LIS_OUT_OF_MEMORY;
    memset(s->schur_lu,0,entries*sizeof(LIS_SCALAR));

    ng=0;
    for(i=0;i<s->n;i++)
    {
        if(s->gamma_of[i]>=0)
        {
            s->gamma_of[i]=ng;
            s->gamma[ng]=i;
            s->pivots[ng]=ng;
            ng++;
        }
        else
        {
            /*
             * diag(A_II)^-1 is meaningful only when every interior row
             * has a usable diagonal relative to its own coefficient
             * scale. Quietly dropping a zero/tiny KKT multiplier pivot
             * would create a misleading approximate Schur operator.
             */
            LIS_REAL row_max=0.0,diag_abs=0.0;
            LIS_SCALAR diag=0.0;
            for(k=A->ptr[i];k<A->ptr[i+1];k++)
            {
                LIS_REAL mag=fabs(A->value[k]);
                if(!isfinite((double)mag))
                {
                    LIS_SETERR(LIS_ERR_ILL_ARG,
                        "ILUT RAS Schur interior contains nonfinite entry");
                    return LIS_ERR_ILL_ARG;
                }
                if(mag>row_max) row_max=mag;
                if(A->index[k]==i) { diag=A->value[k]; diag_abs=mag; }
            }
            if(!(row_max>0.0) || !(diag_abs>1.0e-10*row_max))
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUT RAS Schur interior has missing/weak diagonal");
                return LIS_ERR_ILL_ARG;
            }
            s->diag_inv[i]=1.0/diag;
            if(!isfinite((double)s->diag_inv[i]))
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUT RAS Schur interior reciprocal overflows");
                return LIS_ERR_ILL_ARG;
            }
        }
    }

    for(g=0;g<s->ngamma;g++)
    {
        i=s->gamma[g];
        for(k=A->ptr[i];k<A->ptr[i+1];k++)
        {
            j=A->index[k];
            if(s->gamma_of[j]>=0)
                s->schur_lu[(size_t)g*s->ngamma+s->gamma_of[j]]+=A->value[k];
            else if(s->diag_inv[j]!=0.0)
                for(q=A->ptr[j];q<A->ptr[j+1];q++)
                {
                    t=A->index[q];
                    if(s->gamma_of[t]>=0)
                        s->schur_lu[(size_t)g*s->ngamma+s->gamma_of[t]]
                          -=A->value[k]*s->diag_inv[j]*A->value[q];
                }
        }
    }
    /*
     * Reject nonfinite or dangerously small pivots relative to the
     * overall approximate Schur scale, not a hard-coded absolute unit.
     * This is a safety gate, NOT a convergence/conditioning guarantee.
     */
    for(i=0;i<s->ngamma*s->ngamma;i++)
    {
        val=fabs(s->schur_lu[i]);
        if(!isfinite((double)val))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "ILUT RAS approximate Schur has nonfinite coefficients");
            return LIS_ERR_ILL_ARG;
        }
        if(val>schur_norm) schur_norm=val;
    }
    if(!(schur_norm>0.0))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "ILUT RAS approximate Schur is zero");
        return LIS_ERR_ILL_ARG;
    }

    /* Partial-pivot LU: preserve row-permutation for each RHS. */
    for(k=0;k<s->ngamma;k++)
    {
        best=k;peak=fabs(s->schur_lu[(size_t)k*s->ngamma+k]);
        for(i=k+1;i<s->ngamma;i++)
        {
            val=fabs(s->schur_lu[(size_t)i*s->ngamma+k]);
            if(val>peak) { best=i;peak=val; }
        }
        if(!isfinite((double)peak) || peak<=1.0e-12*schur_norm)
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "ILUT RAS approximate Schur complement has singular pivot");
            return LIS_ERR_ILL_ARG;
        }
        if(best!=k)
        {
            for(j=0;j<s->ngamma;j++)
            {
                tmp=s->schur_lu[(size_t)k*s->ngamma+j];
                s->schur_lu[(size_t)k*s->ngamma+j]=
                    s->schur_lu[(size_t)best*s->ngamma+j];
                s->schur_lu[(size_t)best*s->ngamma+j]=tmp;
            }
            p=s->pivots[k]; s->pivots[k]=s->pivots[best]; s->pivots[best]=p;
        }
        for(i=k+1;i<s->ngamma;i++)
        {
            s->schur_lu[(size_t)i*s->ngamma+k] /=
                s->schur_lu[(size_t)k*s->ngamma+k];
            for(j=k+1;j<s->ngamma;j++)
                s->schur_lu[(size_t)i*s->ngamma+j] -=
                    s->schur_lu[(size_t)i*s->ngamma+k]*
                    s->schur_lu[(size_t)k*s->ngamma+j];
        }
    }

    return lis_vector_duplicate(A,&s->residual);
}

LIS_INT
lis_precon_create_ilut_ras(LIS_SOLVER solver, LIS_PRECON precon)
{
    struct LIS_ILUT_RAS_STRUCT *s=NULL;
    LIS_MATRIX A=solver->A;
    LIS_INT n=A->n,d,i,k,err;
    if(A->nprocs!=1 || A->n!=A->gn || A->is!=0 ||
       A->matrix_type!=LIS_MATRIX_CSR || A->ptr==NULL ||
       A->index==NULL || A->value==NULL ||
       solver->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "ILUT RAS currently supports only serial double CSR");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
    if(solver->ilut_ras_overlap<0 || solver->ilut_ras_overlap>2 ||
       solver->ilut_ras_domains<2 || solver->ilut_ras_domains>16 ||
       solver->ilut_ras_domains>n)
        return LIS_ERR_ILL_ARG;
    if(solver->ilutp_full && !solver->ilutp_enabled)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "ILUT RAS: -ilutp_full requires -ilutp 1");
        return LIS_ERR_ILL_ARG;
    }
    if(solver->ilutp_enabled &&
       solver->precon_half_mode!=LIS_PRECON_HALF_OFF)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "ILUTP local factors do not support half storage");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
    s=(struct LIS_ILUT_RAS_STRUCT *)lis_malloc(sizeof(*s),"ILUT RAS state");
    if(!s) return LIS_OUT_OF_MEMORY;
    memset(s,0,sizeof(*s));
    precon->ilut_ras=s; /* generic preconditioner destruction owns cleanup */
    s->A=A; s->n=n; s->ndomains=solver->ilut_ras_domains;
    s->overlap=solver->ilut_ras_overlap;
    s->owner=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"RAS owners");
    s->parts=(struct LIS_ILUT_RAS_PART *)lis_malloc(
        (size_t)s->ndomains*sizeof(*s->parts),"RAS subdomains");
    /*
     * If one of these allocations fails, the generic preconditioner
     * destructor still visits any allocated parts. Initialize them first.
     */
    if(s->parts)
        memset(s->parts,0,(size_t)s->ndomains*sizeof(*s->parts));
    if(!s->owner || !s->parts) return LIS_OUT_OF_MEMORY;
    /*
     * Use identical integer partition boundaries for owner mapping and
     * extracted core domains. floor(i*ndomains/n) is NOT equivalent
     * to the partition intervals [floor(d*n/ndomains), ...); for
     * example n=9, ndomains=2 mislabels row 4 as owner 0.
     * Quotient/remainder arithmetic also avoids double rounding and
     * multiplication overflow for very large LIS_INT.
     */
    {
        LIS_INT block=n/s->ndomains,rem=n%s->ndomains;
        for(d=0;d<s->ndomains;d++)
        {
            struct LIS_ILUT_RAS_PART *p=&s->parts[d];
            p->start=d*block+(d*rem)/s->ndomains;
            p->end=(d+1)*block+((d+1)*rem)/s->ndomains;
            for(i=p->start;i<p->end;i++) s->owner[i]=d;
            err=lis_ilut_ras_make_subdomain(solver,s,p);
            if(err) return err;
        }
    }
    if(solver->ilut_schur)
    {
        err=lis_ilut_ras_make_schur(s);
        if(err) return err;
    }
    return LIS_SUCCESS;
}

LIS_INT
lis_psolve_ilut_ras(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
    struct LIS_ILUT_RAS_STRUCT *s=solver->precon->ilut_ras;
    LIS_INT d,i,k,gi,err,ng;
    LIS_SCALAR val;
    if(!s || !B || !X || B->n!=s->n || X->n!=s->n)
        return LIS_ERR_ILL_ARG;
    err=lis_vector_set_all((LIS_SCALAR)0.0,X);
    if(err) return err;
    for(d=0;d<s->ndomains;d++)
    {
        struct LIS_ILUT_RAS_PART *p=&s->parts[d];
        for(i=0;i<p->nlocal;i++)
            p->b->value[i]=B->value[p->global[i]];
        err=lis_psolve(p->solver,p->b,p->x);
        if(err) return err;
        /* RAS: inject each solve only into its nonoverlapping owner core. */
        for(i=0;i<p->nlocal;i++)
        {
            gi=p->global[i];
            if(gi>=p->start && gi<p->end)
                X->value[gi]=p->x->value[i];
        }
    }
    if(s->ngamma==0) return LIS_SUCCESS;
    ng=s->ngamma;
    err=lis_matvec(s->A,X,s->residual);
    if(err) return err;
    for(i=0;i<ng;i++)
        s->schur_rhs[i] =
            B->value[s->gamma[s->pivots[i]]] -
            s->residual->value[s->gamma[s->pivots[i]]];
    /* Solve factored approximate Schur via pivoted forward/backward solve. */
    for(i=0;i<ng;i++)
    {
        val=s->schur_rhs[i];
        for(k=0;k<i;k++) val-=s->schur_lu[(size_t)i*ng+k]*
                              s->schur_rhs[k];
        s->schur_rhs[i]=val;
    }
    for(i=ng-1;i>=0;i--)
    {
        val=s->schur_rhs[i];
        for(k=i+1;k<ng;k++) val-=s->schur_lu[(size_t)i*ng+k]*
                              s->schur_delta[k];
        s->schur_delta[i]=val/s->schur_lu[(size_t)i*ng+i];
        if(!isfinite((double)s->schur_delta[i])) return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<ng;i++) X->value[s->gamma[i]]+=s->schur_delta[i];
    /* Diagonal-approximate interior back substitution. */
    for(i=0;i<s->n;i++)
    {
        if(s->gamma_of[i]>=0 || s->diag_inv[i]==0.0) continue;
        val=0.0;
        for(k=s->A->ptr[i];k<s->A->ptr[i+1];k++)
        {
            gi=s->gamma_of[s->A->index[k]];
            if(gi>=0) val+=s->A->value[k]*s->schur_delta[gi];
        }
        X->value[i]-=s->diag_inv[i]*val;
    }
    return LIS_SUCCESS;
}

#endif /* real-valued ILUT RAS implementation */
