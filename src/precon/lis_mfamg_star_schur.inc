/*
 * Experimental MFAMG serial arrowhead terminal solver.
 *
 * This file is included by lis_precon_mfamg.c after the existing dense
 * factorization/solve helpers.  Do not enable it by default.
 *
 * The initial greedy aggregate identifies a small candidate core. The
 * induced outside graph may contain isolated 1x1 and 2x2 blocks. If
 * larger components remain, degree-aware core promotion falls back to
 * the diagonal-leaf strategy.  We then factor S=H-F*D^-1*G with a
 * block-diagonal D, subject to relative pivot and residual checks.
 */

struct LIS_MFAMG_STAR_SOLVER
{
    LIS_INT n, c, l;
    LIS_MATRIX A; /* borrowed terminal matrix, owned by the hierarchy */
    LIS_INT *core_map, *leaf_map, *core_rows, *leaf_rows, *pivots;
    LIS_INT *partner; /* leaf-index partner, or -1 for a scalar leaf */
    LIS_REAL *row_scale;
    LIS_SCALAR *d_inv, *d_cross, *leaf_work;
    LIS_SCALAR *f, *g, *schur, *lu, *rhs, *sol;
};

static void
lis_mfamg_star_destroy(struct LIS_MFAMG_STAR_SOLVER *s)
{
    if(s==NULL) return;
    if(s->core_map) lis_free(s->core_map);
    if(s->leaf_map) lis_free(s->leaf_map);
    if(s->core_rows) lis_free(s->core_rows);
    if(s->leaf_rows) lis_free(s->leaf_rows);
    if(s->pivots) lis_free(s->pivots);
    if(s->partner) lis_free(s->partner);
    if(s->row_scale) lis_free(s->row_scale);
    if(s->d_inv) lis_free(s->d_inv);
    if(s->d_cross) lis_free(s->d_cross);
    if(s->leaf_work) lis_free(s->leaf_work);
    if(s->f) lis_free(s->f);
    if(s->g) lis_free(s->g);
    if(s->schur) lis_free(s->schur);
    if(s->lu) lis_free(s->lu);
    if(s->rhs) lis_free(s->rhs);
    if(s->sol) lis_free(s->sol);
    lis_free(s);
}

static LIS_INT
lis_mfamg_star_prepare(
    LIS_MATRIX A,
    const LIS_INT *ag,
    LIS_INT nag,
    LIS_INT core_limit,
    struct LIS_MFAMG_STAR_SOLVER **out)
{
#ifndef _COMPLEX
    struct LIS_MFAMG_STAR_SOLVER *s=NULL;
    LIS_INT *sizes=NULL,*core=NULL,*degree=NULL;
    LIS_INT i,j,p,a,largest,largest_size,c,l,k,m,err=LIS_SUCCESS;
    LIS_INT is_matching=LIS_TRUE;
    LIS_REAL mag,row_max;
    size_t cl,cc;

    if(out==NULL) return LIS_ERR_ILL_ARG;
    *out=NULL;
    if(A==NULL || ag==NULL || nag<=0 || core_limit<2 ||
       A->nprocs!=1 || A->matrix_type!=LIS_MATRIX_CSR ||
       A->n!=A->gn || A->is!=0 || A->n<=0 || A->n>8192 ||
       A->ptr==NULL || A->index==NULL || A->value==NULL)
        return LIS_SUCCESS;

    sizes=(LIS_INT *)lis_malloc((size_t)nag*sizeof(LIS_INT),
                                "mfamg_star::sizes");
    core=(LIS_INT *)lis_malloc((size_t)A->n*sizeof(LIS_INT),
                               "mfamg_star::core");
    degree=(LIS_INT *)lis_malloc((size_t)A->n*sizeof(LIS_INT),
                                 "mfamg_star::degree");
    if(sizes==NULL || core==NULL || degree==NULL)
    {
        err=LIS_OUT_OF_MEMORY;
        goto done;
    }
    memset(sizes,0,(size_t)nag*sizeof(LIS_INT));
    memset(core,0,(size_t)A->n*sizeof(LIS_INT));
    for(i=0;i<A->n;i++)
    {
        if(ag[i]<0 || ag[i]>=nag) { err=LIS_ERR_ILL_ARG; goto done; }
        sizes[ag[i]]++;
    }
    largest=0; largest_size=sizes[0];
    for(i=1;i<nag;i++)
        if(sizes[i]>largest_size)
        {
            largest=i;
            largest_size=sizes[i];
        }
    if(largest_size<2 || largest_size>core_limit) goto done;
    for(i=0;i<A->n;i++) if(ag[i]==largest) core[i]=1;

    /*
     * A matching (undirected degree <=1) in the induced outside graph
     * means D consists exclusively of 1x1 and 2x2 blocks. Recognize
     * both directed orientations; missing reverse entries are zero.
     * Use degree[] here as a temporary original-row partner array.
     */
    for(i=0;i<A->n;i++) degree[i]=-1;
    for(i=0;i<A->n && is_matching;i++)
    {
        if(core[i]) continue;
        for(p=A->ptr[i];p<A->ptr[i+1];p++)
        {
            j=A->index[p];
            if(j<0 || j>=A->n)
            {
                err=LIS_ERR_ILL_ARG;
                goto done;
            }
            if(j==i || core[j] || A->value[p]==(LIS_SCALAR)0)
                continue;
            if((degree[i]>=0 && degree[i]!=j) ||
               (degree[j]>=0 && degree[j]!=i))
            {
                is_matching=LIS_FALSE;
                break;
            }
            degree[i]=j;
            degree[j]=i;
        }
    }

    c=largest_size;
    if(!is_matching)
    {
        /*
         * Unsupported components retain the previous degree-aware
         * vertex-cover fallback. It may promote a high-degree hub,
         * but must never exceed the core limit.
         */
        while(c<core_limit)
        {
            LIS_INT best=-1;
            LIS_INT best_degree=0;
            memset(degree,0,(size_t)A->n*sizeof(LIS_INT));
            for(i=0;i<A->n;i++)
            {
                if(core[i]) continue;
                for(p=A->ptr[i];p<A->ptr[i+1];p++)
                {
                    j=A->index[p];
                    if(j<0 || j>=A->n)
                    {
                        err=LIS_ERR_ILL_ARG;
                        goto done;
                    }
                    if(j!=i && A->value[p]!=(LIS_SCALAR)0 && !core[j])
                    {
                        degree[i]++;
                        degree[j]++;
                    }
                }
            }
            for(i=0;i<A->n;i++)
            {
                if(!core[i] && degree[i]>best_degree)
                {
                    best=i;
                    best_degree=degree[i];
                }
            }
            if(best<0) break;
            core[best]=1;
            c++;
        }
        /* The fallback must leave a strictly diagonal outside block. */
        for(i=0;i<A->n;i++)
        {
            if(core[i]) continue;
            for(p=A->ptr[i];p<A->ptr[i+1];p++)
            {
                j=A->index[p];
                if(j<0 || j>=A->n)
                {
                    err=LIS_ERR_ILL_ARG;
                    goto done;
                }
                if(j!=i && A->value[p]!=(LIS_SCALAR)0 && !core[j])
                    goto done;
            }
        }
    }

    l=A->n-c;
    if(c>core_limit || c<=0 || l<=c || l<=0) goto done;
    if((size_t)c>((size_t)-1)/(size_t)l ||
       (size_t)c>((size_t)-1)/(size_t)c)
    { err=LIS_OUT_OF_MEMORY; goto done; }
    cl=(size_t)c*(size_t)l;
    cc=(size_t)c*(size_t)c;
    if(cl>((size_t)-1)/sizeof(LIS_SCALAR) ||
       cc>((size_t)-1)/sizeof(LIS_SCALAR))
    { err=LIS_OUT_OF_MEMORY; goto done; }

    /*
     * Limit dense F/G plus S/LU workspace to 64 MiB.  Exceeding the
     * budget is not an allocation error: retain the normal hierarchy.
     */
    if(cl+cc>((size_t)(64U*1024U*1024U))/
              (2*sizeof(LIS_SCALAR)))
        goto done;

    s=(struct LIS_MFAMG_STAR_SOLVER *)lis_malloc(sizeof(*s),
                                                 "mfamg_star::state");
    if(s==NULL) { err=LIS_OUT_OF_MEMORY; goto done; }
    memset(s,0,sizeof(*s));
    s->n=A->n; s->c=c; s->l=l; s->A=A;
#define STAR_ALLOC(field,count,type) do { \
    s->field=(type *)lis_malloc((size_t)(count)*sizeof(type), \
                                "mfamg_star::" #field); \
    if(s->field==NULL) { err=LIS_OUT_OF_MEMORY; goto done; } \
} while(0)
    STAR_ALLOC(core_map,A->n,LIS_INT);
    STAR_ALLOC(leaf_map,A->n,LIS_INT);
    STAR_ALLOC(core_rows,c,LIS_INT);
    STAR_ALLOC(leaf_rows,l,LIS_INT);
    STAR_ALLOC(pivots,c,LIS_INT);
    STAR_ALLOC(partner,l,LIS_INT);
    STAR_ALLOC(row_scale,l,LIS_REAL);
    STAR_ALLOC(d_inv,l,LIS_SCALAR);
    STAR_ALLOC(d_cross,l,LIS_SCALAR);
    STAR_ALLOC(leaf_work,l,LIS_SCALAR);
    STAR_ALLOC(f,cl,LIS_SCALAR);
    STAR_ALLOC(g,cl,LIS_SCALAR);
    STAR_ALLOC(schur,cc,LIS_SCALAR);
    STAR_ALLOC(lu,cc,LIS_SCALAR);
    STAR_ALLOC(rhs,c,LIS_SCALAR);
    STAR_ALLOC(sol,c,LIS_SCALAR);
#undef STAR_ALLOC
    memset(s->f,0,cl*sizeof(LIS_SCALAR));
    memset(s->g,0,cl*sizeof(LIS_SCALAR));
    memset(s->schur,0,cc*sizeof(LIS_SCALAR));
    memset(s->d_inv,0,(size_t)l*sizeof(LIS_SCALAR));
    memset(s->d_cross,0,(size_t)l*sizeof(LIS_SCALAR));
    for(k=0;k<l;k++) s->partner[k]=-1;

    c=0; l=0;
    for(i=0;i<A->n;i++)
    {
        s->core_map[i]=-1; s->leaf_map[i]=-1;
        if(core[i])
        {
            s->core_map[i]=c;
            s->core_rows[c++]=i;
        }
        else
        {
            s->leaf_map[i]=l;
            s->leaf_rows[l++]=i;
        }
    }
    c=s->c; l=s->l;

    /*
     * Reconstruct the matching in compact leaf indices. Every nonzero
     * off-diagonal outside edge must agree with the same reciprocal
     * partner; directed and nonsymmetric matrices are supported.
     */
    for(i=0;i<A->n;i++)
    {
        if(s->core_map[i]>=0) continue;
        k=s->leaf_map[i];
        for(p=A->ptr[i];p<A->ptr[i+1];p++)
        {
            j=A->index[p];
            if(j<0 || j>=A->n) { err=LIS_ERR_ILL_ARG; goto done; }
            if(j==i || s->leaf_map[j]<0 ||
               A->value[p]==(LIS_SCALAR)0)
                continue;
            m=s->leaf_map[j];
            if((s->partner[k]>=0 && s->partner[k]!=m) ||
               (s->partner[m]>=0 && s->partner[m]!=k))
                goto done;
            s->partner[k]=m;
            s->partner[m]=k;
        }
    }

    for(i=0;i<A->n;i++)
    {
        a=s->core_map[i];
        row_max=0;
        for(p=A->ptr[i];p<A->ptr[i+1];p++)
        {
            j=A->index[p];
            if(j<0 || j>=A->n) { err=LIS_ERR_ILL_ARG; goto done; }
            mag=fabs((double)A->value[p]);
            if(!isfinite((double)mag)) goto done;
            if(mag>row_max) row_max=mag;
            if(a>=0)
            {
                if(s->core_map[j]>=0)
                    s->schur[(size_t)a*c+s->core_map[j]]+=A->value[p];
                else
                    s->f[(size_t)a*l+s->leaf_map[j]]+=A->value[p];
            }
            else if(s->core_map[j]>=0)
                s->g[(size_t)s->leaf_map[i]*c+s->core_map[j]]+=A->value[p];
            else if(i==j)
                s->d_inv[s->leaf_map[i]]+=A->value[p];
            else if(s->partner[s->leaf_map[i]]==s->leaf_map[j])
                s->d_cross[s->leaf_map[i]]+=A->value[p];
            else if(A->value[p]!=(LIS_SCALAR)0)
                goto done;
        }
        if(a<0)
            s->row_scale[s->leaf_map[i]]=row_max;
    }

    /* Invert each isolated scalar or nonsymmetric 2x2 block. */
    for(k=0;k<l;k++)
    {
        LIS_REAL scale,det,denom,aa,bb,cc0,dd;
        m=s->partner[k];
        if(m>=0 && m<k) continue;
        if(m<0)
        {
            mag=fabs((double)s->d_inv[k]);
            if(!isfinite((double)mag) || s->row_scale[k]==0 ||
               mag<=1.0e-12*s->row_scale[k])
                goto done;
            s->d_inv[k]=(LIS_SCALAR)1/s->d_inv[k];
            if(!isfinite((double)s->d_inv[k])) goto done;
            continue;
        }
        aa=(LIS_REAL)s->d_inv[k];
        dd=(LIS_REAL)s->d_inv[m];
        bb=(LIS_REAL)s->d_cross[k];
        cc0=(LIS_REAL)s->d_cross[m];
        scale=fabs((double)aa);
        if(fabs((double)bb)>scale) scale=fabs((double)bb);
        if(fabs((double)cc0)>scale) scale=fabs((double)cc0);
        if(fabs((double)dd)>scale) scale=fabs((double)dd);
        if(!isfinite((double)scale) || scale==0.0 ||
           scale<=1.0e-12*s->row_scale[k] ||
           scale<=1.0e-12*s->row_scale[m])
            goto done;
        aa/=scale; bb/=scale; cc0/=scale; dd/=scale;
        det=aa*dd-bb*cc0;
        if(!isfinite((double)det) || fabs((double)det)<=1.0e-12)
            goto done;
        denom=det*scale;
        s->d_inv[k]=(LIS_SCALAR)(dd/denom);
        s->d_inv[m]=(LIS_SCALAR)(aa/denom);
        s->d_cross[k]=(LIS_SCALAR)(-bb/denom);
        s->d_cross[m]=(LIS_SCALAR)(-cc0/denom);
        if(!isfinite((double)s->d_inv[k]) ||
           !isfinite((double)s->d_inv[m]) ||
           !isfinite((double)s->d_cross[k]) ||
           !isfinite((double)s->d_cross[m]))
            goto done;
    }

    for(k=0;k<l;k++)
        for(i=0;i<c;i++)
        {
            LIS_SCALAR t=s->f[(size_t)i*l+k];
            if(t==(LIS_SCALAR)0) continue;
            for(j=0;j<c;j++)
            {
                LIS_SCALAR v=s->d_inv[k]*s->g[(size_t)k*c+j];
                if(s->partner[k]>=0)
                    v+=s->d_cross[k]*
                       s->g[(size_t)s->partner[k]*c+j];
                if(v!=(LIS_SCALAR)0)
                    s->schur[(size_t)i*c+j]-=t*v;
            }
        }
    for(k=0;k<c*c;k++)
        if(!isfinite((double)s->schur[k])) goto done;
    memcpy(s->lu,s->schur,cc*sizeof(LIS_SCALAR));
    err=lis_mfamg_dense_factor(c,s->lu,s->pivots);
    if(err) goto done;
    /*
     * lis_mfamg_dense_factor() already performs scaled partial
     * pivoting and marks relatively small pivots by exact zero.
     * Do not impose an unrelated absolute 1e-14 cutoff: valid
     * matrices may have all entries far below that scale.
     */
    for(i=0;i<c;i++)
        if(!isfinite((double)s->lu[(size_t)i*c+i]) ||
           s->lu[(size_t)i*c+i]==(LIS_SCALAR)0)
            goto done;

    *out=s;
    s=NULL;

done:
    if(sizes) lis_free(sizes);
    if(core) lis_free(core);
    if(degree) lis_free(degree);
    lis_mfamg_star_destroy(s);
    return err;
#else
    (void)A; (void)ag; (void)nag; (void)core_limit;
    if(out==NULL) return LIS_ERR_ILL_ARG;
    *out=NULL;
    return LIS_SUCCESS;
#endif
}

static LIS_INT
lis_mfamg_star_solve(
    struct LIS_MFAMG_STAR_SOLVER *s,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    LIS_INT i,j,k,err;
    LIS_SCALAR t;
    if(s==NULL || b==NULL || x==NULL ||
       b->n!=s->n || x->n!=s->n) return LIS_ERR_ILL_ARG;
    for(k=0;k<s->l;k++)
    {
        t=s->d_inv[k]*b->value[s->leaf_rows[k]];
        if(s->partner[k]>=0)
            t+=s->d_cross[k]*
               b->value[s->leaf_rows[s->partner[k]]];
        s->leaf_work[k]=t;
    }
    for(i=0;i<s->c;i++)
    {
        t=b->value[s->core_rows[i]];
        for(k=0;k<s->l;k++)
            if(s->f[(size_t)i*s->l+k]!=(LIS_SCALAR)0)
                t-=s->f[(size_t)i*s->l+k]*s->leaf_work[k];
        s->rhs[i]=t;
    }
    err=lis_mfamg_dense_solve_values(
        s->c,s->schur,s->lu,s->pivots,s->rhs,s->sol,LIS_FALSE);
    if(err) return err;
    for(k=0;k<s->l;k++)
    {
        t=b->value[s->leaf_rows[k]];
        for(j=0;j<s->c;j++)
            if(s->g[(size_t)k*s->c+j]!=(LIS_SCALAR)0)
                t-=s->g[(size_t)k*s->c+j]*s->sol[j];
        s->leaf_work[k]=t;
    }
    for(k=0;k<s->l;k++)
    {
        t=s->d_inv[k]*s->leaf_work[k];
        if(s->partner[k]>=0)
            t+=s->d_cross[k]*s->leaf_work[s->partner[k]];
        x->value[s->leaf_rows[k]]=t;
    }
    for(i=0;i<s->c;i++)
        x->value[s->core_rows[i]]=s->sol[i];

#ifndef _COMPLEX
    /*
     * Check the original coarse operator, not only the reduced Schur
     * equations.  A poor factorization must never silently return an
     * inaccurate preconditioner correction.
     */
    if(s->A==NULL || s->A->n!=s->n ||
       s->A->ptr==NULL || s->A->index==NULL ||
       s->A->value==NULL)
        return LIS_ERR_ILL_ARG;

    for(i=0;i<s->n;i++)
    {
        LIS_SCALAR ax=(LIS_SCALAR)0;
        LIS_REAL scale=fabs((double)b->value[i]);
        LIS_REAL discrepancy;
        LIS_INT p;
        for(p=s->A->ptr[i];p<s->A->ptr[i+1];p++)
        {
            j=s->A->index[p];
            if(j<0 || j>=s->n) return LIS_ERR_ILL_ARG;
            ax+=s->A->value[p]*x->value[j];
            scale+=fabs((double)s->A->value[p])*
                   fabs((double)x->value[j]);
        }
        discrepancy=fabs((double)(b->value[i]-ax));
        if(!isfinite((double)scale) ||
           !isfinite((double)discrepancy) ||
           (scale==0.0 ? discrepancy!=0.0 :
                           discrepancy>1.0e-9*scale))
        {
            LIS_SETERR(LIS_BREAKDOWN,
                "MFAMG star-Schur correction failed backward-error check\n");
            return LIS_BREAKDOWN;
        }
    }
#endif

    return LIS_SUCCESS;
}

/* Internal regression hook. */
LIS_INT
lis_precon_mfamg_star_schur_for_test(
    LIS_MATRIX A,
    const LIS_INT *aggregate,
    LIS_INT nag,
    LIS_INT core_limit,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_INT *used,
    LIS_INT *core_size)
{
    struct LIS_MFAMG_STAR_SOLVER *s=NULL;
    LIS_INT err;
    if(used==NULL || core_size==NULL) return LIS_ERR_ILL_ARG;
    *used=LIS_FALSE;
    *core_size=0;
    err=lis_mfamg_star_prepare(A,aggregate,nag,core_limit,&s);
    if(err) return err;
    if(s!=NULL)
    {
        *used=LIS_TRUE;
        *core_size=s->c;
        err=lis_mfamg_star_solve(s,b,x);
    }
    lis_mfamg_star_destroy(s);
    return err;
}
