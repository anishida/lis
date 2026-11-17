/*
 * Opt-in serial real-double ILUTP with dynamic COLUMN pivoting.
 *
 * The active permutation maps factor-space columns to source columns.
 * U stores source-column labels while factors are constructed.  This
 * avoids O(nnz(U)) rewrites for each pivot.  Upon completion U indices
 * are remapped once to final factor-space positions.
 *
 * L U approximates A Q, where Q maps factor-space coordinates to
 * original columns (precon->ilutp_colperm).
 *
 * Keep this implementation separate from the existing ILUT kernels;
 * no legacy behavior is changed when -ilutp 0.
 */

typedef struct
{
    LIS_INT col;
    LIS_SCALAR value;
    LIS_REAL magnitude;
} LIS_ILUTP_ENTRY;

static int
lis_ilutp_compare_desc(const void *va,const void *vb)
{
    const LIS_ILUTP_ENTRY *a=(const LIS_ILUTP_ENTRY *)va;
    const LIS_ILUTP_ENTRY *b=(const LIS_ILUTP_ENTRY *)vb;
    if(a->magnitude>b->magnitude) return -1;
    if(a->magnitude<b->magnitude) return 1;
    if(a->col<b->col) return -1;
    if(a->col>b->col) return 1;
    return 0;
}

/* A min-heap avoids O(lower_fill^2) selection scans. Lower column fill
 * may appear during elimination, so a fixed initial sort is insufficient. */
static void lis_ilutp_heap_push(LIS_INT *a,LIS_INT *count,LIS_INT col)
{
    size_t i=(size_t)(*count);
    (*count)++;
    while(i>0 && a[(i-1)/2]>col)
    {
        a[i]=a[(i-1)/2];
        i=(i-1)/2;
    }
    a[i]=col;
}
static LIS_INT lis_ilutp_heap_pop(LIS_INT *a,LIS_INT *count)
{
    LIS_INT first=a[0],last;
    size_t i=0,n,c;
    (*count)--;
    if(*count==0) return first;
    last=a[*count];
    n=(size_t)(*count);
    while(2*i+1<n)
    {
        c=2*i+1;
        if(c+1<n && a[c+1]<a[c]) c++;
        if(last<=a[c]) break;
        a[i]=a[c];
        i=c;
    }
    a[i]=last;
    return first;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_ilutp_csr"
static LIS_INT
lis_precon_create_ilutp_csr(LIS_SOLVER solver,LIS_PRECON precon)
{
#if defined(_COMPLEX) || defined(USE_COMPLEX)
    (void)solver; (void)precon;
    return LIS_ERR_NOT_IMPLEMENTED;
#else
    LIS_MATRIX A=solver->A;
    LIS_MATRIX_ILU L=NULL,U=NULL;
    LIS_VECTOR D=NULL;
    LIS_INT n=A->n,i,j,k,orig,stamp,nlo,nup,processed;
    LIS_INT best,t,col,cap,nnz,m,err=LIS_SUCCESS,kept;
    LIS_INT *perm=NULL,*inv=NULL,*seen=NULL,*lower=NULL,*finished=NULL,*upper=NULL;
    LIS_SCALAR *w=NULL,fact,delta,pivot,swapval;
    LIS_ILUTP_ENTRY *entries=NULL;
    LIS_REAL drop,rate,tau,pivot_tol,rowmax,rowsum,tolnorm,maxval,v;
    LIS_INT swapcol,swaps=0;
    LIS_INT full=solver->ilutp_full;
    uint64_t fill_used=0,fill_limit=0;
#ifdef _OPENMP
    /* Existing ILU triangular kernels partition factors by OpenMP rank. */
    if(omp_get_max_threads()!=1)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "ILUTP requires one solver OpenMP thread");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
#endif
    if(A->matrix_type!=LIS_MATRIX_CSR || A->nprocs!=1 ||
       A->n!=A->gn || A->is!=0 || n<1 ||
       A->ptr==NULL || A->index==NULL || A->value==NULL ||
       solver->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE ||
       solver->precon_half_mode!=LIS_PRECON_HALF_OFF)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "ILUTP supports serial real-double CSR with FP16 disabled");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
    if((size_t)n>SIZE_MAX/sizeof(LIS_ILUTP_ENTRY) ||
       (size_t)n>SIZE_MAX/sizeof(LIS_SCALAR) ||
       (size_t)n>SIZE_MAX/sizeof(LIS_INT))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"ILUTP temporary allocation size overflow");
        return LIS_ERR_ILL_ARG;
    }
    if(A->ptr[0]!=0 || A->ptr[n]!=A->nnz)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"ILUTP CSR row pointers are inconsistent");
        return LIS_ERR_ILL_ARG;
    }

    drop=(LIS_REAL)solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
    rate=(LIS_REAL)solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
    pivot_tol=(LIS_REAL)solver->params[
        LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN];
    tau=solver->ilutp_threshold;
    if(!isfinite((double)drop) || drop<0.0 ||
       !isfinite((double)rate) || rate<0.0 ||
       !isfinite((double)pivot_tol) || pivot_tol<0.0 ||
       !isfinite((double)tau) || tau<0.0 || tau>1.0)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "ILUTP requires finite nonnegative drop/rate/pivot tolerance and threshold in [0,1]");
        return LIS_ERR_ILL_ARG;
    }
    {
        double budget=((double)A->nnz/(2.0*(double)n))*rate;
        cap=(budget>=(double)(n-1)) ? n-1 : (LIS_INT)budget;
    }
    if(full)
    {
        /* No dropping or diagonal perturbation. Pivot against the
         * surviving Schur-complement column when the requested
         * -ilutp_threshold deems the current pivot insufficient.
         * A*Q=L*U in exact arithmetic; Q is carried by LIS ILUTP.
         * Factor nnz is bounded independently of sparse pivot fill. */
        uint64_t imax=sizeof(LIS_INT)==4 ? INT32_MAX : INT64_MAX;
        uint64_t auto_limit=(uint64_t)A->nnz;
        if(auto_limit>imax/20) auto_limit=imax;
        else auto_limit*=20;
        if(auto_limit<(uint64_t)n*4) auto_limit=(uint64_t)n*4;
        fill_limit=solver->ilutp_full_limit ?
            (uint64_t)solver->ilutp_full_limit : auto_limit;
        drop=0.0;
        pivot_tol=0.0;
        /* Respect -ilutp_threshold also in complete-fill mode.
         * A strongest-column pivot at every row can introduce
         * unnecessary permutations and sparse factor fill. */
        cap=n-1;
    }

    /* Never overwrite live factor/permutation pointers on a reused object. */
    if(precon->L || precon->U || precon->D ||
       precon->ilutp_colperm || precon->ilut_perm)
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"ILUTP requires a fresh preconditioner");
        return LIS_ERR_ILL_ARG;
    }
    /* The no-pivot case is precisely legacy ILUT, including drop/fill
     * selection and floating-point evaluation order. */
    if(tau==0.0 && !full)
    {
        err=lis_precon_create_ilut_csr(solver,precon);
        if(err) return err;
        perm=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),
                                    "ILUTP identity permutation");
        if(!perm)
        {
            /* Legacy ILUT succeeded but the identity-permutation
             * allocation failed. Do not expose incomplete setup state. */
            lis_matrix_ilu_destroy(precon->L);
            lis_matrix_ilu_destroy(precon->U);
            lis_vector_destroy(precon->D);
            precon->L=NULL;
            precon->U=NULL;
            precon->D=NULL;
            LIS_SETERR_MEM((size_t)n*sizeof(LIS_INT));
            return LIS_OUT_OF_MEMORY;
        }
        for(i=0;i<n;i++) perm[i]=i;
        precon->ilutp_colperm=perm;
        return LIS_SUCCESS;
    }
    perm=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),
                               "ILUTP column permutation");
    inv=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),
                              "ILUTP inverse permutation");
    seen=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),
                               "ILUTP sparse row markers");
    lower=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),
                                "ILUTP lower min-heap");
    finished=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),
                                   "ILUTP processed lower entries");
    upper=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),
                                "ILUTP upper queue");
    w=(LIS_SCALAR *)lis_malloc((size_t)n*sizeof(LIS_SCALAR),
                               "ILUTP row coefficients");
    entries=(LIS_ILUTP_ENTRY *)lis_malloc(
        (size_t)n*sizeof(LIS_ILUTP_ENTRY),"ILUTP fill ranking");
    if(!perm || !inv || !seen || !lower || !finished || !upper ||
       !w || !entries)
    {
        err=LIS_OUT_OF_MEMORY;
        goto cleanup;
    }
    memset(seen,0,(size_t)n*sizeof(LIS_INT));
    for(i=0;i<n;i++) perm[i]=inv[i]=i;

    err=lis_matrix_ilu_create(n,1,&L);
    if(err) goto cleanup;
    err=lis_matrix_ilu_create(n,1,&U);
    if(err) goto cleanup;
    err=lis_matrix_ilu_setCR(L);
    if(err) goto cleanup;
    err=lis_matrix_ilu_setCR(U);
    if(err) goto cleanup;
    err=lis_vector_duplicate(A,&D);
    if(err) goto cleanup;

    for(i=0;i<n;i++)
    {
        stamp=i+1;
        nlo=nup=processed=0;
        rowmax=rowsum=0.0;
        nnz=A->ptr[i+1]-A->ptr[i];
        seen[i]=stamp;
        w[i]=(LIS_SCALAR)0.0;

        for(k=A->ptr[i];k<A->ptr[i+1];k++)
        {
            orig=A->index[k];
            if(orig<0 || orig>=n ||
               !isfinite((double)A->value[k]))
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUTP CSR contains invalid column indices or nonfinite values");
                err=LIS_ERR_ILL_ARG;
                goto cleanup;
            }
            col=inv[orig];
            if(seen[col]!=stamp)
            {
                seen[col]=stamp;
                w[col]=(LIS_SCALAR)0.0;
                if(col<i) lis_ilutp_heap_push(lower,&nlo,col);
                else if(col>i) upper[nup++]=col;
            }
            w[col]+=A->value[k];
            v=fabs(A->value[k]);
            rowsum+=v;
            if(v>rowmax) rowmax=v;
        }
        tolnorm=nnz>0 ? drop*rowsum/(LIS_REAL)nnz : 0.0;

        /* Lower fill enters the min-heap as it is discovered; pop the
         * smallest remaining column each time. */
        while(nlo>0)
        {
            j=lis_ilutp_heap_pop(lower,&nlo);
            /* Explicit CSR zeros and exact cancellation must not
             * propagate zero multipliers or introduce fake fill. */
            if(w[j]==0.0) continue;
            finished[processed++]=j;
            fact=w[j]*D->value[j];
            if(!isfinite((double)fact))
            {
                err=LIS_BREAKDOWN;
                goto cleanup;
            }
            w[j]=fact;

            for(k=0;k<U->nnz[j];k++)
            {
                /* Earlier U entries retain source-column labels until
                   all column pivot swaps are finished. */
                col=inv[U->index[j][k]];
                delta=-fact*U->value[j][k];
                if(!isfinite((double)delta))
                {
                    err=LIS_BREAKDOWN;
                    goto cleanup;
                }
                if(delta==0.0) continue;
                if(seen[col]!=stamp)
                {
                    if(fabs(delta)<tolnorm) continue;
                    seen[col]=stamp;
                    w[col]=(LIS_SCALAR)0.0;
                    if(col<i) lis_ilutp_heap_push(lower,&nlo,col);
                    else if(col>i) upper[nup++]=col;
                }
                w[col]+=delta;
            }
        }

        /* Select the largest surviving upper entry as column pivot when
           the current diagonal fails the relative threshold. */
        best=i;
        maxval=fabs(w[i]);
        for(t=0;t<nup;t++)
        {
            col=upper[t];
            v=fabs(w[col]);
            if(v>maxval)
            {
                maxval=v;
                best=col;
            }
        }
        if(best!=i && fabs(w[i])<tau*maxval)
        {
            swapcol=perm[i];
            perm[i]=perm[best];
            perm[best]=swapcol;
            inv[perm[i]]=i;
            inv[perm[best]]=best;
            swapval=w[i];
            w[i]=w[best];
            w[best]=swapval;
            ++swaps;
        }

        pivot=w[i];
        if(pivot_tol>0.0)
            pivot=lis_precon_regularize_ilu_pivot(
                pivot,rowmax,pivot_tol);
        if(!isfinite((double)pivot) || pivot==0.0 ||
           !isfinite(1.0/(double)pivot))
        {
            err=LIS_BREAKDOWN;
            goto cleanup;
        }
        D->value[i]=(LIS_SCALAR)(1.0/pivot);

        /* Guard actual L+U storage before allocating the current row.
         * The full mode fails cleanly rather than silently dropping fill
         * or exceeding its explicit/automatic memory budget. */
        if(full)
        {
            uint64_t next=fill_used;
            if(!isfinite((double)w[i]))
            {
                err=LIS_BREAKDOWN;
                goto cleanup;
            }
            for(t=0;t<processed;t++)
            {
                if(!isfinite((double)w[finished[t]]))
                { err=LIS_BREAKDOWN; goto cleanup; }
                if(w[finished[t]]!=0.0) next++;
            }
            for(t=0;t<nup;t++)
            {
                if(!isfinite((double)w[upper[t]]))
                { err=LIS_BREAKDOWN; goto cleanup; }
                if(w[upper[t]]!=0.0) next++;
            }
            if(next>fill_limit)
            {
                LIS_SETERR(LIS_OUT_OF_MEMORY,
                    "ILUTP full factor exceeds L+U nonzero safety limit");
                err=LIS_OUT_OF_MEMORY;
                goto cleanup;
            }
            fill_used=next;
        }
        /* Keep the largest lower multipliers (existing ILUT fill budget). */
        m=0;
        for(t=0;t<processed;t++)
        {
            col=finished[t];
            if(w[col]==0.0) continue;
            entries[m].col=col;
            entries[m].value=w[col];
            entries[m].magnitude=fabs(w[col]);
            m++;
        }
        if(!full)
            qsort(entries,(size_t)m,sizeof(entries[0]),lis_ilutp_compare_desc);
        kept=m<cap?m:cap;
        L->nnz[i]=kept;
        if(kept)
        {
            L->index[i]=(LIS_INT *)malloc((size_t)kept*sizeof(LIS_INT));
            L->value[i]=(LIS_SCALAR *)malloc((size_t)kept*sizeof(LIS_SCALAR));
            if(!L->index[i] || !L->value[i])
            {
                err=LIS_OUT_OF_MEMORY;
                goto cleanup;
            }
            for(t=0;t<kept;t++)
            {
                L->index[i][t]=entries[t].col;
                L->value[i][t]=entries[t].value;
            }
        }

        /* Store original upper column labels for future pivot changes. */
        m=0;
        for(t=0;t<nup;t++)
        {
            col=upper[t];
            if(w[col]==0.0) continue;
            entries[m].col=perm[col];
            entries[m].value=w[col];
            entries[m].magnitude=fabs(w[col]);
            m++;
        }
        if(!full)
            qsort(entries,(size_t)m,sizeof(entries[0]),lis_ilutp_compare_desc);
        kept=m<cap?m:cap;
        U->nnz[i]=kept;
        if(kept)
        {
            U->index[i]=(LIS_INT *)malloc((size_t)kept*sizeof(LIS_INT));
            U->value[i]=(LIS_SCALAR *)malloc((size_t)kept*sizeof(LIS_SCALAR));
            if(!U->index[i] || !U->value[i])
            {
                err=LIS_OUT_OF_MEMORY;
                goto cleanup;
            }
            for(t=0;t<kept;t++)
            {
                U->index[i][t]=entries[t].col;
                U->value[i][t]=entries[t].value;
            }
        }
    }
    /* Convert stored original column IDs into final upper positions. */
    for(i=0;i<n;i++)
        for(j=0;j<U->nnz[i];j++)
        {
            col=inv[U->index[i][j]];
            if(col<=i)
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUTP upper index mapping violates triangular structure");
                err=LIS_ERR_ILL_ARG;
                goto cleanup;
            }
            U->index[i][j]=col;
        }

    precon->L=L; L=NULL;
    precon->U=U; U=NULL;
    precon->D=D; D=NULL;
    precon->ilutp_colperm=perm; perm=NULL;
    precon->ilutp_swap_count=swaps;

cleanup:
    lis_matrix_ilu_destroy(L);
    lis_matrix_ilu_destroy(U);
    lis_vector_destroy(D);
    lis_free(perm);
    lis_free(inv);
    lis_free(seen);
    lis_free(lower);
    lis_free(finished);
    lis_free(upper);
    lis_free(w);
    lis_free(entries);
    return err;
#endif
}
