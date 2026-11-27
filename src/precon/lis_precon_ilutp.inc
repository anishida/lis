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

/* Maintain the cap strongest entries using a worst-at-root heap.
 * Deterministic tie-break (column ascending) matches the old full qsort.
 * Sorting only retained entries reduces limited-fill selection to
 * O(m log(cap) + cap log(cap)), instead of O(m log(m)).
 */
static void
lis_ilutp_worst_heap_down(LIS_ILUTP_ENTRY *entry,LIS_INT pos,LIS_INT n)
{
    LIS_ILUTP_ENTRY held=entry[pos];
    while(pos<n/2)
    {
        LIS_INT child=2*pos+1;
        if(child>=n) break;
        if(child+1<n &&
           lis_ilutp_compare_desc(&entry[child+1],&entry[child])>0)
            child++;
        if(lis_ilutp_compare_desc(&entry[child],&held)<=0) break;
        entry[pos]=entry[child];
        pos=child;
    }
    entry[pos]=held;
}
static LIS_INT
lis_ilutp_keep_largest(LIS_ILUTP_ENTRY *entry,LIS_INT m,LIS_INT cap)
{
    LIS_INT i,kept;
    kept=m<cap?m:cap;
    if(kept<=0) return 0;
    if(m>kept)
    {
        for(i=kept/2;i>0;i--)
            lis_ilutp_worst_heap_down(entry,i-1,kept);
        for(i=kept;i<m;i++)
        {
            if(lis_ilutp_compare_desc(&entry[i],&entry[0])<0)
            {
                entry[0]=entry[i];
                lis_ilutp_worst_heap_down(entry,0,kept);
            }
        }
    }
    qsort(entry,(size_t)kept,sizeof(entry[0]),lis_ilutp_compare_desc);
    return kept;
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

/* Experimental fixed-pattern full LU for symmetric CSR structure.
 * Numeric values need not be symmetric. No pivoting or dropping.
 * Separate symbolic reach allows a marker-free numeric scatter.
 */
#if !defined(_COMPLEX) && !defined(USE_COMPLEX)
static int lis_ilutp_cmp_index(const void *a,const void *b)
{
    LIS_INT x=*(const LIS_INT*)a,y=*(const LIS_INT*)b;
    return (x>y)-(x<y);
}
/* Each row writes its own factors, reads only completed dependency rows. */
static LIS_INT lis_ilutp_numeric_row(LIS_MATRIX A,LIS_MATRIX_ILU L,LIS_MATRIX_ILU U,LIS_VECTOR D,
  LIS_INT i,LIS_SCALAR *w,const unsigned char *pairable,const LIS_INT *sp,const LIS_INT *se,
  int runs,int pairs,int groups,int group8,uint64_t *pairs_used,uint64_t *groups3,uint64_t *groups4,uint64_t *groups8,const LIS_INT *perm,const LIS_INT *inv)
{
    LIS_INT n=A->n,j,k,t,c; LIS_SCALAR f,pivot;
        if(perm) {
            LIS_INT old=perm[i];
            for(k=A->ptr[old];k<A->ptr[old+1];k++) {
                LIS_SCALAR v=A->value[k];
                if(!isfinite((double)v))return LIS_ERR_ILL_ARG;
                w[inv[A->index[k]]]=v;
            }
        } else {
            for(k=A->ptr[i];k<A->ptr[i+1];k++) w[A->index[k]]=A->value[k];
        }
        for(t=0;t<L->nnz[i];t++) {
            j=L->index[i][t]; f=w[j]*D->value[j];
            if(!isfinite((double)f)) {return LIS_BREAKDOWN;}
            L->value[i][t]=f; w[j]=0;
            /* Eight-row fused update, with the same sequential arithmetic. */
            if(group8 && groups && pairs && f!=0 && j+7<n && t+7<L->nnz[i] &&
               pairable[j+0] && pairable[j+1] && pairable[j+2] && pairable[j+3] && pairable[j+4] && pairable[j+5] && pairable[j+6] &&
               L->index[i][t+1]==j+1 && L->index[i][t+2]==j+2 && L->index[i][t+3]==j+3 && L->index[i][t+4]==j+4 && L->index[i][t+5]==j+5 && L->index[i][t+6]==j+6 && L->index[i][t+7]==j+7) {
                const LIS_SCALAR *v0=U->value[j+0]; const LIS_SCALAR *v1=U->value[j+1]; const LIS_SCALAR *v2=U->value[j+2]; const LIS_SCALAR *v3=U->value[j+3]; const LIS_SCALAR *v4=U->value[j+4]; const LIS_SCALAR *v5=U->value[j+5]; const LIS_SCALAR *v6=U->value[j+6]; const LIS_SCALAR *v7=U->value[j+7];
                LIS_SCALAR f1,f2,f3,f4,f5,f6,f7;
                const LIS_INT *idx=U->index[j+7]; LIS_INT len=U->nnz[j+7];
                w[j+1]-=f*v0[0];
                f1=w[j+1]*D->value[j+1];
                if(!isfinite((double)f1)) {return LIS_BREAKDOWN;}
                L->value[i][t+1]=f1; w[j+1]=0;
                w[j+2]-=f*v0[1];
                if(f1!=0) w[j+2]-=f1*v1[0];
                f2=w[j+2]*D->value[j+2];
                if(!isfinite((double)f2)) {return LIS_BREAKDOWN;}
                L->value[i][t+2]=f2; w[j+2]=0;
                w[j+3]-=f*v0[2];
                if(f1!=0) w[j+3]-=f1*v1[1];
                if(f2!=0) w[j+3]-=f2*v2[0];
                f3=w[j+3]*D->value[j+3];
                if(!isfinite((double)f3)) {return LIS_BREAKDOWN;}
                L->value[i][t+3]=f3; w[j+3]=0;
                w[j+4]-=f*v0[3];
                if(f1!=0) w[j+4]-=f1*v1[2];
                if(f2!=0) w[j+4]-=f2*v2[1];
                if(f3!=0) w[j+4]-=f3*v3[0];
                f4=w[j+4]*D->value[j+4];
                if(!isfinite((double)f4)) {return LIS_BREAKDOWN;}
                L->value[i][t+4]=f4; w[j+4]=0;
                w[j+5]-=f*v0[4];
                if(f1!=0) w[j+5]-=f1*v1[3];
                if(f2!=0) w[j+5]-=f2*v2[2];
                if(f3!=0) w[j+5]-=f3*v3[1];
                if(f4!=0) w[j+5]-=f4*v4[0];
                f5=w[j+5]*D->value[j+5];
                if(!isfinite((double)f5)) {return LIS_BREAKDOWN;}
                L->value[i][t+5]=f5; w[j+5]=0;
                w[j+6]-=f*v0[5];
                if(f1!=0) w[j+6]-=f1*v1[4];
                if(f2!=0) w[j+6]-=f2*v2[3];
                if(f3!=0) w[j+6]-=f3*v3[2];
                if(f4!=0) w[j+6]-=f4*v4[1];
                if(f5!=0) w[j+6]-=f5*v5[0];
                f6=w[j+6]*D->value[j+6];
                if(!isfinite((double)f6)) {return LIS_BREAKDOWN;}
                L->value[i][t+6]=f6; w[j+6]=0;
                w[j+7]-=f*v0[6];
                if(f1!=0) w[j+7]-=f1*v1[5];
                if(f2!=0) w[j+7]-=f2*v2[4];
                if(f3!=0) w[j+7]-=f3*v3[3];
                if(f4!=0) w[j+7]-=f4*v4[2];
                if(f5!=0) w[j+7]-=f5*v5[1];
                if(f6!=0) w[j+7]-=f6*v6[0];
                f7=w[j+7]*D->value[j+7];
                if(!isfinite((double)f7)) {return LIS_BREAKDOWN;}
                L->value[i][t+7]=f7; w[j+7]=0;
                t+=7; (*groups8)++;
                if(f1!=0 && f2!=0 && f3!=0 && f4!=0 && f5!=0 && f6!=0 && f7!=0) for(c=0;c<len;c++) { LIS_SCALAR z=w[idx[c]];
                    z-=f*v0[c+7];
                    z-=f1*v1[c+6];
                    z-=f2*v2[c+5];
                    z-=f3*v3[c+4];
                    z-=f4*v4[c+3];
                    z-=f5*v5[c+2];
                    z-=f6*v6[c+1];
                    z-=f7*v7[c+0];
                    w[idx[c]]=z;
                } else for(c=0;c<len;c++) { LIS_SCALAR z=w[idx[c]];
                    z-=f*v0[c+7];
                    if(f1!=0) z-=f1*v1[c+6];
                    if(f2!=0) z-=f2*v2[c+5];
                    if(f3!=0) z-=f3*v3[c+4];
                    if(f4!=0) z-=f4*v4[c+3];
                    if(f5!=0) z-=f5*v5[c+2];
                    if(f6!=0) z-=f6*v6[c+1];
                    if(f7!=0) z-=f7*v7[c+0];
                    w[idx[c]]=z;
                }
                continue;
            }
            /* Fused three/four-row updates: chain of identical trailing
             * patterns, and consecutive entries in this lower row. */
            if(groups && pairs && f!=0 && j+2<n && t+2<L->nnz[i] &&
               pairable[j] && pairable[j+1] && L->index[i][t+2]==j+2 && L->index[i][t+1]==j+1) {
                const LIS_SCALAR *v1=U->value[j],*v2=U->value[j+1],*v3=U->value[j+2];
                const LIS_INT *idx; LIS_INT len; LIS_SCALAR g,h;
                w[j+1]-=f*v1[0]; g=w[j+1]*D->value[j+1];
                w[j+2]-=f*v1[1]; if(g!=0) w[j+2]-=g*v2[0]; h=w[j+2]*D->value[j+2];
                if(!isfinite((double)g)||!isfinite((double)h)) {return LIS_BREAKDOWN;}
                L->value[i][++t]=g; L->value[i][++t]=h; w[j+1]=0; w[j+2]=0;
                if(j+3<n && t+1<L->nnz[i] && pairable[j+2] && L->index[i][t+1]==j+3) {
                    const LIS_SCALAR *v4=U->value[j+3]; LIS_SCALAR q;
                    w[j+3]-=f*v1[2]; if(g!=0) w[j+3]-=g*v2[1]; if(h!=0) w[j+3]-=h*v3[0]; q=w[j+3]*D->value[j+3];
                    if(!isfinite((double)q)) {return LIS_BREAKDOWN;}
                    L->value[i][++t]=q; w[j+3]=0; (*groups4)++;
                    idx=U->index[j+3]; len=U->nnz[j+3];
                    if(g!=0 && h!=0 && q!=0) for(c=0;c<len;c++) {LIS_SCALAR z=w[idx[c]];z-=f*v1[c+3];z-=g*v2[c+2];z-=h*v3[c+1];z-=q*v4[c];w[idx[c]]=z;}
                    else for(c=0;c<len;c++) {LIS_SCALAR z=w[idx[c]];z-=f*v1[c+3];if(g!=0) z-=g*v2[c+2];if(h!=0) z-=h*v3[c+1];if(q!=0) z-=q*v4[c];w[idx[c]]=z;}
                } else {
                    (*groups3)++; idx=U->index[j+2]; len=U->nnz[j+2];
                    if(g!=0 && h!=0) for(c=0;c<len;c++) {LIS_SCALAR z=w[idx[c]];z-=f*v1[c+2];z-=g*v2[c+1];z-=h*v3[c];w[idx[c]]=z;}
                    else for(c=0;c<len;c++) {LIS_SCALAR z=w[idx[c]];z-=f*v1[c+2];if(g!=0) z-=g*v2[c+1];if(h!=0) z-=h*v3[c];w[idx[c]]=z;}
                }
                continue;
            }
            if(pairs && pairable[j] && t+1<L->nnz[i] && L->index[i][t+1]==j+1 && f!=0) {
                const LIS_SCALAR *v1=U->value[j],*v2=U->value[j+1];
                const LIS_INT *idx=U->index[j+1]; LIS_INT len=U->nnz[j+1];
                LIS_SCALAR g;
                w[j+1]-=f*v1[0]; g=w[j+1]*D->value[j+1];
                if(!isfinite((double)g)) {return LIS_BREAKDOWN;}
                L->value[i][++t]=g; w[j+1]=0; (*pairs_used)++;
                if(g!=0) for(c=0;c<len;c++) {LIS_SCALAR z=w[idx[c]]; z-=f*v1[c+1]; z-=g*v2[c]; w[idx[c]]=z;}
                else for(c=0;c<len;c++) w[idx[c]]-=f*v1[c+1];
                continue;
            }
            if(f!=0) {
                LIS_INT len=U->nnz[j]; const LIS_INT *idx=U->index[j]; const LIS_SCALAR *val=U->value[j];
                if(runs) {
                    LIS_INT seg,beg=0,stop;
                    for(seg=sp[j];seg<sp[j+1];seg++) {
                        stop=se[seg];
                        if(stop>0) {
                            LIS_SCALAR *dst=w+idx[beg]; const LIS_SCALAR *src=val+beg;
                            LIS_INT count_run=stop-beg;
                            #pragma omp simd
                            for(c=0;c<count_run;c++) dst[c]-=f*src[c];
                        } else {
                            stop=-stop;
                            for(c=beg;c<stop;c++) w[idx[c]]-=f*val[c];
                        }
                        beg=stop;
                    }
                } else for(c=0;c<len;c++) w[idx[c]]-=f*val[c];
            }
        }
        pivot=w[i]; w[i]=0;
        if(!isfinite((double)pivot)||pivot==0||!isfinite(1.0/(double)pivot)) {return LIS_BREAKDOWN;}
        D->value[i]=1.0/pivot;
        for(t=0;t<U->nnz[i];t++) {
            j=U->index[i][t]; f=w[j]; w[j]=0;
            if(!isfinite((double)f)) {return LIS_BREAKDOWN;}
            U->value[i][t]=f;
        }
    return LIS_SUCCESS;
}
static LIS_INT lis_ilutp_full_symbolic(LIS_MATRIX A,LIS_PRECON precon,uint64_t limit)
{
    LIS_INT n=A->n,i,j,k,t,m,lo,hi,mid,err=LIS_SUCCESS;
    LIS_INT *parent=NULL,*mark=NULL,*reach=NULL,*count=NULL,*perm=NULL;
    LIS_MATRIX_ILU L=NULL,U=NULL;
    LIS_VECTOR D=NULL;
    LIS_SCALAR *w=NULL;
    uint64_t used=0,segments=0,contiguous=0;
    LIS_INT *sp=NULL,*se=NULL,si,start,end,sc,k0;
    double clock0=lis_wtime(),validate_time=0,symbolic_time=0,segment_time=0,numeric_time=0;
    unsigned char *pairable=NULL;
    uint64_t pairs_used=0,pairs_possible=0,groups4=0,groups3=0,groups8=0;
    int group8=getenv("LIS_EXPERIMENT_ILUTP_GROUP8")!=NULL;
    int groups=getenv("LIS_EXPERIMENT_ILUTP_GROUPS")!=NULL;
    int pairs=getenv("LIS_EXPERIMENT_ILUTP_PAIRS")!=NULL;
    int runs=getenv("LIS_EXPERIMENT_ILUTP_RUNS")!=NULL;
    /* Only sorted, duplicate-free, structurally symmetric CSR is supported.
     * Validate all rows before binary-searching transposed entries. */
    for(i=0;i<n;i++) {
        if(A->ptr[i]<0 || A->ptr[i+1]<A->ptr[i] || A->ptr[i+1]>A->nnz)
            return LIS_ERR_ILL_ARG;
        for(k=A->ptr[i];k<A->ptr[i+1];k++) {
            j=A->index[k];
            if(j<0 || j>=n || !isfinite((double)A->value[k])) return LIS_ERR_ILL_ARG;
            if(k>A->ptr[i] && A->index[k-1]>=j) return LIS_ERR_NOT_IMPLEMENTED;
        }
    }
    for(i=0;i<n;i++) for(k=A->ptr[i];k<A->ptr[i+1];k++) {
        j=A->index[k]; if(j==i) continue;
        lo=A->ptr[j]; hi=A->ptr[j+1];
        while(lo<hi) {mid=lo+(hi-lo)/2; if(A->index[mid]<i) lo=mid+1; else hi=mid;}
        if(lo==A->ptr[j+1] || A->index[lo]!=i) return LIS_ERR_NOT_IMPLEMENTED;
    }
    validate_time=lis_wtime()-clock0; clock0=lis_wtime();
    parent=lis_malloc((size_t)n*sizeof(*parent),"symbolic LU parent");
    mark=lis_malloc((size_t)n*sizeof(*mark),"symbolic LU mark");
    reach=lis_malloc((size_t)n*sizeof(*reach),"symbolic LU reach");
    count=lis_malloc((size_t)n*sizeof(*count),"symbolic LU count");
    perm=lis_malloc((size_t)n*sizeof(*perm),"symbolic LU identity");
    w=lis_malloc((size_t)n*sizeof(*w),"symbolic LU workspace");
    if(!parent||!mark||!reach||!count||!perm||!w) {err=LIS_OUT_OF_MEMORY;goto done;}
    for(i=0;i<n;i++) {parent[i]=-1; mark[i]=-1; count[i]=0; perm[i]=i; w[i]=0;}
    if((err=lis_matrix_ilu_create(n,1,&L)) || (err=lis_matrix_ilu_create(n,1,&U)) ||
       (err=lis_matrix_ilu_setCR(L)) || (err=lis_matrix_ilu_setCR(U)) ||
       (err=lis_vector_duplicate(A,&D))) goto done;
    /* Compute elimination tree and each lower row pattern. Every visited
     * node is a structural lower entry. Parent always exceeds child. */
    for(i=0;i<n;i++) {
        m=0; mark[i]=i;
        for(k=A->ptr[i];k<A->ptr[i+1] && A->index[k]<i;k++) {
            j=A->index[k];
            while(mark[j]!=i) {
                mark[j]=i; reach[m++]=j;
                if(parent[j]<0) parent[j]=i;
                j=parent[j];
            }
        }
        if((uint64_t)m>(limit-used)/2) {err=LIS_OUT_OF_MEMORY;goto done;}
        used+=2*(uint64_t)m;
        if(m) {
            qsort(reach,(size_t)m,sizeof(*reach),lis_ilutp_cmp_index);
            L->index[i]=malloc((size_t)m*sizeof(LIS_INT));
            L->value[i]=malloc((size_t)m*sizeof(LIS_SCALAR));
            if(!L->index[i]||!L->value[i]) {L->nnz[i]=m;err=LIS_OUT_OF_MEMORY;goto done;}
            L->nnz[i]=m;
            memcpy(L->index[i],reach,(size_t)m*sizeof(*reach));
            for(t=0;t<m;t++) count[reach[t]]++;
        }
    }
    for(i=0;i<n;i++) if(count[i]) {
        U->index[i]=malloc((size_t)count[i]*sizeof(LIS_INT));
        U->value[i]=malloc((size_t)count[i]*sizeof(LIS_SCALAR));
        U->nnz[i]=count[i];
        if(!U->index[i]||!U->value[i]) {err=LIS_OUT_OF_MEMORY;goto done;}
    }
    memset(count,0,(size_t)n*sizeof(*count));
    for(i=0;i<n;i++) for(t=0;t<L->nnz[i];t++) {
        j=L->index[i][t]; U->index[j][count[j]++]=i;
    }
    symbolic_time=lis_wtime()-clock0; clock0=lis_wtime();
    if(runs) {
        sp=lis_malloc((size_t)(n+1)*sizeof(*sp),"LU segment offsets");
        if(!sp) {err=LIS_OUT_OF_MEMORY;goto done;}
        /* Group adjacent columns into vectorizable runs; combine remaining
         * indices into scatter segments. Positive end=contiguous, negative=scatter. */
        for(sc=0;sc<2;sc++) {
            si=0;
            for(i=0;i<n;i++) {
                const LIS_INT *idx=U->index[i]; LIS_INT len=U->nnz[i];
                sp[i]=si; k0=0; start=0;
                while(start<len) {
                    end=start+1;
                    while(end<len && idx[end]==idx[end-1]+1) end++;
                    if(end-start>=4) {
                        if(k0<start) {if(sc) se[si]=-start;si++;}
                        if(sc) se[si]=end;si++;
                        if(!sc) contiguous+=(uint64_t)(end-start);
                        k0=end;
                    }
                    start=end;
                }
                if(k0<len) {if(sc) se[si]=-len;si++;}
            }
            sp[n]=si;
            if(!sc) {
                segments=si; se=lis_malloc((size_t)(si?si:1)*sizeof(*se),"LU segment boundaries");
                if(!se) {err=LIS_OUT_OF_MEMORY;goto done;}
            }
        }
    }
    if(pairs) {
        pairable=lis_malloc((size_t)n,"LU adjacent row flags");
        if(!pairable) {err=LIS_OUT_OF_MEMORY;goto done;}
        memset(pairable,0,(size_t)n);
        for(j=0;j<n-1;j++) {
            LIS_INT len=U->nnz[j];
            if(len>=2 && U->index[j][0]==j+1 && len==U->nnz[j+1]+1 &&
               !memcmp(U->index[j]+1,U->index[j+1],(size_t)(len-1)*sizeof(LIS_INT))) {
                pairable[j]=1; pairs_possible++;
            }
        }
    }
    segment_time=lis_wtime()-clock0; clock0=lis_wtime();
    /* Numeric LU: symbolic reach guarantees every touched entry belongs
     * to this row's L/diagonal/U. Clear all those entries after use.
     * Check finite coefficients and pivots before exposing the factors. */
#ifdef _OPENMP
    if(getenv("LIS_EXPERIMENT_ILUTP_PARALLEL_NUMERIC") && omp_get_max_threads()>1) {
        LIS_INT *level=parent,*order=reach,*offset=NULL,*cursor=NULL,depth=0,lev,nt=omp_get_max_threads();
        LIS_SCALAR *work=NULL;
        uint64_t work_total=0,work_critical=0;
        /* Symbolic arrays are no longer needed. Compute longest dependency depth. */
        for(i=0;i<n;i++) {
            lev=0;
            for(t=0;t<L->nnz[i];t++) { j=L->index[i][t]; if(lev<level[j]+1) lev=level[j]+1; }
            level[i]=lev; if(depth<lev+1) depth=lev+1;
        }
        offset=calloc((size_t)depth+1,sizeof(*offset)); cursor=calloc((size_t)depth,sizeof(*cursor));
        if(!offset||!cursor || (size_t)nt>SIZE_MAX/sizeof(*work)/(size_t)n || !(work=calloc((size_t)nt*n,sizeof(*work)))) {
            free(offset);free(cursor);free(work);err=LIS_OUT_OF_MEMORY;goto done;
        }
        for(i=0;i<n;i++) offset[level[i]+1]++;
        for(lev=0;lev<depth;lev++) {offset[lev+1]+=offset[lev];cursor[lev]=offset[lev];}
        for(i=0;i<n;i++) order[cursor[level[i]]++]=i;
        for(lev=0;lev<depth;lev++) {
            uint64_t maxcost=0;
            for(k=offset[lev];k<offset[lev+1];k++) {
                uint64_t cost=0; i=order[k];
                for(t=0;t<L->nnz[i];t++) cost+=(uint64_t)U->nnz[L->index[i][t]];
                work_total+=cost;if(cost>maxcost)maxcost=cost;
            } work_critical+=maxcost;
        }
        fprintf(stderr,"LU_LEVELS depth=%d max_threads=%d work=%llu critical=%llu\n",(int)depth,(int)nt,(unsigned long long)work_total,(unsigned long long)work_critical);
        /* Implicit barriers enforce dependencies between adjacent levels.
         * A failure is reduced across the whole level before exiting. */
        #pragma omp parallel reduction(+:pairs_used,groups3,groups4,groups8)
        {
            LIS_SCALAR *local=work+(size_t)omp_get_thread_num()*n;
            for(LIS_INT l=0;l<depth;l++) {
                #pragma omp for schedule(dynamic,4) reduction(max:err)
                for(LIS_INT pos=offset[l];pos<offset[l+1];pos++) {
                    LIS_INT e=lis_ilutp_numeric_row(A,L,U,D,order[pos],local,pairable,sp,se,runs,pairs,groups,group8,&pairs_used,&groups3,&groups4,&groups8,NULL,NULL);
                    if(e>err) err=e;
                }
                if(err!=LIS_SUCCESS) break;
            }
        }
        free(offset);free(cursor);free(work);
        if(err!=LIS_SUCCESS) goto done;
    } else
#endif
    for(i=0;i<n;i++) {
        err=lis_ilutp_numeric_row(A,L,U,D,i,w,pairable,sp,se,runs,pairs,groups,group8,
                                 &pairs_used,&groups3,&groups4,&groups8,NULL,NULL);
        if(err) goto done;
    }
    numeric_time=lis_wtime()-clock0;
    fprintf(stderr,"LU_PROFILE runs=%d validation=%.6f symbolic=%.6f segment=%.6f numeric=%.6f segments=%llu contiguous=%llu upper=%llu\n",runs,validate_time,symbolic_time,segment_time,numeric_time,(unsigned long long)segments,(unsigned long long)contiguous,(unsigned long long)(used/2));
    fprintf(stderr,"LU_PAIRS mode=%d possible=%llu used=%llu\n",pairs,(unsigned long long)pairs_possible,(unsigned long long)pairs_used);
    fprintf(stderr,"LU_GROUPS mode=%d triples=%llu quads=%llu\n",groups,(unsigned long long)groups3,(unsigned long long)groups4);
    fprintf(stderr,"LU_GROUP8 mode=%d used=%llu\n",group8,(unsigned long long)groups8);
    precon->L=L;L=NULL;precon->U=U;U=NULL;precon->D=D;D=NULL;
    precon->ilutp_colperm=perm;perm=NULL;precon->ilutp_swap_count=0;
 done:
    lis_matrix_ilu_destroy(L);lis_matrix_ilu_destroy(U);lis_vector_destroy(D);
    lis_free(sp);lis_free(se);lis_free(pairable);
    lis_free(parent);lis_free(mark);lis_free(reach);lis_free(count);lis_free(perm);lis_free(w);
    return err;
}
#endif

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
    if(omp_get_max_threads()!=1 && !(full && getenv("LIS_EXPERIMENT_ILUTP_SYMBOLIC") && getenv("LIS_EXPERIMENT_ILUTP_SERIAL_APPLY")))
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
    if(full && tau==0.0 && getenv("LIS_EXPERIMENT_ILUTP_SYMBOLIC"))
    {
        err=lis_ilutp_full_symbolic(A,precon,fill_limit);
        if(err==LIS_SUCCESS) fprintf(stderr,"ILUTP_SYMBOLIC full LU selected\n");
        else if(err==LIS_ERR_NOT_IMPLEMENTED) fprintf(stderr,"ILUTP_SYMBOLIC fallback to dynamic LU\n");
        if(err!=LIS_ERR_NOT_IMPLEMENTED) return err;
        /* A general CSR falls back to the original dynamic algorithm. */
    }
#ifdef _OPENMP
    /* Only the fixed-pattern no-pivot path is extended to multi-thread
     * outer solvers. Never enter a legacy block factorization here. */
    if(omp_get_max_threads()!=1) return LIS_ERR_NOT_IMPLEMENTED;
#endif
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
                LIS_SETERR1(LIS_BREAKDOWN,
                    "ILUTP lower multiplier became nonfinite at local factor row %D (0-based)\n",i);
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
                    LIS_SETERR1(LIS_BREAKDOWN,
                        "ILUTP elimination became nonfinite at local factor row %D (0-based)\n",i);
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
            /* A nonzero input diagonal does not guarantee a usable
             * pivot after sparse dropping and column swaps.
             * This is the local, zero-based factorization-order row,
             * not the original/global (or Fortran 1-based) matrix row. */
            LIS_SETERR1(LIS_BREAKDOWN,
                "ILUTP found no usable pivot at local factor row %D (0-based); "
                "consider fill/drop and explicit pivot regularization\n",i);
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
        kept=full ? m : lis_ilutp_keep_largest(entries,m,cap);
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
        kept=full ? m : lis_ilutp_keep_largest(entries,m,cap);
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
