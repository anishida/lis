/* Experimental caller-owned session for fixed-pattern full ILUTP.
 * No change to public LIS structures. No shared cache. The session owns
 * its factors and a byte-exact copy of the original CSR pattern.
 * Any failed refresh invalidates the factors; destroy/recreate session.
 * Solver and matrices must remain serial real-double CSR. No scaling.
 */
#if !defined(USE_MPI) && !defined(_COMPLEX) && !defined(USE_COMPLEX)
typedef struct LIS_ILUTP_REUSE_STRUCT {
    LIS_INT n, nnz;
    LIS_INT *ptr, *index;
    LIS_PRECON factors;
} *LIS_ILUTP_REUSE;

static int reuse_config(LIS_SOLVER s) {
    return s && sizeof(LIS_INT)==4 && sizeof(LIS_SCALAR)==8 &&
        s->options[LIS_OPTIONS_PRECON]==LIS_PRECON_TYPE_ILUT &&
        s->options[LIS_OPTIONS_PRECISION]==LIS_PRECISION_DOUBLE &&
        s->options[LIS_OPTIONS_SCALE]==LIS_SCALE_NONE &&
        !s->options[LIS_OPTIONS_ADDS] && s->ilutp_enabled && s->ilutp_full &&
        s->ilutp_threshold==0 && !s->ilut_rcm && !s->ilut_row_scale &&
        !s->ilut_equil && s->ilut_ras_overlap<0 && !s->ilut_schur &&
        s->precon_half_mode==LIS_PRECON_HALF_OFF &&
        getenv("LIS_EXPERIMENT_ILUTP_MD") && getenv("LIS_EXPERIMENT_ILUTP_SYMBOLIC");
}
static int reuse_matrix(LIS_MATRIX A) {
    return A && A->matrix_type==LIS_MATRIX_CSR && A->nprocs==1 && A->is==0 &&
        A->n==A->gn && A->n>0 && A->nnz>0 && A->ptr && A->index && A->value &&
        A->ptr[0]==0 && A->ptr[A->n]==A->nnz;
}
void lis_ilutp_reuse_destroy(LIS_ILUTP_REUSE c) {
    if(!c)return;
    if(c->factors)lis_precon_destroy(c->factors);
    free(c->ptr);free(c->index);free(c);
}
LIS_INT lis_ilutp_reuse_create(LIS_SOLVER s,LIS_MATRIX A,LIS_ILUTP_REUSE *out) {
    LIS_ILUTP_REUSE c; LIS_INT err; LIS_MATRIX old;
    if(!out)return LIS_ERR_ILL_ARG;
    *out=NULL;
    if(!reuse_config(s)||!reuse_matrix(A))return LIS_ERR_NOT_IMPLEMENTED;
    c=calloc(1,sizeof(*c));if(!c)return LIS_OUT_OF_MEMORY;
    c->n=A->n;c->nnz=A->nnz;
    c->ptr=malloc(((size_t)c->n+1)*sizeof(LIS_INT));
    c->index=malloc((size_t)c->nnz*sizeof(LIS_INT));
    if(!c->ptr||!c->index){lis_ilutp_reuse_destroy(c);return LIS_OUT_OF_MEMORY;}
    memcpy(c->ptr,A->ptr,((size_t)c->n+1)*sizeof(LIS_INT));
    memcpy(c->index,A->index,(size_t)c->nnz*sizeof(LIS_INT));
    old=s->A;s->A=A;err=lis_precon_create(s,&c->factors);s->A=old;
    /* lis_precon_create destroys its output on failure in this version. */
    if(err){c->factors=NULL;lis_ilutp_reuse_destroy(c);return err;}
    if(!c->factors->ilut_perm||!c->factors->ilutp_colperm||c->factors->ilutp_swap_count||!c->factors->L||!c->factors->U||!c->factors->D){lis_ilutp_reuse_destroy(c);return LIS_ERR_NOT_IMPLEMENTED;}
    *out=c;return LIS_SUCCESS;
}
LIS_PRECON lis_ilutp_reuse_factors(LIS_ILUTP_REUSE c){return c?c->factors:NULL;}
static LIS_INT reuse_refresh_impl(LIS_ILUTP_REUSE c,LIS_SOLVER s,LIS_MATRIX A) {
    LIS_INT err=LIS_ERR_ILL_ARG,i,j,n,*inv=NULL;
    LIS_PRECON p;
    LIS_SCALAR *w=NULL;unsigned char *pairable=NULL;
    uint64_t fill=0,pairs_used=0,g3=0,g4=0,g8=0,limit;
    int pairs=getenv("LIS_EXPERIMENT_ILUTP_PAIRS")!=NULL;
    int groups=getenv("LIS_EXPERIMENT_ILUTP_GROUPS")!=NULL;
    int group8=getenv("LIS_EXPERIMENT_ILUTP_GROUP8")!=NULL;
    double t=lis_wtime(),verify=0,permute=0,numeric=0;
    if(!c||!(p=c->factors))return LIS_ERR_ILL_ARG;
    n=c->n;
    if(!reuse_config(s)||!reuse_matrix(A)||A->n!=n||A->nnz!=c->nnz||s->ilutp_full_limit<0)goto done;
    if(memcmp(c->ptr,A->ptr,((size_t)n+1)*sizeof(LIS_INT))||memcmp(c->index,A->index,(size_t)c->nnz*sizeof(LIS_INT)))goto done;
    for(i=0;i<n;i++)fill+=(uint64_t)p->L->nnz[i]+p->U->nnz[i];
    limit=s->ilutp_full_limit?(uint64_t)s->ilutp_full_limit:(uint64_t)c->nnz*20;
    if(!s->ilutp_full_limit){if(limit>INT32_MAX)limit=INT32_MAX;if(limit<(uint64_t)n*4)limit=(uint64_t)n*4;}
    if(fill>limit){err=LIS_OUT_OF_MEMORY;goto done;}
    verify=lis_wtime()-t;t=lis_wtime();
    inv=malloc((size_t)n*sizeof(*inv));w=calloc((size_t)n,sizeof(*w));pairable=calloc((size_t)n,1);
    if(!inv||!w||!pairable){err=LIS_OUT_OF_MEMORY;goto done;}
    for(i=0;i<n;i++)inv[p->ilut_perm[i]]=i;
    if(pairs)for(j=0;j<n-1;j++){
        LIS_INT len=p->U->nnz[j];
        if(len>=2 && p->U->index[j][0]==j+1 && len==p->U->nnz[j+1]+1 &&
           !memcmp(p->U->index[j]+1,p->U->index[j+1],(size_t)(len-1)*sizeof(LIS_INT)))pairable[j]=1;
    }
    permute=lis_wtime()-t;t=lis_wtime();
    /* Reuse complete symbolic pattern; sequential row kernel, no RUNS. */
    for(i=0;i<n;i++)if((err=lis_ilutp_numeric_row(A,p->L,p->U,p->D,i,w,pairable,NULL,NULL,0,pairs,groups,group8,&pairs_used,&g3,&g4,&g8,p->ilut_perm,inv)))goto done;
    numeric=lis_wtime()-t;
    fprintf(stderr,"REUSE verify=%.6f permute=%.6f numeric=%.6f fill=%llu\n",verify,permute,numeric,(unsigned long long)fill);
 done:
    free(inv);free(w);free(pairable);
    return err;
}
LIS_INT lis_ilutp_reuse_refresh(LIS_ILUTP_REUSE c,LIS_SOLVER s,LIS_MATRIX A) {
    LIS_INT err=reuse_refresh_impl(c,s,A);
    if(err && c && c->factors){if(s && s->precon==c->factors)s->precon=NULL;lis_precon_destroy(c->factors);c->factors=NULL;}
    return err;
}
/* The built-in ILUT preconditioner owns its snapshot; no factor ownership
 * cycle. Public user-data access must reject this reserved state. */
static LIS_INT reuse_attached_destroy(LIS_PRECON p) {
    LIS_ILUTP_REUSE c=p->user_data;
    p->user_data=NULL;p->user_destroy=NULL;
    if(c){c->factors=NULL;lis_ilutp_reuse_destroy(c);}
    return LIS_SUCCESS;
}
static LIS_INT reuse_attach(LIS_SOLVER s,LIS_PRECON p) {
    LIS_MATRIX A=s->A; LIS_ILUTP_REUSE c;
    if(!getenv("LIS_EXPERIMENT_ILUTP_REUSE") || !reuse_config(s) || !reuse_matrix(A))return LIS_SUCCESS;
    if(p->precon_type!=LIS_PRECON_TYPE_ILUT || p->user_data || p->user_destroy || !p->ilut_perm || !p->ilutp_colperm || p->ilutp_swap_count || p->half)return LIS_SUCCESS;
    c=calloc(1,sizeof(*c));if(!c)return LIS_OUT_OF_MEMORY;
    c->n=A->n;c->nnz=A->nnz;
    c->ptr=malloc(((size_t)c->n+1)*sizeof(LIS_INT));c->index=malloc((size_t)c->nnz*sizeof(LIS_INT));
    if(!c->ptr||!c->index){lis_ilutp_reuse_destroy(c);return LIS_OUT_OF_MEMORY;}
    memcpy(c->ptr,A->ptr,((size_t)c->n+1)*sizeof(LIS_INT));memcpy(c->index,A->index,(size_t)c->nnz*sizeof(LIS_INT));
    c->factors=p;p->user_data=c;p->user_destroy=reuse_attached_destroy;
    return LIS_SUCCESS;
}
LIS_INT lis_ilutp_refresh_attached(LIS_SOLVER s,LIS_PRECON p) {
    if(!getenv("LIS_EXPERIMENT_ILUTP_REUSE") || !p || p->user_destroy!=reuse_attached_destroy || !p->user_data)return LIS_ERR_NOT_IMPLEMENTED;
    /* Caller must destroy this preconditioner if update fails. */
    return reuse_refresh_impl(p->user_data,s,s->A);
}
#else
LIS_INT lis_ilutp_refresh_attached(LIS_SOLVER s,LIS_PRECON p) { (void)s;(void)p;return LIS_ERR_NOT_IMPLEMENTED; }
#endif
