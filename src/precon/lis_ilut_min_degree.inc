/* Experimental, serial exact minimum degree, 32-bit indices. All state
 * belongs to this call. Input sorted, unique, structurally symmetric CSR.
 * Permutation new -> old; lower-bound heap, quotient graph.
 */
typedef struct { LIS_INT *a,n,cap; } MDVec;
typedef struct { LIS_INT *heap,*pos,hn; int64_t *key; } MDHeap;
static int md_less(MDHeap *h,LIS_INT a,LIS_INT b) {
    return h->key[a]<h->key[b] || (h->key[a]==h->key[b] && a>b);
}
static void md_swap(MDHeap *h,LIS_INT a,LIS_INT b) {
    LIS_INT v=h->heap[a],u=h->heap[b]; h->heap[a]=u;h->heap[b]=v;h->pos[u]=a;h->pos[v]=b;
}
static void md_up(MDHeap *h,LIS_INT i) {
    while(i && md_less(h,h->heap[i],h->heap[(i-1)/2])) {LIS_INT j=(i-1)/2;md_swap(h,i,j);i=j;}
}
static void md_down(MDHeap *h,LIS_INT i) {
    while(i<h->hn/2) {LIS_INT l=2*i+1;if(l+1<h->hn && md_less(h,h->heap[l+1],h->heap[l])) ++l;
        if(!md_less(h,h->heap[l],h->heap[i])) break;md_swap(h,i,l);i=l;}
}
static LIS_INT md_pop(MDHeap *h) {
    LIS_INT v=h->heap[0];--h->hn;if(h->hn){h->heap[0]=h->heap[h->hn];h->pos[h->heap[0]]=0;md_down(h,0);}h->pos[v]=-1;return v;
}
static LIS_INT md_push(MDVec *v,LIS_INT x,MDVec *el) {
    if(v->n==v->cap) { LIS_INT j=0,i,newcap; LIS_INT *p;
        for(i=0;i<v->n;i++) if(el[v->a[i]].a) v->a[j++]=v->a[i];v->n=j;
        if(v->n==v->cap) {if(v->cap>INT32_MAX/2) return LIS_OUT_OF_MEMORY;newcap=v->cap?2*v->cap:4;
            p=realloc(v->a,(size_t)newcap*sizeof(*p));if(!p)return LIS_OUT_OF_MEMORY;v->a=p;v->cap=newcap;}
    }v->a[v->n++]=x;return LIS_SUCCESS;
}
static LIS_INT lis_ilut_min_degree(LIS_MATRIX A,LIS_INT *perm,uint64_t limit) {
    LIS_INT n=A->n,i,k,j,v,u,l,r,m,d,step,q,stamp=0,err=LIS_SUCCESS;
    const LIS_INT *p=A->ptr,*idx=A->index; MDVec *el=NULL,*inc=NULL;MDHeap h={0};
    LIS_INT *seen=NULL,*mark=NULL,*neighbors=NULL;uint64_t lower=0,checks=0;double t=lis_wtime();
    if(sizeof(LIS_INT)!=4 || n<=0 || n>=INT32_MAX || !p || !idx || p[0] || p[n]!=A->nnz) return LIS_ERR_ILL_ARG;
    if((size_t)n>SIZE_MAX/sizeof(MDVec) || (size_t)n>SIZE_MAX/sizeof(int64_t)) return LIS_OUT_OF_MEMORY;
    for(i=0;i<n;i++) {if(p[i]<0||p[i+1]<p[i]||p[i+1]>A->nnz)return LIS_ERR_ILL_ARG;
        for(k=p[i];k<p[i+1];k++) if(idx[k]<0||idx[k]>=n||(k>p[i]&&idx[k-1]>=idx[k]))return LIS_ERR_ILL_ARG;}
    for(i=0;i<n;i++) for(k=p[i];k<p[i+1];k++) {v=idx[k];l=p[v];r=p[v+1];while(l<r){m=l+(r-l)/2;if(idx[m]<i)l=m+1;else r=m;}
        if(l==p[v+1]||idx[l]!=i)return LIS_ERR_NOT_IMPLEMENTED;}
    el=calloc((size_t)n,sizeof(*el));inc=calloc((size_t)n,sizeof(*inc));h.heap=malloc((size_t)n*sizeof(*h.heap));h.pos=malloc((size_t)n*sizeof(*h.pos));h.key=malloc((size_t)n*sizeof(*h.key));seen=calloc((size_t)n,sizeof(*seen));mark=malloc((size_t)n*sizeof(*mark));neighbors=malloc((size_t)n*sizeof(*neighbors));
    if(!el||!inc||!h.heap||!h.pos||!h.key||!seen||!mark||!neighbors){err=LIS_OUT_OF_MEMORY;goto done;}
    h.hn=n;for(i=0;i<n;i++){h.heap[i]=h.pos[i]=i;mark[i]=-1;h.key[i]=0;for(k=p[i];k<p[i+1];k++)if(idx[k]!=i)h.key[i]++;}for(i=n/2;i-->0;)md_down(&h,i);
    for(step=0;step<n;step++) {
        for(;;){v=h.heap[0];d=0;if(++stamp==INT32_MAX){memset(seen,0,(size_t)n*sizeof(*seen));stamp=1;}seen[v]=stamp;
            for(k=p[v];k<p[v+1];k++){u=idx[k];if(h.pos[u]>=0&&seen[u]!=stamp){seen[u]=stamp;d++;}}
            for(k=0;k<inc[v].n;k++){MDVec *e=&el[inc[v].a[k]];if(!e->a)continue;for(q=0;q<e->n;q++){u=e->a[q];if(h.pos[u]>=0&&seen[u]!=stamp){seen[u]=stamp;d++;}}}
            checks++;if(d<h.key[v]){err=LIS_ERR_ILL_ARG;goto done;}if(d==h.key[v])break;h.key[v]=d;md_down(&h,0);
        }
        v=md_pop(&h);d=0;perm[step]=v;
        for(k=p[v];k<p[v+1];k++){u=idx[k];if(h.pos[u]>=0&&mark[u]!=v){mark[u]=v;neighbors[d++]=u;}}
        for(k=0;k<inc[v].n;k++){MDVec *e=&el[inc[v].a[k]];if(!e->a)continue;for(q=0;q<e->n;q++){u=e->a[q];if(h.pos[u]>=0&&mark[u]!=v){mark[u]=v;neighbors[d++]=u;}}free(e->a);e->a=NULL;e->n=e->cap=0;}
        free(inc[v].a);inc[v].a=NULL;inc[v].n=inc[v].cap=0;lower+=(uint64_t)d;
        if(lower>limit/2){err=LIS_OUT_OF_MEMORY;goto done;}
        if(d==h.hn){for(k=step+1;k<n;k++){perm[k]=md_pop(&h);lower+=(uint64_t)(n-k-1);}if(lower>limit/2)err=LIS_OUT_OF_MEMORY;break;}
        if(d){el[v].a=malloc((size_t)d*sizeof(*el[v].a));if(!el[v].a){err=LIS_OUT_OF_MEMORY;goto done;}memcpy(el[v].a,neighbors,(size_t)d*sizeof(*neighbors));el[v].n=el[v].cap=d;}
        for(q=0;q<d;q++){int64_t old;u=neighbors[q];old=h.key[u];h.key[u]=old-1>d-1?old-1:d-1;if((err=md_push(&inc[u],v,el)))goto done;if(h.key[u]<old)md_up(&h,h.pos[u]);else md_down(&h,h.pos[u]);}
    }
    if(!err) fprintf(stderr,"ILUT_MD n=%d fill=%llu checks=%llu wall=%.6f\n",(int)n,(unsigned long long)(2*lower),(unsigned long long)checks,lis_wtime()-t);
done:
    if(el)for(i=0;i<n;i++)free(el[i].a);if(inc)for(i=0;i<n;i++)free(inc[i].a);free(el);free(inc);free(h.heap);free(h.pos);free(h.key);free(seen);free(mark);free(neighbors);return err;
}
/* Construct sorted permuted CSR without the per-entry matrix insertion API. */
typedef struct { LIS_INT j; LIS_SCALAR v; } MDEntry;
static int md_entry_cmp(const void *a,const void *b){LIS_INT x=((const MDEntry*)a)->j,y=((const MDEntry*)b)->j;return (x>y)-(x<y);}
static LIS_INT lis_ilut_md_permute(LIS_MATRIX A,LIS_MATRIX P,const LIS_INT *perm,const LIS_INT *inv,const LIS_SCALAR *weights,const LIS_SCALAR *colweights) {
    LIS_INT i,k,j,m=0,old,len,off=0,err;LIS_INT *rp=NULL,*rj=NULL;LIS_SCALAR *rv=NULL;MDEntry *tmp=NULL;
    for(i=0;i<A->n;i++)if(A->ptr[i+1]-A->ptr[i]>m)m=A->ptr[i+1]-A->ptr[i];
    if((size_t)m>SIZE_MAX/sizeof(*tmp))return LIS_OUT_OF_MEMORY;
    tmp=malloc((size_t)(m?m:1)*sizeof(*tmp));if(!tmp)return LIS_OUT_OF_MEMORY;
    err=lis_matrix_malloc_csr(A->n,A->nnz,&rp,&rj,&rv);if(err){free(tmp);return err;}
    rp[0]=0;for(i=0;i<A->n;i++){old=perm[i];len=A->ptr[old+1]-A->ptr[old];for(k=0;k<len;k++){j=A->ptr[old]+k;tmp[k].j=inv[A->index[j]];tmp[k].v=weights[i]*A->value[j];if(colweights)tmp[k].v*=colweights[tmp[k].j];if(!isfinite((double)tmp[k].v)){err=LIS_ERR_ILL_ARG;goto done;}}
        qsort(tmp,(size_t)len,sizeof(*tmp),md_entry_cmp);for(k=0;k<len;k++){rj[off]=tmp[k].j;rv[off++]=tmp[k].v;}rp[i+1]=off;}
    err=lis_matrix_set_csr(A->nnz,rp,rj,rv,P);if(!err){rp=NULL;rj=NULL;rv=NULL;err=lis_matrix_assemble(P);}
done:free(tmp);lis_free(rp);lis_free(rj);lis_free(rv);return err;
}
