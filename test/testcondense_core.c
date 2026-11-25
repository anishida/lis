#include "lis_condense_core.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int tests=0,errors=0;
#define CHECK(x) do{tests++;if(!(x)){fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#x);errors++;}}while(0)

static void matvec(const lis_sc_csr *a,const double *x,double *b){
    size_t i,k;
    for(i=0;i<a->n;i++){b[i]=0;for(k=a->rowptr[i];k<a->rowptr[i+1];k++)b[i]+=a->values[k]*x[a->colind[k]];}
}
static void csr_from_dense(size_t n,const double *D,lis_sc_csr *a){
    size_t i,j,k=0,nnz=0;
    for(i=0;i<n*n;i++)if(D[i]!=0.0)nnz++;
    a->n=n;a->nnz=nnz;
    size_t *rp=(size_t*)calloc(n+1,sizeof(size_t)),*ci=(size_t*)calloc(nnz+1,sizeof(size_t));
    double *v=(double*)calloc(nnz+1,sizeof(double));
    if(!rp||!ci||!v){fprintf(stderr,"OOM\n");exit(2);}
    for(i=0;i<n;i++){
        rp[i]=k;
        for(j=0;j<n;j++)if(D[i*n+j]!=0.0){ci[k]=j;v[k++]=D[i*n+j];}
    }
    rp[n]=k;
    a->rowptr=rp;a->colind=ci;a->values=v;
}
static void csr_free(lis_sc_csr *a){free((void*)a->rowptr);free((void*)a->colind);free((void*)a->values);memset(a,0,sizeof(*a));}
static double relative(const double *a,const double *b,size_t n){
    double num=0,den=0;size_t i;
    for(i=0;i<n;i++){num=hypot(num,a[i]-b[i]);den=hypot(den,b[i]);}
    return num/(den ? den:1);
}

/* Densely exercise variable-sized local blocks and interleaved DOF numbering. */
static void exercise_mixed(int nonsym, int indef){
    size_t bp[]={0,3,9,18},qd[]={1,3,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36};
    size_t n=43,i,j,k,bl,t,phys[43],np=0;
    int isq[43]={0};
    double D[43*43]={0},x[43]={0},b[43]={0},xr[43]={0},g[43]={0},ur[43]={0},pred[43]={0},Sx[43]={0};
    lis_sc_csr a;lis_sc_plan *p=NULL;
    lis_sc_blocks bs={3,bp,qd};
    lis_sc_status st;
    for(i=0;i<18;i++)isq[qd[i]]=1;
    for(i=0;i<n;i++)if(!isq[i])phys[np++]=i;
    CHECK(np==25);
    for(i=0;i<n;i++){
        x[i]=sin((double)(i+1)*0.37);
        if(!isq[i]){
            D[i*n+i]=6.0;
            for(j=0;j<n;j++)if(!isq[j]&&j!=i&&((i+j)%7==0))
                D[i*n+j]=nonsym ? 0.02*(double)(1+i%3) : 0.05;
        }
    }
    for(bl=0;bl<3;bl++){
        size_t m=bp[bl+1]-bp[bl];
        for(i=0;i<m;i++)for(j=0;j<m;j++){
            size_t di=qd[bp[bl]+i],dj=qd[bp[bl]+j];
            D[di*n+dj]=(i==j)? ((indef&&i==0)?-3.8:3.9+(double)i*.12) : (nonsym ? (0.11+0.02*(double)i):0.11);
        }
        for(t=0;t<4;t++){
            size_t ph=phys[(bl*7+t*3)%np];
            for(k=0;k<m;k++){
                size_t q=qd[bp[bl]+k];
                /* Intentional one-sided edges: K_uq != K_qu^T. */
                D[ph*n+q]=0.04*(double)(k+t+1);
                D[q*n+ph]=nonsym ? (t%2?0.0:0.07*(double)(k+1)) : D[ph*n+q];
            }
        }
    }
    csr_from_dense(n,D,&a);
    st=lis_sc_build(&p,&a,&bs,LIS_SC_FULL,1e-13);
    CHECK(st==LIS_SC_OK);
    if(st==LIS_SC_OK){
        const lis_sc_csr *S=lis_sc_matrix(p);
        CHECK(S->n==25);
        CHECK(lis_sc_internal_dofs(p)==18);
        matvec(&a,x,b);
        st=lis_sc_rhs(p,b,g);CHECK(st==LIS_SC_OK);
        for(i=0;i<25;i++)ur[i]=x[phys[i]];
        matvec(S,ur,Sx);
        CHECK(relative(g,Sx,25)<2e-12);
        st=lis_sc_reconstruct(p,ur,b,xr);CHECK(st==LIS_SC_OK);
        CHECK(relative(xr,x,n)<2e-12);
        double rr=10;
        CHECK(lis_sc_true_residual(p,b,xr,&rr)==LIS_SC_OK);
        CHECK(rr<2e-12);
        /* Reduction preserves labels and arbitrary vectors. */
        long label[43],labred[43];double modes[2*43],mr[2*43];
        for(i=0;i<n;i++){label[i]=(long)(i%4);modes[i]=cos((double)i);modes[n+i]=sin((double)i);}
        CHECK(lis_sc_project_labels(p,label,labred)==LIS_SC_OK);
        CHECK(lis_sc_project_modes(p,2,modes,mr)==LIS_SC_OK);
        for(i=0;i<25;i++){
            CHECK(labred[i]==label[phys[i]]);
            CHECK(mr[i]==modes[phys[i]]);
            CHECK(mr[25+i]==modes[n+phys[i]]);
        }
        /* Different RHS, exactly same operator. */
        for(i=0;i<n;i++)b[i]*=3.2;
        CHECK(lis_sc_rhs(p,b,g)==LIS_SC_OK);
        for(i=0;i<25;i++)ur[i]*=3.2;
        CHECK(lis_sc_reconstruct(p,ur,b,xr)==LIS_SC_OK);
        for(i=0;i<n;i++)pred[i]=x[i]*3.2;
        CHECK(relative(xr,pred,n)<2e-12);
        lis_sc_plan *original=p;
        CHECK(lis_sc_build(&p,&a,&bs,LIS_SC_REUSE_ALL,1e-13)==LIS_SC_OK && p==original);
        ((double*)a.values)[0]+=0.25;
        CHECK(lis_sc_build(&p,&a,&bs,LIS_SC_REUSE_ALL,1e-13)==LIS_SC_PATTERN_MISMATCH && p==original);
        CHECK(lis_sc_build(&p,&a,&bs,LIS_SC_NUMERIC,1e-13)==LIS_SC_OK && p!=NULL);
        ((double*)a.values)[0]-=0.25;
    }
    lis_sc_free(p);csr_free(&a);
}

/* Regression: K_uq has a structural edge absent from K_qu. */
static void test_one_way(void){
    double D[]={4,0,2,  0,4,0,  0,3,2};
    size_t ptr[]={0,1},dof[]={2};
    lis_sc_csr A;lis_sc_plan *p=NULL;lis_sc_blocks B={1,ptr,dof};
    csr_from_dense(3,D,&A);
    CHECK(lis_sc_build(&p,&A,&B,LIS_SC_FULL,1e-14)==LIS_SC_OK);
    if(p){
        const lis_sc_csr *S=lis_sc_matrix(p);
        double u[]={1.3,2.1},full[]={1.3,2.1,0},b[3],g[2],sg[2];
        double expected[]={4,-3,0,4};
        for(size_t i=0;i<2;i++)for(size_t j=0;j<2;j++){
            int found=0;
            for(size_t k=S->rowptr[i];k<S->rowptr[i+1];k++)if(S->colind[k]==j){
                CHECK(fabs(S->values[k]-expected[2*i+j])<1e-14);found=1;break;
            }
            CHECK(found || expected[2*i+j]==0);
        }
        matvec(&A,full,b);
        CHECK(lis_sc_rhs(p,b,g)==LIS_SC_OK);
        matvec(S,u,sg); CHECK(relative(sg,g,2)<1e-13);
    }
    lis_sc_free(p);csr_free(&A);
}
static void test_atomic_errors(void){
    double D[]={4,1,0,  1,3,1,  0,1,2};
    size_t ptr[]={0,1},qd[]={2};lis_sc_csr A;
    lis_sc_blocks B={1,ptr,qd};lis_sc_plan *p=NULL;lis_sc_status st;
    csr_from_dense(3,D,&A);
    CHECK(lis_sc_build(&p,&A,&B,LIS_SC_FULL,1e-12)==LIS_SC_OK);
    lis_sc_plan *before=p;
    size_t wrong[]={9};lis_sc_blocks BAD={1,ptr,wrong};
    st=lis_sc_build(&p,&A,&BAD,LIS_SC_FULL,1e-12);
    CHECK(st==LIS_SC_BAD_MAP && p==before);
    ((double*)A.values)[A.nnz-1]=0.0;
    st=lis_sc_build(&p,&A,&B,LIS_SC_FULL,1e-12);
    /* Last diagonal is zero while q->u remains: q block singular. */
    CHECK(st==LIS_SC_SINGULAR_BLOCK && p==before);
    ((double*)A.values)[A.nnz-1]=2.0;
    lis_sc_free(p);csr_free(&A);
    /* Two separate q-blocks are not independent. */
    double E[]={3,0,0,0, 0,2,0,0, 0,0,4,2, 0,0,2,4};
    size_t ptr2[]={0,1,2},q2[]={2,3};lis_sc_blocks C={2,ptr2,q2};
    csr_from_dense(4,E,&A);p=NULL;
    CHECK(lis_sc_build(&p,&A,&C,LIS_SC_FULL,1e-12)==LIS_SC_COUPLED_BLOCKS);
    CHECK(p==NULL);
    csr_free(&A);
    /* unsorted CSR and nonfinite data rejected before allocations */
    size_t rp[]={0,2,3},ci[]={1,0,1};double vv[]={1,2,3};
    lis_sc_csr badcsr={2,3,rp,ci,vv};lis_sc_blocks Z={0,NULL,NULL};
    CHECK(lis_sc_build(&p,&badcsr,&Z,LIS_SC_FULL,1e-12)==LIS_SC_BAD_CSR);
    ci[0]=0;ci[1]=1;vv[0]=NAN;
    CHECK(lis_sc_build(&p,&badcsr,&Z,LIS_SC_FULL,1e-12)==LIS_SC_NONFINITE);
}
static void test_no_internal(void){
    double D[]={4,1, 2,3};lis_sc_csr A;lis_sc_plan *p=NULL;
    lis_sc_blocks Z={0,NULL,NULL};
    csr_from_dense(2,D,&A);
    CHECK(lis_sc_build(&p,&A,&Z,LIS_SC_FULL,1e-13)==LIS_SC_OK);
    double x[]={1.0,3.0},b[2],red[2],xr[2];
    matvec(&A,x,b);
    CHECK(lis_sc_rhs(p,b,red)==LIS_SC_OK);
    CHECK(relative(red,b,2)==0);
    CHECK(lis_sc_reconstruct(p,x,b,xr)==LIS_SC_OK);
    CHECK(relative(xr,x,2)==0);
    /* Projected modes and RHS remain unchanged if invalid input is detected. */
    double full_modes[]={1.0,2.0,3.0,NAN},reduced_modes[]={-91,-91,-91,-91};
    CHECK(lis_sc_project_modes(p,2,full_modes,reduced_modes)==LIS_SC_NONFINITE);
    CHECK(reduced_modes[0]==-91 && reduced_modes[1]==-91 &&
          reduced_modes[2]==-91 && reduced_modes[3]==-91);
    CHECK(lis_sc_project_modes(p,SIZE_MAX,full_modes,reduced_modes)==LIS_SC_SIZE);
    full_modes[3]=4.0;
    CHECK(lis_sc_project_modes(p,2,full_modes,reduced_modes)==LIS_SC_OK);
    CHECK(reduced_modes[0]==1.0 && reduced_modes[3]==4.0);
    double bad_rhs[]={1.0,NAN},unchanged[]={-70,-70};
    CHECK(lis_sc_rhs(p,bad_rhs,unchanged)==LIS_SC_NONFINITE);
    CHECK(unchanged[0]==-70 && unchanged[1]==-70);
    lis_sc_free(p);csr_free(&A);
}

/* A stored zero may become a live coupling while the CSR pattern stays
 * unchanged. Numeric reuse must preserve the structural Schur fill, and a
 * failed numerical update must leave the prior complete plan untouched. */
static void test_numeric_pattern_reuse(void){
    enum { N=6 };
    size_t rp[N+1],ci[N*N],bp[]={0,1,3},qd[]={3,4,5};
    double v[N*N],x[N],rhs[N],g[N],xr[N],u[3],su[3];
    lis_sc_csr A={N,N*N,rp,ci,v};
    lis_sc_blocks B={2,bp,qd};
    lis_sc_plan *p=NULL,*reference=NULL;
    size_t k=0,i,j;
    for(i=0;i<N;i++){
        rp[i]=k;
        for(j=0;j<N;j++){
            ci[k]=j;
            v[k]=(i==j) ? 5.0+0.1*i : 0.0;
            k++;
        }
        x[i]=sin(0.23*(i+1));
    }
    rp[N]=k;
    /* Initial case: three independent physical DOFs. */
    CHECK(lis_sc_build(&p,&A,&B,LIS_SC_FULL,1e-13)==LIS_SC_OK);
    if(!p)return;
    double old_values[9]={0};
    const lis_sc_csr *initial=lis_sc_matrix(p);
    for(i=0;i<initial->nnz;i++){
        size_t row=0;
        while(row+1<initial->n && i>=initial->rowptr[row+1])row++;
        old_values[row*3+initial->colind[i]]=initial->values[i];
    }
    /* Activate previously zero, structurally existing one-way edges.
     * q3 belongs to block 0; q4,q5 belong to block 1. */
    v[0*N+3]=0.3;
    v[3*N+1]=0.12;
    v[2*N+4]=-0.23;
    v[4*N+0]=0.38;
    v[2*N+5]=-0.08;
    v[5*N+1]=0.05;
    v[4*N+5]=0.17;
    v[5*N+4]=0.11;
    CHECK(lis_sc_build(&reference,&A,&B,LIS_SC_FULL,1e-13)==LIS_SC_OK);
    CHECK(lis_sc_build(&p,&A,&B,LIS_SC_NUMERIC,1e-13)==LIS_SC_OK);
    if(reference && p){
        const lis_sc_csr *R=lis_sc_matrix(reference),*S=lis_sc_matrix(p);
        CHECK(S->n==R->n && S->nnz==R->nnz);
        CHECK(memcmp(S->rowptr,R->rowptr,(S->n+1)*sizeof(size_t))==0);
        CHECK(memcmp(S->colind,R->colind,S->nnz*sizeof(size_t))==0);
        for(k=0;k<S->nnz;k++)CHECK(fabs(S->values[k]-R->values[k])<1e-14);
        matvec(&A,x,rhs);
        CHECK(lis_sc_rhs(p,rhs,g)==LIS_SC_OK);
        const size_t *map=lis_sc_physical_to_full(p);
        for(i=0;i<3;i++)u[i]=x[map[i]];
        matvec(S,u,su);
        CHECK(relative(su,g,3)<1e-13);
        CHECK(lis_sc_reconstruct(p,u,rhs,xr)==LIS_SC_OK);
        CHECK(relative(xr,x,N)<1e-13);
        double rr=0;
        CHECK(lis_sc_true_residual(p,rhs,xr,&rr)==LIS_SC_OK && rr<1e-13);
        /* A second numerical update modifies local block and physical
         * coupling values while preserving the exact CSR structure. */
        v[0*N+3]*=0.5;
        v[4*N+5]*=1.5;
        CHECK(lis_sc_build(&p,&A,&B,LIS_SC_NUMERIC,1e-13)==LIS_SC_OK);
        lis_sc_free(reference);reference=NULL;
        CHECK(lis_sc_build(&reference,&A,&B,LIS_SC_FULL,1e-13)==LIS_SC_OK);
        if(reference){
            R=lis_sc_matrix(reference);S=lis_sc_matrix(p);
            CHECK(S->nnz==R->nnz);
            for(k=0;k<S->nnz;k++)CHECK(fabs(S->values[k]-R->values[k])<1e-14);
        }
        /* Coupling between distinct internal blocks becomes nonzero, even
         * though the pattern was always dense.  Refuse and roll back. */
        lis_sc_plan *kept=p;
        S=lis_sc_matrix(p);
        double *snapshot=(double*)malloc(S->nnz*sizeof(double));
        if(!snapshot){fprintf(stderr,"OOM\n");exit(2);}
        memcpy(snapshot,S->values,S->nnz*sizeof(double));
        v[3*N+4]=0.07;
        CHECK(lis_sc_build(&p,&A,&B,LIS_SC_NUMERIC,1e-13)==LIS_SC_COUPLED_BLOCKS);
        CHECK(p==kept);
        CHECK(memcmp(snapshot,lis_sc_matrix(p)->values,S->nnz*sizeof(double))==0);
        free(snapshot);
        v[3*N+4]=0.0;
        CHECK(lis_sc_build(&p,&A,&B,LIS_SC_REUSE_ALL,1e-13)==LIS_SC_OK);
        v[1*N+5]=0.01;
        CHECK(lis_sc_build(&p,&A,&B,LIS_SC_REUSE_ALL,1e-13)==LIS_SC_PATTERN_MISMATCH);
        CHECK(lis_sc_build(&p,&A,&B,LIS_SC_NUMERIC,1e-13)==LIS_SC_OK);
    }
    /* Original snapshot is not modified by any numeric rebuild. */
    for(i=0;i<3;i++)for(j=0;j<3;j++)if(i==j)
        CHECK(old_values[3*i+j]>0);
    lis_sc_free(reference);lis_sc_free(p);
}

static unsigned rndstate=0xC0FFEEu;
static unsigned rng(void){rndstate=rndstate*1664525u+1013904223u;return rndstate;}
static void randomized_properties(void){
    for(size_t seed=0;seed<100;seed++){
        size_t n=8+(rng()%19),nb=1+(rng()%3),nq=nb+(rng()%3);
        size_t ptr[4]={0},qd[12],used[32]={0},next=0;
        double D[32*32]={0},x[32],b[32],g[32],u[32],recovered[32],Sx[32];
        lis_sc_csr A;lis_sc_plan *P=NULL;
        for(size_t k=0;k<nq;k++){
            size_t q=(k*7+seed*5)%n;
            while(used[q])q=(q+1)%n;
            used[q]=1;qd[k]=q;
        }
        for(size_t k=0;k<nb;k++){
            size_t count=nq/nb+(k<nq%nb);
            ptr[k]=next;next+=count;
        }
        ptr[nb]=nq;
        for(size_t i=0;i<n;i++){
            x[i]=((double)(rng()%2001)-1000.0)/1001.0;
            D[i*n+i]=used[i] ? 2.4 : 4.1;
        }
        for(size_t i=0;i<n;i++)for(size_t j=0;j<n;j++)if(i!=j){
            int bi=-1,bj=-1;
            for(size_t k=0;k<nq;k++){
                if(qd[k]==i)for(size_t r=0;r<nb;r++)if(k>=ptr[r]&&k<ptr[r+1])bi=(int)r;
                if(qd[k]==j)for(size_t r=0;r<nb;r++)if(k>=ptr[r]&&k<ptr[r+1])bj=(int)r;
            }
            if(bi>=0&&bj>=0&&bi!=bj)continue;
            if((rng()%5)!=0)continue;
            D[i*n+j]=((double)(rng()%101)-50.0)/680.0;
        }
        if(seed%11==0){
            /* Negative local pivot while retaining nonsingularity. */
            size_t q=qd[0];D[q*n+q]=-3.0;
        }
        csr_from_dense(n,D,&A);
        lis_sc_blocks B={nb,ptr,qd};
        lis_sc_status st=lis_sc_build(&P,&A,&B,LIS_SC_FULL,1e-13);
        CHECK(st==LIS_SC_OK);
        if(st==LIS_SC_OK){
            const lis_sc_csr *S=lis_sc_matrix(P);
            matvec(&A,x,b);
            CHECK(lis_sc_rhs(P,b,g)==LIS_SC_OK);
            const size_t *map=lis_sc_physical_to_full(P);
            for(size_t i=0;i<S->n;i++)u[i]=x[map[i]];
            matvec(S,u,Sx);
            CHECK(relative(Sx,g,S->n)<4e-12);
            CHECK(lis_sc_reconstruct(P,u,b,recovered)==LIS_SC_OK);
            CHECK(relative(recovered,x,n)<4e-12);
            double rr=0;
            CHECK(lis_sc_true_residual(P,b,recovered,&rr)==LIS_SC_OK && rr<4e-12);
        }
        lis_sc_free(P);csr_free(&A);
    }
}

int main(void){
    randomized_properties();
    test_numeric_pattern_reuse();
    exercise_mixed(0,0);exercise_mixed(1,0);exercise_mixed(1,1);
    test_one_way();test_atomic_errors();test_no_internal();
    printf("condensation core: %d assertions, %d failures\n",tests,errors);
    return errors?1:0;
}
