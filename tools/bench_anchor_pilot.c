/* Native reference benchmark; the same executable runs fresh and pilot modes. */
#include "lislib.h"
#include "lis_anchor_pilot.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(USE_MPI) || defined(_COMPLEX) || defined(USE_COMPLEX) || defined(_LONG__DOUBLE) || defined(USE_QUAD_PRECISION)
int main(void){fputs("Requires a serial real-double LIS build\n",stderr);return 1;}
#else
static LIS_INT load_matrix(char *name,LIS_MATRIX *out)
{
    LIS_INT rc=lis_matrix_create(LIS_COMM_WORLD,out);
    if(!rc)rc=lis_matrix_set_type(*out,LIS_MATRIX_CSR);
    if(!rc)rc=lis_input_matrix(*out,name);
    if(!rc)rc=lis_matrix_sort_csr(*out);
    return rc;
}
/* Standard Matrix Market real general array, one column. */
static LIS_INT load_rhs(const char *name,LIS_VECTOR b)
{
    FILE *f=fopen(name,"r");char line[512];long long n,m;LIS_INT i;
    if(!f)return LIS_ERR_FILE_IO;
    if(!fgets(line,sizeof(line),f) || strncmp(line,"%%MatrixMarket matrix array real general",40)!=0)goto bad;
    do{if(!fgets(line,sizeof(line),f))goto bad;}while(line[0]=='%');
    if(sscanf(line,"%lld%lld",&n,&m)!=2 || n!=b->n || m!=1)goto bad;
    for(i=0;i<b->n;i++)if(fscanf(f,"%lf",&b->value[i])!=1 || !isfinite(b->value[i]))goto bad;
    fclose(f);return LIS_SUCCESS;
bad:
    fclose(f);return LIS_ERR_FILE_IO;
}
/* Opt-in audit only: exact round-trip decimal doubles, no timing claims. */
static LIS_INT save_solution(const char *dir,int step,LIS_VECTOR x)
{
    char path[4096]; FILE *f; LIS_INT i; int len,failed=0;
    len=snprintf(path,sizeof(path),"%s/x-%03d.mtx",dir,step);
    if(len<0 || (size_t)len>=sizeof(path))return LIS_ERR_FILE_IO;
    f=fopen(path,"wx");
    if(!f)return LIS_ERR_FILE_IO;
    if(fprintf(f,"%%%%MatrixMarket matrix array real general\n%lld 1\n",(long long)x->n)<0)failed=1;
    for(i=0;i<x->n && !failed;i++){
        if(!isfinite(x->value[i]) || fprintf(f,"%.17g\n",x->value[i])<0)failed=1;
    }
    if(fclose(f)!=0)failed=1;
    return failed?LIS_ERR_FILE_IO:LIS_SUCCESS;
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,P=NULL;LIS_VECTOR b=NULL,x=NULL;
    lis_anchor_pilot *ctx=NULL;lis_anchor_pilot_report r;
    LIS_INT rc=0;int i,pilot,low_memory,rolling,accepted=0,creates=0,steps;
    double total=0,setup=0;
    const char *audit_dir=getenv("LIS_ANCHOR_AUDIT_DIR");
    if(argc<4 || (argc%2) || (strcmp(argv[1],"fresh") && strcmp(argv[1],"pilot-low-memory") && strcmp(argv[1],"pilot-rolling") && strcmp(argv[1],"pilot-single-pass"))){
        fprintf(stderr,"usage: %s fresh|pilot-single-pass|pilot-low-memory|pilot-rolling A1.mtx b1.mtx [A2.mtx b2.mtx ...]\n",argv[0]);return 2;
    }
    low_memory=!strcmp(argv[1],"pilot-low-memory");
    rolling=!strcmp(argv[1],"pilot-rolling");
    pilot=rolling||low_memory||!strcmp(argv[1],"pilot-single-pass");steps=(argc-2)/2;
    if(lis_initialize(&argc,&argv))return 1;
    for(i=2;i<argc;i+=2){
        double start,initial=0,elapsed;int first=pilot&&!ctx;
        if((rc=load_matrix(argv[i],&A)) || (rc=lis_vector_duplicate(A,&b)) ||
           (rc=lis_vector_duplicate(A,&x)) || (rc=load_rhs(argv[i+1],b)))goto done;
        start=lis_wtime();
        if(first){P=A;if((rc=rolling?lis_anchor_pilot_create_rolling(P,&ctx):low_memory?lis_anchor_pilot_create_low_memory(P,&ctx):lis_anchor_pilot_create(P,&ctx)))goto done;initial=lis_wtime()-start;}
        rc=pilot?lis_anchor_pilot_solve(ctx,A,b,x,&r):lis_anchor_fresh_solve(A,b,x,&r);
        elapsed=lis_wtime()-start;
        accepted+=r.accepted;creates+=r.factorizations+first;total+=elapsed;setup+=r.setup_seconds+initial;
        printf("STEP step=%d accepted=%d api=%lld raw_status=%lld pilot_iter=%lld continuation_iter=%lld fresh_iter=%lld continued=%d fallback=%d factors=%d ratio=%.12e residual=%.12e setup_s=%.6f total_s=%.6f\n",
            i/2,r.accepted,(long long)rc,(long long)r.raw_status,(long long)r.pilot_iterations,
            (long long)r.continuation_iterations,(long long)r.fresh_iterations,r.continued,r.fallback,
            r.factorizations+first,r.ratio,r.relative_residual,r.setup_seconds+initial,elapsed);
        fflush(stdout);
        if(rc)goto done;
        if(audit_dir && r.accepted && (rc=save_solution(audit_dir,i/2,x)))goto done;
        lis_vector_destroy(x);x=NULL;lis_vector_destroy(b);b=NULL;
        if(rolling){
            LIS_MATRIX next=lis_anchor_pilot_current_anchor(ctx);
            if(next!=P){lis_matrix_destroy(P);P=next;}
        }
        if(A!=P)lis_matrix_destroy(A);
        A=NULL;
    }
    printf("SUMMARY accepted=%d/%d factors=%d setup_s=%.6f total_s=%.6f\n",accepted,steps,creates,setup,total);
done:
    if(ctx)lis_anchor_pilot_destroy(ctx);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A && A!=P)lis_matrix_destroy(A);
    if(P)lis_matrix_destroy(P);
    lis_finalize();
    return rc || accepted!=steps;
}
#endif
