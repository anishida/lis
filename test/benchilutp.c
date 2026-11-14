/*
 * Native ILUT versus ILUTP benchmark for original CCX Matrix Market A/b.
 * One process, one OpenMP thread. Never fabricate or replace the given RHS.
 *
 * Usage: ./test/benchilutp A.mtx b.mtx [--maxiter N] [--restart N]
 *        [--rtol V] [--drop V] [--rate V] [--pivot-tol V]
 *        [--rcm 0|1] [--row-scale 0|1] [--equil 0|1]
 *
 * Emits CSV, comparing ILUT against ILUTP with tau=0.5 and 1.0.
 * Factor setup and Krylov CPU times are measured separately; the original
 * relative residual ||b-A*x||_2/||b||_2 is always computed.
 *
 * This intentionally does not use lis_input_vector(): LIS's MM reader has
 * historical vector-format expectations that differ from the CCX
 * '%%MatrixMarket matrix array real general' exports.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <errno.h>
#include <time.h>

#if defined(_COMPLEX) || defined(USE_COMPLEX)
int main(void)
{
    fprintf(stderr,"benchilutp supports real-double LIS only\n");
    return 77;
}
#else
typedef struct
{
    LIS_INT maxiter,restart,rcm,row_scale,equil,full,only_full,full_limit;
    double rtol,drop,rate,pivot_tol,full_tau;
} BENCH_ARGS;

static double cpu_seconds(void)
{
    return (double)clock()/(double)CLOCKS_PER_SEC;
}

/* Read the original array-style Matrix Market RHS exported by CCX. */
static int read_ccx_rhs(const char *filename,LIS_VECTOR b,LIS_INT n)
{
    FILE *f;
    char line[4096],mm[128],object[128],fmt[128],type[128],shape[128];
    long rows,cols,i;
    double value;
    f=fopen(filename,"r");
    if(!f)
    {
        fprintf(stderr,"Cannot open RHS: %s\n",filename);
        return 1;
    }
    if(!fgets(line,sizeof(line),f) ||
       sscanf(line,"%127s %127s %127s %127s %127s",
              mm,object,fmt,type,shape)!=5 ||
       strcmp(mm,"%%MatrixMarket")!=0 ||
       strcmp(object,"matrix")!=0 ||
       strcmp(fmt,"array")!=0 ||
       strcmp(type,"real")!=0 ||
       strcmp(shape,"general")!=0)
    {
        fprintf(stderr,"Expected Matrix Market array real general RHS\n");
        fclose(f);
        return 1;
    }
    do
    {
        if(!fgets(line,sizeof(line),f))
        {
            fclose(f);
            return 1;
        }
    }while(line[0]=='%');
    if(sscanf(line,"%ld %ld",&rows,&cols)!=2 ||
       rows!=(long)n || cols!=1)
    {
        fprintf(stderr,"RHS shape must be %ld x 1\n",(long)n);
        fclose(f);
        return 1;
    }
    for(i=0;i<rows;i++)
    {
        if(!fgets(line,sizeof(line),f) ||
           sscanf(line,"%lf",&value)!=1 || !isfinite(value))
        {
            fprintf(stderr,"Invalid RHS at row %ld\n",i);
            fclose(f);
            return 1;
        }
        b->value[i]=(LIS_SCALAR)value;
    }
    fclose(f);
    return 0;
}

static int run_case(FILE *csv,const char *name,LIS_MATRIX A,LIS_VECTOR b,
                    LIS_VECTOR x,LIS_VECTOR resid,double norm_b,
                    const BENCH_ARGS *args,int use_ilutp,double tau,int full)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err=0,i,iter=0,retcode=0;
    long long nnz_l=0,nnz_u=0,permuted=0;
    char opts[1024];
    double setup_s=0.0,solve_s=0.0,t,relative=-1.0,rnorm;
    const char *status="internal_error";

    err=lis_solver_create(&solver);
    if(err || !solver)goto report;
    if(snprintf(opts,sizeof(opts),
        "-i fgmres -p ilut -scale none -print none "
        "-initx_zeros true -restart %lld -maxiter %lld "
        "-tol %.16g -ilut_drop %.16g -ilut_rate %.16g "
        "-ilu_pivot_tol %.16g -omp_num_threads 1 "
        "-ilut_rcm %lld -ilut_row_scale %lld -ilut_equil %lld "
        "-ilutp %d -ilutp_threshold %.16g "
        "-ilutp_full %d -ilutp_full_limit %lld",
        (long long)args->restart,(long long)args->maxiter,args->rtol,
        args->drop,args->rate,args->pivot_tol,
        (long long)args->rcm,(long long)args->row_scale,
        (long long)args->equil,use_ilutp,tau,full,
        full ? (long long)args->full_limit : 0LL)>=(int)sizeof(opts))
    {
        err=LIS_ERR_ILL_ARG;
        goto report;
    }
    err=lis_solver_set_option(opts,solver);
    if(err)goto report;
    err=lis_solver_set_sequence_policy(solver,LIS_SEQUENCE_POLICY_AUTO);
    if(err)goto report;
    t=cpu_seconds();
    err=lis_solver_sequence_prepare_precon(A,solver,&pc);
    setup_s=cpu_seconds()-t;
    if(err || !pc)
    {
        status="setup_failed";
        goto report;
    }
    if(pc->L && pc->U)
    {
        for(i=0;i<A->n;i++)
        {
            nnz_l+=(long long)pc->L->nnz[i];
            nnz_u+=(long long)pc->U->nnz[i];
        }
    }
    if(pc->ilutp_colperm)
        for(i=0;i<A->n;i++)
            if(pc->ilutp_colperm[i]!=i)permuted++;

    err=lis_vector_set_all((LIS_SCALAR)0.0,x);
    if(err)goto report;
    t=cpu_seconds();
    /* Factor setup was timed above. Borrow exactly those factors for one
     * Krylov run: lis_solve_sequence may retry/rebuild after breakdown and
     * would otherwise attribute another factorization to solve_cpu_s. */
    err=lis_solve_kernel(A,b,x,solver,pc);
    solve_s=cpu_seconds()-t;
    retcode=solver->retcode;
    iter=solver->iter;
    if(!err)
    {
        LIS_INT v_err=lis_matvec(A,x,resid);
        if(!v_err)v_err=lis_vector_axpy((LIS_SCALAR)-1.0,b,resid);
        if(!v_err)v_err=lis_vector_nrm2(resid,&rnorm);
        if(v_err)err=v_err;
        else relative=rnorm/norm_b;
    }
    status=err?"api_failed":retcode==LIS_SUCCESS?"converged":"not_converged";

report:
    if(solver && iter==0)iter=solver->iter;
    fprintf(csv,"%s,%d,%lld,%.6g,%s,%lld,%lld,%lld,%lld,%.8g,%.8g,%.15g,%d,%lld\n",
           name,use_ilutp,(long long)args->equil,tau,status,(long long)iter,
           nnz_l,nnz_u,permuted,setup_s,solve_s,relative,
           (int)err,(long long)retcode);
    fflush(csv);
    if(solver)lis_solver_destroy(solver);
    return 0; /* Solver failure is measured, not a benchmark crash. */
}

/* Unlike floating conversion, this rejects fractional/overflow limits. */
static int parse_nonnegative_lis_int(const char *s,LIS_INT *out)
{
    char *end=NULL;
    unsigned long long v,maxval=(sizeof(LIS_INT)==4 ? INT32_MAX : INT64_MAX);
    if(!s || !*s || *s=='-' || *s=='+')return 1;
    errno=0;
    v=strtoull(s,&end,10);
    if(errno==ERANGE || end==s || *end || v>maxval)return 1;
    *out=(LIS_INT)v;
    return 0;
}

static int parse_num(const char *s,double *value)
{
    char *end=NULL;
    *value=strtod(s,&end);
    return end==s || *end!='\0' || !isfinite(*value);
}

int main(int argc,char **argv)
{
    BENCH_ARGS a;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,resid=NULL;
    LIS_REAL bnorm;
    double t;
    FILE *csv=stdout;
    const char *csv_filename=NULL;
    int err=1,i;
    LIS_INT rc;
    if(argc<3)
    {
        fprintf(stderr,
            "Usage: %s A.mtx b.mtx [--maxiter N] [--restart N] "
            "[--rtol V] [--drop V] [--rate V] [--pivot-tol V] "
            "[--rcm 0|1] [--row-scale 0|1] [--equil 0|1] "
            "[--full 0|1] [--only-full 0|1] [--full-limit N] [--full-tau 0..1] "
            "[--csv FILE]\n",argv[0]);
        return 2;
    }
    a.maxiter=1000;a.restart=120;
    a.rtol=1e-5;a.drop=0.005;a.rate=12.0;
    a.pivot_tol=1e-9;a.rcm=1;a.row_scale=1;a.equil=0;
    a.full=0;a.only_full=0;a.full_limit=0;a.full_tau=1.0;
    for(i=3;i<argc;i+=2)
    {
        double num;
        if(i+1==argc)
        {
            fprintf(stderr,"Missing argument for %s\n",argv[i]);
            return 2;
        }
        if(strcmp(argv[i],"--csv")==0)
        {
            csv_filename=argv[i+1];
            continue;
        }
        if(strcmp(argv[i],"--full-limit")==0)
        {
            if(parse_nonnegative_lis_int(argv[i+1],&a.full_limit))
            {
                fprintf(stderr,"Invalid nonnegative LIS_INT fill limit: %s\n",argv[i+1]);
                return 2;
            }
            continue;
        }
        /* Convert integer options as integers, not via double: a large
         * finite double may overflow LIS_INT or silently truncate. */
        if(strcmp(argv[i],"--maxiter")==0 ||
           strcmp(argv[i],"--restart")==0 ||
           strcmp(argv[i],"--rcm")==0 ||
           strcmp(argv[i],"--row-scale")==0 ||
           strcmp(argv[i],"--equil")==0 ||
           strcmp(argv[i],"--full")==0 ||
           strcmp(argv[i],"--only-full")==0)
        {
            LIS_INT value;
            if(parse_nonnegative_lis_int(argv[i+1],&value) ||
               ((strcmp(argv[i],"--rcm")==0 ||
                 strcmp(argv[i],"--row-scale")==0 ||
                 strcmp(argv[i],"--equil")==0 ||
                 strcmp(argv[i],"--full")==0 ||
                 strcmp(argv[i],"--only-full")==0) && value>1))
            {
                fprintf(stderr,"Invalid integer argument for %s: %s\n",
                        argv[i],argv[i+1]);
                return 2;
            }
            if(strcmp(argv[i],"--maxiter")==0)a.maxiter=value;
            else if(strcmp(argv[i],"--restart")==0)a.restart=value;
            else if(strcmp(argv[i],"--rcm")==0)a.rcm=value;
            else if(strcmp(argv[i],"--row-scale")==0)a.row_scale=value;
            else if(strcmp(argv[i],"--equil")==0)a.equil=value;
            else if(strcmp(argv[i],"--full")==0)a.full=value;
            else a.only_full=value;
            continue;
        }
        if(parse_num(argv[i+1],&num))
        {
            fprintf(stderr,"Invalid argument near %s\n",argv[i]);
            return 2;
        }
        if(strcmp(argv[i],"--rtol")==0)a.rtol=num;
        else if(strcmp(argv[i],"--drop")==0)a.drop=num;
        else if(strcmp(argv[i],"--rate")==0)a.rate=num;
        else if(strcmp(argv[i],"--pivot-tol")==0)a.pivot_tol=num;
        else if(strcmp(argv[i],"--full-tau")==0)a.full_tau=num;
        else
        {
            fprintf(stderr,"Unknown option %s\n",argv[i]);
            return 2;
        }
    }
    if(a.maxiter<1 || a.restart<1 || a.rtol<=0 || a.drop<0 ||
       a.rate<0 || a.pivot_tol<0 ||
       (a.rcm!=0 && a.rcm!=1) ||
       (a.row_scale!=0 && a.row_scale!=1) ||
       (a.equil!=0 && a.equil!=1) ||
       (a.only_full && !a.full) ||
       a.full_tau<0.0 || a.full_tau>1.0)
        return 2;

    rc=lis_initialize(&argc,&argv);
    if(rc)return 1;
    rc=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(rc)goto cleanup;
    rc=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(rc)goto cleanup;
    t=cpu_seconds();
    rc=lis_input_matrix(A,argv[1]);
    if(rc)goto cleanup;
    fprintf(stderr,"BENCH_ILUTP_INPUT n=%lld nnz=%lld load_cpu=%.6g\n",
        (long long)A->n,(long long)A->nnz,cpu_seconds()-t);
    if(A->n<1 || A->nprocs!=1 || A->n!=A->gn)
    {
        fprintf(stderr,"Benchmark requires one MPI rank and serial CSR\n");
        goto cleanup;
    }
    rc=lis_vector_duplicate(A,&b);
    if(rc)goto cleanup;
    rc=lis_vector_duplicate(A,&x);
    if(rc)goto cleanup;
    rc=lis_vector_duplicate(A,&resid);
    if(rc)goto cleanup;
    if(read_ccx_rhs(argv[2],b,A->n))goto cleanup;
    rc=lis_vector_nrm2(b,&bnorm);
    if(rc || !(bnorm>0))goto cleanup;

    if(csv_filename)
    {
        csv=fopen(csv_filename,"w");
        if(!csv)
        {
            fprintf(stderr,"Cannot open CSV output: %s\n",csv_filename);
            goto cleanup;
        }
    }
    fprintf(csv,"variant,ilutp,equil,tau,status,iterations,nnz_L,nnz_U,"
           "moved_columns,setup_cpu_s,solve_cpu_s,physical_relres,"
           "api_error,solver_status\n");
    if(!a.only_full)
    {
        run_case(csv,"legacy_ilut",A,b,x,resid,(double)bnorm,&a,0,0,0);
        run_case(csv,"ilutp_0p5",A,b,x,resid,(double)bnorm,&a,1,0.5,0);
        run_case(csv,"ilutp_1p0",A,b,x,resid,(double)bnorm,&a,1,1.0,0);
    }
    if(a.full)
        run_case(csv,"ilutp_full",A,b,x,resid,(double)bnorm,&a,1,a.full_tau,1);
    err=0;
cleanup:
    if(csv && csv!=stdout)fclose(csv);
    if(resid)lis_vector_destroy(resid);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return err;
}

#endif /* real-double only */
