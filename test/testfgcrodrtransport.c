/* Regression for FGCRO-DR cross-solve TRANSPORT policy.
 * Tests a nonuniform value change in the same assembled CSR matrix object
 * with the same sparsity graph and a reused LIS_SOLVER.
 * Asserts policy counters, matrix generations, AZ=C and true residual.
 * No MFAMG, CalculiX, or private data are needed.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include <math.h>
#include "lis.h"
#include "lis_solver.h"
#include "lislib.h"

#ifdef _COMPLEX
int main(int argc,char **argv)
{
    (void)argc; (void)argv;
    puts("SKIP: real-only transport regression");
    return 77;
}
#else
#define N 120

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,stored=NULL;
    LIS_VECTOR exact=NULL,b=NULL,x=NULL,ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,p,is=0,ie=0,failed=1,dim=0,valid=0;
    LIS_UNSIGNED_INT uses0=0,keeps0=0,transports0=0,rejects0=0;
    LIS_UNSIGNED_INT uses1=0,keeps1=0,transports1=0,rejects1=0;
    LIS_UNSIGNED_INT sg0=0,vg0=0,sg1=0,vg1=0,ssg=0,svg=0;
    LIS_REAL orth=0.0,relation=0.0,rnorm=0.0,bnorm=0.0;

    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    if((err=lis_matrix_create(LIS_COMM_WORLD,&A))) goto done;
    if((err=lis_matrix_set_size(A,0,N))) goto done;
    if((err=lis_matrix_get_range(A,&is,&ie))) goto done;
    for(i=is;i<ie;i++)
    {
        if((err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A)))
            goto done;
        if(i>0 && (err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A)))
            goto done;
        if(i<N-1 && (err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.5,A)))
            goto done;
    }
    if((err=lis_matrix_set_type(A,LIS_MATRIX_CSR))) goto done;
    if((err=lis_matrix_assemble(A))) goto done;
    if((err=lis_vector_duplicate(A,&exact))) goto done;
    if((err=lis_vector_duplicate(A,&b))) goto done;
    if((err=lis_vector_duplicate(A,&x))) goto done;
    if((err=lis_vector_duplicate(A,&ax))) goto done;
    for(i=0;i<exact->n;i++)
        exact->value[i]=(LIS_SCALAR)((((is+i)*37)%23)-11.0);
    if((err=lis_matvec(A,exact,b))) goto done;
    if((err=lis_solver_create(&solver))) goto done;

    /* KEEP=1 requires near-perfect alignment. REJECT=0 never rejects
     * nonnegative quality; a genuine nonparallel update must TRANSPORT. */
    if((err=lis_solver_set_option(
        "-i fgcrodr -p jacobi -scale none -restart 2 -maxiter 600 "
        "-tol 1e-5 -fgcrodr_k 1 -fgcrodr_keep 1 "
        "-fgcrodr_reject 0 -fgcrodr_intra 1 "
        "-fgcrodr_projected_relative 1 -fgcrodr_intra_trace 0 "
        "-print none",solver))) goto done;
    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_policy_stats(
        solver,&uses0,&keeps0,&transports0,&rejects0))) goto done;
    if((err=lis_solver_get_fgcrodr_recycle(
        solver,&dim,NULL,NULL,NULL))) goto done;
    if(dim<1) {err=LIS_ERR_ILL_ARG; goto done;}
    if((err=lis_matrix_get_generation(A,&sg0,&vg0))) goto done;

    /* Nonuniform diagonal change, CSR graph/pointer/column order unchanged.
     * Diagonal remains strictly dominant (min 2.6 vs 1.5 off-diagonal). */
    if(A->matrix_type!=LIS_MATRIX_CSR || A->ptr==NULL ||
       A->index==NULL || A->value==NULL)
    {
        err=LIS_ERR_ILL_ARG; goto done;
    }
    for(i=0;i<A->n;i++)
    {
        LIS_INT global_row=is+i;
        int found=0;
        for(p=A->ptr[i];p<A->ptr[i+1];p++)
        {
            if(A->index[p]==i)
            {
                A->value[p]+=(LIS_SCALAR)(0.35*(((global_row*23)%9)-4));
                found=1;
                break;
            }
        }
        if(!found) {err=LIS_ERR_ILL_ARG; goto done;}
    }
    if((err=lis_matrix_mark_values_changed(A))) goto done;
    if((err=lis_matrix_get_generation(A,&sg1,&vg1))) goto done;
    if(sg0!=sg1 || vg1<=vg0) {err=LIS_ERR_ILL_ARG; goto done;}
    if((err=lis_matvec(A,exact,b))) goto done;
    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_policy_stats(
        solver,&uses1,&keeps1,&transports1,&rejects1))) goto done;
    if((err=lis_fgcrodr_get_lifecycle_state(
        solver,&valid,&stored,&ssg,&svg))) goto done;
    if(!valid || stored!=A || ssg!=sg1 || svg!=vg1 ||
       transports1<=transports0 || uses1<=uses0 ||
       keeps1!=keeps0 || rejects1!=rejects0)
    {
        err=LIS_ERR_ILL_ARG;
        goto done;
    }
    if((err=lis_fgcrodr_validate_recycle(
        solver,&orth,&relation))) goto done;
    if((err=lis_matvec(A,x,ax))) goto done;
    if((err=lis_vector_xpay(b,-1.0,ax))) goto done;
    if((err=lis_vector_nrm2(ax,&rnorm))) goto done;
    if((err=lis_vector_nrm2(b,&bnorm))) goto done;
    if(!isfinite(rnorm) || bnorm<=0.0 || rnorm/bnorm>1e-5 ||
       !isfinite(orth) || !isfinite(relation) ||
       orth>1e-8 || relation>1e-8)
    {
        err=LIS_ERR_ILL_ARG;
        goto done;
    }
    printf("FGCRODR_TRANSPORT_POLICY_OK cross_uses=%llu->%llu "
           "transports=%llu->%llu relative=%.12e orth=%.6e relation=%.6e\n",
           (unsigned long long)uses0,(unsigned long long)uses1,
           (unsigned long long)transports0,
           (unsigned long long)transports1,
           (double)(rnorm/bnorm),(double)orth,(double)relation);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_TRANSPORT_POLICY_FAIL api=%d "
                "status=%d dim=%d trans=%llu->%llu "
                "keeps=%llu->%llu rejects=%llu->%llu "
                "use=%llu->%llu orth=%.6e relation=%.6e\n",
                (int)err,solver?(int)solver->retcode:-1,(int)dim,
                (unsigned long long)transports0,
                (unsigned long long)transports1,
                (unsigned long long)keeps0,
                (unsigned long long)keeps1,
                (unsigned long long)rejects0,
                (unsigned long long)rejects1,
                (unsigned long long)uses0,
                (unsigned long long)uses1,
                (double)orth,(double)relation);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
#endif
