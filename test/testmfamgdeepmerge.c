/* MFAMG deep sibling regression:
 * (1) standalone star grouping, (2) actual multilevel GMRES/MFAMG solves,
 * (3) option transitions and rank-local MPI fallback verification.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <math.h>
#include <stdio.h>
#include "lislib.h"
#include "lis_solver_options.h"

#define STAR_N 121
#define STAR_ARMS 32
#define ARM_LENGTH 32
#define DEEP_N (1 + STAR_ARMS*ARM_LENGTH)
#define MPI_N 2049
#define MAX_LEVELS 32

extern LIS_INT lis_precon_mfamg_merge_siblings_for_test(
    LIS_MATRIX A, LIS_INT target, LIS_INT *aggregate, LIS_INT *nc);
extern LIS_INT lis_precon_mfamg_get_hierarchy_info(
    LIS_PRECON precon, LIS_INT max_levels, LIS_INT *nlevels,
    LIS_INT *coarse_n, LIS_INT *level_n, LIS_INT *owns_A);
extern LIS_INT lis_precon_mfamg_get_deep_merge_info(
    LIS_PRECON precon, LIS_INT *attempts, LIS_INT *events,
    LIS_INT *reduced, LIS_INT *max_depth);

/* This structural case deliberately has 120 equally weighted leaves. */
static int structural_star_test(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT id[STAR_N],sizes[STAR_N],nc=STAR_N,err=0,i;
    int failed=1;
    err=lis_solver_create(&solver);
    if(err) goto done;
    if(solver->mfamg_deep_sibling_merge!=0) goto done;
    err=lis_solver_set_option("-mfamg_deep_sibling_merge 1",solver);
    if(err || solver->mfamg_deep_sibling_merge!=1) goto done;
    if(lis_solver_set_option("-mfamg_deep_sibling_merge 2",solver)==LIS_SUCCESS ||
       solver->mfamg_deep_sibling_merge!=1) goto done;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,STAR_N);
    if(err) goto done;
    for(i=0;i<STAR_N;i++) {
        id[i]=i; sizes[i]=0;
        err=lis_matrix_set_value(
            LIS_INS_VALUE,i,i,(i==0 ? 121.0 : 2.0),A);
        if(err) goto done;
        if(i>0) {
            err=lis_matrix_set_value(LIS_INS_VALUE,0,i,-1.0,A);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,i,0,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_precon_mfamg_merge_siblings_for_test(A,96,id,&nc);
    if(err) goto done;
    if(nc!=STAR_N-95) {
        fprintf(stderr,"MFAMG DEEP MERGE rooted star expected %d groups, got %d\n",
                STAR_N-95,(int)nc);
        goto done;
    }
    /* Rooted aggregate contains hub and 95 of the leaves.
     * No group composed solely of disconnected leaves is permitted. */
    if(id[0]<0 || id[0]>=nc) goto done;

    for(i=0;i<STAR_N;i++) {
        if(id[i]<0 || id[i]>=nc) goto done;
        if(++sizes[id[i]]>96) goto done;
    }
    if(sizes[id[0]]!=96) goto done;
    for(i=0;i<nc;i++)
        if(i!=id[0] && sizes[i]!=1) {
            fprintf(stderr,"MFAMG DEEP MERGE disconnected residual group size=%d\n",
                    (int)sizes[i]);
            goto done;
        }
    printf("MFAMG DEEP MERGE structural rooted star=%d groups=%d PASS\n",
           STAR_N,(int)nc);
    failed=0;
done:
    if(err) fprintf(stderr,"MFAMG DEEP MERGE structural rc=%d\n",(int)err);
    if(A) lis_matrix_destroy(A);
    if(solver) lis_solver_destroy(solver);
    return failed;
}

/* A regular 1-D chain has at most two incoming anchor choices.
 * fan_in>=3 must leave its aggregation unchanged. */
static int structural_line_gate_test(void)
{
    LIS_MATRIX A=NULL;
    LIS_INT ids[9],nc=9,i,err=0;
    int failed=1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,9);
    if(err) goto done;
    for(i=0;i<9;i++) {
        ids[i]=i;
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,2.5,A);
        if(err) goto done;
        if(i>0) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
            if(err) goto done;
        }
        if(i<8) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_precon_mfamg_merge_siblings_for_test(A,8,ids,&nc);
    if(err) goto done;
    if(nc!=9) {
        fprintf(stderr,"MFAMG DEEP MERGE line fan-in gate merged %d groups\n",
                (int)(9-nc));
        goto done;
    }
    failed=0;
done:
    if(A) lis_matrix_destroy(A);
    return failed;
}

/* Leaf edges with larger raw magnitude must not beat a stronger
 * diagonal-normalized anchor. SPD: diagonally dominant 6x6 matrix.
 * Four leaves each connect to hub 0 with magnitude 2 and to hub 1
 * with magnitude 1. The smaller-diagonal hub 1 must win.
 */
static int structural_normalized_anchor_test(void)
{
    LIS_MATRIX A=NULL;
    LIS_INT ids[6]={0,1,2,3,4,5},nc=6,i,err=0;
    int failed=1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,6);
    if(err) goto done;
    err=lis_matrix_set_value(LIS_INS_VALUE,0,0,100.0,A);
    if(err) goto done;
    err=lis_matrix_set_value(LIS_INS_VALUE,1,1,5.0,A);
    if(err) goto done;
    for(i=2;i<6;i++) {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,i,0,-2.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,0,i,-2.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,i,1,-1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,1,i,-1.0,A);
        if(err) goto done;
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_precon_mfamg_merge_siblings_for_test(A,5,ids,&nc);
    if(err) goto done;
    if(nc!=2 || ids[0]==ids[1]) goto done;
    for(i=2;i<6;i++) if(ids[i]!=ids[1]) goto done;
    printf("MFAMG DEEP MERGE normalized anchor selects hub 1 PASS\n");
    failed=0;
done:
    if(failed) fprintf(stderr,
        "MFAMG DEEP MERGE normalized anchor regression failed rc=%d nc=%d\n",
        (int)err,(int)nc);
    if(A) lis_matrix_destroy(A);
    return failed;
}

/* Assemble a diagonally shifted graph Laplacian.
 * Serial: 32 chains of length 32 meeting at node 0.
 * MPI: a distributed SPD line, avoiding a hub-dominated partition.
 */
static LIS_INT build_matrix(int distributed, LIS_MATRIX *A)
{
    LIS_INT g,n,degree,neighbors[3],count,k,err;
    n=distributed ? MPI_N : DEEP_N;
    err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return err;
    err=lis_matrix_set_size(*A,0,n);
    if(err) return err;
    for(g=(*A)->is;g<(*A)->ie;g++) {
        count=0;
        if(distributed) {
            if(g>0) neighbors[count++]=g-1;
            if(g+1<n) neighbors[count++]=g+1;
        } else if(g==0) {
            /* hub row has STAR_ARMS links: handled separately */
            for(k=0;k<STAR_ARMS;k++) {
                err=lis_matrix_set_value(
                    LIS_INS_VALUE,g,1+k*ARM_LENGTH,
                    (LIS_SCALAR)-1.0,*A);
                if(err) return err;
            }
            count=STAR_ARMS;
        } else {
            LIS_INT t=(g-1)%ARM_LENGTH;
            if(t==0) neighbors[count++]=0;
            else neighbors[count++]=g-1;
            if(t+1<ARM_LENGTH) neighbors[count++]=g+1;
        }
        degree=count;
        if(g!=0 || distributed) {
            for(k=0;k<count;k++) {
                err=lis_matrix_set_value(LIS_INS_VALUE,g,neighbors[k],
                                         (LIS_SCALAR)-1.0,*A);
                if(err) return err;
            }
        }
        err=lis_matrix_set_value(LIS_INS_VALUE,g,g,
                                 (LIS_SCALAR)degree+(LIS_SCALAR)1.0e-3,*A);
        if(err) return err;
    }
    return lis_matrix_assemble(*A);
}

static int integration_test(int distributed)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR modes[2]={NULL,NULL},exact=NULL,b=NULL,x=NULL,ax=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_INT levels[MAX_LEVELS],owns[MAX_LEVELS],level_baseline[MAX_LEVELS];
    LIS_INT nlevels=0,coarse_n=0,base_nlevels=0,base_coarse_n=0;
    LIS_INT attempts=0,events=0,reduced=0,max_depth=0;
    LIS_INT i,g,iter=0,baseline_iter=0,err=0,changed=0,phase;
    LIS_REAL norm_b=0.0,norm_r=0.0,physical=0.0,baseline=0.0;
    const int flags[3]={0,1,0};
    char option[200];
    int failed=1;

    err=build_matrix(distributed,&A);
    if(err) goto done;
    for(i=0;i<2;i++) {
        err=lis_vector_duplicate(A,&modes[i]);
        if(err) goto done;
    }
    err=lis_vector_duplicate(A,&exact);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&ax);
    if(err) goto done;
    err=lis_vector_duplicate(A,&r);
    if(err) goto done;

    for(i=0;i<A->n;i++) {
        g=A->is+i;
        modes[0]->value[i]=(LIS_SCALAR)1.0;
        /* Constant along each arm, smooth over the arm index.
         * Distributed line uses blocks of eight unknowns. */
        if(distributed)
            modes[1]->value[i]=(LIS_SCALAR)((g/8)+1)/32.0;
        else if(g==0)
            modes[1]->value[i]=(LIS_SCALAR)0.5;
        else
            modes[1]->value[i]=(LIS_SCALAR)((g-1)/ARM_LENGTH+1)/STAR_ARMS;
        exact->value[i]=(LIS_SCALAR)(1.0+0.02*(g%7));
    }
    err=lis_matvec(A,exact,b);
    if(err) goto done;
    err=lis_vector_nrm2(b,&norm_b);
    if(err || !isfinite((double)norm_b) || norm_b<=0.0) goto done;

    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_near_nullspace(solver,2,modes);
    if(err) goto done;
    for(phase=0;phase<3;phase++) {
        snprintf(option,sizeof(option),
          "-i gmres -p mfamg -scale none -print none "
          "-tol 1.0e-9 -maxiter 600 -restart 60 "
          "-mfamg_deep_sibling_merge %d",flags[phase]);
        err=lis_solver_set_option_replace_if_changed(option,solver,&changed);
        if(err || solver->mfamg_deep_sibling_merge!=flags[phase])
            goto done;
        if(phase>0 && changed!=LIS_TRUE) {
            fprintf(stderr,"MFAMG DEEP MERGE option transition not detected\n");
            goto done;
        }
        err=lis_solver_set_matrix(A,solver);
        if(err) goto done;
        /*
         * Construct the exact same MFAMG preconditioner path as lis_solve,
         * and inspect its internal, per-instance merge counters.
         */
        err=lis_precon_create(solver,&precon);
        if(err) goto done;
        err=lis_precon_mfamg_get_hierarchy_info(
            precon,MAX_LEVELS,&nlevels,&coarse_n,levels,owns);
        if(err) goto done;
        err=lis_precon_mfamg_get_deep_merge_info(
            precon,&attempts,&events,&reduced,&max_depth);
        if(err) goto done;
        if(distributed) {
            /* A shallow MPI hierarchy cannot prove that the serial guard
             * actually suppresses the deep-level merge call. */
            if(nlevels<6) {
                fprintf(stderr,"MFAMG DEEP MERGE MPI insufficient depth=%d (need level 4 active)\n",
                        (int)nlevels);
                goto done;
            }
            if(phase==0) {
                base_nlevels=nlevels;
                base_coarse_n=coarse_n;
                for(i=0;i<nlevels;i++) level_baseline[i]=levels[i];
            } else {
                if(nlevels!=base_nlevels || coarse_n!=base_coarse_n) {
                    fprintf(stderr,"MFAMG DEEP MERGE MPI option altered hierarchy\n");
                    goto done;
                }
                for(i=0;i<nlevels;i++)
                    if(levels[i]!=level_baseline[i]) goto done;
            }
            if(attempts!=0 || events!=0 || reduced!=0 || max_depth!=0) {
                fprintf(stderr,"MFAMG DEEP MERGE MPI must fall back (attempts=%d)\n",
                        (int)attempts);
                goto done;
            }
        } else {
            if(phase==0) {
                if(nlevels<6 || attempts!=0 || events!=0 || reduced!=0 ||
                   max_depth!=0) {
                    fprintf(stderr,"MFAMG DEEP MERGE baseline depth/counters invalid %d\n",
                            (int)nlevels);
                    goto done;
                }
                base_nlevels=nlevels; base_coarse_n=coarse_n;
                for(i=0;i<nlevels;i++) level_baseline[i]=levels[i];
            } else if(phase==1) {
                if(attempts<1 || events<1 || reduced<1 ||
                   max_depth<4 || nlevels<6) {
                    fprintf(stderr,
                       "MFAMG DEEP MERGE not activated: levels=%d attempts=%d events=%d reduced=%d depth=%d\n",
                       (int)nlevels,(int)attempts,(int)events,
                       (int)reduced,(int)max_depth);
                    goto done;
                }
                /* The next level after the deepest merge must shrink.
                 * This is less fragile than hardcoding level 5. */
                if(max_depth+1>=nlevels ||
                   max_depth+1>=base_nlevels ||
                   levels[max_depth+1]>=level_baseline[max_depth+1]) {
                    fprintf(stderr,
                       "MFAMG DEEP MERGE no reduction after depth=%d\n",
                       (int)max_depth);
                    goto done;
                }
            } else if(nlevels!=base_nlevels || coarse_n!=base_coarse_n ||
                      attempts!=0 || events!=0 ||
                      reduced!=0 || max_depth!=0) {
                fprintf(stderr,"MFAMG DEEP MERGE 1->0 did not restore hierarchy\n");
                goto done;
            } else {
                for(i=0;i<nlevels;i++) if(levels[i]!=level_baseline[i])
                    goto done;
            }
        }
        printf("MFAMG DEEP MERGE profile=%d distributed=%d levels=%d level5=%d attempts=%d merged=%d depth=%d\n",
               flags[phase],distributed,(int)nlevels,
               (int)(nlevels>5?levels[5]:-1),
               (int)attempts,(int)reduced,(int)max_depth);
        lis_precon_destroy(precon);
        precon=NULL;

        err=lis_vector_set_all((LIS_SCALAR)0.0,x);
        if(err) goto done;
        err=lis_solve(A,b,x,solver);
        if(err!=LIS_SUCCESS) {
            fprintf(stderr,"MFAMG DEEP MERGE solve phase=%d rc=%d\n",
                    phase,(int)err);
            goto done;
        }
        err=lis_solver_get_iter(solver,&iter);
        if(err || iter<=0 || iter>600) goto done;
        err=lis_matvec(A,x,ax);
        if(err) goto done;
        err=lis_vector_copy(b,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,ax,r);
        if(err) goto done;
        err=lis_vector_nrm2(r,&norm_r);
        if(err) goto done;
        physical=norm_r/norm_b;
        if(!isfinite((double)physical) || physical>=1.0e-8) {
            fprintf(stderr,"MFAMG DEEP MERGE bad physical residual phase=%d rel=%.9e\n",
                    phase,(double)physical);
            goto done;
        }
        if(phase==0) {
            baseline=physical;
            baseline_iter=iter;
        } else if((distributed || phase==2) &&
                  (iter!=baseline_iter || physical!=baseline)) {
            fprintf(stderr,
                "MFAMG DEEP MERGE reproducibility failure distributed=%d phase=%d baseline_iter=%d current_iter=%d baseline_rel=%.17e current_rel=%.17e\n",
                distributed,phase,(int)baseline_iter,(int)iter,
                (double)baseline,(double)physical);
            goto done;
        }
        printf("MFAMG DEEP MERGE SOLVE phase=%d mpi=%d iter=%d rel=%.9e PASS\n",
               phase,distributed,(int)iter,(double)physical);
    }
    failed=0;
done:
    if(err && failed) fprintf(stderr,"MFAMG DEEP MERGE integration rc=%d\n",(int)err);
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    for(i=0;i<2;i++) if(modes[i]) lis_vector_destroy(modes[i]);
    if(A) lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
    int failed=1,nprocs=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#ifdef USE_MPI
    MPI_Comm_size(LIS_COMM_WORLD,&nprocs);
#endif
    if(nprocs==1) {
        if(structural_star_test()!=0 ||
           structural_line_gate_test()!=0 ||
           structural_normalized_anchor_test()!=0) goto done;
    }
    failed=integration_test(nprocs>1);
done:
    lis_finalize();
    return failed;
}
