/* Real native-LIS integration regression: static condensation -> LIS FGMRES
 * -> full-space reconstruction and true residual.
 *
 * Default: native FGMRES / no preconditioner, tested by `make check`.
 * Optional: ./testcondense_lis mfamg | ca-mfamg | ilutp |
 *           mfamg_ilutp | ca-mfamg_ilutp
 *
 * All matrices are generated internally; no CalculiX geometry or external
 * numerical packages are required. The internal blocks are of size 2 and 4,
 * with DOFs interleaved with the retained physical DOFs.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lis.h"
#include "lis_condense_lis.h"
#include "lis_matrix_materialize.h"

#define NP 240
#define NQ 6
#define N (NP+NQ)
#define NB 2
#define NMODES 3
#define NEXT_MODES 4

/* An optional native MFAMG+ILUTP composite is provided by upstream develop.
 * CA mode is opt-in, with partition metadata propagated to the MFAMG child. */
static int uses_mfamg(const char *profile){
    return !strcmp(profile,"mfamg") || !strcmp(profile,"ca-mfamg") ||
           !strcmp(profile,"mfamg_ilutp") || !strcmp(profile,"ca-mfamg_ilutp");
}
static int uses_contact(const char *profile){
    return !strcmp(profile,"ca-mfamg") || !strcmp(profile,"ca-mfamg_ilutp");
}

static LIS_INT qd[NQ]={4,41,82,131,187,224};
static LIS_INT bp[NB+1]={0,2,6};
static LIS_INT phys[NP], full2phys[N], internal[N];

static void make_map(void){
    LIS_INT i,k=0;
    for(i=0;i<N;i++) {internal[i]=0;full2phys[i]=-1;}
    for(i=0;i<NQ;i++)internal[qd[i]]=1;
    for(i=0;i<N;i++)if(!internal[i]){
        phys[k]=i;full2phys[i]=k++;}
}

static LIS_INT put(LIS_MATRIX a,LIS_INT i,LIS_INT j,double v){
    return lis_matrix_set_value(LIS_INS_VALUE,i,j,(LIS_SCALAR)v,a);
}

static LIS_INT build_matrix(LIS_MATRIX *a,double diagonal_shift){
    LIS_INT i,j,k,err;
    err=lis_matrix_create(LIS_COMM_WORLD,a);if(err)return err;
    err=lis_matrix_set_size(*a,0,N);if(err)return err;
    for(i=0;i<N;i++){
        err=put(*a,i,i,5.2+diagonal_shift);if(err)return err;
    }
    /* 1-D algebraic elasticity surrogate with three contact partitions. */
    for(i=0;i<NP-1;i++){
        double w=(i==79 || i==159)?-0.32:-0.75;
        err=put(*a,phys[i],phys[i+1],w);if(err)return err;
        err=put(*a,phys[i+1],phys[i],w);if(err)return err;
    }
    /* Two independent local blocks with different sizes. */
    for(k=0;k<NB;k++)for(i=bp[k];i<bp[k+1];i++)
        for(j=bp[k];j<bp[k+1];j++)if(i!=j){
            err=put(*a,qd[i],qd[j],-0.11);if(err)return err;
        }
    for(i=0;i<NQ;i++){
        /* Coupling is symmetric, but block metadata makes no such assumption. */
        LIS_INT p=(13*i+21)%NP;
        err=put(*a,qd[i],phys[p],-0.13);if(err)return err;
        err=put(*a,phys[p],qd[i],-0.13);if(err)return err;
        err=put(*a,qd[i],phys[(p+3)%NP],-0.07);if(err)return err;
        err=put(*a,phys[(p+3)%NP],qd[i],-0.07);if(err)return err;
    }
    return lis_matrix_assemble(*a);
}

/* A second nonlinear-like coefficient update, with the same CSR sparsity
 * and the same internal-DOF metadata. Crucially this function builds a NEW
 * Krylov solver and new vectors after the session's Schur was replaced;
 * no AMG or ILUTP setup is assumed reusable until LIS exposes that safely. */
static int run_update_step(lis_sc_lis_session *sc, LIS_VECTOR exact,
                            const LIS_INT labels[NP],
                            const char *profile, char *options){
    LIS_MATRIX next_a=NULL, s=NULL, old=NULL;
    LIS_VECTOR next_b=NULL,next_g=NULL,next_u=NULL,next_x=NULL;
    LIS_VECTOR modes[NEXT_MODES]={NULL,NULL,NULL,NULL};
    LIS_VECTOR full_modes[NEXT_MODES]={NULL,NULL,NULL,NULL};
    LIS_SOLVER next_solver=NULL;
    LIS_INT err=0,it=0,k;
    double rr=1;
    int failed=1;
    err=build_matrix(&next_a,0.125);if(err)goto done;
    err=lis_vector_duplicate(next_a,&next_b);if(err)goto done;
    err=lis_vector_duplicate(next_a,&next_x);if(err)goto done;
    err=lis_matvec(next_a,exact,next_b);if(err)goto done;
    old=lis_sc_lis_matrix(sc);
    if(lis_sc_lis_update(sc,next_a,NB,bp,qd,1e-13,LIS_SC_NUMERIC)!=LIS_SC_OK)goto done;
    s=lis_sc_lis_matrix(sc);
    if(!s||s==old||s->n!=NP)goto done;
    err=lis_vector_duplicate(s,&next_g);if(err)goto done;
    err=lis_vector_duplicate(s,&next_u);if(err)goto done;
    if(lis_sc_lis_rhs(sc,next_b,next_g)!=LIS_SC_OK)goto done;
    /* The partition count and near-nullspace dimension both change between
     * steps: 3 -> 4. These are application-supplied algebraic modes, never
     * inferred from DOF number or element type by the preconditioner. */
    for(k=0;k<NEXT_MODES;k++){
        err=lis_vector_duplicate(next_a,&full_modes[k]);if(err)goto done;
        err=lis_vector_duplicate(s,&modes[k]);if(err)goto done;
        for(LIS_INT i=0;i<NP;i++){
            err=lis_vector_set_value(LIS_INS_VALUE,phys[i],
                 (LIS_SCALAR)(labels[i]==k?1.0:0.0),full_modes[k]);
            if(err)goto done;
        }
    }
    if(lis_sc_lis_modes(sc,NEXT_MODES,full_modes,modes)!=LIS_SC_OK)goto done;
    err=lis_solver_create(&next_solver);if(err)goto done;
    err=lis_solver_set_option(options,next_solver);if(err)goto done;
    if(uses_mfamg(profile)){
        err=lis_solver_set_near_nullspace(next_solver,NEXT_MODES,modes);if(err)goto done;
    }
    if(uses_contact(profile)){
        err=lis_solver_set_mfamg_contact_partition(next_solver,NP,labels);if(err)goto done;
    }
    err=lis_solve(s,next_g,next_u,next_solver);if(err)goto done;
    if(lis_sc_lis_reconstruct(sc,next_u,next_b,next_x)!=LIS_SC_OK)goto done;
    if(lis_sc_lis_residual(sc,next_b,next_x,&rr)!=LIS_SC_OK||!isfinite(rr)||rr>=1e-6)goto done;
    (void)lis_solver_get_iter(next_solver,&it);
    printf("LIS updated Schur: %s PASS; n=%d, nphys=%d, true residual=%.3e, iters=%d\n",
           profile,N,NP,rr,(int)it);
    failed=0;
done:
    if(failed)fprintf(stderr,"LIS updated Schur: %s FAIL (err=%d, residual=%.3e)\n",
                      profile,(int)err,rr);
    if(next_solver)lis_solver_destroy(next_solver);
    for(k=0;k<NEXT_MODES;k++){
        if(modes[k])lis_vector_destroy(modes[k]);
        if(full_modes[k])lis_vector_destroy(full_modes[k]);
    }
    if(next_g)lis_vector_destroy(next_g);
    if(next_u)lis_vector_destroy(next_u);
    if(next_x)lis_vector_destroy(next_x);
    if(next_b)lis_vector_destroy(next_b);
    if(next_a)lis_matrix_destroy(next_a);
    return failed;
}

/* Negative contract: CA composite MUST reject absent partition metadata.
 * If the MFAMG child forgets contact-aware flags, this silently succeeds
 * and the regression fails. Valid labels are tested in the normal path. */
static int check_missing_ca_labels(LIS_MATRIX s, LIS_VECTOR g,
                                   LIS_VECTOR const modes[NMODES]){
    LIS_SOLVER probe=NULL;
    LIS_VECTOR x=NULL;
    LIS_INT err;
    int failed=1;
    char options[]= "-i fgmres -p mfamg_ilutp -mfamg_contact_aware 1 "
        "-scale none -print none -maxiter 40 -restart 20 -tol 1e-8 "
        "-ilutp_threshold 0.5 -ilut_drop 0.005 -ilut_rate 8";
    err=lis_vector_duplicate(s,&x);if(err)goto done;
    err=lis_solver_create(&probe);if(err)goto done;
    err=lis_solver_set_option(options,probe);if(err)goto done;
    err=lis_solver_set_near_nullspace(probe,NMODES,(LIS_VECTOR*)modes);
    if(err)goto done;
    /* Deliberately do NOT call lis_solver_set_mfamg_contact_partition(). */
    err=lis_solve(s,g,x,probe);
    if(!err){
        fprintf(stderr,"CA composite ERROR: accepted missing partition labels\n");
        goto done;
    }
    puts("CA composite rejects missing partition labels: PASS");
    failed=0;
done:
    if(probe)lis_solver_destroy(probe);
    if(x)lis_vector_destroy(x);
    return failed;
}

/* Native reference solve on the ORIGINAL matrix, not a second Schur solve.
 * This verifies physical residual independently of the condensation path.
 * Baseline timings in this 246x246 manufactured test are not benchmarks. */
static int run_full_reference(LIS_MATRIX A,LIS_VECTOR b,
                              lis_sc_lis_session *sc){
    LIS_SOLVER solver=NULL;
    LIS_VECTOR fullx=NULL;
    LIS_INT err,it=0;
    double relative=1.0;
    char options[]="-i fgmres -p none -scale none -print none "
                   "-maxiter 500 -restart 120 -tol 1e-8";
    int failed=1;
    err=lis_vector_duplicate(A,&fullx);if(err)goto done;
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(options,solver);if(err)goto done;
    err=lis_solve(A,b,fullx,solver);if(err)goto done;
    if(lis_sc_lis_residual(sc,b,fullx,&relative)!=LIS_SC_OK ||
       !isfinite(relative)||relative>=1e-6)goto done;
    (void)lis_solver_get_iter(solver,&it);
    printf("LIS full reference: no preconditioner PASS; true residual=%.3e, iters=%d\n",
           relative,(int)it);
    failed=0;
done:
    if(failed)fprintf(stderr,"LIS full reference FAIL, error=%d, res=%.3e\n",
                      (int)err,relative);
    if(solver)lis_solver_destroy(solver);
    if(fullx)lis_vector_destroy(fullx);
    return failed;
}

static int check_managed_session(void){
    LIS_MATRIX a=NULL,next=NULL,old_s=NULL;
    LIS_VECTOR exact=NULL,b=NULL,x=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    lis_sc_solver_session *session=NULL;
    lis_sc_solver_report report;
    LIS_REAL nr,nb;
    int failed=1;

    if(build_matrix(&a,0.0) || lis_vector_duplicate(a,&exact) ||
       lis_vector_duplicate(a,&b) || lis_vector_duplicate(a,&x) ||
       lis_vector_duplicate(a,&r) || lis_vector_set_all((LIS_SCALAR)1,exact) ||
       lis_matvec(a,exact,b))goto done;
    if(lis_sc_solver_session_create(a,NB,bp,qd,1e-13,&session)!=LIS_SC_OK)
        goto done;
    solver=lis_sc_solver_session_get_solver(session);
    if(!solver || lis_solver_set_option(
        "-i fgmres -p none -scale none -print none -maxiter 500 -restart 120 -tol 1e-9",solver))
        goto done;
    old_s=lis_sc_lis_matrix(lis_sc_solver_session_get_condensation(session));
    if(lis_sc_solver_session_solve(session,b,x,1e-7)!=LIS_SC_OK ||
       lis_sc_solver_session_get_report(session,&report)!=LIS_SC_OK ||
       !report.solver_called || !report.solution_committed ||
       !report.full_residual_valid || report.full_relative_residual>1e-7)
        goto done;
    if(lis_sc_solver_session_update(session,a,NB,bp,qd,1e-13,LIS_SC_REUSE_ALL)!=LIS_SC_OK ||
       lis_sc_lis_matrix(lis_sc_solver_session_get_condensation(session))!=old_s ||
       lis_sc_solver_session_solve(session,b,x,1e-7)!=LIS_SC_OK)
        goto done;
    if(build_matrix(&next,0.125) || lis_matvec(next,exact,b))goto done;
    if(lis_sc_solver_session_update(session,next,NB,bp,qd,1e-13,LIS_SC_REUSE_ALL)==LIS_SC_OK ||
       lis_sc_lis_matrix(lis_sc_solver_session_get_condensation(session))!=old_s)
        goto done;
    if(lis_sc_solver_session_update(session,next,NB,bp,qd,1e-13,LIS_SC_NUMERIC)!=LIS_SC_OK ||
       lis_sc_lis_matrix(lis_sc_solver_session_get_condensation(session))==old_s ||
       lis_sc_solver_session_get_report(session,&report)==LIS_SC_OK)
        goto done;
    if(lis_vector_set_all((LIS_SCALAR)7,x) ||
       lis_sc_solver_session_solve(session,b,x,-1)!=LIS_SC_BAD_ARGUMENT ||
       lis_sc_solver_session_get_report(session,&report)!=LIS_SC_OK ||
       report.solver_called || report.solution_committed ||
       x->value[0]!=(LIS_SCALAR)7)
        goto done;
    if(lis_sc_solver_session_solve(session,b,x,1e-7)!=LIS_SC_OK ||
       lis_matvec(next,x,r) || lis_vector_axpy((LIS_SCALAR)-1,b,r) ||
       lis_vector_nrm2(r,&nr) || lis_vector_nrm2(b,&nb) || !(nb>0) ||
       !(nr/nb<1e-7))
        goto done;
    failed=0;
done:
    lis_sc_solver_session_destroy(session);
    if(r)lis_vector_destroy(r);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(exact)lis_vector_destroy(exact);
    if(next)lis_matrix_destroy(next);
    if(a)lis_matrix_destroy(a);
    printf("MANAGED_SCHUR_SESSION %s\n",failed?"FAILED":"PASSED");
    return failed;
}

/* Use a genuinely different Schur P while certifying against physical S.
 * Exercise native MFAMG/CA-MFAMG and their ILUTP composite with the
 * explicit interpolation basis, never with mathematical near-nullspace. */
static int check_paired_schur_profile(LIS_MATRIX S,LIS_VECTOR g,
                  LIS_VECTOR modes[],const LIS_INT labels[],const char *profile,
                  char *options)
{
    LIS_MATRIX shifted=NULL,P=NULL;
    lis_sc_lis_session *pc_sc=NULL;
    LIS_MATRIX_AUTO_SESSION session=NULL;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR x=NULL,r=NULL;
    LIS_PRECON kept=NULL;
    LIS_REAL nr,nb;
    LIS_INT status,pp,pa,pass;
    int failed=1;
    if(!uses_mfamg(profile) && strcmp(profile,"ilutp"))return 0;
    if(build_matrix(&shifted,0.125) ||
       lis_sc_lis_create(shifted,NB,bp,qd,1e-13,&pc_sc)!=LIS_SC_OK)
        goto done;
    P=lis_sc_lis_matrix(pc_sc);
    if(!P || P==S || lis_vector_duplicate(S,&x) ||
       lis_vector_duplicate(S,&r) || lis_matrix_auto_session_create(&session) ||
       lis_matrix_auto_session_get_solver(session,&solver) ||
       lis_solver_set_option(options,solver) ||
       lis_solver_set_num_threads(solver,1))goto done;
    if(uses_mfamg(profile) && lis_solver_set_amg_coarse_basis(solver,NMODES,modes))
        goto done;
    if(uses_contact(profile) &&
       lis_solver_set_mfamg_contact_partition(solver,NP,labels))goto done;
    for(pass=0;pass<6;pass++)
    {
        /* Second pair of solves tests independent P with local RAS
         * factors. Standalone MFAMG has no ILUT child. */
        if(pass==2)
        {
            if(strcmp(profile,"ilutp") &&
               strcmp(profile,"mfamg_ilutp") &&
               strcmp(profile,"ca-mfamg_ilutp"))break;
            if(lis_solver_set_option(
                "-ilut_ras_overlap 1 -ilut_ras_domains 2",solver))
                goto done;
            kept=NULL;
        }
        if(pass==4)
        {
            if(lis_solver_set_option("-ilut_schur 1",solver))goto done;
            kept=NULL;
        }
        if(lis_vector_set_all((LIS_SCALAR)0,x) ||
           lis_matrix_auto_session_solve_with_precon_matrix(
                session,S,P,g,x,NULL,NULL,&pp,&pa) ||
           pp!=LIS_MATRIX_POLICY_DIRECT || pa!=LIS_MATRIX_POLICY_ACTION_DIRECT ||
           lis_solver_get_status(solver,&status) || status!=LIS_SUCCESS ||
           !solver->sequence_precon || solver->sequence_precon_matrix!=P ||
           (kept && solver->sequence_precon!=kept) ||
           lis_matvec(S,x,r) || lis_vector_axpy((LIS_SCALAR)-1,g,r) ||
           lis_vector_nrm2(r,&nr) || lis_vector_nrm2(g,&nb) ||
           !(nb>0) || !(nr/nb<1e-7))goto done;
        kept=solver->sequence_precon;
    }
    failed=0;
done:
    lis_matrix_auto_session_destroy(session);
    if(r)lis_vector_destroy(r);
    if(x)lis_vector_destroy(x);
    lis_sc_lis_destroy(pc_sc);
    if(shifted)lis_matrix_destroy(shifted);
    printf("PAIRED_SCHUR_%s %s\n",profile,failed?"FAILED":"PASSED");
    return failed;
}

static int run(int argc,char **argv){
    LIS_MATRIX a=NULL,s=NULL;
    LIS_VECTOR xtrue=NULL,b=NULL,g=NULL,u=NULL,xfull=NULL;
    LIS_VECTOR allmode[NMODES]={NULL,NULL,NULL};
    LIS_VECTOR redmode[NMODES]={NULL,NULL,NULL};
    LIS_SOLVER solver=NULL;
    lis_sc_lis_session *sc=NULL;
    LIS_INT labels[NP],changed_labels[NP],full_labels[N],full_changed[N],i,k,err,iter=0;
    double rr=1.0,rr_exact=1.0;
    int failed=1;
    const char *profile=argc>1?argv[1]:"none";
    char options[220];

    make_map();
    if(!strcmp(profile,"none") && check_managed_session())return 1;
    if(strcmp(profile,"none") && strcmp(profile,"mfamg") &&
       strcmp(profile,"ca-mfamg") && strcmp(profile,"ilutp") &&
       strcmp(profile,"mfamg_ilutp") && strcmp(profile,"ca-mfamg_ilutp")){
        fprintf(stderr,"Usage: testcondense_lis [none|mfamg|ca-mfamg|ilutp|mfamg_ilutp|ca-mfamg_ilutp]\n");
        return 2;
    }
    err=build_matrix(&a,0.0);if(err)goto done;
    err=lis_vector_duplicate(a,&xtrue);if(err)goto done;
    err=lis_vector_duplicate(a,&b);if(err)goto done;
    err=lis_vector_duplicate(a,&xfull);if(err)goto done;
    for(i=0;i<N;i++){
        err=lis_vector_set_value(LIS_INS_VALUE,i,
               (LIS_SCALAR)(1.2+0.3*sin((double)(i+1)*0.023)),xtrue);
        if(err)goto done;
    }
    err=lis_matvec(a,xtrue,b);if(err)goto done;
    if(lis_sc_lis_create(a,NB,bp,qd,1e-13,&sc)!=LIS_SC_OK)goto done;
    s=lis_sc_lis_matrix(sc);
    if(!s || s->n!=NP)goto done;
    if(run_full_reference(a,b,sc))goto done;
    err=lis_vector_duplicate(s,&g);if(err)goto done;
    err=lis_vector_duplicate(s,&u);if(err)goto done;
    if(lis_sc_lis_rhs(sc,b,g)!=LIS_SC_OK)goto done;

    /* Independent exact solution reconstruction checks the operator before
     * evaluating a Krylov convergence claim. */
    for(i=0;i<NP;i++)u->value[i]=xtrue->value[phys[i]];
    if(lis_sc_lis_reconstruct(sc,u,b,xfull)!=LIS_SC_OK)goto done;
    if(lis_sc_lis_residual(sc,b,xfull,&rr_exact)!=LIS_SC_OK || rr_exact>2e-12)goto done;

    for(i=0;i<N;i++){
        full_labels[i]=-1;
        full_changed[i]=-1;
    }
    for(i=0;i<NP;i++){
        full_labels[phys[i]]=i/80; /* first step: three contact groups */
        full_changed[phys[i]]=i/60; /* second step: four groups */
    }
    /* Short labels buffers must fail without touching the destination. */
    labels[0]=-7777;
    if(lis_sc_lis_labels(sc,full_labels,N-1,labels,NP)!=LIS_SC_BAD_ARGUMENT ||
       labels[0]!=-7777 ||
       lis_sc_lis_labels(sc,full_labels,N,labels,NP-1)!=LIS_SC_BAD_ARGUMENT ||
       labels[0]!=-7777)goto done;
    /* Exercise the actual full->reduced partition projection API. */
    if(lis_sc_lis_labels(sc,full_labels,N,labels,NP)!=LIS_SC_OK ||
       lis_sc_lis_labels(sc,full_changed,N,changed_labels,NP)!=LIS_SC_OK)goto done;
    if(labels[0]!=0 || labels[NP-1]!=2 ||
       changed_labels[0]!=0 || changed_labels[NP-1]!=3)goto done;
    for(k=0;k<NMODES;k++){
        err=lis_vector_duplicate(a,&allmode[k]);if(err)goto done;
        err=lis_vector_duplicate(s,&redmode[k]);if(err)goto done;
        for(i=0;i<NP;i++){
            err=lis_vector_set_value(LIS_INS_VALUE,phys[i],
                (LIS_SCALAR)((i/80==k)?1.0:0.0),allmode[k]);
            if(err)goto done;
        }
    }
    if(lis_sc_lis_modes(sc,NMODES,allmode,redmode)!=LIS_SC_OK)goto done;
    for(i=0;i<NP;i++){
        for(k=0;k<NMODES;k++){
            double want=(i/80==k)?1.0:0.0;
            if(fabs(redmode[k]->value[i]-want)>1e-15)goto done;
        }
    }

    if(!strcmp(profile,"ca-mfamg_ilutp") &&
       check_missing_ca_labels(s,g,redmode))goto done;
    err=lis_solver_create(&solver);if(err)goto done;
    if(!strcmp(profile,"ilutp")){
        strcpy(options,"-i fgmres -p ilut -ilutp 1 -scale none -print none -maxiter 500 -restart 120 -tol 1e-8");
    }else if(!strcmp(profile,"mfamg_ilutp") || !strcmp(profile,"ca-mfamg_ilutp")){
        const char *extra=uses_contact(profile)?" -mfamg_contact_aware 1":"";
        if(snprintf(options,sizeof(options),
            "-i fgmres -p mfamg_ilutp -scale none -print none -maxiter 500 "
            "-restart 120 -tol 1e-8 -ilutp_threshold 0.5 -ilut_drop 0.005 -ilut_rate 8%s",
            extra)>=(int)sizeof(options))goto done;
    }else{
        if(snprintf(options,sizeof(options),
            "-i fgmres -p %s -scale none -print none -maxiter 500 -restart 120 -tol 1e-8",
            profile)>=(int)sizeof(options))goto done;
    }
    err=lis_solver_set_option(options,solver);if(err)goto done;
    if(uses_mfamg(profile)){
        err=lis_solver_set_near_nullspace(solver,NMODES,redmode);if(err)goto done;
    }
    if(uses_contact(profile)){
        err=lis_solver_set_mfamg_contact_partition(solver,NP,labels);if(err)goto done;
    }
    err=lis_solve(s,g,u,solver);if(err){
        fprintf(stderr,"native LIS solve failed, status=%d, profile=%s\n",(int)err,profile);
        goto done;
    }
    if(lis_sc_lis_reconstruct(sc,u,b,xfull)!=LIS_SC_OK)goto done;
    if(lis_sc_lis_residual(sc,b,xfull,&rr)!=LIS_SC_OK || !isfinite(rr) || rr>=1e-6)goto done;
    (void)lis_solver_get_iter(solver,&iter);
    printf("LIS condensation integration: %s PASS; n=%d, nphys=%d, ",profile,N,NP);
    printf("exact reconstruction=%.3e, true residual=%.3e, iters=%d\n",rr_exact,rr,(int)iter);
    /* All objects referring to the previous Schur are destroyed before update.
     * The original full matrix, b and exact solution remain caller-owned. */
    if(lis_sc_lis_update(sc,a,NB,bp,qd,1e-13,LIS_SC_REUSE_ALL)!=LIS_SC_OK ||
       lis_sc_lis_matrix(sc)!=s)goto done;
    lis_solver_destroy(solver);solver=NULL;
    if(check_paired_schur_profile(s,g,redmode,labels,profile,options))goto done;
    lis_vector_destroy(g);g=NULL;
    lis_vector_destroy(u);u=NULL;
    for(k=0;k<NMODES;k++) {lis_vector_destroy(redmode[k]);redmode[k]=NULL;}
    if(run_update_step(sc,xtrue,changed_labels,profile,options))goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"LIS condensation integration: %s FAIL (last error=%d, true_res=%.3e)\n",profile,(int)err,rr);
    if(solver)lis_solver_destroy(solver);
    for(k=0;k<NMODES;k++){
        if(allmode[k])lis_vector_destroy(allmode[k]);
        if(redmode[k])lis_vector_destroy(redmode[k]);
    }
    if(g)lis_vector_destroy(g);
    if(u)lis_vector_destroy(u);
    if(xfull)lis_vector_destroy(xfull);
    if(b)lis_vector_destroy(b);
    if(xtrue)lis_vector_destroy(xtrue);
    lis_sc_lis_destroy(sc);
    if(a)lis_matrix_destroy(a);
    return failed;
}
int main(int argc,char **argv){
    int status;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 1;
#if defined(USE_MPI) || defined(USE_COMPLEX) || defined(_COMPLEX)
    (void)run;
    status=77; /* Automake treats exit code 77 as a genuine SKIP. */
    puts("LIS condensation integration: SKIP (MPI / complex)");
#else
    status=run(argc,argv);
#endif
    lis_finalize();
    return status;
}
