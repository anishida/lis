/* Generic pair-session regression: matrix lifetime, copied results,
 * value/structure update, callback failure, and repeat solve.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "lis.h"
#include "lis_matrix_materialize.h"

typedef struct {
    LIS_INT n,fail_row,dense,reads;
    LIS_REAL scale;
} ROW_CONTEXT;

static LIS_INT get_row(void *p,LIS_INT row,LIS_INT *nnz,
    LIS_INT *cols,LIS_SCALAR *vals)
{
    ROW_CONTEXT *c=(ROW_CONTEXT *)p;
    LIS_INT required,k=0,cap;
    if(!c || !nnz || row<0 || row>=c->n) return LIS_ERR_ILL_ARG;
    c->reads++;
    if(c->fail_row==row) return LIS_ERR_ILL_ARG;
    required=1+(c->dense && row>0 ? 1:0);
    cap=*nnz;
    if(!cols && !vals){*nnz=required;return LIS_SUCCESS;}
    if(!cols || !vals || cap<required)return LIS_ERR_ILL_ARG;
    if(c->dense && row>0){cols[k]=row-1;vals[k++]=(LIS_SCALAR)0.0;}
    cols[k]=row;
    vals[k++]=(LIS_SCALAR)(c->scale*(row+2));
    *nnz=k;
    return LIS_SUCCESS;
}

static LIS_INT matvec(void *p,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    ROW_CONTEXT *c=(ROW_CONTEXT *)p;
    LIS_INT i;
    if(!c || !x || !y)return LIS_ERR_ILL_ARG;
    for(i=0;i<c->n;i++)
        y[i]=(LIS_SCALAR)(c->scale*(i+2))*x[i];
    return LIS_SUCCESS;
}

static LIS_INT make_user(ROW_CONTEXT *c,LIS_MATRIX *A)
{
    LIS_INT err;
    *A=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,A);if(err)return err;
    err=lis_matrix_set_size(*A,0,c->n);if(err)return err;
    err=lis_matrix_set_user(*A,c,matvec,matvec);if(err)return err;
    return lis_matrix_set_user_row(*A,get_row);
}

static LIS_INT make_identity(LIS_INT n,LIS_MATRIX *B)
{
    LIS_INT i,err;
    *B=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,B);if(err)return err;
    err=lis_matrix_set_size(*B,0,n);if(err)return err;
    for(i=0;i<n;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,*B);
        if(err)return err;
    }
    err=lis_matrix_set_type(*B,LIS_MATRIX_CSR);if(err)return err;
    return lis_matrix_assemble(*B);
}

int main(int argc,char **argv)
{
    ROW_CONTEXT ctx={6,-1,0,0},bctx={6,-1,0,0};
    LIS_MATRIX A=NULL,B=NULL,C=NULL;
    LIS_VECTOR x=NULL,out=NULL,kept=NULL;
    LIS_GESOLVE_MATRIX_SESSION session=NULL;
    LIS_ESOLVER e=NULL;
    LIS_SCALAR lambda,oldvalue;
    LIS_INT err,pa=-123,pb=-456,nev=0;
    int failed=1;
    ctx.scale=(LIS_REAL)1;
    bctx.scale=(LIS_REAL)0.5;
    if(lis_initialize(&argc,&argv))return 1;
    err=make_user(&ctx,&A);if(err)goto done;
    err=make_identity(ctx.n,&B);if(err)goto done;
    err=lis_vector_duplicate(B,&x);if(err)goto done;
    err=lis_vector_duplicate(B,&out);if(err)goto done;
    err=lis_vector_duplicate(B,&kept);if(err)goto done;
    err=lis_gesolve_matrix_session_create(&session);if(err)goto done;
    err=lis_gesolve_matrix_session_get_esolver(session,&e);if(err)goto done;
    err=lis_esolver_set_option(
        "-e gii -emaxiter 300 -etol 1e-8 -eprint none",e);
    if(err)goto done;
    err=lis_vector_set_all((LIS_SCALAR)1,x);if(err)goto done;
    err=lis_gesolve_matrix_session_solve(session,A,B,x,&pa,&pb);
    if(err || pa!=LIS_MATRIX_POLICY_CACHED_CSR ||
       pb!=LIS_MATRIX_POLICY_DIRECT)goto done;
    err=lis_gesolve_matrix_session_get_nev(session,&nev);
    if(err || nev!=1)goto done;
    err=lis_gesolve_matrix_session_get_eigenpair(
        session,0,&lambda,kept);
    if(err)goto done;
    memcpy(&oldvalue,&lambda,sizeof(lambda));
    err=lis_gesolve_matrix_session_get_eigenpair(
        session,-1,&lambda,out);
    if(err==LIS_SUCCESS)goto done;

    /* Structural and value changes must not expose dangling CSR pointers. */
    ctx.scale=(LIS_REAL)2;
    ctx.dense=1;
    err=lis_matrix_mark_values_changed(A);if(err)goto done;
    err=lis_vector_set_all((LIS_SCALAR)1,x);if(err)goto done;
    err=lis_gesolve_matrix_session_solve(session,A,B,x,&pa,&pb);
    if(err || pa!=LIS_MATRIX_POLICY_CACHED_CSR)goto done;
    err=lis_gesolve_matrix_session_get_nev(session,&nev);
    if(err || nev!=1)goto done;

    /* B-side USER callback failure must invalidate BOTH caches/results
     * after successfully preparing the A-side cache. */
    err=make_user(&bctx,&C);if(err)goto done;
    err=lis_vector_set_all((LIS_SCALAR)1,x);if(err)goto done;
    err=lis_gesolve_matrix_session_solve(session,A,C,x,&pa,&pb);
    if(err || pa!=LIS_MATRIX_POLICY_CACHED_CSR ||
       pb!=LIS_MATRIX_POLICY_CACHED_CSR)goto done;
    bctx.fail_row=1;
    err=lis_matrix_mark_values_changed(C);if(err)goto done;
    pa=101;pb=202;
    err=lis_gesolve_matrix_session_solve(session,A,C,x,&pa,&pb);
    if(err==LIS_SUCCESS || pa!=101 || pb!=202)goto done;
    if(lis_gesolve_matrix_session_get_nev(session,&nev)==LIS_SUCCESS)
        goto done;
    bctx.fail_row=-1;
    err=lis_vector_set_all((LIS_SCALAR)1,x);if(err)goto done;
    err=lis_gesolve_matrix_session_solve(session,A,C,x,&pa,&pb);
    if(err)goto done;

    /* Subspace mode: publish a genuinely indexed collection rather than
     * treating every generalized eigensolver as single-pair only. */
    err=lis_gesolve_matrix_session_get_esolver(session,&e);
    if(err)goto done;
    err=lis_esolver_set_option(
        "-e gsi -ss 3 -nev 2 -emaxiter 400 -etol 1e-7 -eprint none",e);
    if(err)goto done;
    err=lis_vector_set_all((LIS_SCALAR)1,x);if(err)goto done;
    err=lis_gesolve_matrix_session_solve(session,A,B,x,&pa,&pb);
    if(err)goto done;
    err=lis_gesolve_matrix_session_get_nev(session,&nev);
    if(err || nev<2)goto done;
    err=lis_gesolve_matrix_session_get_eigenpair(
        session,1,&lambda,out);
    if(err)goto done;

    err=lis_gesolve_matrix_session_reset(session);if(err)goto done;
    if(lis_gesolve_matrix_session_get_nev(session,&nev)==LIS_SUCCESS)
        goto done;
    /* Caller-owned eigenvector remains valid after reset. */
    if(kept->value==NULL || kept->n!=ctx.n)goto done;
    failed=0;
done:
    printf("GENERIC_PAIR_SESSION=%s\n",failed?"FAIL":"PASS");
    if(session)lis_gesolve_matrix_session_destroy(session);
    if(kept)lis_vector_destroy(kept);
    if(out)lis_vector_destroy(out);
    if(x)lis_vector_destroy(x);
    if(C)lis_matrix_destroy(C);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
