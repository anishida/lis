/* Numerical and OpenMP regression for the legacy BSR kernels.
 * Uses stack-owned LIS descriptors to exercise the kernels in isolation.
 */
#include "lis.h"
#include "lis_matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern LIS_INT lis_matrix_normf_bsr(LIS_MATRIX A,LIS_SCALAR *norm);

static int close_to(LIS_SCALAR got,double expected)
{
    return fabs((double)got-expected) <= 1e-10*(1.0+fabs(expected));
}

static int check_block_scaling(int nblocks)
{
    struct LIS_MATRIX_STRUCT a;
    struct LIS_MATRIX_CORE_STRUCT lower,upper;
    struct LIS_MATRIX_DIAG_STRUCT diag,scale;
    LIS_INT *lp=NULL,*up=NULL,*li=NULL,*ui=NULL;
    LIS_SCALAR *lv=NULL,*uv=NULL,*dv=NULL,*s=NULL;
    int i,k,failed=1;
    memset(&a,0,sizeof(a));
    memset(&lower,0,sizeof(lower));
    memset(&upper,0,sizeof(upper));
    memset(&diag,0,sizeof(diag));
    memset(&scale,0,sizeof(scale));
    lp=(LIS_INT *)calloc((size_t)nblocks+1,sizeof(*lp));
    up=(LIS_INT *)calloc((size_t)nblocks+1,sizeof(*up));
    li=(LIS_INT *)calloc((size_t)nblocks,sizeof(*li));
    ui=(LIS_INT *)calloc((size_t)nblocks,sizeof(*ui));
    lv=(LIS_SCALAR *)calloc((size_t)nblocks*9,sizeof(*lv));
    uv=(LIS_SCALAR *)calloc((size_t)nblocks*9,sizeof(*uv));
    dv=(LIS_SCALAR *)calloc((size_t)nblocks*9,sizeof(*dv));
    s=(LIS_SCALAR *)calloc((size_t)nblocks*9,sizeof(*s));
    if(!lp || !up || !li || !ui || !lv || !uv || !dv || !s) goto done;
    for(i=0;i<nblocks;i++)
    {
        lp[i]=up[i]=i;
        li[i]=(i+1)%nblocks; ui[i]=(i+2)%nblocks;
        for(k=0;k<9;k++)
        {
            lv[9*i+k]=(LIS_SCALAR)((k+1)*(i%13+1));
            uv[9*i+k]=(LIS_SCALAR)((k+2)*(i%11+1));
        }
        s[9*i+0]=2; s[9*i+4]=3; s[9*i+8]=4;
    }
    lp[nblocks]=up[nblocks]=nblocks;
    lower.bptr=lp; lower.bindex=li; lower.value=lv;
    upper.bptr=up; upper.bindex=ui; upper.value=uv;
    diag.value=dv; scale.value=s;
    a.nr=nblocks; a.bnr=a.bnc=3; a.is_splited=1;
    a.L=&lower; a.U=&upper; a.D=&diag;
    if(lis_matrix_bscale_bsr(&a,&scale)!=LIS_SUCCESS) goto done;
    for(i=0;i<nblocks;i++)
    {
        for(k=0;k<9;k++)
        {
            int row=k%3;
            double fac=(double)(row+2);
            double expected_l=(double)((k+1)*(i%13+1))*fac;
            double expected_u=(double)((k+2)*(i%11+1))*fac;
            if(!close_to(lv[9*i+k],expected_l) ||
               !close_to(uv[9*i+k],expected_u))
            {
                fprintf(stderr,"BSR scale error block %d entry %d\n",i,k);
                goto done;
            }
        }
        for(k=0;k<9;k++)
        {
            if(!close_to(dv[9*i+k],(k==0||k==4||k==8)?1.0:0.0))
                goto done;
        }
    }
    failed=0;
done:
    free(s); free(dv); free(uv); free(lv);
    free(ui); free(li); free(up); free(lp);
    return failed;
}

static int check_symm_and_norm(int nblocks)
{
    struct LIS_MATRIX_STRUCT a;
    struct LIS_MATRIX_CORE_STRUCT lower,upper;
    struct LIS_MATRIX_DIAG_STRUCT diag;
    LIS_INT *ptr=NULL,*cols=NULL;
    LIS_SCALAR *lv=NULL,*uv=NULL,*dv=NULL,*scale=NULL;
    LIS_SCALAR norm;
    double expected_norm=0.0;
    int i,k,failed=1;
    memset(&a,0,sizeof(a));
    memset(&lower,0,sizeof(lower));
    memset(&upper,0,sizeof(upper));
    memset(&diag,0,sizeof(diag));
    ptr=(LIS_INT *)calloc((size_t)nblocks+1,sizeof(*ptr));
    cols=(LIS_INT *)calloc((size_t)nblocks,sizeof(*cols));
    lv=(LIS_SCALAR *)calloc((size_t)nblocks*9,sizeof(*lv));
    uv=(LIS_SCALAR *)calloc((size_t)nblocks*9,sizeof(*uv));
    dv=(LIS_SCALAR *)calloc((size_t)nblocks*9,sizeof(*dv));
    scale=(LIS_SCALAR *)calloc((size_t)nblocks*3,sizeof(*scale));
    if(!ptr || !cols || !lv || !uv || !dv || !scale) goto done;
    for(i=0;i<nblocks;i++)
    {
        ptr[i]=i;
        cols[i]=(i*17+1)%nblocks;
        for(k=0;k<9;k++)
        {
            lv[9*i+k]=(LIS_SCALAR)(k+1);
            uv[9*i+k]=(LIS_SCALAR)(k+2);
            dv[9*i+k]=(LIS_SCALAR)(k+3);
            expected_norm+=(k+1)*(k+1)+(k+2)*(k+2)+(k+3)*(k+3);
        }
        scale[3*i+0]=2+(i%3);
        scale[3*i+1]=3+(i%4);
        scale[3*i+2]=4+(i%5);
    }
    ptr[nblocks]=nblocks;
    lower.bptr=upper.bptr=ptr;
    lower.bindex=upper.bindex=cols;
    lower.value=lv; upper.value=uv; diag.value=dv;
    a.nr=nblocks; a.bnr=a.bnc=3; a.is_splited=1;
    a.L=&lower; a.U=&upper; a.D=&diag;
    if(lis_matrix_normf_bsr(&a,&norm)!=LIS_SUCCESS ||
       !close_to(norm,sqrt(expected_norm)))
    {
        fprintf(stderr,"Split BSR norm mismatch\n");
        goto done;
    }
    if(lis_matrix_scale_symm_bsr(&a,scale)!=LIS_SUCCESS) goto done;
    for(i=0;i<nblocks;i++)
    {
        int col=cols[i];
        for(k=0;k<9;k++)
        {
            int row=k%3,inner_col=k/3;
            double diag_expected=(k+3)*scale[3*i+row]*scale[3*i+inner_col];
            double l_expected=(k+1)*scale[3*i+row]*scale[3*col+inner_col];
            double u_expected=(k+2)*scale[3*i+row]*scale[3*col+inner_col];
            if(!close_to(dv[9*i+k],diag_expected) ||
               !close_to(lv[9*i+k],l_expected) ||
               !close_to(uv[9*i+k],u_expected))
            {
                fprintf(stderr,"Symmetric BSR scale mismatch block=%d entry=%d\n",i,k);
                goto done;
            }
        }
    }
    failed=0;
done:
    free(scale); free(dv); free(uv); free(lv); free(cols); free(ptr);
    return failed;
}

static int check_unsplit_norm(void)
{
    struct LIS_MATRIX_STRUCT a;
    LIS_INT ptr[3]={0,1,2};
    LIS_SCALAR values[18];
    LIS_SCALAR norm;
    double expected=0.0;
    int i;
    memset(&a,0,sizeof(a));
    a.nr=2; a.bnr=a.bnc=3; a.is_splited=0;
    a.bptr=ptr; a.value=values;
    for(i=0;i<18;i++)
    {
        values[i]=(LIS_SCALAR)(i+1);
        expected+=(i+1)*(i+1);
    }
    return lis_matrix_normf_bsr(&a,&norm)!=LIS_SUCCESS ||
           !close_to(norm,sqrt(expected));
}

int main(int argc,char **argv)
{
    int failed=0;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 2;
    failed|=check_block_scaling(512);
    failed|=check_symm_and_norm(512);
    failed|=check_unsplit_norm();
    lis_finalize();
    if(failed) { fprintf(stderr,"BSR audit regression FAIL\n"); return 1; }
    puts("BSR audit regression PASS");
    return 0;
}
