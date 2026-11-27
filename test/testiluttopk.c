#include "lislib.h"
#include <stdio.h>
#include <stdlib.h>

static int retained(LIS_MATRIX_ILU f,LIS_INT row,LIS_INT want)
{
    if(!f || f->nnz[row]!=1 || f->index[row][0]!=want)
    {
        fprintf(stderr,"row=%lld want=%lld count=%lld first=%lld\n",
            (long long)row,(long long)want,
            (long long)(f?f->nnz[row]:-1),
            (long long)(f && f->nnz[row]>0?f->index[row][0]:-1));
        return 0;
    }
    return 1;
}

/* Triangular fixtures isolate the first U and last L rows. There are
 * sixteen entries across eight rows, so lfil=1 in serial and OpenMP. */
static int run(int lower)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER s=NULL;
    LIS_PRECON p=NULL;
    LIS_INT i,j,rc;
    int fail=1;
    const double weights[3]={-1.,3.,-7.};
    if(lis_matrix_create(LIS_COMM_WORLD,&A)||lis_matrix_set_size(A,0,8))goto done;
    for(i=0;i<8;i++)if(lis_matrix_set_value(LIS_INS_VALUE,i,i,32.,A))goto done;
    for(i=0;i<8;i+=4)
    {
        for(j=0;j<3;j++)
        {
            LIS_INT row=lower?i+3:i;
            LIS_INT col=lower?i+j:i+j+1;
            if(lis_matrix_set_value(LIS_INS_VALUE,row,col,weights[j],A))goto done;
        }
        if(lis_matrix_set_value(LIS_INS_VALUE,i+1,lower?i:i+2,0.25,A))goto done;
    }
    if(lis_matrix_set_type(A,LIS_MATRIX_CSR)||lis_matrix_assemble(A)||lis_solver_create(&s))goto done;
    s->A=A;
    if(lis_solver_set_option("-p ilut -ilutp 0 -ilut_drop 0 -ilut_rate 1 -scale none",s))goto done;
    rc=lis_precon_create(s,&p);
    if(rc){p=NULL;goto done;}
    if(lower)
    {
        if(!retained(p->L,3,2)||!retained(p->L,7,6))goto done;
    }
    else if(!retained(p->U,0,3)||!retained(p->U,4,7))goto done;
    fail=0;
done:
    if(p)lis_precon_destroy(p);
    if(s)lis_solver_destroy(s);
    if(A)lis_matrix_destroy(A);
    printf("ILUT_TOPK_%s=%s\n",lower?"LOWER":"UPPER",fail?"FAIL":"PASS");
    return fail;
}
int main(int argc,char **argv)
{
    int fail;
    if(lis_initialize(&argc,&argv))return 1;
    fail=run(0);fail|=run(1);
    lis_finalize();
    return fail;
}
