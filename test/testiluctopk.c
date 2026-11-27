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

/* Two disjoint four-node stars; each worker owns a full star with two
 * workers. The largest magnitude is at a noncontiguous candidate index.
 * Blocks are scalar multiples of I, so the same ordering holds for BSR. */
static int run(int bs)
{
    LIS_MATRIX A=NULL,B=NULL;
    LIS_SOLVER s=NULL;
    LIS_PRECON p=NULL;
    LIS_INT i,j,k,rc;
    int fail=1;
    const double weights[3]={-1.,3.,-7.};
    if(lis_matrix_create(LIS_COMM_WORLD,&A) ||
       lis_matrix_set_size(A,0,8*bs)) goto done;
    for(i=0;i<8*bs;i++)
        if(lis_matrix_set_value(LIS_INS_VALUE,i,i,32.,A)) goto done;
    for(i=0;i<8;i+=4)
        for(j=1;j<4;j++)
            for(k=0;k<bs;k++)
            {
                if(lis_matrix_set_value(LIS_INS_VALUE,i*bs+k,(i+j)*bs+k,weights[j-1],A) ||
                   lis_matrix_set_value(LIS_INS_VALUE,(i+j)*bs+k,i*bs+k,weights[j-1],A)) goto done;
            }
    if(lis_matrix_set_type(A,LIS_MATRIX_CSR)||lis_matrix_assemble(A)) goto done;
    if(bs>1)
    {
        if(lis_matrix_duplicate(A,&B)||lis_matrix_set_blocksize(B,bs,bs,NULL,NULL)||
           lis_matrix_set_type(B,LIS_MATRIX_BSR)||lis_matrix_convert(A,B)) goto done;
    }
    if(lis_solver_create(&s)) goto done;
    s->A=B?B:A;
    if(lis_solver_set_option(bs==1?
       "-p iluc -iluc_drop 0 -iluc_rate 1 -scale none":
       bs==2?"-p iluc -iluc_drop 0 -iluc_rate 1 -scale none -storage bsr -storage_block 2":
             "-p iluc -iluc_drop 0 -iluc_rate 1 -scale none -storage bsr -storage_block 3",s)) goto done;
    rc=lis_precon_create(s,&p);
    if(rc){p=NULL;goto done;}
    if(!retained(p->U,0,3)||!retained(p->L,0,3)||
       !retained(p->U,4,7)||!retained(p->L,4,7)) goto done;
    fail=0;
done:
    if(p)lis_precon_destroy(p);
    if(s)lis_solver_destroy(s);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    printf("ILUC_TOPK_BLOCK_%d=%s\n",bs,fail?"FAIL":"PASS");
    return fail;
}
int main(int argc,char **argv)
{
    int fail;
    if(lis_initialize(&argc,&argv))return 1;
    /* BSR uses a different factor representation; run after CSR. */
    fail=run(1);fail|=run(2);fail|=run(3);
    lis_finalize();
    return fail;
}
