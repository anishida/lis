/* Generalized block eigensolver metric with a massless degree of freedom.
 * B=diag(2,0,5) has rank 2; the massless direction must be deflated.
 */
#include "lis.h"
#include "lis_esolver.h"
#include <math.h>
#include <stdio.h>

int main(int argc,char **argv)
{
    LIS_MATRIX B=NULL;
    LIS_VECTOR W[3]={NULL,NULL,NULL},BW[3]={NULL,NULL,NULL};
    LIS_VECTOR Q[3]={NULL,NULL,NULL},BQ[3]={NULL,NULL,NULL};
    LIS_SCALAR gram[9],eval[3],evec[9],transform[9];
    LIS_INT i,j,k,rank=-1,err=0;
    int failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 2;
    err=lis_matrix_create(LIS_COMM_WORLD,&B);if(err)goto done;
    err=lis_matrix_set_size(B,0,3);if(err)goto done;
    err=lis_matrix_set_value(LIS_INS_VALUE,0,0,2.0,B);if(err)goto done;
    err=lis_matrix_set_value(LIS_INS_VALUE,2,2,5.0,B);if(err)goto done;
    err=lis_matrix_assemble(B);if(err)goto done;
    for(i=0;i<3;i++)
    {
        err=lis_vector_duplicate(B,&W[i]);if(err)goto done;
        err=lis_vector_duplicate(B,&BW[i]);if(err)goto done;
        err=lis_vector_duplicate(B,&Q[i]);if(err)goto done;
        err=lis_vector_duplicate(B,&BQ[i]);if(err)goto done;
        err=lis_vector_set_all(0.0,W[i]);if(err)goto done;
        W[i]->value[i]=1.0;
    }
    err=lis_eigen_block_brank_deflate(
        B,3,3,W,BW,Q,BQ,1e-12,&rank,gram,eval,evec,transform);
    if(err || rank!=2)
    {
        fprintf(stderr,"SEMIDEFINITE_B rank=%d err=%d (expected 2)\n",
                (int)rank,(int)err);
        goto done;
    }
    for(i=0;i<rank;i++)
    {
        for(j=0;j<rank;j++)
        {
            LIS_SCALAR value=0.0;
            for(k=0;k<3;k++)
                value+=Q[i]->value[k]*BQ[j]->value[k];
            if(!isfinite((double)value) ||
               fabs((double)value-(i==j?1.0:0.0))>1e-9)
            {
                fprintf(stderr,"SEMIDEFINITE_B B-orthogonality failure\n");
                goto done;
            }
        }
    }
    failed=0;
done:
    for(i=0;i<3;i++)
    {
        if(BQ[i])lis_vector_destroy(BQ[i]);
        if(Q[i])lis_vector_destroy(Q[i]);
        if(BW[i])lis_vector_destroy(BW[i]);
        if(W[i])lis_vector_destroy(W[i]);
    }
    if(B)lis_matrix_destroy(B);
    lis_finalize();
    if(failed)return 1;
    puts("SEMIDEFINITE_B B-orthogonality PASS");
    return 0;
}
