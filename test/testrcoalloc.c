/* Generic RCO allocation-size validation and storage ownership regression. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include "lis.h"
#include "lis_matrix.h"
#include <stdint.h>
#include <limits.h>
#include <stdio.h>

static void cleanup(LIS_INT n,LIS_INT *row,LIS_INT **index,LIS_SCALAR **value)
{
    LIS_INT i;
    if(index!=NULL || value!=NULL)
        for(i=0;i<n;i++)
        {
            if(index!=NULL && index[i]!=NULL) lis_free(index[i]);
            if(value!=NULL && value[i]!=NULL) lis_free(value[i]);
        }
    if(row!=NULL) lis_free(row);
    if(index!=NULL) lis_free(index);
    if(value!=NULL) lis_free(value);
}

static int check_inputs_and_rows(void)
{
    LIS_INT counts[3]={1,0,2},invalid[3]={1,-1,2};
    LIS_INT *row=NULL,**indices=NULL;
    LIS_SCALAR **values=NULL;
    LIS_INT err;
    int failed=0;

    err=lis_matrix_malloc_rco(-1,NULL,&row,&indices,&values);
    if(err!=LIS_ERR_ILL_ARG || row!=NULL || indices!=NULL || values!=NULL)
        failed=1;
    err=lis_matrix_malloc_rco(3,invalid,&row,&indices,&values);
    if(err!=LIS_ERR_ILL_ARG || row!=NULL || indices!=NULL || values!=NULL)
        failed=1;

    err=lis_matrix_malloc_rco(3,counts,&row,&indices,&values);
    if(err!=LIS_SUCCESS || row==NULL || indices==NULL || values==NULL)
    {
        fprintf(stderr,"RCO positive allocation failed: %d\n",(int)err);
        /* When allocation failed, ownership is retained by the allocator. */
        return 1;
    }
    if(indices[0]==NULL || values[0]==NULL ||
       indices[1]!=NULL || values[1]!=NULL ||
       indices[2]==NULL || values[2]==NULL ||
       row[0]!=0 || row[1]!=0 || row[2]!=0)
        failed=1;

    indices[0][0]=17;
    err=lis_matrix_realloc_rco(0,-1,&indices,&values);
    if(err!=LIS_ERR_ILL_ARG || indices[0][0]!=17) failed=1;
    err=lis_matrix_realloc_rco(0,0,&indices,&values);
    if(err!=LIS_ERR_ILL_ARG || indices[0][0]!=17) failed=1;
    err=lis_matrix_realloc_rco(-1,2,&indices,&values);
    if(err!=LIS_ERR_ILL_ARG || indices[0][0]!=17) failed=1;

    err=lis_matrix_realloc_rco(0,3,&indices,&values);
    if(err!=LIS_SUCCESS || indices[0]==NULL || values[0]==NULL)
        failed=1;
    else if(indices[0][0]!=17) failed=1;

    cleanup(3,row,indices,values);
    return failed;
}

static int check_size_overflow(void)
{
    LIS_INT limit;
    LIS_INT huge[1];
    LIS_INT *row=NULL,**indices=NULL;
    LIS_SCALAR **values=NULL;
    LIS_INT err;
    int failed=0;

    limit=(sizeof(LIS_INT)==sizeof(long long))
          ?(LIS_INT)LLONG_MAX:(LIS_INT)INT_MAX;

    /* On common 32-bit LIS_INT / 64-bit size_t builds, valid LIS_INT
     * values cannot overflow size_t; avoid requesting huge allocations.
     */
    if((uintmax_t)limit > (uintmax_t)SIZE_MAX / sizeof(LIS_SCALAR))
    {
        huge[0]=limit;
        err=lis_matrix_malloc_rco(1,huge,&row,&indices,&values);
        if(err!=LIS_OUT_OF_MEMORY || row!=NULL ||
           indices!=NULL || values!=NULL)
            failed=1;
        err=lis_matrix_realloc_rco(0,limit,&indices,&values);
        /* Invalid pointer arrays are reported before sizing checks. */
        if(err!=LIS_ERR_ILL_ARG) failed=1;
    }
    if((uintmax_t)limit > (uintmax_t)SIZE_MAX / sizeof(LIS_INT *))
    {
        err=lis_matrix_malloc_rco(limit,NULL,&row,&indices,&values);
        if(err!=LIS_OUT_OF_MEMORY || row!=NULL ||
           indices!=NULL || values!=NULL)
            failed=1;
    }
    return failed;
}

int main(int argc,char **argv)
{
    LIS_INT err=lis_initialize(&argc,&argv);
    int failed;
    if(err) return 1;
    failed=check_inputs_and_rows() || check_size_overflow();
    lis_finalize();
    if(failed) fprintf(stderr,"RCO_ALLOC_SAFETY=FAILED\n");
    else printf("RCO_ALLOC_SAFETY=PASSED\n");
    return failed;
}
