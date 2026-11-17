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
#include <string.h>

#define EXPECT(expr) do { \
    if(!(expr)) { \
        fprintf(stderr,"testrcoalloc:%d: %s failed\n",__LINE__,#expr); \
        failed=1; \
    } \
} while(0)

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
    LIS_INT *row=NULL,**indices=NULL,**null_indices=NULL;
    LIS_SCALAR **values=NULL,**null_values=NULL;
    unsigned char saved_scalar[sizeof(LIS_SCALAR)];
    LIS_INT err;
    int failed=0;

    EXPECT(lis_matrix_malloc_rco(1,NULL,NULL,&indices,&values)
           ==LIS_ERR_ILL_ARG);
    EXPECT(lis_matrix_malloc_rco(1,NULL,&row,NULL,&values)
           ==LIS_ERR_ILL_ARG);
    EXPECT(lis_matrix_malloc_rco(1,NULL,&row,&indices,NULL)
           ==LIS_ERR_ILL_ARG);

    err=lis_matrix_malloc_rco(-1,NULL,&row,&indices,&values);
    EXPECT(err==LIS_ERR_ILL_ARG &&
           row==NULL && indices==NULL && values==NULL);
    err=lis_matrix_malloc_rco(3,invalid,&row,&indices,&values);
    EXPECT(err==LIS_ERR_ILL_ARG &&
           row==NULL && indices==NULL && values==NULL);

    err=lis_matrix_malloc_rco(3,counts,&row,&indices,&values);
    if(err!=LIS_SUCCESS || row==NULL || indices==NULL || values==NULL)
    {
        fprintf(stderr,"testrcoalloc: positive allocation failed: %d\n",(int)err);
        cleanup(3,row,indices,values);
        return 1;
    }
    EXPECT(indices[0]!=NULL && values[0]!=NULL &&
           indices[1]==NULL && values[1]==NULL &&
           indices[2]!=NULL && values[2]!=NULL &&
           row[0]==0 && row[1]==0 && row[2]==0);
    if(indices[0]!=NULL && values[0]!=NULL)
    {
        indices[0][0]=17;
        memset(&values[0][0],0x5a,sizeof(LIS_SCALAR));
        memcpy(saved_scalar,&values[0][0],sizeof(saved_scalar));

        err=lis_matrix_realloc_rco(0,-1,&indices,&values);
        EXPECT(err==LIS_ERR_ILL_ARG && indices[0][0]==17);
        err=lis_matrix_realloc_rco(0,0,&indices,&values);
        EXPECT(err==LIS_ERR_ILL_ARG && indices[0][0]==17);
        err=lis_matrix_realloc_rco(-1,2,&indices,&values);
        EXPECT(err==LIS_ERR_ILL_ARG && indices[0][0]==17);
        err=lis_matrix_realloc_rco(0,2,NULL,&values);
        EXPECT(err==LIS_ERR_ILL_ARG);
        err=lis_matrix_realloc_rco(0,2,&indices,NULL);
        EXPECT(err==LIS_ERR_ILL_ARG);
        err=lis_matrix_realloc_rco(0,2,&null_indices,&values);
        EXPECT(err==LIS_ERR_ILL_ARG);
        err=lis_matrix_realloc_rco(0,2,&indices,&null_values);
        EXPECT(err==LIS_ERR_ILL_ARG);

        EXPECT(memcmp(&values[0][0],saved_scalar,sizeof(saved_scalar))==0);

        err=lis_matrix_realloc_rco(0,3,&indices,&values);
        EXPECT(err==LIS_SUCCESS && indices[0]!=NULL && values[0]!=NULL);
        if(err==LIS_SUCCESS)
        {
            EXPECT(indices[0][0]==17);
            EXPECT(memcmp(&values[0][0],saved_scalar,
                          sizeof(saved_scalar))==0);
        }
    }

    /* A row with zero initial capacity must grow from NULL pointers. */
    err=lis_matrix_realloc_rco(1,2,&indices,&values);
    EXPECT(err==LIS_SUCCESS && indices[1]!=NULL && values[1]!=NULL);
    if(err==LIS_SUCCESS)
    {
        indices[1][0]=23;
        EXPECT(indices[1][0]==23);
    }

    cleanup(3,row,indices,values);
    return failed;
}

static int check_size_overflow(void)
{
    /* Compute signed LIS_INT_MAX without assuming long or long long. */
    const uintmax_t unsigned_limit=
        ((uintmax_t)1 << (sizeof(LIS_INT)*CHAR_BIT-1))-1;
    const LIS_INT limit=(LIS_INT)unsigned_limit;
    LIS_INT counts[1]={1},huge[1];
    LIS_INT *row=NULL,**indices=NULL;
    LIS_SCALAR **values=NULL;
    LIS_INT *other_row=NULL,**other_indices=NULL;
    LIS_SCALAR **other_values=NULL;
    LIS_INT err;
    int failed=0;
    int tested=0;

    err=lis_matrix_malloc_rco(1,counts,&row,&indices,&values);
    if(err!=LIS_SUCCESS || row==NULL || indices==NULL || values==NULL)
    {
        fprintf(stderr,"testrcoalloc: overflow fixture allocation failed\n");
        cleanup(1,row,indices,values);
        return 1;
    }
    indices[0][0]=17;
    if(unsigned_limit > (uintmax_t)SIZE_MAX/sizeof(LIS_SCALAR) ||
       unsigned_limit > (uintmax_t)SIZE_MAX/sizeof(LIS_INT))
    {
        LIS_INT *old_index=indices[0];
        LIS_SCALAR *old_value=values[0];
        err=lis_matrix_realloc_rco(0,limit,&indices,&values);
        EXPECT(err==LIS_OUT_OF_MEMORY);
        EXPECT(indices[0]==old_index && values[0]==old_value);
        EXPECT(indices[0][0]==17);
        tested=1;

        huge[0]=limit;
        err=lis_matrix_malloc_rco(1,huge,&other_row,
                                  &other_indices,&other_values);
        EXPECT(err==LIS_OUT_OF_MEMORY &&
               other_row==NULL &&
               other_indices==NULL &&
               other_values==NULL);
    }

    if(unsigned_limit > (uintmax_t)SIZE_MAX/sizeof(LIS_INT *) ||
       unsigned_limit > (uintmax_t)SIZE_MAX/sizeof(LIS_SCALAR *))
    {
        err=lis_matrix_malloc_rco(limit,NULL,&other_row,
                                  &other_indices,&other_values);
        EXPECT(err==LIS_OUT_OF_MEMORY &&
               other_row==NULL &&
               other_indices==NULL &&
               other_values==NULL);
        tested=1;
    }

    if(!tested)
        puts("RCO_ALLOC_OVERFLOW: SKIPPED (LIS_INT range fits size_t)");
    else
        puts("RCO_ALLOC_OVERFLOW: exercised real size_t guards");

    cleanup(1,row,indices,values);
    return failed;
}

int main(int argc,char **argv)
{
    LIS_INT err=lis_initialize(&argc,&argv);
    int failed=0;
    if(err) return 1;
    failed |= check_inputs_and_rows();
    failed |= check_size_overflow();
    if(lis_finalize()!=LIS_SUCCESS) failed=1;
    if(failed) fprintf(stderr,"RCO_ALLOC_SAFETY=FAILED\n");
    else printf("RCO_ALLOC_SAFETY=PASSED\n");
    return failed;
}
