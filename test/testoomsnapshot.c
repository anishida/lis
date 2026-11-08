/* Native OOM and cleanup regression for the opt-in generic snapshot APIs.
 * This is intentionally limited to file lifecycle and ownership safety.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_snapshot.h"
#include <stdio.h>
#include <string.h>

#ifndef LIS_TEST_OOM
#error "testoomsnapshot requires --enable-oom-tests"
#endif

#define CHECK(expr) do { \
    if(!(expr)) { \
        fprintf(stderr,"testoomsnapshot:%d: %s failed\n",__LINE__,#expr); \
        return 1; \
    } \
} while(0)

#define PREFIX "lis-oom-snapshot-regression"

static const char *const outputs[] = {
    PREFIX "_A.mtx", PREFIX "_A.mtx.tmp",
    PREFIX "_b.mtx", PREFIX "_b.mtx.tmp",
    PREFIX "_dofs.csv", PREFIX "_dofs.csv.tmp"
};

static int no_outputs(void)
{
    size_t i;
    FILE *f;
    for(i=0;i<sizeof(outputs)/sizeof(outputs[0]);i++)
    {
        f=fopen(outputs[i],"rb");
        if(f!=NULL)
        {
            fclose(f);
            fprintf(stderr,"Unexpected snapshot file: %s\n",outputs[i]);
            return 0;
        }
    }
    return 1;
}

static LIS_INT label_equation(void *context,LIS_INT eq,
                              const char **node,const char **dof,
                              const char **body)
{
    (void)context;
    if(eq<0 || eq>1) return LIS_ERR_ILL_ARG;
    *node=(eq==0)?"node,1":"node-2";
    *dof=(eq==0)?"UX":"TEMP";
    *body=NULL;
    return LIS_SUCCESS;
}

static LIS_INT fail_second_label(void *context,LIS_INT eq,
                                 const char **node,const char **dof,
                                 const char **body)
{
    if(eq==1) return LIS_ERR_ILL_ARG;
    return label_equation(context,eq,node,dof,body);
}

static int check_csv(LIS_MATRIX A,size_t baseline)
{
    size_t fail_at;
    LIS_INT err;
    FILE *f;

    CHECK(no_outputs());
    for(fail_at=1;fail_at<=2;fail_at++)
    {
        lis_test_oom_begin(fail_at);
        err=lis_snapshot_export_dof_map_csv(A,label_equation,NULL,PREFIX);
        lis_test_oom_end();
        CHECK(err==LIS_OUT_OF_MEMORY);
        /* The function allocates both buffers before checking either;
         * cleanup paths may make further tracked allocation attempts.
         */
        CHECK(lis_test_oom_calls()>=fail_at);
        CHECK(lis_test_oom_live()==baseline);
        CHECK(no_outputs());
    }

    /* Failure after one callback must remove the partially written CSV. */
    err=lis_snapshot_export_dof_map_csv(A,fail_second_label,NULL,PREFIX);
    CHECK(err==LIS_ERR_ILL_ARG);
    CHECK(lis_test_oom_live()==baseline);
    CHECK(no_outputs());

    lis_test_oom_begin(0);
    err=lis_snapshot_export_dof_map_csv(A,label_equation,NULL,PREFIX);
    lis_test_oom_end();
    CHECK(err==LIS_SUCCESS);
    CHECK(lis_test_oom_calls()==2);
    CHECK(lis_test_oom_live()==baseline);
    f=fopen(PREFIX "_dofs.csv","rb");
    CHECK(f!=NULL);
    CHECK(fclose(f)==0);
    CHECK(remove(PREFIX "_dofs.csv")==0);
    CHECK(no_outputs());
    return 0;
}

#ifndef _COMPLEX
static int check_mm(LIS_MATRIX A,LIS_VECTOR b,size_t baseline)
{
    size_t fail_at,allocation_count;
    LIS_INT err;
    FILE *f;

    CHECK(no_outputs());

    /* Probe the actual number of LIS allocations on the successful path.
     * Explicit and convertible CSR inputs may require different counts.
     */
    lis_test_oom_begin(0);
    err=lis_snapshot_export_mm(A,b,NULL,NULL,0,NULL,PREFIX);
    lis_test_oom_end();
    allocation_count=lis_test_oom_calls();
    CHECK(err==LIS_SUCCESS);
    CHECK(allocation_count>=2);
    CHECK(lis_test_oom_live()==baseline);
    f=fopen(PREFIX "_A.mtx","rb");
    CHECK(f!=NULL);
    CHECK(fclose(f)==0);
    f=fopen(PREFIX "_b.mtx","rb");
    CHECK(f!=NULL);
    CHECK(fclose(f)==0);
    CHECK(remove(PREFIX "_A.mtx")==0);
    CHECK(remove(PREFIX "_b.mtx")==0);
    CHECK(no_outputs());

    for(fail_at=1;fail_at<=allocation_count;fail_at++)
    {
        size_t calls;
        lis_test_oom_begin(fail_at);
        err=lis_snapshot_export_mm(A,b,NULL,NULL,0,NULL,PREFIX);
        lis_test_oom_end();
        calls=lis_test_oom_calls();
        if(err==LIS_SUCCESS || lis_test_oom_live()!=baseline ||
           !no_outputs())
        {
            fprintf(stderr,
                    "MM OOM injection fail_at=%lu, calls=%lu, rc=%d, "
                    "live=%lu baseline=%lu\n",
                    (unsigned long)fail_at,(unsigned long)calls,(int)err,
                    (unsigned long)lis_test_oom_live(),
                    (unsigned long)baseline);
            return 1;
        }
    }
    return 0;
}
#endif

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL;
    LIS_INT err;
    size_t baseline;
    int failed=0;

    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    CHECK(no_outputs());
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) { failed=1; goto done; }
    err=lis_matrix_set_size(A,0,2);
    if(err) { failed=1; goto done; }
    err=lis_matrix_set_value(LIS_INS_VALUE,0,0,(LIS_SCALAR)4,A);
    if(err) { failed=1; goto done; }
    err=lis_matrix_set_value(LIS_INS_VALUE,1,1,(LIS_SCALAR)3,A);
    if(err) { failed=1; goto done; }
    err=lis_matrix_assemble(A);
    if(err) { failed=1; goto done; }
    err=lis_vector_duplicate(A,&b);
    if(err) { failed=1; goto done; }
    err=lis_vector_set_all((LIS_SCALAR)1,b);
    if(err) { failed=1; goto done; }
    baseline=lis_test_oom_live();

#ifndef _COMPLEX
    if(check_mm(A,b,baseline)) failed=1;
#endif
    if(!failed && check_csv(A,baseline)) failed=1;
    if(lis_test_oom_live()!=baseline) failed=1;

done:
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    if(lis_finalize()!=LIS_SUCCESS) failed=1;
    if(failed) fprintf(stderr,"SNAPSHOT_NATIVE_OOM_REGRESSION=FAILED\n");
    else printf("SNAPSHOT_NATIVE_OOM_REGRESSION=PASSED\n");
    return failed;
}
