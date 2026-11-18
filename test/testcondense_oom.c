/* Deterministic test-only allocation failures in transactional Schur setup.
 * No global malloc overrides, no injection code in production compilation.
 */
#include "lis_condense_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef LIS_SC_TEST_ALLOC_FAIL
#error Test must compile its own copy of lis_condense_core.c with LIS_SC_TEST_ALLOC_FAIL
#endif

static size_t rp[]={0,2,3,6,8};
static size_t ci[]={0,2,1,0,2,3,2,3};
static double old_values[]={4,0.5,5,0.25,3,0.35,0.4,6};
static double new_values[]={4,0.75,5,0.25,3.125,0.35,0.4,6};
static size_t bp[]={0,1};
static size_t qd[]={2};
static lis_sc_blocks blocks={1,bp,qd};
static int checked=0;

#define CHECK(expr) do {checked++; if(!(expr)){ \
    fprintf(stderr,"OOM assertion %d at line %d: %s\n",checked,__LINE__,#expr);return 1; \
}} while(0)

static int equal_matrix(const lis_sc_csr *a,const lis_sc_csr *b){
    return a && b && a->n==b->n && a->nnz==b->nnz &&
        !memcmp(a->rowptr,b->rowptr,(a->n+1)*sizeof(size_t)) &&
        !memcmp(a->colind,b->colind,a->nnz*sizeof(size_t)) &&
        !memcmp(a->values,b->values,a->nnz*sizeof(double));
}

static int exhaust_fault_positions(lis_sc_plan **slot,const lis_sc_csr *input,
                                    lis_sc_reuse_mode mode){
    long failure,failures=0;
    for(failure=0;failure<300;failure++){
        lis_sc_plan *old=*slot;
        const lis_sc_csr *old_s=lis_sc_matrix(old);
        double *before=(double*)malloc(old_s->nnz*sizeof(double));
        lis_sc_status st;
        if(!before){fprintf(stderr,"test fixture OOM\n");return 1;}
        memcpy(before,old_s->values,old_s->nnz*sizeof(double));
        lis_sc_test_fail_after(failure);
        st=lis_sc_build(slot,input,&blocks,mode,1e-13);
        lis_sc_test_fail_after(-1);
        if(st==LIS_SC_OK){
            free(before);
            CHECK(*slot!=old);
            CHECK(failures>=8);
            printf("OOM %s: %ld allocation failures rejected atomically\n",
                   mode==LIS_SC_FULL?"FULL":"NUMERIC",failures);
            return 0;
        }
        CHECK(st==LIS_SC_OUT_OF_MEMORY);
        CHECK(*slot==old);
        CHECK(!memcmp(lis_sc_matrix(*slot)->values,before,old_s->nnz*sizeof(double)));
        failures++;
        free(before);
    }
    fprintf(stderr,"Failed to exhaust allocator within 300 positions\n");
    return 1;
}

int main(void){
    lis_sc_csr A={4,8,rp,ci,old_values};
    lis_sc_csr B={4,8,rp,ci,new_values};
    lis_sc_plan *slot=NULL,*ref=NULL;
    lis_sc_test_fail_after(-1);
    CHECK(lis_sc_build(&slot,&A,&blocks,LIS_SC_FULL,1e-13)==LIS_SC_OK);
    CHECK(slot!=NULL);
    CHECK(exhaust_fault_positions(&slot,&B,LIS_SC_NUMERIC)==0);
    CHECK(lis_sc_build(&ref,&B,&blocks,LIS_SC_FULL,1e-13)==LIS_SC_OK);
    CHECK(equal_matrix(lis_sc_matrix(slot),lis_sc_matrix(ref)));
    CHECK(exhaust_fault_positions(&slot,&A,LIS_SC_FULL)==0);
    lis_sc_free(ref);ref=NULL;
    CHECK(lis_sc_build(&ref,&A,&blocks,LIS_SC_FULL,1e-13)==LIS_SC_OK);
    CHECK(equal_matrix(lis_sc_matrix(slot),lis_sc_matrix(ref)));
    lis_sc_free(ref);lis_sc_free(slot);
    printf("condensation allocator faults: %d assertions, 0 failures\n",checked);
    return 0;
}
