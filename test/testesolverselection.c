/* Focused regression for converged eigensolver-pair selection. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lis.h"
#include "lis_esolver_selection.h"

static int fail(const char *msg, int err)
{
    fprintf(stderr,"FAIL: %s (err=%d)\n",msg,err);
    return 1;
}

int main(int argc, char *argv[])
{
    LIS_ESOLVER esolver = NULL;
    LIS_ESOLVER_SELECTED_PAIR pairs[4];
    LIS_INT available = 0;
    LIS_INT err;

    err = lis_initialize(&argc,&argv);
    if(err) return fail("lis_initialize",err);

    err = lis_esolver_create(&esolver);
    if(err) return fail("lis_esolver_create",err);

    esolver->options[LIS_EOPTIONS_ESOLVER] = LIS_ESOLVER_LI;
    esolver->options[LIS_EOPTIONS_SUBSPACE] = 4;
    esolver->nev = 4;

    esolver->evalue = (LIS_SCALAR *)calloc(6,sizeof(LIS_SCALAR));
    esolver->resid = (LIS_REAL *)calloc(6,sizeof(LIS_REAL));
    if(esolver->evalue==NULL || esolver->resid==NULL)
        return fail("allocation",LIS_OUT_OF_MEMORY);

    esolver->evalue[0] = (LIS_SCALAR)3.0;
    esolver->evalue[1] = (LIS_SCALAR)(-1.0);
    esolver->evalue[2] = (LIS_SCALAR)2.0;
    esolver->evalue[3] = (LIS_SCALAR)1.0;

    esolver->resid[0] = (LIS_REAL)1.0e-8;
    esolver->resid[1] = (LIS_REAL)1.0e-8;
    esolver->resid[2] = (LIS_REAL)1.0e-3;
    esolver->resid[3] = (LIS_REAL)1.0e-9;

    err = lis_esolver_select_converged_pairs(
        esolver,4,2,(LIS_REAL)1.0e-6,
        LIS_ESOLVER_SELECT_ASCENDING_REAL,
        pairs,&available);
    if(err) return fail("algebraic selection",err);

    if(available!=3 ||
       (double)pairs[0].value!=-1.0 ||
       (double)pairs[1].value!=1.0 ||
       pairs[0].mode!=1 || pairs[1].mode!=3)
        return fail("algebraic ordering",LIS_BREAKDOWN);

    err = lis_esolver_select_converged_pairs(
        esolver,4,2,(LIS_REAL)1.0e-6,
        LIS_ESOLVER_SELECT_ASCENDING_REAL |
        LIS_ESOLVER_SELECT_POSITIVE_REAL,
        pairs,&available);
    if(err) return fail("positive selection",err);

    if(available!=2 ||
       (double)pairs[0].value!=1.0 ||
       (double)pairs[1].value!=3.0 ||
       pairs[0].mode!=3 || pairs[1].mode!=0)
        return fail("positive ordering",LIS_BREAKDOWN);

    available = -1;
    err = lis_esolver_select_converged_pairs(
        esolver,6,1,(LIS_REAL)1.0e-6,
        LIS_ESOLVER_SELECT_ASCENDING_REAL,
        pairs,&available);
    if(err==LIS_SUCCESS)
        return fail("out-of-range scan error propagation",LIS_BREAKDOWN);

    printf("PASS: eigensolver converged-pair selection\n");

    lis_esolver_destroy(esolver);
    lis_finalize();
    return 0;
}
