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

static double scalar_real(LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (double)creal(value);
#else
    return (double)value;
#endif
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_ESOLVER esolver=NULL;
    LIS_VECTOR out0=NULL,out1=NULL;
    LIS_VECTOR outputs[2];
    LIS_SCALAR values[2];
    LIS_INT available=-1,err,i;
    int rc=1;

    err=lis_initialize(&argc,&argv); if(err) return 1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A); if(err) goto done;
    err=lis_matrix_set_size(A,0,3); if(err) goto done;
    err=lis_esolver_create(&esolver); if(err) goto done;

    esolver->options[LIS_EOPTIONS_ESOLVER]=LIS_ESOLVER_LI;
    esolver->options[LIS_EOPTIONS_SUBSPACE]=4;
    esolver->nev=4;
    esolver->evalue=(LIS_SCALAR *)calloc(4,sizeof(LIS_SCALAR));
    esolver->resid=(LIS_REAL *)calloc(4,sizeof(LIS_REAL));
    esolver->evector=(LIS_VECTOR *)calloc(4,sizeof(LIS_VECTOR));
    if(!esolver->evalue || !esolver->resid || !esolver->evector) goto done;
    esolver->evector_owned_count=4;

    esolver->evalue[0]=(LIS_SCALAR)3.0;
    esolver->evalue[1]=(LIS_SCALAR)(-1.0);
    esolver->evalue[2]=(LIS_SCALAR)2.0;
    esolver->evalue[3]=(LIS_SCALAR)1.0;
    esolver->resid[0]=(LIS_REAL)1.0e-8;
    esolver->resid[1]=(LIS_REAL)1.0e-8;
    esolver->resid[2]=(LIS_REAL)1.0e-3;
    esolver->resid[3]=(LIS_REAL)1.0e-9;

    for(i=0;i<4;i++){
        err=lis_vector_duplicate(A,&esolver->evector[i]); if(err) goto done;
        lis_vector_set_all((LIS_SCALAR)(10+i),esolver->evector[i]);
    }
    err=lis_vector_duplicate(A,&out0); if(err) goto done;
    err=lis_vector_duplicate(A,&out1); if(err) goto done;
    outputs[0]=out0; outputs[1]=out1;

    err=lis_esolver_copy_selected_eigenpairs(
        esolver,4,2,(LIS_REAL)1.0e-6,
        LIS_ESOLVER_SELECT_ASCENDING_REAL,
        values,outputs,&available);
    if(err) goto done;

    if(available!=3 ||
       fabs(scalar_real(values[0])+1.0)>1.0e-12 ||
       fabs(scalar_real(values[1])-1.0)>1.0e-12 ||
       fabs(scalar_real(out0->value[0])-11.0)>1.0e-12 ||
       fabs(scalar_real(out1->value[0])-13.0)>1.0e-12)
        goto done;

    printf("LIS_ESOLVER_COPY_SELECTED PASSED\n");
    rc=0;

done:
    if(out1) lis_vector_destroy(out1);
    if(out0) lis_vector_destroy(out0);
    if(esolver) lis_esolver_destroy(esolver);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return rc;
}
