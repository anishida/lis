#include "lis_camfamg_composite_core.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
static int tests=0,errors=0;
#define CHECK(x) do{tests++;if(!(x)){fprintf(stderr,"FAIL: %s at %d\n",#x,__LINE__);errors++;}}while(0)
typedef struct {double a00,a01,a10,a11;int calls;} matrix;
static int matvec(void *v,const double*r,double*z){matrix*m=(matrix*)v;m->calls++;z[0]=m->a00*r[0]+m->a01*r[1];z[1]=m->a10*r[0]+m->a11*r[1];return 0;}
static int ca(void *v,const double*r,double*z){int *n=(int*)v;(*n)++;z[0]=r[0]/4.0;z[1]=r[1]/4.0;return 0;}
static int ilut(void*v,const double*r,double*z){int*n=(int*)v;(*n)++;z[0]=r[0]/3.0;z[1]=r[1]/5.0;return 0;}
static int fails(void*v,const double*r,double*z){(void)v;(void)r;(void)z;return 1;}
int main(void){
    matrix m={4,1,1,3,0};int ca_calls=0,ilut_calls=0;
    lis_cam_composite op={2,matvec,&m,ca,&ca_calls,ilut,&ilut_calls,0};
    double r[]={1.5,2.5},z[]={1234,5678};
    CHECK(lis_cam_apply_multiplicative(&op,r,z)==LIS_CAM_OK);
    CHECK(ca_calls==1&&ilut_calls==1&&m.calls==1);
    /* Explicit two-stage formula. */
    double z1[]={r[0]/4,r[1]/4};double az[]={4*z1[0]+z1[1],z1[0]+3*z1[1]};
    CHECK(fabs(z[0]-(z1[0]+(r[0]-az[0])/3.0))<1e-14);
    CHECK(fabs(z[1]-(z1[1]+(r[1]-az[1])/5.0))<1e-14);
    op.safe_damping=1; CHECK(lis_cam_apply_multiplicative(&op,r,z)==LIS_CAM_OK);
    CHECK(isfinite(z[0])&&isfinite(z[1]));
    op.ca_apply=fails;z[0]=11;z[1]=12;
    CHECK(lis_cam_apply_multiplicative(&op,r,z)==LIS_CAM_CALLBACK_ERROR);
    CHECK(z[0]==11&&z[1]==12); /* atomic output on callback failure */
    printf("composite core: %d assertions, %d failures\n",tests,errors);
    return errors?1:0;
}
