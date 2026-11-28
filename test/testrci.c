/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include "lis_rci.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"RCI test line %d: %s\n",__LINE__,#c); return 1; } } while (0)
int main(void)
{
    const double exact[3]={1,2,3};
    double b[3],x0[3]={0,0,0};
    lis_rci_request q;
    lis_rci_info info;
    CHECK(lis_rci_create(NULL)==NULL);
    lis_rci_options bad={SIZE_MAX,20,3,1e-12,0,LIS_RCI_CG};
    CHECK(lis_rci_create(&bad)==NULL);
    for(int m=LIS_RCI_CG;m<=LIS_RCI_FGMRES;m++){
        lis_rci_options o={3,100,3,1e-12,1e-14,(lis_rci_method)m};
        lis_rci_context *c=lis_rci_create(&o);
        CHECK(c!=NULL);
        CHECK(lis_rci_get_info(c,&info)==0 && info.status==LIS_RCI_STATE && isnan(info.residual));
        CHECK(lis_rci_advance(c,0,&q)==LIS_RCI_STATE);
        for(int rep=0;rep<2;rep++){
            /* SPD for CG; nonsymmetric for the other methods. */
            double a[3][3]={{4,1,0},{1,5,1},{0,1,6}};
            if(m!=LIS_RCI_CG)a[1][0]=-1;
            for(size_t i=0;i<3;i++){b[i]=0;for(size_t j=0;j<3;j++)b[i]+=a[i][j]*exact[j]*(rep+1);}
            int s=lis_rci_start(c,b,x0,&q),calls=0;
            CHECK(s==LIS_RCI_MATVEC);
            CHECK(lis_rci_start(c,b,x0,&q)==LIS_RCI_STATE);
            while(s>0){
                CHECK(q.n==3 && q.in && q.out && ++calls<1000);
                for(size_t i=0;i<3;i++){
                    if(s==LIS_RCI_MATVEC){q.out[i]=0;for(size_t j=0;j<3;j++)q.out[i]+=a[i][j]*q.in[j];}
                    else {CHECK(s==LIS_RCI_PRECON);q.out[i]=q.in[i]/(a[i][i]*(m==LIS_RCI_FGMRES?1+.1*((i+calls)%3):1));}
                }
                s=lis_rci_advance(c,0,&q);
            }
            CHECK(s==LIS_RCI_DONE && q.in==NULL && q.out==NULL && q.n==0);
            CHECK(lis_rci_get_info(c,&info)==0 && info.status==s && info.matvecs>0);
            const double *x=lis_rci_solution(c);
            double r2=0,b2=0;
            for(size_t i=0;i<3;i++){
                double r=b[i];for(size_t j=0;j<3;j++)r-=a[i][j]*x[j];
                r2+=r*r;b2+=b[i]*b[i];
                CHECK(fabs(x[i]-exact[i]*(rep+1))<1e-9);
            }
            CHECK(sqrt(r2)<=1.05e-12*sqrt(b2));
        }
        CHECK(lis_rci_start(c,b,NULL,&q)==LIS_RCI_MATVEC);
        CHECK(lis_rci_advance(c,73,&q)==LIS_RCI_CALLBACK);
        CHECK(lis_rci_get_info(c,&info)==0 && info.callback_code==73);
        CHECK(lis_rci_start(c,b,NULL,&q)==LIS_RCI_MATVEC);
        q.out[0]=0;
        CHECK(lis_rci_advance(c,0,&q)==LIS_RCI_NONFINITE);
        CHECK(lis_rci_start(c,b,NULL,&q)==LIS_RCI_MATVEC);
        for(size_t i=0;i<3;i++)q.out[i]=0;
        CHECK(lis_rci_advance(c,0,&q)==LIS_RCI_PRECON);
        CHECK(lis_rci_advance(c,23,&q)==LIS_RCI_CALLBACK);
        lis_rci_destroy(c);
    }
    puts("RCI core, reuse and protocol errors PASS");
    return 0;
}
