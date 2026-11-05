/* White-box BE hard-lock side accounting regression.
 * Include the implementation only in this standalone test translation
 * unit so the static Ritz selection helper can be exercised without
 * depending on which end happens to converge first in a synthetic solve.
 * Do NOT also link a separate lis_krylov.c object.
 */
#include "../src/esolver/lis_krylov.c"
#include <stdio.h>

static int check(int nev,int high_locked,int low_locked,
                 const int *expected,const int *side,int ncheck)
{
    const double val[6]={-6.0,-5.0,2.0,10.0,11.0,12.0};
    int order[6],work[6],i;
    sort_ritz(6,val,order,"BE",high_locked,low_locked,nev,work);
    for(i=0;i<ncheck;i++){
        if(order[i]!=expected[i] || work[i]!=side[i]){
            fprintf(stderr,"BE failure nev=%d high_locked=%d low_locked=%d index=%d actual=(%d,%d) expected=(%d,%d)\n",
                    nev,high_locked,low_locked,i,order[i],work[i],
                    expected[i],side[i]);
            return 1;
        }
    }
    return 0;
}
int main(void)
{
    static const int e0[]={5,0,4,1},s0[]={1,0,1,0};
    static const int e1[]={5},s1[]={1};
    static const int e2[]={0},s2[]={0};
    static const int e3[]={5,4},s3[]={1,1};
    static const int e4[]={0,1},s4[]={0,0};
    static const int e5[]={0,5},s5[]={0,1};
    int i,order[6],work[6];
    const double val[6]={-6.0,-5.0,2.0,10.0,11.0,12.0};
    if(check(4,0,0,e0,s0,4))return 1;
    if(check(2,0,1,e1,s1,1))return 2;
    if(check(2,1,0,e2,s2,1))return 3;
    if(check(4,0,2,e3,s3,2))return 4;
    if(check(4,2,0,e4,s4,2))return 5;
    if(check(4,1,1,e0,s0,2))return 6;
    if(check(3,0,1,e3,s3,2))return 7;
    if(check(3,1,0,e5,s5,2))return 8;
    sort_ritz(6,val,order,"LA",0,0,4,work);
    for(i=0;i<6;i++)if(order[i]!=5-i)return 9;
    puts("PASS BE side quotas for out-of-order hard locks and unchanged LA");
    return 0;
}
