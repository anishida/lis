/* Test-only Matrix Market driver (coordinate real general, 1-based indexes).
 * Usage: ./check_real_mtx matrix_A.mtx blocks.txt
 * blocks.txt: first line N_BLOCKS, then one row per block: m i0 ... i(m-1) (0-based).
 * Only tests Schur identity/reconstruction; does NOT solve with LIS or SuperLU.
 */
#include "lis_condense_core.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef struct {size_t r,c;double v;} triple;
static int cmp(const void*ap,const void*bp){const triple*a=ap,*b=bp;
 if(a->r!=b->r)return (a->r>b->r)-(a->r<b->r);
 return (a->c>b->c)-(a->c<b->c);}
static void err(const char*m){fprintf(stderr,"ERROR %s\n",m);exit(2);}
static void mv(const lis_sc_csr*a,const double*x,double*b){
 for(size_t i=0;i<a->n;i++){long double z=0;
  for(size_t k=a->rowptr[i];k<a->rowptr[i+1];k++)z+=(long double)a->values[k]*x[a->colind[k]];
  b[i]=(double)z;}
}
static double norm_rel(const double*a,const double*b,size_t n){double num=0,den=0;
 for(size_t i=0;i<n;i++){num=hypot(num,a[i]-b[i]);den=hypot(den,b[i]);}
 return num/(den>0?den:1);}
int main(int argc,char**argv){
 if(argc!=3){fprintf(stderr,"usage: %s A.mtx blocks.txt\n",argv[0]);return 2;}
 FILE*f=fopen(argv[1],"r");if(!f)err("cannot read mtx");
 char line[4096];if(!fgets(line,sizeof(line),f)||strstr(line,"%%MatrixMarket matrix coordinate real general")!=line)
   err("expected Matrix Market coordinate real general");
 do{if(!fgets(line,sizeof(line),f))err("missing dims");}while(line[0]=='%');
 size_t n,ncol,nnz;
 if(sscanf(line,"%zu%zu%zu",&n,&ncol,&nnz)!=3||n!=ncol)err("bad dimensions");
 triple *t=calloc(nnz?nnz:1,sizeof(triple));if(!t)err("out of memory triples");
 for(size_t k=0;k<nnz;k++){
   size_t ri,ci;double v;
   if(fscanf(f,"%zu %zu %lf",&ri,&ci,&v)!=3||ri==0||ci==0||ri>n||ci>n)err("bad matrix triplet");
   t[k]=(triple){ri-1,ci-1,v};
 }
 fclose(f);qsort(t,nnz,sizeof(triple),cmp);
 size_t unique=0;for(size_t k=0;k<nnz;k++)if(!k||cmp(&t[k],&t[k-1])!=0)unique++;
 size_t *ptr=calloc(n+1,sizeof(size_t)),*idx=calloc(unique+1,sizeof(size_t));
 double *val=calloc(unique+1,sizeof(double));if(!ptr||!idx||!val)err("out of memory csr");
 size_t at=0;for(size_t k=0;k<nnz;k++){
   if(k && cmp(&t[k],&t[k-1])==0){val[at-1]+=t[k].v;continue;}
   ptr[t[k].r+1]++;idx[at]=t[k].c;val[at++]=t[k].v;
 }
 for(size_t i=0;i<n;i++)ptr[i+1]+=ptr[i];
 free(t);
 lis_sc_csr A={n,unique,ptr,idx,val};
 FILE*bf=fopen(argv[2],"r");if(!bf)err("cannot read blocks");
 size_t nb=0;
 if(fscanf(bf,"%zu",&nb)!=1||!nb)err("bad blocks header");
 size_t *bp=calloc(nb+1,sizeof(size_t)),*qd=calloc(n,sizeof(size_t));
 if(!bp||!qd)err("out of memory blocks");
 size_t nq=0;
 for(size_t b=0;b<nb;b++){
   size_t m;if(fscanf(bf,"%zu",&m)!=1||!m||m>n-nq)err("bad block count");
   bp[b]=nq;
   for(size_t k=0;k<m;k++)if(fscanf(bf,"%zu",&qd[nq++])!=1)err("bad block index");
 }
 bp[nb]=nq;fclose(bf);
 lis_sc_blocks B={nb,bp,qd};lis_sc_plan *P=NULL;
 clock_t t0=clock();lis_sc_status st=lis_sc_build(&P,&A,&B,LIS_SC_FULL,1e-13);
 double t_setup=(double)(clock()-t0)/CLOCKS_PER_SEC;
 if(st!=LIS_SC_OK){fprintf(stderr,"CONDENSATION FAILED %s\n",lis_sc_status_string(st));return 3;}
 size_t nu=lis_sc_physical_dofs(P);
 double *x=calloc(n,sizeof(double)),*b=calloc(n,sizeof(double)),*xr=calloc(n,sizeof(double));
 double *u=calloc(nu,sizeof(double)),*g=calloc(nu,sizeof(double)),*Su=calloc(nu,sizeof(double));
 if(!x||!b||!xr||!u||!g||!Su)err("out of memory vectors");
 for(size_t i=0;i<n;i++)x[i]=sin((double)(i+1)*0.025)+cos((double)(i+1)*0.006);
 mv(&A,x,b);
 const size_t*u2=lis_sc_physical_to_full(P);
 for(size_t i=0;i<nu;i++)u[i]=x[u2[i]];
 st=lis_sc_rhs(P,b,g);if(st!=LIS_SC_OK)err(lis_sc_status_string(st));
 mv(lis_sc_matrix(P),u,Su);
 double eS=norm_rel(Su,g,nu);
 st=lis_sc_reconstruct(P,u,b,xr);if(st!=LIS_SC_OK)err(lis_sc_status_string(st));
 double ex=norm_rel(xr,x,n),rr=0;
 st=lis_sc_true_residual(P,b,xr,&rr);if(st!=LIS_SC_OK)err(lis_sc_status_string(st));
 printf("n=%zu nphys=%zu nq=%zu blocks=%zu nnzA=%zu nnzS=%zu setup_cpu_s=%.6f\n",n,nu,nq,nb,unique,lis_sc_matrix(P)->nnz,t_setup);
 printf("Schur identity relerr=%.12g, reconstruction relerr=%.12g, true residual=%.12g\n",eS,ex,rr);
 /* A new Newton tangent with the same CSR sparsity and block metadata:
  * change one retained diagonal, then compare true numerical reuse against
  * a completely rebuilt Schur.  This is not a benchmark of native LIS. */
 size_t physical0=lis_sc_physical_to_full(P)[0], diagpos=SIZE_MAX;
 for(size_t z=ptr[physical0];z<ptr[physical0+1];z++)
     if(idx[z]==physical0){diagpos=z;break;}
 if(diagpos==SIZE_MAX)err("missing retained diagonal");
 val[diagpos]+=0.000123*(fabs(val[diagpos])+1.0);
 clock_t t1=clock();
 st=lis_sc_build(&P,&A,&B,LIS_SC_NUMERIC,1e-13);
 double t_numeric=(double)(clock()-t1)/CLOCKS_PER_SEC;
 if(st!=LIS_SC_OK){fprintf(stderr,"NUMERIC REUSE FAILED %s\n",lis_sc_status_string(st));return 5;}
 lis_sc_plan *Fresh=NULL;
 t1=clock();st=lis_sc_build(&Fresh,&A,&B,LIS_SC_FULL,1e-13);
 double t_fresh=(double)(clock()-t1)/CLOCKS_PER_SEC;
 if(st!=LIS_SC_OK){fprintf(stderr,"FRESH REBUILD FAILED %s\n",lis_sc_status_string(st));return 6;}
 const lis_sc_csr *SN=lis_sc_matrix(P),*SF=lis_sc_matrix(Fresh);
 if(SN->n!=SF->n||SN->nnz!=SF->nnz||
    memcmp(SN->rowptr,SF->rowptr,(SN->n+1)*sizeof(size_t))!=0||
    memcmp(SN->colind,SF->colind,SN->nnz*sizeof(size_t))!=0)
     err("Schur pattern mismatch between FULL and NUMERIC");
 double schur_num_error=norm_rel(SN->values,SF->values,SN->nnz);
 mv(&A,x,b);
 st=lis_sc_rhs(P,b,g);if(st!=LIS_SC_OK)err(lis_sc_status_string(st));
 mv(SN,u,Su);
 double updated_error=norm_rel(Su,g,nu);
 st=lis_sc_reconstruct(P,u,b,xr);if(st!=LIS_SC_OK)err(lis_sc_status_string(st));
 double new_rr=0;
 st=lis_sc_true_residual(P,b,xr,&new_rr);if(st!=LIS_SC_OK)err(lis_sc_status_string(st));
 printf("numeric_cpu_s=%.6f fresh_cpu_s=%.6f numerical_schur_err=%.12g updated_identity_err=%.12g updated_residual=%.12g\n",
        t_numeric,t_fresh,schur_num_error,updated_error,new_rr);
 free(x);free(b);free(xr);free(u);free(g);free(Su);
 lis_sc_free(P);lis_sc_free(Fresh);free(ptr);free(idx);free(val);free(bp);free(qd);
 return (eS<1e-9 && ex<1e-9 && rr<1e-9 && schur_num_error<1e-10 &&
         updated_error<1e-9 && new_rr<1e-9)?0:4;
}
