/* Generic explicit, opt-in Matrix Market snapshot exporter for LIS.
 * No coupling to any FE application or element/DOF numbering.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lislib.h"
#include "lis_matrix_materialize.h"
#include "lis_snapshot.h"

#ifdef _COMPLEX
LIS_INT
lis_snapshot_export_mm(
    LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x0, LIS_VECTOR x,
    LIS_INT basis_count, LIS_VECTOR basis[], const char *prefix)
{
    (void)A; (void)b; (void)x0; (void)x;
    (void)basis_count; (void)basis; (void)prefix;
    return LIS_ERR_NOT_IMPLEMENTED;
}
#else

static LIS_INT
lis_snapshot_check_vector(LIS_VECTOR v, LIS_INT n)
{
    LIS_INT i;
    if(v==NULL) return LIS_SUCCESS;
    if(v->n!=n || v->gn!=n || v->is!=0 ||
       v->nprocs!=1 || v->value==NULL)
        return LIS_ERR_ILL_ARG;
    for(i=0;i<n;i++)
        if(!isfinite((double)v->value[i]))
            return LIS_ERR_ILL_ARG;
    return LIS_SUCCESS;
}

static LIS_INT
lis_snapshot_check_csr(LIS_MATRIX A)
{
    LIS_INT i,k,nnz;
    if(A==NULL || A->ptr==NULL || A->n<=0 ||
       A->nnz<0 || A->n!=A->gn ||
       A->nprocs!=1 || A->is!=0 || A->ptr[0]!=0)
        return LIS_ERR_ILL_ARG;
    nnz=A->ptr[A->n];
    if(nnz<0 || nnz!=A->nnz ||
       (nnz>0 && (A->index==NULL || A->value==NULL)))
        return LIS_ERR_ILL_ARG;
    for(i=0;i<A->n;i++)
        if(A->ptr[i]<0 || A->ptr[i]>A->ptr[i+1] ||
           A->ptr[i+1]>nnz)
            return LIS_ERR_ILL_ARG;
    for(k=0;k<nnz;k++)
        if(A->index[k]<0 || A->index[k]>=A->n ||
           !isfinite((double)A->value[k]))
            return LIS_ERR_ILL_ARG;
    return LIS_SUCCESS;
}

static LIS_INT
lis_snapshot_write_matrix(FILE *out, LIS_MATRIX A)
{
    LIS_INT i,k;
    if(fprintf(out,"%%%%MatrixMarket matrix coordinate real general\n"
                   "%lld %lld %lld\n",
               (long long)A->n,(long long)A->n,
               (long long)A->ptr[A->n])<0)
        return LIS_ERR_FILE_IO;
    for(i=0;i<A->n;i++)
        for(k=A->ptr[i];k<A->ptr[i+1];k++)
            if(fprintf(out,"%lld %lld %.17g\n",
                       (long long)i+1,
                       (long long)A->index[k]+1,
                       (double)A->value[k])<0)
                return LIS_ERR_FILE_IO;
    return LIS_SUCCESS;
}

static LIS_INT
lis_snapshot_write_vector(FILE *out, LIS_VECTOR v)
{
    LIS_INT i;
    if(fprintf(out,"%%%%MatrixMarket matrix array real general\n"
                   "%lld 1\n",(long long)v->n)<0)
        return LIS_ERR_FILE_IO;
    for(i=0;i<v->n;i++)
        if(fprintf(out,"%.17g\n",(double)v->value[i])<0)
            return LIS_ERR_FILE_IO;
    return LIS_SUCCESS;
}

static LIS_INT
lis_snapshot_write_basis(FILE *out,LIS_INT n,LIS_INT count,LIS_VECTOR basis[])
{
    LIS_INT m,i;
    if(fprintf(out,"%%%%MatrixMarket matrix array real general\n"
                   "%lld %lld\n",(long long)n,(long long)count)<0)
        return LIS_ERR_FILE_IO;
    for(m=0;m<count;m++)
        for(i=0;i<n;i++)
            if(fprintf(out,"%.17g\n",(double)basis[m]->value[i])<0)
                return LIS_ERR_FILE_IO;
    return LIS_SUCCESS;
}

/* All file writes go through a single close-and-check path. */
static LIS_INT
lis_snapshot_emit(char *path, size_t cap, const char *prefix,
                  const char *label, LIS_MATRIX A, LIS_VECTOR v,
                  LIS_INT count, LIS_VECTOR basis[])
{
    FILE *file;
    LIS_INT err;
    if(snprintf(path,cap,"%s_%s.mtx",prefix,label)<0)
        return LIS_ERR_ILL_ARG;
    file=fopen(path,"w");
    if(file==NULL) return LIS_ERR_FILE_IO;
    if(A!=NULL)
        err=lis_snapshot_write_matrix(file,A);
    else if(v!=NULL)
        err=lis_snapshot_write_vector(file,v);
    else
        err=lis_snapshot_write_basis(file,basis[0]->n,count,basis);
    if(fclose(file)!=0 && err==LIS_SUCCESS)
        err=LIS_ERR_FILE_IO;
    if(err!=LIS_SUCCESS)
        remove(path);
    return err;
}

LIS_INT
lis_snapshot_export_mm(
    LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x0, LIS_VECTOR x,
    LIS_INT basis_count, LIS_VECTOR basis[], const char *prefix)
{
    LIS_MATRIX csr=NULL,owned=NULL;
    const char *labels[4]={"A","b","x0","x"};
    LIS_VECTOR vectors[4];
    LIS_INT m,i,err;
    size_t n;
    char *path=NULL;

    if(A==NULL || prefix==NULL || prefix[0]=='\0' ||
       A->n<=0 || A->gn!=A->n || A->nprocs!=1 || A->is!=0 ||
       basis_count<0 || (basis_count>0 && basis==NULL))
        return LIS_ERR_ILL_ARG;

    /* Explicit opt-in: no output, allocation or side effect in lis_solve(). */
    err=lis_snapshot_check_vector(b,A->n);
    if(err) return err;
    err=lis_snapshot_check_vector(x0,A->n);
    if(err) return err;
    err=lis_snapshot_check_vector(x,A->n);
    if(err) return err;
    for(m=0;m<basis_count;m++)
    {
        if(basis[m]==NULL) return LIS_ERR_ILL_ARG;
        err=lis_snapshot_check_vector(basis[m],A->n);
        if(err) return err;
    }

    if(A->matrix_type==LIS_MATRIX_CSR)
    {
        /* Split CSR may contain entries in separately managed arrays. */
        if(A->is_splited) return LIS_ERR_NOT_IMPLEMENTED;
        csr=A;
    }
    else if(A->matrix_type==LIS_MATRIX_USER ||
            A->matrix_type==LIS_MATRIX_OPERATOR)
    {
        err=lis_matrix_materialize_csr_owned(A,&owned);
        if(err) return err;
        csr=owned;
    }
    else
    {
        err=lis_matrix_duplicate(A,&owned);
        if(err) return err;
        err=lis_matrix_set_type(owned,LIS_MATRIX_CSR);
        if(!err) err=lis_matrix_convert(A,owned);
        if(err) goto cleanup;
        csr=owned;
    }
    err=lis_snapshot_check_csr(csr);
    if(err) goto cleanup;

    n=strlen(prefix);
    if(n>((size_t)-1)-32)
    {
        err=LIS_ERR_ILL_ARG;
        goto cleanup;
    }
    path=(char *)lis_malloc(n+32,"lis_snapshot_export_mm::path");
    if(path==NULL)
    {
        err=LIS_OUT_OF_MEMORY;
        goto cleanup;
    }
    vectors[0]=NULL;
    vectors[1]=b;
    vectors[2]=x0;
    vectors[3]=x;
    for(i=0;i<4;i++)
    {
        if(i>0 && vectors[i]==NULL)
            continue;
        err=lis_snapshot_emit(path,n+32,prefix,labels[i],
                              i==0?csr:NULL,vectors[i],0,NULL);
        if(err) goto cleanup;
    }
    if(basis_count>0)
        err=lis_snapshot_emit(path,n+32,prefix,"modes",NULL,NULL,
                              basis_count,basis);

cleanup:
    if(path) lis_free(path);
    if(owned) lis_matrix_destroy(owned);
    return err;
}
#endif

/*
 * Write application-owned physical DOF labels, not inferred FE metadata.
 * This routine is intentionally available even in complex LIS builds:
 * the output is independent of numerical scalar type.
 */
static LIS_INT
lis_snapshot_csv_field(FILE *out, const char *value)
{
    const unsigned char *p=(const unsigned char *)(value ? value : "");
    if(fputc('"',out)==EOF) return LIS_ERR_FILE_IO;
    for(;*p;p++)
    {
        if(*p=='"' && fputc('"',out)==EOF)
            return LIS_ERR_FILE_IO;
        if(fputc((int)*p,out)==EOF)
            return LIS_ERR_FILE_IO;
    }
    if(fputc('"',out)==EOF) return LIS_ERR_FILE_IO;
    return LIS_SUCCESS;
}

LIS_INT
lis_snapshot_export_dof_map_csv(
    LIS_MATRIX A,
    LIS_SNAPSHOT_DOF_MAP_CALLBACK callback,
    void *user_data,
    const char *prefix)
{
    size_t length;
    char *path=NULL;
    FILE *out=NULL;
    LIS_INT i,err=LIS_SUCCESS;
    int needed;

    if(A==NULL || callback==NULL || prefix==NULL || prefix[0]=='\0' ||
       A->n<=0 || A->n!=A->gn || A->is!=0 || A->nprocs!=1)
        return LIS_ERR_ILL_ARG;

    length=strlen(prefix);
    if(length>((size_t)-1)-32)
        return LIS_ERR_ILL_ARG;
    path=(char *)lis_malloc(length+32,"lis_snapshot_export_dof_map_csv::path");
    if(path==NULL) return LIS_OUT_OF_MEMORY;

    needed=snprintf(path,length+32,"%s_dofs.csv",prefix);
    if(needed<0 || (size_t)needed>=length+32)
    {
        err=LIS_ERR_ILL_ARG;
        goto done;
    }
    out=fopen(path,"w");
    if(out==NULL)
    {
        err=LIS_ERR_FILE_IO;
        goto done;
    }
    if(fputs("eq,node,dof,body\n",out)==EOF)
    {
        err=LIS_ERR_FILE_IO;
        goto close_file;
    }

    for(i=0;i<A->n;i++)
    {
        const char *node=NULL,*dof=NULL,*body=NULL;
        err=callback(user_data,i,&node,&dof,&body);
        if(err) goto close_file;
        if(node==NULL || node[0]=='\0' || dof==NULL || dof[0]=='\0')
        {
            err=LIS_ERR_ILL_ARG;
            goto close_file;
        }
        if(fprintf(out,"%lld,",(long long)i)<0)
        {
            err=LIS_ERR_FILE_IO;
            goto close_file;
        }
        err=lis_snapshot_csv_field(out,node);
        if(err) goto close_file;
        if(fputc(',',out)==EOF) { err=LIS_ERR_FILE_IO; goto close_file; }
        err=lis_snapshot_csv_field(out,dof);
        if(err) goto close_file;
        if(fputc(',',out)==EOF) { err=LIS_ERR_FILE_IO; goto close_file; }
        err=lis_snapshot_csv_field(out,body);
        if(err) goto close_file;
        if(fputc('\n',out)==EOF) { err=LIS_ERR_FILE_IO; goto close_file; }
    }

close_file:
    if(fclose(out)!=0 && err==LIS_SUCCESS)
        err=LIS_ERR_FILE_IO;
    out=NULL;
    if(err!=LIS_SUCCESS)
        remove(path);
done:
    if(path) lis_free(path);
    return err;
}
