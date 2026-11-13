/* Generic graph queries derived from matrix entries. */

#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#include <stdlib.h>
#include "lislib.h"
#include "lis_precon.h"
#include "lis_matrix_graph.h"

typedef struct
{
    LIS_INT col;
    LIS_SCALAR value;
} LIS_MATRIX_GRAPH_ENTRY;

static int
lis_matrix_graph_entry_compare(
    const void *a,
    const void *b)
{
    const LIS_MATRIX_GRAPH_ENTRY *ea;
    const LIS_MATRIX_GRAPH_ENTRY *eb;

    ea = (const LIS_MATRIX_GRAPH_ENTRY *)a;
    eb = (const LIS_MATRIX_GRAPH_ENTRY *)b;
    if(ea->col<eb->col) return -1;
    if(ea->col>eb->col) return 1;
    return 0;
}

static LIS_REAL
lis_matrix_graph_scalar_abs(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return (LIS_REAL)cabs(value);
#else
    return (LIS_REAL)fabs(value);
#endif
}

static LIS_INT
lis_matrix_graph_find(
    LIS_INT *parent,
    LIS_INT x)
{
    LIS_INT root;
    LIS_INT next;

    root = x;
    while(parent[root]!=root) root = parent[root];

    while(parent[x]!=x)
    {
        next = parent[x];
        parent[x] = root;
        x = next;
    }

    return root;
}

static void
lis_matrix_graph_union(
    LIS_INT *parent,
    LIS_INT *rank,
    LIS_INT a,
    LIS_INT b)
{
    LIS_INT ra;
    LIS_INT rb;

    ra = lis_matrix_graph_find(parent,a);
    rb = lis_matrix_graph_find(parent,b);
    if(ra==rb) return;

    if(rank[ra]<rank[rb]) parent[ra] = rb;
    else if(rank[ra]>rank[rb]) parent[rb] = ra;
    else
    {
        parent[rb] = ra;
        rank[ra]++;
    }
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_connected_components"
LIS_INT
lis_matrix_connected_components(
    LIS_MATRIX A,
    LIS_REAL drop_tolerance,
    LIS_INT *labels,
    LIS_INT *component_count)
{
    LIS_MATRIX C;
    LIS_MATRIX active;
    LIS_PRECON_ROW_SOURCE row;
    LIS_MATRIX_GRAPH_ENTRY *entries;
    LIS_INT entry_capacity;
    LIS_INT *parent;
    LIS_INT *rank;
    LIS_INT *root_label;
    LIS_INT i,k,j,root,count,err;
    LIS_INT row_source_ready;

    LIS_DEBUG_FUNC_IN;

    if(A==NULL || labels==NULL || component_count==NULL ||
       !(drop_tolerance>=(LIS_REAL)0.0))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,"invalid matrix connected-components arguments\n");
        return LIS_ERR_ILL_ARG;
    }

    err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
    if(err) return err;

    if(A->nprocs!=1 || A->n!=A->gn)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "distributed matrix connected components are not implemented\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }

    if(A->n==0)
    {
        *component_count = 0;
        return LIS_SUCCESS;
    }

    if(A->matrix_type==LIS_MATRIX_OPERATOR)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "matrix connected components require explicit or row-backed matrix entries\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }

    C = NULL;
    active = A;
    if(A->matrix_type!=LIS_MATRIX_CSR && A->matrix_type!=LIS_MATRIX_USER)
    {
        err = lis_matrix_duplicate(A,&C);
        if(err) return err;
        err = lis_matrix_set_type(C,LIS_MATRIX_CSR);
        if(!err) err = lis_matrix_convert(A,C);
        if(err)
        {
            lis_matrix_destroy(C);
            return err;
        }
        active = C;
    }

    parent = NULL;
    rank = NULL;
    root_label = NULL;
    entries = NULL;
    entry_capacity = 0;
    row_source_ready = LIS_FALSE;

    parent = (LIS_INT *)lis_malloc(active->n*sizeof(LIS_INT),
        "lis_matrix_connected_components::parent");
    rank = (LIS_INT *)lis_malloc(active->n*sizeof(LIS_INT),
        "lis_matrix_connected_components::rank");
    root_label = (LIS_INT *)lis_malloc(active->n*sizeof(LIS_INT),
        "lis_matrix_connected_components::root_label");
    if(parent==NULL || rank==NULL || root_label==NULL)
    {
        err = LIS_ERR_OUT_OF_MEMORY;
        LIS_SETERR_MEM(3*active->n*sizeof(LIS_INT));
        goto cleanup;
    }

    for(i=0;i<active->n;i++)
    {
        parent[i] = i;
        rank[i] = 0;
        root_label[i] = -1;
    }

    err = lis_precon_row_source_init(active,&row);
    if(err) goto cleanup;
    row_source_ready = LIS_TRUE;

    for(i=0;i<active->n;i++)
    {
        err = lis_precon_row_source_get(&row,i);
        if(err) goto cleanup;

        if(row.nnz>entry_capacity)
        {
            LIS_MATRIX_GRAPH_ENTRY *tmp;
            tmp = (LIS_MATRIX_GRAPH_ENTRY *)lis_realloc(
                entries,row.nnz*sizeof(LIS_MATRIX_GRAPH_ENTRY));
            if(tmp==NULL)
            {
                err = LIS_ERR_OUT_OF_MEMORY;
                LIS_SETERR_MEM(row.nnz*sizeof(LIS_MATRIX_GRAPH_ENTRY));
                goto cleanup;
            }
            entries = tmp;
            entry_capacity = row.nnz;
        }

        for(k=0;k<row.nnz;k++)
        {
            if(row.index[k]<0 || row.index[k]>=active->n)
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "matrix graph row contains an out-of-range column\n");
                err = LIS_ERR_ILL_ARG;
                goto cleanup;
            }
            entries[k].col = row.index[k];
            entries[k].value = row.value[k];
        }

        if(row.nnz>1)
            qsort(entries,(size_t)row.nnz,sizeof(*entries),
                  lis_matrix_graph_entry_compare);

        k = 0;
        while(k<row.nnz)
        {
            LIS_INT col;
            LIS_SCALAR combined;

            col = entries[k].col;
            combined = entries[k].value;
            j = k+1;
            while(j<row.nnz && entries[j].col==col)
            {
                combined += entries[j].value;
                j++;
            }

            if(col!=i && lis_matrix_graph_scalar_abs(combined)>drop_tolerance)
                lis_matrix_graph_union(parent,rank,i,col);

            k = j;
        }
    }

    count = 0;
    for(i=0;i<active->n;i++)
    {
        root = lis_matrix_graph_find(parent,i);
        if(root_label[root]<0) root_label[root] = count++;
        labels[i] = root_label[root];
    }

    *component_count = count;
    err = LIS_SUCCESS;

cleanup:
    if(row_source_ready) lis_precon_row_source_destroy(&row);
    lis_free(entries);
    lis_free(root_label);
    lis_free(rank);
    lis_free(parent);
    if(C) lis_matrix_destroy(C);

    LIS_DEBUG_FUNC_OUT;
    return err;
}
