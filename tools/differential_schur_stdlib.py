#!/usr/bin/env python3
"""Deterministic, stdlib-only differential test of the portable C Schur core.

Builds temporary C shared library with system C compiler, compares against an
independently written Python dense Gaussian elimination reference. No NumPy,
SciPy, SuperLU or external numerical package is used, in build or execution.
Does NOT exercise native LIS, MFAMG or FGMRES.
"""
from __future__ import annotations
import ctypes as C
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
KIND = C.c_size_t
DBL = C.c_double

class CSR(C.Structure):
    _fields_ = [('n', KIND), ('nnz', KIND), ('rowptr', C.POINTER(KIND)),
                ('colind', C.POINTER(KIND)), ('values', C.POINTER(DBL))]

class Blocks(C.Structure):
    _fields_ = [('nblocks', KIND), ('block_ptr', C.POINTER(KIND)), ('block_dof', C.POINTER(KIND))]


def ref_solve(mat, rhs):
    """Independent small dense partial-pivot reference."""
    n=len(rhs)
    a=[list(mat[i])+[float(rhs[i])] for i in range(n)]
    for k in range(n):
        pivot=max(range(k,n),key=lambda i:abs(a[i][k]))
        if abs(a[pivot][k])<1e-12:
            raise ValueError('near singular reference matrix')
        a[k],a[pivot]=a[pivot],a[k]
        for i in range(k+1,n):
            factor=a[i][k]/a[k][k]
            for j in range(k+1,n+1):
                a[i][j]-=factor*a[k][j]
            a[i][k]=0
    x=[0.0]*n
    for i in range(n-1,-1,-1):
        x[i]=(a[i][n]-sum(a[i][j]*x[j] for j in range(i+1,n)))/a[i][i]
    return x


def schur_reference(mat,blocks):
    n=len(mat)
    q=[d for block in blocks for d in block]
    qs=set(q)
    u=[i for i in range(n) if i not in qs]
    Kqq=[[mat[i][j] for j in q] for i in q]
    Kuu=[[mat[i][j] for j in u] for i in u]
    if not q:
        return u, Kuu
    Kqu=[[mat[i][j] for j in u] for i in q]
    Kuq=[[mat[i][j] for j in q] for i in u]
    reduced=[row[:] for row in Kuu]
    for j in range(len(u)):
        sol=ref_solve(Kqq,[row[j] for row in Kqu])
        for i in range(len(u)):
            reduced[i][j]-=sum(Kuq[i][k]*sol[k] for k in range(len(q)))
    return u, reduced


def setup(lib):
    pp=C.POINTER(C.c_void_p)
    lib.lis_sc_build.argtypes=[pp,C.POINTER(CSR),C.POINTER(Blocks),C.c_int,DBL]
    lib.lis_sc_build.restype=C.c_int
    lib.lis_sc_free.argtypes=[C.c_void_p]
    lib.lis_sc_matrix.argtypes=[C.c_void_p]
    lib.lis_sc_matrix.restype=C.POINTER(CSR)
    lib.lis_sc_rhs.argtypes=[C.c_void_p,C.POINTER(DBL),C.POINTER(DBL)]
    lib.lis_sc_rhs.restype=C.c_int
    lib.lis_sc_reconstruct.argtypes=[C.c_void_p,C.POINTER(DBL),C.POINTER(DBL),C.POINTER(DBL)]
    lib.lis_sc_reconstruct.restype=C.c_int
    lib.lis_sc_true_residual.argtypes=[C.c_void_p,C.POINTER(DBL),C.POINTER(DBL),C.POINTER(DBL)]
    lib.lis_sc_true_residual.restype=C.c_int


def c_csr(rows):
    n=len(rows)
    ptr=[0]
    col=[]
    val=[]
    for row in rows:
        col.extend(sorted(row))
        val.extend(row[j] for j in sorted(row))
        ptr.append(len(col))
    arptr=(KIND*len(ptr))(*ptr)
    arcol=(KIND*len(col))(*col)
    arval=(DBL*len(val))(*val)
    return CSR(n,len(col),arptr,arcol,arval),(arptr,arcol,arval)


def c_blocks(blocks):
    ptr=[0]
    q=[]
    for b in blocks:
        q.extend(b)
        ptr.append(len(q))
    bp=(KIND*len(ptr))(*ptr)
    dof=(KIND*len(q))(*q)
    return Blocks(len(blocks),bp,dof),(bp,dof)


def extract_mat(c):
    return [[float(c.values[k]) for k in range(int(c.rowptr[i]),int(c.rowptr[i+1]))]
            for i in range(int(c.n))]


def dense_from_rows(rows):
    n=len(rows)
    return [[rows[i].get(j,0.0) for j in range(n)] for i in range(n)]


def verify_plan(lib,plan,rows,blocks,rand):
    mat=dense_from_rows(rows)
    u,S=schur_reference(mat,blocks)
    sv=lib.lis_sc_matrix(plan).contents
    assert sv.n==len(u), (sv.n,len(u))
    for i in range(len(u)):
        got={int(sv.colind[k]):float(sv.values[k])
             for k in range(int(sv.rowptr[i]),int(sv.rowptr[i+1]))}
        for j in range(len(u)):
            err=abs(got.get(j,0.0)-S[i][j])
            assert err<2e-10,(i,j,err,got.get(j,0),S[i][j])
    xtrue=[rand.uniform(-2,2) for _ in range(len(mat))]
    b=[sum(aij*xtrue[j] for j,aij in enumerate(row)) for row in mat]
    B=(DBL*len(b))(*b)
    g=(DBL*len(u))()
    assert lib.lis_sc_rhs(plan,B,g)==0
    for i in range(len(u)):
        check=sum(S[i][j]*xtrue[u[j]] for j in range(len(u)))
        assert abs(g[i]-check)<2e-9,(i,g[i],check)
    U=(DBL*len(u))(*(xtrue[i] for i in u))
    out=(DBL*len(mat))()
    assert lib.lis_sc_reconstruct(plan,U,B,out)==0
    assert max(abs(out[i]-xtrue[i]) for i in range(len(mat)))<2e-9
    resid=DBL()
    assert lib.lis_sc_true_residual(plan,B,out,C.byref(resid))==0
    assert resid.value<2e-9
    return sv


def gen_case(rand):
    nphys=rand.randint(2,28)
    nb=rand.randint(0,4)
    blocksizes=[rand.randint(1,6) for _ in range(nb)]
    n=nphys+sum(blocksizes)
    dofs=list(range(n))
    rand.shuffle(dofs)
    q=dofs[:sum(blocksizes)]
    blocks=[]
    off=0
    for sz in blocksizes:
        blocks.append(q[off:off+sz]);off+=sz
    block_of={d:ib for ib,block in enumerate(blocks) for d in block}
    rows=[{} for _ in range(n)]
    for i in range(n):
        sign=-1.0 if i in block_of and rand.random()<0.25 else 1.0
        rows[i][i]=sign*(8.0+rand.random()*2)
        for j in range(n):
            if i==j:continue
            if i in block_of and j in block_of and block_of[i]!=block_of[j]:
                if rand.random()<0.04:rows[i][j]=0.0 # structural zero, must not break locality
                continue
            p=0.75 if (i in block_of and j in block_of) else (0.35 if (i in block_of or j in block_of) else 0.08)
            if rand.random()<p:
                rows[i][j]=0.0 if rand.random()<0.2 else rand.uniform(-0.35,0.35)
    return rows,blocks


def update_same_pattern(rand,rows,blocks):
    changed=[dict(r) for r in rows]
    block_of={d:ib for ib,block in enumerate(blocks) for d in block}
    for i,r in enumerate(changed):
        for j,v in list(r.items()):
            if i==j:
                r[j]=v+math.copysign(rand.uniform(.001,.15),v)
            elif (i not in block_of or j not in block_of or block_of[i]==block_of[j]) and rand.random()<0.3:
                r[j]=rand.uniform(-.3,.3) # includes activating structural zeros
    # Preserve all cross-block zero constraints for valid numerical updates.
    return changed


def run_cases(lib, count=500):
    rand=random.Random(0xCA21F00D)
    verified=0
    rejected_cross=0
    rejected_reuse=0
    zero_block=0
    for case in range(count):
        rows,blocks=gen_case(rand)
        if not blocks:zero_block+=1
        A,ka=c_csr(rows)
        BP,kb=c_blocks(blocks)
        slot=C.c_void_p(None)
        assert lib.lis_sc_build(C.byref(slot),C.byref(A),C.byref(BP),0,1e-13)==0,case
        try:
            verify_plan(lib,slot,rows,blocks,rand)
            updated=update_same_pattern(rand,rows,blocks)
            B,kc=c_csr(updated)
            assert lib.lis_sc_build(C.byref(slot),C.byref(B),C.byref(BP),1,1e-13)==0,case
            sv=verify_plan(lib,slot,updated,blocks,rand)
            values=(float(sv.values[k]) for k in range(int(sv.nnz)))
            old_values=list(values)
            old_ptr=slot.value
            # REUSE_ALL must reject changed coefficients even with the same CSR.
            bad_reuse=lib.lis_sc_build(C.byref(slot),C.byref(A),C.byref(BP),2,1e-13)
            assert bad_reuse==9 and slot.value==old_ptr,(case,bad_reuse)
            rejected_reuse+=1
            assert lib.lis_sc_build(C.byref(slot),C.byref(B),C.byref(BP),2,1e-13)==0
            assert slot.value==old_ptr # REUSE_ALL must be true no-op
            assert [float(sv.values[k]) for k in range(int(sv.nnz))]==old_values
            # FULL and NUMERIC must be equivalent even with activated zeros.
            full=C.c_void_p(None)
            assert lib.lis_sc_build(C.byref(full),C.byref(B),C.byref(BP),0,1e-13)==0
            try:
                s2=lib.lis_sc_matrix(full).contents
                assert int(sv.nnz)==int(s2.nnz)
                for k in range(int(sv.nnz)):
                    assert int(sv.colind[k])==int(s2.colind[k])
                    assert abs(float(sv.values[k])-float(s2.values[k]))<2e-10
            finally:
                lib.lis_sc_free(full)
            if len(blocks)>=2:
                # Activating explicitly stored cross-block zero must be rejected.
                i,j=blocks[0][0],blocks[1][0]
                if j in rows[i]:
                    invalid=[dict(r) for r in updated]
                    invalid[i][j]=0.4
                    BAD,kbad=c_csr(invalid)
                    preserved=slot.value
                    code=lib.lis_sc_build(C.byref(slot),C.byref(BAD),C.byref(BP),1,1e-13)
                    assert code==4 and slot.value==preserved,(case,code)
                    rejected_cross+=1
            verified+=1
        finally:
            lib.lis_sc_free(slot)
    assert rejected_cross>0 and zero_block>0
    print(f'STDLIB DIFFERENTIAL: {verified}/{count} cases PASS, zero-block={zero_block}, cross-block rejected={rejected_cross}, REUSE_ALL mismatch rejected={rejected_reuse}; independent dense reference, zero external numerical dependencies')


def main():
    with tempfile.TemporaryDirectory(prefix='ca_mfamg_diff_') as t:
        lib=Path(t)/'libcondense.so'
        subprocess.run([os.environ.get('CC','cc'),'-std=c99','-O2','-Wall','-Wextra','-Werror',
                        '-shared','-fPIC','-I'+str(ROOT/'include'),
                        str(ROOT/'src/matrix/lis_condense_core.c'),'-o',str(lib),'-lm'],check=True)
        dll=C.CDLL(str(lib))
        setup(dll)
        run_cases(dll,500)

if __name__=='__main__':
    main()
