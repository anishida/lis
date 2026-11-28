#!/usr/bin/env python3
"""Qualification-only native positive-k FGCRO-DR continuations.
Source-specific, real double, non-MPI. No numerical recurrence is reimplemented.
The zero-recycle classical FGMRES path stays in its original translation unit.
"""
import re
import generate_linear_continuation as g

p = "src/solver/lis_solver_fgcrodr.c"
sha = "a0cccf9503b7f0571c3eafa2bdf96e1df5719a10" if "--complex" in __import__("sys").argv else "fea73bea4d85b65d1204c0715ae7e8660a76bbdc"
specs = [
    ("install", "install_recycle", ["solver", "v", "z", "jdim"], 1),
    ("harmonic_build", "try_harmonic_build", ["solver", "v", "z", "jdim", "hraw", "hld"], 1),
    ("replace", "replace_recycle_basis", ["solver", "neww", "newz", "nsel", "count_update", "action", "quality", "quality_min"], 1),
    ("commit", "commit_recycle_update", ["solver", "neww", "newz", "nsel"], 1),
    ("harmonic_update", "try_harmonic_update", ["solver", "v", "z", "jdim", "hraw", "hld", "bcoef", "bld"], 1),
    ("build", "build_recycle", ["solver", "v", "z", "jdim", "hraw", "hld", "bcoef", "bld"], 3),
    ("prepare", "prepare_recycle", ["solver"], 2),
    ("core", "flexible_core", ["solver"], 5),
    ("fgcro", "", ["solver"], 2),
]
def native(suffix):
    return "lis_fgcrodr" + ("_" + suffix if suffix else "")

g.SPECS = [g.SPECS[0]] + [(name, native(suffix), p, sha) for name, suffix, _, _ in specs]
g.IDS = {name:i+1 for i,(name,*_) in enumerate(g.SPECS)}
import sys
g.COMPLEX_ONLY = '--complex' in sys.argv
g.REAL_ONLY = not g.COMPLEX_ONLY
g.EXTRA_CALLS = {native(suffix):(name, params) for name,suffix,params,_ in specs if suffix}
g.EXPECTED_COUNTS = {"initial":2, **{name:count for name,_,_,count in specs}}
original_extract = g.extract

def extract(src, name):
    params, body = original_extract(src, name)
    if name == "lis_fgcrodr_commit_recycle_update":
        body = "\n    LIS_INT err;\n" + body
    children = "|".join(g.EXTRA_CALLS)
    body, n = re.subn(r"\breturn\s+((?:"+children+r")\s*\([^;]*?\));",
                      r"err = \1;\n    return err;", body)
    expected = {"lis_fgcrodr_commit_recycle_update":1,"lis_fgcrodr_build_recycle":2}.get(name,0)
    assert n == expected, (name,n)
    if name == "lis_fgcrodr":
        body,n = re.subn(r"return\s+lis_fgmres\s*\(\s*solver\s*\);",
                         "return LIS_ERR_ILL_ARG;",body)
        assert n == 1
    return params, body

g.extract = extract
out, counts = g.generate()
(g.DIR/"fgcro_continuation.inc").write_text(out.replace("lis_rci_native_linear_restore","lis_rci_native_fgcro_restore"))
print("NATIVE FGCRO GENERATED", counts)
src = g.ROOT/p
s = src.read_text()
a = "LIS_INT\nlis_fgcrodr(\n    LIS_SOLVER solver)"
assert s.count(a) == 1
s = s.replace(a,"LIS_INT\nlis_fgcrodr_native_reference(\n    LIS_SOLVER solver)")
a = "#undef FGCRODR_DEFAULT_RECYCLE_DIM"
assert s.count(a) == 1
s = s.replace(a,'#include "../../qualification/rci-restore/lis_rci_native_fgcro_restore.c"\n\n'+a)
src.write_text(s)
print("NATIVE FGCRO QUALIFICATION PATCH READY")
