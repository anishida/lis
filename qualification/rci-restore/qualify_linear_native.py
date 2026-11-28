#!/usr/bin/env python3
"""Patch only the fresh qualification checkout, not upstream sources."""
from pathlib import Path
from generate_linear_continuation import ROOT, DIR, SPECS, generate
output, counts = generate()
(DIR/"linear_continuation.inc").write_text(output)
print("NATIVE LINEAR GENERATED", counts)
for kind,func,path,sha in SPECS[1:]:
    p=ROOT/path
    s=p.read_text()
    a="LIS_INT "+func+"(LIS_SOLVER solver)"
    assert s.count(a)==1,a
    p.write_text(s.replace(a,"LIS_INT "+func+"_native_reference(LIS_SOLVER solver)"))
p=ROOT/"src/solver/lis_solver_cg.c"
p.write_text(p.read_text()+'\n#include "../../qualification/rci-restore/lis_rci_native_linear_restore.c"\n')
print("NATIVE LINEAR QUALIFICATION PATCH READY")
