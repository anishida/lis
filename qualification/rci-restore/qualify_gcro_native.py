#!/usr/bin/env python3
"""Generate native real GCRO-DR, including all recycle callback sites."""
import generate_linear_continuation as g
p="src/solver/lis_solver_gcrodr.c"
sha="4801f40031c48c14ea359e67f5b4aeb5a8bd6c24"
g.SPECS=[g.SPECS[0],
    ("refresh","lis_gcrodr_refresh_recycle",p,sha),
    ("assess","lis_gcrodr_assess_refresh_recycle",p,sha),
    ("update","lis_gcrodr_update_recycle",p,sha),
    ("gcro","lis_gcrodr",p,sha)]
g.IDS={name:i+1 for i,(name,*_) in enumerate(g.SPECS)}
import sys
g.COMPLEX_ONLY = '--complex' in sys.argv
g.REAL_ONLY = not g.COMPLEX_ONLY
g.EXTRA_CALLS={
    "lis_gcrodr_refresh_recycle":("refresh",["solver"]),
    "lis_gcrodr_assess_refresh_recycle":("assess",["solver","scratch_z","scratch_c"]),
    "lis_gcrodr_update_recycle":("update",["solver","v","jdim","hraw","hld","bcoef","bld"])}
g.EXPECTED_COUNTS={"initial":2,"refresh":2,"assess":3,"update":1,"gcro":6}
out,counts=g.generate()
(g.DIR/"gcro_continuation.inc").write_text(out.replace("lis_rci_native_linear_restore","lis_rci_native_gcro_restore"))
print("NATIVE GCRO GENERATED",counts)
src=g.ROOT/p
s=src.read_text()
a="LIS_INT lis_gcrodr(LIS_SOLVER solver)"
assert s.count(a)==1
s=s.replace(a,"LIS_INT lis_gcrodr_native_reference(LIS_SOLVER solver)")
a="#undef GCRODR_NWORK\n#undef GCRODR_TINY\n#undef GCRODR_REORTH_TINY"
assert s.count(a)==1
s=s.replace(a,'#include "../../qualification/rci-restore/lis_rci_native_gcro_restore.c"\n\n'+a)
src.write_text(s)
print("NATIVE GCRO QUALIFICATION PATCH READY")
