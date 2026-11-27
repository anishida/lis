#!/usr/bin/env python3
"""Qualification only: patch a disposable native source checkout.

Never commit the patched file. All original algorithms and ABI symbols stay
available in this build; the public classical Krylov entry point uses the
same generated continuation as the external RCI. The original core remains
available only as a direct differential reference.
"""
from pathlib import Path
import hashlib
import subprocess

root = Path(__file__).resolve().parents[2]
p = root / "src/esolver/lis_krylov.c"
src = p.read_bytes()
blob = hashlib.sha1(b"blob " + str(len(src)).encode() + b"\0" + src).hexdigest()
assert blob == "f275ffbb5cf3ab0a6757d5e59eaefd9dd9f9b74c", blob
module = root / "qualification/rci-restore"
gen = module / "krylov_continuation.inc"
import tempfile
with tempfile.TemporaryDirectory() as tmp:
    check = Path(tmp) / "generated.inc"
    subprocess.run(["node", str(module / "generate_krylov_continuation.js"), str(p), str(check)], check=True)
    assert check.read_bytes() == gen.read_bytes()
s = src.decode("utf-8")
def replace_once(a, b):
    global s
    assert s.count(a) == 1, a
    s = s.replace(a, b)
replace_once('#include "lis_krylov.h"', '''#include "lis_krylov.h"
#define LIS_RCI_KRYLOV_NATIVE_LIBRARY 1
#include "../../qualification/rci-restore/lis_rci_krylov_restore.h"''')
replace_once('''    return lis_krylov_symmetric_core(workspace,p,o,ext,values,vectors,
                                      residuals,stats,ext_stats);''',
'''    return lis_rci_krylov_restore_classic64(p,o,ext,workspace,NULL,values,
                                           vectors,residuals,stats,ext_stats,NULL);''')
s += '''
/* Qualification-only integration, not a proposed permanent include path. */
#include "../../qualification/rci-restore/lis_rci_krylov_restore.c"
int lis_rci_krylov_restore_reference64(const lis_krylov_problem64 *p,
 const lis_krylov_options *o,const lis_krylov_ext_options *ext,
 lis_krylov_workspace *workspace,double *values,double *vectors,
 double *residuals,lis_krylov_stats *stats,lis_krylov_ext_stats *ext_stats)
{
 return lis_krylov_symmetric_core(workspace,p,o,ext,values,vectors,residuals,stats,ext_stats);
}
'''
p.write_text(s,encoding="utf-8")
print("KRYLOV NATIVE QUALIFICATION PATCH READY")
