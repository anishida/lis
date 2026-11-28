#!/usr/bin/env python3
"""Materialize an opt-in installed-API qualification checkout.
Run once on pristine pinned sources. Do not use on a working tree with edits.
No generated source is committed by this script.
"""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parents[2]
d=Path(__file__).resolve().parent
paths=list((root/"src").rglob("*.c"))
old={p:p.read_bytes() for p in paths}
for script in ("qualify_all_native.py","qualify_native_lifecycle.py","qualify_native_feast_lifecycle.py"):
    subprocess.run([sys.executable,str(d/script)],check=True,cwd=root)
for path,includes in {
 "src/esolver/lis_esolver_feast.c":["lis_rci_native_session.c"],
 "src/solver/lis_solver.c":["lis_rci_native_lifecycle.c"],
 "src/esolver/lis_esolver.c":["lis_rci_native_feast_lifecycle.c","lis_rci_native_public.c"]
}.items():
    p=root/path
    with p.open("a") as f:
        for name in includes:f.write('\n#include "../../qualification/rci-restore/'+name+'"\n')
for p,original in old.items():
    changed=p.read_bytes()
    if changed!=original:
        p.write_bytes(b'#ifdef HAVE_CONFIG_H\n#include "lis_config.h"\n#endif\n#ifdef USE_RCI\n'+changed+b'\n#else\n'+original+b'\n#endif\n')
for name in ("Makefile.am","Makefile.in"):
    p=root/"include"/name
    s=p.read_text()
    if name=="Makefile.am":
        a="nobase_include_HEADERS += lis_rci.h lis_rci_bridge.h"
    else:
        a="@ENABLE_RCI_TRUE@am__append_1 = lis_rci.h lis_rci_bridge.h"
    assert s.count(a)==1
    s=s.replace(a,a+" lis_rci_native.h")
    if name=="Makefile.in":
        a="\tlis_rci_bridge.h\n"
        assert s.count(a)==1
        s=s.replace(a,"\tlis_rci_bridge.h lis_rci_native.h\n")
    p.write_text(s)
print("INSTALLED NATIVE API QUALIFICATION SOURCE READY")
