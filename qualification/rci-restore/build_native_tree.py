#!/usr/bin/env python3
"""Generate native RCI only in an independent recursive build directory."""
from pathlib import Path
import re, shutil, subprocess, sys, tempfile
src=Path(sys.argv[1]).resolve()
dest=Path(sys.argv[2]).resolve()
if dest==src or dest in src.parents or dest.name!="rci-generated":
    raise SystemExit("invalid generated directory")
module=Path("qualification/rci-restore")
config=Path(sys.argv[3]).read_text()
complex_scalar=bool(re.search(r"^#define _COMPLEX 1$",config,re.M))
sources=["src/solver/lis_solver.c","src/solver/lis_solver_cg.c","src/solver/lis_solver_bicgstab.c","src/solver/lis_solver_gmres.c","src/solver/lis_solver_gcrodr.c","src/solver/lis_solver_fgcrodr.c","src/esolver/lis_esolver.c","src/esolver/lis_esolver_block.c","src/esolver/lis_esolver_feast.c","src/esolver/lis_krylov.c"]
public={"lis_rci_native_"+x for x in ("create","destroy","start_linear","start_feast","start_krylov","advance","cancel","pending","status","execute_default")}
pattern=re.compile(r"(?m)^(?P<prefix>(?:extern\s+)?(?:LIS_INT|int|void|lis_rci_[A-Za-z0-9_]+)\s*\*?\s*)(?P<name>lis_[A-Za-z0-9_]+)(?=\s*\()")
def hide(m):
    name=m["name"]
    if name in public: return m[0]
    if name.startswith(("lis_rci_native_","lis_rci_krylov_restore_")) or name.endswith("_native_reference") or name in {"lis_eigen_block_rank_transform","lis_eigen_block_zero_outputs","lis_rci_restore_harmonic_native"}:
        return '__attribute__((visibility("hidden"))) '+m[0]
    return m[0]
with tempfile.TemporaryDirectory(prefix=".rci-build-",dir=dest.parent) as tmp:
    root=Path(tmp)/"tree"
    for name in sources:
        target=root/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(src/name,target)
    shutil.copytree(src/module,root/module,ignore=shutil.ignore_patterns("__pycache__"))
    if complex_scalar:
        p=root/"src/solver/lis_solver_fgcrodr.c"
        s=p.read_text()
        for old,new in [("state->c[keep],\n                    state->c[j],\n                    &alpha","state->c[j],\n                    state->c[keep],\n                    &alpha"),("v[i1v],\n                    v[k],\n                    &t","v[k],\n                    v[i1v],\n                    &t")]:
            assert s.count(old)==1,old
            s=s.replace(old,new)
        p.write_text(s)
    for script in ("qualify_all_native.py","qualify_native_lifecycle.py","qualify_native_feast_lifecycle.py"):
        subprocess.run([sys.executable,str(root/module/script)]+(["--complex"] if complex_scalar else []),cwd=root,check=True)
    for path,names in {
       "src/esolver/lis_esolver_feast.c":["lis_rci_native_session.c"],
       "src/solver/lis_solver.c":["lis_rci_native_lifecycle.c"],
       "src/esolver/lis_esolver.c":["lis_rci_native_feast_lifecycle.c","lis_rci_native_public.c"]
    }.items():
        with (root/path).open("a") as f:
            for name in names:f.write('\n#include "../../qualification/rci-restore/'+name+'"\n')
    for file in list((root/"src").rglob("*.c"))+list((root/module).glob("*.c"))+list((root/module).glob("*.h"))+list((root/module).glob("*.inc")):
        file.write_text(pattern.sub(hide,file.read_text()))
    if dest.exists():shutil.rmtree(dest)
    root.rename(dest)
print("NATIVE RCI BUILD TREE READY")
