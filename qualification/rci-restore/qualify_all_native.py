#!/usr/bin/env python3
"""Apply all eight qualified native dispatches to one disposable checkout."""
from pathlib import Path
import subprocess
import sys
root=Path(__file__).resolve().parents[2]
module=Path(__file__).resolve().parent
# Separate interpreter per generator: GCRO/FGCRO configure shared module globals.
# The initial-residual source is read, never rewritten.
for name in ("linear","gcro","fgcro","krylov","feast"):
    subprocess.run([sys.executable,str(module/("qualify_"+name+"_native.py"))]+(["--complex"] if "--complex" in sys.argv and name in ("gcro","fgcro","feast") else []),cwd=root,check=True)
print("ALL EIGHT NATIVE CORE PATCHES READY")
