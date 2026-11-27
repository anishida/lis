#!/usr/bin/env python3
"""Qualification only: lift the actual native linear solve kernel's lifecycle.
The classical function is untouched. Limit the external candidate to the
configuration explicitly checked by its entry point; no unsupported fallback.
"""
import re
import generate_linear_continuation as g
p="src/solver/lis_solver.c"
src=g.get_source(p,"4feb18f8dbe399fee5cc3848b4edb77ce9a060d8")
_,body=g.extract(g.strip_comments(src),"lis_solve_kernel")
import sys
g.COMPLEX_ONLY = '--complex' in sys.argv
g.REAL_ONLY = not g.COMPLEX_ONLY
body=g.no_quad(body)
# Original top-level locals are explicit persistent fields. The loop index
# is only live within non-suspending loops and remains an OpenMP local.
fields={
 "A":"LIS_MATRIX","b":"LIS_VECTOR","x":"LIS_VECTOR","solver":"LIS_SOLVER","precon":"LIS_PRECON",
 "comm":"LIS_Comm","nsolver":"LIS_INT","precon_type":"LIS_INT","maxiter":"LIS_INT","maxiter_noimp":"LIS_INT",
 "err":"LIS_INT","rhistory":"LIS_REAL *","xx":"LIS_VECTOR","output":"LIS_INT","scale":"LIS_INT",
 "conv_cond":"LIS_INT","precision":"LIS_INT","is_use_at":"LIS_INT","storage":"LIS_INT","block":"LIS_INT",
 "n":"LIS_INT","shell_scale":"LIS_INT","p_c_time":"double","p_i_time":"double","itime":"double",
 "nrm2":"LIS_REAL","tol":"LIS_REAL","tol_w":"LIS_REAL","t":"LIS_VECTOR","residual_err":"LIS_INT",
 "bb":"LIS_VECTOR","shell_b":"LIS_VECTOR","AA":"LIS_MATRIX","B":"LIS_MATRIX","shell_A":"LIS_MATRIX",
 "Ah":"LIS_MATRIX","buf":"char","r":"LIS_VECTOR","z":"LIS_VECTOR"}
assert "LIS_DEBUG_FUNC_IN;" in body
body=body.split("LIS_DEBUG_FUNC_IN;",1)[1].replace("LIS_DEBUG_FUNC_OUT;","")
# Avoid holding the native debug stack across an external suspension.
a="err = lis_solver_execute[nsolver](solver);"
assert body.count(a)==1
body=body.replace(a,"return life_core(K);\nresume_core: err=K->completion;")
a=re.compile(r"residual_err\s*=\s*shell_scale\s*\?\s*lis_matvec\(A,x,t\)\s*:\s*lis_matvec\(A,xx,t\);")
body,n=a.subn("return life_residual(K);\nresume_residual: residual_err=K->completion;",body)
assert n==1
assert not re.search(r"\blis_(matvec|psolve)\s*\(",body)
# On external cancellation/error do not launch certification callbacks.
a="xx->precision = LIS_PRECISION_DEFAULT;"
assert body.count(a)==1
body=body.replace(a,a+"\n if(K->forced)goto residual_done;")
lex=re.compile(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r"->|\.|[A-Za-z_]\w*|[^\s]")
prev="";pos=0;out=[]
for m in lex.finditer(body):
    out.append(body[pos:m.start()])
    tok=m[0]
    out.append("F->"+tok if tok in fields and prev not in (".","->") else tok)
    prev=tok;pos=m.end()
body="".join(out)+body[pos:]
out=["/* Qualification-only generated native solve lifecycle. */","struct life_frame {"]
out += [" "+t+" "+n+("[64]" if n=="buf" else "")+";" for n,t in fields.items()]
out += ["};"]
(g.DIR/"native_lifecycle_frame.inc").write_text("\n".join(out)+"\n")
(g.DIR/"native_lifecycle_step.inc").write_text(
    "static LIS_INT life_step(lis_rci_native_lifecycle *K){\n struct life_frame *F=&K->f; LIS_INT i;\n"
    " if(K->phase==1)goto resume_core;\n if(K->phase==2)goto resume_residual;\n"+body+"\n}\n")
print("NATIVE FULL LINEAR LIFECYCLE GENERATED: core and final certification suspension")
