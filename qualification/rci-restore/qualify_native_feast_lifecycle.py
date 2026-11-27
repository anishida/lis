#!/usr/bin/env python3
"""Lift the native single-attempt FEAST setup/finalization around the RCI core."""
import re
import generate_linear_continuation as g
p="src/esolver/lis_esolver.c"
src=g.get_source(p,"3ba869be3d881ce93ff46d418b3cdcbf06978243")
_,body=g.extract(g.strip_comments(src),"lis_gesolve_once")
import sys
g.COMPLEX_ONLY = '--complex' in sys.argv
g.REAL_ONLY = not g.COMPLEX_ONLY
body=g.no_quad(body)
fields={"A":"LIS_MATRIX","B":"LIS_MATRIX","x":"LIS_VECTOR","evalue0":"LIS_SCALAR *","esolver":"LIS_ESOLVER",
 "comm":"LIS_Comm","nesolver":"LIS_INT","niesolver":"LIS_INT","nigesolver":"LIS_INT","emaxiter":"LIS_INT",
 "evalue":"LIS_SCALAR *","evector":"LIS_VECTOR *","resid":"LIS_REAL *","rhistory":"LIS_REAL *","iter":"LIS_INT *","iter2":"LIS_INT *",
 "err":"LIS_INT","output":"LIS_INT","ss":"LIS_INT","mode":"LIS_INT","time":"double","oshift":"LIS_SCALAR","ishift":"LIS_SCALAR",
 "estorage":"LIS_INT","eblock":"LIS_INT","A0":"LIS_MATRIX","eprecision":"LIS_INT","rval":"LIS_INT","xx":"LIS_VECTOR"}
body=body.split("LIS_DEBUG_FUNC_IN;",1)[1].replace("LIS_DEBUG_FUNC_OUT;","")
a="err = lis_esolver_execute[nesolver](esolver);"
assert body.count(a)==1
body=body.replace(a,"return ef_core(K);\nresume_core: err=K->completion;")
# The external driver owns callback completion status. Do not borrow the
# synchronous callback-scope TLS across suspension or inspect an outer scope.
body=body.replace("if(lis_callback_scope_status()) err=lis_callback_scope_status();","")
body=body.replace("lis_callback_scope_status() ||","K->forced ||")
assert "lis_callback_scope_status" not in body
assert not re.search(r"\blis_(matvec|psolve|gesolve|esolve)\s*\(",body)
lex=re.compile(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r"->|\.|[A-Za-z_]\w*|[^\s]")
prev="";pos=0;out=[]
for m in lex.finditer(body):
    out.append(body[pos:m.start()]);tok=m[0]
    out.append("F->"+tok if tok in fields and prev not in (".","->") else tok)
    prev=tok;pos=m.end()
body="".join(out)+body[pos:]
(g.DIR/"native_feast_lifecycle_frame.inc").write_text("struct ef_frame {\n"+"\n".join(" "+t+" "+n+";" for n,t in fields.items())+"\n};\n")
(g.DIR/"native_feast_lifecycle_step.inc").write_text(
 "static LIS_INT ef_step(lis_rci_native_feast_lifecycle *K){\n struct ef_frame *F=&K->f; LIS_INT i;\n if(K->phase)goto resume_core;\n"+body+"\n}\n")
print("NATIVE FEAST FULL LIFECYCLE GENERATED")
