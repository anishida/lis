#!/usr/bin/env python3
"""Pinned continuation of the native split-complex FEAST projector only."""
import re
import generate_linear_continuation as g

path = "src/esolver/lis_esolver_feast.c"
sha = "906c28dd648156682a80b49d2cbc4ca4d0c1bf0e"
name = "lis_efeast_projector_apply_split"
src = g.get_source(path,sha)
params,body = g.extract(g.strip_comments(src),name)
fields = {}
for p in g.args(params):
    m = re.fullmatch(r"(.*?)\s*([A-Za-z_]\w*)",p)
    assert m,p
    fields[m[2]] = (m[1].strip(),"")
decl = re.compile(r"^[ \t]*(LIS_REAL|LIS_VECTOR|LIS_INT)\s+(\w+)(\[\d+\])?\s*;",re.M)
def hoist(m):
    assert m[2] not in fields,m[2]
    fields[m[2]] = (m[1],m[3] or "")
    return ""
body = decl.sub(hoist,body)
assert len(fields)==21,fields
lex=re.compile(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r"->|\.|[A-Za-z_]\w*|[^\s]")
prev=""
pos=0
out=[]
for m in lex.finditer(body):
    out.append(body[pos:m.start()])
    tok=m[0]
    out.append("F->"+tok if tok in fields and prev not in (".","->") else tok)
    prev=tok
    pos=m.end()
body="".join(out)+body[pos:]
pat=re.compile(r"\bF->err\s*=\s*(lis_matvec|lis_efeast_shifted_solve_split)\s*\(")
pos=0
out=[]
sites=[]
for m in pat.finditer(body):
    end=g.closing(body,m.end()-1,"(",")")
    args=g.args(body[m.end():end])
    tail=end+1
    while body[tail].isspace():tail+=1
    assert body[tail]==";"
    i=len(sites)+1
    if m[1]=="lis_matvec":
        assert len(args)==3
        req="fp_suspend(F,FP_B,"+args[0]+",NULL,0,0,"+args[1]+",NULL,"+args[2]+",NULL,NULL)"
    else:
        assert len(args)==9
        req="fp_suspend(F,FP_SHIFT,"+",".join(args)+")"
    out.append(body[pos:m.start()])
    out.append("F->pc="+str(i)+"; return "+req+";\nresume_"+str(i)+": F->pc=0; F->err=F->completion;")
    pos=tail+1
    sites.append(m[1])
body="".join(out)+body[pos:]
assert sites==["lis_matvec","lis_matvec","lis_efeast_shifted_solve_split","lis_efeast_shifted_solve_split"],sites
assert not pat.search(body)
assert not re.search(r"\b(lis_matvec|lis_efeast_shifted_solve_split)\s*\(",body)
body=re.sub(r"\breturn\s+([^;]+);",lambda m:m[0] if m[1].startswith("fp_suspend(") else "return fp_finish(F,"+m[1]+");",body)
out=["/* Generated from pinned native LIS FEAST projector source. */",
"struct lis_rci_native_feast_projector { int pc,active,status,strict,forced; LIS_INT completion,native; uint64_t token; lis_rci_native_feast_projector_request request;"]
out += ["    "+t+" "+n+a+";" for n,(t,a) in fields.items()]
out += ["};","static int fp_finish(lis_rci_native_feast_projector *,LIS_INT);",
"static int fp_suspend(lis_rci_native_feast_projector *,int,LIS_MATRIX,LIS_MATRIX,LIS_REAL,LIS_REAL,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER);",
"#undef __FUNC__",'#define __FUNC__ "lis_efeast_projector_apply_split"',
"static int fp_step(lis_rci_native_feast_projector *F){",
"switch(F->pc){case 0:break;"+"".join("case "+str(i)+":goto resume_"+str(i)+";" for i in range(1,5))+"default:return fp_finish(F,LIS_BREAKDOWN);}",body,"}"]
(g.DIR/"feast_projector_continuation.inc").write_text("\n".join(out)+"\n")
a="LIS_INT "+name+"("
assert src.count(a)==1
src=src.replace(a,"LIS_INT lis_efeast_projector_native_reference(")
# Insert after all definitions so the private vector compatibility helper
# and the original shift solver remain in the same translation unit.
(g.ROOT/path).write_text(src+'\n#include "../../qualification/rci-restore/lis_rci_native_feast_projector.c"\n')
print("NATIVE FEAST PROJECTOR GENERATED",sites)
