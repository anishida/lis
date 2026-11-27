#!/usr/bin/env python3
"""Generate an entire native real-double FEAST continuation, from pinned LIS sources.
No hidden operator/shift callbacks in the generated outer iteration.
Private synchronous block helpers are exposed only in the disposable CI checkout.
"""
import re
import generate_linear_continuation as g
import sys
g.COMPLEX_ONLY = '--complex' in sys.argv
g.REAL_ONLY = not g.COMPLEX_ONLY
p="src/esolver/lis_esolver_feast.c"
b="src/esolver/lis_esolver_block.c"
e="src/esolver/lis_esolver.c"
sources={
 p:g.get_source(p,"906c28dd648156682a80b49d2cbc4ca4d0c1bf0e"),
 b:g.get_source(b,"b8812f250cc3027ec8fd7c913e2fce306a0b7d94"),
 e:g.get_source(e,"3ba869be3d881ce93ff46d418b3cdcbf06978243")}
specs=[
 ("projector","lis_efeast_projector_apply_split",p,4),
 ("bgram","lis_eigen_block_bgram",b,1),
 ("brank","lis_eigen_block_brank_deflate",b,1),
 ("filter","lis_efeast_filtered_block_basis",p,2),
 ("rr","lis_efeast_projected_rayleigh_ritz",p,1),
 ("standard","lis_esolver_standard_residual",e,1),
 ("generalized","lis_esolver_generalized_residual",e,2),
 ("cert","lis_efeast_certify_ritz_pairs",p,2),
 ("feast","lis_efeast",p,3)]
ids={name:i+1 for i,(name,*_) in enumerate(specs)}
children={func:name for name,func,_,_ in specs}
primitive={"lis_matvec","lis_matvec_block","lis_efeast_shifted_solve_split"}
params={}
frames={}
for name,func,path,_ in specs:
    par,body=g.extract(g.strip_comments(sources[path]),func)
    fields={}
    for arg in g.args(par):
        m=re.fullmatch(r"(.*?)\s*([A-Za-z_]\w*)(\[\])?",arg)
        assert m,(name,arg)
        fields[m[2]]=(m[1].strip()+(" *" if m[3] else ""),"")
    params[name]=list(fields)
    body=g.no_quad(body)
    decl=re.compile(r"^[ \t]*(LIS_MATRIX|LIS_VECTOR|LIS_SCALAR|LIS_REAL|LIS_INT|LIS_SOLVER|LIS_ESOLVER|size_t)\s+([^;]+);",re.M)
    def hoist(m):
        out=[]
        for arg in g.args(m[2]):
            n=re.fullmatch(r"(\**)\s*(\w+)(\[\d+\])?\s*(?:=\s*(.*))?",arg,re.S)
            assert n,(name,arg)
            assert n[2] not in fields,(name,n[2])
            fields[n[2]]=(m[1]+(" "+n[1] if n[1] else ""),n[3] or "")
            if n[4]:out.append(n[2]+" = "+n[4]+";")
        return "\n".join(out)
    body=decl.sub(hoist,body)
    lex=re.compile(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r"->|\.|[A-Za-z_]\w*|[^\s]")
    prev="";pos=0;out=[]
    for m in lex.finditer(body):
        out.append(body[pos:m.start()])
        tok=m[0]
        out.append("F->"+tok if tok in fields and prev not in (".","->") else tok)
        prev=tok;pos=m.end()
    body="".join(out)+body[pos:]
    frames[name]=(fields,body)
# These synchronous private helpers remain shared with the original block
# routines; export only inside the qualification build, never in upstream API.
helpers=[]
block=sources[b]
for helper in ("lis_eigen_block_rank_transform","lis_eigen_block_zero_outputs"):
    par,_=g.extract(g.strip_comments(block),helper)
    helpers.append("extern LIS_INT "+helper+"("+par+");")
    block,n=re.subn(r"\bstatic\s+LIS_INT\s+"+helper+r"\s*\(","LIS_INT "+helper+"(",block)
    assert n==1,(helper,n)

calls="|".join(sorted(set(children)|primitive))
pat=re.compile(r"\bF->err\s*=\s*("+calls+r")\s*\(")
out=["/* Generated, qualification only. Native source copyrights retained. */",*helpers]
for name,(fields,_) in frames.items():
    out.append("struct nf_"+name+" { int pc,parent;")
    out.extend(" "+t+" "+n+a+";" for n,(t,a) in fields.items())
    out.append("};")
out.append("struct lis_rci_native_feast { int active,status,strict,forced; LIS_INT completion,native; uint64_t token; LIS_ESOLVER esolver; LIS_MATRIX owner_a,owner_b; LIS_UNSIGNED_INT as, av, bs, bv; lis_rci_native_feast_request request; LIS_VECTOR one_in,one_out;")
out.extend(" struct nf_"+name+" "+name+";" for name in frames)
out+=["};",
      "static int nf_return(lis_rci_native_feast *,int,LIS_INT);",
      "static int nf_block(lis_rci_native_feast *,int,LIS_MATRIX,LIS_INT,LIS_VECTOR *,LIS_VECTOR *);",
      "static int nf_one(lis_rci_native_feast *,int,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR);",
      "static int nf_shift(lis_rci_native_feast *,LIS_MATRIX,LIS_MATRIX,LIS_REAL,LIS_REAL,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER);"]
counts={}
for name,func,path,expected in specs:
    fields,body=frames[name]
    pos=0;parts=[];sites=[]
    for m in pat.finditer(body):
        if m.start()<pos:continue
        end=g.closing(body,m.end()-1,"(",")")
        a=g.args(body[m.end():end])
        tail=end+1
        while body[tail].isspace():tail+=1
        assert body[tail]==";"
        site=len(sites)+1;callee=m[1]
        code=["F->pc="+str(site)+";"]
        if callee in children:
            child=children[callee]
            assert len(a)==len(params[child]),(name,callee,a)
            code += ["K->"+child+"."+field+" = "+val+";" for field,val in zip(params[child],a)]
            code += ["K->"+child+".pc=0; K->"+child+".parent="+str(ids[name])+";",
                     "K->active="+str(ids[child])+"; return 3;"]
        elif callee=="lis_efeast_shifted_solve_split":
            assert len(a)==9
            code+=["return nf_shift(K,"+",".join(a)+");"]
        else:
            op="NF_B" if a[0]=="F->B" else "NF_A"
            if callee=="lis_matvec":
                assert len(a)==3
                code+=["return nf_one(K,"+op+","+",".join(a)+");"]
            else:
                assert len(a)==4
                code+=["return nf_block(K,"+op+","+",".join(a)+");"]
        parts.append(body[pos:m.start()])
        parts.append("\n".join(code)+"\nresume_"+str(site)+": F->pc=0; F->err=K->completion;")
        sites.append(callee);pos=tail+1
    body="".join(parts)+body[pos:]
    assert len(sites)==expected,(name,sites)
    assert not re.search(r"\b("+calls+r")\s*\(",body),(name,"hidden call")
    body=re.sub(r"\breturn\s+([^;]+);",lambda m:m[0] if m[1]=="3" or m[1].startswith(("nf_shift(","nf_one(","nf_block(")) else "return nf_return(K,F->parent,"+m[1]+");",body)
    counts[name]=len(sites)
    out += ["#undef __FUNC__", '#define __FUNC__ "'+func+'"',
        "static int nf_step_"+name+"(lis_rci_native_feast *K){",
        "struct nf_"+name+" *F=&K->"+name+";",
        "switch(F->pc){case 0:break;"+"".join("case "+str(i)+":goto resume_"+str(i)+";" for i in range(1,len(sites)+1))+"default:return nf_return(K,F->parent,LIS_BREAKDOWN);}",body,"}"]
out.append("static int nf_step(lis_rci_native_feast *K){switch(K->active){")
out.extend("case "+str(ids[name])+":return nf_step_"+name+"(K);" for name in frames)
out.append("default:return nf_return(K,0,LIS_BREAKDOWN);}}")
(g.DIR/"feast_native_continuation.inc").write_text("\n".join(out)+"\n")
src=sources[p]
a="LIS_INT lis_efeast(\n        LIS_ESOLVER esolver)"
assert src.count(a)==1
src=src.replace(a,"LIS_INT lis_efeast_native_reference(LIS_ESOLVER esolver)")
(g.ROOT/p).write_text(src+'\n#include "../../qualification/rci-restore/lis_rci_native_feast.c"\n')
(g.ROOT/b).write_text(block)
print("NATIVE FULL FEAST GENERATED",counts)
