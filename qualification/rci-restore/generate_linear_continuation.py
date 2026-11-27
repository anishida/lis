#!/usr/bin/env python3
"""Pinned, source-specific conversion of native LIS linear recurrences.
The generator leaves all recurrence, residual, history and error statements
unchanged except for frame-local references and explicit yield sites.
It is deliberately not a general C parser. Refuse unreviewed source changes.
"""
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIR = Path(__file__).resolve().parent
SPECS = [
    ("initial", "lis_solver_get_initial_residual", "src/solver/lis_solver.c", "4feb18f8dbe399fee5cc3848b4edb77ce9a060d8"),
    ("cg", "lis_cg", "src/solver/lis_solver_cg.c", "18aa1489bba0b7be4a5c3d14aeb6fe88c1fef888"),
    ("bi", "lis_bicgstab", "src/solver/lis_solver_bicgstab.c", "89b228f34bb868eb356bbfefaaa1a2cb46ae7cce"),
    ("gm", "lis_gmres", "src/solver/lis_solver_gmres.c", "be206b6d00588b60dbc2d6ab4e67925ecbdf0995"),
    ("fgm", "lis_fgmres", "src/solver/lis_solver_gmres.c", "be206b6d00588b60dbc2d6ab4e67925ecbdf0995"),
]
IDS = {name:i+1 for i,(name,*_) in enumerate(SPECS)}
REAL_ONLY = False
COMPLEX_ONLY = False
EXTRA_CALLS = {}
EXPECTED_COUNTS = {"initial":2,"cg":3,"bi":5,"gm":5,"fgm":4}
def digest(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
def get_source(path, sha):
    data = (ROOT/path).read_bytes()
    assert digest(data)==sha, (path, digest(data))
    return data.decode()
def closing(s, start, op, cl):
    depth=0
    for match in re.finditer(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r'[\s\S]',s[start:]):
        c=match[0]
        if c==op: depth+=1
        elif c==cl:
            depth-=1
            if not depth: return start+match.start()
    raise ValueError("unbalanced")
def extract(src,name):
    m=re.search(r"LIS_INT\s+"+name+r"\s*\(([^;]*?)\)\s*\{",src)
    assert m, name
    start=m.end()-1
    return m[1],src[start+1:closing(src,start,"{","}")]

def strip_comments(s):
    return re.sub(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r"/\*[\s\S]*?\*/|//[^\n]*",
        lambda m: re.sub(r"[^\n]"," ",m[0]) if m[0].startswith("/") else m[0],s)

def no_quad(s):
    # Remove only known disabled preprocessor branches; preserve unrelated
    # conditionals, including any nested OpenMP directive inside live code.
    stack=[]
    out=[]
    for line in s.splitlines(True):
        m=re.match(r"\s*#\s*(ifdef|ifndef|if|else|endif)\b(.*)",line)
        parent=all(active for active,_ in stack)
        if m:
            kind,arg=m[1],m[2].strip()
            if kind in ("ifdef","ifndef","if"):
                known=(arg=="USE_QUAD_PRECISION" or (REAL_ONLY and arg in ("_COMPLEX","USE_MPI")) or (COMPLEX_ONLY and arg in ("_COMPLEX","USE_MPI"))) and kind in ("ifdef","ifndef")
                active=((kind=="ifdef") if COMPLEX_ONLY and arg=="_COMPLEX" else (kind=="ifndef")) if known else True
                stack.append((active,known))
                if parent and not known: out.append(line)
            elif kind=="else":
                active,known=stack.pop()
                stack.append((not active if known else True,known))
                if all(x for x,_ in stack[:-1]) and not known: out.append(line)
            else:
                active,known=stack.pop()
                if all(x for x,_ in stack) and not known: out.append(line)
        elif parent: out.append(line)
    assert not stack
    return "".join(out)

def args(s):
    depth=0
    result=[]
    start=0
    for i,c in enumerate(s):
        if c in "([{": depth+=1
        elif c in ")]}": depth-=1
        elif c=="," and not depth:
            result.append(s[start:i].strip()); start=i+1
    result.append(s[start:].strip())
    return result

def frame(name,func,path,sha):
    param,body=extract(strip_comments(get_source(path,sha)),func)
    body=no_quad(body)
    fields={}
    for p in args(param):
        m=re.fullmatch(r"(.*?)\s*([A-Za-z_]\w*)",p)
        assert m
        fields[m[2]]=m[1].strip()
    # Local declarations occupy the beginning of each native function.
    # No initialized automatic object is live only on the C call stack.
    decl=re.compile(r"^([ \t]*)(LIS_Comm|LIS_MATRIX|LIS_VECTOR|LIS_PRECON|LIS_SCALAR|LIS_REAL|LIS_INT|LIS_GCRODR_EIG|LIS_FGCRODR_EIG|LIS_FGCRODR_STATE|size_t|double)[ \t]+([^;]+);",re.M)
    def hoist(m):
        out=[]
        for v in args(m[3]):
            d=re.fullmatch(r"(\**)\s*(\w+)\s*(?:=\s*(.*))?",v,re.S)
            assert d,(name,v)
            assert d[2] not in fields,(name,d[2])
            fields[d[2]]=m[2]+(" "+d[1] if d[1] else "")
            if d[3]: out.append(d[2]+" = "+d[3]+";")
        return m[1]+("\n"+m[1]).join(out)
    body=decl.sub(hoist,body)
    # GMRES's original OpenMP loop variable is local to a non-suspending
    # loop. Keep it private, rather than incorrectly privatizing a field.
    if name=="gm":
        body,n=re.subn(r"#pragma omp parallel for private\(k\)(\s*#endif\s*for)\(k=0;k<n;k\+\+\)(\s*\{\s*z->value\[)k(\]\s*=\s*s->value\[0\]\*v\[0\]->value\[)k(\];\s*\})",
            r"#pragma omp parallel for\1(LIS_INT local_k=0;local_k<n;local_k++)\2local_k\3local_k\4",body)
        assert n==1
    # Rename tokens only; structure member names and string contents stay.
    lex=re.compile(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r"->|\.|[A-Za-z_]\w*|[^\s]")
    prev=""
    out=[]
    pos=0
    for m in lex.finditer(body):
        out.append(body[pos:m.start()])
        tok=m[0]
        out.append("F->"+tok if tok in fields and prev not in (".","->") else tok)
        prev=tok;pos=m.end()
    body="".join(out)+body[pos:]
    # Normal callback statements, including the native initial residual
    # helper, are lowered into continuations with the exact original return.
    calls = "|".join(["lis_matvec","lis_psolve","lis_solver_get_initial_residual",*EXTRA_CALLS])
    pattern=re.compile(r"^([ \t]*)(F->err\s*=\s*|if\s*\(\s*)("+calls+r")\s*\(",re.M)
    pos=0
    result=[]
    sites=[]
    for m in pattern.finditer(body):
        if m.start()<pos: continue
        end=closing(body,m.end()-1,"(",")")
        a=args(body[m.end():end])
        site=len(sites)+1
        indent=m[1]
        code=[]
        if m[3]=="lis_solver_get_initial_residual":
            assert len(a)==5
            for field,val in zip(("solver","M","t","r","bnrm2"),a):
                code.append("K->initial."+field+" = "+val+";")
            code += ["K->initial.pc=0; K->initial.parent="+str(IDS[name])+";",
                     "F->pc="+str(site)+"; K->active=1; return L_STEP_AGAIN;"]
            tail=end+1
            while body[tail].isspace(): tail+=1
            assert body[tail]==")"
            result.append(body[pos:m.start()])
            result.append(indent+("\n"+indent).join(code)+"\n"+indent+"resume_"+str(site)+": F->pc=0;\n"+indent+"if(K->retval)")
            pos=tail+1
        elif m[3] in EXTRA_CALLS:
            child,param_names = EXTRA_CALLS[m[3]]
            assert len(param_names)==len(a),(name,m[3],a)
            for field,val in zip(param_names,a):
                code.append("K->"+child+"."+field+" = "+val+";")
            code += ["K->"+child+".pc=0; K->"+child+".parent="+str(IDS[name])+";",
                     "F->pc="+str(site)+"; K->active="+str(IDS[child])+"; return L_STEP_AGAIN;"]
            tail=end+1
            while body[tail].isspace(): tail+=1
            assert body[tail]==";"
            result.append(body[pos:m.start()])
            result.append(indent+("\n"+indent).join(code)+"\n"+indent+"resume_"+str(site)+": F->pc=0; F->err=K->retval;")
            pos=tail+1
        else:
            assert len(a)==3
            op="LIS_RCI_OP_A" if m[3]=="lis_matvec" else "LIS_RCI_OP_PRECON"
            matrix=a[0] if m[3]=="lis_matvec" else "F->solver->A"
            code += ["F->pc="+str(site)+";",
                     "return linear_suspend(K,"+op+","+matrix+","+a[1]+","+a[2]+");"]
            tail=end+1
            while body[tail].isspace(): tail+=1
            assert body[tail]==";"
            result.append(body[pos:m.start()])
            result.append(indent+("\n"+indent).join(code)+"\n"+indent+"resume_"+str(site)+": F->pc=0; F->err=K->retval;")
            pos=tail+1
        sites.append(site)
    body="".join(result)+body[pos:]
    assert not re.search(r"\b("+calls+r")\s*\(",body),name
    body=re.sub(r"\breturn\s+([^;]+);",lambda m:m[0] if m[1]=="L_STEP_AGAIN" or m[1].startswith("linear_suspend(") else "return linear_return(K,F->parent,"+m[1]+");",body)
    # Clear pointers on original terminal cleanup to avoid stale ownership.
    body=body.replace("lis_free(F->h);","lis_free(F->h); F->h=NULL;")
    return fields,body,sites

def generate():
    frames={}
    for spec in SPECS:
        frames[spec[0]]=frame(*spec)
    counts={name:len(f[2]) for name,f in frames.items()}
    assert counts==EXPECTED_COUNTS,counts
    out=["/* Generated from native LIS solver functions. Copyright notices remain in their source files; LIS BSD terms apply. */"]
    for name,(fields,body,sites) in frames.items():
        out.append("struct linear_frame_"+name+" { int pc,parent;")
        out.extend("    "+t+" "+n+";" for n,t in fields.items())
        out.append("};")
    out.append("struct lis_rci_native_linear_restore { int active,method,status,strict; LIS_INT retval,native_code; uint64_t token,matvecs,precons; LIS_SOLVER solver; lis_rci_native_linear_request request;")
    out.extend("    struct linear_frame_"+n+" "+n+";" for n in frames)
    out.append("};")
    out.extend(["#define L_STEP_AGAIN 1","#define L_STEP_YIELD 2","#define L_STEP_STOP 3",
        "static int linear_suspend(lis_rci_native_linear_restore *,uint32_t,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR);",
        "static int linear_return(lis_rci_native_linear_restore *,int,LIS_INT);"])
    for name,func,*_ in SPECS:
        fields,body,sites=frames[name]
        out.extend(["#undef __FUNC__", '#define __FUNC__ "'+func+'"',
                    "static int linear_step_"+name+"(lis_rci_native_linear_restore *K) {",
                    "struct linear_frame_"+name+" *F=&K->"+name+";",
                    "switch(F->pc){case 0:break;"+ "".join("case "+str(i)+":goto resume_"+str(i)+";" for i in sites)+"default:return linear_return(K,F->parent,LIS_BREAKDOWN);}",
                    body,"}"])
    out.append("static int linear_step(lis_rci_native_linear_restore *K){switch(K->active){")
    out.extend("case "+str(IDS[n])+":return linear_step_"+n+"(K);" for n in frames)
    out.append("default:return L_STEP_STOP;}}")
    return "\n".join(x.rstrip() for x in "\n".join(out).splitlines())+"\n",counts

if __name__=="__main__":
    output,counts=generate()
    Path(sys.argv[1]).write_text(output)
    print("NATIVE LINEAR GENERATED",counts)
