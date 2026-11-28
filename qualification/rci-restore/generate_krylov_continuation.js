#!/usr/bin/env node
/* Source-specific continuation generator; deliberately refuse unknown upstream blobs. */
const fs=require("fs"),crypto=require("crypto");
function generate(src){
 const names=["apply_metric","verify_inverse","apply_operator","orthogonalize","append_range_seed","append_vector","append_extra_seeds","eigenpair_residual","lis_krylov_symmetric_core"];
 const ids=Object.fromEntries(names.map((n,i)=>[n,i+1]));
 const fail=(x)=>{throw Error(x)};
 function split(s){let depth=0,parts=[],start=0;for(let i=0;i<s.length;i++){let c=s[i];if("([{".includes(c))depth++;if(")]}".includes(c))depth--;if(c===","&&!depth){parts.push(s.slice(start,i).trim());start=i+1;}}parts.push(s.slice(start).trim());return parts;}
 function closing(s,start,open,close){let depth=0,quote="";for(let i=start;i<s.length;i++){let c=s[i];if(quote){if(c==="\\")i++;else if(c===quote)quote="";continue;}if(c==="'"||c==='"'){quote=c;continue;}if(c===open)depth++;if(c===close&&!--depth)return i;}fail("unbalanced "+start);}
 src=src.replace(/"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|\/\*[\s\S]*?\*\/|\/\/[^\n]*/g,m=>m.startsWith("/")?m.replace(/[^\n]/g," "):m);
 let frames=[];
 for(const name of names){
  let match=new RegExp("static int "+name+"\\s*\\(([^;]*?)\\)\\s*\\{").exec(src);
  if(!match)fail("function "+name);
  let start=match.index+match[0].lastIndexOf("{"),end=closing(src,start,"{","}");
  let body=src.slice(start+1,end);
  let params=split(match[1]).map(p=>{let m=/^(.*?)([A-Za-z_]\w*)$/.exec(p.trim());return {type:m[1].trim(),name:m[2]};});
  let fields=new Map(params.map(x=>[x.name,x.type]));
  const field=(n,t)=>{t=t.replace(/\s+/g," ").trim();if(fields.has(n)&&fields.get(n).replace(/\s+/g," ").trim()!==t)fail("shadow type "+name+" "+n);fields.set(n,t);};
  if(name==="eigenpair_residual"){
   let a=/\} else if \(p->apply_b\(p->ctx, s->n, s->x, s->rhs\)\) \{\s*return -1;\s*\}/;
   if(!a.test(body))fail("else-if anchor");
   body=body.replace(a,"} else {\n        if (p->apply_b(p->ctx, s->n, s->x, s->rhs)) return -1;\n    }");
  }
  if(name==="apply_operator"){
   let a=/if \(!err && all_finite\(s->n, y\)\)\s*err = verify_inverse\(s, metric_x, y\);/;
   if(!a.test(body))fail("conditional anchor");
   body=body.replace(a,"if (!err && all_finite(s->n, y)) {\n            err = verify_inverse(s, metric_x, y);\n        }");
  }
  const decl=/^([ \t]*)(const[ \t]+)?(struct krylov_workspace|lis_krylov_problem64|lis_krylov_matvec64|int64_t|size_t|unsigned|double|int|char)[ \t]+([^;]+);/gm;
  body=body.replace(decl,(full,indent,qual,type,rest)=>{
   let out=[];let base=(qual||"")+type;
   for(let item of split(rest)){
    let m=/^(\**)\s*([A-Za-z_]\w*)\s*(?:=\s*([\s\S]+))?$/.exec(item);
    if(!m)fail("declaration "+name+" "+item);
    let t=base+(m[1]?" "+m[1]:"");
    field(m[2],t);
    let init=m[3];
    if(init){if(init.startsWith("{"))init="("+t+")"+init;out.push(m[2]+" = "+init+";");}
   }
   return indent+(out.length?out.join("\n"+indent):";");
  });
  body=body.replace(/for\s*\(int64_t row\s*=/g,()=>{field("row","int64_t");return "for (row =";});
  // Exact pinned source: disjoint scopes with equal names/types share a slot.
  // Rename identifier tokens, but not structure members or literal contents.
  let prev="",out="",pos=0;
  let lex=/"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*'|->|\.|[A-Za-z_]\w*|[^\s]/g;
  for(let m of body.matchAll(lex)){out+=body.slice(pos,m.index);let tok=m[0];out+=(fields.has(tok)&&prev!=="."&&prev!=="->")?"F->"+tok:tok;pos=m.index+tok.length;prev=tok;}body=out+body.slice(pos);
  frames.push({name,params,fields,body});
 }
 let total=0, counts={};
 for(const frame of frames){
  let b=frame.body, sites=[],pos=0,output="";
  let callee="("+names.join("|")+"|F->p->apply_a|F->p->apply_b|F->p->solve|F->fn)";
  let re=new RegExp("^([ \\t]*)(if\\s*\\(\\s*|return\\s+|F->[A-Za-z_]\\w*\\s*=\\s*)"+callee+"\\s*\\(","gm");
  for(;;){
   let m=re.exec(b);if(!m)break;
   let from=m.index,start=re.lastIndex-1,end=closing(b,start,"(",")");
   let args=split(b.slice(start+1,end)),prefix=m[2],name=m[3],label=sites.length+1;
   let indent=m[1],code=[];
   if(ids[name]){
    const params=frames.find(f=>f.name===name).params;
    if(params.length!==args.length)fail("args "+name);
    for(let i=0;i<args.length;i++)code.push("K->f_"+name+"."+params[i].name+" = "+args[i]+";");
    code.push("K->f_"+name+".pc = 0; K->f_"+name+".parent = "+ids[frame.name]+";");
    code.push("F->pc = "+label+"; K->active = "+ids[name]+"; return K_STEP_AGAIN;");
   }else{
    if(args.length!==4)fail("callback args "+name);
    let op=name==="F->p->apply_a"?"LIS_RCI_OP_A":name==="F->p->apply_b"?"LIS_RCI_OP_B":name==="F->p->solve"?"LIS_RCI_OP_SHIFTED_SOLVE":"(F->p->transform == LIS_KRYLOV_SHIFT_INVERT_A ? LIS_RCI_OP_A : LIS_RCI_OP_B)";
    code.push("F->pc = "+label+"; return k_suspend(K, "+op+", "+args[1]+", "+args[2]+", "+args[3]+");");
   }
   code.push("resume_"+label+": F->pc = 0;");
   if(prefix.trim().startsWith("if")){
    let e=end+1;while(/\s/.test(b[e]))e++;
    if(b[e]!==")")fail("if condition "+frame.name);
    output+=b.slice(pos,from)+indent+code.join("\n"+indent)+"\n"+indent+"if (K->retval)";
    pos=e+1;
   }else{
    let e=end+1;while(/\s/.test(b[e]))e++;
    if(b[e]!==";")fail("assignment/return "+frame.name+" "+b.slice(end,e+10));
    output+=b.slice(pos,from)+indent+code.join("\n"+indent)+"\n"+indent+prefix+"K->retval;";
    pos=e+1;
   }
   re.lastIndex=pos;sites.push(label);total++;
  }
  b=output+b.slice(pos);
  if(new RegExp("\\b(?:"+names.join("|")+")\\s*\\(|F->(?:p->(?:apply_a|apply_b|solve)|fn)\\s*\\(").test(b))fail("unlowered call "+frame.name);
  b=b.replace(/\breturn\s+([^;]+);/g,(m,e)=>{
    if(e==="K_STEP_AGAIN"||e.startsWith("k_suspend("))return m;
    return "return k_return(K,F->parent,"+e+");";
  });
  frame.body=b;frame.sites=sites;counts[frame.name]=sites.length;
 }
 if(total!==28)fail("unexpected suspension count "+total+" "+JSON.stringify(counts));
 let inc="/* Generated from pinned native lis_krylov.c; Copyright (C) 2026, contributors to LIS.\n * Redistribution and use are governed by the LIS project BSD license. */\n";
 for(let f of frames){
  inc+="struct frame_"+f.name+" { int pc,parent;\n";
  for(let [n,t] of f.fields)inc+="    "+t+" "+n+";\n";
  inc+="};\n";
 }
 inc+="struct lis_rci_krylov_restore {\n int active,retval,status,native_code; uint64_t token; lis_rci_krylov_restore_request request; lis_krylov_budget budget; lis_krylov_budget_stats counts;\n";
 for(let f of frames)inc+=" struct frame_"+f.name+" f_"+f.name+";\n";
 inc+="};\n#define K_STEP_AGAIN 1\n#define K_STEP_YIELD 2\n#define K_STEP_STOP 3\n";
 inc+="static int k_return(lis_rci_krylov_restore *,int,int);\nstatic int k_suspend(lis_rci_krylov_restore *,uint32_t,int64_t,const double *,double *);\n";
 for(let f of frames){
  inc+="static int step_"+f.name+"(lis_rci_krylov_restore *K)\n{\n struct frame_"+f.name+" *F=&K->f_"+f.name+";\n switch(F->pc){case 0:break;";
  for(let label of f.sites)inc+="case "+label+":goto resume_"+label+";";
  inc+="default:return k_return(K,F->parent,LIS_KRYLOV_BREAKDOWN);}\n"+f.body+"\n}\n";
 }
 inc+="static int k_step(lis_rci_krylov_restore *K)\n{ switch(K->active){\n";
 for(let f of frames)inc+="case "+ids[f.name]+":return step_"+f.name+"(K);\n";
 inc+="default:return K_STEP_STOP;\n}}\n";
 return {inc:inc.split("\n").map(x=>x.trimEnd()).join("\n"),counts,coreId:ids.lis_krylov_symmetric_core};
}
const src=fs.readFileSync(process.argv[2]);
const gitBlob=crypto.createHash("sha1").update("blob "+src.length+"\0").update(src).digest("hex");
if(gitBlob!=="f275ffbb5cf3ab0a6757d5e59eaefd9dd9f9b74c")throw Error("Native Krylov source changed; re-audit continuation sites: "+gitBlob);
const g=generate(src.toString("utf8"));
fs.writeFileSync(process.argv[3],g.inc);
console.log(JSON.stringify(g.counts));
