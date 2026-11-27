import numpy as np, pathlib, sys
from scipy.sparse import csr_matrix
out=pathlib.Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True)
rng=np.random.default_rng(72117)
def write(a,path):
 a=csr_matrix(a);a.sort_indices(); n,z=a.shape[0],a.nnz
 with open(path,'wb') as f:
  np.array([n,z],np.int32).tofile(f);a.indptr.astype(np.int32).tofile(f);a.indices.astype(np.int32).tofile(f);a.data.astype(np.float64).tofile(f)
for case in range(12):
 n=[1,2,3,7,12,16,24,31,48,64,93,120][case]; prob=[0,.2,.6,.2,.35,.5,.05,.15,.1,.15,.08,.05][case]
 g=rng.random((n,n))<prob;g=np.triu(g,1);g=g|g.T
 for step in range(3):
  a=rng.uniform(-2,2,(n,n))*g; a[np.diag_indices(n)]=np.sum(np.abs(a),axis=1)+rng.uniform(0.5,3,n)
  if step==2: a[::3,:]*=-1 # indefinite in diagonal sign; nonsymmetric, nonsingular via row dominance
  write(a,out/f'case{case:02d}_{step}.csr')
 # final case with changed pattern, and deliberately singular same pattern retained in existing reuse test
 if n>=3:
  g[0,1]=g[1,0]=not g[0,1];a=rng.uniform(-2,2,(n,n))*g; a[np.diag_indices(n)]=np.sum(np.abs(a),axis=1)+1
  write(a,out/f'case{case:02d}_changed.csr')
