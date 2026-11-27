"""Independent certificate of original CCX pencil, using the private scaled inputs."""
import sys
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
binary,data,reference,output=map(Path,sys.argv[1:5])
def read(name):
    with (data/(name+'.csr')).open('rb') as f:
        n,z=np.fromfile(f,np.int32,2);p=np.fromfile(f,np.int32,n+1);j=np.fromfile(f,np.int32,z);v=np.fromfile(f,np.float64,z)
    return csr_matrix((v,j,p),shape=(n,n))
A=read('A');B=read('B');n=A.shape[0];raw=np.fromfile(binary).reshape(-1,n+1)[:6];lam=raw[:,0];v=raw[:,1:].T;ix=np.argsort(lam);lam=lam[ix];v=v[:,ix];av=A@v;bv=B@v;d=np.load(data/'d.npy');av0=av/d[:,None];bv0=bv/d[:,None];res=np.linalg.norm(av0-bv0*lam,axis=0)/(np.linalg.norm(av0,axis=0)+abs(lam)*np.linalg.norm(bv0,axis=0));C=v.T@av;D=np.sqrt(np.diag(C));C=C/D[:,None]/D[None,:];ortho=np.max(abs(C-np.eye(len(lam))));ref=np.loadtxt(reference,delimiter=',',skiprows=1)[:,0];delta=abs(lam*1e9-ref[:len(lam)])/ref[:len(lam)];np.savetxt(output,np.c_[lam*1e9,res,delta],delimiter=',',header='lambda,original_res,relative_error_vs_reference',comments='');ok=len(lam)==6 and max(res)<=1e-8 and max(delta)<=1e-8 and ortho<=1e-8;print('certificate',ok,'pairs',len(lam),'max_original_residual',max(res),'K_orthogonality',ortho,'max_relative_eigenvalue_error',max(delta));sys.exit(0 if ok else 1)
