"""Minimize measured residual within the computed trial span at each physical gradient cap."""
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from physical_features import pressure_residual

freq=0.41491
root=Path(sys.argv[1])
degrees=list(map(int,sys.argv[2:]))
budgets=(1790,1800,2490,2500,4000,8000)

def constrained(B,H,cap):
    dmin=eigh(H,eigvals_only=True)[0]
    if dmin>cap:
        return None
    def vector(tau):
        _,y=eigh(B+tau*H,subset_by_index=(0,0))
        return y[:,0]
    def cost(y):return y@H@y
    lo,hi=0.,1e-15
    if cost(vector(0))<=cap:return vector(0)
    while cost(vector(hi))>cap:
        hi*=4
    for _ in range(55):
        mid=(hi+lo)/2
        if cost(vector(mid))>cap:lo=mid
        else:hi=mid
    return vector(hi)

for p in degrees:
    paths=sorted(root.glob(f'v_{p}_*.npz'))
    assert paths
    vdata=[];pdata=[];vs=[];hs=[];rs=[]
    for path in paths:
        v=np.load(path)['v']
        V,H,R,pr=pressure_residual(p,1,v,freq**2)
        vdata.append(v);pdata.append(pr)
        vs.append(V);hs.append(H);rs.append(R)
    def gram(items):
        mat=np.stack(items)
        return mat@mat.T
    mass=gram(vs);E=gram(hs);B=gram(rs)
    vals,w=eigh(mass)
    ix=vals>1e-13*max(vals)
    w=w[:,ix]/np.sqrt(vals[ix])
    for cap in budgets:
        y=constrained(w.T@B@w,w.T@E@w,cap)
        if y is None:continue
        coeff=w@y
        v=np.array(vdata).T@coeff
        q=np.array(pdata).T@coeff
        residual=np.linalg.norm(np.array(rs).T@coeff)
        energy=np.linalg.norm(np.array(hs).T@coeff)**2
        print('bounded',p,cap,residual,energy,flush=True)
        np.savez(root/f'cap_{p}_{cap}.npz',v=v,q=q,pressure_degree=p+10,degree=p)
