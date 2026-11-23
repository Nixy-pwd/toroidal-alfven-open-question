"""Necessary polynomial moment conditions for a true circular-torus eigenfunction."""
import argparse
import time
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from continue_modes import DiskPolynomials, pressure_map

a=.3
mu=.41491**2
b_weight=.7

def ycoord(out,source):
    def entries():
        for j,(m,d) in enumerate(source.modes):
            for side in (-1,1):
                k=m+side
                beta=m if abs(k)>m else -m
                sign=-side if source.sine else side
                yield k,d+1,j,sign*(d+beta+2)/(4*(d+1))
                yield k,d-1,j,sign*(d-beta)/(4*(d+1))
    return out.matrix(entries(),source.size)


def tangent(degree, first, second, odd):
    x,z=[],[]
    col=0
    def add(m,terms):
        nonlocal col
        for spin,j,w in terms:
            k=abs(m+spin)
            x.append((k,k+2*j,col,w))
            z.append((k,k+2*j,col,(-1 if odd else 1)*spin*w))
        col+=1
    for m in range(degree+2):
        spins=(1,-1) if m else (1,)
        ids=[(s,j) for s in spins for j in range((degree-abs(m+s))//2+1)]
        for s,j in ids:
            if j: add(m,[(s,j,1),(s,j-1,-1)])
        if (1,0) in ids and (-1,0) in ids:add(m,[(1,0,1),(-1,0,-1)])
        elif not m and odd and (1,0) in ids:add(m,[(1,0,1)])
    return first.matrix(x,col),second.matrix(z,col)

def build(P,odd,freq):
    s,c=DiskPolynomials(P,odd),DiskPolynomials(P,not odd)
    rs,rc=s.radius(),c.radius()
    dx=s.derivative(s,False);dz=s.derivative(c,True)
    dxt=c.derivative(c,False);dzt=c.derivative(s,True)
    x,z=tangent(P-4,s,c,odd)
    tw=4*x+rs@(dx@x+dz@z)
    rs4=(rs@rs)@(rs@rs);rc4=(rc@rc)@(rc@rc)
    q=pressure_map(P,s)
    ev=sp.vstack((5*x-2*tw-freq**2*rs4@x,tw-2*x-freq**2*rs4@tw,z-freq**2*rc4@z))
    gr=sp.vstack((rs@(dx@q),q,rc@(dzt@q)))
    weights=np.r_[s.weights,s.weights,c.weights]
    root=np.sqrt(weights)
    B=(sp.diags(root)@sp.hstack((ev,gr))).tocsc()
    colnorm=np.sqrt(np.asarray(B.power(2).sum(axis=0)).ravel())
    B=B@sp.diags(1/colnorm)
    Bs=B.tocsc()
    degs=np.array([d for space in (s,s,c) for m,d in space.modes])
    physical=sp.diags(root)@sp.block_diag((rs,rs,rc))@sp.diags(1/root)
    M=(physical+sp.diags(a*(degs==P))).tocsc()
    def flux(space,other,odd,dx,dz,rad,rad2):
        tx,tz=tangent(P-1,space,other,odd)
        eye=sp.eye(space.size,format='csc');eye2=sp.eye(other.size,format='csc')
        X=(rad-eye)/a
        Y=ycoord(space,other)
        YD=other.derivative(space,True)
        dot=X@tx+Y@tz
        div=tx+(rad-b_weight*eye)@(dx@tx+dz@tz)+b_weight*((2/a)*dot+X@(dx@dot)+Y@(YD@dot))
        C=sp.diags(np.sqrt(space.weights))@div
        W=(tx.T@sp.diags(space.weights)@(rad-b_weight*eye)@tx+
           tz.T@sp.diags(other.weights)@(rad2-b_weight*eye2)@tz+
           b_weight*dot.T@sp.diags(space.weights)@dot)
        return C.tocsc(),sla.splu(W.tocsc())
    return B,M,degs,(s,c),flux(s,c,odd,dx,dz,rs,rc),flux(c,s,not odd,dxt,dzt,rc,rs)

def null_basis(B,P):
    n,k=B.shape;dim=n-k
    # The projector is obtained by solving exact moment constraints, not normal equations.
    mat=sp.bmat([[sp.eye(n), B],[B.T,None]],format='csc')
    lu=sla.splu(mat,permc_spec='COLAMD')
    rng=np.random.default_rng(4177)
    blocks=[]
    for j in range(0,dim+10,20):
        b=min(dim+10-j,20)
        rhs=np.zeros((n+k,b),order='F')
        rhs[:n]=rng.standard_normal((n,b))
        sol=lu.solve(rhs)
        # One correction resolves near dependencies of the ideal and pressure columns.
        sol+=lu.solve(rhs-mat@sol)
        blocks.append(sol[:n].copy())
    proj=np.hstack(blocks)
    Q,S,_=la.svd(proj,full_matrices=False,lapack_driver='gesdd')
    print('null',P,n,k,dim,S[dim-1],S[dim],np.max(np.abs(B.T@Q[:,:dim])),flush=True)
    return Q[:,:dim].copy()

def compute(P,odd,freq):
    B,M,degs,spaces,fs,fc=build(P,odd,freq)
    t=time.monotonic()
    Q=null_basis(B,P)
    s,c=spaces
    u,tor,z=np.split(Q,[s.size,2*s.size])
    az=(2*(u-tor).T@(u-tor)+z.T@z)/(1+a)
    D=az
    for part,(C,factor) in [(u,fs),(tor,fs),(z,fc)]:
        rhs=(C.T@part)
        D+=rhs.T@factor.solve(rhs)
    MQ=Q.T@(M@Q)
    tail=b_weight*(P+1)*(P+3)/((1+2*a)*a**2)
    Z=Q.T@((b_weight*degs*(degs+2)/a**2)[:,None]*Q)
    D+=Z
    theta,y=la.eigh((D+D.T)/2,(MQ+MQ.T)/2,subset_by_index=[0,0])
    lower=min(theta[0],tail)
    print('energy',P,odd,theta[0],lower,time.monotonic()-t,flush=True)
    return lower
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('degrees',nargs='+',type=int)
    ap.add_argument('--even',action='store_true')
    ap.add_argument('--freq',type=float,default=.41491)
    o=ap.parse_args()
    for p in o.degrees:compute(p,not o.even,o.freq)
