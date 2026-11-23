"""Produce dual witnesses for excluding normalized eigenfunctions under a gradient budget."""
import argparse
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from strengthened import build,null_basis,a,b_weight,tangent,ycoord
from pathlib import Path

def run(P,odd,C0,dest):
    freq=.41491
    B,M,degs,(s,c),fs,fc=build(P,odd,freq)
    Q=null_basis(B,P)
    L=26 if C0<=3000 else int(np.sqrt((1+2*a)*C0*a*a/b_weight))+4
    low=np.flatnonzero(degs<=L)
    high=np.flatnonzero(degs>L)
    assert P>L
    Z=b_weight*degs*(degs+2)/a**2
    # M was stabilized only at shell P. Low/high split instead at shell L.
    ML=M[low][:,low].toarray()+np.diag(a*(degs[low]==L))
    MQ=Q[low].T@ML@Q[low]+(1+2*a)*Q[high].T@Q[high]
    D=Q.T@(Z[:,None]*Q)
    parts=np.split(Q,[s.size,2*s.size])
    for part,(cm,fac) in zip(parts,[fs,fs,fc]):
        rhs=cm.T@part
        D+=rhs.T@fac.solve(rhs)
    eig,y=la.eigh((D+D.T)/2,(MQ+MQ.T)/2,subset_by_index=(0,0))
    print('relaxed',P,odd,L,eig[0],flush=True)
    vc=Q@y[:,0]
    us=np.split(vc,[s.size,2*s.size])
    grad=[];gg=[]
    for v,(cm,fac) in zip(us,[fs,fs,fc]):
        w=fac.solve(cm.T@v)
        grad.append(w);gg.append(cm@w)
    g=np.concatenate(gg)
    grnorm=np.sqrt(g@vc)
    g/=grnorm
    grad=[v/grnorm for v in grad]
    del Q,D,MQ,parts
    delta=C0/2
    H=C0*ML+delta*np.eye(len(low))
    d=Z.copy();d[low]+=delta;d[high]-=(1+2*a)*C0
    assert min(d)>0
    n,k=B.shape
    br=sp.diags(1/np.sqrt(d))@B
    # Equilibrate the constraint tests; their scale has no mathematical significance.
    cs=1/np.sqrt(np.asarray(br.power(2).sum(axis=0)).ravel())
    br=br@sp.diags(cs)
    system=sp.bmat([[sp.eye(n),br],[br.T,None]],format='csc')
    lu=sla.splu(system)
    def project(u):
        b=np.vstack((u/np.sqrt(d)[:,None],np.zeros((k,u.shape[1]))))
        sol=lu.solve(b)
        sol+=lu.solve(b-system@sol)
        return sol[:n],cs[:,None]*sol[n:]
    gp,zg=project(g[:,None])
    coeff=np.empty((k,len(low)));eta=np.empty(len(low))
    cols=np.empty((n,len(low)))
    for i in range(0,len(low),16):
        sl=slice(i,min(i+16,len(low)))
        rhs=np.zeros((n,min(len(low)-i,16)));rhs[low]=H[:,sl]
        u,z=project(rhs)
        yi=(gp.T@u)/(1+gp.T@gp)
        coeff[:,sl]=z-zg@yi
        eta[sl]=yi.ravel()
        cols[:,sl]=(rhs-B@coeff[:,sl]-g[:,None]*yi)/np.sqrt(d)[:,None]
    gram=cols.T@cols+np.outer(eta,eta)
    margin=la.eigh(H-(gram+gram.T)/2,subset_by_index=(0,0),eigvals_only=True)[0]
    print('dual',P,odd,margin,'relative',la.eigh(gram,H,subset_by_index=(len(low)-1,len(low)-1),eigvals_only=True)[0],'coeff',np.max(np.abs(coeff)),flush=True)
    np.savez(dest,coeff=coeff,eta=eta,g0=grad[0],g1=grad[1],g2=grad[2],P=P,odd=odd,C0=C0,L=L,margin=margin)
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('P',type=int);ap.add_argument('C',type=int);ap.add_argument('dest',type=Path)
    ap.add_argument('--even',action='store_true')
    opt=ap.parse_args();run(opt.P,not opt.even,opt.C,opt.dest)
