"""Constrained residual search for the circular torus, in a Zernike coefficient basis."""
import argparse
import os
import csv
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from numpy.polynomial import chebyshev
from alf_disk import assemble
from alf_res import pressure_residual

A = 0.3
FREQUENCY = 0.41491


class DiskPolynomials:
    """Real sin/cos Zernike families, with unnormalized Fourier and radial factors."""
    def __init__(self, degree, sine):
        self.sine = sine
        self.modes = [(m, d) for m in range(int(sine), degree + 1)
                      for d in range(m, degree + 1, 2)]
        self.indices = {key: j for j, key in enumerate(self.modes)}
        self.weights = np.array([A*A*(0.5 if m else 1)/(2*(d+1))
                                 for m, d in self.modes])
        self.size = len(self.modes)

    def matrix(self, entries, columns):
        rows, cols, values = [], [], []
        for m, d, column, value in entries:
            idx = self.indices.get((abs(m), d))
            if idx is not None and value:
                rows.append(idx)
                cols.append(column)
                values.append(-value if m < 0 and self.sine else value)
        mat = sp.coo_matrix((values, (rows, cols)), shape=(self.size, columns)).tocsc()
        mat.eliminate_zeros()
        return mat

    def radius(self):
        def entries():
            for j, (m,d) in enumerate(self.modes):
                yield m,d,j,1
                for side in (-1,1):
                    target = m+side
                    signed = m*(1 if abs(target)>m else -1)
                    yield target,d+1,j,A*(d+signed+2)/(4*(d+1))
                    yield target,d-1,j,A*(d-signed)/(4*(d+1))
        return self.matrix(entries(), self.size)

    def derivative(self, source, vertical):
        def entries():
            for j, (m,d) in enumerate(source.modes):
                for side in (-1,1):
                    sign = (-side if source.sine else side) if vertical else 1
                    for lower in range(abs(m+side),d,2):
                        yield m+side,lower,j,sign*(lower+1)/A
        return self.matrix(entries(), source.size)


def velocity_maps(p, sine, cosine):
    x, z, keys = [], [], []
    def include(m, terms, key):
        column = len(keys)
        keys.append((m,*key))
        for spin,j,weight in terms:
            k = abs(m+spin)
            x.append((k,k+2*j,column,weight))
            z.append((k,k+2*j,column,-spin*weight))
    for m in range(p+2):
        spins = (1,-1) if m else (1,)
        ids = [(s,j) for s in spins for j in range((p-abs(m+s))//2+1)]
        for s,j in ids:
            if j:
                include(m,[(s,j,1),(s,j-1,-1)],(s,j))
        if (1,0) in ids and (-1,0) in ids:
            include(m,[(1,0,1),(-1,0,-1)],(0,0))
        elif not m and (1,0) in ids:
            include(m,[(1,0,1)],(0,0))
    return sine.matrix(x,len(keys)), cosine.matrix(z,len(keys)), keys


def pressure_map(p, sine):
    terms = []
    cols = 0
    for m,deg in sine.modes:
        if deg<=p:
            terms.append((m,deg,cols,1))
            if deg>m:
                terms.append((m,deg-2,cols,-1))
            cols+=1
    return sine.matrix(terms,cols)


def weight(matrix, radius, space, exponent, degree):
    # Polynomial weights are used only in the search. Final energies and residuals
    # are recomputed with physical R weights on an independent integration grid.
    coeff = chebyshev.chebinterpolate(lambda x:(1+A*x)**exponent,degree)
    unit = sp.eye(radius.shape[0],format='csc')
    scaled = (radius-unit)/A
    last = matrix
    total = coeff[0]*matrix
    current=scaled@matrix
    for i in range(1,degree+1):
        total += coeff[i]*current
        if i != degree:
            last,current=current,2*scaled@current-last
    return sp.diags(np.sqrt(space.weights))@total


def model(p):
    pd = p+4
    sine,cosine = (DiskPolynomials(pd+10,s) for s in (True,False))
    rs,rc = sine.radius(),cosine.radius()
    dxs,dxc = sine.derivative(sine,False),cosine.derivative(cosine,False)
    dzs,dzc = cosine.derivative(sine,True),sine.derivative(cosine,True)
    fx,fz,keys = velocity_maps(p,sine,cosine)
    t = dxs@fx+dzc@fz
    t.eliminate_zeros()
    q = pressure_map(pd,sine)
    s2,c2 = rs@rs,rc@rc
    s4,c4 = s2@s2,c2@c2
    mu = FREQUENCY**2
    targets = [
        (5*fx-mu*(s4@fx)-2*rs@t, -s4@(rs@(dxs@q)),rs,sine),
        (fz-mu*(c4@fz), -c4@(rc@(dzs@q)),rc,cosine),
        (rs@t-mu*(s4@(rs@t))-2*fx, -s4@q,rs,sine)]
    nr = len(keys)
    def gram(items):
        out = None
        for x in items:
            local = x.T@x
            out = local if out is None else out+local
        return out.tocsc()
    normal = gram(weight(sp.hstack((v,pr),format='csc'),r,s,-4.5,4)
                  for v,pr,r,s in targets)
    mass = gram(weight(v,r,s,-.5,2) for v,r,s in
                [(fx,rs,sine),(fz,rc,cosine),(rs@t,rs,sine)])
    grads = [(rs@(dxs@fx)-fx,rs,sine),(rc@(dzs@fx),rc,cosine),
             (rc@(dxc@fz)-fz,rc,cosine),(rs@(dzc@fz),rs,sine),
             (s2@(dxs@t),rs,sine),(c2@(dzs@t),rc,cosine),
             (fx-rs@t,rs,sine),(fx-rs@t,rs,sine),(fz,rc,cosine)]
    energy = gram(weight(v,r,s,-1.5,3) for v,r,s in grads)
    def extend(mat):
        return sp.block_diag((mat,sp.csc_matrix((q.shape[1],q.shape[1]))),format='csc')
    return normal,extend(energy),extend(mass),keys


def anchor():
    p=16
    k,m=assemble(p,1)
    vals,vec=sla.eigsh(k,k=1,M=m,sigma=FREQUENCY**2,
                       v0=np.random.default_rng(14).standard_normal(m.shape[0]))
    spaces=[DiskPolynomials(p+2,b) for b in (True,False)]
    _,_,keys=velocity_maps(p,*spaces)
    return dict(zip(keys,vec[:,0]))


def search(p, penalties, base, dest):
    normal,energy,mass,keys=model(p)
    unit=np.zeros(normal.shape[0])
    unit[:len(keys)]=[base.get(k,0.) for k in keys]
    print('model',p,normal.shape,normal.nnz,flush=True)
    with (dest/'bounded.csv').open('a',newline='') as f:
        writer=csv.writer(f)
        for tau in penalties:
            system=normal+tau*energy
            scales=1/np.sqrt(system.diagonal())
            diag=sp.diags(scales)
            scaled_system=(diag@system@diag).tocsc()
            threshold=float(os.getenv("ALF_DROP", "8e-16"))
            if threshold:
                scaled_system.data[np.abs(scaled_system.data)<threshold]=0
                scaled_system.eliminate_zeros()
            del system
            factor=sla.splu(scaled_system,permc_spec='MMD_AT_PLUS_A',
                             diag_pivot_thresh=.001)
            u=unit.copy()
            for _ in range(3):
                u=scales*factor.solve(scales*(mass@u))
                u/=np.sqrt(u@(mass@u))
            velocity=u[:len(keys)]
            bound,h,c,norm=pressure_residual(p,1,velocity,FREQUENCY**2)
            overlap=abs(unit@(mass@u))/np.sqrt(unit@(mass@unit))
            row=(p,tau,bound,h,overlap)
            print(*row,flush=True)
            writer.writerow(row);f.flush()
            np.savez(dest/f'v_{p}_{tau:.3e}.npz',v=velocity/np.sqrt(norm))
            del factor

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('degrees',type=int,nargs='+')
    ap.add_argument('--tau',type=float,nargs='+',default=[1e-11,3e-12,1e-12,3e-13,1e-13])
    ap.add_argument('--out',type=Path,default=Path('.'))
    opt=ap.parse_args()
    opt.out.mkdir(exist_ok=True,parents=True)
    ref=anchor()
    for degree in opt.degrees:
        search(degree,opt.tau,ref,opt.out)
