"""Rigorous lower-bound checker. Floats choose witnesses; inequalities are verified on integers."""
import argparse,ctypes,subprocess,tempfile
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
from fractions import Fraction as F
from collections import defaultdict
from math import isqrt,ceil,lcm
from verify import Field

a=F(3,10); b=F(7,10); freq=F(41491,100000)
ONE=F(1)
BF,BZ,BR,BH=57,20,15,24
def shift(n):return 1<<n
def floorfrac(x):return x.numerator//x.denominator
def upperroot(x,scale=1<<48):
    num=x.numerator*scale**2
    n=isqrt(num//x.denominator)
    return F(n+int(n*n*x.denominator<num),scale)
def zeta(w):
    return F(isqrt(w.numerator*(1<<120)//w.denominator),1<<60)
def yy(f):
    target=not f.sine
    mult=lcm(*(d+1 for m,d in f.coeff))
    t=defaultdict(int)
    for (m,d),v in f.coeff.items():
        base=v*(mult//(d+1))
        for side in (-1,1):
            h=m+side;k=abs(h)
            if target and not k:continue
            sign=(-side if f.sine else side)*(-1 if target and h<0 else 1)
            beta=m if k>m else -m
            t[k,d+1]+=sign*base*(d+beta+2)
            if d-1>=k:t[k,d-1]+=sign*base*(d-beta)
    return Field(t,4*mult*f.den,target)
def txbasis(p,odd):
    def make(m,terms):
        x,z=defaultdict(int),defaultdict(int)
        for s,j,val in terms:
            k=abs(m+s); d=k+2*j
            if k or not odd:x[k,d]+=val
            if k or odd:z[k,d]+=(-1 if odd else 1)*s*val
        return Field(x,1,odd),Field(z,1,not odd)
    for m in range(p+2):
        spins=(1,-1) if m else (1,)
        ids=[(s,j) for s in spins for j in range((p-abs(m+s))//2+1)]
        for s,j in ids:
            if j:yield make(m,[(s,j,1),(s,j-1,-1)])
        if (1,0) in ids and (-1,0) in ids:yield make(m,[(1,0,1),(-1,0,-1)])
        elif not m and odd and (1,0) in ids:yield make(m,[(1,0,1)])
def grad_test(p,odd,coeff):
    den=1<<48
    xx,zz=defaultdict(int),defaultdict(int)
    for (x,z),val in zip(txbasis(p-1,odd),coeff,strict=True):
        w=round(F(float(val))*den)
        for k,v in x.coeff.items():xx[k]+=w*v
        for k,v in z.coeff.items():zz[k]+=w*v
    x,z=Field(xx,den,odd),Field(zz,den,not odd)
    dot=x.x().times(1/a)+yy(z)
    div0=x.derivative(False)+z.derivative(True)
    div=x+div0.R()-div0.times(b)+(dot.times(2/a)+dot.derivative(False).x().times(1/a)+yy(dot.derivative(True))).times(b)
    def inner(f,g):
        return ((f+g).norm_squared()-(f-g).norm_squared())/4
    W=a*a*(inner(x,x.R()-x.times(b))+inner(z,z.R()-z.times(b))+b*dot.norm_squared())
    return div,W

def check(path,work,frequency=freq):
    src=np.load(path)
    P,odd,C,L=int(src['P']),bool(src['odd']),int(src['C0']),int(src['L'])
    assert C==2492 and L==26 and P>L
    masks=(odd,odd,not odd)
    modes=[[(m,d) for m in range(int(s),P+1) for d in range(m,P+1,2)] for s in masks]
    rows=[(i,m,d) for i,items in enumerate(modes) for m,d in items]
    ri={k:i for i,k in enumerate(rows)}
    n=len(rows);low=[i for i,(_,m,d) in enumerate(rows) if d<=L]
    small={i:l for l,i in enumerate(low)};h=len(low)
    weights=[a*a/F(2*(d+1)*(2 if m else 1)) for _,m,d in rows]
    zt=[zeta(w) for w in weights]
    alpha=[z*z/w for z,w in zip(zt,weights)]
    delta=F(C,2)
    D=[v*(b*d*(d+2)/(a*a)+(delta if d<=L else -(1+2*a)*C)) for (_,m,d),v in zip(rows,alpha)]
    assert min(D)>0 and min(alpha)>F(99,100) and max(alpha)<=1
    hi=[upperroot(1/v,1<<40) for v in D]
    # H_L and its preweighted coordinate embedding, both formed from exact rationals.
    H0=np.zeros((h,h),dtype=np.int64)
    target=np.zeros((n,h),dtype=np.int64)
    for lj,j in enumerate(low):
        component,m,d=rows[j]; f=Field({(m,d):1},1,masks[component]).R()
        for (k,e),v in f.coeff.items():
            i=ri[component,k,e]
            if i not in small:continue
            entry=C*zt[i]*zt[j]*F(v,f.den)/weights[j]
            if i==j:entry+=alpha[i]*(delta+(C*a if d==L else 0))
            H0[small[i],lj]=round(entry*shift(BH))
            target[i,lj]=round(hi[i]*entry*shift(BH))
    assert np.array_equal(H0,H0.T)
    urows=[];ucols=[];uvals=[]
    upper_zeta=[upperroot(w,1<<60) for w in weights]
    iv=[];jv=[];av=[]
    count=np.zeros(n,dtype=np.int64);j=0;test_mass=F(0)
    def add(fields,scale):
        nonlocal j
        for component,f in enumerate(fields):
            for (m,d),v in f.coeff.items():
                i=ri[component,m,d]
                value=round(hi[i]*zt[i]*scale*F(v,f.den)*shift(BF))
                count[i]+=1
                if value: iv.append(i);jv.append(j);av.append(value)
        j+=1
    def normalized(fields):
        # The rounded scale is a free choice of a nonzero rational test factor.
        normsq=0.
        for i,f in enumerate(fields):
            for (m,d),v in f.coeff.items():
                normsq+=float(weights[ri[i,m,d]])*(v/f.den)**2
        return F(1/np.sqrt(normsq))
    for x,z in txbasis(P-4,odd):
        t=x.times(4)+(x.derivative(False)+z.derivative(True)).R()
        fields=(x.times(5)-t.times(2)-x.R(4).times(freq**2),
                t-x.times(2)-t.R(4).times(freq**2),z-z.R(4).times(freq**2))
        scale=normalized(fields)
        if frequency != freq:
            delta_mu = frequency**2 - freq**2
            fields = tuple(f - v.R(4).times(delta_mu) for f,v in zip(fields,(x,t,z)))
        # Overestimate each test norm by an integer before summation.
        test_mass+=ceil((1+a)**7*a*a*scale**2*(x.norm_squared()+t.norm_squared()+z.norm_squared()))
        for component,poly in enumerate((x,t,z)):
            for (mode,degree),v in poly.coeff.items():
                i=ri[component,mode,degree]
                urows.append(i);ucols.append(j)
                uvals.append(round(upper_zeta[i]*scale*F(v,poly.den)*shift(BF)))
        add(fields,scale)
    eigen_columns=j
    for m,d in modes[0]:
        terms={(m,d):1}
        if d>m:terms[m,d-2]=-1
        q=Field(terms,1,odd)
        fields=(q.derivative(False).R(),q,q.derivative(True).R())
        add(fields,normalized(fields))
    ds=[];W=F(0)
    for i,parity in enumerate(masks):
        g,w=grad_test(P,parity,src['g'+str(i)])
        ds.append(g);W+=w
    add(ds,F(1))
    assert W>0
    ug=upperroot(W,1<<40)
    Fmat=sp.csr_matrix((np.asarray(av,dtype=np.int64),(iv,jv)),shape=(n,j))
    coeff=np.vstack([src['coeff'],src['eta']])
    assert coeff.shape==(j,h) and np.isfinite(coeff).all()
    assert np.max(np.abs(coeff))*shift(BZ)<2**62
    witness=np.ascontiguousarray(np.rint(coeff*shift(BZ)),dtype=np.int64)
    # These *integers* define the witness, irrespective of how they were chosen.
    maxz=int(np.max(np.abs(witness)));maxt=int(np.max(np.abs(target)))
    maxa=int(np.max(np.abs(Fmat.data)));terms=int(max(count))
    tshift=BF+BZ-BH;rshift=BF+BZ-BR
    worst=(maxt<<tshift)+terms*maxz*maxa
    assert worst < 2**126
    # Check the conversion after exact accumulation inside C++, before narrowing.
    out=np.zeros((n,h),dtype=np.int64)
    lib=ctypes.CDLL(str(work));fun=lib.rounded_residual
    ptr=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS')
    fun.argtypes=[ctypes.c_int,ctypes.c_int,ptr,ptr,ptr,ptr,ptr,ctypes.c_int,ctypes.c_int,ptr]
    fun.restype=ctypes.c_int
    ret=fun(n,h,Fmat.indptr.astype('int64'),Fmat.indices.astype('int64'),Fmat.data,witness,target,tshift,rshift,out)
    assert ret==0
    extra=np.array([round(ug*F(int(z),shift(BZ))*shift(BR)) for z in witness[-1]],dtype=np.int64)
    Rint=np.vstack((out,extra))
    err=F(1,shift(BR))+F(1,shift(BH))+F(terms*maxz,shift(BF+BZ))
    err2=(n+1)*h*err**2
    guard=ceil(101*err2+F(h,shift(BH)))+1
    gram=exact_gram(Rint)
    # Common integer denominator is 100*2^(2BR); HD is a conservative minorant.
    den=100*shift(2*BR)
    HD=H0*(den//shift(BH))-101*gram
    assert max(int(np.max(np.abs(H0)))*(den//shift(BH)),101*int(np.max(np.abs(gram))))<2**60
    assert guard*den<2**60
    HD[np.diag_indices(h)]-=guard*den
    G=HD.astype(float)/den
    mine=la.eigh(G,subset_by_index=(0,0),eigvals_only=True)[0]
    print('finite verification',P,odd,'cols',j,'largest row',terms,'guard',guard,'trial pivot margin',mine,flush=True)
    assert mine>2
    s=max(1,int(.97*mine))
    lc=la.cholesky(G-s*np.eye(h),lower=True)
    li=np.asarray(np.rint(lc*shift(BR)),dtype=np.int64)
    gg=exact_gram(li.T)
    assert int(abs(HD).max())+100*int(abs(gg).max())<2**62
    rest=HD-100*gg
    assert int(abs(rest).max())*h<2**62
    gaps=2*np.diag(rest)-np.abs(rest).sum(axis=1)
    assert np.min(gaps)>0
    slack=F(int(min(gaps)),den)
    assert slack>=1
    # Conservative witness norm for a nonzero original-equation residual.
    # Bound the composed physical velocity test map, retaining cancellations
    # between different test columns before taking its Frobenius norm.
    umat=sp.csr_matrix((np.asarray(uvals,dtype=np.int64),(urows,ucols)),shape=(n,j))
    unz=int(max(np.diff(umat.indptr)))
    assert unz*int(abs(umat.data).max())*maxz<2**126
    zero=np.zeros_like(target)
    uout=np.empty_like(zero)
    assert fun(n,h,umat.indptr.astype('int64'),umat.indices.astype('int64'),umat.data,witness,zero,0,BF+BZ,uout)==0
    square_sum=0
    for col in uout.T:
        square_sum+=sum(int(k)**2 for k in col)
    eu=1+F(unz*maxz,shift(BF+BZ))
    WZ=upperroot((1+a)**7)*(upperroot(F(square_sum))+upperroot(F(n*h))*eu)
    print('certified',P,odd,'H-T >=',floorfrac(slack),'I','physical test norm <=',ceil(WZ),flush=True)
    return floorfrac(slack),ceil(WZ)

def exact_gram(R):
    R=np.ascontiguousarray(R)
    assert R.shape[0]*int(abs(R).max())**2<2**62
    col2=(R*R).sum(axis=0)
    assert int(max(col2))<2**52
    # Each dot product is a sum of integers. By Cauchy-Schwarz every partial sum
    # and product has modulus < 2^52, so IEEE binary64 BLAS is exact here.
    f=R.astype(float)
    return (f.T@f).astype(np.int64)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('witness',nargs='+',type=Path)
    opt=ap.parse_args()
    with tempfile.TemporaryDirectory() as folder:
        lib=Path(folder)/'dots.so'
        subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(Path(__file__).with_name('dots.cpp')),'-o',str(lib)],check=True)
        results=[check(p,lib) for p in opt.witness]
    print('proof margins and test norms:',results)
