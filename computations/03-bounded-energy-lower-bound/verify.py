"""Verify the circular-torus trial certificate with integer and rational arithmetic only."""
import argparse
import json
import math
from collections import defaultdict
from fractions import Fraction as F
from functools import reduce
from math import gcd,lcm
from pathlib import Path

radius=F(3,10)
mu=F(41491,100000)**2

class Field:
    """A real Zernike polynomial with a single exact common denominator."""
    def __init__(self,terms,den,sine):
        self.coeff={k:v for k,v in terms.items() if v}
        div=reduce(gcd,self.coeff.values(),den)
        self.den=den//div
        self.sine=sine
        if div!=1:
            self.coeff={k:v//div for k,v in self.coeff.items()}
    def times(self,c):
        c=F(c)
        return Field({k:v*c.numerator for k,v in self.coeff.items()},self.den*c.denominator,self.sine)
    def __add__(self,other):
        assert self.sine==other.sine
        den=lcm(self.den,other.den)
        first,second=den//self.den,den//other.den
        terms=defaultdict(int,{k:v*first for k,v in self.coeff.items()})
        for k,v in other.coeff.items(): terms[k]+=second*v
        return Field(terms,den,self.sine)
    def __sub__(self,other): return self+other.times(-1)
    def x(self):
        # x = R-1 = (3/10) rho cos(theta), with cos/sin normalized at rho=1.
        mult=lcm(*(d+1 for m,d in self.coeff))
        terms=defaultdict(int)
        for (m,d),v in self.coeff.items():
            base=3*v*(mult//(d+1))
            for side in (-1,1):
                k=abs(m+side)
                if self.sine and k==0: continue
                sign=-1 if self.sine and m+side<0 else 1
                beta=m if k>m else -m
                terms[k,d+1]+=sign*base*(d+beta+2)
                if d-1>=k: terms[k,d-1]+=sign*base*(d-beta)
        return Field(terms,40*mult*self.den,self.sine)
    def R(self,power=1):
        out=self
        for _ in range(power):
            out=out+out.x()
        return out
    def derivative(self,vertical):
        sine=self.sine!=vertical
        terms=defaultdict(int)
        for (m,d),v in self.coeff.items():
            for side in (-1,1):
                h=m+side;k=abs(h)
                if sine and not k:continue
                sign=(-side if self.sine else side) if vertical else 1
                if h<0 and sine:sign=-sign
                for ell in range(k,d,2):
                    terms[k,ell]+=10*sign*(ell+1)*v
        return Field(terms,3*self.den,sine)
    def norm_squared(self):
        mult=lcm(*(d+1 for m,d in self.coeff))
        value=sum(v*v*(1 if m else 2)*(mult//(d+1)) for (m,d),v in self.coeff.items())
        return F(value,4*mult*self.den**2)


def polynomial_weighted_norm(fields,alpha,order):
    """Bounds on integral sum(fields^2) R^(-2 alpha), omitting one common area factor."""
    coeff=[F(1)]
    for j in range(order+1):
        coeff.append(-coeff[-1]*(alpha+j)/(j+1))
    error=abs(coeff[-1])*radius**(order+1)
    ratio=radius*max(F(1),(alpha+order+1)/(order+2))
    relative=(1+radius)**math.ceil(alpha)*error/(1-ratio)
    assert relative < 1
    value=F(0)
    for f in fields:
        g=f.times(coeff[order])
        for c in reversed(coeff[:order]):
            g=g.x()+f.times(c)
        value+=g.norm_squared()
    return value/(1+relative)**2,value/(1-relative)**2


def read_trial(path):
    obj=json.loads(Path(path).read_text())
    den=obj['denominator']
    fields=[]
    for name,sine in [('Fx',True),('Fz',False),('q',True)]:
        vals={(m,d):v for m,d,v in obj[name]}
        fields.append(Field(vals,den,sine))
    return fields


def verify(path):
    x,z,q=read_trial(path)
    # At rho=1 each radial Zernike equals 1. Check cos(theta) Fx + sin(theta) Fz.
    flux=defaultdict(int)
    for (m,d),v in x.coeff.items():
        for k in (m+1,m-1):
            if k:flux[abs(k)]+=v*z.den*(1 if k>0 else -1)
    for (m,d),v in z.coeff.items():
        for k in (1+m,1-m):
            if k:flux[abs(k)]+=v*x.den*(1 if k>0 else -1)
    assert not any(flux.values()), 'Boundary flux is nonzero'
    t=x.derivative(False)+z.derivative(True)
    rt=t.R()
    mlo,mhi=polynomial_weighted_norm([x,z,rt],F(1,2),6)
    assert mlo>0
    dx=x.derivative(False).R()-x
    dz=x.derivative(True).R()
    ex=z.derivative(False).R()-z
    ez=z.derivative(True).R()
    gt=t.derivative(False).R(2)
    gz=t.derivative(True).R(2)
    cyl=x-rt
    _,ehi=polynomial_weighted_norm([dx,dz,ex,ez,gt,gz,cyl,cyl,z],F(3,2),8)
    hx=x.times(5)-x.R(4).times(mu)-rt.times(2)-q.derivative(False).R(5)
    hz=z-z.R(4).times(mu)-q.derivative(True).R(5)
    ht=rt-rt.R(4).times(mu)-x.times(2)-q.R(4)
    _,bhi=polynomial_weighted_norm([hx,hz,ht],F(9,2),12)
    energy=ehi/mlo
    squared_error=bhi/mlo
    # Integer round-up gives rigorous decimal enclosures, independent of any float.
    upper_energy=(energy.numerator+energy.denominator-1)//energy.denominator
    places=10**10
    nr=squared_error.numerator*places**2
    k=math.isqrt(nr//squared_error.denominator)
    if k*k*squared_error.denominator<nr:k+=1
    print('gradient energy <=',upper_energy,'; normalized residual <=',f'{k}/{places}',flush=True)
    return upper_energy,F(k,places)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('trial',type=Path)
    opt=ap.parse_args()
    verify(opt.trial)
