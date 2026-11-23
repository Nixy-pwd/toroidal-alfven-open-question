"""Optional finite-time bound for the full viscous and resistive linear model."""
from fractions import Fraction as F
from math import isqrt
import sys
from verify import read_trial, polynomial_weighted_norm, verify


def curl_upper(x, toroidal, z, power):
    # Vector components are (x, i*toroidal, z)/R^power, n=1.
    c1=z-toroidal.derivative(True).R()
    c2=toroidal.derivative(False).R()+toroidal.times(1-power)-x
    c3=(x.derivative(True)-z.derivative(False)).R()+z.times(power)
    _,upper=polynomial_weighted_norm([c1,c2,c3],F(2*power+1,2),12)
    return upper


def sqrt_upper(value):
    scale=10**12
    target=value.numerator*scale**2
    k=isqrt(target//value.denominator)
    return F(k+(k*k*value.denominator<target),scale)


def check(path):
    _,error=verify(path)
    x,z,_=read_trial(path)
    lam=F(41491,100000)
    rt=(x.derivative(False)+z.derivative(True)).R()
    mass,_=polynomial_weighted_norm([x,z,rt],F(1,2),6)
    numerator=curl_upper(x,rt,z,1)+curl_upper(x,rt-x.times(2),z,3)/lam**2
    form=numerator/(mass*(2-error/lam**2))
    form_ceil=(form.numerator+form.denominator-1)//form.denominator
    print('Normalized full-state curl energy <=',form_ceil,flush=True)
    rate=sqrt_upper(error**2/(2*lam**2-error))
    for eps in [F(1,10**8),F(1,10**9),F(1,10**10)]:
        def allowed(periods):
            time=F(44*periods,7)/lam
            return rate*time+sqrt_upper(eps*time*form_ceil/2)<=F(1,10)
        n=0
        while allowed(n+1):n+=1
        print('epsilon',eps,': full dissipative evolution within 10% through',n,'whole periods',flush=True)


if __name__=='__main__':check(sys.argv[1])
