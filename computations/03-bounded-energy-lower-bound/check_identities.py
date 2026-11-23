"""Verify the radial multiplication and derivative identities as exact polynomials."""
from functools import lru_cache
from math import comb
import sys

@lru_cache(maxsize=None)
def radial(d,m):
    if d<m or (d-m)%2:
        return {}
    return {d-2*j:(-1)**j*comb(d-j,j)*comb(d-2*j,(d-m)//2-j)
            for j in range((d-m)//2+1)}

def check(degree):
    for m in range(degree+1):
        for step in (-1,1):
            k=abs(m+step)
            beta=m if k>m else -m
            lower_sum={}
            for d in range(m,degree+1,2):
                current=radial(d,m)
                lower=radial(d-1,k)
                higher=radial(d+1,k)
                for exponent,val in lower.items():
                    lower_sum[exponent]=lower_sum.get(exponent,0)+2*d*val
                for power in range((k%2),d+2,2):
                    derivative=(power+1-step*m)*current.get(power+1,0)
                    assert derivative==lower_sum.get(power,0)
                    product=2*(d+1)*current.get(power-1,0)
                    assert product==(d+beta+2)*higher.get(power,0)+(d-beta)*lower.get(power,0)
    print('Exact radial coefficient identities checked through degree',degree,flush=True)

if __name__=='__main__':check(int(sys.argv[1]))
