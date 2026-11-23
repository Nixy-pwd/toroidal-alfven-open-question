"""Write the exact rational trial; round free coefficients before enforcing tangency."""
import sys,json
from pathlib import Path
from collections import defaultdict
import numpy as np
from continue_modes import DiskPolynomials,velocity_maps
p,cap=map(int,sys.argv[1:3])
path=Path(sys.argv[3])
u=np.load(path/f'cap_{p}_{cap}.npz')
den=2**72
def rounded(x):
    num,div=float(x).as_integer_ratio()
    return round(__import__('fractions').Fraction(num*den,div))
_,_,keys=velocity_maps(p,DiskPolynomials(p+1,True),DiskPolynomials(p+1,False))
x,z,q=defaultdict(int),defaultdict(int),defaultdict(int)
for (m,s,j),v in zip(keys,u['v'],strict=True):
    a=rounded(v)
    if s: terms=[(s,j,1),(s,j-1,-1)]
    elif m: terms=[(1,0,1),(-1,0,-1)]
    else:terms=[(1,0,1)]
    for ss,jj,w in terms:
        k=abs(m+ss)
        if k:x[k,k+2*jj]+=w*a
        z[k,k+2*jj]-=ss*w*a
qr=iter(u['q'])
for m in range(1,int(u['pressure_degree'])+1):
    for d in range(m,int(u['pressure_degree'])+1,2):
        a=rounded(next(qr))
        q[m,d]+=a
        if d>m:q[m,d-2]-=a
assert next(qr,None) is None
obj={'degree':p,'denominator':den}
for key,table in [('Fx',x),('Fz',z),('q',q)]:
    obj[key]=[[m,d,a] for (m,d),a in sorted(table.items()) if a]
(path/f'trial_{p}_{cap}.json').write_text(json.dumps(obj,separators=(',',':')))
