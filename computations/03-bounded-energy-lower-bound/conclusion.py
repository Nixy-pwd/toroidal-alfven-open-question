"""Recheck the two parities, then derive the certified continuous residual floor."""
from fractions import Fraction as F
from pathlib import Path
import tempfile
import subprocess
from certify import check

def main():
    root=Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory() as tmp:
        lib=Path(tmp)/'dots.so'
        subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(root/'dots.cpp'),'-o',str(lib)],check=True)
        results=[check(root/f'witness-{p}.npz',lib) for p in ('even','odd')]
        endpoints=(F(41490985,10**8),F(41491200,10**8))
        endpoint_results=[[check(root/f'witness-{p}.npz',lib,frequency=f) for p in ('even','odd')] for f in endpoints]
    slack=min(v[0] for v in results)
    operator=max(v[1] for v in results)
    C=2492; a=F(3,10); b=F(7,10); L=26
    assert C*(1+2*a+F(1,2))<5234
    assert F(99,100)*C*(1-a+F(1,2))>2960
    assert 1/(1+a)-C*a*a/(b*(L+1)*(L+3))>F(9,25)
    bound=F(slack,2*5234)*2960*F(3,5)/operator
    print('continuous lower bound:',bound,'approximately',float(bound),flush=True)
    rounded=F(455,10**10)
    assert bound>rounded
    width=F(5,10**8);lam=F(41491,100000)
    assert (2*lam*width+width**2)<rounded
    print('uniform exclusion at energy <= 2492 for |frequency - 0.41491| <= 0.00000005')
    edge_slack=min(result[0] for edge in endpoint_results for result in edge)
    edge_norm=max(result[1] for edge in endpoint_results for result in edge)
    uniform=F(edge_slack,2*5234)*2960*F(3,5)/edge_norm
    assert uniform>0
    print('uniform lower bound on [0.41490985, 0.41491200]:',uniform,'approximately',float(uniform))
    gamma=F(slack,5234)
    k=2960*gamma/(1-gamma)
    energy=(C+k/(1+a))/(1+k*a*a/(b*(L+1)*(L+3)))
    print('necessary exact mode energy at center >=',float(energy))
    assert energy>2580

if __name__=='__main__': main()
