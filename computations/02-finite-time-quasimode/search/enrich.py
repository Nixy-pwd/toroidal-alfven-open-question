"""Embed previously computed trials exactly in the next disk polynomial space."""
from pathlib import Path
import argparse
import numpy as np
from continue_modes import DiskPolynomials,velocity_maps

def velocity_keys(p):
    return velocity_maps(p,DiskPolynomials(p+1,True),DiskPolynomials(p+1,False))[2]

def enrich(root,source,target):
    target_idx={key:i for i,key in enumerate(velocity_keys(target))}
    indices=[target_idx[key] for key in velocity_keys(source)]
    for trial in sorted(root.glob(f'v_{source}_*.npz')):
        vector=np.zeros(len(target_idx))
        vector[indices]=np.load(trial)['v']
        np.savez(root/f'v_{target}_from_{trial.stem}.npz',v=vector)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('directory',type=Path)
    ap.add_argument('source',type=int)
    ap.add_argument('target',type=int)
    args=ap.parse_args()
    enrich(args.directory,args.source,args.target)
