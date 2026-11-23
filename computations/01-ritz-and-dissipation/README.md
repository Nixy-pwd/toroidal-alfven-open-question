# Reproducing the thin circular torus computation

Run from this directory. The calculations used Python 3.12, NumPy 2.3.5, SciPy 1.17.0, double precision and one OpenBLAS thread. The ideal operator and dissipative generator are assembled independently except for the common exactly solenoidal trial basis. No mesh or external dataset is required.

To reproduce the main ideal-operator table, run

```sh
export OPENBLAS_NUM_THREADS=1
python alf_cluster.py 72 96 120 144
```

To reproduce the dissipative table and individual Ritz data, run

```sh
export OPENBLAS_NUM_THREADS=1
for p in 60 72 84 96; do
  NEIG=6 python alf_disk.py "$p" 1
  python alf_decay.py "$p" 1e-7 1e-8 1e-9 1e-10 1e-11
done
python alf_res.py 60 72 84 96
```

For the individual ideal eigenvector quoted at p=144 use `NEIG=4 python alf_disk.py 144 1`. Array output files are written into the current directory. The second argument to alf_disk.py is 1 for odd radial component under vertical reflection and 0 for even. The coupled solver reports a local branch. The last field is the smallest damping ratio among the six returned eigenvalues only.

The cluster solver computes the 20 nearest ideal Ritz vectors. For each indicated frequency window it selects the linear combination minimizing the curl energy, then evaluates the normalized gradient energy and the residual against an independently fitted scalar pressure. The pressure trial space is ten degrees larger than the velocity truncation. Unrounded run output is provided in cluster.txt and decay.txt.

alf_disk.py constructs the trial basis, ideal K and mass, and evaluates physical gradients. alf_res.py performs the pressure least-squares projection. alf_decay.py constructs the full two-field generator with curl dissipation. alf_cluster.py implements the local subspace calculation.

Run `python check_polynomial.py` for an exact, independent polynomial-boundary check at the rational target frequency. The printed ranks certify the absence of finite-polynomial pressure solutions through the tested degrees over the rationals; the general proof is in the accompanying note.

At p=108, generate the ideal reference using `NEIG=4 python alf_disk.py 108 1`, then run `OPENBLAS_NUM_THREADS=2 python alf_decay_banded.py 108 1e-9`. This solver permutes by polar angular order, uses a seven-order band as a preconditioner, corrects with GMRES on the full shifted matrix, and checks the final eigenpair against the full matrix. Output is in decay_108.txt.
