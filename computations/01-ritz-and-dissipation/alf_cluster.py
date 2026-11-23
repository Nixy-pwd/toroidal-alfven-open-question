"""Measure the smoothest combination of ideal modes in a prescribed frequency window."""
import sys
import numpy as np
import scipy.linalg as la
import scipy.sparse.linalg as sla
from alf_disk import assemble
from alf_decay import mhd
from alf_res import pressure_residual


def run(degree, eigenvalue_count=20):
    target = 0.41491
    stiffness, mass = assemble(degree, 1)
    initial = np.random.default_rng(1234).standard_normal(mass.shape[0])
    eigenvalues, vectors = sla.eigsh(stiffness, k=eigenvalue_count, M=mass,
                                     sigma=target**2, v0=initial, tol=1e-10)
    del stiffness, mass
    curl_form = mhd(degree, 1)[2]
    local_curl = vectors.T @ (curl_form @ vectors)
    del curl_form
    for half_width in (6e-5, 2e-4):
        inside = np.abs(np.sqrt(eigenvalues) - target) < half_width
        _, coefficients = la.eigh(local_curl[np.ix_(inside, inside)])
        smoothest = coefficients[:, 0]
        velocity = vectors[:, inside] @ smoothest
        frequency = np.sqrt(eigenvalues[inside] @ (smoothest**2))
        upper, gradient_energy, _, _ = pressure_residual(degree, 1, velocity, target**2)
        print('cluster', degree, half_width, inside.sum(), frequency, upper, gradient_energy, flush=True)
        np.savez(f'alf_cluster_{degree}_{half_width:g}.npz',
                 v=velocity, mu=target**2, upper=upper, h=gradient_energy)


if __name__ == '__main__':
    for argument in sys.argv[1:]:
        run(int(argument))
