"""Check one finer dissipative run with a banded angular preconditioner and full-pencil residuals."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from alf_decay import mhd
from alf_disk import basis
from numpy.linalg import norm

def run(p, epses):
    G, D, C = mhd(p, 1)
    z = np.load('alf_result_%d_1.npz' % p)
    v = z['y'][:, int(z['selected'])]
    vref = np.r_[G @ v, np.zeros_like(v)]
    size = G.shape[0]
    mass = sp.block_diag((G, G), format='csc')
    dissip = sp.block_diag((C, C), format='csc')
    ideal = sp.bmat([[None, D.T], [D, None]], format='csc')
    del G, D, C
    sizes = [b[1] for b in basis(p, 1, np.array([0.4]))]
    ends = np.cumsum([0] + sizes)
    perm = np.concatenate([np.r_[np.arange(a, b), np.arange(a, b) + size] for a, b in zip(ends, ends[1:])])
    angular = np.repeat(np.arange(len(sizes)), 2 * np.array(sizes))

    def permute(mat):
        return mat[perm][:, perm].tocsc()
    mass = permute(mass)
    ideal = permute(ideal)
    dissip = permute(dissip)
    vref = vref[perm]
    scale = np.sqrt(mass.diagonal())
    for eps in epses:
        shift = 0.41491 - 1j * eps * 2000
        A = (ideal - shift * mass - 1j * eps * dissip).tocsc()
        coo = A.tocoo()
        good = np.abs(angular[coo.row] - angular[coo.col]) <= 7
        good &= np.abs(coo.data) > 1e-11 * scale[coo.row] * scale[coo.col]
        approx = sp.coo_matrix((coo.data[good], (coo.row[good], coo.col[good])), shape=A.shape).tocsc()
        del coo, good
        lu = sla.splu(approx, permc_spec='NATURAL', diag_pivot_thresh=0)
        del approx
        print('factor', p, eps, lu.nnz, flush=True)
        pre = sla.LinearOperator(A.shape, matvec=lu.solve, dtype=complex)

        def inverse(x):
            rhs = mass @ x
            ans = lu.solve(rhs)
            ans, info = sla.gmres(A, rhs, x0=ans, M=pre, rtol=1e-10, atol=0, restart=16, maxiter=4)
            if info and norm(rhs - A @ ans) > 1e-07 * norm(rhs):
                raise RuntimeError('inverse residual exceeds tolerance')
            return ans
        inv = sla.LinearOperator(A.shape, matvec=inverse, dtype=complex)
        e, y = sla.eigs(inv, k=6, which='LM', tol=2e-10, v0=np.random.default_rng(7).standard_normal(len(vref)))
        fr = shift + 1 / e
        ol = np.abs(vref @ y) / np.sqrt(np.real(np.sum(np.conj(y) * (mass @ y), axis=0)))
        j = int(np.argmax(ol))
        err = norm(A @ y[:, j] + (shift - fr[j]) * (mass @ y[:, j])) / norm(mass @ y[:, j])
        print('residual', p, eps, err, flush=True)
        print('branch', p, eps, fr[j].real, -fr[j].imag / eps, ol[j], 'nearby min', min(-fr.imag / eps), flush=True)
        del inv, pre, lu, A
if __name__ == '__main__':
    run(int(sys.argv[1]), list(map(float, sys.argv[2:])))
