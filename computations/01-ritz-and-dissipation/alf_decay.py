"""Assemble and solve the coupled velocity–magnetic-field problem with equal diffusivities."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.special import roots_legendre
from alf_disk import a, n, basis

def mhd(p, odd):
    """Return mass, ideal velocity-to-field coupling and the quadratic form of curl dissipation."""
    q, w = roots_legendre(p + 10)
    r = (q + 1) / 2
    rw = w / 2 * r * a * a
    wt = (-1.0) ** np.arange(len(w)) * np.sqrt((1 - q * q) * w)
    dif = r[:, None] - r[None, :]
    np.fill_diagonal(dif, 1)
    dr = wt[None, :] / wt[:, None] / dif
    np.fill_diagonal(dr, 0)
    dr -= np.diag(dr.sum(axis=1))
    dr /= a
    blocks = basis(p, odd, r)
    nang = 4 * (p + 12)
    th = np.arange(nang) * 2 * np.pi / nang
    R = 1 + a * r[:, None] * np.cos(th)
    W = {k: np.fft.rfft(R ** k, axis=1).real / nang * rw[:, None] for k in range(-3, 2)}

    def inn(k, i, j, cos):
        z = W[k]
        si = -1 if i < 0 and (not cos) else 1
        sj = -1 if j < 0 and (not cos) else 1
        i, j = (abs(i), abs(j))
        return (si * sj * (z[:, abs(i - j)] + (1 if cos else -1) * z[:, i + j]) * 0.5)[:, None]
    parts = []
    ez = -1 if odd else 1
    for m, dm, vm, tm in blocks:
        der = dr @ tm
        pm = m / (a * r[:, None]) * tm
        up = (der - pm) * 0.5
        lo = (der + pm) * 0.5
        cr = [(0, m + 1, up), (0, m - 1, lo), (-1, m, tm)]
        cz = [(0, m + 1, ez * up), (0, m - 1, -ez * lo)]
        cp = []
        rot = np.zeros_like(tm)
        for j, sg, v in vm:
            cr.append((-2, j, -n * v))
            cz.append((-2, j, -n * ez * sg * v))
            cp.append((-2, j, ez * sg * v))
            rot += ez * sg * (dr @ v + (1 + m * sg) / (a * r[:, None]) * v)
        cp.append((-1, m, -rot))
        parts.append((cr, cz, cp))
    dims = [b[1] for b in blocks]
    offs = np.cumsum([0] + dims)
    rows = []
    cols = []
    mats = [[] for j in range(3)]
    for l, (m, dm, vm, tm) in enumerate(blocks):
        for ll in range(l, min(l + 25, len(blocks))):
            k, dk, vk, tk = blocks[ll]
            mass = tm.T @ (inn(1, m, k, not odd) * tk)
            sym = n * (tm.T @ (inn(-1, m, k, not odd) * tk))
            for i, si, vi in vm:
                for j, sj, vj in vk:
                    mass += vi.T @ ((inn(-1, i, j, not odd) + si * sj * inn(-1, i, j, odd)) * vj)
                    sym += n * vi.T @ ((inn(-3, i, j, not odd) + si * sj * inn(-3, i, j, odd)) * vj)
            dij = sym.copy()
            dji = sym.copy()
            for i, si, vi in vm:
                dji -= 2 * vi.T @ (inn(-2, i, k, not odd) * tk)
            for j, sj, vj in vk:
                dij -= 2 * tm.T @ (inn(-2, m, j, not odd) * vj)
            damp = np.zeros_like(mass)
            for h, cos in enumerate((not odd, odd, odd)):
                for power, i, va in parts[l][h]:
                    temp = np.zeros((len(r), dk))
                    for power2, j, vb in parts[ll][h]:
                        temp += inn(1 + power + power2, i, j, cos) * vb
                    damp += va.T @ temp
            row = np.repeat(np.arange(offs[l], offs[l + 1]), dk)
            col = np.tile(np.arange(offs[ll], offs[ll + 1]), dm)
            rows.append(row)
            cols.append(col)
            for li, mat in zip(mats, (mass, dij, damp)):
                li.append(mat.ravel())
            if ll != l:
                rows.append(col)
                cols.append(row)
                for li, mat in zip(mats, (mass, dji, damp)):
                    li.append(mat.ravel())
    pos = (np.concatenate(rows), np.concatenate(cols))
    return [sp.coo_matrix((np.concatenate(z), pos), shape=(offs[-1], offs[-1])).tocsc() for z in mats]

def main():
    p = int(sys.argv[1])
    epses = list(map(float, sys.argv[2:]))
    G, D, C = mhd(p, 1)
    check = np.load('alf_result_%d_1.npz' % p)
    j = int(check['selected'])
    v = check['y'][:, j]
    print('formcheck', p, v @ (G @ v), v @ (C @ v), check['rs'][j, 1], flush=True)
    mass = sp.block_diag((G, G), format='csc')
    cdiag = sp.block_diag((C, C), format='csc')
    D2 = sp.bmat([[None, D.T], [D, None]], format='csc')
    vref = np.r_[G @ v, np.zeros_like(v)]
    for eps in epses:
        op = D2 - 1j * eps * cdiag
        shift = 0.41491 - 1j * eps * 2000
        lu = sla.splu((op - shift * mass).tocsc())
        inv = sla.LinearOperator(mass.shape, matvec=lambda x: lu.solve(mass @ x), dtype=complex)
        e, y = sla.eigs(inv, k=6, which='LM', tol=2e-11, v0=np.random.default_rng(7).standard_normal(len(vref)))
        fr = shift + 1 / e
        ol = np.abs(vref @ y) / np.sqrt(np.real(np.sum(np.conj(y) * (mass @ y), axis=0)))
        j = int(np.argmax(ol))
        print('branch', p, eps, fr[j].real, -fr[j].imag / eps, ol[j], 'nearby min', min(-fr.imag / eps), flush=True)
        del lu, inv
if __name__ == '__main__':
    main()
