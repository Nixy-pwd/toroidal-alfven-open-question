"""Compute an upper estimate on the continuum residual using an enriched pressure trial space."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.special import eval_jacobi, roots_legendre
from alf_disk import a, n, evaluate

def pressure_residual(p, odd, y, mu, ext=10):
    """Embed velocity coefficients unchanged in an enriched quadrature calculation and fit a scalar pressure."""
    m = p + ext

    def index(deg, reflection):
        ans = []
        for k in range(deg + 2):
            spins = [1, -1] if k else [1]
            ids = [(s, j) for s in spins for j in range((deg - abs(k + s)) // 2 + 1)]
            if not ids:
                break
            ans.extend(((k, s, j) for s, j in ids if j))
            if k and (1, 0) in ids and ((-1, 0) in ids) or (not k and reflection):
                ans.append((k, 0, 0))
        return ans
    ind = {k: i for i, k in enumerate(index(m, odd))}
    u = np.zeros(len(ind))
    u[[ind[k] for k in index(p, odd)]] = y
    U, V, T, R, grad, integ, norm, curl, h = evaluate(m, odd, u)
    q, w = roots_legendre(m + 18)
    r = (q + 1) / 2
    rw = w / 2 * r * a * a
    nang = R.shape[1]
    th = np.arange(nang) * 2 * np.pi / nang
    c = np.cos(th)
    s = np.sin(th)
    B = n * n * R ** (-4) - mu
    hr = (B + 4 * R ** (-4)) * U - 2 * n * R ** (-4) * T
    hz = B * V
    ht = B * T - 2 * n * R ** (-4) * U
    fr = R * (hr * c + hz * s)
    fa = R * (-hr * s + hz * c)
    to_four = lambda z: np.fft.rfft(z, axis=1) / nang
    fr = to_four(fr)
    fa = to_four(fa)
    ft = to_four(ht)
    fr = -fr.imag if odd else fr.real
    fa = fa.real if odd else -fa.imag
    ft = -ft.imag if odd else ft.real
    W = {k: to_four(R ** k).real * rw[:, None] for k in (1, -1)}

    def inn(k, i, j, cos):
        z = W[k]
        return ((z[:, abs(i - j)] + (1 if cos else -1) * z[:, i + j]) * 0.5)[:, None]
    bases = []
    for k in range(int(odd), m + 1):
        j = np.arange((m - k) // 2 + 1)[None, :]
        t = r[:, None]
        J = eval_jacobi(j, 0, k, 2 * t * t - 1)
        Jp = eval_jacobi(np.maximum(j - 1, 0), 1, k + 1, 2 * t * t - 1) * (j > 0)
        val = t ** k * J
        der = (k * t ** (k - 1) * J + 2 * (j + k + 1) * t ** (k + 1) * Jp) / a
        val[:, 1:] -= val[:, :-1].copy()
        der[:, 1:] -= der[:, :-1].copy()
        ang = val * (1 if odd else -1) * k / (a * t)
        bases.append((k, val, der, ang))
    sizes = [z[1].shape[1] for z in bases]
    offs = np.cumsum([0] + sizes)
    rows = []
    cols = []
    dat = []
    rhs = []
    for i, (k, v, d, b) in enumerate(bases):
        rhs.extend(d.T @ (rw * fr[:, k]) + b.T @ (rw * fa[:, k]) + n * v.T @ (rw * ft[:, k]))
        for l in range(i, min(i + 25, len(bases))):
            j, Vj, Dj, Bj = bases[l]
            mat = n * n * v.T @ (inn(-1, k, j, not odd) * Vj)
            if l - i <= 1:
                mat += d.T @ (inn(1, k, j, not odd) * Dj) + b.T @ (inn(1, k, j, odd) * Bj)
            row = np.repeat(np.arange(offs[i], offs[i + 1]), sizes[l])
            col = np.tile(np.arange(offs[l], offs[l + 1]), sizes[i])
            rows.append(row)
            cols.append(col)
            dat.append(mat.ravel())
            if l != i:
                rows.append(col)
                cols.append(row)
                dat.append(mat.ravel())
    A = sp.coo_matrix((np.concatenate(dat), (np.concatenate(rows), np.concatenate(cols))), shape=(offs[-1], offs[-1])).tocsc()
    sol = sla.spsolve(A, np.array(rhs))
    err1 = hr * c + hz * s
    err2 = -hr * s + hz * c
    err3 = ht.copy()
    for i, (k, v, d, b) in enumerate(bases):
        coeff = sol[offs[i]:offs[i + 1]]
        C = np.sin(k * th) if odd else np.cos(k * th)
        S = np.cos(k * th) if odd else np.sin(k * th)
        err1 -= (d @ coeff)[:, None] * C
        err2 -= (b @ coeff)[:, None] * S
        err3 -= n * (v @ coeff)[:, None] * C / R
    return (np.sqrt(integ(err1 ** 2 + err2 ** 2 + err3 ** 2) / norm), h, curl, norm)
if __name__ == '__main__':
    for sp0 in sys.argv[1:]:
        p = int(sp0)
        z = np.load('alf_result_%d_1.npz' % p)
        j = int(z['selected'])
        print(p, pressure_residual(p, 1, z['y'][:, j], z['e'][j]), flush=True)
