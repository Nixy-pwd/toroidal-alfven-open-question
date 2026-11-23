"""Solenoidal Jacobi disk basis and Rayleigh–Ritz calculation for the ideal operator."""
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.special import eval_jacobi, roots_legendre
a = 0.3
n = 1

def basis(p, odd, r):
    """Return normal-flux-conforming radial families, grouped by polar angle index."""
    blocks = []
    for m in range(p + 2):
        spins = [1, -1] if m else [1]
        ids = [(s, j) for s in spins for j in range((p - abs(m + s)) // 2 + 1)]
        if not ids:
            break
        c = len(ids)
        F = np.zeros((c, c))
        count = 0
        for i, (s, j) in enumerate(ids):
            if j:
                F[i, count] = 1
                F[i - 1, count] = -1
                count += 1
        if m and (1, 0) in ids and ((-1, 0) in ids):
            F[ids.index((1, 0)), count] = 1
            F[ids.index((-1, 0)), count] = -1
            count += 1
        if not m and odd:
            F[0, count] = 1
            count += 1
        if not count:
            continue
        F = F[:, :count]
        vals = []
        div = np.zeros((len(r), count))
        for s in spins:
            ix = [i for i, (ss, j) in enumerate(ids) if ss == s]
            if not ix:
                continue
            k = abs(m + s)
            t = r[:, None]
            j = np.arange(len(ix))[None, :]
            J = eval_jacobi(j, 0, k, 2 * t * t - 1)
            Jp = eval_jacobi(np.maximum(j - 1, 0), 1, k + 1, 2 * t * t - 1) * (j > 0)
            raw = t ** k * J
            d = ((k + 1 + m * s) * t ** (k - 1) * J + 2 * (j + k + 1) * t ** (k + 1) * Jp) / (a * n)
            div += d @ F[ix]
            vals.append((k, s, raw @ F[ix]))
        if odd and (not m):
            div[:] = 0
        blocks.append((m, count, vals, div))
    return blocks

def assemble(p, odd):
    """Assemble the ideal stiffness and mass using radial Gauss and angular Fourier moments."""
    q, w = roots_legendre(p + 10)
    r = (q + 1) / 2
    rw = w / 2 * r * a * a
    blocks = basis(p, odd, r)
    nang = 4 * (p + 12)
    th = np.arange(nang) * 2 * np.pi / nang
    R = 1 + a * r[:, None] * np.cos(th)
    W = {k: np.fft.rfft(R ** k, axis=1).real / nang * rw[:, None] for k in (1, -1, -3, -4, -5)}

    def inn(k, i, j, cos):
        z = W[k]
        return ((z[:, abs(i - j)] + (1 if cos else -1) * z[:, i + j]) * 0.5)[:, None]
    dims = [b[1] for b in blocks]
    offs = np.cumsum([0] + dims)
    rows = []
    cols = []
    gg = []
    hh = []
    for l, (m, dm, vm, tm) in enumerate(blocks):
        for ll in range(l, len(blocks)):
            k, dk, vk, tk = blocks[ll]
            if k - m > 26:
                break
            g = tm.T @ (inn(1, m, k, not odd) * tk)
            h = n * n * (tm.T @ (inn(-3, m, k, not odd) * tk))
            for i, si, vi in vm:
                h -= 2 * n * vi.T @ (inn(-4, i, k, not odd) * tk)
                for j, sj, vj in vk:
                    a0 = inn(-1, i, j, not odd) + si * sj * inn(-1, i, j, odd)
                    a1 = (n * n + 4) * inn(-5, i, j, not odd) + n * n * si * sj * inn(-5, i, j, odd)
                    g += vi.T @ (a0 * vj)
                    h += vi.T @ (a1 * vj)
            for j, sj, vj in vk:
                h -= 2 * n * tm.T @ (inn(-4, m, j, not odd) * vj)
            row = np.repeat(np.arange(offs[l], offs[l + 1]), dk)
            col = np.tile(np.arange(offs[ll], offs[ll + 1]), dm)
            rows.append(row)
            cols.append(col)
            gg.append(g.ravel())
            hh.append(h.ravel())
            if ll != l:
                rows.append(col)
                cols.append(row)
                gg.append(g.ravel())
                hh.append(h.ravel())
    pos = (np.concatenate(rows), np.concatenate(cols))
    return [sp.coo_matrix((np.concatenate(z), pos), shape=(offs[-1], offs[-1])).tocsc() for z in (hh, gg)]

def evaluate(p, odd, y):
    """Reconstruct a velocity on an oversampled grid, and evaluate its mass, curl and gradient energies."""
    q, w = roots_legendre(p + 18)
    r = (q + 1) / 2
    wt = (-1.0) ** np.arange(len(w)) * np.sqrt((1 - q * q) * w)
    dif = r[:, None] - r[None, :]
    np.fill_diagonal(dif, 1)
    D = wt[None, :] / wt[:, None] / dif
    np.fill_diagonal(D, 0)
    D -= np.diag(D.sum(axis=1))
    nang = 4 * (p + 15)
    th = np.arange(nang) * 2 * np.pi / nang
    c = np.cos(th)
    s = np.sin(th)
    R = 1 + a * r[:, None] * c
    fx = np.zeros_like(R)
    fz = fx.copy()
    T = fx.copy()
    off = 0
    for m, dm, v, d in basis(p, odd, r):
        b = y[off:off + dm]
        off += dm
        for k, sign, rad in v:
            b0 = (rad @ b)[:, None]
            fx += b0 * (np.sin(k * th) if odd else np.cos(k * th))
            fz += b0 * sign * (-np.cos(k * th) if odd else np.sin(k * th))
        T += (d @ b)[:, None] * (np.sin(m * th) if odd else np.cos(m * th))
    U = fx / R
    V = fz / R

    def grad(z):
        dr = D @ z / a
        da = np.fft.irfft(np.fft.rfft(z, axis=1) * (1j * np.arange(nang // 2 + 1)), n=nang, axis=1) / (a * r[:, None])
        return (c * dr - s * da, s * dr + c * da)

    def integ(z):
        return np.sum(w[:, None] / 2 * r[:, None] * a * a * R * z) / nang
    Ur, Uz = grad(U)
    Vr, Vz = grad(V)
    Tr, Tz = grad(T)
    norm = integ(U * U + V * V + T * T)
    curl = integ((n * V / R - Tz) ** 2 + (Uz - Vr) ** 2 + (Tr + T / R - n * U / R) ** 2) / norm
    H = integ(Ur ** 2 + Uz ** 2 + Vr ** 2 + Vz ** 2 + Tr ** 2 + Tz ** 2 + ((n * U - T) ** 2 + (U - n * T) ** 2 + (n * V) ** 2) / R ** 2) / norm
    return (U, V, T, R, grad, integ, norm, curl, H)

def main():
    p = int(sys.argv[1])
    odd = int(sys.argv[2])
    lam = float(sys.argv[3]) if len(sys.argv) > 3 else 0.41491
    H, G = assemble(p, odd)
    mu = lam * lam
    num = int(os.getenv('NEIG', '14'))
    e, y = sla.eigsh(H, k=num, M=G, sigma=mu, v0=np.random.default_rng(123).standard_normal(G.shape[0]), tol=1e-11)
    rs = []
    for j in range(num):
        U, V, T, R, grad, integr, norm, curl, h = evaluate(p, odd, y[:, j])
        B = n * n * R ** (-4) - e[j]
        q = R / n * (B * T - 2 * n * R ** (-4) * U)
        qr, qz = grad(q)
        res = np.sqrt(integr(((B + 4 * R ** (-4)) * U - 2 * n * R ** (-4) * T - qr) ** 2 + (B * V - qz) ** 2) / norm)
        rs.append((np.sqrt(e[j]), curl, h, res))
    rows = sorted(enumerate(rs), key=lambda kv: kv[1][1])
    print(p, odd, G.shape, [(j, tuple((round(x, 9) for x in row))) for j, row in rows[:3]], flush=True)
    j = rows[0][0]
    np.savez('alf_result_%d_%d.npz' % (p, odd), e=e, y=y, selected=j, rs=rs)
if __name__ == '__main__':
    main()
