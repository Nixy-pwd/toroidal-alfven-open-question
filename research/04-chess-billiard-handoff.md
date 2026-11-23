# Ray-dynamics continuation — Claude reply, summarized from chat

**Source:** A later Claude response supplied by the user; `john.py`, `j7.py` and `j5.py` were **not** in the available Library artifacts. Numerical claims here are NOT reproducible from this repository; treat this as a hypothesis / handoff only.

The turning chord in the fixed circular section is claimed at

\[
R_T=[\lambda_0(1+\sqrt 2)]^{-1/2}\approx0.99916.
\]

The hyperbolic-side characteristic curves were modeled as \(Z\mp G(R)=\text{constant}\), \(G'(R)=\sqrt{-d(R)}\), and a chess-billiard map switches families on reflection from the wall and turning chord.

Weighted rotation-number estimates near exact λ₀ were reported as `0.5528347408146` (20,000 iterates) and `0.5528347408193` (40,000 iterates); a nearby 68/123 locking tongue was also reported. These are **floating-point**, not interval-certified, and cannot prove irrationality, a Diophantine bound or absence of an H¹ eigenmode.

The cited analogue is Colin de Verdière–Li, *Internal Waves in 2D Domains with Ergodic Classical Dynamics* (https://arxiv.org/abs/2306.13834), whose no-point-spectrum mechanism under irrational rotation number does **not directly apply** to this variable-coefficient, mixed elliptic/hyperbolic Alfvén system. In particular, the turning chord transports amplitude and flux, not only ray paths. An asserted transfer relation `(I − b*) f = K_λ f` with a smoothing remainder was not proved, and a compact perturbation alone is not enough to supply a Fredholm determinant in the presence of small divisors.

See [turning-chord.md](03-turning-chord.md) for a later analytic note on elliptic-side energy and the transformed singular potential `7/(36s²)`.
