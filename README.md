# Does a regular Alfvén eigenmode exist in a circular torus?

> **Open question, not a claimed mathematical breakthrough.**
> This repository makes AI-assisted, unreviewed mathematical notes and reproducible computational artifacts available to specialists. **Neither the claimed proofs nor novelty have received independent expert review or Lean 4 certification.**

## The question

Consider **linearized incompressible MHD** about the static current-free toroidal magnetic field \(B_0=\Phi\mathbf e_\phi/R\) in an axisymmetric solid torus, with boundary conditions defined in the [original note](research/01-operator-spectrum-original.md). Fix the *circular* meridional cross-section

\[
D=\{(R,Z):(R-1)^2+Z^2<(3/10)^2\},\qquad n=1,
\]

and the **exact** frequency \(\lambda_0=41491/100000\), \(\mu_0=\lambda_0^2\). The candidate operator on divergence-free wall-tangent fields in the \(e^{i\phi}\) Fourier sector is

\[
K=P\left(R^{-4}\begin{pmatrix}5&2i&0\\-2i&1&0\\0&0&1\end{pmatrix}\right)P,
\]

where \(P\) is the Leray projection in physical-volume \(L^2\).

**Does there exist a nonzero \(v\in H^1\) with \(Kv=\mu_0v\)?** Equivalently, is there a *fixed finite* \(C\) and a sequence of unit admissible fields \(v_j\) such that

\[
\|\nabla v_j\|_{L^2}^2\le C,\quad\|(K-\mu_0)v_j\|_{L^2}\to0?
\]

The roughness measure is the **full 3D cylindrical-coordinate gradient energy**, including coordinate-basis dependence. The equivalence depends on compactness and the stated operator/domain setup.

**We do not know the yes/no answer.** This is the *fixed circular torus*, **not** the other smooth boundaries specially chosen to admit an eigenmode. It does not solve global well-posedness of nonlinear 3D MHD.

## What is available?

| Material | Scope | Verification status |
|---|---|---|
| [Initial operator and spectrum](research/01-operator-spectrum-original.md) | Linearized MHD operator, spectral bands, dissipation | Candidate analysis; proofs partially sketched |
| [Embedded-mode note](research/02-embedded-eigenmodes.md) | Smooth eigenmodes on **tailored noncircular** domains; finite-expansion exclusions on the *circle* | Proposed analytical arguments, **not** full H¹ solution |
| [Turning-chord note](research/03-turning-chord.md) | Elliptic/hyperbolic interface and \(k^{4/3}\) trace-energy scaling | Proposed analytic bound, not externally checked |
| [Ray-dynamics handoff](research/04-chess-billiard-handoff.md) | Numerical reflection-map rotation and its limitations | **Floating-point only; original scripts not archived** |
| [Ritz/dissipation scripts](computations/01-ritz-and-dissipation/README.md) | Ideal/discrete spectral numerics | Numerical evidence |
| [Finite-time quasimode certificate](computations/02-finite-time-quasimode/README.md) | At \(\lambda_0\), claimed \(E\le2492\) and residual \(\le8.9983\cdot10^{-6}\), 430-cycle finite-time bound | Supplied exact-arithmetic checks, independent proof audit pending |
| [Bounded-energy lower-bound certificate](computations/03-bounded-energy-lower-bound/README.md) | For \(E\le2492\), claimed minimum residual \(>4.55\cdot10^{-8}\); necessary true-mode energy \(\ge2584.30\) | Supplied integer-checking program/witnesses, independent proof audit pending |
| [Prior-art research](prior-art/2026-11-23-audit.md) | Comparison with classic/modern MHD, inertial-wave and PDE work | **Novelty not established**; key full texts missing |

The residual upper/lower bounds concern a **fixed energy budget**; a small residual and a long-lived finite-time oscillation are **not** proofs of an exact eigenfunction.

## Reproducing a subset of the evidence

Python 3 is needed; the large lower-bound certificate additionally requires NumPy, SciPy and a C++ compiler supporting 128-bit integers. Follow each source README for full dependencies and caveats.

\`\`\`sh
cd computations/02-finite-time-quasimode
python3 check_identities.py 210
python3 consequences.py trial.json
\`\`\`

\`\`\`sh
cd computations/03-bounded-energy-lower-bound
python3 check_identities.py 105
OPENBLAS_NUM_THREADS=2 python3 conclusion.py
\`\`\`

The latter uses large NPZ witness arrays and may be expensive to run. Exact original ZIP SHA-256 hashes are in [computations/README.md](computations/README.md). The distributed files are **extracted** from the original source archives to permit line-by-line inspection. The two large binary witness arrays are included in the full project bundle; see [DATA-AVAILABILITY.md](DATA-AVAILABILITY.md) for the status of the initial GitHub mirror and the one-step completion procedure.

## Feedback welcome

We are hoping to learn whether any specific result is **correct, previously known, or worth investigating**. Particularly useful responses would identify a flaw, an exact prior theorem, a boundary/domain issue, or a rigorous result settling the circular H¹ question at this precise frequency. Corrections are as valuable as confirmations. Please open an Issue; there is **no request to referee the entire collection**.

## Origin and scope

This work began as a curiosity-driven attempt to explore a difficult question in three-dimensional resistive MHD. With AI tools assisting the derivations, computations, code, and literature exploration, the investigation developed into the narrower spectral question documented here. The materials are shared in the hope that they may be useful to researchers working in this area. No solution, priority, or novelty claim is made; the mathematical arguments and computations remain open to independent expert scrutiny. See [ORIGIN.md](ORIGIN.md) and [STATUS.md](STATUS.md).

**Caution:** Some original notes use labels like \`VERIFIED\` for *within-session* checks. These are not independent specialist review, and the reader should inspect the actual arguments rather than trust the labels.

## Reuse

The repository contents are offered under the permissive [MIT License](LICENSE), to the extent they are copyrightable. Researchers may critique, adapt, improve, and use these materials without further involvement from the initiator. The mathematical ideas are not claimed as proprietary.
