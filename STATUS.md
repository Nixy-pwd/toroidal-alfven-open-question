# Verification status and open questions

**All assertions are provisional as of November 23, 2026.** The notes and programs have *not* undergone full independent mathematical review or Lean kernel verification.

| Claim | Evidence provided | Key limitation |
|---|---|---|
| Proposed linearized projected MHD operator | Derivation / algebraic checks in initial note | Operator domain and boundary conditions require expert audit |
| Smooth embedded eigenmode on certain tailored cross-sections | Analytic construction in `research/02-embedded-eigenmodes.md` | **Does not establish such a mode on the given circle** |
| No finite vertical polynomial / finite vertical Fourier ansatz eigenmode on the circle | Analytic exclusion argument | Does not exclude arbitrary infinite expansions |
| On the circle, quasimode residual ≤8.9983×10⁻⁶ and E≤2492 | Rational trial coefficients and executable verifier | Does not establish exact point spectrum; proof steps still need external verification |
| On the circle, inf residual >4.55×10⁻⁸ for E≤2492 | Integer-assisted certificates / witnesses | Bound depends on fixed budget; source computation has not been formally proved |
| A true circle eigenmode would require E≥2584.30 | Computer-assisted necessary-condition calculation | **Necessary, not sufficient** |
| Reflection rotation number ~0.5528347408 | User-supplied Claude floating-point account | Cannot certify irrationality or full elliptic–hyperbolic transfer; original scripts absent |
| Turning-chord trace extension energy ~k^(4/3) | Proposed analytic argument | A restricted elliptic-side extension problem, not the full spectral problem |
| Novelty of above model-specific claims | Prior-art survey | **Unknown**; numerous related historical papers |

**Central open question.** At the exact parameter \(\mu_0=(41491/100000)^2\), is there a nonzero H¹ eigenfunction of the proposed operator K in the fixed circular torus at n=1? Equivalently, can the normalized residual be taken to zero while the full three-dimensional gradient energy remains bounded by *some* finite constant, fixed across the sequence?

There is no claimed yes/no result in this archive.


**Distribution note.** The initial GitHub commit contains the source and text artifacts but not the two large `.npz` witness arrays. The full bundle contains them; see `DATA-AVAILABILITY.md` and `PUBLISH.md` for the completion step.
