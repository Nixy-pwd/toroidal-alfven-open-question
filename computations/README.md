# Computation archive

These three directories contain extracted data from original AI-produced ZIP files. Each subfolder README explains which claims its code attempts to verify. No original ZIP duplication is included.

| Original ZIP | Path | SHA-256 of the original ZIP |
|---|---|---|
| `alfven-reproduction.zip` | `computations/01-ritz-and-dissipation` | `7ef436455c0051bb76fda2809889ac7fff06e1f54cd6b0fd79b06c966c47addd` |
| `alfven-bounded-mode-certificate.zip` | `computations/02-finite-time-quasimode` | `bde420e96ea3e65de80953a5e222e8188976d923e109aeea58063f1e3a0c8bb1` |
| `alfven-residual-lower-bound.zip` | `computations/03-bounded-energy-lower-bound` | `8de7c5d6ad73ab60a562670d83bcc321b1af97d5b65be726b2c5007898605aee` |

Code or certificates being executable does **not** guarantee that underlying PDE/domain reasoning is mathematically sound. Source links to private ChatGPT Library files were replaced by local links where possible. For the bounded energy certificate the two NPZ witnesses are large (~36 MB combined); use git push, rather than GitHub's browser file-upload, if encountering an individual-file size limit.
