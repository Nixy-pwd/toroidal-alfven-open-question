# Data availability and upload status

The GitHub repository contains the extracted research notes, source code, textual outputs, trial JSON, and checksums for the original archives.

The bounded-energy certificate also depends on two binary NumPy witness arrays:

- `computations/03-bounded-energy-lower-bound/witness-even.npz` (7,445,565 bytes)
- `computations/03-bounded-energy-lower-bound/witness-odd.npz` (29,874,741 bytes)

These two arrays are present in the full project bundle and their SHA-256 digests are recorded in `SHA256SUMS`. They were not included in the initial connector-based GitHub commit; the large binary arrays still need to be pushed from the full local bundle. Until `publish.sh` has been run successfully, reproduction of the bounded-energy certificate is incomplete from the GitHub repository alone.

To complete the mirror, download and extract the full project bundle, then run `bash publish.sh` from the project directory. The script clones the existing repository, synchronizes the complete project tree, and pushes any missing files. It requires Git, rsync, and GitHub authentication on the local machine. The largest individual witness file is below GitHub's usual 100 MB Git-push limit.
