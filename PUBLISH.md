# Completing the GitHub mirror

The public repository already exists at https://github.com/Nixy-pwd/toroidal-alfven-open-question and contains the initial text/source upload.

To add the two large binary witness arrays that were left out of the connector-based upload:

1. Download and unzip the full project bundle.
2. Open Terminal in the extracted `toroidal-alfven-open-question` directory.
3. Run `bash publish.sh`.
4. Complete GitHub authentication if requested.

The script clones the existing `main` branch into a temporary directory, synchronizes the complete project tree (including removals), commits any differences, and pushes them. The script itself uses an ordinary push and never force-pushes or rewrites branch history. It requires Git, rsync, and GitHub authentication on the local machine.

The largest individual witness file is approximately 29 MB, below GitHub's usual 100 MB Git-push limit. The two witness files are necessary to reproduce the bounded-energy certificate; see `DATA-AVAILABILITY.md`.
