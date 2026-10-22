# Publication checklist — P014/E3 referee code

1. Create a **new public repository** under `hnasreddine` named `P014-Tetrahedral-Quantum-Certificate`. Do not use or change the existing FortiSec/ARI repositories. GitHub UI: <https://github.com/new>.
2. For connector-assisted upload, initialize it with a README to create a `main` branch; then request an upload of the prepared package via the connected GitHub integration. (If uploading manually, unzip this archive and add **all its contents**, not the ZIP itself.)
3. Review file visibility, especially that no unpublished manuscript, personal contacts, internal audit, environment files or company code are present.
4. Review [`LICENSE`](LICENSE), [`LICENSE-DATA`](LICENSE-DATA), and [`docs/LICENSE_SCOPE.md`](docs/LICENSE_SCOPE.md).
5. Verify GitHub Actions `Verify E3 exact certificate` completes on the public `main` branch; you may also run `bash RUN_PRECOMPUTED.sh` locally.
6. Make an immutable tagged release (e.g., `v0.1.0`) referencing its commit SHA. If depositing on Zenodo, wait for a real DOI before updating `CITATION.cff` and the manuscript.
7. Only after the public URL and release exist, update the private Overleaf manuscript's Data Availability Statement with the exact URL, tag, and (if available) DOI.

The precomputed check verifies the supplied mathematical certificate using exact arithmetic. The optional full reconstruction is slower; the physical example scripts are supporting checks, not proofs of all resource asymptotics.
