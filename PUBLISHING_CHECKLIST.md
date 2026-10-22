# Referee reproduction checklist — immutable release v0.1.1

This document refers to the public `v0.1.1` release. All instructions below operate on the versioned archive rather than a moving `main` branch.

1. Read `docs/PROOF_MAP.md` and understand the mathematical assumptions connecting the tracial SOS certificate to the published E3 theorem.
2. Install GMP/C++17 and Python dependencies as documented in `README.md` and `requirements-tested.txt`.
3. Check release integrity: `sha256sum -c MANIFEST_SHA256.txt`.
4. Reproduce the fast exact proof-data verification: `bash RUN_PRECOMPUTED.sh`.
5. Run `bash RUN_EXACT.sh` for the slower reconstruction of all 356 coefficient equations from **shared** canonicalization rules, independent field arithmetic, and exact certificates. Check `sha256sum -c MANIFEST_SHA256.txt` again afterward.
6. The independent 2-to-3 saturation is checked symbolically; physical scripts provide examples, not universal proofs.
7. Citation: use the immutable GitHub release `v0.1.1` and its pinned commit until a Zenodo DOI is actually available. Do not cite unpublished internal drafts from this repository.

No mathematical coefficient, E3 bound, or published input certificate has been altered by the v0.1.1 maintenance patch. The v0.1.0 release is preserved unchanged.
