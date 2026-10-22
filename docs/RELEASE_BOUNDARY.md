# Public/private release boundary

## Public: this repository
Only theorem E3 verification data, source code, optional certificate-generator source and finite constructed examples. Outputs are mathematical/data artifacts for referee reproducibility. The exact archive identifier is the version tag plus git commit and `MANIFEST_SHA256.txt`.

## Private: separate Overleaf manuscript project
`main.tex`, `references.bib`, journal figures, author disclosures, proof-draft notes, old frozen baseline, and internal referee/validation reports. **Never commit or copy** the Overleaf ZIP, personal email, historical drafts, or internal reports into this repository.

## Before making the repository public
Verify rights to distribute all included materials, including the original Claude-assisted certificate and any dependencies; verify the license scope in `docs/LICENSE_SCOPE.md` and the two included license texts. Inspect `MANIFEST_SHA256.txt` and the final git diff. Create the repository using a new EMPTY GitHub repository rather than pushing to any existing company/product codebase. Do not claim that the code is a substitute for expert mathematical review.
