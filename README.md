# Exact tetrahedral compression certificate — P014 reproducibility

**Public computational supplement for referees** · prepared release candidate 2026-10-23

This repository accompanies the PRA manuscript *Sharp Completeness Thresholds and Reservoir Costs Across Temperature Regimes*. It contains the machine-readable **E3 tracial sum-of-squares (SOS) certificate**, two independently implemented exact checking routes, exact equality checks, and non-proof numerical tests of the thermal implementation.

**Scientific claim certified by the E3 code:** for four contractions \(K_i:\mathbb C^m\to\mathbb C^n\),

\[
 m-\frac14\sum_{i,j=1}^4(P_\triangle)_{ji}\operatorname{Tr}(K_iK_j^\dagger)
 \ge \frac{3-2\sqrt2}{3}(2m-n)_+.
\]

The manuscript appendix supplies the **mathematical reduction** from physical contractions to the finite SOS identity. The code verifies its explicit algebraic data; running it alone is not a substitute for checking the mathematical reduction. Finite-bath numerical tests do not establish the general thermal resource laws by themselves.

## For a referee: quick exact verification (no FLINT)

Prerequisites: Python 3.10+, NumPy, SymPy, C++17 compiler, GMP and GMPXX development libraries (for example `libgmp-dev` on Ubuntu).

```bash
python -m pip install -r requirements.txt
bash RUN_PRECOMPUTED.sh
```

This checks the supplied exact coefficient rows, the three dependent identities (rows 0, 225, 260), rigorous pivot enclosures, exact positive definiteness of the seven rational Gram approximants, the free-variable roundoff, and the exact 2-to-3 saturation. The 353-pivot error is bounded using integer/rational inequalities, not floating eigenvalues.

## Full independent reconstruction (without FLINT)

```bash
bash RUN_EXACT.sh
```

This first reconstructs the 356 SOS coefficient equations from the original canonical-word and moment rules, regenerates the rational integer-check inputs (NumPy is used only to *propose* an approximate inverse), then validates them with exact GMP arithmetic. It can be considerably slower than the precomputed check. In the packaging audit, the reconstruction, rational enclosure, seven-block positivity and fixed-free checks all completed separately; the combined command hit the runtime limit before its final symbolic-saturation step. That symbolic step passed independently under `RUN_PRECOMPUTED.sh`.

## Original certificate verifier (optional, needs python-flint)

```bash
python -m pip install python-flint
cd E3_certificate && python verify_certificate.py
```

Expected success marker: `CERTIFICATE VERIFIED`. To regenerate the candidate certificate from its optimisation pipeline, additional NumPy, SciPy, Clarabel and FLINT dependencies are needed; see `E3_certificate/generate/run_all.sh`.

## Independent physical sanity checks (not universal proofs)

```bash
python physical_checks/check_explicit_energy_unitary.py
python physical_checks/check_taper_and_e2.py
```

The first constructs a finite 40×40 global unitary that conserves total energy and reproduces the binary Born-rule probability. The second tests integer-degeneracy tapered baths and outer-crossover finite examples. These support, but do not replace, the physical arguments in the manuscript.

## Structure

- `E3_certificate/`: original SOS data, Pauli-word algebra, moment rules, FLINT verifier, equality examples, optional generation pipeline.
- `ind/`: independent number-field coefficient reconstruction, exact omitted-row relations, rigorous rational-inverse enclosure and exact saturation.
- `prior/`: precomputed rounded Gram inputs and exact rational LDLᵀ positivity checker.
- `physical_checks/`: small numerical/constructive tests, not proof of global optimality.
- `docs/PROOF_MAP.md`: correspondence with E3 certificate proof obligations and limitations.
- `MANIFEST_SHA256.txt`: per-file release hashes, excluding the manifest itself.

This repository deliberately does **not** contain the unpublished manuscript source/PDF, private research notes, editorial reports, personal contact information, generated binaries, or legacy attempts. There is no DOI yet; cite a tagged release/commit until an archive DOI is genuinely registered.

**Licensing:** the original software and documentation are licensed under **BSD-3-Clause** (`LICENSE`), while specified mathematical certificate datasets are dedicated under **CC0-1.0** (`LICENSE-DATA`). See [`docs/LICENSE_SCOPE.md`](docs/LICENSE_SCOPE.md) for the precise file-level boundary. Reusing CC0 data does not legally require attribution, but citing the scientific work is requested.
