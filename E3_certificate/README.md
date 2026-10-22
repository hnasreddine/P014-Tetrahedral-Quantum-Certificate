# E3 / P014 — exact certificate for mu_triangle = (3 - 2 sqrt2)/3

Requirements: Python 3, numpy, scipy, sympy, python-flint (>= 0.6); clarabel only for regeneration.

## Proof-critical check (exact arithmetic)
    python3 verify_certificate.py        # ~1 min; must print "CERTIFICATE VERIFIED"
Recomputes from scratch: word algebra and canonical cyclic forms (algebra.py), the A4-symmetric/adjoint
moment relations and their exact real parametrization (moments.py), the two sector moment matrices, and the
identity (eq. sos of the note) coefficient by coefficient over Q(sqrt2, sqrt3); re-solves the 353 pivot Gram
unknowns exactly (FLINT); proves positive definiteness of all seven Gram blocks (exact rational LDL^T +
ball-arithmetic perturbation bound + Weyl).  Data: certificate_data.json.gz (exact strings).

## Supporting exact checks
    python3 verify_2to3_exact.py         # SymPy: 2->3 block isometric, F(2,3) >= (3+2sqrt2)/3, cost = g
    python3 configs_exact.py             # exact 4-in-C^2(x)C^3 optimal labelled configuration, cost = 2g

## Supplementary numerical checks (NOT part of the proof)
    python3 crosscheck_canon.py          # canonical forms vs explicit matrix traces (20000 random words)
    python3 sanity_random_configs.py     # both sides of the SOS identity on explicit symmetrized configurations
    python3 check_reduction_R.py         # covariant amplification map (Proposition 3.2) on random isometries

## Regeneration (optional)
    cd generate && sh run_all.sh         # exact data -> exact faces -> Clarabel interior point -> exact FLINT correction
