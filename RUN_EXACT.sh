#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/ind"
python build_fast_exact.py
python verify_omitted_from_exact_data.py
OPENBLAS_NUM_THREADS=1 python make_certified_integer_data.py
g++ -O3 -std=c++17 certify_inverse_rational.cpp -lgmpxx -lgmp -o certify_inverse_rational
./certify_inverse_rational integer_certificate_data.txt
cd "$ROOT/prior"
python check_rounded_symmetry.py
g++ -O3 -std=c++17 check_rounded_grams.cpp -lgmpxx -lgmp -o check_rounded_grams
./check_rounded_grams rounded_grams.txt
cd "$ROOT"
python check_free_roundoff.py
cd "$ROOT/ind"
python verify_2to3_sympy_independent.py
echo 'UNIVERSAL E3 SOS CERTIFICATE: ALL EXACT CHECKS PASS'
