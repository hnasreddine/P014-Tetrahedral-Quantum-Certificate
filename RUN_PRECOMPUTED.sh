#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/ind"
python verify_omitted_from_exact_data.py
g++ -O2 -std=c++17 certify_inverse_rational.cpp -lgmpxx -lgmp -o certify_inverse_rational
./certify_inverse_rational integer_certificate_data.txt
cd "$ROOT/prior"
python check_rounded_symmetry.py
g++ -O2 -std=c++17 check_rounded_grams.cpp -lgmpxx -lgmp -o check_rounded_grams
./check_rounded_grams rounded_grams.txt
cd "$ROOT"
python check_free_roundoff.py
python ind/verify_2to3_sympy_independent.py
echo 'E3 PRECOMPUTED EXACT VERIFICATION: PASS'
