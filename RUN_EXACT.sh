#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
TMP_DIR_E3="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR_E3"' EXIT
export E3_GENERATED_ROWS_PATH="$TMP_DIR_E3/fast_exact_rows.json.gz"
export E3_INTEGER_OUTPUT="$TMP_DIR_E3/integer_certificate_data.txt"
cd "$ROOT/ind"
python build_fast_exact.py
# Compare mathematical payload, not incidental gzip header metadata.
cmp <(gzip -cd "$E3_GENERATED_ROWS_PATH") <(gzip -cd "$ROOT/ind/fast_exact_rows.json.gz")
python verify_omitted_from_exact_data.py
OPENBLAS_NUM_THREADS=1 python make_certified_integer_data.py
cmp "$E3_INTEGER_OUTPUT" "$ROOT/ind/integer_certificate_data.txt"
g++ -O3 -std=c++17 certify_inverse_rational.cpp -lgmpxx -lgmp -o certify_inverse_rational
./certify_inverse_rational "$E3_INTEGER_OUTPUT"
cd "$ROOT/prior"
python check_rounded_symmetry.py
g++ -O3 -std=c++17 check_rounded_grams.cpp -lgmpxx -lgmp -o check_rounded_grams
./check_rounded_grams rounded_grams.txt
cd "$ROOT"
python check_free_roundoff.py
cd "$ROOT/ind"
python verify_2to3_sympy_independent.py
echo 'UNIVERSAL E3 SOS CERTIFICATE: ALL EXACT CHECKS PASS'
