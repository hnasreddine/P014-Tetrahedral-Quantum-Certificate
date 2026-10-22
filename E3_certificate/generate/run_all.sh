#!/bin/sh
# Regenerates the certificate from scratch (requires numpy, scipy, clarabel, python-flint).
set -e
export PYTHONPATH=..
python3 build_cert.py      # exact moment data, symmetry-adapted bases   -> cert_data.pkl
python3 face_exact.py      # exact optimal-configuration ranges, face bases -> face_data.pkl
python3 reduce_blocks.py   # exact reduced block coefficients              -> reduced.pkl
python3 numsolve.py        # interior-point (Clarabel) solution on the face  -> xnum.pkl
python3 exactify.py        # exact rational/algebraic correction (FLINT)     -> xexact.pkl
python3 export_cert.py     # portable certificate file                       -> certificate_data_regenerated.json.gz
