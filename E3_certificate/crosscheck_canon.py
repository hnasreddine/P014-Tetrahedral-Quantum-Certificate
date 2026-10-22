#!/usr/bin/env python3
"""Numerical cross-check of the exact canonical-form routine: for 20000 random words in sigma_j and Q_c,
tau(word) computed from explicit matrices (random orthogonal projections, N=3) must equal i^k tau(canonical word),
and words declared zero must have zero trace."""
import numpy as np, random
from algebra import canon, keytoks
random.seed(1); rng=np.random.default_rng(0)
letters=[('P',1),('P',2),('P',3)]+[('Q',c) for c in range(4)]
sig=[np.eye(2),np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1]).astype(complex)]
N=3; U=np.linalg.qr(rng.standard_normal((6,6))+1j*rng.standard_normal((6,6)))[0]
Pis=[U[:,0:2]@U[:,0:2].conj().T,U[:,2:3]@U[:,2:3].conj().T,U[:,3:5]@U[:,3:5].conj().T,U[:,5:6]@U[:,5:6].conj().T]
def op(toks):
    M=np.eye(6,dtype=complex)
    for t,x in toks: M=M@(np.kron(sig[x],np.eye(N)) if t=='P' else Pis[x])
    return M
bad=0
for trial in range(20000):
    w=[random.choice(letters) for _ in range(random.randint(0,9))]
    true=np.trace(op(w))/6; r=canon(0,w)
    val=0 if r is None else (1j**r[0])*(1 if r[1]=='I' else np.trace(op(keytoks(r[1])))/6)
    bad+= abs(val-true)>1e-10
print('mismatches:',bad,'of 20000')
