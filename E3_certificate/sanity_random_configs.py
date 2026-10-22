#!/usr/bin/env python3
"""Independent NUMERICAL sanity check (not part of the proof): evaluates both sides of the certificate identity
   g - sum_c tau(D_c Pi_c) = sum_rho Tr(Y_rho W_rho^* M_sector W_rho)
on explicit A4-symmetrized configurations (random ones and the two optimal ones), using matrices only.
Run verify_certificate.py first (it writes Y_float.npz)."""
import json, gzip, numpy as np
from algebra import GROUP, SIGNS
from field import K
from configs_exact import config43, ideal_config
cert=json.load(gzip.open('certificate_data.json.gz','rt')); ORDER=cert['blocks']
Y=dict(np.load('Y_float.npz'))
def fl(l):
    v=[float(eval(s)) if '/' in s else float(int(s)) for s in l]
    return v[0]+v[1]*2**0.5+v[2]*3**0.5+v[3]*6**0.5
W={n:np.array([[fl(e[:4])+1j*fl(e[4:]) for e in col] for col in cert['W'][n]]).T for n in ORDER}
def words(which):
    L=[()]
    if which=='Q': L+=[(('P',j),) for j in (1,2,3)]
    L+=[(('P',a),('Q',c),('P',b)) for c in range(4) for a in (1,2,3) for b in (1,2,3) if a!=b]
    L+=[(('Q',c),('P',a),('Q',d),('P',b)) for c in range(4) for d in range(4) if c!=d for a in (1,2,3) for b in (1,2,3) if a!=b]
    return L
L={'P':words('P'),'Q':words('Q')}
g=(3-2*np.sqrt(2))/3
sig=[np.eye(2),np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1]).astype(complex)]
def U_of(R):
    # qubit unitary implementing the rotation R (adjoint action on Pauli vector); search signed products
    import itertools
    cands=[np.eye(2)]+[ -1j*s for s in sig[1:]]
    for s in itertools.product([1,-1],repeat=3):
        cands.append((np.eye(2)-1j*(s[0]*sig[1]+s[1]*sig[2]+s[2]*sig[3]))/2)
        cands.append((np.eye(2)+1j*(s[0]*sig[1]+s[1]*sig[2]+s[2]*sig[3]))/2)
    for U in cands:
        ok=all(np.allclose(U@sig[j+1]@U.conj().T, sum(R[k][j]*sig[k+1] for k in range(3))) for j in range(3))
        if ok: return U
    raise ValueError
UT=[U_of(R) for R,pm in GROUP]
def symmetrize(Pis,N):
    out=[np.zeros((2*N*12,2*N*12),complex) for _ in range(4)]
    for t,(R,pm) in enumerate(GROUP):
        U=np.kron(UT[t],np.eye(N)); E=np.zeros((12,12)); E[t,t]=1
        for c in range(4):
            blk=(U@Pis[c]@U.conj().T).reshape(2,N,2,N)
            out[pm[c]]+=np.einsum('aibj,st->aisbjt',blk,E).reshape(2*N*12,2*N*12)
    return out,12*N
def word_op(w,P,M):
    O=np.eye(2*M,dtype=complex)
    for typ,v in w: O=O@(np.kron(sig[v],np.eye(M)) if typ=='P' else P[v])
    return O
def both_sides(Pis,N):
    S,M=symmetrize(Pis,N); Pi=sum(S); tau=lambda A: np.trace(A)/(2*M)
    lhs=g-sum(tau(np.kron((4*g-1)/2*np.eye(2)+0.5*sum(SIGNS[c][j]/np.sqrt(3)*sig[j+1] for j in range(3)),np.eye(M))@S[c]) for c in range(4))
    rhs=0
    for sec,proj in (('P',S[0]),('Q',np.eye(2*M)-Pi)):
        V=np.array([(word_op(w,S,M)@proj).reshape(-1) for w in L[sec]])
        Ms=(V.conj()@V.T)/(2*M)
        for n in ORDER:
            if n[0]==sec: rhs+=np.trace(Y[n]@W[n].conj().T@Ms@W[n])
    return lhs.real,rhs
rng=np.random.default_rng(2026)
worst=0
for trial in range(8):
    N=int(rng.integers(1,4)); ranks=list(rng.multinomial(int(rng.integers(1,2*N+1)),[0.25]*4))
    U=np.linalg.qr(rng.standard_normal((2*N,2*N))+1j*rng.standard_normal((2*N,2*N)))[0]
    Pis=[];o=0
    for d in ranks: Pis.append(U[:,o:o+d]@U[:,o:o+d].conj().T); o+=d
    lhs,rhs=both_sides(Pis,N); worst=max(worst,abs(lhs-rhs))
    print(f'random N={N} ranks={[int(r) for r in ranks]}:  g - obj = {lhs:+.13f}   SOS value = {rhs.real:+.13f}')
def tofl(Ps): return [np.array([[complex(x) for x in row] for row in P]) for P in Ps]
P,N=config43(1); lhs,rhs=both_sides(tofl(P),N); print(f'optimal 4-in-3 configuration: g - obj = {lhs:+.2e}, SOS value = {rhs.real:+.2e}')
P,N=ideal_config(0); lhs,rhs=both_sides(tofl(P),N); print(f'optimal ideal piece:          g - obj = {lhs:+.2e}, SOS value = {rhs.real:+.2e}')
print('max |difference| on random configurations:',f'{worst:.1e}')
