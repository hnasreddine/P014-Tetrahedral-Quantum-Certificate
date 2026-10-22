"""Numerical sanity check (NOT a proof) of the covariant reduction E3 -> (R):
For isometries K_a : C^m -> C^n (a in {0,1,2,3} ~ Paulis), with Pauli-orbit SIC v_a = sigma_a xi,
   F = (1/8) || sum_a v_a (x) K_a ||_2^2 .
Define L_c = (1/2) sum_a chi_c(a) [(I (x) sigma_a)|Phi>] (x) K_a : C^m -> C^2 (x) C^2 (x) C^n,
chi_c(a) = +1 if sigma_a, sigma_c commute, -1 otherwise.  Claims (proved in the note):
  (i) L_c^dag L_d = delta_cd I_m ;
  (ii) || (eta_c^dag (x) I) L_c ||^2 = F for every c, where eta_c = sigma_c conj(xi).
"""
import numpy as np
rng=np.random.default_rng(0)
s=[np.eye(2),np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1]).astype(complex)]
n_=np.ones(3)/np.sqrt(3); H=(np.eye(2)+sum(n_[k]*s[k+1] for k in range(3)))/2
w,U=np.linalg.eigh(H); xi=U[:,1]
v=[s[a]@xi for a in range(4)]
def comm(a,c): return 1 if (a==0 or c==0 or a==c) else -1
Phi=np.array([1,0,0,1])/np.sqrt(2)
for trial in range(5):
    m=int(rng.integers(1,4)); n=int(rng.integers(m,2*m+1))
    K=[np.linalg.qr(rng.standard_normal((n,m))+1j*rng.standard_normal((n,m)))[0] for _ in range(4)]
    T=sum(np.kron(v[a].reshape(2,1),K[a]) for a in range(4)); F=np.linalg.norm(T)**2/8
    L=[0.5*sum(comm(a,c)*np.kron((np.kron(np.eye(2),s[a])@Phi).reshape(4,1),K[a]) for a in range(4)) for c in range(4)]
    err_iso=max(np.abs(L[c].conj().T@L[d]-(c==d)*np.eye(m)).max() for c in range(4) for d in range(4))
    eta=[s[c]@np.conj(xi) for c in range(4)]
    vals=[np.linalg.norm(np.kron(eta[c].conj().reshape(1,2),np.eye(2*n))@L[c])**2 for c in range(4)]
    print(f'm={m} n={n}: isometry/orthogonality error {err_iso:.1e};  F={F:.12f};  ||(eta_c^* x I)L_c||^2 = {np.round(vals,12)}')
