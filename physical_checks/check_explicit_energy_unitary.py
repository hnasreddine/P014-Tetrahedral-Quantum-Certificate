"""Independent explicit full-system test of a 2->3 tetrahedral bath corner.
System: 8 nondegenerate levels, bath: 2 zero-energy and 3 unit-energy eigenstates.
All computations apart from the reference exact value use NumPy matrices built here.
"""
import numpy as np
from math import sqrt, log
np.set_printoptions(precision=12)
Delta=1.0
El=np.array([0.,3.,10.,24.]); Eu=El+Delta
Es=np.r_[El,Eu]
Eb=np.array([0.,0.,1.,1.,1.]); dimB=5; dimS=8
Etotal=(Es[:,None]+Eb[None,:]).reshape(-1)
D=dimS*dimB
omega=np.exp(2j*np.pi/3)
v=np.vstack([np.array([1.,0.]), *[np.array([1/sqrt(3),sqrt(2/3)*omega**j]) for j in range(3)]])
P=v@v.conj().T/2
I4=np.eye(4,dtype=complex)
S=np.column_stack((I4[:,0],(I4[:,1]+I4[:,2])/sqrt(2),I4[:,3]))
K=[]
for z in v:
    A=S.conj().T@np.kron(z[:,None],np.eye(2))
    R=sqrt(2)*np.eye(2)+(1-sqrt(2))*np.outer(z,z.conj())
    k=A@R
    assert np.allclose(k.conj().T@k,np.eye(2),atol=1e-12)
    K.append(k)
G=np.array([[np.trace(k@l.conj().T) for l in K] for k in K])
F=(np.trace(P@G)/4).real
g=(3-2*sqrt(2))/3
assert abs(F-(2-g))<1e-12,(F,2-g)
U=np.eye(D,dtype=complex)
for i,k in enumerate(K):
    incoming=[(i+4)*dimB+b for b in [0,1]]
    outgoing=[i*dimB+b for b in [2,3,4]]
    indices=incoming+outgoing
    bottom=np.eye(3)-k@k.conj().T
    M=np.block([[np.zeros((2,2)), -k.conj().T], [k,bottom]])
    assert np.allclose(M.conj().T@M,np.eye(5),atol=1e-12)
    U[np.ix_(indices,indices)]=M
Ue=np.max(np.abs(U.conj().T@U-np.eye(D)))
comm=np.max(np.abs((Etotal[:,None]-Etotal[None,:])*U))
psiS=np.r_[np.zeros(4),np.ones(4)/2].astype(complex)
I=np.eye(D)
# Projector on tetrahedral lower band, tensored with bath identity.
PS=np.zeros((8,8),dtype=complex); PS[:4,:4]=P
PM=np.kron(PS,np.eye(dimB))
Z=2+3/2
prob=0
for b in range(dimB):
    psiB=np.zeros(dimB);psiB[b]=1
    psi=np.kron(psiS,psiB)
    out=U@psi
    w=(1 if b<2 else 0.5)/Z
    prob+=w*np.vdot(out,PM@out).real
fail=1-prob
pred=(2-F+1.5)/Z
pred2=(g+1.5)/Z
print('U unitarity max error',Ue)
print('[U,Htotal] max error',comm)
print('tetrahedral score F(2,3)',F)
print('direct Born failure',fail)
print('energy shell formula',pred)
print('formula using sharp g',pred2)
assert Ue<1e-12 and comm<1e-12
assert abs(fail-pred)<1e-12 and abs(fail-pred2)<1e-12
print('PASS: full energy-conserving 40x40 unitary and binary Born probability')
