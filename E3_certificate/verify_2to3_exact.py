"""Exact symbolic verification (SymPy) of the symmetric-square 2->3 block.
Checks: K_i isometries C^2 -> Sym^2(C^2) ~ C^3, and F = (1/4) sum_ij P_ji Tr(K_i K_j^dag) = (3+2 sqrt2)/3 exactly,
for the tetrahedral frame v_1=(1,0), v_{k+2}=(1/sqrt3, sqrt(2/3) w^k), P_ij = <v_i,v_j>/2 (antilinear in first slot).
Also checks the conjugate convention P -> conj(P) gives the same value."""
import sympy as sp
r2,r3=sp.sqrt(2),sp.sqrt(3)
w=sp.Rational(-1,2)+sp.I*r3/2
v=[sp.Matrix([1,0])]+[sp.Matrix([1/r3, sp.sqrt(sp.Rational(2,3))*w**k]) for k in range(3)]
def ip(a,b): return sp.simplify((a.H*b)[0])
P=sp.Matrix(4,4,lambda i,j: ip(v[i],v[j])/2)
assert sp.simplify(P*P-P)==sp.zeros(4,4) and sp.simplify(P.trace())==2
# inclusion S: Sym^2 -> C^2 (x) C^2 with basis |00>, (|01>+|10>)/sqrt2, |11>
S=sp.Matrix([[1,0,0],[0,1/r2,0],[0,1/r2,0],[0,0,1]])
K=[]
for z in v:
    A=S.H*sp.kronecker_product(z,sp.eye(2))          # 3x2
    AA=sp.simplify(A.H*A)
    assert sp.simplify(AA-(sp.eye(2)+z*z.H)/2)==sp.zeros(2,2)
    Rm=r2*sp.eye(2)+(1-r2)*z*z.H                       # (A^dag A)^(-1/2)
    assert sp.simplify(Rm*Rm*AA-sp.eye(2))==sp.zeros(2,2)
    Ki=sp.simplify(A*Rm)
    assert sp.simplify(Ki.H*Ki-sp.eye(2))==sp.zeros(2,2)
    K.append(Ki)
def score(PP,KK): return sp.simplify(sum(PP[j,i]*(KK[i]*KK[j].H).trace() for i in range(4) for j in range(4))/4)
F_conjP=score(P.conjugate(),K)          # convention P_ij=<v_j,v_i>/2 with K from v
F_P=score(P,[Ki.conjugate() for Ki in K])  # convention P_ij=<v_i,v_j>/2 with conjugated block
F_mismatch=score(P,K)                   # mismatched convention (not optimal)
print('score, convention P_ij=<v_j,v_i>/2, K from v       :',F_conjP)
print('score, convention P_ij=<v_i,v_j>/2, K from conj(v) :',F_P)
print('mismatched pairing (for information)               :',F_mismatch)
F=F_conjP
assert sp.simplify(F_conjP-(3+2*r2)/3)==0 and sp.simplify(F_P-(3+2*r2)/3)==0
g=(3-2*r2)/3
print('cost 2-F =',sp.simplify(2-F),'= g:',sp.simplify(2-F-g)==0)
print('sum_i K_i K_i^dag =',sp.simplify(sum((Ki*Ki.H for Ki in K),sp.zeros(3,3))))
