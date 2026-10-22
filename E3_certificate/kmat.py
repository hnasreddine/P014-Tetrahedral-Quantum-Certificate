"""Small dense matrices over K (lists of lists of K)."""
from field import K, F, KZ, KO, KI, R2, R3, R6
from flint import fmpq
def zeros(n,m): return [[KZ for _ in range(m)] for _ in range(n)]
def eye(n): return [[KO if i==j else KZ for j in range(n)] for i in range(n)]
def mm(A,B):
    n=len(A); m=len(B[0]); l=len(B)
    out=zeros(n,m)
    for i in range(n):
        Ai=A[i]
        for k in range(l):
            a=Ai[k]
            if a.iszero(): continue
            Bk=B[k]
            row=out[i]
            for j in range(m):
                b=Bk[j]
                if b.iszero(): continue
                row[j]=row[j]+a*b
    return out
def add(A,B): return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def scal(c,A): return [[c*A[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def dag(A): return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]
def kron(A,B):
    n,m=len(A),len(A[0]); p,q=len(B),len(B[0])
    out=zeros(n*p,m*q)
    for i in range(n):
        for j in range(m):
            a=A[i][j]
            if a.iszero(): continue
            for k in range(p):
                for l in range(q):
                    if not B[k][l].iszero(): out[i*p+k][j*q+l]=a*B[k][l]
    return out
def tr(A): 
    s=KZ
    for i in range(len(A)): s=s+A[i][i]
    return s
def iszero(A): return all(x.iszero() for r in A for x in r)
def tofloat(A):
    import numpy as np
    return np.array([[complex(x) for x in r] for r in A])
SIG=[ [[KO,KZ],[KZ,KO]], [[KZ,KO],[KO,KZ]], [[KZ,K(0,-1)],[KI,KZ]], [[KO,KZ],[KZ,K(-1)]] ]
def rref_rows(vecs):
    """exact row echelon of a list of K-vectors; returns independent basis (list of vectors)"""
    basis=[]; pivots=[]
    for v in vecs:
        v=list(v)
        for (p,b) in zip(pivots,basis):
            if not v[p].iszero():
                c=v[p]
                v=[v[i]-c*b[i] for i in range(len(v))]
        nzi=[i for i,x in enumerate(v) if not x.iszero()]
        if not nzi: continue
        p=nzi[0]; inv=v[p].inv(); v=[x*inv for x in v]
        # eliminate p from previous basis vectors
        for idx in range(len(basis)):
            b=basis[idx]
            if not b[p].iszero():
                c=b[p]; basis[idx]=[b[i]-c*v[i] for i in range(len(v))]
        basis.append(v); pivots.append(p)
    return basis,pivots
def nullspace_of_rows(rows,n):
    """basis of {x in K^n : sum_i r_i x_i = 0 for all rows r} (bilinear, no conjugation)"""
    basis,piv=rref_rows(rows)
    free=[i for i in range(n) if i not in piv]
    out=[]
    for f in free:
        x=[KZ]*n; x[f]=KO
        for (p,b) in zip(piv,basis): x[p]=-b[f]
        out.append(x)
    return out
