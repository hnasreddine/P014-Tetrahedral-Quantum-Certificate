"""Exact optimal configurations (ideal piece and the 4-in-3 configuration) in the Pauli frame."""
from kmat import *
from algebra import SIGNS
from flint import fmpq
half=fmpq(1,2)
def rho_label(c):
    # (I + (s.sigma)/sqrt3)/2 ; 1/sqrt3 = sqrt3/3
    s=SIGNS[c]; inv3=F(0,0,fmpq(1,3),0)
    M=scal(K(half),SIG[0])
    for j in range(3):
        M=add(M,scal(K(F(0,0,fmpq(s[j],6),0)),SIG[j+1]))
    return M
def ideal_config(c=0):
    P=[zeros(2,2) for _ in range(4)]; P[c]=rho_label(c); return P,1
def config43(s=1):
    p=F(fmpq(1,2),fmpq(1,3),0,0)            # (3+2sqrt2)/6
    q=F(fmpq(1,2),fmpq(-1,3),0,0)           # 1-p
    rho=rho_label(0); I2=SIG[0]; rhop=add(I2,scal(K(-1),rho))
    w=K(fmpq(-1,2),F(0))+K(0,F(0,0,fmpq(1,2),0))   # omega
    w2=w*w
    third=K(fmpq(1,3))
    a=[KO,KO,KO]; b=[KO,w,w2] if s==1 else [KO,w2,w]
    aa=[[third*a[i]*a[j].conj() for j in range(3)] for i in range(3)]
    bb=[[third*b[i]*b[j].conj() for j in range(3)] for i in range(3)]
    ab=[[third*a[i]*b[j].conj() for j in range(3)] for i in range(3)]
    ba=[[third*b[i]*a[j].conj() for j in range(3)] for i in range(3)]
    coef=K(F(0,0,0,fmpq(1,12)))             # (1/6)*sqrt(3/2) = sqrt6/12
    X=mm(mm(rhop,SIG[1]),rho)               # rho^perp sigma_x rho
    Xd=mm(mm(rho,SIG[1]),rhop)
    P0=add(add(kron(scal(K(p),rho),aa),kron(scal(K(q),rhop),bb)),add(kron(scal(coef,X),ba),kron(scal(coef,Xd),ab)))
    Rs=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]]
    Ps=[]
    for c in range(4):
        Rc=[[K(Rs[c][i]) if i==j else KZ for j in range(3)] for i in range(3)]
        Uc=kron(SIG[c],Rc)
        Ps.append(mm(mm(Uc,P0),dag(Uc)))
    return Ps,3
def check(Ps,N):
    ok=True
    for c in range(4):
        if not iszero(add(mm(Ps[c],Ps[c]),scal(K(-1),Ps[c]))): ok=False
        for d in range(4):
            if c!=d and not iszero(mm(Ps[c],Ps[d])): ok=False
    return ok
if __name__=='__main__':
    for s in (1,2):
        Ps,N=config43(s); print('s',s,'projections & orthogonal exactly:',check(Ps,N))
        # cost = sum_c Tr[(rho_c^perp x I) Pi_c] ; expect 2g = 2(3-2sqrt2)/3
        cost=KZ
        for c in range(4):
            rp=add(SIG[0],scal(K(-1),rho_label(c)))
            cost=cost+tr(mm(kron(rp,eye(3)),Ps[c]))
        print('   exact cost =',cost, ' (2g = 2 - (4/3)sqrt2)')
    P,N=ideal_config(0); print('ideal projection:',iszero(add(mm(P[0],P[0]),scal(K(-1),P[0]))))
