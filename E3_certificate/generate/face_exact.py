"""Exact ranges of optimal-configuration moment blocks, and exact complement bases B_rho."""
import pickle, time, itertools
from kmat import *
from algebra import GROUP, SIGNS
from configs_exact import ideal_config, config43
from field import K, F, KZ, KO, KI
from flint import fmpq
half=K(fmpq(1,2))
def su2_elements():
    els=[]
    for sgn in (1,-1):
        for idx in range(4):
            q=[0,0,0,0]; q[idx]=sgn; els.append(q)
    for s in itertools.product([1,-1],repeat=4): els.append([fmpq(x,2) for x in s])
    return els
def quat_to_U(q):
    a,b,c,d=[K(fmpq(x)) for x in q]
    # U = a I - i (b X + c Y + d Z)
    U=scal(a,SIG[0])
    for co,S in ((b,SIG[1]),(c,SIG[2]),(d,SIG[3])):
        U=add(U,scal(K(0,-1)*co,S))
    return U
def rotation_of(U):
    # R_kj with U sigma_j U^dag = sum_k R_kj sigma_k
    Ud=dag(U); R=[[0]*3 for _ in range(3)]
    for j in range(3):
        M=mm(mm(U,SIG[j+1]),Ud)
        for k in range(3):
            v=tr(mm(SIG[k+1],M))*K(fmpq(1,2))
            assert v.im.iszero()
            R[k][j]=v.re
    return R
UT={}
for q in su2_elements():
    U=quat_to_U(q); R=rotation_of(U)
    Rint=[[int(float(x)) if True else None for x in row] for row in R]
    if not all(F(Rint[i][j])==R[i][j] for i in range(3) for j in range(3)): continue
    for t,(Rt,pm) in enumerate(GROUP):
        if Rt==Rint and t not in UT: UT[t]=U
assert len(UT)==12
def image(Ps,N,t):
    U=kron(UT[t],eye(N)); Ud=dag(U); pm=GROUP[t][1]
    new=[None]*4
    for c in range(4): new[pm[c]]=mm(mm(U,Ps[c]),Ud)
    return new
def word_op(w,Ps,N):
    M=eye(2*N)
    for typ,x in w:
        M=mm(M,kron(SIG[x],eye(N)) if typ=='P' else Ps[x])
    return M
def range_vectors(Ps,N,L,sec):
    """vectors x_{kl}[u] = conj((w_u P)_{kl}) for all images; returns list of K-vectors of length len(L)"""
    vecs=[]
    for t in range(12):
        P_=image(Ps,N,t)
        if sec=='P': proj=P_[0]
        else:
            Pi=P_[0]
            for c in (1,2,3): Pi=add(Pi,P_[c])
            proj=add(eye(2*N),scal(K(-1),Pi))
        ops=[mm(word_op(w,P_,N),proj) for w in L]
        for k in range(2*N):
            for l in range(2*N):
                v=[ops[u][k][l].conj() for u in range(len(L))]
                if any(not x.iszero() for x in v): vecs.append(v)
    return vecs
if __name__=='__main__':
    D=pickle.load(open('cert_data.pkl','rb'))
    LP,LQ=D['data']['LP'],D['data']['LQ']; SB=D['SB']
    confs=[ideal_config(0),config43(1)]
    rv={}
    for sec,L in (('P',LP),('Q',LQ)):
        t0=time.time(); vs=[]
        for Ps,N in confs: vs+=range_vectors(Ps,N,L,sec)
        rv[sec]=vs; print('sector',sec,'range vectors',len(vs),round(time.time()-t0,1),'s',flush=True)
    B={}; RNG={}
    for name,Qb in SB.items():
        sec=name[0]; n=len(Qb)
        proj=[[sum((Qb[a][u].conj()*x[u] for u in range(len(x)) if not Qb[a][u].iszero()),KZ) for a in range(n)] for x in rv[sec]]
        basis,piv=rref_rows(proj)
        RNG[name]=basis
        # complement: vectors y with <r, y> = sum conj(r_a) y_a = 0 for all r in range -> nullspace of rows conj(r)
        comp=nullspace_of_rows([[x.conj() for x in r] for r in basis],n)
        B[name]=comp
        print(name,'size',n,'exact optimal-range dim',len(basis),'complement dim',len(comp),flush=True)
    pickle.dump(dict(B=B,RNG=RNG),open('face_data.pkl','wb'))
