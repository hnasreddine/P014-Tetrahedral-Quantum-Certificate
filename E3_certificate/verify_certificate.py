#!/usr/bin/env python3
"""
Exact verification of the sector sum-of-squares certificate for the tracial inequality (R*):

  For mutually orthogonal projections Pi_0,...,Pi_3 on C^2 (x) C^N (any N), tau = Tr/(2N),
        sum_c tau(D_c Pi_c) <= g = (3 - 2 sqrt2)/3,
  D_c = (4g-1)/2 * 1 + (1/2) r_c . sigma,   r_c = SIGNS[c]/sqrt3 (regular tetrahedron).

Certificate (identity on all A4-symmetric tracial moment functionals y):
  g*tau(1) - sum_c tau(D_c Pi_c) = sum_{rho in P-blocks} Tr( Y_rho W_rho^* M_P(y) W_rho )
                                  + sum_{rho in Q-blocks} Tr( Y_rho W_rho^* M_Q(y) W_rho ),
  M_P(y)[u,v] = tau(u^* v Pi_0)            (u,v in word list L_P),
  M_Q(y)[u,v] = tau(u^* v (1 - Pi))        (u,v in word list L_Q),
with Y_rho Hermitian positive definite.  M_P, M_Q are Gram matrices for every actual configuration,
hence PSD, which gives the inequality.

This script recomputes everything from scratch with exact arithmetic in Q(sqrt2, sqrt3, i):
  1. word algebra & canonical cyclic forms, A4-symmetric/adjoint relations among moments;
  2. the exact linear identity (coefficient matching for every free real moment parameter);
  3. the exact Gram unknowns: stored rational 'free' values + exact re-solution of the 'pivot' values;
  4. rigorous positive definiteness of every Y_rho (exact rational LDL^T + ball-arithmetic perturbation bound).
"""
import json, gzip, sys, time
import flint
from flint import fmpq, fmpq_mat, arb, ctx
from field import F, K, KZ, KO
from algebra import canon, adjoint, SIGNS
from moments import Moments, IK
t_start=time.time()
def fparse(l): return F(*[fmpq(*map(int,s.split('/'))) if '/' in s else fmpq(int(s)) for s in l])
def kparse(l): return K(fparse(l[:4]),fparse(l[4:]))
cert=json.load(gzip.open('certificate_data.json.gz','rt'))
ORDER=cert['blocks']
# ---------------- 1. words and exact moment structure ----------------
def words(which):
    L=[()]
    if which=='Q': L+=[(('P',j),) for j in (1,2,3)]
    L+=[(('P',a),('Q',c),('P',b)) for c in range(4) for a in (1,2,3) for b in (1,2,3) if a!=b]
    L+=[(('Q',c),('P',a),('Q',d),('P',b)) for c in range(4) for d in range(4) if c!=d for a in (1,2,3) for b in (1,2,3) if a!=b]
    return L
L={'P':words('P'),'Q':words('Q')}
def sector_terms(u,v,sec):
    _,ua=adjoint(0,list(u)); base=ua+list(v)
    if sec=='P': return [(KO,base+[('Q',0)])]
    return [(KO,base)]+[(K(-1),base+[('Q',c)]) for c in range(4)]
raw={}; keys=set()
for sec in ('P','Q'):
    for i,u in enumerate(L[sec]):
        for j,v in enumerate(L[sec]):
            lst=[]
            for co,toks in sector_terms(u,v,sec):
                r=canon(0,toks)
                if r is None: continue
                k,key=r; lst.append((co*IK[k],key))
                if key!='I': keys.add(key)
            raw[(sec,i,j)]=lst
g=F(1,fmpq(-2,3),0,0)                         # (3-2sqrt2)/3
objraw=[]
for c in range(4):
    r=canon(0,[('Q',c)]); objraw.append((K(F(fmpq(3,2),fmpq(-4,3),0,0)),r))          # (4g-1)/2
    for j in (1,2,3):
        r=canon(0,[('P',j),('Q',c)]); objraw.append((K(F(0,0,fmpq(SIGNS[c][j-1],6),0)),r))   # r_cj/2
objraw=[(co*IK[r[0]],r[1]) for co,r in objraw if r is not None]
for _,key in objraw:
    if key!='I': keys.add(key)
mom=Moments(); mom.close(keys); nz=mom.parametrize()
def lin(lst):
    d={}
    for co,key in lst:
        if key=='I': d['c']=d.get('c',KZ)+co; continue
        for zi,b in mom.expr(key): d[zi]=d.get(zi,KZ)+co*b
    return {k:v for k,v in d.items() if not v.iszero()}
C={k:lin(v) for k,v in raw.items()}; OBJ=lin(objraw)
assert all(v.im.iszero() for v in OBJ.values())
print(f'[1] words |L_P|={len(L["P"])} |L_Q|={len(L["Q"])}; moment classes {len(mom.par)}; free real moment parameters {nz}')
# ---------------- 2. reduced blocks W^* M W ----------------
W={n:[[kparse(e) for e in col] for col in cert['W'][n]] for n in ORDER}
def reduce(sec,Wc):
    n=len(L[sec]); r=len(Wc)
    T=[[{} for _ in range(r)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            e=C.get((sec,i,j))
            if not e: continue
            for b in range(r):
                x=Wc[b][j]
                if x.iszero(): continue
                acc=T[i][b]
                for k,v in e.items(): acc[k]=acc.get(k,KZ)+v*x
    Cr=[[{} for _ in range(r)] for _ in range(r)]
    for a in range(r):
        for i in range(n):
            y=Wc[a][i]
            if y.iszero(): continue
            yc=y.conj()
            for b in range(r):
                for k,v in T[i][b].items(): Cr[a][b][k]=Cr[a][b].get(k,KZ)+yc*v
    for a in range(r):
        for b in range(r): Cr[a][b]={k:v for k,v in Cr[a][b].items() if not v.iszero()}
    return Cr
red={n:reduce(n[0],W[n]) for n in ORDER}
for n in ORDER:
    Cr=red[n]; r=len(Cr)
    for a in range(r):
        for b in range(r):
            for k in set(Cr[a][b])|set(Cr[b][a]):
                assert (Cr[a][b].get(k,KZ)-Cr[b][a].get(k,KZ).conj()).iszero(), 'non-Hermitian block'
print('[2] reduced blocks:',{n:len(red[n]) for n in ORDER},' (exactly Hermitian)')
# ---------------- 3. exact linear identity ----------------
lay={}; nu=0
for n in ORDER:
    r=len(red[n]); idx={}
    for a in range(r): idx[('d',a)]=nu; nu+=1
    for a in range(r):
        for b in range(a+1,r): idx[('re',a,b)]=nu; nu+=1; idx[('im',a,b)]=nu; nu+=1
    lay[n]=(r,idx)
rows=[dict() for _ in range(nz+1)]
rid=lambda k: 0 if k=='c' else 1+k
for n in ORDER:
    Cr=red[n]; r,idx=lay[n]
    for a in range(r):
        for k,v in Cr[a][a].items():
            R=rows[rid(k)]; R[idx[('d',a)]]=R.get(idx[('d',a)],F(0))+v.re
        for b in range(a+1,r):
            for k,v in Cr[b][a].items():
                R=rows[rid(k)]; i1=idx[('re',a,b)]; i2=idx[('im',a,b)]
                R[i1]=R.get(i1,F(0))+v.re*2; R[i2]=R.get(i2,F(0))-v.im*2
rhs=[F(0)]*(nz+1); rhs[0]=g-OBJ.get('c',KZ).re
for k,v in OBJ.items():
    if k!='c': rhs[1+k]=rhs[1+k]-v.re
piv=cert['pivots']; pivset=set(piv); x=[None]*nu
for j,s in cert['free'].items(): x[int(j)]=F(fmpq(*map(int,s.split('/'))))
assert all(x[j] is None for j in piv) and sum(1 for v in x if v is None)==len(piv)
keep=cert['rows_used']; npv=len(piv); pidx={j:t for t,j in enumerate(piv)}
Mx=fmpq_mat(4*npv,4*npv); bx=fmpq_mat(4*npv,1)
for rr,i in enumerate(keep):
    acc=rhs[i]
    for j,v in rows[i].items():
        if j in pivset:
            c=pidx[j]; a0,a1,a2,a3=v.a,v.b,v.c,v.d
            blk=[[a0,2*a1,3*a2,6*a3],[a1,a0,3*a3,3*a2],[a2,2*a3,a0,2*a1],[a3,a2,a1,a0]]
            for p in range(4):
                for q in range(4):
                    if blk[p][q]!=0: Mx[4*rr+p,4*c+q]+=blk[p][q]
        else: acc=acc-v*x[j]
    for p,val in enumerate((acc.a,acc.b,acc.c,acc.d)): bx[4*rr+p,0]=val
sol=Mx.solve(bx)
for t,j in enumerate(piv): x[j]=F(sol[4*t,0],sol[4*t+1,0],sol[4*t+2,0],sol[4*t+3,0])
viol=0
for i,R in enumerate(rows):
    acc=F(0)
    for j,v in R.items(): acc=acc+v*x[j]
    if not (acc-rhs[i]).iszero(): viol+=1
print(f'[3] exact identity: {len(rows)} coefficient equations over Q(sqrt2,sqrt3), violated: {viol}')
# ---------------- 4. rigorous positive definiteness ----------------
DEN=10**14
def rr_(v): return fmpq(int(round(float(v)*DEN)),DEN)
def embed(n,val,zero):
    r,idx=lay[n]; m=2*r; E=[[zero for _ in range(m)] for _ in range(m)]
    for a in range(r):
        d=val(x[idx[('d',a)]]); E[a][a]=d; E[a+r][a+r]=d
        for b in range(a+1,r):
            re=val(x[idx[('re',a,b)]]); im=val(x[idx[('im',a,b)]])
            E[a][b]=re; E[b][a]=re; E[a+r][b+r]=re; E[b+r][a+r]=re
            E[a+r][b]=im; E[b][a+r]=im; E[b+r][a]=-im; E[a][b+r]=-im
    return E
def ldl_pd(E):
    n=len(E); A=[row[:] for row in E]
    for k in range(n):
        p=A[k][k]
        if p<=0: return False
        for i in range(k+1,n):
            if A[i][k]==0: continue
            f=A[i][k]/p
            for j in range(k+1,n):
                if A[k][j]!=0: A[i][j]-=f*A[k][j]
    return True
ctx.prec=300; allpd=True
for n in ORDER:
    Er=embed(n,rr_,fmpq(0))
    Ed=embed(n,lambda v: v.arb(300)-arb(rr_(v)),arb(0))
    fro=arb(0)
    for row in Ed:
        for e in row: fro+=e*e
    fro=fro.sqrt(); ub=fro.mid()+fro.rad()
    mu=fmpq(int(float(ub)*1e18)+2,10**18)*2          # mu > 2*||Delta||_F (rigorous upper bound)
    assert arb(mu)>2*fro
    ok=ldl_pd([[Er[i][j]-(mu if i==j else 0) for j in range(len(Er))] for i in range(len(Er))])
    allpd&=ok
    print(f'    block {n}: complex size {len(Er)//2}: E(Y_rat) - mu*I positive definite (exact LDL^T): {ok};  ||E(Y)-E(Y_rat)||_F <= {float(ub):.2e} < mu = {float(mu):.2e}')
print(f'[4] all Gram blocks positive definite: {allpd}')
ok_all=(viol==0) and allpd
# export floating-point copies of the exact Gram blocks for the independent numerical sanity check
import numpy as np
Yf={}
for n in ORDER:
    r,idx=lay[n]; Y=np.zeros((r,r),complex)
    for a in range(r):
        Y[a,a]=float(x[idx[('d',a)]])
        for b in range(a+1,r):
            Y[a,b]=float(x[idx[('re',a,b)]])+1j*float(x[idx[('im',a,b)]]); Y[b,a]=np.conj(Y[a,b])
    Yf[n]=Y
np.savez('Y_float.npz',**Yf)
print('CERTIFICATE VERIFIED' if ok_all else 'VERIFICATION FAILED', f'(total {time.time()-t_start:.0f}s)')
sys.exit(0 if ok_all else 1)
