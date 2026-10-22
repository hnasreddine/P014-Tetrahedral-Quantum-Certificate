"""Exact construction of the sector moment data for the certificate.
Sector P (projection Q_0):  M_P[u,v] = tau(u^* v Q_0),  words L_P
Sector Q (complement 1-Pi): M_Q[u,v] = tau(u^* v) - sum_c tau(u^* v Q_c),  words L_Q
Target identity (on A4-symmetric tracial moment functionals):
   g*tau(1) - sum_c tau(D_c Q_c) = <Z_P, M_P> + <Z_Q, M_Q>,
   D_c = (4g-1)/2 * 1 + (1/2) r_c.sigma,  r_c = SIGNS[c]/sqrt3,  g=(3-2sqrt2)/3.
"""
import pickle, time
from algebra import *
from moments import Moments, IK
from field import K, F, KZ, KO, KI, kgauss
from kmat import *
from configs_exact import ideal_config, config43
from flint import fmpq
def words(which):
    L=[()]
    if which=='Q': L+=[(('P',j),) for j in (1,2,3)]
    L+=[(('P',a),('Q',c),('P',b)) for c in range(4) for a in (1,2,3) for b in (1,2,3) if a!=b]
    L+=[(('Q',c),('P',a),('Q',d),('P',b)) for c in range(4) for d in range(4) if c!=d for a in (1,2,3) for b in (1,2,3) if a!=b]
    return L
def sector_terms(u,v,sec):
    k,ua=adjoint(0,list(u)); base=ua+list(v)
    if sec=='P': return [(K(1),base+[('Q',0)])]
    return [(K(1),base)]+[(K(-1),base+[('Q',c)]) for c in range(4)]
G_SEC={'P':stabilizer0(),'Q':list(range(12))}
def irreps(sec):
    from field import OMEGA
    w=OMEGA; w2=w*w; pw=[KO,w,w2]
    if sec=='P':
        return [('Z0',1,lambda t:[[KO]]),('Z1',1,lambda t:[[pw[z3exp(t)]]]),('Z2',1,lambda t:[[pw[(2*z3exp(t))%3]]])]
    return [('A0',1,lambda t:[[KO]]),('A1',1,lambda t:[[pw[z3exp(t)]]]),('A2',1,lambda t:[[pw[(2*z3exp(t))%3]]]),
            ('T',3,lambda t:[[K(GROUP[t][0][i][j]) for j in range(3)] for i in range(3)])]
def build():
    t0=time.time()
    LP,LQ=words('P'),words('Q')
    mom=Moments()
    raw={}
    keys=set()
    for sec,L in (('P',LP),('Q',LQ)):
        for i,u in enumerate(L):
            for j,v in enumerate(L):
                lst=[]
                for coef,toks in sector_terms(u,v,sec):
                    r=canon(0,toks)
                    if r is None: continue
                    k,key=r
                    lst.append((coef*IK[k],key))
                    if key!='I': keys.add(key)
                raw[(sec,i,j)]=lst
    # objective sum_c tau(D_c Q_c)
    g=F(1,fmpq(-2,3),0,0)
    objraw=[]
    for c in range(4):
        r=canon(0,[('Q',c)]); objraw.append((K(F(fmpq(3,2),fmpq(-4,3),0,0)) ,r))   # (4g-1)/2 = (9-8sqrt2)/6 = 3/2 - (4/3)sqrt2
        for j in (1,2,3):
            r=canon(0,[('P',j),('Q',c)])
            objraw.append((K(F(0,0,fmpq(SIGNS[c][j-1],6),0)),r))                # s_cj/(2 sqrt3) = s_cj sqrt3/6
    objraw=[(co*IK[r[0]],r[1]) for co,r in objraw if r is not None]
    for co,key in objraw:
        if key!='I': keys.add(key)
    mom.close(keys); nz=mom.parametrize()
    def lin(lst):
        d={}
        for co,key in lst:
            if key=='I': d['c']=d.get('c',KZ)+co; continue
            for zi,b in mom.expr(key): d[zi]=d.get(zi,KZ)+co*b
        return {k:v for k,v in d.items() if not v.iszero()}
    C={k:lin(v) for k,v in raw.items()}
    OBJ=lin(objraw)
    print('words',len(LP),len(LQ),'keys',len(mom.par),'real moment params',nz,'build',round(time.time()-t0,1),'s',flush=True)
    return dict(LP=LP,LQ=LQ,C=C,OBJ=OBJ,nz=nz,mom=mom,g=g)
def sym_bases(data):
    """exact symmetry-adapted 'first copy' bases Q_rho (columns as K-vectors) for each sector"""
    out={}
    for sec,L in (('P',data['LP']),('Q',data['LQ'])):
        widx={w:i for i,w in enumerate(L)}; n=len(L); grp=G_SEC[sec]
        act_t={}
        for t in grp:
            act_t[t]=[]
            for w in L:
                s,w2=act(t,list(w)); act_t[t].append((widx[tuple(w2)],s))
        for name,d,D in irreps(sec):
            cols=[]
            seen=set()
            for i0 in range(n):
                if i0 in seen: continue
                orb=set(act_t[t][i0][0] for t in grp); seen|=orb
                vecs=[]
                for i in sorted(orb):
                    v=[KZ]*n
                    for t in grp:
                        j,s=act_t[t][i]
                        co=D(t)[0][0].conj()*K(fmpq(d*s,len(grp)))
                        v[j]=v[j]+co
                    vecs.append(v)
                basis,_=rref_rows(vecs)
                cols+=basis
            out[sec+name]=cols
    return out
def block_coeffs(data,Qb,sec):
    """C_rho[a][b] = dict over z (and 'c') of K coefficients, for basis columns Qb (list of K-vectors)"""
    n=len(Qb[0]); d=len(Qb)
    supp=[[(i,q[i]) for i in range(n) if not q[i].iszero()] for q in Qb]
    Cr=[[{} for _ in range(d)] for _ in range(d)]
    C=data['C']
    for a in range(d):
        for b in range(d):
            acc={}
            for (i,qa) in supp[a]:
                qac=qa.conj()
                for (j,qb) in supp[b]:
                    e=C.get((sec,i,j))
                    if not e: continue
                    f=qac*qb
                    for k,v in e.items(): acc[k]=acc.get(k,KZ)+f*v
            Cr[a][b]={k:v for k,v in acc.items() if not v.iszero()}
    return Cr
if __name__=='__main__':
    data=build()
    SB=sym_bases(data)
    print({k:len(v) for k,v in SB.items()})
    blocks={}
    for name,Qb in SB.items():
        t0=time.time(); blocks[name]=block_coeffs(data,Qb,name[0]); print(name,len(Qb),'coeffs',round(time.time()-t0,1),'s',flush=True)
    pickle.dump(dict(data={k:v for k,v in data.items() if k!='mom'},SB=SB,blocks=blocks),open('cert_data.pkl','wb'))
