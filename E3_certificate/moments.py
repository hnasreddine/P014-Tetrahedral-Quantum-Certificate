"""Exact parametrization of A4-symmetric tracial moment functionals (union-find with Gaussian-unit phases)."""
from algebra import *
from field import K, F, kgauss, KZ, KO, KI
from flint import fmpq
IK=[K(1),K(0,1),K(-1),K(0,-1)]   # i^k
class Moments:
    def __init__(s):
        s.par={}; s.tr={}; s.self={}; s.zero=set()
    def add(s,c):
        if c not in s.par: s.par[c]=c; s.tr[c]=(0,0); s.self[c]=[]; return True
        return False
    def find(s,c):
        # y_c = i^k * conj^f(y_root)
        if s.par[c]==c: return c,(0,0)
        r,(k2,f2)=s.find(s.par[c]); k1,f1=s.tr[c]
        # y_c = i^k1 conj^f1( i^k2 conj^f2 y_r ) = i^(k1 + (-1)^f1 k2) conj^(f1+f2) y_r
        k=(k1+(k2 if f1==0 else -k2))%4; f=(f1+f2)%2
        s.par[c]=r; s.tr[c]=(k,f); return r,(k,f)
    def relate(s,a,b,k,f):
        """impose y_a = i^k conj^f(y_b)"""
        ra,(ka,fa)=s.find(a); rb,(kb,fb)=s.find(b)
        # i^ka conj^fa(y_ra) = i^k conj^f( i^kb conj^fb(y_rb) ) = i^(k + (-1)^f kb) conj^(f+fb)(y_rb)
        q=(k+(kb if f==0 else -kb))%4; ff=(f+fb)%2
        # conj^fa(y_ra) = i^(q-ka) conj^ff y_rb  => y_ra = i^((-1)^fa (q-ka)) conj^(ff+fa) y_rb
        kk=((q-ka) if fa==0 else -(q-ka))%4; f2=(ff+fa)%2
        if ra==rb: s.self[ra].append((kk,f2)); return
        s.par[ra]=rb; s.tr[ra]=(kk,f2)
        for (k3,f3) in s.self[ra]:
            # y_ra = i^k3 conj^f3 y_ra with y_ra = i^kk conj^f2 y_rb
            # i^kk conj^f2 y_rb = i^k3 conj^f3( i^kk conj^f2 y_rb ) = i^(k3 + (-1)^f3 kk) conj^(f3+f2) y_rb
            # => conj^f2 y_rb = i^(k3 + (-1)^f3 kk - kk) conj^(f3+f2) y_rb ; apply conj^f2
            e=(k3+(kk if f3==0 else -kk)-kk)%4
            e=(e if f2==0 else -e)%4
            s.self[rb].append((e,f3))
        s.self[ra]=[]
    def close(s,keys):
        """add keys, close under adjoint and the A4 action, record relations"""
        fr=[c for c in keys if s.add(c)]
        while fr:
            c=fr.pop()
            toks=keytoks(c)
            r=canon(*adjoint(0,toks))
            if r is None: s.zero.add(c)
            else:
                k2,c2=r
                assert c2!='I'
                if s.add(c2): fr.append(c2)
                # tau(c^*) = conj tau(c) ;  tau(c^*) = i^k2 y_c2  => y_c2 = i^-k2 conj(y_c)
                s.relate(c2,c,(-k2)%4,1)
            for t in range(12):
                sg,w2=act(t,toks)
                r=canon(0 if sg==1 else 2,w2)
                if r is None: s.zero.add(c); continue
                k3,c3=r
                if s.add(c3): fr.append(c3)
                # tau(t.w) = tau(w): i^k3 y_c3 = y_c  => y_c3 = i^-k3 y_c
                s.relate(c3,c,(-k3)%4,0)
    def parametrize(s):
        """exact real parametrization: each root r gets real parameters z_j with y_r = sum_j z_j beta_j (beta Gaussian ints)"""
        roots=sorted(set(s.find(c)[0] for c in s.par), key=repr)
        zroots=set(s.find(c)[0] for c in s.zero)
        s.param={}; n=0
        for r in roots:
            if r in zroots: s.param[r]=[]; continue
            # allowed set: y with y = i^k conj^f y for all self relations
            cands=[KO,KI]   # basis of C over R
            basis=None
            lin=[(k,f) for (k,f) in s.self[r]]
            # solve over R: y = u + i v ; relation y - i^k conj^f(y) = 0  (real 2x2 linear system)
            rows=[]
            for (k,f) in lin:
                ph=[(1,0),(0,1),(-1,0),(0,-1)][k]
                # i^k*(u+iv) or i^k*(u-iv)
                if f==0: M=[[1-ph[0], ph[1]],[-ph[1],1-ph[0]]]   # (u,v) -> real,imag of y - ph*y
                else:    M=[[1-ph[0], -ph[1]],[-ph[1],1+ph[0]]]
                rows+=M
            # nullspace of rows (integer 2-column) exactly
            ns=_null2(rows)
            s.param[r]=[]
            for (u,v) in ns:
                s.param[r].append((n,K(u,v))); n+=1
        s.nz=n
        return n
    def expr(s,c):
        """y_c as list of (z index, K coefficient)"""
        if c in s.zero: return []
        r,(k,f)=s.find(c)
        out=[]
        for zi,beta in s.param[r]:
            b=beta.conj() if f else beta
            out.append((zi,IK[k]*b))
        return out
def _null2(rows):
    if not rows: return [(1,0),(0,1)]
    # rank of integer matrix with 2 columns
    import itertools
    nz=[r for r in rows if r[0]!=0 or r[1]!=0]
    if not nz: return [(1,0),(0,1)]
    a,b=nz[0]
    v=(-b,a)   # candidate null vector
    for r in nz:
        if r[0]*v[0]+r[1]*v[1]!=0: return []
    return [v]
