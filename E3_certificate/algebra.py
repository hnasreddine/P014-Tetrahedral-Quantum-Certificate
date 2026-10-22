"""Exact word algebra for the tracial moment certificate.
Generators: sigma_1, sigma_2, sigma_3 (Pauli matrices tensor identity) and four mutually orthogonal projections Q_0..Q_3.
Normalized trace tau on M_2(C) (x) M_N(C).  Phases are powers of i, stored as integers mod 4.
"""
import itertools
# Pauli products: PMUL[(a,b)] = (k, c) meaning sigma_a sigma_b = i^k sigma_c   (sigma_0 = identity)
PMUL={}
for a in range(4):
    PMUL[(0,a)]=(0,a); PMUL[(a,0)]=(0,a)
for a in (1,2,3): PMUL[(a,a)]=(0,0)
for (a,b,c) in ((1,2,3),(2,3,1),(3,1,2)):
    PMUL[(a,b)]=(1,c); PMUL[(b,a)]=(3,c)
def reduce_linear(k,toks):
    """toks: list of ('P',a) / ('Q',c). Returns (k,toks) reduced, or None if the word is 0."""
    out=[]
    for t in toks:
        if t[0]=='P':
            if t[1]==0: continue
            if out and out[-1][0]=='P':
                kk,c=PMUL[(out[-1][1],t[1])]; k=(k+kk)%4; out.pop()
                if c!=0: out.append(('P',c))
                # identity created: may make Q Q adjacent
                if c==0 and len(out)>=2 and out[-1][0]=='Q' and False: pass
            else: out.append(t)
        else:
            if out and out[-1][0]=='Q':
                if out[-1][1]==t[1]: continue
                return None
            out.append(t)
        # fix Q Q adjacency created by cancellation
        while len(out)>=2 and out[-1][0]=='Q' and out[-2][0]=='Q':
            if out[-1][1]==out[-2][1]: out.pop()
            else: return None
    return k,out
def canon(k,toks):
    """Canonical cyclic representative for tau(i^k * word). Returns (k, key) with key a tuple of (c,a) pairs,
    or (k,'I') for a multiple of the identity, or None if tau vanishes identically by the algebra."""
    r=reduce_linear(k,list(toks))
    if r is None: return None
    k,t=r
    qs=[i for i,x in enumerate(t) if x[0]=='Q']
    if not qs:
        if not t: return (k,'I')
        return None                       # tau(sigma_j) = 0
    i=qs[0]; t=t[i:]+t[:i]
    r=reduce_linear(k,t)
    if r is None: return None
    k,t=r
    if t[-1][0]=='Q' and len(t)>1:        # wrap-around Q...Q
        if t[-1][1]==t[0][1]: t=t[:-1]
        else: return None
    pairs=[];j=0
    while j<len(t):
        c=t[j][1]; a=0
        if j+1<len(t) and t[j+1][0]=='P': a=t[j+1][1]; j+=2
        else: j+=1
        pairs.append([c,a])
    changed=True
    while changed and len(pairs)>1:
        changed=False
        for q in range(len(pairs)):
            if pairs[q][1]==0:
                nq=(q+1)%len(pairs)
                if pairs[nq][0]!=pairs[q][0]: return None
                pairs[q][1]=pairs[nq][1]; pairs.pop(nq); changed=True; break
    pairs=[tuple(p) for p in pairs]
    rots=[tuple(pairs[q:]+pairs[:q]) for q in range(len(pairs))]
    return (k,min(rots))
def keytoks(key):
    toks=[]
    for c,a in key:
        toks.append(('Q',c))
        if a: toks.append(('P',a))
    return toks
def adjoint(k,toks): return ((-k)%4, list(reversed(toks)))
# ---- symmetry group: rotations of the tetrahedron r_c = s_c/sqrt3, s_0=(1,1,1), s_1=(1,-1,-1), s_2=(-1,1,-1), s_3=(-1,-1,1)
SIGNS=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
def _group():
    G=[]
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product([1,-1],repeat=3):
            R=[[0]*3 for _ in range(3)]
            for i in range(3): R[perm[i]][i]=sg[i]
            det=sg[0]*sg[1]*sg[2]*_permsign(perm)
            if det!=1: continue
            pm=[]
            for c in range(4):
                img=tuple(sum(R[r][j]*SIGNS[c][j] for j in range(3)) for r in range(3))
                if img not in SIGNS: break
                pm.append(SIGNS.index(img))
            if len(pm)==4: G.append((R,pm))
    return G
def _permsign(p):
    s=1
    for i in range(3):
        for j in range(i+1,3):
            if p[i]>p[j]: s=-s
    return s
GROUP=_group(); assert len(GROUP)==12
def act(t,toks):
    """apply group element t to a word: sigma_j -> sum_k R_kj sigma_k = s sigma_{pi j}, Q_c -> Q_{pm(c)}. returns (sign, toks)"""
    R,pm=GROUP[t]; s=1; out=[]
    for typ,x in toks:
        if typ=='P':
            col=[R[r][x-1] for r in range(3)]; kk=[r for r in range(3) if col[r]!=0][0]
            s*=col[kk]; out.append(('P',kk+1))
        else: out.append(('Q',pm[x]))
    return s,out
def stabilizer0(): return [t for t in range(12) if GROUP[t][1][0]==0]
def z3exp(t):
    R,_=GROUP[t]
    P=[[1 if R[i][j]!=0 else 0 for j in range(3)] for i in range(3)]
    if P==[[1,0,0],[0,1,0],[0,0,1]]: return 0
    if P[1][0]==1: return 1   # e0 -> e1
    return 2
