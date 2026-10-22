import pickle, numpy as np
from field import K, F, KZ, KO
from flint import fmpq
ORDER=['PZ0','PZ1','PZ2','QA0','QA1','QA2','QT']
def layout(red):
    lay={}; off=0
    for name in ORDER:
        r=len(red[name]); idx={}
        for a in range(r): idx[('d',a)]=off; off+=1
        for a in range(r):
            for b in range(a+1,r):
                idx[('re',a,b)]=off; off+=1; idx[('im',a,b)]=off; off+=1
        lay[name]=(r,idx)
    return lay,off
def build_system(red,OBJ,g,nz):
    """rows: 0 -> constant ('c'), 1+z -> moment parameter z.  Entries in F (exact)."""
    lay,nu=layout(red)
    rows=[dict() for _ in range(nz+1)]
    def rowid(k): return 0 if k=='c' else 1+k
    for name in ORDER:
        Cr=red[name]; r,idx=lay[name]
        for a in range(r):
            for k,v in Cr[a][a].items():
                assert v.im.iszero()
                R=rows[rowid(k)]; R[idx[('d',a)]]=R.get(idx[('d',a)],F(0))+v.re
            for b in range(a+1,r):
                for k,v in Cr[b][a].items():
                    R=rows[rowid(k)]
                    i1=idx[('re',a,b)]; i2=idx[('im',a,b)]
                    R[i1]=R.get(i1,F(0))+v.re*2
                    R[i2]=R.get(i2,F(0))-v.im*2
    rhs=[F(0)]*(nz+1)
    rhs[0]=g-OBJ.get('c',KZ).re
    for k,v in OBJ.items():
        if k!='c': rhs[1+k]=rhs[1+k]-v.re
    return rows,rhs,lay,nu
def Y_from_x(x,lay,name,conv=lambda t:t):
    r,idx=lay[name]
    Y=[[None]*r for _ in range(r)]
    for a in range(r):
        Y[a][a]=(conv(x[idx[('d',a)]]),conv(0*x[idx[('d',a)]]) if False else None)
    return Y
