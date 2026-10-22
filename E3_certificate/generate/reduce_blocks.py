import pickle, time
from field import K, F, KZ, KO
def reduce_block(Cb,Bv):
    d=len(Cb); r=len(Bv)
    # T[i][b] = sum_j C[i][j] B[b][j]
    T=[[{} for _ in range(r)] for _ in range(d)]
    for i in range(d):
        for j in range(d):
            e=Cb[i][j]
            if not e: continue
            for b in range(r):
                x=Bv[b][j]
                if x.iszero(): continue
                acc=T[i][b]
                for k,v in e.items(): acc[k]=acc.get(k,KZ)+v*x
    Cr=[[{} for _ in range(r)] for _ in range(r)]
    for a in range(r):
        for i in range(d):
            y=Bv[a][i]
            if y.iszero(): continue
            yc=y.conj()
            for b in range(r):
                e=T[i][b]
                if not e: continue
                acc=Cr[a][b]
                for k,v in e.items(): acc[k]=acc.get(k,KZ)+yc*v
    for a in range(r):
        for b in range(r):
            Cr[a][b]={k:v for k,v in Cr[a][b].items() if not v.iszero()}
    return Cr
if __name__=='__main__':
    D=pickle.load(open('cert_data.pkl','rb')); FD=pickle.load(open('face_data.pkl','rb'))
    red={}
    for name in D['blocks']:
        t0=time.time(); Cr=reduce_block(D['blocks'][name],FD['B'][name])
        # exact Hermiticity check
        r=len(Cr); ok=True
        for a in range(r):
            for b in range(r):
                ka=set(Cr[a][b])|set(Cr[b][a])
                for k in ka:
                    if not (Cr[a][b].get(k,KZ)-Cr[b][a].get(k,KZ).conj()).iszero(): ok=False
        red[name]=Cr; print(name,'reduced size',r,'hermitian exact',ok,round(time.time()-t0,1),'s',flush=True)
    OBJ=D['data']['OBJ']
    print('objective coefficients real:',all(v.im.iszero() for v in OBJ.values()))
    pickle.dump(red,open('reduced.pkl','wb'))
