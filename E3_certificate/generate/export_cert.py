"""Export the certificate in a portable exact format (JSON, rationals as strings).
Stored: for each block rho: W_rho = Q_rho B_rho (word-coordinate basis, exact entries in Q(sqrt2,sqrt3,i)),
the rounded rational values of the 'free' Gram unknowns, and the list of 'pivot' unknowns that the verifier
re-solves exactly."""
import pickle, json, gzip, numpy as np, scipy.linalg as sl
from field import K, F, KZ
from linsys import ORDER, build_system
D=pickle.load(open('cert_data.pkl','rb')); FD=pickle.load(open('face_data.pkl','rb'))
red=pickle.load(open('reduced.pkl','rb')); XN=pickle.load(open('xnum.pkl','rb'))
def fstr(f): return [str(f.a),str(f.b),str(f.c),str(f.d)]
def kstr(k): return fstr(k.re)+fstr(k.im)
W={}
for name in ORDER:
    Q=D['SB'][name]; B=FD['B'][name]; n=len(Q[0]); d=len(Q); r=len(B)
    # W[:,b] = sum_a B[b][a] Q[a]
    cols=[]
    for b in range(r):
        col=[KZ]*n
        for a in range(d):
            if B[b][a].iszero(): continue
            for i in range(n):
                if not Q[a][i].iszero(): col[i]=col[i]+B[b][a]*Q[a][i]
        cols.append([kstr(v) for v in col])
    W[name]=cols
OBJ=D['data']['OBJ']; g=D['data']['g']; nz=D['data']['nz']
rows,rhs,lay,nu=build_system(red,OBJ,g,nz)
A=np.zeros((len(rows),nu))
for i,R in enumerate(rows):
    for j,v in R.items(): A[i,j]=float(v)
keep=list(XN['keep']); Qr,Rr,Pc=sl.qr(A[keep],pivoting=True,mode='economic'); piv=sorted(int(p) for p in Pc[:len(keep)])
DEN=10**12
free={}
for j in range(nu):
    if j not in set(piv): free[j]=str(int(round(XN['x'][j]*DEN)))+'/'+str(DEN)
cert=dict(description='E3 sector SOS certificate (see note)',blocks=ORDER,
          block_sizes={n:len(W[n]) for n in ORDER},W=W,pivots=piv,free=free,rows_used=[int(k) for k in keep])
with gzip.open('certificate_data_regenerated.json.gz','wt') as f: json.dump(cert,f)
print('exported; unknowns',nu,'pivots',len(piv),'free',len(free))
