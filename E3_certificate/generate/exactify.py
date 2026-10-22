import pickle, numpy as np, time, scipy.linalg as sl
from fractions import Fraction
import flint
from flint import fmpq, fmpq_mat
from field import F, K, KZ
from linsys import *
D=pickle.load(open('cert_data.pkl','rb')); red=pickle.load(open('reduced.pkl','rb')); XN=pickle.load(open('xnum.pkl','rb'))
OBJ=D['data']['OBJ']; g=D['data']['g']; nz=D['data']['nz']
rows,rhs,lay,nu=build_system(red,OBJ,g,nz)
keep=list(XN['keep']); xnum=XN['x']
A=np.zeros((len(rows),nu))
for i,R in enumerate(rows):
    for j,v in R.items(): A[i,j]=float(v)
Ak=A[keep]
Qr,Rr,Pc=sl.qr(Ak,pivoting=True,mode='economic')
piv=sorted(Pc[:len(keep)]); pivset=set(piv); free=[j for j in range(nu) if j not in pivset]
print('pivot cols',len(piv),'free',len(free),'cond(A_piv)~',np.linalg.cond(Ak[:,piv]))
def rat(v,den=10**12): return fmpq(int(round(v*den)),den)
x=[None]*nu
for j in free: x[j]=F(rat(xnum[j]))
# rhs' = b - A_free x_free  (exact over F)
t0=time.time()
rhs2=[]
for i in keep:
    acc=rhs[i]
    for j,v in rows[i].items():
        if j not in pivset: acc=acc-v*x[j]
    rhs2.append(acc)
# regular representation
n=len(piv); pidx={j:k for k,j in enumerate(piv)}
M=fmpq_mat(4*n,4*n); bvec=fmpq_mat(4*n,1)
for r_,i in enumerate(keep):
    for j,v in rows[i].items():
        if j not in pivset: continue
        c=pidx[j]; a0,a1,a2,a3=v.a,v.b,v.c,v.d
        blk=[[a0,2*a1,3*a2,6*a3],[a1,a0,3*a3,3*a2],[a2,2*a3,a0,2*a1],[a3,a2,a1,a0]]
        for p in range(4):
            for q in range(4):
                if blk[p][q]!=0: M[4*r_+p,4*c+q]=M[4*r_+p,4*c+q]+blk[p][q]
    rr=rhs2[r_]
    for p,val in enumerate((rr.a,rr.b,rr.c,rr.d)): bvec[4*r_+p,0]=val
print('assembled',round(time.time()-t0,1),'s; solving',4*n,'x',4*n,'rational system',flush=True)
t0=time.time(); sol=M.solve(bvec); print('solved',round(time.time()-t0,1),'s',flush=True)
for c,j in enumerate(piv):
    x[j]=F(sol[4*c,0],sol[4*c+1,0],sol[4*c+2,0],sol[4*c+3,0])
dev=max(abs(float(x[j])-xnum[j]) for j in piv)
print('max |exact pivot - numerical| =',dev)
# exact verification of ALL rows
t0=time.time(); bad=0
for i,R in enumerate(rows):
    acc=F(0)
    for j,v in R.items(): acc=acc+v*x[j]
    if not (acc-rhs[i]).iszero(): bad+=1
print('exact identity rows violated:',bad,'of',len(rows),round(time.time()-t0,1),'s')
pickle.dump(dict(x=x,lay=lay),open('xexact.pkl','wb'))
