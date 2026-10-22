import pickle, numpy as np, scipy.sparse as sp, clarabel, time
from linsys import *
D=pickle.load(open('cert_data.pkl','rb')); red=pickle.load(open('reduced.pkl','rb'))
OBJ=D['data']['OBJ']; g=D['data']['g']; nz=D['data']['nz']
rows,rhs,lay,nu=build_system(red,OBJ,g,nz)
A=np.zeros((nz+1,nu)); b=np.array([float(x) for x in rhs])
for i,R in enumerate(rows):
    for j,v in R.items(): A[i,j]=float(v)
print('system',A.shape,'rank',np.linalg.matrix_rank(A))
import scipy.linalg as sl
Qr,Rr,Pr=sl.qr(A.T,pivoting=True,mode='economic')
rk=int(np.sum(np.abs(np.diag(Rr))>1e-9*abs(Rr[0,0])))
keep=np.sort(Pr[:rk])
xl=np.linalg.lstsq(A[keep],b[keep],rcond=None)[0]
print('independent rows',rk,' consistency residual of all rows at lstsq solution',np.abs(A@xl-b).max())
A_full,b_full=A,b
A=A[keep]; b=b[keep]; nrow=len(keep)
# Clarabel: variables [x (nu), lam]; maximize lam ; A x = b ; E(Y_rho) - lam I in PSD
nv=nu+1
Arows=[sp.csr_matrix(np.hstack([A,np.zeros((nrow,1))]))]; bv=[b]; cones=[clarabel.ZeroConeT(nrow)]
for name in ORDER:
    r,idx=lay[name]; D2=2*r
    pos={};tt=0
    for col in range(D2):
        for row in range(col+1): pos[(row,col)]=tt;tt+=1
    G=sp.lil_matrix((tt,nv)); h=np.zeros(tt)
    def put(i,j,var,co):
        if i>j: i,j=j,i
        sc=1.0 if i==j else np.sqrt(2)
        G[pos[(i,j)],var]+= -co*sc      # s = h - G x  => entry = co*var
    for a in range(r):
        v=idx[('d',a)]; put(a,a,v,1.0); put(a+r,a+r,v,1.0)
        put(a,a,nu,-1.0); put(a+r,a+r,nu,-1.0)
        for bb in range(a+1,r):
            vr=idx[('re',a,bb)]; vi=idx[('im',a,bb)]
            # Y_ab = vr + i vi ; E = [[ReY,-ImY],[ImY,ReY]]
            put(a,bb,vr,1.0); put(a+r,bb+r,vr,1.0)
            put(a+r,bb,vi,1.0); put(bb+r,a,vi,-1.0)
    Arows.append(G.tocsr()); bv.append(h); cones.append(clarabel.PSDTriangleConeT(D2))
AA=sp.vstack(Arows).tocsc(); bb_=np.concatenate(bv)
q=np.zeros(nv); q[nu]=-1.0
st=clarabel.DefaultSettings(); st.verbose=False; st.max_iter=500
for k in ('tol_gap_abs','tol_gap_rel','tol_feas'): setattr(st,k,1e-12)
t0=time.time(); sol=clarabel.DefaultSolver(sp.csc_matrix((nv,nv)),q,AA,bb_,cones,st).solve()
x=np.array(sol.x); print('status',sol.status,'lambda* =',x[nu],'time',round(time.time()-t0,1))
print('equality residual (all rows)',np.abs(A_full@x[:nu]-b_full).max())
pickle.dump(dict(x=x[:nu],lam=x[nu],keep=keep),open('xnum.pkl','wb'))
