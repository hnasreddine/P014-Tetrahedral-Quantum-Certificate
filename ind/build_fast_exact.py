"""Pure Python integer exact coefficient builder for Claude E3 SOS.
Stores all entries in Z[sqrt2,sqrt3] using a uniform W denominator.
This is independent of FLINT. No float arithmetic occurs.
"""
import sys,os,types,json,gzip,time,math
from fractions import Fraction as Fr
SRC=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','E3_certificate'))
sys.path.insert(0,SRC)
class E:
 __slots__=('v',)
 def __init__(self,a=0,b=0,c=0,d=0):
  if isinstance(a,E) and b==c==d==0:self.v=a.v
  else:self.v=tuple(x if isinstance(x,int) else x.v[0] for x in (a,b,c,d))
 @classmethod
 def raw(cls,v):
  z=object.__new__(cls);z.v=tuple(v);return z
 def __add__(a,b):
  if not isinstance(b,E):b=E(b)
  return E.raw([a.v[i]+b.v[i] for i in range(4)])
 __radd__=__add__
 def __neg__(a):return E.raw([-e for e in a.v])
 def __sub__(a,b):return a+(-b if isinstance(b,E) else E(-b))
 def __rsub__(a,b):return -a+b
 def __mul__(a,b):
  if not isinstance(b,E):return E.raw([e*b for e in a.v])
  u,v=a.v,b.v
  return E.raw([u[0]*v[0]+2*u[1]*v[1]+3*u[2]*v[2]+6*u[3]*v[3],
   u[0]*v[1]+u[1]*v[0]+3*u[2]*v[3]+3*u[3]*v[2],
   u[0]*v[2]+u[2]*v[0]+2*u[1]*v[3]+2*u[3]*v[1],
   u[0]*v[3]+u[3]*v[0]+u[1]*v[2]+u[2]*v[1]])
 __rmul__=__mul__
 def __eq__(a,b):return a.v==(b.v if isinstance(b,E) else E(b).v)
 def iszero(a):return a.v==(0,0,0,0)
class C:
 __slots__=('re','im')
 def __init__(self,re=0,im=0):self.re=re if isinstance(re,E) else E(re);self.im=im if isinstance(im,E) else E(im)
 def __add__(a,b):
  if not isinstance(b,C):b=C(b)
  return C(a.re+b.re,a.im+b.im)
 __radd__=__add__
 def __neg__(a):return C(-a.re,-a.im)
 def __sub__(a,b):return a+(-b if isinstance(b,C) else C(-b))
 def __mul__(a,b):
  if isinstance(b,C):return C(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re)
  return C(a.re*b,a.im*b)
 __rmul__=__mul__
 def conj(a):return C(a.re,-a.im)
 def iszero(a):return a.re.iszero() and a.im.iszero()
 def __eq__(a,b):return (a-b).iszero()
FO=E(1);KO=C(1);KZ=C(0);KI=C(0,1)
def fmpq(num=0,den=1):
 assert den==1
 return int(num)
fakeflint=types.ModuleType('flint');fakeflint.fmpq=fmpq;sys.modules['flint']=fakeflint
modfield=types.ModuleType('field');modfield.F=E;modfield.K=C;modfield.KZ=KZ;modfield.KO=KO;modfield.KI=KI;modfield.kgauss=lambda z:C(round(z.real),round(z.imag));sys.modules['field']=modfield
from algebra import canon,adjoint,SIGNS
from moments import Moments,IK
cert=json.load(gzip.open(os.path.join(SRC,'certificate_data.json.gz'),'rt'));ORDER=cert['blocks']
dens=set()
for n in ORDER:
 for col in cert['W'][n]:
  for el in col:
   for x in el:dens.add(Fr(x).denominator)
DEN=math.lcm(*dens)
print('Common W denominator',DEN,flush=True)
def fparse(l):return E(*[int(Fr(x)*DEN) for x in l])
def kparse(l):return C(fparse(l[:4]),fparse(l[4:]))
def words(which):
 L=[()]
 if which=='Q':L += [(('P',j),) for j in (1,2,3)]
 L += [(('P',a),('Q',c),('P',b)) for c in range(4) for a in (1,2,3) for b in (1,2,3) if a!=b]
 L += [(('Q',c),('P',a),('Q',d),('P',b)) for c in range(4) for d in range(4) if c!=d for a in (1,2,3) for b in (1,2,3) if a!=b]
 return L
L={'P':words('P'),'Q':words('Q')}
def sector_terms(u,v,sec):
 _,ua=adjoint(0,list(u));base=ua+list(v)
 if sec=='P':return [(KO,base+[('Q',0)])]
 return [(KO,base)]+[(C(-1),base+[('Q',c)]) for c in range(4)]
print('Starting raw',flush=True);raw={};keys=set()
for sec in ('P','Q'):
 for i,u in enumerate(L[sec]):
  for j,v in enumerate(L[sec]):
   lst=[]
   for co,toks in sector_terms(u,v,sec):
    r=canon(0,toks)
    if r is None:continue
    k,key=r;lst.append((co*IK[k],key))
    if key!='I':keys.add(key)
   raw[(sec,i,j)]=lst
objraw=[]
for c in range(4):
 r=canon(0,[('Q',c)]);objraw.append((C(E(9,-8)),r))  # 6*(4g-1)/2
 for j in (1,2,3):
  r=canon(0,[('P',j),('Q',c)]);objraw.append((C(E(0,0,SIGNS[c][j-1])),r))  # 6*r_cj/2
objraw=[(co*IK[r[0]],r[1]) for co,r in objraw if r is not None]
for _,key in objraw:
 if key!='I':keys.add(key)
print('Moments',flush=True);mom=Moments();mom.close(keys);nz=mom.parametrize()
def lin(lst):
 d={}
 for co,key in lst:
  if key=='I':d['c']=d.get('c',KZ)+co;continue
  for zi,b in mom.expr(key):d[zi]=d.get(zi,KZ)+co*b
 return {k:v for k,v in d.items() if not v.iszero()}
print('Starting lin',len(raw),flush=True);Crows={k:lin(v) for k,v in raw.items()};OBJ6=lin(objraw)
assert all(v.im.iszero() for v in OBJ6.values())
W={n:[[kparse(e) for e in col] for col in cert['W'][n]] for n in ORDER}
def reduce(sec,Wc):
 n=len(L[sec]);r=len(Wc)
 T=[[{} for _ in range(r)] for _ in range(n)]
 for i in range(n):
  for j in range(n):
   ent=Crows.get((sec,i,j))
   if not ent:continue
   for b in range(r):
    x=Wc[b][j]
    if x.iszero():continue
    acc=T[i][b]
    for k,v in ent.items():acc[k]=acc.get(k,KZ)+v*x
 Cr=[[{} for _ in range(r)] for _ in range(r)]
 for a in range(r):
  for i in range(n):
   y=Wc[a][i]
   if y.iszero():continue
   yc=y.conj()
   for b in range(r):
    for k,v in T[i][b].items():Cr[a][b][k]=Cr[a][b].get(k,KZ)+yc*v
 for a in range(r):
  for b in range(r):Cr[a][b]={k:v for k,v in Cr[a][b].items() if not v.iszero()}
 return Cr
print('Starting reduced block coefficients',flush=True);red={}
for n in ORDER:
 t=time.time();red[n]=reduce(n[0],W[n]);print('block',n,'seconds',round(time.time()-t,1),flush=True)
lay={};nu=0
for n in ORDER:
 r=len(red[n]);idx={}
 for a in range(r):idx[('d',a)]=nu;nu+=1
 for a in range(r):
  for b in range(a+1,r):idx[('re',a,b)]=nu;nu+=1;idx[('im',a,b)]=nu;nu+=1
 lay[n]=(r,idx)
rows=[{} for _ in range(nz+1)];rid=lambda k:0 if k=='c' else 1+k
for n in ORDER:
 Cr=red[n];r,idx=lay[n]
 for a in range(r):
  for k,v in Cr[a][a].items():
   R=rows[rid(k)];j=idx[('d',a)];R[j]=R.get(j,E(0))+v.re
  for b in range(a+1,r):
   for k,v in Cr[b][a].items():
    R=rows[rid(k)];i1=idx[('re',a,b)];i2=idx[('im',a,b)]
    R[i1]=R.get(i1,E(0))+v.re*2;R[i2]=R.get(i2,E(0))-v.im*2
rhs6=[E(0)]*(nz+1);rhs6[0]=E(6,-4)-OBJ6.get('c',KZ).re
for k,v in OBJ6.items():
 if k!='c':rhs6[1+k]=rhs6[1+k]-v.re
piv=cert['pivots'];pivset=set(piv);keep=cert['rows_used'];npv=len(piv);pidx={j:t for t,j in enumerate(piv)}
free={int(j):Fr(s) for j,s in cert['free'].items()}
# Equations: 6* A_N / DEN^2 * Y = RHS6, yielding 6 A_N* pivot = DEN^2 * RHS6 - 6 A_N*free.
full=[]
for i,R in enumerate(rows):
 acc=[Fr(k*DEN**2) for k in rhs6[i].v]
 for j,v in R.items():
  if j not in pivset:
   val=free[j];acc=[x - 6*Fr(y)*val for x,y in zip(acc,v.v)]
 full.append(([(pidx[j], [6*k for k in v.v]) for j,v in R.items() if j in pivset],acc))
out={'W_DEN':DEN,'pivots':piv,'used':keep,'blocks':ORDER,'sizes':cert['block_sizes'],
     'allrows':[[[[int(j),v] for j,v in row], [str(z) for z in b]] for row,b in full]}
p=os.path.join(os.path.dirname(__file__),'fast_exact_rows.json.gz')
with gzip.open(p,'wt') as f:json.dump(out,f)
print('DONE',len(full),'rows',npv,'pivots',sum(len(r) for r,b in full),'entries',os.stat(p).st_size,flush=True)
