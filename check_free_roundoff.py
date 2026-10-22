from fractions import Fraction as Q
import json,gzip, pathlib
base=pathlib.Path(__file__).parent
D=json.load(gzip.open(base/'E3_certificate'/'certificate_data.json.gz','rt'))
raw=(base/'prior'/'rounded_grams.txt').read_text().split()
ix=1; got={}; index=0; shift=Q(5,100000)
for name in D['blocks']:
 label=raw[ix]; n=int(raw[ix+1]);ix+=2; r=n//2
 assert name==label and r==D['block_sizes'][name]
 aa=list(map(int,raw[ix:ix+n*n]));ix+=n*n
 mat=[aa[k*n:(k+1)*n] for k in range(n)]
 for a in range(r):
  got[index]=Q(mat[a][a],10**11)+shift;index+=1
 for a in range(r):
  for b in range(a+1,r):
   got[index]=Q(mat[a][b],10**11);index+=1
   got[index]=Q(mat[a+r][b],10**11);index+=1
assert index==len(D['pivots'])+len(D['free'])
err= max((abs(got[int(j)]-Q(s)) for j,s in D['free'].items()))
print('Number gram variables',len(got));print('Exact max rounding gap fixed-free variables',err,'=',float(err))
print('5e-5/(2*max_r)',float(Q(5,100000)/(2*max(D['block_sizes'].values()))))
print('extra operator norm from free rounding at max block 30:',float(2*max(D['block_sizes'].values())*err))
