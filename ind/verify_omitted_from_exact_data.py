"""Fully exact check of three omitted SOS equations using a supplied Z[sqrt2,sqrt3] row-combination certificate.
No numerical solver, numpy, scipy or flint required.
"""
import gzip,json,os
from fractions import Fraction as Fr
P=os.path.dirname(__file__)
z=json.load(gzip.open(os.path.join(P,'fast_exact_rows.json.gz'),'rt'))
rel=json.load(open(os.path.join(P,'omitted_row_relations_exact.json')))

def mul(u,v):
 a,b,c,d=u;e,f,g,h=v
 return (a*e+2*b*f+3*c*g+6*d*h,a*f+b*e+3*c*h+3*d*g,a*g+c*e+2*b*h+2*d*f,a*h+d*e+b*g+c*f)
# note above last component is a*h+d*e+b*g+c*f, not an approximation.
R=[{j:tuple(v) for j,v in row} for row,b in z['allrows']]
B=[[Fr(s) for s in b] for row,b in z['allrows']]
for stridx,cert in rel.items():
 idx=int(stridx);scale=cert['multiplier'];lhs={};rhs=[Fr(0)]*4
 for use_i,u in cert['terms']:
  for col,v in R[use_i].items():
   term=mul(u,v);p=lhs.get(col,(0,0,0,0));lhs[col]=tuple(x+y for x,y in zip(p,term))
  term=mul(u,B[use_i]);rhs=[x+y for x,y in zip(rhs,term)]
 for col in set(lhs)|set(R[idx]):
  assert lhs.get(col,(0,0,0,0))==tuple(scale*x for x in R[idx].get(col,(0,0,0,0))),(idx,col)
 assert rhs==[scale*x for x in B[idx]],(idx,'rhs')
 print(f'row {idx}: EXACT PASS with {len(cert["terms"])} terms, multiplier {scale}')
print('ALL OMITTED SOS ROWS EXACTLY IMPLIED')
