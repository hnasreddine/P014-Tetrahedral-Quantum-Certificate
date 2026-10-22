import sympy as s
r2=s.sqrt(2);r3=s.sqrt(3)
w=s.Rational(-1,2)+s.I*r3/2
v=[s.Matrix([1,0])]+[s.Matrix([1/r3,s.sqrt(s.Rational(2,3))*w**j]) for j in range(3)]
S=s.Matrix.hstack(s.Matrix([1,0,0,0]),s.Matrix([0,1/r2,1/r2,0]),s.Matrix([0,0,0,1]))
I=s.eye(2);K=[]
for vi in v:
 A=S.H*s.kronecker_product(vi,I)
 R=r2*I+(1-r2)*vi*vi.H
 Ki=s.simplify(A*R)
 assert s.simplify(Ki.H*Ki-I)==s.zeros(2)
 K.append(Ki)
P=s.Matrix(4,4,lambda i,j: s.simplify((v[j].H*v[i])[0]/2))
F=s.simplify(sum(P[j,i]*s.trace(K[i]*K[j].H) for i in range(4) for j in range(4))/4)
g=(3-2*r2)/3
print('F(2,3) =',F,'expected:',s.simplify((3+2*r2)/3))
assert s.simplify(F-(3+2*r2)/3)==0
assert s.simplify(2-F-g)==0
print('SYMPY 2->3 SATURATION VERIFIED EXACTLY')
