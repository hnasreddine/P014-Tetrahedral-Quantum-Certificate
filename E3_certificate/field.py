"""Exact arithmetic in F = Q(sqrt2, sqrt3) and K = F(i), using FLINT rationals.
F element: (a,b,c,d) meaning a + b*sqrt2 + c*sqrt3 + d*sqrt6.
K element: (re, im) with re, im in F.
"""
import flint
from flint import fmpq, arb, ctx
Q0=fmpq(0); Q1=fmpq(1)
class F:
    __slots__=('a','b','c','d')
    def __init__(s,a=0,b=0,c=0,d=0):
        s.a=fmpq(a); s.b=fmpq(b); s.c=fmpq(c); s.d=fmpq(d)
    def __add__(s,o):
        if not isinstance(o,F): o=F(o)
        return F(s.a+o.a,s.b+o.b,s.c+o.c,s.d+o.d)
    __radd__=__add__
    def __neg__(s): return F(-s.a,-s.b,-s.c,-s.d)
    def __sub__(s,o):
        if not isinstance(o,F): o=F(o)
        return F(s.a-o.a,s.b-o.b,s.c-o.c,s.d-o.d)
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if not isinstance(o,F):
            o=fmpq(o); return F(s.a*o,s.b*o,s.c*o,s.d*o)
        a,b,c,d=s.a,s.b,s.c,s.d; e,f,g,h=o.a,o.b,o.c,o.d
        # (a + b r2 + c r3 + d r6)(e + f r2 + g r3 + h r6)
        A=a*e+2*b*f+3*c*g+6*d*h
        B=a*f+b*e+3*c*h+3*d*g
        C=a*g+c*e+2*b*h+2*d*f
        D=a*h+d*e+b*g+c*f
        return F(A,B,C,D)
    __rmul__=__mul__
    def iszero(s): return s.a==0 and s.b==0 and s.c==0 and s.d==0
    def conj2(s): return F(s.a,-s.b,s.c,-s.d)   # sqrt2 -> -sqrt2
    def conj3(s): return F(s.a,s.b,-s.c,-s.d)   # sqrt3 -> -sqrt3
    def inv(s):
        # multiply by Galois conjugates: N = s*s2*s3*s23 is rational
        s2=s.conj2(); s3=s.conj3(); s23=s2.conj3()
        p=s*s2*s3*s23
        assert p.b==0 and p.c==0 and p.d==0 and p.a!=0
        return (s2*s3*s23)*(Q1/p.a)
    def __truediv__(s,o):
        if not isinstance(o,F): return s*(Q1/fmpq(o))
        return s*o.inv()
    def arb(s,prec=200):
        ctx.prec=prec
        return arb(s.a)+arb(s.b)*arb(2).sqrt()+arb(s.c)*arb(3).sqrt()+arb(s.d)*arb(6).sqrt()
    def sign(s):
        if s.iszero(): return 0
        for prec in (128,512,2048,8192):
            v=s.arb(prec)
            if v>0: return 1
            if v<0: return -1
        raise ValueError('sign undetermined')
    def __float__(s): return float(s.a)+float(s.b)*2**0.5+float(s.c)*3**0.5+float(s.d)*6**0.5
    def __repr__(s): return f'F({s.a},{s.b},{s.c},{s.d})'
    def __eq__(s,o):
        if not isinstance(o,F): o=F(o)
        return (s-o).iszero()
    def comps(s): return (s.a,s.b,s.c,s.d)
FZ=F(0); FO=F(1)
R2=F(0,1,0,0); R3=F(0,0,1,0); R6=F(0,0,0,1)
class K:
    __slots__=('re','im')
    def __init__(s,re=0,im=0):
        s.re=re if isinstance(re,F) else F(re); s.im=im if isinstance(im,F) else F(im)
    def __add__(s,o):
        if not isinstance(o,K): o=K(o)
        return K(s.re+o.re,s.im+o.im)
    __radd__=__add__
    def __neg__(s): return K(-s.re,-s.im)
    def __sub__(s,o):
        if not isinstance(o,K): o=K(o)
        return K(s.re-o.re,s.im-o.im)
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if isinstance(o,K): return K(s.re*o.re-s.im*o.im, s.re*o.im+s.im*o.re)
        if isinstance(o,F): return K(s.re*o,s.im*o)
        o=fmpq(o); return K(s.re*o,s.im*o)
    __rmul__=__mul__
    def conj(s): return K(s.re,-s.im)
    def iszero(s): return s.re.iszero() and s.im.iszero()
    def inv(s):
        n=s.re*s.re+s.im*s.im; ni=n.inv(); return K(s.re*ni,-(s.im*ni))
    def __truediv__(s,o):
        if isinstance(o,K): return s*o.inv()
        if isinstance(o,F): return s*o.inv()
        return s*(Q1/fmpq(o))
    def __complex__(s): return complex(float(s.re),float(s.im))
    def __repr__(s): return f'K({s.re},{s.im})'
    def __eq__(s,o):
        if not isinstance(o,K): o=K(o)
        return (s-o).iszero()
KZ=K(0); KO=K(1); KI=K(0,1)
def kgauss(z):
    """exact Gaussian integer/unit from a python complex with integer parts"""
    re=round(z.real); im=round(z.imag)
    assert abs(z.real-re)<1e-12 and abs(z.imag-im)<1e-12, z
    return K(re,im)
OMEGA=K(fmpq(-1,2),F(0,0,fmpq(0),0))+K(0,F(0,0,fmpq(1,2),0))   # -1/2 + i sqrt3/2
assert (OMEGA*OMEGA*OMEGA)==KO
