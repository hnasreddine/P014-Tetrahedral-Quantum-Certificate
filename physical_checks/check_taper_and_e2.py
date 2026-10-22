"""Independent integer-bath and uniform-crossover checks, no external dependencies."""
import math
from fractions import Fraction
SQ2=math.sqrt(2); g=(3-2*SQ2)/3

def Phi(x): return 2*g if x==0 else math.log1p(2*g*x)/x

def weighted_taper(eps,x):
    q=2*(1+x*eps); theta=q/2; rho=eps**.25
    c=Phi(x)+eps**.125
    J=math.ceil(c/eps); H=math.ceil(eps**(-.5)); L=J+H
    d=1; t=1.
    rho_frac=Fraction(rho).limit_denominator(10**12)
    dn=[]; tvals=[]
    failure_numerator=0.;Z=0.; M1=M2=0.
    for j in range(L+1):
        if j<=J:
            if j>0: d*=2
        else:
            k=(rho_frac.numerator*d+rho_frac.denominator-1)//rho_frac.denominator
            assert 0<k<=d//2
            # j from J+1: use preceding d_j (save in loop)
            d=2*d-k
            t=t*(1-k/(2*d_prev))
        weight=t*theta**(-j)
        tvals.append(t);dn.append(d)
        Z+=weight;M1+=j*weight;M2+=j*j*weight
        if J<=j<L:
            d_prev=d
            k=(rho_frac.numerator*d+rho_frac.denominator-1)//rho_frac.denominator
            failure_numerator+=g*t*theta**(-j)*(k/d)
    failure_numerator+=tvals[-1]*theta**(-L)
    p=failure_numerator/Z
    B=sum(dn)
    logB=math.log2(B)
    mean=M1/Z; var=M2/Z-mean**2
    # weighted telescope equality check at q=2 uses exact relation, finite weighted for x !=0
    print(f'x={x:>3}, eps={eps:.4g},  phi={Phi(x):.9g}, eps log2B={eps*logB:.9g}, eps W/Delta={eps*L:.9g}, p/eps={p/eps:.9g}, eps E/Delta={eps*mean:.9g}, eps^2 C={eps*eps*math.log(q)**2*var:.9g}')
    assert p<eps,(p,eps)
    if x==0:
        assert abs(p-(2*g+(1-2*g)*tvals[-1])/Z)<1e-12
    return p

print('FINITE INTEGER BATHS:')
for x in [-1.,0.,1.]:
 for eps in [0.01,0.005,0.002]:
    weighted_taper(eps,x)

print('\nUNIFORM E2 LADDERS:')
for r in [2,3]:
 for delta in [1e-4,1e-6,1e-8]:
    eta=delta**.5
    h=math.log1p(eta); ell=math.ceil(math.log1p(eta/delta)/h)
    t=(1+eta)**(-1)
    error=1.0/math.exp(math.log(math.expm1(h*ell)/eta)) if h*ell<700 else eta/math.expm1(h*ell)
    lnB=math.log((r**1)/(r-1))+ell*math.log(r)+math.log1p(-math.exp(-ell*math.log(r))) # log((r^ell-1)/(r-1))
    leading=math.log(r)/h * math.log(eta/delta)
    print(f'r={r}, δ={delta:.0e}, eta={eta:.4g}, ell={ell}, error/δ={error/delta:.8g}, lnB/leading={lnB/leading:.8g}')
    assert error<=delta*1.000000001
print('PASS: finite integer taper and E2 geometrical ladder')
