import mpmath as mp
mp.mp.dps = 60
from mst import MST, I
from polar import E_spheroidal
f = mp.factorial

def Hodd_full(s, l, m, cw, eps):
    cw, eps = mp.mpf(cw), mp.mpf(eps)
    q = 2*cw/eps; kabs = mp.sqrt(q**2-1)
    lam = E_spheroidal(s, l, m, cw) - 2*m*cw + cw**2 - s*(s+1)
    H = []
    for k in (I*kabs, -I*kabs):
        P = MST(eps, q, m, s, lam, k, N=30); P.solve_nu(l - mp.mpf('0.5')*eps**2)
        H.append(P.amplitudes()['H'])
    return (H[0]-H[1])/2

def H0(s, l, m, cw):        # closed-form leading eta-odd element, per unit eps
    return (-mp.mpf(1)/2*(-1)**m*(4*I*cw)**(2*l+1)*f(l+m)*f(l-m)*(f(l-s)*f(l+s))**2/(f(2*l)*f(2*l+1))**2)

def rich(func, eps0):       # 3-point Richardson in eps: removes O(eps) and O(eps^2)
    e = [eps0, eps0/2, eps0/4]; v = [func(x)/x for x in e]
    return (8*v[2] - 6*v[1] + v[0])/3

print("== l=2 (s=-2): H_odd/H0 vs 1 + (10m/9)c + 5(475 m^2+261)c^2/3969, residual/c^3 should be ~const ==")
for m in (2, 0, -2):
    for cw in ('0.01', '0.02', '0.03'):
        r = rich(lambda e: Hodd_full(-2, 2, m, cw, e), mp.mpf('4e-6'))/H0(-2, 2, m, mp.mpf(cw))
        c = mp.mpf(cw); ser = 1 + mp.mpf(10*m)/9*c + 5*(475*m**2+261)*c**2/3969
        print(f"  m={m:+d} c={cw}: ratio={mp.nstr(r, 14)}  series={mp.nstr(ser, 14)}  (ratio-series)/c^3={mp.nstr((r-ser)/c**3, 6)}")

print("\n== l=3 (s=-2): leading closed form  H0(3,m)/(eps c^7) = (8/893025) i (-1)^m (3+m)!(3-m)! ==")
for m in (2, 0, -1, 3):
    vals = []
    for cw in ('0.01', '0.02'):
        c = mp.mpf(cw)
        vals.append(rich(lambda e: Hodd_full(-2, 3, m, cw, e), mp.mpf('4e-6'))/c**7)
    lead = 2*vals[0] - vals[1]                         # remove O(c)
    print(f"  m={m:+d}: MST -> {mp.nstr(lead, 10)}   closed form {mp.nstr(H0(-2, 3, m, 1), 10)}")
