import mpmath as mp
mp.mp.dps = 60
from mst import MST, I
from polar import E_spheroidal          # reuse the spectral eigenvalue routine (module guard below)

def H_odd(s, l, m, aw, eps):
    aw, eps = mp.mpf(aw), mp.mpf(eps)
    q = 2*aw/eps; kabs = mp.sqrt(q**2-1)
    lam = E_spheroidal(s, l, m, aw) - 2*m*aw + aw**2 - s*(s+1)
    H = []
    for k in (I*kabs, -I*kabs):
        P = MST(eps, q, m, s, lam, k, N=25); P.solve_nu(l - mp.mpf('0.5')*eps**2)
        H.append(P.amplitudes()['H'])
    return (H[0]-H[1])/2

def coeffs(s, l, m):
    r = []
    for aw in ('0.004', '0.008'):
        e1, e2 = mp.mpf('2e-5'), mp.mpf('1e-5')
        h0 = 2*H_odd(s, l, m, aw, e2)/e2 - H_odd(s, l, m, aw, e1)/e1
        r.append(h0/mp.mpf(aw)**(2*l+1))
    C5 = 2*r[0] - r[1]; C6 = (r[1]-r[0])/mp.mpf('0.004')
    return C5, C6/C5

def closed_form(s, l, m):
    f = mp.factorial
    return (-mp.mpf(1)/2*(-1)**m*(4*I)**(2*l+1)*f(l+m)*f(l-m)*(f(l-s)*f(l+s))**2
            /(f(2*l)*f(2*l+1))**2)

# Expected slopes C6/C5: 10m/9 for s = -2 (Eq. (S7)). For s = 0 the slope vanishes: a three-point eps
# extrapolation at a*w = 0.005...0.025 gives C6/C5 = 0 to 1e-9 for m = 2, 1, -1.
CASES = [(-2, 2, 1, mp.mpf(10)/9), (-2, 2, 0, 0), (-2, 2, -1, -mp.mpf(10)/9), (-2, 2, -2, -mp.mpf(20)/9),
         (0, 2, 2, 0), (0, 2, -1, 0)]

if __name__ == '__main__':
    for (s, l, m, slope) in CASES:
        C5, k6 = coeffs(s, l, m)
        print(f"s={s:+d} l={l} m={m:+d}: C5 = {mp.nstr(C5, 10)}  closed form {mp.nstr(closed_form(s, l, m), 10)} | "
              f"C6/C5 = {mp.nstr(k6, 8)}  expected {mp.nstr(slope, 8)}", flush=True)
