"""Absorption content of the exchange-odd part: Gamma_pm = 1 - |S_pm|^2 for the physical solution and its image.
Check (Gamma_+ - Gamma_-)/eps^6 -> -4 A6^-  and the q = 0 values (Page: Gamma_+ = 4 eps^6/225)."""
import mpmath as mp
mp.mp.dps = 50
from mst import MST, I
from phys import spheroidal_c
from csm import secant, D_TS
def S(m, q, eps, sign):
    q, eps = mp.mpf(q), mp.mpf(eps); c = q*eps/2
    E, _ = spheroidal_c(-2, 2, m, c); lam = E - 2*m*c + c**2 - 2
    k0 = mp.sqrt(1 - q**2)
    k, lk = (k0, mp.log(k0)) if sign > 0 else (-k0, mp.log(k0) + I*mp.pi)
    P = MST(eps, q, m, -2, lam, k, logkappa=lk, N=40)
    x = secant(lambda x: P.g(2 + x*eps**2)/eps**2, mp.mpf(-107)/210, mp.mpf(10)**(-35))
    P.nu = 2 + x*eps**2
    A = P.amplitudes()
    return (-1)**3*D_TS(lam, c, m)/(16*(eps/2)**4)*A['S']
def A6m(m, q):
    q = mp.mpf(q); k = mp.sqrt(1 - q**2)
    if m == 0 or q == 0: return -k**5/225
    w = m*q/k; return -(mp.pi/900)*k**5*w*(1 + w**2)*(4 + w**2)*(1/mp.tanh(mp.pi*w) - 1)
for m, q in ((2, '0'), (0, '0.3'), (2, '0.3'), (-1, '0.5')):
    vals = {}
    for e in ('2e-3', '1e-3'):
        ee = mp.mpf(e)
        Gp, Gm = 1 - abs(S(m, q, e, 1))**2, 1 - abs(S(m, q, e, -1))**2
        vals[e] = (Gp - Gm)/ee**6
        print(f"m={m:+d} q={q} eps={e}: Gamma+/eps^6 = {mp.nstr(Gp/ee**6, 10):>14}  Gamma-/eps^6 = {mp.nstr(Gm/ee**6, 10):>14}"
              f"  (G+ - G-)/eps^6 = {mp.nstr((Gp - Gm)/ee**6, 10):>14}   -4 A6^- = {mp.nstr(-4*A6m(m, q), 10)}")
    rich = 2*vals['1e-3'] - vals['2e-3']             # removes the O(eps) correction
    print(f"   linear extrapolation eps -> 0: {mp.nstr(rich, 10)}   |difference| = {mp.nstr(abs(rich + 4*A6m(m, q)), 3)}")
print("Page: 4/225 =", mp.nstr(mp.mpf(4)/225, 10))
