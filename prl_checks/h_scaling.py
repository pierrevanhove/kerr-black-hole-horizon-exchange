"""Scaling of the near-zone factor H = (1+i e^{i pi nu}X)/(1-i e^{-i pi nu} rho X) - 1 at fixed q on both branches,
and of the exchange-even/odd parts: is there an eps^4/nu2 piece in the even part that the far zone must cancel?"""
import mpmath as mp
mp.mp.dps = 60
from mst import MST, I
from phys import spheroidal_c
from csm import secant
def H(m, q, eps, sign):
    q, eps = mp.mpf(q), mp.mpf(eps); c = q*eps/2
    E, _ = spheroidal_c(-2, 2, m, c); lam = E - 2*m*c + c**2 - 2
    k0 = mp.sqrt(1 - q**2)
    k, lk = (k0, mp.log(k0)) if sign > 0 else (-k0, mp.log(k0) + I*mp.pi)
    P = MST(eps, q, m, -2, lam, k, logkappa=lk, N=40)
    x = secant(lambda x: P.g(2 + x*eps**2)/eps**2, mp.mpf(-107)/210, mp.mpf(10)**(-40))
    P.nu = 2 + x*eps**2
    A = P.amplitudes()
    return A['H'], A['Sfar']
m, q = 2, '0.3'
for e in ('1e-3', '5e-4', '2.5e-4'):
    Hp, Sp = H(m, q, e, 1); Hm, Sm = H(m, q, e, -1)
    ev, od = (Hp + Hm)/2, (Hp - Hm)/2
    ee = mp.mpf(e)
    print(f"eps={e}: |H_even|/eps^4 = {mp.nstr(abs(ev)/ee**4, 8)}  |H_even|/eps^5 = {mp.nstr(abs(ev)/ee**5, 8)}  |H_odd|/eps^6 = {mp.nstr(abs(od)/ee**6, 8)}"
          f"  |Sfar_+ - Sfar_-| = {mp.nstr(abs(Sp - Sm), 3)}  Sfar-1 ~ {mp.nstr(Sp - 1, 6)}")
