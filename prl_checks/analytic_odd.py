"""Analytic control of the exchange-odd response.  Checks, with full MST (physical branch, real eps):
 (a) delta^2 X / [(i/1800) eps^7 kappa^5 p(tau)] = 1 + c1 eps + O(eps^2),  p(tau) = tau(1+tau^2)(4+tau^2), tau exact,
     X = K_{-nu-1}/K_nu, delta = nu - 2;  expect c1 = 5 m q/9.
 (b) (e^{i pi nu} + e^{-i pi nu} rho) eps/(2 i delta) = 1 + (pi - i nu2) eps + O(eps^2),  nu2 = -107/210.
 (c) the law X_-/X_+ = e^{i pi eta(2nu+1)} sin pi(nu - i tau)/sin pi(nu + i tau) (eta = +1) on the same data."""
import mpmath as mp
mp.mp.dps = 60
from mst import MST, I
from phys import spheroidal_c
from csm import secant
def run(m, q, eps):
    q, eps = mp.mpf(q), mp.mpf(eps)
    c = q*eps/2
    E, _ = spheroidal_c(-2, 2, m, c)
    lam = E - 2*m*c + c**2 - 2
    k = mp.sqrt(1 - q**2)
    P = MST(eps, q, m, -2, lam, k, N=40)
    x = secant(lambda x: P.g(2 + x*eps**2)/eps**2, mp.mpf(-107)/210, mp.mpf(10)**(-40))
    P.nu = 2 + x*eps**2
    A = P.amplitudes()
    d = P.nu - 2; t = P.tau
    R = d**2*A['X']/(1j/1800*eps**7*k**5*t*(1 + t**2)*(4 + t**2))
    Efac = (mp.exp(I*mp.pi*P.nu) + mp.exp(-I*mp.pi*P.nu)*A['rho'])*eps/(2j*d)
    # law: dual branch
    Pd = MST(eps, q, m, -2, lam, -k, logkappa=mp.log(k) + I*mp.pi, N=40); Pd.nu = P.nu
    Ad = Pd.amplitudes()
    law = mp.exp(I*mp.pi*(2*P.nu + 1))*mp.sin(mp.pi*(P.nu - I*t))/mp.sin(mp.pi*(P.nu + I*t))
    return R, Efac, abs(Ad['X']/A['X']/law - 1), x
nu2 = -mp.mpf(107)/210
for m, q in ((2, '0.3'), (1, '0.3'), (-1, '0.5'), (2, '0.5')):
    out = [(e,) + run(m, q, e) for e in ('2e-3', '1e-3')]
    (e1, R1, E1, L1, x1), (e2, R2, E2, L2, x2) = out
    e1, e2 = mp.mpf(e1), mp.mpf(e2)
    c1 = (2*(R2 - 1)/e2 - (R1 - 1)/e1)       # Richardson on (R-1)/eps
    cE = (2*(E2 - 1)/e2 - (E1 - 1)/e1)
    print(f"m={m:+d} q={q}: c1 = {mp.nstr(c1, 10)}  (5mq/9 = {mp.nstr(mp.mpf(5*m)/9*mp.mpf(q), 10)}) | "
          f"E-coefficient = {mp.nstr(cE, 10)}  (pi - i nu2 = {mp.nstr(mp.pi - 1j*nu2, 10)}) | law residual {mp.nstr(max(L1, L2), 3)} | x = {mp.nstr(x2, 8)}")
