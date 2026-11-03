"""(c) X / X0 = 1 + (5/9) m q eps - (5/9) eps^2 + ..., X0 = leading-term MST ratio (all Gamma functions exact):
X0 = (2 eps kappa)^{2nu+1} G(-2nu)G(-2nu-1)/(G(2nu+1)G(2nu+2)) * G(nu+1-i tau)G(nu+1-s-i eps)G(nu+1-s+i eps)
     / (G(-nu-i tau)G(-nu-s-i eps)G(-nu-s+i eps)).   Analytic derivation: a_1 = (5i/18) eps kappa (3 - i tau), a_{-1} = O(eps^3)."""
import mpmath as mp
mp.mp.dps = 60
from mst import MST, I
from phys import spheroidal_c
from csm import secant
def X_and_X0(m, q, eps):
    q, eps = mp.mpf(q), mp.mpf(eps); c = q*eps/2; s = -2
    E, _ = spheroidal_c(s, 2, m, c); lam = E - 2*m*c + c**2 - s*(s+1)
    k = mp.sqrt(1 - q**2)
    P = MST(eps, q, m, s, lam, k, N=40)
    x = secant(lambda x: P.g(2 + x*eps**2)/eps**2, mp.mpf(-107)/210, mp.mpf(10)**(-45)); nu = 2 + x*eps**2; P.nu = nu
    A = P.amplitudes(); t = P.tau; G = mp.gamma
    X0 = (mp.exp((2*nu + 1)*mp.log(2*eps*k))*G(-2*nu)*G(-2*nu - 1)/(G(2*nu + 1)*G(2*nu + 2))
          *G(nu + 1 - I*t)*G(nu + 1 - s - I*eps)*G(nu + 1 - s + I*eps)/(G(-nu - I*t)*G(-nu - s - I*eps)*G(-nu - s + I*eps)))
    return A['X']/X0, x
for m, q in ((2, '0.3'), (1, '0.5'), (-2, '0.4'), (0, '0.6')):
    vals = []
    for e in ('4e-3', '2e-3', '1e-3'):
        r, x = X_and_X0(m, q, e); vals.append((mp.mpf(e), r, x))
    # fit r = 1 + c1 e + c2 e^2
    M = mp.matrix([[e, e**2] for e, _, _ in vals[:2]]); sol = mp.lu_solve(M, mp.matrix([r - 1 for _, r, _ in vals[:2]]))
    pred = mp.mpf(5*m)/9*mp.mpf(q)
    print(f"m={m:+d} q={q}: c1 = {mp.nstr(sol[0], 12)}  (5mq/9 = {mp.nstr(pred, 12)})   c2 = {mp.nstr(sol[1], 6)}  (-5/9 + ... )"
          f"   (nu-2)/eps^2 at eps=1e-3: {mp.nstr(vals[2][2], 10)}  -> nu3 ~ {mp.nstr((vals[2][2] + mp.mpf(107)/210)/mp.mpf('1e-3'), 6)} (mq/12 = {mp.nstr(m*mp.mpf(q)/12, 6)})")
