"""
Physical-region (|a| < M) double expansion of the Teukolsky partial-wave S-matrix in (eps, q),
eps = 2 M omega, q = a/M, spheroidal parameter c = a omega = q eps/2.  M = 1.
Branch: sign=+1 -> kappa = +sqrt(1-q^2) (outer horizon, physical);
        sign=-1 -> kappa = -sqrt(1-q^2) (inner-horizon boundary condition = irreducible-mass dual),
        with log kappa = log sqrt(1-q^2) + i pi (continuous on the q-circle).
Everything is analytic in (eps, q) near (0,0) after removing exp(2 i eps ln eps), up to O(eps^8 ln eps).
"""
import mpmath as mp
from mst import MST, I
from polar import xmat
from csm import D_TS, secant


def spheroidal_c(s, l, m, c, Lextra=10):
    """Complex-c spin-weighted spheroidal eigenvalue and spherical-basis vector (sum b^2 = 1, b_l ~ +1)."""
    lmin = max(abs(s), abs(m))
    ls, X = xmat(s, m, lmin, l + Lextra)
    n = len(ls)
    X2 = X*X
    H = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            H[i, j] = (ls[i]*(ls[i]+1) if i == j else 0) - c**2*X2[i, j] + 2*c*s*X[i, j]
    E, ER = mp.eig(H)
    k = min(range(n), key=lambda i: abs(E[i] - l*(l+1)))
    v = [ER[i, k] for i in range(n)]
    nrm = mp.sqrt(mp.fsum(x**2 for x in v))
    v = [x/nrm for x in v]
    if mp.re(v[ls.index(l)]) < 0:
        v = [-x for x in v]
    return E[k], {L: v[i] for i, L in enumerate(ls)}


def S_tilde(l, m, s, q, eps, sign, x0=None, N=20, rtol=None):
    """exp(-2 i eps ln eps) * (-1)^(l+1) D/(16 w^4) B_ref/B_inc  and the spheroidal vector."""
    q, eps = mp.mpmathify(q), mp.mpmathify(eps)
    c = q*eps/2
    E, b = spheroidal_c(s, l, m, c)
    lam = E - 2*m*c + c**2 - s*(s+1)
    k0 = mp.sqrt(1 - q**2)
    if sign > 0:
        kappa, logk = k0, mp.log(k0)
    else:
        kappa, logk = -k0, mp.log(k0) + I*mp.pi
    P = MST(eps, q, m, s, lam, kappa, logkappa=logk, N=N)
    if x0 is None:
        lm = mp.mpf(l)
        x0 = (-2 - s**2/(lm*(lm+1)) + ((lm+1)**2 - s**2)**2/((2*lm+1)*(2*lm+2)*(2*lm+3))
              - (lm**2 - s**2)**2/((2*lm-1)*2*lm*(2*lm+1)))/(2*lm+1)
    rtol = rtol or mp.mpf(10)**(-mp.mp.dps//2)
    Gf = lambda x: P.g(l + x*eps**2)/eps**2
    x = secant(Gf, x0, rtol)
    P.nu = l + x*eps**2
    A = P.amplitudes()
    om = eps/2
    S = (-1)**(l+1)*D_TS(lam, c, m)/(16*om**4)*A['S']
    return S*mp.exp(-2*I*eps*mp.log(eps)), b, x, abs(Gf(x))
