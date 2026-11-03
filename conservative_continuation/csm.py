"""
Classical-limit (eps -> 0 at fixed c = a*omega) partial-wave S-matrix of super-extremal Kerr,
full far+near MST, both horizon branches.  M = 1, eps = 2 omega, q = a = 2c/eps.
Branch eta: kappa_eta = eta * i * q * sqrt(1 - 1/q^2)   (analytic in q; = eta*i*sqrt(q^2-1) for q>0).
Physical (parity-averaged, helicity-preserving) phase:  S = (-1)^(l+1) D/(16 om^4) * B_ref/B_inc.
"""
import mpmath as mp
from mst import MST, I
from polar import xmat


def spheroidal(s, l, m, c, Lextra=25):
    """Eigenvalue E (-> l(l+1) at c=0) and spherical-basis coefficients b_L (sum b_L^2 = 1, b_l > 0)."""
    lmin = max(abs(s), abs(m))
    ls, X = xmat(s, m, lmin, l + Lextra)
    n = len(ls)
    X2 = X*X
    H = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            H[i, j] = (ls[i]*(ls[i]+1) if i == j else 0) - c**2*X2[i, j] + 2*c*s*X[i, j]
    ev, Q = mp.eigsy(H)
    k = min(range(n), key=lambda i: abs(ev[i] - l*(l+1)))
    v = [Q[i, k] for i in range(n)]
    sg = 1 if v[ls.index(l)] > 0 else -1
    return ev[k], {L: sg*v[i] for i, L in enumerate(ls)}


def D_TS(lam, c, m):
    return mp.sqrt(lam**2*(lam+2)**2 - 8*lam*((c**2 - m*c)*(5*lam + 6) - 12*c**2) + 144*(c**2 - m*c)**2)


def secant(G, x0, rtol, maxit=40):
    xa = mp.mpmathify(x0); xb = xa*mp.mpf('1.0001') + mp.mpf('1e-6')
    fa, fb = G(xa), G(xb)
    best = min([(abs(fa), xa), (abs(fb), xb)], key=lambda t: t[0])
    for _ in range(maxit):
        if fb == fa:
            break
        xc = xb - fb*(xb - xa)/(fb - fa)
        fc = G(xc)
        if abs(fc) < best[0]:
            best = (abs(fc), xc)
        if abs(xc - xb) < rtol*max(1, abs(xc)):
            break
        xa, fa, xb, fb = xb, fb, xc, fc
    return best[1]


def solve_x(G, x0, rtol):
    """Robust root of the scaled characteristic function: bracket Re G on a grid around x0, then
    Illinois on Re G, then complex secant polish.  Returns (x, |G(x)|, |G| scale)."""
    x0 = mp.re(mp.mpmathify(x0))
    h = mp.mpf('0.0125')
    grid = [x0 + (j + mp.mpf('0.37'))*h for j in range(-48, 48)]
    grid = [x for x in grid if abs(x) > mp.mpf('1e-3')]
    vals = [mp.re(G(x)) for x in grid]
    scale = max(abs(v) for v in vals)
    brs = [(grid[j], grid[j+1], vals[j], vals[j+1]) for j in range(len(grid)-1) if vals[j]*vals[j+1] < 0
           and not (grid[j] < 0 < grid[j+1])]
    brs = [b for b in brs if max(abs(b[2]), abs(b[3])) < 0.2*scale]      # discard sign flips across poles
    a, b, fa, fb = min(brs, key=lambda t: abs((t[0]+t[1])/2 - x0))
    side = 0
    for _ in range(200):
        xc = (a*fb - b*fa)/(fb - fa)
        fc = mp.re(G(xc))
        if fc*fb < 0:
            a, fa = b, fb
            if side == -1: pass
            side = -1
        else:
            fa = fa/2
            side = 1
        b, fb = xc, fc
        if abs(b - a) < rtol*1e3 or fc == 0:
            break
    x = b
    # complex secant polish (imaginary parts are at the noise level)
    xa, xb = x, x*(1 + mp.mpf('1e-12'))
    ga, gb = G(xa), G(xb)
    for _ in range(6):
        if gb == ga: break
        xc = xb - gb*(xb - xa)/(gb - ga)
        xa, ga, xb, gb = xb, gb, xc, G(xc)
        if abs(xb - xa) < rtol: break
    if abs(gb) > abs(ga): xb, gb = xa, ga
    return xb, abs(gb), scale


def make(l, m, s, c, eta, eps, E, N=40):
    c, eps = mp.mpf(c), mp.mpf(eps)
    lam = E - 2*m*c + c**2 - s*(s+1)
    if c == 0:
        q, kappa = mp.mpf(0), mp.mpf(1)
    else:
        q = 2*c/eps
        kappa = eta*I*q*mp.sqrt(1 - 1/q**2)
    return MST(eps, q, m, s, lam, kappa, N=N), lam


def S_phys(l, m, s, c, eta, eps, x0, N=40, E=None, rtol=mp.mpf('1e-24')):
    """Return (S, x=(nu-l)/eps^2, MST amplitude dict)."""
    c, eps = mp.mpf(c), mp.mpf(eps)
    if E is None:
        E, _ = spheroidal(s, l, m, c)
    P, lam = make(l, m, s, c, eta, eps, E, N)
    Gf = lambda x: P.g(l + x*eps**2)/eps**2
    x = secant(Gf, x0, rtol)
    gres, gscale = abs(Gf(x)), abs(Gf(x + mp.mpf('0.01')))
    if not gres < mp.mpf(10)**(-mp.mp.dps//3)*gscale:          # fall back to the bracketing solver
        x, gres, gscale = solve_x(Gf, x0, rtol)
    P.nu = l + x*eps**2
    A = P.amplitudes()
    A['gres'], A['gscale'] = gres, gscale
    om = eps/2
    S = (-1)**(l+1)*D_TS(lam, c, m)/(16*om**4)*A['S']
    return S, x, A


def s1(S, eps):
    eps = mp.mpf(eps)
    return (S*mp.exp(-2*I*eps*mp.log(eps)) - 1)/eps


def nu2_schw(l, s=-2):
    l = mp.mpf(l)
    return (-2 - s**2/(l*(l+1)) + ((l+1)**2 - s**2)**2/((2*l+1)*(2*l+2)*(2*l+3))
            - (l**2 - s**2)**2/((2*l-1)*2*l*(2*l+1)))/(2*l+1)
