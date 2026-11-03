"""Segment representation of the eta-odd classical partial waves, Supplement Eq. (S16) (s = -2, m = 2).
Classical limit = first order in M at fixed a, omega.  At M = 0 the radial Teukolsky equation is the
spheroidal equation: R0(r) = Delta0(r) * S(i r/a), S = {}_{-2}S_{l2}(x; c) (entire), regular at both
horizons r = +-ia.  The continuation eta puts the ingoing condition at r = eta*ia (exponents {3,-1} at +ia,
degenerate {1,1} at -ia).  The Green identity of Delta^s T_M between the horizons gives
   e^{2i delta}|_odd = eps * G / (8 i omega^2 B_inc B_ref),
   G = FP Int_{-ia}^{ia} Delta0^{-2} R0 T1 R0 dr  -  (i/2a) R0'(-ia)^2 [1/(2ia) + R0''/R0' - i omega]_{r=-ia}
     = FP Int ...  +  S(1) [(2c - 1) S(1) - 4 S'(1)]      (endpoint term in angular form; checked below),
with T1 = dT/dM at M = 0, Eq. (S15).  Checks printed: radial ODE residual of R0, cancellation of the double
pole (C2) against the boundary Wronskian, vanishing residue (C1), S = 1 at M = 0, and the endpoint identity.
Run directly: compares with full MST at eps = 1e-12 (nonlinear_check.py) for l = 2.
vs_production.py compares l = 2..6 with the 110-digit data of ../conservative_test/polar_run2.json."""
import sys
import mpmath as mp
mp.mp.dps = 40
from csm import spheroidal
I = mp.mpc(0, 1)
s, m, a = -2, 2, mp.mpf(1)

def make_S(l, c, Lextra=34):
    """S(x) = sum_L b_L Y_L(x), Y_L the orthonormal {}_{-2}Y_{L2} built by the three-term recurrence of
    the same x-matrix that defines b_L (so conventions match); Y_2 = sqrt(5/64pi) (1+x)^2."""
    from polar import xmat
    E, b = spheroidal(s, l, m, c, Lextra)
    ls, X = xmat(s, m, max(abs(s), abs(m)), l + Lextra)
    def S(x, d=0):
        Y = [mp.sqrt(mp.mpf(5)/(64*mp.pi))*(1+x)**2]
        Y1 = [mp.sqrt(mp.mpf(5)/(64*mp.pi))*2*(1+x)]
        Y2 = [mp.sqrt(mp.mpf(5)/(64*mp.pi))*2]
        for i in range(len(ls)-1):
            prev = (Y[i-1], Y1[i-1], Y2[i-1]) if i > 0 else (0, 0, 0)
            off = X[i, i-1] if i > 0 else 0
            nx = X[i, i+1]
            Y.append(((x - X[i, i])*Y[i] - off*prev[0])/nx)
            Y1.append((Y[i] + (x - X[i, i])*Y1[i] - off*prev[1])/nx)
            Y2.append((2*Y1[i] + (x - X[i, i])*Y2[i] - off*prev[2])/nx)
        arr = (Y, Y1, Y2)[d]
        return mp.fsum(b[L]*arr[i] for i, L in enumerate(ls))
    return E, S

def setup(l, c):
    c = mp.mpf(c); om = c/a
    E, S = make_S(l, c)
    A = E - s*(s+1)
    lam = A + c**2 - 2*m*c
    # R0 = Delta0 * S(i r/a) and derivatives
    def R0(r, d=0):
        x = I*r/a; D0 = r**2 + a**2
        S0, S1, S2 = S(x), S(x, 1)*(I/a), S(x, 2)*(I/a)**2
        if d == 0: return D0*S0
        if d == 1: return 2*r*S0 + D0*S1
        return 2*S0 + 4*r*S1 + D0*S2
    K = lambda r: om*(r**2 + a**2) - a*m
    def T0(r):
        D0 = r**2 + a**2
        return (D0*R0(r, 2) + 2*(s+1)*r*R0(r, 1)
                + ((K(r)**2 - 2*I*s*r*K(r))/D0 + 4*I*s*om*r - lam)*R0(r))
    def T1R0(r):
        D0 = r**2 + a**2
        return (-2*r*R0(r, 2) - 2*(s+1)*R0(r, 1)
                + (2*I*s*K(r)/D0 + 2*r*(K(r)**2 - 2*I*s*r*K(r))/D0**2)*R0(r))
    F = lambda r: (r**2 + a**2)**s*R0(r)*T1R0(r)
    return dict(c=c, om=om, E=E, S=S, lam=lam, R0=R0, T0=T0, F=F, K=K)

# ---------------- asymptotic amplitudes of R0 by Taylor integration of the M = 0 ODE ----------------
def polys(om, lam):
    # P2 R'' + P1 R' + P0 R = 0  (multiplied by Delta0), coefficients as lists in powers of r
    a2 = a**2
    P2 = [a2**2, 0, 2*a2, 0, 1]
    P1 = [0, 2*(s+1)*a2, 0, 2*(s+1)]
    P0 = [om**2*a2**2 - 2*a**3*m*om + a2*m**2 - lam*a2,
          2*I*s*om*a2 + 2*I*s*a*m,
          2*om**2*a2 - 2*a*m*om - lam,
          2*I*s*om,
          om**2]
    return P2, P1, P0

def shift(P, rc):
    # coefficients of P(rc + t) in powers of t
    n = len(P); out = [mp.mpf(0)]*n
    for j, pj in enumerate(P):
        for k in range(j+1):
            out[k] += pj*mp.binomial(j, k)*rc**(j-k)
    return out

def taylor_step(Ps, rc, y0, y1, h, nterm=90):
    P2, P1, P0 = [shift(P, rc) for P in Ps]
    R = [y0, y1]
    for n in range(nterm):
        acc = 0
        for k in range(1, len(P2)):
            if n-k+2 >= 0: acc += P2[k]*(n-k+2)*(n-k+1)*R[n-k+2]
        for k in range(len(P1)):
            if n-k+1 >= 0: acc += P1[k]*(n-k+1)*R[n-k+1]
        for k in range(len(P0)):
            if n-k >= 0: acc += P0[k]*R[n-k]
        R.append(-acc/(P2[0]*(n+2)*(n+1)))
    v = mp.fsum(R[n]*h**n for n in range(len(R)))
    d = mp.fsum(n*R[n]*h**(n-1) for n in range(1, len(R)))
    return v, d

def asym_series(om, lam, sigma, K=40):
    p = -2*s - 1 if sigma == 1 else -1
    P2, P1, P0 = polys(om, lam)
    P2 = P2 + [0]*(5-len(P2)); P1 = P1 + [0]*(5-len(P1)); P0 = P0 + [0]*(5-len(P0))
    g = [mp.mpf(1)]
    def En(n, g):   # coefficient of r^{p+4-n}
        tot = 0
        for j in range(5):
            k = j - 6 + n
            if 0 <= k < len(g): tot += P2[j]*g[k]*(p-k)*(p-k-1)
            k = j - 5 + n
            if 0 <= k < len(g): tot += (2*I*sigma*om*P2[j] + P1[j])*g[k]*(p-k)
            k = j - 4 + n
            if 0 <= k < len(g): tot += (-om**2*P2[j] + I*sigma*om*P1[j] + P0[j])*g[k]
        return tot
    for n in range(1, K):
        g.append(mp.mpf(0))
        rest = En(n+1, g)                      # with g_n = 0
        g[n] = rest/(2*I*sigma*om*n)
        assert abs(En(n+1, g)) < mp.mpf(10)**(-mp.mp.dps+12)*max(1, abs(rest))
    def f(r):
        v = mp.exp(I*sigma*om*r)*mp.fsum(g[k]*r**(p-k) for k in range(K))
        dv = mp.exp(I*sigma*om*r)*mp.fsum(g[k]*((p-k)*r**(p-k-1) + I*sigma*om*r**(p-k)) for k in range(K))
        return v, dv
    return f

def amplitudes(P, r1=mp.mpf(1), rend=None):
    om, lam = P['om'], P['lam']
    rend = rend or mp.mpf(60)/om
    Ps = polys(om, lam)
    r, y0, y1 = r1, P['R0'](r1), P['R0'](r1, 1)
    while r < rend:
        h = min(mp.mpf('0.35')*mp.sqrt(r**2 + a**2), mp.mpf(2.5)/om, rend - r)
        y0, y1 = taylor_step(Ps, r, y0, y1, h)
        r += h
    fin, fref = asym_series(om, lam, -1), asym_series(om, lam, +1)
    vi, di = fin(r); vr, dr = fref(r)
    det = vi*dr - vr*di
    Binc = (y0*dr - vr*y1)/det
    Bref = (vi*y1 - y0*di)/det
    return Binc, Bref

def D_TS(lam, c):
    return mp.sqrt(lam**2*(lam+2)**2 - 8*lam*((c**2 - m*c)*(5*lam + 6) - 12*c**2) + 144*(c**2 - m*c)**2)

def S_odd(l, c, verbose=True):
    P = setup(l, c)
    c, om = P['c'], P['om']
    # checks: radial ODE residual at a few points
    res = max(abs(P['T0'](r))/abs(P['R0'](r)) for r in (mp.mpf('0.3'), I*mp.mpf('0.5'), mp.mpf('2.0')))
    # local data at -ia
    rm = -I*a
    Am = P['R0'](rm, 1); c1 = P['R0'](rm, 2)/(2*Am)
    # Laurent data of F at -ia: g(y) = y^2 F(-ia + y)
    from fp_integral import fp_segment
    FP, C2, C1 = fp_segment(P['F'], a)
    G = FP - Am**2*(I/(2*a))*(1/(2*I*a) + 2*c1 - I*om)
    S1, Sp1 = P['S'](mp.mpf(1)), P['S'](mp.mpf(1), 1)
    endpoint_identity = abs(-Am**2*(I/(2*a))*(1/(2*I*a) + 2*c1 - I*om) - S1*((2*c - 1)*S1 - 4*Sp1))
    Binc, Bref = amplitudes(P)
    S0 = (-1)**(l+1)*D_TS(P['lam'], c)/(16*om**4)*Bref/Binc
    val = G/(8*I*om**2*Binc*Bref)
    if verbose:
        print(f"l={l} c={mp.nstr(c,4)}: ODE residual {mp.nstr(res,3)}; C2/(A^2 i/2a) = {mp.nstr(C2/(Am**2*I/(2*a)),12)}; "
              f"C1 = {mp.nstr(C1,3)}; S0(M=0) = {mp.nstr(S0,15)}; endpoint identity {mp.nstr(endpoint_identity,3)}")
        print(f"   S_odd/eps (segment formula) = {mp.nstr(val, 15)}")
    return val

if __name__ == "__main__":
    mst = {'0.025': mp.mpc(0, '-8.81079077887e-9'), '0.05': mp.mpc(0, '-2.98193418097e-7'),
           '0.1': mp.mpc(0, '-1.0683752673e-5'), '0.2': mp.mpc(0, '-4.30179592639e-4')}
    for c in ('0.1', '0.2', '0.05'):
        v = S_odd(2, mp.mpf(c))
        print(f"   full MST (eps=1e-12):          {mp.nstr(mst[c], 12)}    ratio-1 = {mp.nstr(v/mst[c]-1, 3)}")
        sys.stdout.flush()
