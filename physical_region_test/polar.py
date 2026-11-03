"""
Sharp check: eta-odd part of the O(G) polar Compton amplitude for super-extremal Kerr,
s=-2, l=2, m=2 (only channel contributing to the eta-odd part through O(a^6)).
Prediction from the swap law (leading order):  H_odd = -(64/75) i (a w)^5 eps [1 + kappa6 a w + ...]
BGKV/Bautista polar result requires kappa6 - 2/9 = 2  (i.e. kappa6 = 20/9).
"""
import mpmath as mp
if __name__ == "__main__":
    mp.mp.dps = 60
from mst import MST, I

s, l, m = -2, 2, 2

def xmat(s, m, lmin, L):
    ls = list(range(lmin, L+1)); n = len(ls); X = mp.matrix(n, n)
    for i, li in enumerate(ls):
        X[i, i] = mp.mpf(-m*s)/(li*(li+1))
        if i+1 < n:
            lp = li+1
            X[i, i+1] = X[i+1, i] = mp.sqrt(mp.mpf((lp**2-m**2)*(lp**2-s**2)))/(lp*mp.sqrt(mp.mpf((2*li+1)*(2*li+3))))
    return ls, X

def E_spheroidal(s, l, m, c, L=40):
    lmin = max(abs(s), abs(m))
    ls, X = xmat(s, m, lmin, L)
    n = len(ls)
    H = mp.matrix(n, n)
    X2 = X*X
    for i in range(n):
        for j in range(n):
            H[i, j] = (ls[i]*(ls[i]+1) if i == j else 0) - c**2*X2[i, j] + 2*c*s*X[i, j]
    ev = mp.eigsy(H)[0]
    return min((ev[i] for i in range(n)), key=lambda e: abs(e - l*(l+1)))

def H_odd(aw, eps):
    aw, eps = mp.mpf(aw), mp.mpf(eps)
    q = 2*aw/eps; kabs = mp.sqrt(q**2-1)
    E = E_spheroidal(s, l, m, aw)
    lam = E - 2*m*aw + aw**2 - s*(s+1)
    out = []
    for k in (I*kabs, -I*kabs):
        P = MST(eps, q, m, s, lam, k, N=25)
        P.solve_nu(l - mp.mpf('0.5')*eps**2)
        out.append(P.amplitudes()['H'])
    return (out[0]-out[1])/2

if __name__ == "__main__":
    pred5 = -mp.mpf(64)/75*I
    rows = []
    for aw in ['0.004', '0.008', '0.012', '0.016', '0.02']:
        e1, e2 = mp.mpf('2e-5'), mp.mpf('1e-5')
        h1, h2 = H_odd(aw, e1)/e1, H_odd(aw, e2)/e2
        h0 = 2*h2 - h1                        # Richardson: remove O(eps) correction
        r = h0/mp.mpf(aw)**5
        rows.append((mp.mpf(aw), r))
        print(f"a*w={aw}:  H_odd/(eps (aw)^5) -> {mp.nstr(r, 14)}   (prediction at leading order: {mp.nstr(pred5, 10)})")

    # polynomial fit r(aw) = C5 + C6 aw + C7 aw^2 + C8 aw^3 + C9 aw^4 (exact interpolation through 5 points)
    A = mp.matrix([[x**k for k in range(5)] for x, _ in rows]); b = mp.matrix([r for _, r in rows])
    C = mp.lu_solve(A, b)
    print("\nC5 =", mp.nstr(C[0], 12), "  expected", mp.nstr(pred5, 12))
    print("kappa6 = C6/C5 =", mp.nstr(C[1]/C[0], 10), "  needed for BGKV polar formula: 20/9 =", mp.nstr(mp.mpf(20)/9, 10))
    print("C7/C5 =", mp.nstr(C[2]/C[0], 8))
