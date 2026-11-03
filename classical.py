"""Classical limit (eps -> 0 at fixed c = a omega) of the near-zone response, s = 0, l = 2, m = 2 and 1:
on the two super-extremal branches kappa = +-i sqrt(q^2-1) the exchange-odd/even ratio of the response is
X_odd/X_even = delta/u, with nu = l + delta and i tau = -m + u (structure of Supplement Eq. (S3)).
Also compares the scalar spheroidal eigenvalue with its O(c^2) formula."""
import mpmath as mp
mp.mp.dps = 40
from mst import MST, I

def E_exact(l, m, c, nmax=40):
    """Scalar (s=0) spheroidal eigenvalue for [(1-x^2)S']' + [c^2 x^2 + E - m^2/(1-x^2)]S = 0."""
    ls = [L for L in range(abs(m), abs(m)+2*nmax) if (L-l) % 2 == 0]
    n = len(ls); A = mp.matrix(n, n)
    for i, L in enumerate(ls):
        A[i, i] = L*(L+1) - c**2*mp.mpf(2*L*(L+1)-2*m*m-1)/((2*L-1)*(2*L+3))
        if i+1 < n:
            off = -c**2*mp.sqrt(mp.mpf(((L+1)**2-m*m)*((L+2)**2-m*m)))/((2*L+3)*mp.sqrt(mp.mpf((2*L+1)*(2*L+5))))
            A[i, i+1] = A[i+1, i] = off
    ev = mp.eigsy(A)[0]
    return min((ev[i] for i in range(n)), key=lambda e: abs(e-l*(l+1)))

aw = mp.mpf('0.3'); l = 2
for m in (2, 1):
    E = E_exact(l, m, aw)
    print(f"m={m}: E(c=0.3) = {mp.nstr(E,15)}  (O(c^2) formula: {mp.nstr(l*(l+1)-aw**2*(2*l*(l+1)-2*m*m-1)/((2*l-1)*(2*l+3)),15)})")
    for eps in ['0.04', '0.02', '0.01', '0.005']:
        e = mp.mpf(eps); q = 2*aw/e; kabs = mp.sqrt(q**2-1)
        lam = E - 2*m*aw + aw**2
        P = MST(e, q, m, 0, lam, I*kabs, N=60); P.solve_nu(l - 0.3*float(eps)**2); A = P.amplitudes()
        Q = MST(e, q, m, 0, lam, -I*kabs, N=60); Q.solve_nu(P.nu); B = Q.amplitudes()
        d = (P.nu-l).real; u = (I*P.tau + m).real
        r = (A['X']-B['X'])/(A['X']+B['X'])
        print(f"   eps={eps}: delta/eps^2={mp.nstr(d/e**2,8)}  u/eps^2={mp.nstr(u/e**2,8)}  X_odd/X_even={mp.nstr(r,8)}  delta/u={mp.nstr(d/u,8)}")
