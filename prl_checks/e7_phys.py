"""Self-dual eps^7 function for l=2, m=2 in the physical region.  Log bookkeeping on the discrete eps-circle:
c_n(r) = A_n + B_n ln r + sum_{k!=n} B_k r^(k-n) i W_(k-n)   (principal ln eps on the nodes),
solved jointly for (A_n, B_n), n = 6,7,8 from two radii.  Then test
E7 = (1/kappa) d/dtau0 T6 + (u0 + u1 q) T6 + poly(q),   T6 = -(i/450) P(tau0) [ln kappa + Re psi(1 + i tau0)],
B7 = -(i/450) [P'(tau0)/kappa + (b0 + b1 q) P(tau0)] + poly(q)."""
import json, mpmath as mp
mp.mp.dps = 60
NE, NQ, RQ, m = 12, 40, mp.mpf('0.3'), 2
qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
th = [(mp.pi*(2*j + 1)/NE if mp.pi*(2*j + 1)/NE <= mp.pi else mp.pi*(2*j + 1)/NE - 2*mp.pi) for j in range(NE)]
W = {d: mp.fsum(t*mp.expj(d*t) for t in th)/NE for d in range(-4, 5)}
radii = ('2e-5', '1e-5')
C = {}
for RE in radii:
    data = json.load(open(f'phys_l2_{RE}.json')); D = {(r['sign'], r['j']): r['pts'] for r in data}
    eps = [mp.mpf(RE)*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
    for k in range(NQ):
        F = [(mp.mpc(*map(mp.mpf, D[(1, j)][k]['St'])) + mp.mpc(*map(mp.mpf, D[(-1, j)][k]['St'])))/2 for j in range(NE)]
        for n in (6, 7, 8):
            C[(RE, n, k)] = mp.fsum(F[j]*eps[j]**(-n) for j in range(NE))/NE
ns = (6, 7, 8)
A = {n: [] for n in ns}; B = {n: [] for n in ns}
for k in range(NQ):
    rows, rhs = [], []
    for RE in radii:
        r = mp.mpf(RE)
        for n in ns:
            row = []
            for kk in ns:       # unknown order: A6,B6,A7,B7,A8,B8
                row += [1 if kk == n else 0, (mp.log(r) if kk == n else r**(kk - n)*1j*W[kk - n])]
            rows.append(row); rhs.append(C[(RE, n, k)])
    sol = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
    for i, n in enumerate(ns):
        A[n].append(sol[2*i]); B[n].append(sol[2*i + 1])
tq = lambda vals: [mp.fsum(v*q**(-p) for v, q in zip(vals, qs))/NQ for p in range(NQ)]
A6, B6, A7, B7 = tq(A[6]), tq(B[6]), tq(A[7]), tq(B[7])
print("B6 check (p=1,3): ", mp.nstr(B6[1], 15), mp.nstr(B6[3], 15), "  expected 4i/225 =", mp.nstr(mp.mpf(4)/225, 15), " 12i/225 =", mp.nstr(mp.mpf(12)/225, 15))
kap = lambda q: mp.sqrt(1 - q**2)
def Pf(q):
    k = kap(q); t = -m*q/k; return k**5*t*(1 + t**2)*(4 + t**2)
def dPf(q):         # d/dtau0 of P at fixed kappa
    k = kap(q); t = -m*q/k; return k**5*(1 + 6*t**2 + 3*0 + 5*t**4 - 0)*1 if False else k**5*((1 + t**2)*(4 + t**2) + 2*t**2*(4 + t**2) + 2*t**2*(1 + t**2))
def Repsi(q):
    k = kap(q); t = -m*q/k; return (mp.digamma(1 + 1j*t) + mp.digamma(1 - 1j*t))/2
def dRepsi(q):      # d/dtau0 Re psi(1 + i tau0)
    k = kap(q); t = -m*q/k; return 1j*(mp.psi(1, 1 + 1j*t) - mp.psi(1, 1 - 1j*t))/2
T6f = lambda q: -1j/450*Pf(q)*(mp.log(kap(q)) + Repsi(q))
dT6f = lambda q: -1j/450*(dPf(q)*(mp.log(kap(q)) + Repsi(q)) + Pf(q)*dRepsi(q))/kap(q)
Nn, rr = 160, mp.mpf('0.3')
zs = [rr*mp.expjpi(2*(mp.mpf(k) + mp.mpf(1)/3)/Nn) for k in range(Nn)]
def tay(f, P=40):
    v = [f(z) for z in zs]; return [mp.fsum(vv*z**(-p) for vv, z in zip(v, zs))/Nn for p in range(P)]
T6, dT6, Pt, dPt = tay(T6f), tay(dT6f), tay(Pf), tay(lambda q: dPf(q)/kap(q))
P_ = 24
def fit(target, base, extra_cols, deg, labels):
    cols = list(extra_cols) + [[1 if p == k else 0 for p in range(P_)] for k in range(deg + 1)]
    M = mp.matrix([[c[p] for c in cols] for p in range(P_)])
    sol, res = mp.qr_solve(M, mp.matrix([target[p] - base[p] for p in range(P_)]))
    return sol, res
shift = lambda arr, j: [arr[p - j] if p >= j else 0 for p in range(P_)]
sol, res = fit(A7, dT6, [T6[:P_], shift(T6, 1), shift(T6, 2)], 8, None)
print("\nA7 fit:  u0 =", mp.nstr(sol[0], 12), "  u1 =", mp.nstr(sol[1], 12), "  u2 =", mp.nstr(sol[2], 4), "  residual", mp.nstr(res, 3))
print("        poly coeffs:", [mp.nstr(s, 6) for s in sol[3:]])
u0x = mp.pi + 1j*(2*mp.euler + 2*mp.log(2) - mp.mpf(37)/12 + mp.mpf(107)/210)
print("u0 closed form pi + i(2 gamma_E + 2 ln2 - 37/12 + 107/210) =", mp.nstr(u0x, 14), "  fit - closed:", mp.nstr(sol[0] - u0x, 3))
# refit with u0, u1 fixed to closed forms -> polynomial part
fixed = [dT6[p] + u0x*T6[p] + mp.mpf(5*m)/9*(T6[p-1] if p >= 1 else 0) for p in range(P_)]
solp, resp = fit(A7, fixed, [], 6, None)
print("A7 with u0,u1 fixed: residual", mp.nstr(resp, 3), " poly:", [mp.nstr(s, 12) for s in solp])
print("64*[q^6]poly =", mp.nstr(64*solp[6], 12), "  classical eta=+1 c^6 (m=2) = 0.000185311568809566i")
predB7 = [-1j/450*(dPt[p] + u0x*Pt[p] + mp.mpf(5*m)/9*(Pt[p-1] if p >= 1 else 0)) for p in range(P_)]
print("B7 data vs -(i/450)[P'/kappa + (u0 + 5m q/9) P]:")
for p in range(0, 10):
    print(f"   p={p}: {mp.nstr(B7[p], 12):>36}   {mp.nstr(predB7[p], 12):>36}   diff {mp.nstr(B7[p]-predB7[p], 4)}")
json.dump({'poly': [str(s) for s in solp], 'u0': str(u0x)}, open('e7_fit.json', 'w'))
print("(5/9) m =", mp.nstr(mp.mpf(5*m)/9, 12), "  anti-self-dual u0 = pi - 0.0330838328901 i")
json.dump({'A7': [str(x) for x in A7], 'B7': [str(x) for x in B7]}, open('e7_phys.json', 'w'))
