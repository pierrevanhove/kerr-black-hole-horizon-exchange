"""Uniqueness test for the exchange-even eps^6 closed form (m = 2): the polynomial part is fixed from the |q| = 0.3
Taylor data (joint log bookkeeping), then the closed form is continued to 1.5 <= q <= 5000 and compared with direct
super-extremal MST.  A non-polynomial remainder invisible in the Taylor data would show up here."""
import json, mpmath as mp
from collections import defaultdict
mp.mp.dps = 50
NE, NQ, m = 12, 40, 2
th = [(mp.pi*(2*j + 1)/NE if mp.pi*(2*j + 1)/NE <= mp.pi else mp.pi*(2*j + 1)/NE - 2*mp.pi) for j in range(NE)]
W = {d: mp.fsum(t*mp.expj(d*t) for t in th)/NE for d in range(-4, 5)}
qs = [mp.mpf('0.3')*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
F = {}
for RE in ('2e-5', '1e-5'):
    data = json.load(open(f'phys_l2_{RE}.json')); D = {(r['sign'], r['j']): r['pts'] for r in data}
    F[RE] = {sg: [[mp.mpc(*map(mp.mpf, D[(sg, j)][k]['St'])) for k in range(NQ)] for j in range(NE)] for sg in (1, -1)}
A6 = {}
for sg in (1, -1):
    vals = []
    for k in range(NQ):
        rows, rhs = [], []
        for RE in ('2e-5', '1e-5'):
            r = mp.mpf(RE); eps = [r*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
            for n in (6, 7, 8):
                rows.append(sum(([1 if kk == n else 0, (mp.log(r) if kk == n else r**(kk - n)*1j*W[kk - n])] for kk in (6, 7, 8)), []))
                rhs.append(mp.fsum(F[RE][sg][j][k]*eps[j]**(-n) for j in range(NE))/NE)
        vals.append(mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))[0])
    A6[sg] = vals
E6 = [(a + b)/2 for a, b in zip(A6[1], A6[-1])]
def T6(q, k, lk):
    w = m*q/k; Pi = m*q*(1 + (m**2 - 1)*q**2)*(4 + (m**2 - 4)*q**2)
    return 1j/450*Pi*(lk + (mp.digamma(1 + 1j*w) + mp.digamma(1 - 1j*w))/2)
R = [e - T6(q, mp.sqrt(1 - q**2), mp.log(mp.sqrt(1 - q**2))) for e, q in zip(E6, qs)]
Rt = [mp.fsum(x*q**(-p) for x, q in zip(R, qs))/NQ for p in range(NQ)]
print("remainder Taylor coefficients |p| for p = 6..20 (should vanish):", mp.nstr(max(abs(Rt[p]) for p in range(6, 21)), 3))
poly = Rt[:6]
def E6_cont(q):
    kq = mp.sqrt(q**2 - 1); k = -1j*kq
    return mp.fsum(c*q**p for p, c in enumerate(poly)) + T6(q, k, mp.log(kq) - 1j*mp.pi/2)
res = []
for fn in ('e6_se.json', 'e6_se2.json', 'e6_se3.json'):
    res += json.load(open(fn))
E = defaultdict(dict)
for q, r, evd, od in res:
    E[q][r] = {int(k): mp.mpmathify(v) for k, v in evd.items()}
print("   q     rel. difference (continued closed form vs direct super-extremal MST)")
rd = {}
for qs_ in sorted(E, key=float):
    rs = sorted(E[qs_], key=float); n = min(3, len(rs))
    M = mp.matrix([[1, mp.log(mp.mpf(r))] + ([mp.mpf(r)] if n == 3 else []) for r in rs[:n]])
    a6 = mp.lu_solve(M, mp.matrix([E[qs_][r][6] for r in rs[:n]]))[0]
    cf = E6_cont(mp.mpf(qs_))
    rd[float(qs_)] = abs(a6 - cf)/abs(cf)
    print(f"{qs_:>6}   {mp.nstr(abs(a6 - cf)/abs(cf), 2)}")

# an extra term c_p q^p (p >= 6) would grow faster than the q^5 behaviour of the continued function:
# bound from the agreement at the largest q (Supplement S5 (c))
qmax = max(rd); qM = mp.mpf(qmax)
print(f"\nbounds on an additional term c_p q^p from the agreement at q = {qmax:g} (|E6(q)|/q^5 = {mp.nstr(abs(E6_cont(qM))/qM**5, 5)}):")
for p in (6, 7, 8):
    print(f"   |c_{p}| <~ {mp.nstr(rd[qmax]*abs(E6_cont(qM))/qM**p, 3)}")

# compare with the polynomial from the naive two-radius extraction (phys_e6_hi.json + e6_fit.json)
fit = {k: mp.mpf(v) for k, v in json.load(open('e6_fit.json')).items()}
raw = json.load(open('phys_e6_hi.json'))
Aold = {tuple(int(t) for t in k.split('_')): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in raw.items()}
ev = [(Aold[(1, p)] + Aold[(-1, p)])/2 for p in range(25)]
Tt = [mp.fsum(T6(q, mp.sqrt(1 - q**2), mp.log(mp.sqrt(1 - q**2)))*q**(-p) for q in qs)/NQ for p in range(6)]
print("\n p   joint-bookkeeping poly coeff                 naive poly coeff (ev - T6)                 difference")
for p in range(6):
    old = ev[p] - Tt[p]
    print(f" {p}   {mp.nstr(poly[p], 14):>40}   {mp.nstr(old, 14):>40}   {mp.nstr(abs(poly[p] - old), 3)}")

# apply the same joint log bookkeeping to the two-radius super-extremal data (e6_se.json, NE = 16 nodes)
NE2 = 16
th2 = [(mp.pi*(2*j + 1)/NE2 if mp.pi*(2*j + 1)/NE2 <= mp.pi else mp.pi*(2*j + 1)/NE2 - 2*mp.pi) for j in range(NE2)]
W2 = {d: mp.fsum(t*mp.expj(d*t) for t in th2)/NE2 for d in range(-4, 5)}
print("\n   q     rel. difference with joint bookkeeping on both sides")
for qs_ in sorted(E, key=float):
    rs = sorted(E[qs_], key=float)
    if len(rs) != 2: continue
    rows, rhs = [], []
    for r in rs:
        rr = mp.mpf(r)
        for n in (6, 7, 8):
            rows.append(sum(([1 if kk == n else 0, (mp.log(rr) if kk == n else rr**(kk - n)*1j*W2[kk - n])] for kk in (6, 7, 8)), []))
            rhs.append(E[qs_][r][n])
    a6 = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))[0]
    cf = E6_cont(mp.mpf(qs_))
    print(f"{qs_:>6}   {mp.nstr(abs(a6 - cf)/abs(cf), 2)}")
