"""Self-dual eps^7 closed-form test for m = 1, -1, -2, 0 (phys_allm.json, NQ=24):
A7 = (1/kappa) d/dtau0 T6 + (u0 + u1 q) T6 + poly6,  B7 = -(i/450)[P'/kappa + (u0 + u1 q) P]  (expect u1 = 5m/9)."""
import json, mpmath as mp
mp.mp.dps = 60
NE, NQ, RQ = 12, 24, mp.mpf('0.3')
qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
th = [(mp.pi*(2*j + 1)/NE if mp.pi*(2*j + 1)/NE <= mp.pi else mp.pi*(2*j + 1)/NE - 2*mp.pi) for j in range(NE)]
W = {d: mp.fsum(t*mp.expj(d*t) for t in th)/NE for d in range(-4, 5)}
data = json.load(open('phys_allm.json')); D = {(r['m'], r['sign'], r['RE'], r['j']): r['St'] for r in data}
u0x = mp.pi + 1j*(2*mp.euler + 2*mp.log(2) - mp.mpf(37)/12 + mp.mpf(107)/210)
kap = lambda q: mp.sqrt(1 - q**2)
Nn, rr = 128, mp.mpf('0.3')
zs = [rr*mp.expjpi(2*(mp.mpf(k) + mp.mpf(1)/3)/Nn) for k in range(Nn)]
def tay(f, P):
    v = [f(z) for z in zs]; return [mp.fsum(vv*z**(-p) for vv, z in zip(v, zs))/Nn for p in range(P)]
for m in (1, -1, -2, 0):
    C = {}
    for RE in ('2e-5', '1e-5'):
        eps = [mp.mpf(RE)*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
        for k in range(NQ):
            F = [(mp.mpc(*map(mp.mpf, D[(m, 1, RE, j)][k])) + mp.mpc(*map(mp.mpf, D[(m, -1, RE, j)][k])))/2 for j in range(NE)]
            for n in (6, 7, 8):
                C[(RE, n, k)] = mp.fsum(F[j]*eps[j]**(-n) for j in range(NE))/NE
    A7v, B7v = [], []
    for k in range(NQ):
        rows, rhs = [], []
        for RE in ('2e-5', '1e-5'):
            r = mp.mpf(RE)
            for n in (6, 7, 8):
                row = []
                for kk in (6, 7, 8):
                    row += [1 if kk == n else 0, (mp.log(r) if kk == n else r**(kk - n)*1j*W[kk - n])]
                rows.append(row); rhs.append(C[(RE, n, k)])
        sol = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
        A7v.append(sol[2]); B7v.append(sol[3])
    tq = lambda vals: [mp.fsum(v*q**(-p) for v, q in zip(vals, qs))/NQ for p in range(NQ)]
    A7, B7 = tq(A7v), tq(B7v)
    def Pf(q):
        k = kap(q); t = -m*q/k; return k**5*t*(1 + t**2)*(4 + t**2)
    def dPk(q):
        k = kap(q); t = -m*q/k; return k**4*((1 + t**2)*(4 + t**2) + 2*t**2*(4 + t**2) + 2*t**2*(1 + t**2))
    def Rp(q):
        k = kap(q); t = -m*q/k; return (mp.digamma(1 + 1j*t) + mp.digamma(1 - 1j*t))/2
    def dRp(q):
        k = kap(q); t = -m*q/k; return 1j*(mp.psi(1, 1 + 1j*t) - mp.psi(1, 1 - 1j*t))/2
    P_ = 14
    T6 = tay(lambda q: -1j/450*Pf(q)*(mp.log(kap(q)) + Rp(q)), P_)
    dT6 = tay(lambda q: -1j/450*(dPk(q)*(mp.log(kap(q)) + Rp(q)) + Pf(q)*dRp(q)/kap(q)), P_)
    Pt, dPt = tay(Pf, P_), tay(dPk, P_)
    cols = ([T6, [T6[p-1] if p else 0 for p in range(P_)]] if m != 0 else []) + [[1 if p == j else 0 for p in range(P_)] for j in range(7)]
    M = mp.matrix([[c[p] for c in cols] for p in range(P_)])
    sol, res = mp.qr_solve(M, mp.matrix([A7[p] - dT6[p] for p in range(P_)]))
    predB = [-1j/450*(dPt[p] + u0x*Pt[p] + mp.mpf(5*m)/9*(Pt[p-1] if p else 0)) for p in range(P_)]
    errB = max(abs(B7[p] - predB[p]) for p in range(P_))
    if m != 0:
        print(f"m={m:+d}: u0 = {mp.nstr(sol[0], 11)} (closed {mp.nstr(u0x, 11)})  u1 = {mp.nstr(sol[1], 11)} (5m/9 = {mp.nstr(mp.mpf(5*m)/9, 11)})"
              f"  resid {mp.nstr(res, 2)} | max|B7 - pred| = {mp.nstr(errB, 2)}  (|B7| ~ {mp.nstr(max(abs(x) for x in B7[:6]), 3)})")
    else:
        print(f"m=+0: A7 - (1/kappa) dT6/dtau0 = poly6: resid {mp.nstr(res, 2)};  max|B7 - pred| = {mp.nstr(errB, 2)};  B7(p=0,2,4) = "
              f"{[mp.nstr(B7[p], 10) for p in (0, 2, 4)]}  vs -(2i/225) kappa^4 -> {[mp.nstr(x, 10) for x in tay(lambda q: -2j/225*kap(q)**4, 5)[0:5:2]]}")
