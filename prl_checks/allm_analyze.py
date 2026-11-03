"""All-m analysis of the l=2 eps^6 physical-region coefficients (self-dual / anti-self-dual split) and
comparison of the continued closed form with the super-extremal classical eta-even c^5 coefficient."""
import json, mpmath as mp
mp.mp.dps = 90
NE, NQ, RQ = 12, 24, mp.mpf('0.3')
data = json.load(open('phys_allm.json'))
D = {(r['m'], r['sign'], r['RE'], r['j']): r['St'] for r in data}
qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
def coef(m, sg, RE, n):
    r = mp.mpf(RE)
    eps = [r*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
    F = [[mp.mpc(mp.mpf(D[(m, sg, RE, j)][k][0]), mp.mpf(D[(m, sg, RE, j)][k][1])) for k in range(NQ)] for j in range(NE)]
    G = [mp.fsum(F[j][k]*eps[j]**(-n) for j in range(NE))/NE for k in range(NQ)]
    return [mp.fsum(G[k]*qs[k]**(-p) for k in range(NQ))/NQ for p in range(NQ)]
def Pm(m, q):
    kap = mp.sqrt(1 - q**2); t0 = -m*q/kap
    return kap**5*t0*(1 + t0**2)*(4 + t0**2)
def RePsi(m, q):
    kap = mp.sqrt(1 - q**2); t0 = -m*q/kap
    return (mp.digamma(1 + 1j*t0) + mp.digamma(1 - 1j*t0))/2
def tay(f, P=24, N=128, r=mp.mpf('0.35')):
    zs = [r*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
    v = [f(z) for z in zs]
    return [mp.fsum(vv*z**(-p) for vv, z in zip(v, zs))/N for p in range(P)]
res = {}
for m in (1, -1, -2, 0):
    A = {}
    for sg in (1, -1):
        c2, c1 = coef(m, sg, '2e-5', 6), coef(m, sg, '1e-5', 6)
        B = [(a - b)/mp.log(2) for a, b in zip(c2, c1)]
        A0 = [b - BB*mp.log(mp.mpf('1e-5')) for b, BB in zip(c1, B)]
        A[sg] = (A0, B)
    sd = [(A[1][0][p] + A[-1][0][p])/2 for p in range(NQ)]
    asd = [(A[1][0][p] - A[-1][0][p])/2 for p in range(NQ)]
    Bsd = [(A[1][1][p] + A[-1][1][p])/2 for p in range(NQ)]
    Basd = [(A[1][1][p] - A[-1][1][p])/2 for p in range(NQ)]
    res[m] = dict(sd=sd, asd=asd, Bsd=Bsd, Basd=Basd)
    Pt = tay(lambda q: Pm(m, q)) if m != 0 else [0]*24
    print(f"\n===== m = {m:+d} =====")
    print(" log coefficient B6 (self-dual) vs -(i/450) P^(m):")
    for p in range(0, 8):
        print(f"   p={p}: {mp.nstr(Bsd[p], 14):>40}   pred {mp.nstr(-1j*Pt[p]/450, 14) if m else 0}   (anti-self-dual log: {mp.nstr(Basd[p], 3)})")
    if m != 0:
        Lt = tay(lambda q: Pm(m, q)*mp.log(mp.sqrt(1 - q**2)))
        Pst = tay(lambda q: Pm(m, q)*RePsi(m, q))
        odd = list(range(1, 16, 2))
        d = {p: mp.im(sd[p]) for p in range(NQ)}
        rows = [[1 if p == 1 else 0, 1 if p == 3 else 0, 1 if p == 5 else 0, mp.re(Lt[p]), mp.re(Pst[p])] for p in odd]
        sol = mp.qr_solve(mp.matrix(rows), mp.matrix([d[p] for p in odd]))[0]
        print(f" free fit (odd p<=15): a1={mp.nstr(sol[0],15)} a3={mp.nstr(sol[1],15)} a5={mp.nstr(sol[2],15)}  "
              f"450*lam_log={mp.nstr(450*sol[3],12)}  450*lam_psi={mp.nstr(450*sol[4],12)}")
        lam = -mp.mpf(1)/450
        a = {p: d[p] - lam*(mp.re(Lt[p]) + mp.re(Pst[p])) for p in (1, 3, 5)}
        print(f" fixed lam=-1/450: a1={mp.nstr(a[1],15)} a3={mp.nstr(a[3],15)} a5={mp.nstr(a[5],15)}")
        print("   p   Im(self-dual)              closed-form (p>=7)        residual     Re(self-dual)")
        for p in range(0, 16):
            model = lam*(mp.re(Lt[p]) + mp.re(Pst[p])) + (a[p] if p in a else 0)
            print(f"  {p:2d}  {mp.nstr(d[p], 16):>24}  {mp.nstr(model, 16):>24}  {mp.nstr(d[p]-model, 3):>10}   {mp.nstr(mp.re(sd[p]), 6)}")
        res[m]['a5'] = a[5]; res[m]['a'] = a
    else:
        print("   p   self-dual coefficient (should terminate: polynomial)")
        for p in range(0, 16):
            print(f"  {p:2d}  {mp.nstr(sd[p], 16)}")
        res[m]['a5'] = mp.im(sd[5])
json.dump({str(m): {k: ([[str(mp.re(x)), str(mp.im(x))] for x in v] if isinstance(v, list) else
                         (str(v) if not isinstance(v, dict) else {str(kk): str(vv) for kk, vv in v.items()}))
                    for k, v in res[m].items()} for m in res}, open('allm_phys_result.json', 'w'))
