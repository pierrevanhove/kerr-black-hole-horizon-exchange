import json, mpmath as mp
mp.mp.dps = 40
fit = {k: mp.mpf(v) for k, v in json.load(open('e6_fit.json')).items()}
raw = json.load(open('phys_e6_hi.json'))
A = {tuple(int(t) for t in k.split('_')): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in raw.items()}
ev = [(A[(1, p)] + A[(-1, p)])/2 for p in range(25)]
def E6_closed(q):
    """physical-region closed form of the self-dual eps^6 coefficient, continued to q>1 along Im q > 0
       (kappa -> -i sqrt(q^2-1), log kappa -> log sqrt(q^2-1) - i pi/2)"""
    k = mp.sqrt(q**2 - 1); x = 2*q/k
    re_poly = mp.fsum(mp.re(ev[p])*q**p for p in range(5))
    im_poly = mp.im(ev[0]) + mp.im(ev[2])*q**2 + mp.im(ev[4])*q**4 + fit['a1']*q + fit['a3']*q**3 + fit['a5']*q**5
    F = mp.mpf(4)/225*q*(1 + 3*q**2)*((mp.log(k) - mp.j*mp.pi/2) + (mp.digamma(1 + x) + mp.digamma(1 - x))/2)
    return re_poly + mp.j*(im_poly + F)
res = []
for fn in ('e6_se.json', 'e6_se2.json', 'e6_se3.json'):
    res += json.load(open(fn))
E = {}
for q, r, evd, od in res:
    E.setdefault(q, {})[r] = {int(k): mp.mpmathify(v.replace(' ', '')) for k, v in evd.items()}
print("   q       E6 from MST at super-extremal spin (log-subtracted)        continued physical closed form                     rel.diff")
for qs in sorted(E, key=float):
    rs = sorted(E[qs], key=float)
    n = min(3, len(rs))
    M = mp.matrix([[1, mp.log(mp.mpf(r))] + ([mp.mpf(r)**2] if n == 3 else []) for r in rs[:n]])
    A6 = mp.lu_solve(M, mp.matrix([E[qs][r][6] for r in rs[:n]]))[0]
    cf = E6_closed(mp.mpf(qs))
    print(f"{qs:>6}  {mp.nstr(A6, 18):>50}  {mp.nstr(cf, 18):>50}  {mp.nstr(abs(A6 - cf)/abs(cf), 2)}")
