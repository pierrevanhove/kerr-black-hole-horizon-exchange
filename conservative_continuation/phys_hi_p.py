# extend the physical-region l=2 coefficients to high powers of q (both branches, Richardson in r_eps with log handling)
import json, mpmath as mp
mp.mp.dps = 90
NE, RQ, NQ = 12, mp.mpf('0.3'), 40
def load(fn, RE):
    data = json.load(open(fn)); D = {(r['sign'], r['j']): r['pts'] for r in data}
    eps = [RE*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
    qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
    out = {}
    for sg in (1, -1):
        F = [[mp.mpc(mp.mpf(D[(sg, j)][k]['St'][0]), mp.mpf(D[(sg, j)][k]['St'][1])) for k in range(NQ)] for j in range(NE)]
        G = [mp.fsum(F[j][kk]*eps[j]**(-6) for j in range(NE))/NE for kk in range(NQ)]     # eps^6 coefficient on the q-circle
        for p in range(0, 25):
            out[(sg, p)] = mp.fsum(G[kk]*qs[kk]**(-p) for kk in range(NQ))/NQ
    return out
A1, A2 = load('../prl_checks/phys_l2_2e-5.json', mp.mpf('2e-5')), load('../prl_checks/phys_l2_1e-5.json', mp.mpf('1e-5'))
# model per (sign,p): f(r) = A + B ln r + C r^2 ; with two radii we can only remove the log if C is negligible -> use A = f(r2) - B ln r2, B = (f(r1)-f(r2))/ln 2
res = {}
for key in A1:
    B = (A1[key] - A2[key])/mp.log(2)
    res[key] = (A2[key] - B*mp.log(mp.mpf('1e-5')), B)
json.dump({f"{s}_{p}": [str(mp.re(a)), str(mp.im(a)), str(mp.re(b)), str(mp.im(b))] for (s, p), (a, b) in res.items()}, open('phys_e6_hi.json', 'w'))
mp.mp.dps = 25
for p in range(0, 25):
    e = (res[(1, p)][0] + res[(-1, p)][0])/2; d = (res[(1, p)][0] - res[(-1, p)][0])/2
    be = (res[(1, p)][1] + res[(-1, p)][1])/2
    print(f"p={p:2d}: even A6 = {mp.nstr(e, 16):>40}  odd A6 = {mp.nstr(d, 12):>32}  even log-coeff B6 = {mp.nstr(be, 10)}")
