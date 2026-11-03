"""Background driver: (1) physical-region 2D Cauchy for l=2, m in {1,-1,-2,0}, both branches, two radii;
                      (2) classical super-extremal s1(2,m;c) at 28 Chebyshev nodes, both eta."""
import json, time, mpmath as mp
from multiprocessing import Pool
NE, RQ, NQ = 12, '0.3', 24
def phys_task(args):
    m, sign, j, RE = args
    mp.mp.dps = 90
    from phys import spheroidal_c
    from mst import MST, I
    from csm import D_TS, secant
    eps = mp.mpf(RE)*mp.expjpi(2*(j + mp.mpf(1)/2)/NE)
    out = []
    for k in range(NQ):
        q = mp.mpf(RQ)*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ)
        c = q*eps/2
        E, b = spheroidal_c(-2, 2, m, c)
        lam = E - 2*m*c + c**2 - 2
        k0 = mp.sqrt(1 - q**2)
        kappa, logk = (k0, mp.log(k0)) if sign > 0 else (-k0, mp.log(k0) + I*mp.pi)
        P = MST(eps, q, m, -2, lam, kappa, logkappa=logk, N=20)
        x = secant(lambda x: P.g(2 + x*eps**2)/eps**2, -0.5095, mp.mpf(10)**(-45))
        P.nu = 2 + x*eps**2
        A = P.amplitudes()
        S = (-1)**3*D_TS(lam, c, m)/(16*(eps/2)**4)*A['S']*mp.exp(-2*I*eps*mp.log(eps))
        out.append([str(mp.re(S)), str(mp.im(S))])
    return dict(m=m, sign=sign, j=j, RE=RE, St=out)
def se_task(args):
    m, k = args
    mp.mp.dps = 110
    from csm import S_phys, s1, spheroidal, nu2_schw
    c = mp.mpf('0.2')*mp.cos(mp.pi*(k + mp.mpf(1)/2)/28)
    E, b = spheroidal(-2, 2, m, c, Lextra=24)
    out = dict(m=m, k=k, c=str(c))
    for eta in (1, -1):
        S, x, A = S_phys(2, m, -2, c, eta, '1e-16', nu2_schw(2), E=E, rtol=mp.mpf('1e-45'))
        v = s1(S, '1e-16'); out[f's1_{eta}'] = [str(mp.re(v)), str(mp.im(v))]
    return out
if __name__ == '__main__':
    t0 = time.time()
    with Pool(2) as p:
        se = p.map(se_task, [(m, k) for m in (1, 0, -1, -2) for k in range(28)])
        json.dump(se, open('se_allm.json', 'w')); print(f"SE done {time.time()-t0:.0f}s", flush=True)
        jobs = [(m, sg, j, RE) for m in (1, -1, -2, 0) for sg in (1, -1) for RE in ('2e-5', '1e-5') for j in range(NE)]
        res = []
        for i, r in enumerate(p.imap_unordered(phys_task, jobs)):
            res.append(r)
            if i % 24 == 0: print(f"phys {i+1}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
        json.dump(res, open('phys_allm.json', 'w')); print(f"all done {time.time()-t0:.0f}s", flush=True)
