"""Direct MST check of the analytic continuation of D6(q) to super-extremal q (fixed q, expansion in eps)."""
import sys, mpmath as mp
from multiprocessing import Pool
mp.mp.dps = 90
NE = 16

def S_tilde_k(q, eps, kappa, logk, l=2, m=2, s=-2):
    from mst import MST, I
    from phys import spheroidal_c
    from csm import D_TS, secant
    c = q*eps/2
    E, b = spheroidal_c(s, l, m, c)
    lam = E - 2*m*c + c**2 - s*(s+1)
    P = MST(eps, q, m, s, lam, kappa, logkappa=logk, N=20)
    x = secant(lambda x: P.g(l + x*eps**2)/eps**2, -0.5095, mp.mpf(10)**(-45))
    P.nu = l + x*eps**2
    A = P.amplitudes()
    S = (-1)**(l+1)*D_TS(lam, c, m)/(16*(eps/2)**4)*A['S']
    return S*mp.exp(-2*I*eps*mp.log(eps))

def A6(args):
    qs, rs = args
    mp.mp.dps = 90
    from mst import I
    q, r = mp.mpf(qs), mp.mpf(rs)
    kq = mp.sqrt(q**2 - 1)
    ka, la = -I*kq, mp.log(kq) - I*mp.pi/2       # continuation of the physical branch along Im q > 0
    kb, lb = I*kq, mp.log(kq) + I*mp.pi/2        # its dual (log kappa shifted by +i pi)
    out = {}
    for k in (5, 6):
        out[k] = 0
    vals = []
    for j in range(NE):
        e = r*mp.expjpi(2*(j + mp.mpf(1)/2)/NE)
        d = (S_tilde_k(q, e, ka, la) - S_tilde_k(q, e, kb, lb))/2
        vals.append((e, d))
    res = {k: mp.fsum(d*e**(-k) for e, d in vals)/NE for k in range(0, 9)}
    return qs, rs, {k: str(v) for k, v in res.items()}

if __name__ == '__main__':
    qlist = ['1.5', '2', '3', '5', '8', '12']
    jobs = [(q, r) for q in qlist for r in ('2e-5', '1e-5')]
    with Pool(2) as p:
        res = p.map(A6, jobs)
    import json
    json.dump(res, open('d6_se.json', 'w'))
    print("done")
