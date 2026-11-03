"""Self-dual (duality-even) eps^6 coefficient of the l=2, m=2 S-matrix at super-extremal q (fixed q, expansion in eps)."""
import json, mpmath as mp
from multiprocessing import Pool
from d6_se import S_tilde_k
NE = 16
def job(args):
    qs, rs = args
    mp.mp.dps = 90
    from mst import I
    q, r = mp.mpf(qs), mp.mpf(rs)
    kq = mp.sqrt(q**2 - 1)
    ka, la = -I*kq, mp.log(kq) - I*mp.pi/2
    kb, lb = I*kq, mp.log(kq) + I*mp.pi/2
    vals = []
    for j in range(NE):
        e = r*mp.expjpi(2*(j + mp.mpf(1)/2)/NE)
        Sa, Sb = S_tilde_k(q, e, ka, la), S_tilde_k(q, e, kb, lb)
        vals.append((e, (Sa + Sb)/2, (Sa - Sb)/2))
    ev = {k: mp.fsum(s*e**(-k) for e, s, d in vals)/NE for k in range(0, 9)}
    od = {k: mp.fsum(d*e**(-k) for e, s, d in vals)/NE for k in range(0, 9)}
    return qs, rs, {k: str(v) for k, v in ev.items()}, {k: str(v) for k, v in od.items()}
if __name__ == '__main__':
    qlist = ['1.5', '2', '3', '5', '8', '12', '20', '30', '50', '80']
    jobs = [(q, r) for q in qlist for r in ('2e-5', '1e-5')]
    with Pool(2) as p:
        res = p.map(job, jobs)
    json.dump(res, open('e6_se.json', 'w'))
    print("done")
