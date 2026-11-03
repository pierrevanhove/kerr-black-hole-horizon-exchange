import json, time, mpmath as mp
from multiprocessing import Pool
DPS, RE, NE, RQ, NQ = 90, '2e-5', 12, '0.35', 28
JOBS = [(2, 1), (2, -1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 1)]

def task(args):
    l, sign, j = args
    mp.mp.dps = DPS
    from phys import S_tilde
    eps = mp.mpf(RE)*mp.expjpi(2*(j + mp.mpf(1)/2)/NE)
    out = []
    for k in range(NQ):
        q = mp.mpf(RQ)*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ)
        St, b, x, gr = S_tilde(l, 2, -2, q, eps, sign)
        out.append(dict(k=k, St=[str(mp.re(St)), str(mp.im(St))],
                        b={str(L): [str(mp.re(v)), str(mp.im(v))] for L, v in b.items() if L <= 9}, g=str(gr)))
    return dict(l=l, sign=sign, j=j, pts=out)

if __name__ == '__main__':
    jobs = [(l, sg, j) for (l, sg) in JOBS for j in range(NE)]
    res, t0 = [], time.time()
    with Pool(2) as p:
        for i, r in enumerate(p.imap_unordered(task, jobs)):
            res.append(r)
            print(f"{i+1}/{len(jobs)} l={r['l']} sign={r['sign']} j={r['j']}  {time.time()-t0:.0f}s", flush=True)
    json.dump(res, open('phys_run.json', 'w'))
    print("done", flush=True)
