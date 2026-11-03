import sys, json, time, mpmath as mp
from multiprocessing import Pool
DPS, EPS, R, NC, LMAX = 110, '1e-16', '0.2', 28, 19

def task(args):
    l, k = args
    mp.mp.dps = DPS
    from csm import S_phys, s1, spheroidal, nu2_schw
    c = mp.mpf(R)*mp.cos(mp.pi*(k + mp.mpf(1)/2)/NC)
    t = time.time()
    E, b = spheroidal(-2, l, 2, c, Lextra=24)
    out = dict(l=l, k=k, c=str(c), E=str(E), b={str(L): str(v) for L, v in b.items()})
    for eta in (1, -1):
        S, x, A = S_phys(l, 2, -2, c, eta, EPS, nu2_schw(l), E=E, rtol=mp.mpf('1e-45'))
        v = s1(S, EPS)
        out[f's1_{eta}'] = [str(mp.re(v)), str(mp.im(v))]
        out[f'x_{eta}'] = [str(mp.re(x)), str(mp.im(x))]; out[f'g_{eta}'] = [str(A['gres']), str(A['gscale'])]
    out['time'] = time.time() - t
    return out

if __name__ == '__main__':
    jobs = [(l, k) for l in range(2, LMAX + 1) for k in range(NC)]
    res = []
    with Pool(2) as p:
        for i, r in enumerate(p.imap_unordered(task, jobs)):
            res.append(r)
            if i % 20 == 0:
                print(f"{i+1}/{len(jobs)} l={r['l']} k={r['k']} t={r['time']:.1f}s", flush=True)
                json.dump(res, open('polar_run2.json', 'w'))
    json.dump(res, open('polar_run2.json', 'w'))
    print("done", flush=True)
