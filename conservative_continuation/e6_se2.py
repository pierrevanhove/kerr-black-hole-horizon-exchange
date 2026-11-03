import json, mpmath as mp
from multiprocessing import Pool
from e6_se import job
if __name__ == '__main__':
    jobs = [(q, r) for q in ('50', '80', '120', '200') for r in ('1e-5', '5e-6', '2.5e-6')]
    with Pool(2) as p:
        res = p.map(job, jobs)
    json.dump(res, open('e6_se2.json', 'w'))
    mp.mp.dps = 40
    E = {(q, r): {int(k): mp.mpmathify(v.replace(' ', '')) for k, v in ev.items()} for q, r, ev, od in res}
    for qs in ('50', '80', '120', '200'):
        q = mp.mpf(qs)
        a, b, c = E[(qs, '1e-5')][6], E[(qs, '5e-6')][6], E[(qs, '2.5e-6')][6]
        r1 = (4*b - a)/3; r2 = (4*c - b)/3; rr = (8*r2 - r1)/7
        print(f"q={qs}: raw 32E6/q^5 at r=1e-5,5e-6,2.5e-6: {[mp.nstr(32*x/q**5, 12) for x in (a, b, c)]}")
        print(f"        Richardson r^2: {mp.nstr(32*r1/q**5, 12)}, {mp.nstr(32*r2/q**5, 12)};  r^2+r^3: {mp.nstr(32*rr/q**5, 12)}")
