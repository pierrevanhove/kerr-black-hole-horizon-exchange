import json, mpmath as mp
from multiprocessing import Pool
from e6_se import job
if __name__ == '__main__':
    jobs = []
    for qs in ('300', '600', '1000', '2000', '5000'):
        for c0 in ('2e-3', '1e-3', '5e-4'):
            jobs.append((qs, str(mp.mpf(c0)/mp.mpf(qs))))
    with Pool(2) as p:
        res = p.map(job, jobs)
    json.dump(res, open('e6_se3.json', 'w'))
    print("done")
