import json, mpmath as mp
mp.mp.dps = 90
NE, RQ, NQ = 12, mp.mpf('0.3'), 40
def load(fn, RE):
    data = json.load(open(fn))
    D = {(r['sign'], r['j']): r['pts'] for r in data}
    eps = [RE*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
    qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
    out = {}
    for sg in (1, -1):
        F = [[mp.mpc(mp.mpf(D[(sg, j)][k]['St'][0]), mp.mpf(D[(sg, j)][k]['St'][1])) for k in range(NQ)] for j in range(NE)]
        for k in range(9):
            for p in range(11):
                out[(sg, k, p)] = mp.fsum(F[j][kk]*eps[j]**(-k)*qs[kk]**(-p) for j in range(NE) for kk in range(NQ))/(NE*NQ)
    return out
A1 = load('../prl_checks/phys_l2_2e-5.json', mp.mpf('2e-5'))
A2 = load('../prl_checks/phys_l2_1e-5.json', mp.mpf('1e-5'))
R = {key: (2**(8-key[1])*A2[key] - A1[key])/(2**(8-key[1]) - 1) if key[1] < 8 else A2[key] for key in A1}
print("Richardson change |R - A(1e-5)| for k=6, p=0..8:", [mp.nstr(abs(R[(1,6,p)] - A2[(1,6,p)]), 2) for p in range(9)])
json.dump({f"{s}_{k}_{p}": [str(mp.re(v)), str(mp.im(v))] for (s, k, p), v in R.items()}, open('l2_rich.json', 'w'))
# absorption at q^0 (Schwarzschild limit), order eps^6: Gamma = -(2 Re a6 + sum_{k=1}^{5} a_k conj(a_{6-k}))
for sg in (1, -1):
    a = {k: R[(sg, k, 0)] for k in range(7)}
    G6 = -(2*mp.re(a[6]) + mp.re(mp.fsum(a[k]*mp.conj(a[6-k]) for k in range(1, 6))))
    print(f"sign {sg:+d}: Gamma(q=0) at eps^6 = {mp.nstr(G6, 15)}   (Page/Teukolsky-Press l=2: 4/225 = {mp.nstr(mp.mpf(4)/225, 15)})")
print("\nduality-odd (A+ - A-)/2 at eps^6 [Re, Im]:")
for p in range(11):
    d = (R[(1, 6, p)] - R[(-1, 6, p)])/2
    print(f"  p={p:2d}: {mp.nstr(mp.re(d), 18):>26} {mp.nstr(mp.im(d), 3):>10}")
print("\nduality-even (A+ + A-)/2 at eps^6, Im part:")
for p in range(11):
    e = (R[(1, 6, p)] + R[(-1, 6, p)])/2
    print(f"  p={p:2d}: Im {mp.nstr(mp.im(e), 18):>26}   Re {mp.nstr(mp.re(e), 18):>26}")
