"""Does the physical S-matrix contain the pole term (i/900) eps^4 Pi_m/nu2 of H^+ ?  Taylor coefficients in q of the
eps^4 and eps^5 coefficients (physical branch, m = 2), and the spheroidal eigenvalue coefficients E1, E2."""
import json, mpmath as mp
mp.mp.dps = 60
NE, NQ, RQ = 12, 40, mp.mpf('0.3')
qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
data = json.load(open('phys_l2_1e-5.json')); D = {(r['sign'], r['j']): r['pts'] for r in data}
eps = [mp.mpf('1e-5')*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
for n in (3, 4, 5):
    for sg in (1, -1):
        vals = [mp.fsum(mp.mpc(*map(mp.mpf, D[(sg, j)][k]['St']))*eps[j]**(-n) for j in range(NE))/NE for k in range(NQ)]
        tq = [mp.fsum(v*q**(-p) for v, q in zip(vals, qs))/NQ for p in range(8)]
        print(f"eps^{n} sign {sg:+d}: q^0..q^5 =", [mp.nstr(t, 8) for t in tq[:6]])
nu2 = -mp.mpf(107)/210
print("pole-term prediction (i/900) Pi_2/nu2: q^1 coeff =", mp.nstr(1j/900*8/nu2, 8), "  q^3 coeff =", mp.nstr(1j/900*24/nu2, 8))
print("superradiant eps^5 prediction Pi_2/900: q^1 =", mp.nstr(mp.mpf(8)/900, 8), " q^3 =", mp.nstr(mp.mpf(24)/900, 8))
from phys import spheroidal_c
for m in (2, 1, 0, -1):
    cs = [mp.mpf('1e-4')*k for k in (1, 2, 3)]
    Es = [mp.re(spheroidal_c(-2, 2, m, c)[0]) for c in cs]
    # E = 6 + E1 c + E2 c^2 + E3 c^3
    M = mp.matrix([[c, c**2, c**3] for c in cs]); sol = mp.lu_solve(M, mp.matrix([e - 6 for e in Es]))
    print(f"m={m:+d}: E1 = {mp.nstr(sol[0], 12)} (-4m/3 = {mp.nstr(-mp.mpf(4*m)/3, 12)})   E2 = {mp.nstr(sol[1], 10)}  ((10/189)(9-m^2)-1 = {mp.nstr(mp.mpf(10)/189*(9 - m**2) - 1, 10)})")
