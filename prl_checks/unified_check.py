"""Physical branch only (kappa = +sqrt(1-q^2)): A6 - [(i/450) Pi (ln kappa + psi(1 + i w)) - kappa^5 (1+w^2)(4+w^2)/900]
must be a polynomial of degree <= 5 in q;  B6 = (i/450) Pi.   w = m q/kappa = m Omega_H/(2 pi T_H)."""
import json, mpmath as mp
mp.mp.dps = 50
NE = 12
def eps_nodes(RE): return [mp.mpf(RE)*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
def branch_coeffs(m):
    if m == 2:
        NQ = 40; out = {}
        for RE in ('2e-5', '1e-5'):
            data = json.load(open(f'phys_l2_{RE}.json')); D = {(r['sign'], r['j']): r['pts'] for r in data}
            e = eps_nodes(RE)
            out[RE] = [mp.fsum(mp.mpc(*map(mp.mpf, D[(1, j)][k]['St']))*e[j]**(-6) for j in range(NE))/NE for k in range(NQ)]
    else:
        NQ = 24; out = {}
        data = json.load(open('phys_allm.json')); D = {(r['m'], r['sign'], r['RE'], r['j']): r['St'] for r in data}
        for RE in ('2e-5', '1e-5'):
            e = eps_nodes(RE)
            out[RE] = [mp.fsum(mp.mpc(*map(mp.mpf, D[(m, 1, RE, j)][k]))*e[j]**(-6) for j in range(NE))/NE for k in range(NQ)]
    return NQ, out
for m in (2, 1, -1, -2, 0):
    NQ, C = branch_coeffs(m)
    qs = [mp.mpf('0.3')*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
    B = [(a - b)/mp.log(2) for a, b in zip(C['2e-5'], C['1e-5'])]
    A = [b - BB*mp.log(mp.mpf('1e-5')) for b, BB in zip(C['1e-5'], B)]
    def model(q):
        k = mp.sqrt(1 - q**2); w = m*q/k; Pi = k**5*w*(1 + w**2)*(4 + w**2)
        if m == 0: return -k**5*4/900
        return 1j/450*Pi*(mp.log(k) + mp.digamma(1 + 1j*w)) - k**5*(1 + w**2)*(4 + w**2)/900
    R = [a - model(q) for a, q in zip(A, qs)]
    tq = lambda v: [mp.fsum(x*q**(-p) for x, q in zip(v, qs))/NQ for p in range(NQ)]
    Rt, Bt = tq(R), tq(B)
    Pit = tq([m*q*(1 + (m**2 - 1)*q**2)*(4 + (m**2 - 4)*q**2) for q in qs])
    top = 14 if NQ == 24 else 20
    print(f"m={m:+d}: |remainder Taylor coeff| p=0..5: {[mp.nstr(abs(Rt[p]), 3) for p in range(6)]}  max p=6..{top}: {mp.nstr(max(abs(Rt[p]) for p in range(6, top)), 3)}"
          f"  | max|B6 - (i/450)Pi| = {mp.nstr(max(abs(Bt[p] - 1j*Pit[p]/450) for p in range(top)), 3)}")
