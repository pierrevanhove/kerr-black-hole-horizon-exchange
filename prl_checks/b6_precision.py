"""B6 (eps^6 ln eps coefficient) for every m with the joint log bookkeeping of Supplement S5:
c_n(r) = A_n + B_n ln r + sum_{k != n} B_k r^(k-n) i W_(k-n), solved for n = 6,7,8 from two radii, both branches separately."""
import json, mpmath as mp
mp.mp.dps = 60
NE = 12
th = [(mp.pi*(2*j + 1)/NE if mp.pi*(2*j + 1)/NE <= mp.pi else mp.pi*(2*j + 1)/NE - 2*mp.pi) for j in range(NE)]
W = {d: mp.fsum(t*mp.expj(d*t) for t in th)/NE for d in range(-4, 5)}
def load(m):
    if m == 2:
        out = {}
        for RE in ('2e-5', '1e-5'):
            data = json.load(open(f'phys_l2_{RE}.json')); D = {(r['sign'], r['j']): r['pts'] for r in data}
            out[RE] = {sg: [[mp.mpc(*map(mp.mpf, D[(sg, j)][k]['St'])) for k in range(40)] for j in range(NE)] for sg in (1, -1)}
        return 40, out
    data = json.load(open('phys_allm.json')); D = {(r['m'], r['sign'], r['RE'], r['j']): r['St'] for r in data}
    return 24, {RE: {sg: [[mp.mpc(*map(mp.mpf, D[(m, sg, RE, j)][k])) for k in range(24)] for j in range(NE)] for sg in (1, -1)} for RE in ('2e-5', '1e-5')}
for m in (2, 1, 0, -1, -2):
    NQ, F = load(m)
    qs = [mp.mpf('0.3')*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
    worst_naive, worst_joint = 0, 0
    for sg in (1, -1):
        Bj, Bn = [], []
        for k in range(NQ):
            C = {}
            for RE in ('2e-5', '1e-5'):
                eps = [mp.mpf(RE)*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
                for n in (6, 7, 8):
                    C[(RE, n)] = mp.fsum(F[RE][sg][j][k]*eps[j]**(-n) for j in range(NE))/NE
            rows, rhs = [], []
            for RE in ('2e-5', '1e-5'):
                r = mp.mpf(RE)
                for n in (6, 7, 8):
                    row = []
                    for kk in (6, 7, 8):
                        row += [1 if kk == n else 0, (mp.log(r) if kk == n else r**(kk - n)*1j*W[kk - n])]
                    rows.append(row); rhs.append(C[(RE, n)])
            sol = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
            Bj.append(sol[1]); Bn.append((C[('2e-5', 6)] - C[('1e-5', 6)])/mp.log(2))
        tq = lambda v: [mp.fsum(x*q**(-p) for x, q in zip(v, qs))/NQ for p in range(10)]
        Pi = tq([m*q*(1 + (m**2 - 1)*q**2)*(4 + (m**2 - 4)*q**2) for q in qs])
        Bjt, Bnt = tq(Bj), tq(Bn)
        scale = max(abs(1j*Pi[p]/450) for p in range(10)) or mp.mpf('1e-3')
        worst_joint = max(worst_joint, max(abs(Bjt[p] - 1j*Pi[p]/450) for p in range(10))/scale)
        worst_naive = max(worst_naive, max(abs(Bnt[p] - 1j*Pi[p]/450) for p in range(10))/scale)
    print(f"m={m:+d}: max relative deviation of B6 from i Pi_m/450 (p <= 9, both branches): naive two-radius {mp.nstr(worst_naive, 3)}   joint bookkeeping {mp.nstr(worst_joint, 3)}")
