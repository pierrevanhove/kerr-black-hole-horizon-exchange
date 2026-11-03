"""Anti-self-dual eps^7 closed-form test for all m (Richardson-corrected in the eps-circle radius)."""
import json, mpmath as mp
mp.mp.dps = 60
def coefs(load_fn, NQ, RQ):
    qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
    G = load_fn()
    out = {}
    for key, vals in G.items():
        out[key] = [mp.fsum(v*q**(-p) for v, q in zip(vals, qs))/NQ for p in range(NQ)]
    return out
NE = 12
def eps_nodes(RE): return [mp.mpf(RE)*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
def load_m2():
    out = {}
    for RE in ('2e-5', '1e-5'):
        data = json.load(open(f'phys_l2_{RE}.json')); D = {(r['sign'], r['j']): r['pts'] for r in data}
        eps = eps_nodes(RE)
        for n in (6, 7):
            vals = []
            for k in range(40):
                F = {sg: [mp.mpc(mp.mpf(D[(sg, j)][k]['St'][0]), mp.mpf(D[(sg, j)][k]['St'][1])) for j in range(NE)] for sg in (1, -1)}
                vals.append(mp.fsum((F[1][j] - F[-1][j])/2*eps[j]**(-n) for j in range(NE))/NE)
            out[(RE, n)] = vals
    return out
def load_m(m):
    data = json.load(open('phys_allm.json')); D = {(r['m'], r['sign'], r['RE'], r['j']): r['St'] for r in data}
    out = {}
    for RE in ('2e-5', '1e-5'):
        eps = eps_nodes(RE)
        for n in (6, 7):
            vals = []
            for k in range(24):
                F = {sg: [mp.mpc(mp.mpf(D[(m, sg, RE, j)][k][0]), mp.mpf(D[(m, sg, RE, j)][k][1])) for j in range(NE)] for sg in (1, -1)}
                vals.append(mp.fsum((F[1][j] - F[-1][j])/2*eps[j]**(-n) for j in range(NE))/NE)
            out[(RE, n)] = vals
    return out
kap = lambda q: mp.sqrt(1 - q**2)
def D6f(m, q):
    k = kap(q)
    if m == 0: return -k**5/225
    t = -m*q/k; return -(mp.pi/900)*k**5*t*(1 + t**2)*(4 + t**2)*(1 + mp.coth(mp.pi*t))
def dlog(m, q):
    k = kap(q)
    if m == 0: return mp.pi/k
    t = -m*q/k; return (1/t + 2*t/(1 + t**2) + 2*t/(4 + t**2) + mp.pi*(1 - mp.coth(mp.pi*t)))/k
def tay(f, P, N=160, r=mp.mpf('0.3')):
    zs = [r*mp.expjpi(2*(mp.mpf(k) + mp.mpf(1)/3)/N) for k in range(N)]
    v = [f(z) for z in zs]; return [mp.fsum(vv*z**(-p) for vv, z in zip(v, zs))/N for p in range(P)]
results = {}
for m, loader, NQ in ((2, load_m2, 40), (1, lambda: load_m(1), 24), (-1, lambda: load_m(-1), 24), (-2, lambda: load_m(-2), 24), (0, lambda: load_m(0), 24)):
    C = coefs(loader, NQ, mp.mpf('0.3'))
    D7 = [2*C[('1e-5', 7)][p] - C[('2e-5', 7)][p] for p in range(NQ)]       # Richardson: remove O(r) leakage
    P = 22 if NQ == 40 else 14
    T6, Tsh = tay(lambda q: D6f(m, q), P), tay(lambda q: D6f(m, q)*dlog(m, q), P)
    cols = [T6, [T6[p-1] if p >= 1 else 0 for p in range(P)], [T6[p-2] if p >= 2 else 0 for p in range(P)]]
    M = mp.matrix([[cols[c][p] for c in range(3)] for p in range(P)])
    sol, res = mp.qr_solve(M, mp.matrix([D7[p] - Tsh[p] for p in range(P)]))
    results[m] = sol
    print(f"m={m:+d}: D7 - D6*dlog/dtau0/kappa = (u0 + u1 q + u2 q^2) D6 :  u0 = {mp.nstr(sol[0], 13)}  u1 = {mp.nstr(sol[1], 13)}"
          f"  u2 = {mp.nstr(sol[2], 3)}   resid {mp.nstr(res, 3)}   [(5/9)m = {mp.nstr(mp.mpf(5*m)/9, 10)}]")
u0 = results[2][0]
print("\nu0 - pi =", mp.nstr(u0 - mp.pi, 15))
print("identify Im u0:", mp.identify(mp.im(u0), ['log(2)', 'euler', 'pi']), mp.nstr(mp.im(u0), 15))
