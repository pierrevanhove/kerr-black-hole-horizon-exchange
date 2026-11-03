import json, mpmath as mp
mp.mp.dps = 110
R, NC = mp.mpf('0.2'), 28
data = json.load(open('polar_run2.json'))
ls = sorted(set(d['l'] for d in data)); ks = list(range(NC))
D = {(d['l'], d['k']): d for d in data}
cs = [mp.mpf(D[(2, k)]['c']) for k in ks]
Y0 = lambda L: mp.sqrt(mp.mpf(2*L+1)/(4*mp.pi))
Ls = list(range(2, 10))
# --- sanity: root residuals and branch-independence of nu
worst = max(mp.mpf(D[(l, k)][f'g_{e}'][0])/mp.mpf(D[(l, k)][f'g_{e}'][1]) for l in ls for k in ks for e in (1, -1))
dx = max(abs(mp.mpf(D[(l, k)]['x_1'][0]) - mp.mpf(D[(l, k)]['x_-1'][0])) for l in ls for k in ks)
print(f"worst relative |G| at root: {mp.nstr(worst, 3)};  max |x(eta=+1)-x(eta=-1)| = {mp.nstr(dx, 3)}")
def s1(l, k, eta):
    re, im = D[(l, k)][f's1_{eta}']; return mp.mpc(mp.mpf(re), mp.mpf(im))
print("max |Re s1| =", mp.nstr(max(abs(mp.re(s1(l, k, e))) for l in ls for k in ks for e in (1, -1)), 3))
F, last = {}, {}
for eta in (1, -1):
    for L in Ls:
        vals = []
        for k in ks:
            tot = mp.mpc(0)
            for l in ls:
                b = D[(l, k)]['b']
                S0 = mp.fsum(mp.mpf(v)*Y0(int(Lp)) for Lp, v in b.items())
                term = s1(l, k, eta)*S0*mp.mpf(b[str(L)])
                tot += term
                if l == ls[-1]: last[L] = max(last.get(L, 0), abs(term))
            vals.append(tot)
        F[(eta, L)] = vals
print("max |l=16 term| per L:", {L: mp.nstr(last[L], 2) for L in Ls})
V = mp.matrix([[c**j for j in range(NC)] for c in cs])
T = {}
for key, vals in F.items():
    p = mp.lu_solve(V, mp.matrix(vals)); T[key] = [p[j] for j in range(NC)]
json.dump({f"{e}_{L}": [[str(mp.re(v)), str(mp.im(v))] for v in T[(e, L)]] for (e, L) in T}, open('taylor2.json', 'w'))
for n in range(0, 9):
    print(f"n={n}: " + "  ".join(f"L{L}:{mp.nstr(mp.im(T[(1,L)][n]), 12)}" for L in Ls[:5]))
print("max |Re T_n| (n<=8):", mp.nstr(max(abs(mp.re(T[(e, L)][n])) for e in (1, -1) for L in Ls for n in range(9)), 3))
print("tail |T_n| for L=2 (n=10..19):", [mp.nstr(abs(T[(1, 2)][n]), 3) for n in range(10, NC)])
