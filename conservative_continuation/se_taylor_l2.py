# Taylor coefficients of the super-extremal classical s1(l=2; c) for both branches (from the 28-node, 110-digit run)
import json, mpmath as mp
mp.mp.dps = 110
data = json.load(open('../conservative_test/polar_run2.json'))
D = {(d['l'], d['k']): d for d in data}
NC = 28
cs = [mp.mpf(D[(2, k)]['c']) for k in range(NC)]
V = mp.matrix([[c**j for j in range(NC)] for c in cs])
out = {}
for eta in (1, -1):
    vals = [mp.mpc(mp.mpf(D[(2, k)][f's1_{eta}'][0]), mp.mpf(D[(2, k)][f's1_{eta}'][1])) for k in range(NC)]
    p = mp.lu_solve(V, mp.matrix(vals))
    out[eta] = [p[j] for j in range(12)]
mp.mp.dps = 30
for n in range(9):
    e, o = (out[1][n] + out[-1][n])/2, (out[1][n] - out[-1][n])/2
    print(f"c^{n}: eta-even {mp.nstr(e, 20):>48}   eta-odd {mp.nstr(o, 16)}")
json.dump({str(eta): [[str(mp.re(v)), str(mp.im(v))] for v in out[eta]] for eta in out}, open('se_taylor_l2.json', 'w'))
