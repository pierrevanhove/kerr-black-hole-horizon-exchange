"""eps^7 coefficient of the l=2, m=2 S-matrix at super-extremal q (both continued branches): remove the eps^6 ln eps leakage
exactly, separate A7 + B7 ln r (+ C r), and compare 64 A7/q^6 with the classical a^6 coefficients."""
import json, mpmath as mp
from collections import defaultdict
mp.mp.dps = 60
NE = 16
th = []
for j in range(NE):
    t = mp.pi*(2*j + 1)/NE
    th.append(t if t <= mp.pi else t - 2*mp.pi)          # principal arguments of the nodes
Wm1 = mp.fsum(t*mp.expj(-t) for t in th)/NE            # (1/NE) sum theta_j e^{-i theta_j}
data = defaultdict(list)
for fn in ('e6_se.json', 'e6_se2.json', 'e6_se3.json'):
    for qs, rs, ev, od in json.load(open(fn)):
        data[qs].append((mp.mpf(rs), {int(k): mp.mpc(v.replace(' ', '').replace('j', 'j')) if False else mp.mpmathify(v.replace('j', 'j')) for k, v in ev.items()},
                         {int(k): mp.mpmathify(v) for k, v in od.items()}))
tgt_even, tgt_odd = mp.mpf('1.8964816078651064584'), -mp.mpf(256)/135
print(f"targets (classical c^6, Im): eta-even {mp.nstr(tgt_even, 12)}   eta-odd {mp.nstr(tgt_odd, 12)}")
print("   q     nr   64*A7_even/q^6 (Im)        64*A7_odd/q^6 (Im)       B7_even/q^6            B7_odd/q^6")
rows = []
for qs in sorted(data, key=lambda s: float(s)):
    q = mp.mpf(qs)
    pts = sorted(data[qs], key=lambda t: t[0])
    B6 = 4j*q*(1 + 3*q**2)/225
    res = {}
    for nm, idx in (('even', 1), ('odd', 2)):
        ys = []
        for r, ev, od in pts:
            c7 = (ev if idx == 1 else od)[7]
            if idx == 1:
                c7 = c7 - B6/r*1j*Wm1
            ys.append((r, c7))
        if len(ys) >= 3:
            M = mp.matrix([[1, mp.log(r), r] for r, _ in ys[:3]]); sol = mp.lu_solve(M, mp.matrix([y for _, y in ys[:3]]))
            A, B = sol[0], sol[1]
        else:
            (r1, y1), (r2, y2) = ys[0], ys[1]
            B = (y2 - y1)/mp.log(r2/r1); A = y1 - B*mp.log(r1)
        res[nm] = (A, B)
    rows.append((q, res))
    print(f"{qs:>6} {len(pts):3d}   {mp.nstr(mp.im(64*res['even'][0]/q**6), 14):>22}   {mp.nstr(mp.im(64*res['odd'][0]/q**6), 14):>22}   "
          f"{mp.nstr(res['even'][1]/q**6, 4):>22}   {mp.nstr(res['odd'][1]/q**6, 4):>22}")
json.dump([[str(q), {k: [str(v[0]), str(v[1])] for k, v in r.items()}] for q, r in rows], open('se7_A7.json', 'w'))
# large-q extrapolation of 64 A7/q^6 (polynomial in 1/q through the largest-q points), Supplement S9(h)
big = [(q, r) for q, r in rows if q >= 300]
for npt in (3, 4, len(big)):
    sel = big[-npt:]
    M = mp.matrix([[q**(-k) for k in range(npt)] for q, _ in sel])
    lim = {nm: mp.lu_solve(M, mp.matrix([mp.im(64*r[nm][0]/q**6) for q, r in sel]))[0] for nm in ('even', 'odd')}
    print(f"extrapolation q -> infinity ({npt} points, q >= {mp.nstr(sel[0][0], 4)}):  eta-even {mp.nstr(lim['even'], 10)}"
          f"   eta-odd {mp.nstr(lim['odd'], 10)}   (targets {mp.nstr(tgt_even, 10)}, {mp.nstr(-tgt_odd, 10)} = 256/135)")
