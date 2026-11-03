"""Super-extremal classical l=2 Taylor coefficients for m in {1,0,-1,-2} (plus m=2 from se_taylor_l2.json) and the
per-m test  eta-even c^5 = 32 i a5^(m) + pole_m  with a5^(m) from the physical-region fit (allm_phys_result.json)."""
import json, mpmath as mp
mp.mp.dps = 110
data = json.load(open('se_allm.json'))
D = {(d['m'], d['k']): d for d in data}
NC = 28
T = {}
for m in (1, 0, -1, -2):
    cs = [mp.mpf(D[(m, k)]['c']) for k in range(NC)]
    V = mp.matrix([[c**j for j in range(NC)] for c in cs])
    for eta in (1, -1):
        vals = [mp.mpc(mp.mpf(D[(m, k)][f's1_{eta}'][0]), mp.mpf(D[(m, k)][f's1_{eta}'][1])) for k in range(NC)]
        p = mp.lu_solve(V, mp.matrix(vals))
        T[(m, eta)] = [p[j] for j in range(12)]
old = json.load(open('se_taylor_l2.json'))
for eta in (1, -1):
    T[(2, eta)] = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in old[str(eta)]]
json.dump({f"{m}_{eta}": [[str(mp.re(v)), str(mp.im(v))] for v in T[(m, eta)]] for (m, eta) in T}, open('se_taylor_allm.json', 'w'))
mp.mp.dps = 30
ph = json.load(open('allm_phys_result.json'))
a5 = {m: mp.mpf(ph[str(m)]['a5']) for m in (1, -1, -2, 0)}
a5[2] = mp.mpf('-0.000023458036042842932649012808215965098753284051787938')     # e6_fit.json (m=2 physical fit)
def pole(m, q=mp.mpf(10)**12):
    """32 i x lim_{q->inf} (-1/450) P^(m)(q) Re psi(1+i tau0) / q^5 on the continued branch kappa = -i sqrt(q^2-1)."""
    kap = -1j*mp.sqrt(q**2 - 1); t0 = -m*q/kap
    P = kap**5*t0*(1 + t0**2)*(4 + t0**2)
    return 32j*(-mp.mpf(1)/450)*P*(mp.digamma(1 + 1j*t0) + mp.digamma(1 - 1j*t0))/2/q**5
print(" m   SE eta-even c^5 (MST, classical)        32 i a5 (phys fit) + pole_m                 diff      | pole_m      | eta-odd c^5   (-1)^m (8/225) i (2+m)!(2-m)!")
for m in (2, 1, 0, -1, -2):
    ev = (T[(m, 1)][5] + T[(m, -1)][5])/2; od = (T[(m, 1)][5] - T[(m, -1)][5])/2
    pm = pole(m) if m != 0 else mp.mpf(0)
    pred = 32j*a5[m] + pm
    odd_t = (-1)**m*mp.mpf(8)/225*mp.factorial(2 + m)*mp.factorial(2 - m)
    print(f"{m:+d}  {mp.nstr(mp.im(ev), 22):>26}i   {mp.nstr(mp.im(pred), 22):>26}i   {mp.nstr(mp.im(ev - pred), 3):>9}   {mp.nstr(mp.im(pm), 12):>14}i  {mp.nstr(mp.im(od), 14):>18}i  {mp.nstr(odd_t, 14)}")
print("\nper-branch c^5 (Im):")
for m in (2, 1, 0, -1, -2):
    print(f"  m={m:+d}: eta=+1 {mp.nstr(mp.im(T[(m,1)][5]), 18):>24}   eta=-1 {mp.nstr(mp.im(T[(m,-1)][5]), 18):>24}   32 a5 = {mp.nstr(32*a5[m], 18)}")
