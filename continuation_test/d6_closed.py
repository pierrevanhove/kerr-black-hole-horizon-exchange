import json, mpmath as mp
mp.mp.dps = 50
R = {tuple(int(t) for t in k.split('_')): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in json.load(open('l2_rich.json')).items()}
num = [mp.re((R[(1, 6, p)] - R[(-1, 6, p)])/2) for p in range(11)]
eta = 1
D6 = lambda q: (2*mp.pi/225)*q*(1 + 3*q**2)*(eta - mp.coth(2*mp.pi*q/mp.sqrt(1 - q**2)))
tay = mp.taylor(D6, mp.mpf('1e-40'), 10) if False else None
# Taylor coefficients of D6 at q=0 via Cauchy on a small circle (D6 is analytic at q=0: q*coth(...) is even and regular)
N, r = 64, mp.mpf('0.2')
qs = [r*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
vals = [D6(q) for q in qs]
coef = [mp.fsum(v*q**(-p) for v, q in zip(vals, qs))/N for p in range(11)]
z2, z4, z6 = mp.zeta(2), mp.zeta(4), mp.zeta(6)
rho2 = mp.mpf(4)/45 + 64*mp.j*(mp.j/120*(120*z2 - 1152*z4 + 1024*z6 - 1))/15
print(" p   closed form                  full-MST duality-odd part        difference")
for p in range(11):
    print(f"{p:2d}  {mp.nstr(mp.re(coef[p]), 18):>26}  {mp.nstr(num[p], 18):>26}  {mp.nstr(mp.re(coef[p]) - num[p], 3):>10}")
print("\nq^6 coefficient vs BBCPSS rho_2/64:", mp.nstr(mp.re(coef[6]), 18), mp.nstr(mp.re(rho2)/64, 18))
