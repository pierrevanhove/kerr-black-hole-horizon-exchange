import mpmath as mp
mp.mp.dps = 40
F = lambda q: mp.mpf(4)/225*q*(1 + 3*q**2)*(mp.log(mp.sqrt(1 - q**2)) + (mp.digamma(1 + 2j*q/mp.sqrt(1 - q**2)) + mp.digamma(1 - 2j*q/mp.sqrt(1 - q**2)))/2)
N, r = 96, mp.mpf('0.3')
qs = [r*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
v = [F(q) for q in qs]
c5 = mp.re(mp.fsum(vv*q**(-5) for vv, q in zip(v, qs))/N)
b2 = mp.mpf(128)/45*mp.zeta(3) - mp.mpf(128)/45*mp.zeta(5) - mp.mpf(14)/45
print("160 * [q^5] F (Taylor, physical)   =", mp.nstr(160*c5, 25))
print("16 b2 (2610.40101, eq. 4.2)        =", mp.nstr(16*b2, 25))
Fc = lambda q: mp.mpf(4)/225*q*(1 + 3*q**2)*(mp.log(mp.sqrt(q**2 - 1)) + (mp.digamma(1 + 2*q/mp.sqrt(q**2 - 1)) + mp.digamma(1 - 2*q/mp.sqrt(q**2 - 1)))/2)
for q in ('1e3', '1e5', '1e7'):
    q = mp.mpf(q); print(f"160 * F(q)/q^5 continued, q={mp.nstr(q,2)}: {mp.nstr(160*Fc(q)/q**5, 15)}")
print("64/15 =", mp.nstr(mp.mpf(64)/15, 15))
