"""Re-derive every number quoted in the updated draft (closed forms, identities, unified physical-branch form)."""
import json, mpmath as mp
mp.mp.dps = 40
Z = mp.zeta
# (1) far-zone O(eps) phase coefficient at q=0, and u0
se = json.load(open('se_taylor_l2.json'))
s10 = mp.mpc(*map(mp.mpf, se['1'][0]))
f1 = 2*mp.euler + 2*mp.log(2) - mp.mpf(37)/12
print("(1) s1(c=0) =", mp.nstr(s10, 20), "  i(2gE+2ln2-37/12) =", mp.nstr(1j*f1, 20))
u0 = mp.pi + 1j*(f1 + mp.mpf(107)/210)
print("    u0 = pi + i(2gE + 2ln2 - 37/12 - nu2) =", mp.nstr(u0, 15))
# (2) Pi_m and Gamma form of the duality-odd function
for m in (2, 1, 0, -1, -2):
    q = mp.mpf('0.37'); k = mp.sqrt(1 - q**2); w = m*q/k
    Pi = k**5*w*(1 + w**2)*(4 + w**2); Pi2 = m*q*(1 + (m**2 - 1)*q**2)*(4 + (m**2 - 4)*q**2)
    for eta in (1, -1):
        if m == 0:
            Dc = -k**5/225
        else:
            Dc = -(mp.pi/900)*Pi*(1/mp.tanh(mp.pi*w) - eta)
        Dg = -(k**5/900)*mp.gamma(3 + 1j*w)*mp.gamma(3 - 1j*w)*mp.exp(-eta*mp.pi*w)
        t0 = -m*q/k
        Dold = -k**5/225 if m == 0 else -(mp.pi/900)*k**5*t0*(1 + t0**2)*(4 + t0**2)*(eta + 1/mp.tanh(mp.pi*t0))
        assert abs(Dc - Dg) < 1e-30 and abs(Dc - Dold) < 1e-30
    assert abs(Pi - Pi2) < 1e-30
print("(2) Pi_m = m q [1+(m^2-1)q^2][4+(m^2-4)q^2];  D6 = -(pi/900)Pi(coth(pi w)-eta) = -(k^5/900)|G(3+iw)|^2 e^{-eta pi w}: OK")
# (3) Boltzmann identity pi(tau+eps) = (omega - m Omega_H)/(2 T_H), M = 1
q, eps, m = mp.mpf('0.6'), mp.mpf('0.013'), 2
k = mp.sqrt(1 - q**2); rp = 1 + k; TH = k/(4*mp.pi*rp); OH = q/(2*rp)
print("(3) pi(tau+eps) =", mp.nstr(mp.pi*((eps - m*q)/k + eps), 15), "  (w - m OmegaH)/(2TH) =", mp.nstr((eps/2 - m*OH)/(2*TH), 15))
# (4) one function, two limits (m = 2): C(q) = 160 F(q)
F = lambda q: mp.mpf(4)/225*q*(1 + 3*q**2)*(mp.log(mp.sqrt(1 - q**2)) + (mp.digamma(1 + 2j*q/mp.sqrt(1 - q**2)) + mp.digamma(1 - 2j*q/mp.sqrt(1 - q**2)))/2)
N, r = 96, mp.mpf('0.3'); qs = [r*mp.expjpi(2*mp.mpf(j)/N) for j in range(N)]
c5 = mp.re(mp.fsum(F(z)*z**(-5) for z in qs)/N)
b2 = mp.mpf(128)/45*(Z(3) - Z(5)) - mp.mpf(14)/45
print("(4) 160[q^5]F =", mp.nstr(160*c5, 20), "  16 b2 =", mp.nstr(16*b2, 20))
Fc = lambda q: mp.mpf(4)/225*q*(1 + 3*q**2)*(mp.log(mp.sqrt(q**2 - 1)) + (mp.digamma(1 + 2*q/mp.sqrt(q**2 - 1)) + mp.digamma(1 - 2*q/mp.sqrt(q**2 - 1)))/2)
print("    160 F/q^5 at q = 1e8 (continued):", mp.nstr(160*Fc(mp.mpf(10)**8)/mp.mpf(10)**40, 15), "  64/15 =", mp.nstr(mp.mpf(64)/15, 15))
# (5) duality-odd: q^6 Taylor coefficient = rho_|m|/64 (BCPSS 4.27) and large-q limit
c1 = 1j/15*(10*Z(2) - 92*Z(4) + 80*Z(6) - 1); c2 = 1j/120*(120*Z(2) - 1152*Z(4) + 1024*Z(6) - 1)
rho = {2: mp.mpf(4)/45 + 64j*c2/15, 1: -8j*c1/5 + 16j*c2/15 + mp.mpf(4)/45, 0: mp.mpf(4)/45}
def D6(m, q, k, eta=1):
    if m == 0: return -k**5/225
    w = m*q/k; return -(mp.pi/900)*k**5*w*(1 + w**2)*(4 + w**2)*(1/mp.tanh(mp.pi*w) - eta)
for m in (2, 1, 0, -1, -2):
    c6 = mp.fsum(D6(m, z, mp.sqrt(1 - z**2))*z**(-6) for z in qs)/N
    qq = mp.mpf(10)**6; lim = 32*D6(m, qq, -1j*mp.sqrt(qq**2 - 1))/qq**5
    tgt = (-1)**m*mp.mpf(8)/225*1j*mp.factorial(2 + m)*mp.factorial(2 - m)
    print(f"(5) m={m:+d}: [q^6]D6 = {mp.nstr(mp.re(c6), 16):>22}  rho/64 = {mp.nstr(mp.re(rho[abs(m)])/64, 16):>22} | 32 D6/q^5 -> {mp.nstr(lim, 10)}  target {mp.nstr(tgt, 10)}")
# (6) pole terms of the duality-even function
for m in (2, 1, -1, -2):
    qq = mp.mpf(10)**10; k = -1j*mp.sqrt(qq**2 - 1); w = m*qq/k
    Pi = k**5*w*(1 + w**2)*(4 + w**2)
    pole = 32*(1j/450)*Pi*(mp.digamma(1 + 1j*w) + mp.digamma(1 - 1j*w))/2/qq**5
    print(f"(6) m={m:+d}: pole_m = {mp.nstr(pole, 12)}   x75/i = {mp.nstr(pole*75/1j, 10)}")
