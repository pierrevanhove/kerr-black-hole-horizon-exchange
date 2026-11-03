import mpmath as mp
mp.mp.dps = 40
z = {n: mp.zeta(n) for n in range(2, 8)}
c1 = mp.j/15*(10*z[2] - 92*z[4] + 80*z[6] - 1); c2 = mp.j/120*(120*z[2] - 1152*z[4] + 1024*z[6] - 1)
rho = {2: mp.mpf(4)/45 + 64*mp.j*c2/15, 1: -8*mp.j*c1/5 + 16*mp.j*c2/15 + mp.mpf(4)/45, 0: mp.mpf(4)/45}
def D6(m, q, kap, eta=1):
    t0 = -m*q/kap
    if m == 0:
        return -kap**5/225                              # tau0 -> 0 limit: tau0*coth(pi tau0) -> 1/pi
    return -(mp.pi/900)*kap**5*t0*(1 + t0**2)*(4 + t0**2)*(eta + mp.coth(mp.pi*t0))
# (a) physical-region q^6 coefficients vs BBCPSS rho_m/64 (eqs. 4.26-4.27)
N, r = 64, mp.mpf('0.25')
print("q^6 Taylor coefficient of the coth-part vs BBCPSS rho_m/64:")
for m in (2, 1, 0, -1, -2):
    qs = [r*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
    c6 = mp.fsum(D6(m, q, mp.sqrt(1 - q**2))*q**(-6) for q in qs)/N
    print(f"  m={m:+d}: {mp.nstr(mp.re(c6), 20):>26}   rho_|m|/64 = {mp.nstr(mp.re(rho[abs(m)])/64, 20)}")
# (b) continuation to super-extremal spin along Im q > 0: kappa -> -i sqrt(q^2-1); classical eta-odd a^5 coefficient
print("\nlarge-spin limit 32*D6/q^5 (upper path) vs closed-form eta-odd a^5 coefficient (-1)^m (8/225) i (2+m)!(2-m)!:")
for m in (2, 1, 0, -1, -2):
    vals = []
    for qs in ('1e2', '1e3', '1e4'):
        q = mp.mpf(qs)
        vals.append(32*D6(m, q, -mp.j*mp.sqrt(q**2 - 1))/q**5)
    target = (-1)**m*mp.mpf(8)/225*mp.j*mp.factorial(2 + m)*mp.factorial(2 - m)
    print(f"  m={m:+d}: q=1e2,1e3,1e4 -> {[mp.nstr(v, 10) for v in vals]}   target {mp.nstr(target, 10)}")
