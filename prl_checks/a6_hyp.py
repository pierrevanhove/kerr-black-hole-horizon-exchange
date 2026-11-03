"""Candidate a^6 pole mechanism: keep the polynomial prefactor at tau0 = -m q/kappa but shift the argument of the
singular function to the exact tau = (eps - m q)/kappa;  eps-derivative -> trigamma / csch^2 double pole."""
import mpmath as mp
mp.mp.dps = 40
def P(m, q, kap):
    t0 = -m*q/kap; return kap**5*t0*(1 + t0**2)*(4 + t0**2)
for m in (2, 1, -1, -2):
    out = []
    for qs in ('1e3', '1e4', '1e5'):
        q = mp.mpf(qs); kap = -1j*mp.sqrt(q**2 - 1); t0 = -m*q/kap
        # self-dual: -(i/450) P(tau0) d/deps Re psi(1+i tau)
        Xe = -1j/450*P(m, q, kap)*(1j/kap)*(mp.psi(1, 1 + 1j*t0) - mp.psi(1, 1 - 1j*t0))/2
        # anti-self-dual: -(pi/900) kappa^5 p(tau0) d/deps coth(pi tau)
        Xo = -(mp.pi/900)*P(m, q, kap)*(mp.pi/kap)*(-mp.csch(mp.pi*t0)**2)
        out.append((64*Xe/q**6, 64*Xo/q**6))
    print(f"m={m:+d}:  64 X7_even/q^6 -> {[mp.nstr(v[0], 10) for v in out]}")
    print(f"        64 X7_odd/q^6  -> {[mp.nstr(v[1], 10) for v in out]}")
f = mp.factorial
print("targets: pole part at a^6 = (10m/9) x a^5 pole part;  m=2: 256/135 =", mp.nstr(mp.mpf(256)/135, 10),
      "; m=1: (10/9)(16/75) =", mp.nstr(mp.mpf(10)/9*mp.mpf(16)/75, 10))
