"""Continue the physical-region eps^7 closed forms (l=2, m=2) to super-extremal q along Im q > 0
(kappa -> -i sqrt(q^2-1), ln kappa -> ln sqrt(q^2-1) - i pi/2) and compare with the direct MST data."""
import json, mpmath as mp
mp.mp.dps = 40
m = 2
fit = json.load(open('e7_fit.json')); poly = [mp.mpmathify(s) for s in fit['poly']]; u0 = mp.mpmathify(fit['u0'])
def parts(q):
    kq = mp.sqrt(q**2 - 1); k = -1j*kq; lk = mp.log(kq) - 1j*mp.pi/2
    t = -m*q/k
    p = t*(1 + t**2)*(4 + t**2); dp = (1 + t**2)*(4 + t**2) + 2*t**2*(4 + t**2) + 2*t**2*(1 + t**2)
    P, dP = k**5*p, k**5*dp
    Re = (mp.digamma(1 + 1j*t) + mp.digamma(1 - 1j*t))/2
    dRe = 1j*(mp.psi(1, 1 + 1j*t) - mp.psi(1, 1 - 1j*t))/2
    T6 = -1j/450*P*(lk + Re); dT6 = -1j/450*(dP*(lk + Re) + P*dRe)/k
    cth = mp.coth(mp.pi*t); dcth = -mp.pi*mp.csch(mp.pi*t)**2
    D6 = -(mp.pi/900)*P*(1 + cth); dD6 = -(mp.pi/900)*(dP*(1 + cth) + P*dcth)/k
    E7 = dT6 + (u0 + mp.mpf(5*m)/9*q)*T6 + mp.fsum(c*q**j for j, c in enumerate(poly))
    D7 = dD6 + (u0 + mp.mpf(5*m)/9*q)*D6
    return E7, D7
data = json.load(open('se7_A7.json'))
print("   q      E7 closed-form (continued)          MST eps^7 (even)            rel.diff  |  D7 closed-form            MST (odd)            rel.diff")
for qs, r in data:
    q = mp.mpf(qs); E7, D7 = parts(q)
    Ae, Ao = mp.mpmathify(r['even'][0]), mp.mpmathify(r['odd'][0])
    print(f"{mp.nstr(q, 5):>6}  {mp.nstr(E7, 12):>34} {mp.nstr(Ae, 12):>34}  {mp.nstr(abs(E7-Ae)/abs(Ae), 2):>8} | "
          f"{mp.nstr(D7, 10):>30} {mp.nstr(Ao, 10):>30}  {mp.nstr(abs(D7-Ao)/abs(Ao), 2):>8}")
