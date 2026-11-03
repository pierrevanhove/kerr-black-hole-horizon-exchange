"""Flux balance at M omega = 0.3 (eps = 0.6), including the superradiant regime, with nu followed by continuation
in eps from the low-frequency root.  Beyond eps ~ 0.7 nu reaches 3/2 and moves into the complex plane:
see complexnu_check.py (real-axis continuation stalls at the spurious point nu = 3/2)."""
import mpmath as mp
from mst import MST, I
from phys import spheroidal_c
from csm import D_TS
import flux_check as fc
mp.mp.dps = 60


def follow(l, m, q, eps_targets, N=80, step=mp.mpf('0.05')):
    s = -2; q = mp.mpf(q)
    d0 = fc.setup(l, m, q, mp.mpf('0.2'), N=N)
    nu, e = d0['nu'], mp.mpf('0.2')
    out = {}
    for et in sorted(mp.mpf(x) for x in eps_targets):
        while e < et - mp.mpf('1e-30'):
            e = min(e + step, et)
            c = q*e/2
            E, _ = spheroidal_c(s, l, m, c)
            lam = E - 2*m*c + c**2 - s*(s + 1)
            k0 = mp.sqrt(1 - q**2)
            P = MST(e, q, m, s, lam, k0, logkappa=mp.log(k0), N=N)
            nu = mp.findroot(P.g, nu, tol=mp.mpf(10)**(-55))
        d = dict(s=s, q=q, eps=e, c=c, lam=lam, k0=k0, nu=nu, P=P, gres=abs(P.g(nu)), N=N, m=m, l=l)
        out[et] = d
    return out


def check(d):
    l, m = d['l'], d['m']
    P, A = fc.branch(d, +1)
    om, k0, lam, c = d['eps']/2, d['k0'], d['lam'], d['c']
    rp, rm = 1 + k0, 1 - k0
    kH = om - m*d['q']/(2*rp); kap = (rp - rm)/(4*rp)
    PH = (2*rp)**5*kH*(kH**2 + kap**2)*(kH**2 + 4*kap**2)*abs(A['Btr']/A['Binc'])**2/om**5
    D = D_TS(lam, c, m)
    R = (D**2 + 144*om**2)/(256*om**8)*abs(A['Bref']/A['Binc'])**2
    a = A['a']; tail = max(abs(a[n]) for n in (d['N'], -d['N']))
    return d['nu'], d['gres'], R, PH, R + PH - 1, tail, kH


if __name__ == '__main__':
  for (m, q) in [(2, '0.95'), (2, 0), (-2, '0.95')]:
    res = follow(2, m, q, ['0.6'])
    for e, d in res.items():
        nu, g, R, PH, bal, tail, kH = check(d)
        print(f"m={m:+d} q={mp.nstr(mp.mpf(q),3):>4} eps={mp.nstr(e,3):>3}: nu={mp.nstr(nu,12):>32} |g|={mp.nstr(g,2):>7}"
              f"  R={mp.nstr(R,12):>15}  P_H={mp.nstr(PH,12):>16}{' (superradiant)' if kH < 0 else '':14}"
              f"  R+P_H-1={mp.nstr(bal,2):>8}  |a_N|={mp.nstr(tail,2)}")
