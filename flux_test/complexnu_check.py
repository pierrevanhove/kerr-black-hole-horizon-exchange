"""Flux balance on the complex-nu branch (nu = 3/2 + i y) at M omega = 0.4 - 0.6."""
import mpmath as mp
from mst import MST, I
from phys import spheroidal_c
from csm import D_TS
import flux_check as fc
from highfreq_check import check
mp.mp.dps = 60
for (m, q) in [(2, '0.95'), (2, 0), (-2, '0.95')]:
    q = mp.mpf(q)
    for e in ('0.8', '1.0', '1.2'):
        e = mp.mpf(e); s = -2; l = 2; N = 80
        c = q*e/2
        E, _ = spheroidal_c(s, l, m, c)
        lam = E - 2*m*c + c**2 - s*(s + 1)
        k0 = mp.sqrt(1 - q**2)
        P = MST(e, q, m, s, lam, k0, logkappa=mp.log(k0), N=N)
        best = None
        for guess in (mp.mpc('1.5', '0.05'), mp.mpc('1.5', '0.2'), mp.mpc('1.5', '0.5'), mp.mpc('-0.5', '0.3')):
            try:
                nu = mp.findroot(P.g, guess, tol=mp.mpf(10)**(-55))
            except Exception:
                continue
            if abs(mp.im(nu)) > mp.mpf('1e-6') and abs(P.g(nu)) < mp.mpf('1e-45'):
                best = nu; break
        if best is None:
            print(m, q, e, 'no complex root found'); continue
        d = dict(s=s, q=q, eps=e, c=c, lam=lam, k0=k0, nu=best, P=P, gres=abs(P.g(best)), N=N, m=m, l=l)
        nu, g, R, PH, bal, tail, kH = check(d)
        print(f"m={m:+d} q={mp.nstr(q,3):>4} eps={mp.nstr(e,3):>3}: nu={mp.nstr(nu,12):>34} |g|={mp.nstr(g,2):>7}"
              f"  R={mp.nstr(R,12):>15}  P_H={mp.nstr(PH,12):>16}{' (superradiant)' if kH < 0 else '':14}"
              f"  R+P_H-1={mp.nstr(bal,2):>8}")
