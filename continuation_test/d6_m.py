"""Physical-region check of the general-m conjecture
   D6^(m)(q) = -(pi/900) kappa^5 tau0 (1+tau0^2)(4+tau0^2) [eta + coth(pi tau0)],  tau0 = -m q/kappa, kappa = sqrt(1-q^2),
   for l=2, s=-2 (eta=+1 with log kappa_dual = log kappa + i pi)."""
import mpmath as mp
from multiprocessing import Pool
NE = 16
def job(args):
    m, qs = args
    mp.mp.dps = 90
    from mst import MST, I
    from phys import spheroidal_c
    from csm import D_TS, secant
    q = mp.mpf(qs)
    k0 = mp.sqrt(1 - q**2)
    def St(e, kappa, logk, l=2, s=-2):
        c = q*e/2
        E, b = spheroidal_c(s, l, m, c)
        lam = E - 2*m*c + c**2 - s*(s+1)
        P = MST(e, q, m, s, lam, kappa, logkappa=logk, N=20)
        x = secant(lambda x: P.g(l + x*e**2)/e**2, -0.5095, mp.mpf(10)**(-45))
        P.nu = l + x*e**2
        A = P.amplitudes()
        return (-1)**(l+1)*D_TS(lam, c, m)/(16*(e/2)**4)*A['S']*mp.exp(-2*I*e*mp.log(e))
    out = {}
    for rs in ('2e-5', '1e-5'):
        r = mp.mpf(rs)
        vals = []
        for j in range(NE):
            e = r*mp.expjpi(2*(j + mp.mpf(1)/2)/NE)
            vals.append((e, (St(e, k0, mp.log(k0)) - St(e, -k0, mp.log(k0) + I*mp.pi))/2))
        out[rs] = {k: mp.fsum(d*e**(-k) for e, d in vals)/NE for k in range(0, 8)}
    d6 = (4*out['1e-5'][6] - out['2e-5'][6])/3
    low = max(abs(out['1e-5'][k]) for k in range(6))
    t0 = -m*q/k0
    pred = -(mp.pi/900)*k0**5*t0*(1 + t0**2)*(4 + t0**2)*(1 + mp.coth(mp.pi*t0)) if m != 0 else mp.nan
    return m, qs, str(d6), str(low), str(pred)
if __name__ == '__main__':
    jobs = [(m, q) for m in (2, 1, -1, -2) for q in ('0.3', '0.7')]
    with Pool(2) as p:
        for m, qs, d6, low, pred in p.map(job, jobs):
            mp.mp.dps = 30
            d6, pred = mp.mpmathify(d6), mp.mpmathify(pred)
            print(f"m={m:+d} q={qs}: D6 MST = {mp.nstr(d6, 15):>36}   conjecture = {mp.nstr(pred, 15):>22}   rel.diff {mp.nstr(abs(d6-pred)/abs(pred), 2)}   (max |D_k<6| {mp.nstr(mp.mpmathify(low), 2)})")
