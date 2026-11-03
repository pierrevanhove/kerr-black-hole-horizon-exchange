import mpmath as mp
from multiprocessing import Pool
from d6_m import NE
def job(qs):
    mp.mp.dps = 90
    from mst import MST, I
    from phys import spheroidal_c
    from csm import D_TS, secant
    m, q = 0, mp.mpf(qs)
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
        r = mp.mpf(rs); vals = []
        for j in range(NE):
            e = r*mp.expjpi(2*(j + mp.mpf(1)/2)/NE)
            vals.append((e, (St(e, k0, mp.log(k0)) - St(e, -k0, mp.log(k0) + I*mp.pi))/2))
        out[rs] = {k: mp.fsum(d*e**(-k) for e, d in vals)/NE for k in range(0, 8)}
    rich = {k: (2**(8-k)*out['1e-5'][k] - out['2e-5'][k])/(2**(8-k) - 1) for k in range(8)}
    return qs, {k: str(v) for k, v in rich.items()}
if __name__ == '__main__':
    with Pool(2) as p:
        for qs, d in p.map(job, ['0.3', '0.7']):
            mp.mp.dps = 25
            q = mp.mpf(qs); k0 = mp.sqrt(1 - q**2)
            print(f"m=0 q={qs}: " + "  ".join(f"D{k}={mp.nstr(mp.mpmathify(v), 12)}" for k, v in d.items() if k >= 4)
                  + f"   | -kappa^5/225 = {mp.nstr(-k0**5/225, 12)}")
