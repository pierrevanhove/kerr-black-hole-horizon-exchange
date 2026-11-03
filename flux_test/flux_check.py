"""Unit-flux normalization of the MST far-zone factor F_lm (s = -2), cf. Supplement S10.
M = 1, eps = 2 omega, q = a, kappa = sqrt(1-q^2), r_pm = 1 pm kappa.
For each (l, m, q, eps), with the same MST code (continued fraction for nu, coefficients a_n) as all other checks:
  1. asymptotic coefficients of the two Floquet solutions, computed separately:
       B^inc[R_C^nu] = A_+^nu e^{-i phi}/omega,  B^ref[R_C^nu] = A_-^nu e^{+i phi}/omega^(1+2s),
       and the same with nu -> -nu-1, a_n -> a_{-n}  (phi = eps ln eps - (1-kappa) eps/2)
  2. F_P = (-1)^(l+1) C_P/(16 omega^4) B^ref[R_C^nu]/B^inc[R_C^nu],  C_P = D + 12 i P M omega,
     |C_P|^2 = D^2 + 144 M^2 omega^2  (= Eq. for |C|^2 in S10, checked below);  P = +-1 the parity.
  3. S_P(factorized) = F_P (1 + i e^{i pi nu} K)/(1 - i e^{-i pi nu} rho K),  K = K_{-nu-1}/K_nu,
     S_P(direct)     = (-1)^(l+1) C_P/(16 omega^4) B^ref/B^inc of R_in = K_nu R_C^nu + K_{-nu-1} R_C^{-nu-1}.
  4. Flux balance |S_P|^2 + P_H = 1 with the horizon flux from B^trans (Teukolsky-Press / S10):
       P_H = (2 M r_+)^5 k (k^2 + kap^2)(k^2 + 4 kap^2) |B^trans/B^inc|^2 / omega^5,
       k = omega - m Omega_H, kap = 2 pi T_H = (r_+ - r_-)/(4 M r_+).
  5. The normalization used in all previous checks, S_old = (-1)^(l+1) D/(16 omega^4) B^ref/B^inc (helicity-
     preserving, parity-averaged), equals (S_+ + S_-)/2 and is not the unit-flux one: |S_old/S_P| = D/|C|.
"""
import sys, mpmath as mp
from mst import MST, I
from phys import spheroidal_c
from csm import D_TS, secant

mp.mp.dps = 60


def C2_S10(lam, a, om, m, M=1):
    """|C|^2 exactly as written in Supplement S10."""
    return (((lam + 2)**2 + 4*a*m*om - 4*a**2*om**2)*(lam**2 + 36*a*m*om - 36*a**2*om**2)
            + (2*lam + 3)*(96*a**2*om**2 - 48*a*m*om) + 144*om**2*(M**2 - a**2))


def setup(l, m, q, eps, N=60):
    s = -2
    q, eps = mp.mpf(q), mp.mpf(eps)
    c = q*eps/2
    E, _ = spheroidal_c(s, l, m, c)
    lam = E - 2*m*c + c**2 - s*(s + 1)
    k0 = mp.sqrt(1 - q**2)
    P = MST(eps, q, m, s, lam, k0, logkappa=mp.log(k0), N=N)
    lm = mp.mpf(l)
    x0 = (-2 - s**2/(lm*(lm+1)) + ((lm+1)**2 - s**2)**2/((2*lm+1)*(2*lm+2)*(2*lm+3))
          - (lm**2 - s**2)**2/((2*lm-1)*2*lm*(2*lm+1)))/(2*lm+1)
    Gf = lambda x: P.g(l + x*eps**2)/eps**2
    x = secant(Gf, x0, mp.mpf(10)**(-50))
    nu = l + x*eps**2
    return dict(s=s, q=q, eps=eps, c=c, lam=lam, k0=k0, nu=nu, P=P, gres=abs(P.g(nu)), N=N, m=m, l=l)


def floquet_B(P, nu, a, s, e, k):
    """incoming / outgoing asymptotic coefficients of R_C^nu built from coefficients a_n."""
    om = e/2
    ph = e*mp.log(e) - (1 - k)*e/2
    sum_a = mp.fsum(a[n] for n in a)
    sum_b = mp.fsum((-1)**n*mp.rf(nu+1+s-I*e, n)/mp.rf(nu+1-s+I*e, n)*a[n] for n in a)
    Ap = (mp.exp(-mp.pi*e/2)*mp.exp(I*mp.pi*(nu+1-s)/2)*mp.power(2, -1+s-I*e)
          * mp.gamma(nu+1-s+I*e)/mp.gamma(nu+1+s-I*e)*sum_a)
    Am = mp.power(2, -1-s+I*e)*mp.exp(-I*mp.pi*(nu+1+s)/2)*mp.exp(-mp.pi*e/2)*sum_b
    return Ap*mp.exp(-I*ph)/om, Am*mp.exp(I*ph)/om**(1+2*s)


def branch(d, sign=+1, eta=+1):
    """MST data for the physical solution (sign=+1) or its image (sign=-1, path eta)."""
    q, eps, m, s, l = d['q'], d['eps'], d['m'], d['s'], d['l']
    k0 = d['k0']
    kappa, logk = (k0, mp.log(k0)) if sign > 0 else (-k0, mp.log(k0) + I*mp.pi*eta)
    P = MST(eps, q, m, s, d['lam'], kappa, logkappa=logk, N=d['N'])
    P.nu = d['nu']
    A = P.amplitudes()
    return P, A


def analyse(l, m, q, eps, verbose=True):
    d = setup(l, m, q, eps)
    s, e, nu, lam, c = d['s'], d['eps'], d['nu'], d['lam'], d['c']
    om, k0 = e/2, d['k0']
    P, A = branch(d, +1)
    a = A['a']
    aneg = {n: a[-n] for n in a}
    Binc_nu, Bref_nu = floquet_B(P, nu, a, s, e, k0)
    Binc_mn, Bref_mn = floquet_B(P, -nu-1, aneg, s, e, k0)
    rho = A['rho']
    id_inc = abs(Binc_mn/Binc_nu/(-I*mp.exp(-I*mp.pi*nu)*rho) - 1)
    id_ref = abs(Bref_mn/Bref_nu/(I*mp.exp(I*mp.pi*nu)) - 1)
    Kp, Km = A['Kp'], A['Km']
    Binc = Kp*Binc_nu + Km*Binc_mn
    Bref = Kp*Bref_nu + Km*Bref_mn
    D = D_TS(lam, c, m)
    C2 = C2_S10(lam, d['q'], om, m)
    out = dict(l=l, m=m, q=q, eps=eps, nu=nu, gres=d['gres'],
               id_inc=id_inc, id_ref=id_ref, C2_check=abs(C2/(D**2 + 144*om**2) - 1))
    X = Km/Kp
    frac = (1 + I*mp.exp(I*mp.pi*nu)*X)/(1 - I*mp.exp(-I*mp.pi*nu)*rho*X)
    rp, rm = 1 + k0, 1 - k0
    OmH = d['q']/(2*rp)
    kH = om - m*OmH
    kap = (rp - rm)/(4*rp)
    PH = (2*rp)**5*kH*(kH**2 + kap**2)*(kH**2 + 4*kap**2)*abs(A['Btr']/Binc)**2/om**5
    res = {}
    for Par in (+1, -1):
        CP = D + 12*I*Par*om
        F = (-1)**(l+1)*CP/(16*om**4)*Bref_nu/Binc_nu
        S_fact = F*frac
        S_dir = (-1)**(l+1)*CP/(16*om**4)*Bref/Binc
        res[Par] = dict(F=F, S=S_fact, fact_vs_dir=abs(S_fact/S_dir - 1),
                        flux=abs(S_fact)**2 + PH - 1)
    S_old = (-1)**(l+1)*D/(16*om**4)*Bref/Binc
    out.update(PH=PH, kH=kH, S_P=res[1]['S'], absS=abs(res[1]['S']),
               fact_vs_dir=max(res[1]['fact_vs_dir'], res[-1]['fact_vs_dir']),
               flux_res=max(abs(res[1]['flux']), abs(res[-1]['flux'])),
               flux_old=abs(S_old)**2 + PH - 1,
               old_over_unit=abs(S_old/res[1]['S']), D_over_C=D/mp.sqrt(D**2 + 144*om**2),
               avg_check=abs((res[1]['S'] + res[-1]['S'])/2/S_old - 1),
               flip_check=abs(abs(S_old)**2 + abs((res[1]['S'] - res[-1]['S'])/2)**2 - abs(res[1]['S'])**2),
               F_abs=abs(res[1]['F']), d=d, P=P, A=A, D=D)
    return out


if __name__ == '__main__':
    l = 2
    print("Flux balance and factorization, s=-2, l=2 (60 digits).  Columns:")
    print(" m    q     eps   |  g(nu)   A-identities  S_fact/S_dir-1 | |C|^2(S10)/(D^2+144w^2)-1 |"
          " |S_P|^2+P_H-1   P_H            | old norm: |S_old|^2+P_H-1, |S_old/S_P| - D/|C|")
    # eps <= 0.4: nu is found from the low-frequency guess; larger eps: highfreq_check.py, complexnu_check.py
    cases = [(m, q, e) for (q, e) in [(0, '0.01'), (0, '0.3'), ('0.3', '0.05'), ('0.7', '0.2'),
                                      ('0.95', '0.4')] for m in (2, 1, 0, -1, -2)]
    for m, q, e in cases:
        if q == 0 and m != 2:
            continue
        o = analyse(l, m, q, e)
        print(f"{m:+d} {mp.nstr(mp.mpf(q),3):>5} {mp.nstr(mp.mpf(e),3):>5} | {mp.nstr(o['gres'],2):>8} "
              f"{mp.nstr(max(o['id_inc'], o['id_ref']),2):>8} {mp.nstr(o['fact_vs_dir'],2):>8} | "
              f"{mp.nstr(o['C2_check'],2):>8} | {mp.nstr(o['flux_res'],2):>8}  {mp.nstr(o['PH'],8):>15}"
              f" {'(superradiant)' if o['kH'] < 0 else '':14} | {mp.nstr(o['flux_old'],5):>12}"
              f"  {mp.nstr(o['old_over_unit'] - o['D_over_C'],2):>8}")
        sys.stdout.flush()
