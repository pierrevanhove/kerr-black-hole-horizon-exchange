import mpmath as mp
mp.mp.dps = 50
from mst import MST, X_leading, I

def lam_s0(l, m, q, eps):
    c = mp.mpf(q)*mp.mpf(eps)/2
    E = l*(l+1) - c**2*(2*l*(l+1)-2*m**2-1)/((2*l-1)*(2*l+3))
    return E - 2*mp.mpf(q)*m*mp.mpf(eps)/2 + c**2

def run(case, l, m, s, q, eps):
    kap = mp.sqrt(1-mp.mpf(q)**2)
    lam = lam_s0(l, m, q, eps) if s == 0 else l*(l+1) - s*(s+1)   # s!=0 only used with q=0
    out = {}
    P = MST(eps, q, m, s, lam, kap); P.solve_nu(l - 0.3*float(eps)**2); A = P.amplitudes()
    for eta in (+1, -1):
        D = MST(eps, q, m, s, lam, -kap, logkappa=mp.log(kap)+I*mp.pi*eta)
        D.solve_nu(P.nu)                      # start from the physical nu
        out[eta] = (D, D.amplitudes())
    d = P.nu - l
    print(f"\n[{case}] l={l} m={m} s={s} q={q} eps={eps}")
    print("  nu_phys  =", mp.nstr(P.nu.real, 20), "   nu_dual - nu_phys =", mp.nstr(abs(out[1][0].nu-P.nu), 3))
    print("  |Sfar_dual/Sfar - 1| =", mp.nstr(abs(out[1][1]['Sfar']/A['Sfar']-1), 3))
    for eta in (+1, -1):
        D, B = out[eta]
        rX = B['X']/A['X'] - 1
        lead = -mp.exp(2*I*mp.pi*eta*d)*mp.sin(mp.pi*(d+I*P.tau))/mp.sin(mp.pi*(d-I*P.tau)) - 1
        rH = B['H']/A['H'] - 1
        print(f"  eta={eta:+d}: X_dual/X-1 = {mp.nstr(rX,8)} | leading formula {mp.nstr(lead,8)} | H_dual/H-1 = {mp.nstr(rH,8)}")
    return P, A, out

for eps in ['0.1', '0.05', '0.025']:
    run("Schwarzschild s=0", 2, 0, 0, 0, eps)
for eps in ['0.1', '0.05', '0.025']:
    run("Schwarzschild s=-2", 2, 2, -2, 0, eps)
for eps in ['0.1', '0.05', '0.025']:
    run("Kerr s=0", 2, 2, 0, '0.6', eps)
