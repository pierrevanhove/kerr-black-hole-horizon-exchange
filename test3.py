import mpmath as mp
mp.mp.dps = 50
from mst import MST, I

def lam_s0(l, m, q, eps):
    c = mp.mpf(q)*mp.mpf(eps)/2
    E = l*(l+1) - c**2*(2*l*(l+1)-2*m**2-1)/((2*l-1)*(2*l+3))
    return E - 2*mp.mpf(q)*m*mp.mpf(eps)/2 + c**2

def branches(l, m, q, eps):
    kabs = mp.sqrt(mp.mpf(q)**2-1)
    lam = lam_s0(l, m, q, eps)
    P = MST(eps, q, m, 0, lam, I*kabs); P.solve_nu(l - 0.3*float(eps)**2); A = P.amplitudes()
    Q = MST(eps, q, m, 0, lam, -I*kabs); Q.solve_nu(P.nu); B = Q.amplitudes()
    return P, A, Q, B

print("== exactness of the swap factor for super-extremal spin (s=0) ==")
for (l, m, q, eps) in [(2, 2, '1.5', '0.2'), (2, -1, '3', '0.1'), (1, 1, '2', '0.3')]:
    P, A, Q, B = branches(l, m, q, eps)
    y = I*P.tau                               # real for super-extremal
    Phi = mp.exp(-I*mp.pi*(2*P.nu+1))*mp.sin(mp.pi*(P.nu-y))/mp.sin(mp.pi*(P.nu+y))
    print(f"l={l} m={m} q={q} eps={eps}: nu={mp.nstr(P.nu,12)} (|Im nu|={mp.nstr(abs(P.nu.imag),2)}), "
          f"|nu(+)-nu(-)|={mp.nstr(abs(P.nu-Q.nu),2)}, |X(-)/X(+)/Phi - 1| = {mp.nstr(abs(B['X']/A['X']/Phi-1),3)}, "
          f"|Sfar(-)/Sfar(+)-1| = {mp.nstr(abs(B['Sfar']/A['Sfar']-1),3)}")

print("\n== classical-like regime: a*omega fixed, eps -> 0 (q = 2 a omega/eps) ==")
aw = mp.mpf('0.3')
for (l, m) in [(2, 2), (2, 1), (2, -1), (2, -2)]:
    for eps in ['0.04', '0.02', '0.01']:
        q = 2*aw/mp.mpf(eps)
        P, A, Q, B = branches(l, m, q, eps)
        d = P.nu - l
        u = I*P.tau + m                       # detuning of y = i*tau from the integer -m
        odd, even = (A['X']-B['X'])/2, (A['X']+B['X'])/2
        Hodd, Heven = (A['H']-B['H'])/2, (A['H']+B['H'])/2
        print(f"l={l} m={m:+d} eps={eps} q={mp.nstr(q,4)}: delta={mp.nstr(d.real,6)}  u={mp.nstr(u.real,6)}  "
              f"X_odd/X_even = {mp.nstr(odd/even,6)}  delta/u = {mp.nstr((d/u).real,6)}  H_odd/H_even = {mp.nstr(Hodd/Heven,6)}")
