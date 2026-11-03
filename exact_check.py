import mpmath as mp
mp.mp.dps = 60
from mst import MST, I

def check(l, m, s, q, eps, lam, N=60):
    kap = mp.sqrt(1-mp.mpf(q)**2)
    P = MST(eps, q, m, s, lam, kap, N=N); P.solve_nu(l - 0.3*float(eps)**2); A = P.amplitudes()
    for eta in (+1, -1):
        D = MST(eps, q, m, s, lam, -kap, logkappa=mp.log(kap)+I*mp.pi*eta, N=N); D.solve_nu(P.nu)
        B = D.amplitudes()
        pred = mp.exp(I*mp.pi*eta*(2*P.nu+1))*mp.sin(mp.pi*(P.nu-I*P.tau))/mp.sin(mp.pi*(P.nu+I*P.tau))
        print(f"l={l} m={m} s={s} q={q} eps={eps} eta={eta:+d}:  |X_dual/X / formula - 1| = {mp.nstr(abs(B['X']/A['X']/pred-1),3)}"
              f"   (|X_dual/X - 1| = {mp.nstr(abs(B['X']/A['X']-1),3)})")

check(2, 0, 0, 0, '0.3', 6)
check(2, 2, -2, 0, '0.5', 4)
check(2, 2, -2, '0.7', '0.3', 3.7)     # Kerr, arbitrary (real) separation constant
check(3, 1, -1, '0.4', '0.6', 10.3)
check(1, -1, 0, '0.9', '0.2', 2.1)
