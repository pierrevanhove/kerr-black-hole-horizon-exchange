import json, mpmath as mp
mp.mp.dps = 50
Tj = json.load(open('taylor2.json'))
T = {(int(k.split('_')[0]), int(k.split('_')[1])): [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in v] for k, v in Tj.items()}
Ls = sorted(set(L for (_, L) in T))
Y0 = lambda L: mp.sqrt(mp.mpf(2*L+1)/(4*mp.pi))
def YL(L, x): return mp.sqrt(mp.mpf(2*L+1)/(4*mp.pi))*((1+x)/2)**2*mp.jacobi(L-2, 0, 4, x)
J = {(L, k): 2*mp.pi*mp.quad(lambda x: YL(L, x)*((1+x)/2)**2*((1-x)/2)**k, [-1, 1]) for L in Ls for k in range(0, 13)}
# normalization from n=1:  T_L^(1) = (i zeta/4pi)(-2 sigma) J_L^(0) + C1 Y0(L)
for eta in (1, -1):
    C1 = mp.im(T[(eta, 3)][1])/Y0(3)
    zs = mp.im(T[(eta, 2)][1] - 1j*C1*Y0(2))/(-2*J[(2, 0)]/(4*mp.pi))
    print(f"eta={eta:+d}: zeta*sigma = {mp.nstr(zs, 15)};  C1 = {mp.nstr(C1, 5)};  L>=4 residuals:",
          [mp.nstr(mp.im(T[(eta, L)][1]) - C1*Y0(L), 3) for L in Ls[2:]])

def fit(eta, n, deg, zeta=1, sigma=1, Luse=None):
    """least squares for q_{n,0..deg} and C_n:  Im T_L^(n) = (zeta/4pi)[E_L^(n) + sum_j q_j J_L^(j)] + C_n Y0(L)"""
    Luse = Luse or Ls
    rows, rhs = [], []
    for L in Luse:
        E = (-2*sigma)**n/mp.factorial(n)*J[(L, n-1)]
        rows.append([zeta*J[(L, j)]/(4*mp.pi) for j in range(deg+1)] + [Y0(L)])
        rhs.append(mp.im(T[(eta, L)][n]) - zeta*E/(4*mp.pi))
    A, b = mp.matrix(rows), mp.matrix(rhs)
    sol = mp.qr_solve(A, b)[0]
    res = A*sol - b
    return [sol[j] for j in range(deg+1)], sol[deg+1], max(abs(r) for r in res)

if __name__ == "__main__":
    for n in range(2, 9):
        for eta in (1, -1):
            deg = max(0, min(n-1, len(Ls)-2))
            q, C, r = fit(eta, n, deg)
            print(f"n={n} eta={eta:+d}: C_n={mp.nstr(C, 8):>14}  maxres={mp.nstr(r, 2):>8}  q=", [mp.nstr(v, 10) for v in q])

    # direct test of P(eta=+1)=0: residual after subtracting the exponential and a single forward constant
    print("\n eta=+1 direct test: max_L |Im T_L^(n) - E_L^(n)/4pi - C_n Y0(L)| relative to max_L |eta=-1 contact part|")
    for n in range(5, 12):
        q, C, r = fit(1, n, -1) if False else (None, None, None)
        rows = [[Y0(L)] for L in Ls]
        rp = [mp.im(T[(1, L)][n]) - (-2)**n/mp.factorial(n)*J[(L, n-1)]/(4*mp.pi) for L in Ls]
        rm = [mp.im(T[(-1, L)][n]) - (-2)**n/mp.factorial(n)*J[(L, n-1)]/(4*mp.pi) for L in Ls]
        Cp = mp.qr_solve(mp.matrix(rows), mp.matrix(rp))[0][0]
        Cm = mp.qr_solve(mp.matrix(rows), mp.matrix(rm))[0][0]
        resp = max(abs(rp[i] - Cp*Y0(L)) for i, L in enumerate(Ls))
        resm = max(abs(rm[i] - Cm*Y0(L)) for i, L in enumerate(Ls))
        print(f"  n={n:2d}: eta=+1 residual {mp.nstr(resp, 3):>10}   eta=-1 contact size {mp.nstr(resm, 3):>10}   ratio {mp.nstr(resp/resm, 3)}")
