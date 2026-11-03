import json, mpmath as mp
mp.mp.dps = 50
Tj = json.load(open('taylor2.json'))
T = {(int(k.split('_')[0]), int(k.split('_')[1])): [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in v] for k, v in Tj.items()}
Ls = sorted(set(L for (_, L) in T))
Y0 = lambda L: mp.sqrt(mp.mpf(2*L+1)/(4*mp.pi))
def YL(L, x): return mp.sqrt(mp.mpf(2*L+1)/(4*mp.pi))*((1+x)/2)**2*mp.jacobi(L-2, 0, 4, x)
J = {(L, k): 2*mp.pi*mp.quad(lambda x: YL(L, x)*((1+x)/2)**2*((1-x)/2)**k, [-1, 1]) for L in Ls for k in range(0, 13)}
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

