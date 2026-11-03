# Derive (not just verify) the nine eta-odd coefficients: least-squares fit of the Compton-ansatz structures
# to the BHPT eta-odd amplitude at many random kinematic points (helicity +1).
import mpmath as mp, random
mp.mp.dps = 30
from kinematics import bhpt, ansatz, table
I = mp.mpc(0, 1)
random.seed(11)
pts = [(mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0, 6.2))) for _ in range(12)]

def basis(order):
    names = ['c4'] if order == 5 else ['c9', 'c10']
    cols = []
    for nm in names:
        for mm in range(3):
            c = {k: [0, 0, 0] for k in ('c4', 'c9', 'c10')}; c[nm][mm] = 1
            cols.append((f"{nm}^({mm})", c))
    return cols

for order, sign in ((5, 1), (6, -1)):       # a^6 carries the dot-product sign convention found above
    cols = basis(order)
    rows, rhs = [], []
    for (g, t, p) in pts:
        y = bhpt(g, t, p, order)[0]*sign
        xs = [ansatz(g, t, p, order, +1, c) for _, c in cols]
        for part in (mp.re, mp.im):           # real least squares on Re and Im parts
            rows.append([part(x) for x in xs]); rhs.append(part(y))
    A = mp.matrix(rows); bvec = mp.matrix(rhs)
    sol = mp.lu_solve(A.T*A, A.T*bvec)
    resid = mp.norm(A*sol - bvec)
    print(f"O(a^{order}): residual {mp.nstr(resid, 3)}")
    for (nm, _), v in zip(cols, sol):
        key, mm = nm.split('^')[0], int(nm[-2])
        print(f"   {nm:8s} fitted = {mp.nstr(v, 20):>24s}   table = {mp.nstr(table[key][mm], 20)}")

# non-triviality: perturb one table coefficient and show agreement breaks
pert = {k: list(v) for k, v in table.items()}; pert['c4'][1] += mp.mpf('0.01')
g, t, p = pts[0]
print("\nperturbed c4^(1) by 0.01 -> BHPT/ansatz =", mp.nstr(bhpt(g, t, p, 5)[0]/ansatz(g, t, p, 5, +1, pert), 12))
