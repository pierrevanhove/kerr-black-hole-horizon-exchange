"""Anti-self-dual eps^7 function D7(q) for l=2, m=2 in the physical region, and a closed-form test:
   D(eps,q) = -(eps^6/900) kappa^5 Gamma(3+i tau)Gamma(3-i tau) e^{pi tau} [1 + u(q) eps + ...],  tau = (eps - m q)/kappa
=> D7 = D6 [u(q) + (1/kappa) dlog/dtau0],  dlog/dtau = 1/tau + 2tau/(1+tau^2) + 2tau/(4+tau^2) + pi(1 - coth(pi tau))."""
import json, mpmath as mp
mp.mp.dps = 60
NE, NQ, RQ, m = 12, 40, mp.mpf('0.3'), 2
qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
def load(fn, RE):
    data = json.load(open(fn)); D = {(r['sign'], r['j']): r['pts'] for r in data}
    eps = [RE*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
    out = {}
    for sg in (1, -1):
        F = [[mp.mpc(mp.mpf(D[(sg, j)][k]['St'][0]), mp.mpf(D[(sg, j)][k]['St'][1])) for k in range(NQ)] for j in range(NE)]
        out[sg] = {n: [mp.fsum(F[j][k]*eps[j]**(-n) for j in range(NE))/NE for k in range(NQ)] for n in (6, 7)}
    return out
A, B = load('phys_l2_2e-5.json', mp.mpf('2e-5')), load('phys_l2_1e-5.json', mp.mpf('1e-5'))
def taylor(vals):
    return [mp.fsum(v*q**(-p) for v, q in zip(vals, qs))/NQ for p in range(NQ)]
D7a = taylor([(A[1][7][k] - A[-1][7][k])/2 for k in range(NQ)])
D7b = taylor([(B[1][7][k] - B[-1][7][k])/2 for k in range(NQ)])
D6b = taylor([(B[1][6][k] - B[-1][6][k])/2 for k in range(NQ)])
print("radius stability of the anti-self-dual eps^7 Taylor coefficients (r=2e-5 vs 1e-5):")
for p in range(0, 12):
    print(f"  p={p:2d}  {mp.nstr(D7b[p], 14):>40}   diff {mp.nstr(D7a[p]-D7b[p], 3)}")
# closed-form pieces (eta = +1 in the physical-region convention)
kap = lambda q: mp.sqrt(1 - q**2)
def D6f(q):
    k = kap(q); t = -m*q/k
    return -(mp.pi/900)*k**5*t*(1 + t**2)*(4 + t**2)*(1 + mp.coth(mp.pi*t))
def shiftf(q):
    k = kap(q); t = -m*q/k
    return D6f(q)/k*(1/t + 2*t/(1 + t**2) + 2*t/(4 + t**2) + mp.pi*(1 - mp.coth(mp.pi*t)))
N, r = 160, mp.mpf('0.3')
zs = [r*mp.expjpi(2*(mp.mpf(k) + mp.mpf(1)/3)/N) for k in range(N)]
def tay(f, P=30):
    v = [f(z) for z in zs]; return [mp.fsum(vv*z**(-p) for vv, z in zip(v, zs))/N for p in range(P)]
T6, Tsh = tay(D6f), tay(shiftf)
print("\ncheck D6 closed form vs data (p=0..8):", [mp.nstr(D6b[p] - T6[p], 3) for p in range(9)])
# fit  D7 = Tsh + sum_j u_j q^j D6 + sum_k v_k q^k   (u: degree <= 3, v: degree <= 7)
P = 26
qD6 = [[T6[p - j] if p >= j else 0 for p in range(P)] for j in range(4)]
cols = [qD6[j] for j in range(4)] + [[1 if p == k else 0 for p in range(P)] for k in range(8)]
names = [f"u{j}" for j in range(4)] + [f"v{k}" for k in range(8)]
M = mp.matrix([[cols[c][p] for c in range(len(cols))] for p in range(P)])
rhs = mp.matrix([D7b[p] - Tsh[p] for p in range(P)])
sol, res = mp.qr_solve(M, rhs)
print("\nfit D7 - (1/kappa) dD6/dtau0 = (u0+u1 q+u2 q^2+u3 q^3) D6 + poly_7(q):  residual", mp.nstr(res, 4))
for nm, v in zip(names, sol):
    print(f"   {nm} = {mp.nstr(v, 16)}")
print("   (5/9) m =", mp.nstr(mp.mpf(5)/9*m, 16))
