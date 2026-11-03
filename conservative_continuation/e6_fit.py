import json, mpmath as mp
mp.mp.dps = 50
raw = json.load(open('phys_e6_hi.json'))
A = {tuple(int(t) for t in k.split('_')): (mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])), mp.mpc(mp.mpf(v[2]), mp.mpf(v[3]))) for k, v in raw.items()}
d = {p: mp.im((A[(1, p)][0] + A[(-1, p)][0])/2) for p in range(25)}       # Im of the self-dual eps^6 coefficient
# Taylor coefficients of the candidate functions (Cauchy on a small circle)
N, r = 128, mp.mpf('0.35')
qs = [r*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
def tay(f, P=25):
    v = [f(q) for q in qs]
    return [mp.re(mp.fsum(vv*q**(-p) for vv, q in zip(v, qs))/N) for p in range(P)]
kap = lambda q: mp.sqrt(1 - q**2)
L = tay(lambda q: mp.mpf(2)/225*q*(1 + 3*q**2)*mp.log(1 - q**2))
Psi = tay(lambda q: -8*q*(1 + 3*q**2)*(mp.digamma(1 - 2j*q/kap(q)) + mp.digamma(1 + 2j*q/kap(q)))/2)
# fit d_p (odd p = 1..23) = a1 [p=1] + a3 [p=3] + a5 [p=5] + L_p + lam*Psi_p
odd = list(range(1, 24, 2))
rows = [[1 if p == 1 else 0, 1 if p == 3 else 0, 1 if p == 5 else 0, Psi[p]] for p in odd]
rhs = [d[p] - L[p] for p in odd]
M, b = mp.matrix(rows), mp.matrix(rhs)
sol = mp.qr_solve(M, b)[0]
a1, a3, a5, lam = [sol[i] for i in range(4)]
print("fit: a1 =", mp.nstr(a1, 15), " a3 =", mp.nstr(a3, 15), " a5 =", mp.nstr(a5, 15), " lambda =", mp.nstr(lam, 15))
print(" p   data Im E6               model                      residual")
for p in odd:
    model = (a1 if p == 1 else 0) + (a3 if p == 3 else 0) + (a5 if p == 5 else 0) + L[p] + lam*Psi[p]
    print(f"{p:2d}  {mp.nstr(d[p], 15):>22}  {mp.nstr(model, 15):>22}  {mp.nstr(d[p] - model, 3):>10}")
print("lambda * 450 =", mp.nstr(lam*450, 12), "   lambda*900 =", mp.nstr(lam*900, 12))
for nm, val in (('a1', a1), ('a3', a3), ('a5', a5)):
    print(f"{nm} = {mp.nstr(val, 15)}  identify: {mp.identify(val, ['zeta(3)','euler','log(2)'])}")
json.dump({'a1': str(a1), 'a3': str(a3), 'a5': str(a5), 'lam': str(lam)}, open('e6_fit.json', 'w'))
