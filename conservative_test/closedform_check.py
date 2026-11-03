"""Test the conjectured closed form against the full-MST polar Taylor data, order by order.
   Conjecture:  f^(eta) = f0 e^{-2cv} { 1 + (1-eta) (64/3) c^4 int_0^1 dt t^3 e^{4ct} (e^{4cv(1-t)} - 1) },
   u = cos^2(Th/2), v = sin^2(Th/2), c = a omega (orientation with the exponential e^{-2cv})."""
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, '.')
import sympy as sp, mpmath as mp
from fitlib import T, J, Ls, Y0
mp.mp.dps = 50
c, v, t = sp.symbols('c v t')
NMAX = 13
# series of (64/3) c^4 e^{-2cv} int t^3 e^{4ct}(e^{4cv(1-t)}-1) dt / v   (= sum_n sum_j qcons_{n,j} v^j c^n)
integrand = sp.series(t**3*sp.exp(4*c*t)*(sp.exp(4*c*v*(1-t)) - 1), c, 0, NMAX - 3).removeO()
I = sp.expand(sp.integrate(sp.expand(integrand), (t, 0, 1)))
G = sp.expand(sp.series(sp.Rational(64, 3)*c**4*sp.exp(-2*c*v)*I/v, c, 0, NMAX + 1).removeO())
qcons = {n: [G.coeff(c, n).coeff(v, j) for j in range(n)] for n in range(5, NMAX + 1)}
for n in range(5, 9):
    print(f"n={n}: 2*qcons =", [2*q for q in qcons[n] if q != 0])
print()
for n in range(5, NMAX + 1):
    Luse = [L for L in Ls if L <= 19 - n + 1] if n >= 11 else Ls
    for eta in (1, -1):
        res = []
        for L in Luse:
            E = (-2)**n/mp.factorial(n)*J[(L, n-1)]
            Pc = (1 - eta)*sum(mp.mpf(sp.Rational(q).p)/sp.Rational(q).q*J[(L, j)] for j, q in enumerate(qcons[n]))
            res.append(mp.im(T[(eta, L)][n]) - (E + Pc)/(4*mp.pi))
        Cn = mp.fsum(r*Y0(L) for r, L in zip(res, Luse))/mp.fsum(Y0(L)**2 for L in Luse)
        dev = max(abs(r - Cn*Y0(L)) for r, L in zip(res, Luse))
        scale = max(abs(mp.im(T[(eta, L)][n])) for L in Luse)
        print(f"n={n:2d} eta={eta:+d}  L<={Luse[-1]}:  max deviation {mp.nstr(dev, 3):>9}  (data scale {mp.nstr(scale, 3)}, rel {mp.nstr(dev/scale, 3)})")
