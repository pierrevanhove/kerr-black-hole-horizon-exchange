"""Higher orders of Eq. (17): expansion of the closed form of the polar eta-odd amplitude, Eqs. (18), (S12),
   f^odd/f0 = -(64/3) c^4 e^{-2cv} Int_0^1 dt t^3 e^{4ct} (e^{4cv(1-t)} - 1),   c = a omega, v = sin^2(Theta/2), x = cos Theta,
written as  f^odd/f0 = -(64/15) v c^5 sum_k b_k(x) c^k  (b_0..b_3: Eq. (17); b_4..b_8: Eq. (S14)).
Also checks the elementary form, Eq. (S13), with Q(b) = b^3 - 3b^2 + 6b - 6, that the eta-odd part of Eq. (S12) is
Eq. (18), and that the odd-order b_k vanish at x = -3 (factor 3 + cos Theta).  The closed form itself is validated
against the full-MST data by closedform_check.py (through a^13)."""
import sympy as sp
c, v, t, x = sp.symbols('c v t x')
N = 14
integrand = sp.series(t**3*sp.exp(4*c*t)*(sp.exp(4*c*v*(1 - t)) - 1), c, 0, N).removeO()
G = sp.Rational(64, 3)*c**4*sp.integrate(sp.expand(integrand), (t, 0, 1))
fodd = sp.expand(sp.series(-sp.exp(-2*c*v)*G, c, 0, N + 1).removeO())
lead = -sp.Rational(64, 15)*v*c**5
br = sp.expand(sp.simplify(fodd/lead))
bs = []
for k in range(0, N - 4):
    bk = sp.factor(sp.expand(br.coeff(c, k).subs(v, (1 - x)/2)))
    bs.append(bk)
    print(f"b_{k}(x) =", bk, "   =", sp.expand(bk))
# odd-order coefficients vanish at x = -3 (v = 2), where the closed form is odd in c
print("odd b_k at x=-3:", [sp.simplify(bs[k].subs(x, -3)) for k in range(1, len(bs), 2)])
# elementary closed form and its check
u = 1 - v
Q = lambda b: b**3 - 3*b**2 + 6*b - 6
closed = -sp.exp(-2*c*v)/12*((sp.exp(4*c)*Q(4*c*u) + 6*sp.exp(4*c*v))/u**4 - sp.exp(4*c)*Q(4*c) - 6)
diff = sp.simplify(sp.expand(sp.series(closed, c, 0, N + 1).removeO()) - fodd)
print("elementary form - series:", diff)
# eta = +1 amplitude is e^{-2cv}; full eta form
eta = sp.symbols('eta')
feta = sp.exp(-2*c*v)*(1 + (1 - eta)*G)
fodd2 = sp.expand(sp.series((feta.subs(eta, 1) - feta.subs(eta, -1))/2, c, 0, N + 1).removeO())
print("odd part of f^(eta) - fodd:", sp.simplify(fodd2 - fodd))
# Taylor coefficients in v (as quoted in S4)
for n in range(5, 9):
    print(f"[c^{n}] f^odd/f0 =", sp.factor(fodd.coeff(c, n)))
