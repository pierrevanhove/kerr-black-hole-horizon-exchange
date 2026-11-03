"""Independent check of the O(a^7) eta-odd polar prediction, Eq. (prediction), against the full-MST polar amplitude
(closed form validated against 110-digit MST, l = 2..19, through O(a^13); see conservative_test/):
   f^(eta)/f0 = e^{-2cv} [1 + (1-eta) G],  G = (64/3) c^4 int_0^1 t^3 e^{4ct} (e^{4cv(1-t)} - 1) dt,  v = sin^2(Theta/2).
eta-odd part: f^odd/f0 = -e^{-2cv} G."""
import sympy as sp
c, v, t, x = sp.symbols('c v t x')
N = 8
integrand = sp.series(t**3*sp.exp(4*c*t)*(sp.exp(4*c*v*(1 - t)) - 1), c, 0, N).removeO()
G = sp.Rational(64, 3)*c**4*sp.integrate(sp.expand(integrand), (t, 0, 1))
fodd = sp.expand(sp.series(-sp.exp(-2*c*v)*G, c, 0, N).removeO())
pred = -sp.Rational(64, 15)*v*c**5*(1 + sp.Rational(2, 3)*c*(3 + x) + c**2*(sp.Rational(13, 6) + sp.Rational(29, 21)*x + sp.Rational(11, 42)*x**2))
pred = sp.expand(pred.subs(x, 1 - 2*v))
for n in (5, 6, 7):
    a, b = sp.factor(fodd.coeff(c, n)), sp.factor(pred.coeff(c, n))
    print(f"c^{n}:  MST closed form: {a}    prediction: {b}    difference: {sp.simplify(a - b)}")
