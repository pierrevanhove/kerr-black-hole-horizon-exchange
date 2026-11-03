"""Exact O(c^2) ingredients for the a^7 eta-odd prediction (s=-2).
1) R(c,m) = X_{-nu-1}/X_nu in the classical limit (l=2): series from the limiting MST recurrence.
2) Spin-weighted spheroidal harmonic S_22(theta;c) to O(c^2) (exact perturbation theory, basis l'=2..6).
3) Polar bracket: H_odd(2,2) S_22(0)S_22(Theta) + H_odd(3,2) Y_32(0)Y_32(Theta), normalised to the a^5 term."""
import sympy as sp
c, x, m = sp.symbols('c x m')
s, l = -2, 2
E = l*(l+1) - 2*c*m*s**2/(l*(l+1))                     # spheroidal eigenvalue to O(c) (enough for O(c^2) here)
def al(n): return -2*c*(n+l+1+s)**2*(n+l+1-m)/((n+l+1)*(2*n+2*l+3))
def be(n): return (n+l)*(n+l+1) - E - c**2 - 2*m*s**2*c/((n+l)*(n+l+1))
def ga(n): return 2*c*(n+l-s)**2*(n+l+m)/((n+l)*(2*n+2*l-1))
# minimal solution on the positive side: R_n = -ga_n/(be_n + al_n R_{n+1}), truncated deep enough for O(c^2)
R = {6: 0}
for n in range(5, 0, -1):
    R[n] = -ga(n)/(be(n) + al(n)*R[n+1])
a = {0: sp.Integer(1)}
for n in range(1, 5):
    a[n] = sp.series(R[n]*a[n-1], c, 0, 3).removeO()
Nplus = sum((-1)**n*sp.rf(3-m, n)/sp.rf(3+m, n)*a[n] for n in range(0, 5))   # (2l+1)_n/n! (l+1+s)_n/(l+1-s)_n = 1 for l=2,s=-2
Dminus = sum(a[n] for n in range(0, 5))                                     # weights = 1 for l=2, s=-2
Rser = sp.simplify(sp.series(1/(Nplus*Dminus), c, 0, 3).removeO())
print("R(c,m) =", sp.collect(sp.expand(Rser), c))
print("   check vs 1 + (10m/9)c + 5(475m^2+261)c^2/3969 :", sp.simplify(Rser - (1 + sp.Rational(10,9)*m*c + 5*(475*m**2+261)*c**2/3969)))

# ---- spheroidal harmonic S_{l=2,m=2}(theta;c) to O(c^2), s=-2 ----
def Yl(L):   # spin -2, m=2 harmonic as polynomial in x=cos(theta): sqrt((2L+1)/4pi) ((1+x)/2)^2 P^{(0,4)}_{L-2}(x)
    return sp.sqrt(sp.Rational(2*L+1, 4)/sp.pi)*((1+x)/2)**2*sp.jacobi(L-2, 0, 4, x)
Ls = list(range(2, 7))
def mel(f, L1, L2): return sp.simplify(2*sp.pi*sp.integrate(sp.expand(Yl(L1)*Yl(L2)*f), (x, -1, 1)))
X1 = sp.Matrix(len(Ls), len(Ls), lambda i, j: mel(x, Ls[i], Ls[j]))
X2 = sp.Matrix(len(Ls), len(Ls), lambda i, j: mel(x**2, Ls[i], Ls[j]))
print("norm check:", [mel(1, L, L) for L in Ls])
E0 = [L*(L+1) for L in Ls]
V1, V2 = 2*s*X1, -X2                                    # H = diag(L(L+1)) + c V1 + c^2 V2 ; eigenvalue E (see text)
n0 = 0
coef = {}
for k in range(len(Ls)):
    if k == n0: continue
    d = E0[n0]-E0[k]
    first = V1[k, n0]/d
    second = (V2[k, n0]/d + sum(V1[k, j]*V1[j, n0]/(d*(E0[n0]-E0[j])) for j in range(len(Ls)) if j != n0)
              - V1[n0, n0]*V1[k, n0]/d**2)
    coef[k] = c*first + c**2*second
norm2 = sum((V1[k, n0]/(E0[n0]-E0[k]))**2 for k in range(len(Ls)) if k != n0)
S = (1 - c**2*norm2/2)*Yl(2) + sum(coef[k]*Yl(Ls[k]) for k in coef)
ratio = sp.series(sp.simplify(S*S.subs(x, 1)/(Yl(2)*Yl(2).subs(x, 1))), c, 0, 3).removeO()
ratio = sp.expand(sp.simplify(ratio))
print("S22(0;c)S22(Theta;c)/(Y22 Y22) =", sp.collect(ratio, c))

# ---- assemble polar bracket ----
R2 = (1 + sp.Rational(10,9)*m*c + 5*(475*m**2+261)*c**2/3969).subs(m, 2)
l3 = -sp.Rational(1, 567)*c**2*(3*x - 2)                  # l=3 leading term relative to -(64/15) c^5 sin^2(Theta/2)
br = sp.expand(sp.series(R2*ratio, c, 0, 3).removeO() + l3)
print("\npolar bracket  [1 + c(...) + c^2(...)] with x = cos(Theta):")
print("  O(c):  ", sp.factor(br.coeff(c, 1)))
print("  O(c^2):", sp.collect(sp.expand(br.coeff(c, 2)), x))
