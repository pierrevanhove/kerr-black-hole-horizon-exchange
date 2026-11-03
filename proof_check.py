import sympy as sp
lam, a, b, X, Y = sp.symbols('lambda a b X Y', nonzero=True)
D = sp.diag(lam, 1/lam)                                   # big-loop monodromy in the Floquet basis (R_C^nu, R_C^{-nu-1})
P = sp.Matrix([[1, 1], [X, Y]])                           # columns: in-going (eig a) and out-going (eig 1/a) at r+
Mp = P*sp.diag(a, 1/a)*P.inv()
for name, Mm in [("M- = M_inf M+^-1", D*Mp.inv()), ("M- = M+^-1 M_inf", Mp.inv()*D)]:
    Ysol = sp.solve(sp.Eq(sp.simplify(Mm.trace()), b + 1/b), Y)          # spectrum of M- fixed by local exponents at r-
    for Ys in Ysol:
        if sp.simplify(Ys - X) == 0:                     # degenerate branch (P singular) -> discard
            continue
        M = sp.simplify(Mm.subs(Y, Ys))
        # eigenvector (1, Xt) of M- with eigenvalue b  (first row of (M - b) v = 0)
        Xt = sp.solve(sp.Eq(M[0, 0] - b + M[0, 1]*sp.Symbol('t'), 0), sp.Symbol('t'))[0]
        ratio = sp.factor(sp.simplify(Xt/X))
        print(f"{name}:  X_-/X_+ = {ratio}")
        # check it equals -lam^eta (lam a - b)/(lam b - a) for eta = -1 or +1
        for eta in (-1, 1):
            if sp.simplify(ratio + lam**eta*(lam*a - b)/(lam*b - a)) == 0:
                print(f"   = -lam^{eta:+d} (lam a - b)/(lam b - a)")
# translate: a = e^{2 pi eps+}, b = e^{2 pi eps-}, lam = e^{2 pi i nu}, tau = eps+ - eps-
nu, ep, em = sp.symbols('nu epsilon_p epsilon_m', real=True)
expr = (sp.exp(2*sp.pi*sp.I*nu)*sp.exp(2*sp.pi*ep) - sp.exp(2*sp.pi*em))/(sp.exp(2*sp.pi*sp.I*nu)*sp.exp(2*sp.pi*em) - sp.exp(2*sp.pi*ep))
target = sp.sin(sp.pi*(nu - sp.I*(ep-em)))/sp.sin(sp.pi*(nu + sp.I*(ep-em)))
vals = {nu: sp.Rational(37,17), ep: sp.Rational(3,10), em: -sp.Rational(7,13)}
print("check (lam a - b)/(lam b - a) == sin pi(nu - i tau)/sin pi(nu + i tau):", sp.N((expr - target).subs(vals), 30))
