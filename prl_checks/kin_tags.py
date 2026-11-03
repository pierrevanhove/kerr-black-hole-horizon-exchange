"""Which per-m (l=2) patterns do BGKV's alpha-tagged and untagged a^5 terms correspond to?  (helicity +1, physical sign of w.a)"""
import io, contextlib, mpmath as mp, random
with contextlib.redirect_stdout(io.StringIO()):
    import kin_even as KE
mp.mp.dps = 30
f = mp.factorial
random.seed(5)
pts = [(mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0, 6.2))) for _ in range(10)]
def fit(w, kind):
    rows, rhs = [], []
    for (g, t, p) in pts:
        pref, wa = KE.structures(g, t, p, 1)
        x = (-wa) if kind == 'even' else 1          # physical w.a = -(code w.a); |a| w = 1
        rows.append([pref[j]*x for j in range(3)]); rhs.append(KE.bhpt(g, t, p, w))
    sol, res = mp.qr_solve(mp.matrix(rows), mp.matrix(rhs))
    return [mp.nstr(mp.re(s), 12) for s in sol], mp.nstr(res, 3)
wodd = {m: (-1)**m*f(2+m)*f(2-m) for m in range(-2, 3)}
print("eta-odd  full (-1)^m(2+m)!(2-m)!      ->", fit(wodd, 'odd'))
print("eta-odd  m != 0 part (coth pole)     ->", fit({m: (wodd[m] if m else 0) for m in range(-2, 3)}, 'odd'))
print("eta-odd  m = 0 part (polynomial eta) ->", fit({m: (wodd[m] if m == 0 else 0) for m in range(-2, 3)}, 'odd'))
weven = {m: -mp.sign(m)*wodd[m] for m in range(-2, 3)}
print("eta-even psi-pole part               ->", fit(weven, 'even'))
print("BGKV: c4 = eta*(alpha 64/15, alpha 16/5, 4/15);  c3 = (alpha 64/15, alpha 16/3, 4/15 + alpha 16/15)")
# is the untagged c3 = (0,0,4/15) structure a pure l=2 pattern?  fit it with the two antisymmetric l=2 patterns
rows, rhs = [], []
for (g, t, p) in pts:
    pref, wa = KE.structures(g, t, p, 1)
    rhs.append(pref[2]*(-wa)*mp.mpf(4)/15)
    rows.append([KE.bhpt(g, t, p, {m: (1 if m == 2 else -1 if m == -2 else 0) for m in range(-2, 3)}),
                 KE.bhpt(g, t, p, {m: (1 if m == 1 else -1 if m == -1 else 0) for m in range(-2, 3)}),
                 KE.bhpt(g, t, p, {m: (1 if m == 0 else 0) for m in range(-2, 3)})])
sol, res = mp.qr_solve(mp.matrix(rows), mp.matrix(rhs))
print("untagged c3^(2)=4/15 as l=2 pattern (m=2 antisym, m=1 antisym, m=0):", [mp.nstr(s, 10) for s in sol], " residual", mp.nstr(res, 3))
