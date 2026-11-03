# conservative_test: polar classical amplitude and the closed form (18)

Classical polar Compton amplitude of super-extremal Kerr (incidence along the spin), s = −2,
helicity preserving, first order in G. The classical limit is taken directly: ε = 2GMω = 10⁻¹⁶ at
fixed c = aω, full MST (near and far zone) at 110 digits, both continuations
κ = η i q √(1 − 1/q²) with q = 2c/ε, partial waves ℓ = 2…19, 28 Chebyshev nodes |c| ≤ 0.2.

Paper: Eqs. (17), (18); Supplemental S4, Eqs. (S12)–(S14).

## Result

Through O(a¹³) the data agree with

```math
f^{(\eta)}(\Theta) = f_0\, e^{-2cv}\Big[1 + (1-\eta)\,\frac{64c^4}{3}\int_0^1 dt\, t^3 e^{4ct}\big(e^{4cv(1-t)}-1\big)\Big],
```

with f₀ = M cos⁴(Θ/2)/sin²(Θ/2) and v = sin²(Θ/2). The relative deviation is 3×10⁻¹⁶ through O(a¹⁰)
and 10⁻¹⁵ through O(a¹³), the size of the finite-ε systematic. For η = +1 the amplitude is the pure
exponential exp(−2cv); the exchange-odd part (η-odd) is Eq. (18).

## Scripts

| Script | Purpose | Time |
|---|---|---|
| `run_polar2.py` | Production run, two processes → `polar_run2.json` (included) | ~15 min |
| `analyze3.py` | Projections F_L(c) onto ₋₂Y_L2 and their Taylor coefficients in c → `taylor2.json` (included) | seconds |
| `fitlib.py` | Projections of the exponential and of the contact terms; least-squares fit per order (module) | — |
| `fit_contact2.py` | Order-by-order contact polynomials for η = +1 and η = −1; the η = +1 residual is at the noise level (pure exponential) | seconds |
| `closedform_check.py` | The closed form above against the data, order by order through O(a¹³) | seconds |
| `higher_orders_eq17.py` | Expansion of the closed form (SymPy): the a⁸ term of Eq. (17), the coefficients b₀–b₉ of Eq. (S14), the elementary form (S13), the factor 3 + cos Θ of the odd orders | ~30 s |
| `mst.py`, `polar.py`, `csm.py` | Core modules (copies) | — |

## Expected output

`closedform_check.py` prints, for n = 5…13 and both η, the maximal deviation over the projections,
for example

```
n= 7 eta=-1  L<=9:  max deviation  1.02e-15  (data scale 3.47, rel 2.95e-16)
n=13 eta=-1  L<=7:  max deviation  3.45e-17  (data scale 0.128, rel 2.69e-16)
```

The η = +1 data at n ≥ 11 are at the noise level of the computation, as expected for the pure exponential.

`higher_orders_eq17.py` ends with

```
odd b_k at x=-3: [0, 0, 0, 0, 0]
elementary form - series: 0
odd part of f^(eta) - fodd: 0
```

`polar_run2.json` is also read by `../segment_test/vs_production.py` and
`../conservative_continuation/se_taylor_l2.py`, and `taylor2.json` by
`../physical_region_test/analyze_phys.py`.
