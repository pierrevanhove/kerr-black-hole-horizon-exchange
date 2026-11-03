# continuation_test: the exchange-odd ε⁶ function and its continuation

Exchange-odd ℓ = 2 coefficient at order ε⁶, from the physical region to super-extremal spin.

Paper: Eq. (10); Supplemental S9(d), S9(f).

## Closed form

For s = −2, ℓ = 2 and all m (verified against full MST), the coefficient D₆⁽ᵐ⁾(q) of ε⁶ in
(S_outer − S_inner)/2 is

D₆⁽ᵐ⁾(q) = (π/900) κ⁵ τ₀(1 + τ₀²)(4 + τ₀²) [coth(πϖ) − η]

with q = a/GM, κ = √(1 − q²), τ₀ = −mq/κ and ϖ = mq/κ = mΩ_H/(2πT_H); for m = 0, D₆ = −κ⁵/225.
Written with the Bose factor 1/(exp(2πϖ) − 1), this is Eq. (10).

## Scripts

| Script | Purpose | Time |
|---|---|---|
| `d6_closed.py` | m = 2: closed form vs the eleven physical-region Taylor coefficients (`l2_rich.json`); the q⁶ coefficient equals ρ₂/64 of arXiv:2610.40101 | seconds |
| `d6_m.py` | m = ±1, ±2: closed form vs MST at q = 0.3 and 0.7 (two worker processes) | ~2 min |
| `d6_m0.py` | m = 0 check (two worker processes) | ~30 s |
| `d6_se.py` | Production: direct MST at super-extremal q = 1.5…12 (fixed q, expansion in ε) vs the continued closed form → `d6_se.json` (included) | minutes |
| `d6_general.py` | q⁶ coefficients vs ρ_m/64 for all m; large-q limit 32 D₆/q⁵ vs the exchange-odd a⁵ coefficients (−1)^m (8/225) i (2+m)!(2−m)! | seconds |
| `phys.py`, `mst.py`, `polar.py`, `csm.py` | Core modules (copies) | — |

## Expected output

`d6_closed.py` lists the Taylor coefficients of the data and of the closed form, with differences of order
10⁻¹⁸–10⁻¹⁹, and ends with

```
q^6 coefficient vs BBCPSS rho_2/64: 0.00621387874520450288 0.00621387874520450288
```

`d6_general.py` shows 32 D₆/q⁵ approaching the targets ±0.8533…i, ∓0.2133…i, 0.1422…i as
q = 10², 10³, 10⁴.
