# conservative_continuation: the exchange-even ε⁶ function and 𝒞(q)

Exchange-even analogue of `../continuation_test/`: ℓ = 2, m = 2 (polar), s = −2.

Paper: Eqs. (12), (16); Supplemental S5, S9(e).

## Closed form

The exchange-even ε⁶ coefficient of the Teukolsky S-matrix, E₆(q), the coefficient of ε⁶ in
(S_outer + S_inner)/2 with q = a/GM, is

E₆(q) = polynomial(q, degree ≤ 5) + (4i/225) q (1 + 3q²) [ln(εκ) + Re ψ(1 + 2iq/κ)],   κ = √(1 − q²).

The ln ε piece is the ε⁶ ln ε term; everything else is the ε⁶ coefficient with that logarithm removed.
With F(q) = (4/225) q (1 + 3q²) [ln κ + Re ψ(1 + 2iq/κ)], the contact coefficient of Eq. (16) is
𝒞 = 160 F.

## Scripts

| Script | Purpose | Time |
|---|---|---|
| `phys_hi_p.py` | Physical-region coefficients up to q²⁴ with the ε⁶ ln ε term separated, from `../prl_checks/phys_l2_*.json` → `phys_e6_hi.json` (included) | seconds |
| `e6_fit.py` | Fit: the odd coefficients q⁷…q²³ are reproduced by the ψ + log ansatz with λ = −1/450 (naive log bookkeeping; the refined test is `../prl_checks/uniqueness_check.py`) | seconds |
| `e6_se.py`, `e6_se2.py`, `e6_se3.py` | Production: direct MST at super-extremal q = 1.5…5000 (fixed q, expansion in ε, three radii for the log term) → `e6_se*.json` (included) | not timed |
| `e6_check_se2.py` | Continued closed form vs direct super-extremal MST (10⁻⁶…10⁻⁸ with this bookkeeping) | seconds |
| `se_taylor_l2.py` | Classical (ε → 0 at fixed aω) ℓ = 2 Taylor coefficients, both branches, from `../conservative_test/polar_run2.json` → `se_taylor_l2.json` (included) | seconds |
| `e6_identities.py` | 160 [q⁵]F (Taylor) = 16b₂ of arXiv:2610.40101 exactly; 160 F/q⁵ → 64/15 at large q | seconds |
| `mst.py`, `polar.py`, `csm.py`, `phys.py`, `d6_se.py` | Core modules (copies; `d6_se.py` provides the MST driver used by `e6_se.py`) | — |

## Expected output

`e6_identities.py` prints

```
160 * [q^5] F (Taylor, physical)   = 2.537433225271721941600605
16 b2 (2610.40101, eq. 4.2)        = 2.537433225271721941600605
160 * F(q)/q^5 continued, q=1.0e+7: 4.26666666666808
64/15 = 4.26666666666667
```
