# prl_checks: all-m functions, the ε⁷ operator, general kinematics and analytic control

ℓ = 2, s = −2 Kerr Teukolsky equation, M = 1, ε = 2Mω, q = a/M, κ = √(1 − q²), τ₀ = −mq/κ and

P⁽ᵐ⁾ = κ⁵ τ₀(1 + τ₀²)(4 + τ₀²)  [= −8q(1 + 3q²) for m = 2, −q(4 − 3q²) for m = 1, 0 for m = 0],

so that Π_m of Eq. (12) is −P⁽ᵐ⁾. All scripts run from this folder with Python 3 and mpmath (SymPy
for `a7_polar_check.py`). Except for the production driver `run_allm.py` (and the copy `d6_se.py`), they
read the included JSON files and finish in seconds.

## All-m exchange-even function and its continuation (Eq. (12); S5; S9(f))

| Script | Purpose and expected output |
|---|---|
| `run_allm.py` | Production (~20 min on two cores): physical-region 2D Cauchy integrals for m = 1, −1, −2, 0, both κ branches, two ε radii, and the super-extremal classical data for both η → `phys_allm.json`, `se_allm.json` |
| `allm_analyze.py` | Exchange-even ε⁶ function for every m, E₆ = poly₅(q) − (i/450)P⁽ᵐ⁾[ln(εκ) + Re ψ(1 + iτ₀)]: a free fit gives 450λ = −0.999999997 for both the ln κ and the Re ψ coefficients; for m = 0, E₆ is a degree-four polynomial. Writes `allm_phys_result.json` |
| `se_allm_analyze.py` | Continuation test for every m: the super-extremal exchange-even c⁵ coefficient equals 32i a₅⁽ᵐ⁾ plus the pole contribution pole_m = (64/75, −16/75, 0, 16/75, −64/75)i, to ~10⁻⁹. Also prints the per-branch c⁵ coefficients used in S5(d). Writes `se_taylor_allm.json` |
| `unified_check.py` | Physical branch only: after subtracting the closed form, Eq. (S18), the Taylor coefficients 6 ≤ p ≤ 20 are below 2×10⁻¹¹ for every m (S5, uniqueness (b)) |

## General kinematics (S6, Table S1)

| Script | Purpose and expected output |
|---|---|
| `kinematics.py` | Copy of the top-level script: exchange-odd BHPT amplitude vs the contact-term ansatz of Bautista, Guevara, Kavanagh and Vines at random kinematic points; ratios 1 to 10⁻³⁰ for helicity +1, failure for helicity −1 (negative control) |
| `kin_even.py` | Exchange-even a⁵ contact term in general kinematics: the ψ-pole contributions are the α-tagged terms c₃ = α(64/15, 16/3, 16/15), residual 5×10⁻³⁰ (helicity −1: negative control) |
| `kin_tags.py` | Exchange-odd: the m ≠ 0 (coth-pole) part gives the α-tagged c₄ = ηα(64/15, 16/5, 0); the m = 0 part gives the untagged c₄ = η(0, 0, 4/15). The untagged c₃⁽²⁾ = 4/15 is not a pure ℓ = 2 pattern |
| `kin_a6.py` | a⁶: the m ≠ 0 (pole) part equals the α-tagged c₉, c₁₀ (exchange-odd) and c₆, c₈, f₁ (exchange-even); the m = 0 part equals the untagged c₉⁽²⁾ = −8/45. All ratios are ±1 to 10⁻³⁰ (the sign is the code's a → −a convention for odd-degree structures) |

## The same mechanism at a⁶ (ε⁷); S9(g), S9(h); Eq. (S23)

| Script | Purpose and expected output |
|---|---|
| `se7_analyze.py` | ε⁷ coefficient at super-extremal q = 1.5…5000 (the ε⁶ ln ε leakage removed exactly; A + B ln r + C r with three radii), 64A₇/q⁶ for both sectors, and its extrapolation to q → ∞: 1.8964816 (exchange-even) and 1.8962963 (exchange-odd), against the classical values 1.89648161 and 256/135. The limits commute at a⁶ |
| `se7_compare.py` | The physical-region ε⁷ closed forms continued to super-extremal q reproduce the direct MST coefficients; for 50 ≤ q ≤ 5000 the agreement is 5×10⁻⁹ (exchange-even) and from 3×10⁻¹⁰ down to 5×10⁻¹³ (exchange-odd). For q ≤ 30 the direct data are less precise (10⁻⁴ to 10⁻⁵) |
| `d7_phys.py`, `d7_fit2.py` | Exchange-odd ε⁷ closed form for all m (residual ~10⁻¹²): D₇ = [(1/κ)∂_τ₀ + u₀ + (5m/9)q] D₆ with u₀ = π + i(2γ_E + 2 ln 2 − 37/12 − ν₂), ν₂ = −107/210; no q² term |
| `e7_phys.py`, `e7_allm.py` | Exchange-even ε⁷, all m: the same operator acting on T₆ = −(i/450)P[ln κ + Re ψ], plus a degree-six polynomial; π and the imaginary part of u₀ to 3×10⁻⁸, the coefficient of q equal to 5m/9 to 2×10⁻⁸. For m = 0 the ln ε coefficient is B₇ = −(2i/225)κ⁴ (the m = 0 ω⁷ logarithm of Saketh, Zhou and Ivanov, arXiv:2307.10391). Writes `e7_fit.json`, `e7_phys.json` |
| `a6_hyp.py` | Negative test: shifting only the singular function (trigamma or csch²) gives the wrong ratio |

The only term of the ε⁷ function that grows like q⁶ at large spin is (5m/9)q times the singular ε⁶
function; the τ-derivative term is regular at the pole. The a⁶ α-tagged terms are therefore the a⁵ pole
dressed by (1 + 10mc/9) plus the ℓ = 2/ℓ = 3 spheroidal mixing, with no new transcendental structure.

## Analytic control of the exchange-odd response (S7; S9(j))

| Script | Purpose and expected output |
|---|---|
| `analytic_odd.py` | Full MST: δ²𝒦₊/[(i/1800)ε⁷κ⁵τ(1 + τ²)(4 + τ²)] → 1 + (5mq/9)ε; (exp(iπν) + exp(−iπν)ρ)ε/(2iδ) → 1 + (π − iν₂)ε (to ~10⁻⁵) |
| `series_ratio_check.py` | Full MST at 60 digits: 𝒳 ratio = 1 + (5mq/9)ε + … (to ~10⁻⁶), ν = 2 − (107/210)ε² + (mq/12)ε³ for (m, q) = (2, 0.3), (1, 0.5), (−2, 0.4), (0, 0.6): ν₂ does not depend on q or m (derivation in S7(c)); ν₃ to 0.5% |
| `eps4_check.py` | The pole term (i/900)ε⁴Π_m/ν₂ of the near-zone factor is absent from S (its ε⁴ coefficient is cubic in q); spheroidal coefficients E₁ = −4m/3 and E₂ = (10/189)(9 − m²) − 1 (to 10⁻⁹) |
| `h_scaling.py` | Near-zone factor: exchange-even part ~ε⁴ (the pole term, cancelled by the far zone), exchange-odd part ~ε⁶; the far-zone factor is exchange-invariant to 10⁻⁷⁴ |
| `absorption_check.py` | Γ₊ − Γ₋ = −4ε⁶A₆⁻ (physical solution vs its image), including q = 0 (Page, 4/225); the linear extrapolation in ε agrees to ≤ 2.5×10⁻⁶ |

## Logarithms, uniqueness and the numbers of the paper (S5)

| Script | Purpose and expected output |
|---|---|
| `b6_precision.py` | B₆ = iΠ_m/450 for every m: naive two-radius fit 5×10⁻⁵ to 3×10⁻⁴, joint log bookkeeping 10⁻¹⁰ |
| `uniqueness_check.py` | m = 2 exchange-even ε⁶ polynomial fixed at \|q\| ≤ 0.3 and continued to 1.5 ≤ q ≤ 5000: agreement 4×10⁻¹² for q ≤ 30 (joint bookkeeping on both sides) and 3×10⁻¹⁰ to 6×10⁻¹³ for 50 ≤ q ≤ 5000; bounds \|c₆\| ≲ 3×10⁻¹⁸, \|c₇\| ≲ 6×10⁻²²; comparison of the naive and joint polynomials. The q⁵ coefficient −2.3457997302×10⁻⁵ i equals the η = +1 classical c⁵ coefficient printed by `se_allm_analyze.py`, divided by 32, to 1.5×10⁻¹⁵ (S5(d)) |
| `a7_polar_check.py` | SymPy: the exchange-odd part of the closed form (S12) has the c⁵, c⁶, c⁷ coefficients of Eq. (17) |
| `paper_numbers.py` | Re-derives the closed forms and identities quoted in the main text: the ε-constant i(2γ_E + 2 ln 2 − 37/12), u₀, Eq. (10) in its Bose and Boltzmann forms, π(τ + ε) = (ω − mΩ_H)/2T_H, 160[q⁵]F = 16b₂, [q⁶]D₆ = ρ_\|m\|/64 and the large-q limits, and the pole terms pole_m |

## Files

- `mst.py`, `csm.py`, `polar.py`, `phys.py`: core modules (copies). `d6_se.py`: copy of the production
  driver of `../continuation_test/`.
- Data: `phys_l2_1e-5.json`, `phys_l2_2e-5.json` (from `../physical_region_test/run_phys2.py`; also read by
  `../physical_region_test/rich_phys.py` and `../conservative_continuation/phys_hi_p.py`),
  `phys_allm.json`, `se_allm.json` (`run_allm.py`), `e6_se*.json`, `e6_fit.json`, `phys_e6_hi.json` (from
  `../conservative_continuation/`), `se_taylor_l2.json`, and the outputs of the analysis scripts
  (`allm_phys_result.json`, `se_taylor_allm.json`, `e7_fit.json`, `e7_phys.json`, `se7_A7.json`).
