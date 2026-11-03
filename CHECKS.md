# Checks: map, layout, and notes

Complete reference for the numerical checks in this repository.
All equation, section and item numbers refer to the **final paper** (LAPTH-050/26,
"Horizon exchange and the universal all-spin response of Kerr"):
main-text equations (1)–(23), Supplemental equations (S1)–(S48), sections S1–S8.

## Map of the checks

### Main text

| Paper | What is checked | Script(s) | Expected output |
|---|---|---|---|
| Eq. (6), Supplemental S1 | Monodromy algebra for both orderings of the loop product; identity with sin π(ν−iτ)/sin π(ν+iτ) | [`proof_check.py`](proof_check.py) | symbolic, residual 0 |
| Eq. (6) | Exchange law vs full MST (image computed directly), s = 0, −1, −2, ℓ ≤ 3, q ≤ 0.9, ε ≤ 0.6, 60 digits | [`exact_check.py`](exact_check.py) | relative deviation ≤ 2×10⁻⁶⁰ |
| Eq. (6), super-extremal | Same for s = 0 at q = 1.5, 2, 3 | [`test3.py`](test3.py) | ≤ 10⁻⁴⁹ (50 digits) |
| Fig. 1 | Loops of the figure by direct integration of the radial equation (s = −2, ℓ = m = 2, q = 0.6, ε = 0.3) | [`figure_monodromy/monodromy_loops_check.py`](figure_monodromy/) | M∞ = M₊M₋⁽⁺⁾ = M₋⁽⁻⁾M₊; eigenvalues exp(±2πε±), exp(±2πiν); the lasso above (below) r₊ gives the MST image with η = +1 (−1) to 8 digits |
| Eq. (13) | Exchange-odd ε⁶ coefficient, all m: closed form vs physical-region MST data; q⁶ coefficients equal ρ_\|m\|/64 of arXiv:2610.40101 | [`continuation_test/`](continuation_test/), [`prl_checks/paper_numbers.py`](prl_checks/) | Taylor coefficients to ~10⁻¹⁸; q⁶ coefficients exact |
| Eq. (14) | Interference term Δ_odd 𝒫_H for both η; one half of Page's 4ε⁶/225 at q = 0; fraction κ/(1+κ) for m = 0 | [`flux_test/decomp_check.py`](flux_test/) | 2×10⁻⁶; 5×10⁻⁷ |
| Eq. (15) | Exchange-even ε⁶ function, all m: free fit of the coefficients of Π_m ln κ and Π_m Re ψ(1+iϖ) | [`prl_checks/allm_analyze.py`](prl_checks/) | both −1/450 to 3×10⁻⁹ |
| Eq. (19) | Taylor coefficient [q⁵]𝒞 = 16b₂ of arXiv:2610.40101; 𝒞/q⁵ → 64/15 on the super-extremal branch | [`conservative_continuation/e6_identities.py`](conservative_continuation/) | 25 digits; 64/15 to 10⁻¹² at q = 10⁷ |
| Eq. (22), a⁵–a⁷ | Exact O(c²) ingredients of the law-based prediction (R₂ₘ, S₂₂, the ℓ = 3 term) | [`a7_exact.py`](a7_exact.py) | the bracket coefficients of Eq. (22), symbolically |
| Eq. (22), a⁵–a⁷ | Closed form (23) expanded and compared with the prediction | [`prl_checks/a7_polar_check.py`](prl_checks/) | difference 0 at c⁵, c⁶, c⁷ |
| Eq. (22), a⁸; Eq. (S14) | Expansion of the closed form: b₀–b₉, elementary form (S13), factor 3 + cos Θ of the odd orders | [`conservative_test/higher_orders_eq17.py`](conservative_test/) | symbolic |
| Eq. (23), Eq. (S12) | Polar classical amplitude from full MST (ℓ = 2…19, both η, 110 digits, ε = 10⁻¹⁶) against the closed form, order by order | [`conservative_test/`](conservative_test/) | 3×10⁻¹⁶ through O(a¹⁰), 10⁻¹⁵ through O(a¹³); η = +1 is the pure exponential exp(−2cv) |

### Supplemental Material

| Paper | What is checked | Script(s) | Expected output |
|---|---|---|---|
| S2, Eq. (S3) | Classical-limit structure of the response: X_odd/X_even = δ/u (s = 0) | [`classical.py`](classical.py) | ratio → δ/u as ε → 0 (8 digits at ε = 0.005) |
| S2, Eqs. (S4), (S7); S7(i) | Leading exchange-odd element H^odd_ℓm (s = −2 and s = 0) and the slopes 10m/9 of R₂ₘ (s = −2; the slope vanishes for s = 0) | [`polar.py`](polar.py), [`polar_m.py`](polar_m.py), [`slopes.py`](slopes.py) | leading coefficients to 10⁻⁵–10⁻⁴ (5×10⁻⁶ for ℓ = m = 2); slopes to 0.5% |
| S7(i) | H^odd vs full MST: ℓ = 2 remainder of Eq. (S7) ∝ c³; ℓ = 3 leading coefficient (m = 3, 2, 0, −1) | [`verify_a7.py`](verify_a7.py) | (ratio − series)/c³ constant; ℓ = 3 to 2×10⁻⁴ |
| S3, Eqs. (S16)–(S18) | Segment representation of the exchange-odd classical partial waves vs MST | [`segment_test/`](segment_test/) | ℓ = 2…6 to 4×10⁻¹⁶ |
| S3 | Size of the terms nonlinear in 𝒦 in the fixed-aω limit | [`segment_test/nonlinear_check.py`](segment_test/) | relative correction to the ℓ = 2 element −0.85c³ at c = 0.025 (≈ −0.82c³ as c → 0) |
| S4, logarithms; S7(c) | B₆ = iΠ_m/450 with naive and joint log bookkeeping | [`prl_checks/b6_precision.py`](prl_checks/) | naive 5×10⁻⁵ to 3×10⁻⁴; joint 10⁻¹⁰ |
| S4, uniqueness (a) | Free fit of the transcendental coefficients | [`prl_checks/allm_analyze.py`](prl_checks/) | −1/450 to 3×10⁻⁹ |
| S4, uniqueness (b) | Taylor coefficients 6 ≤ p ≤ 20 after subtracting Eq. (S24) | [`prl_checks/unified_check.py`](prl_checks/) | below 2×10⁻¹¹ |
| S4, uniqueness (c) | Polynomial fixed at \|q\| ≤ 0.3, continued to 1.5 ≤ q ≤ 5000, vs direct MST; bounds on extra terms c_p q^p | [`prl_checks/uniqueness_check.py`](prl_checks/) | 4×10⁻¹² (q ≤ 30), 3×10⁻¹⁰ to 6×10⁻¹³ (50 ≤ q ≤ 5000); \|c₆\| ≲ 3×10⁻¹⁸ |
| S4, uniqueness (d) | q⁵ coefficient of p₅⁽²⁾ vs the η = +1 classical c⁵ coefficient / 32 | [`prl_checks/uniqueness_check.py`](prl_checks/) and [`prl_checks/se_allm_analyze.py`](prl_checks/) | −2.3457997302×10⁻⁵ i; difference 1.5×10⁻¹⁵ |
| S5, Table S1 | General kinematics: BHPT exchange-odd amplitude vs the contact-term ansatz; fit of the nine exchange-odd coefficients; tagged and untagged exchange-even terms; a⁶ | [`kinematics.py`](kinematics.py), [`fitcoeffs.py`](fitcoeffs.py), [`prl_checks/kin_even.py`, `kin_tags.py`, `kin_a6.py`](prl_checks/) | ratios 1 to 10⁻³⁰ for helicity +1 (helicity −1 fails, as it should); fit residual 10⁻²³ |
| S6 (c)–(e); S7(j) | δ²𝒦₊ normalization and its 1 + (5mq/9)ε; E coefficient π − iν₂; ν₂ = −107/210 for all (m, q); ν₃ = mq/12 | [`prl_checks/analytic_odd.py`, `series_ratio_check.py`](prl_checks/) | X-ratio coefficient to ~10⁻⁶; E coefficient to 10⁻⁵; ν₃ to 0.5% |
| S6 (f); S7(j) | Pole term (i/900)ε⁴Π_m/ν₂ absent from S; spheroidal coefficients E₁, E₂; scaling of the near-zone factor | [`prl_checks/eps4_check.py`, `h_scaling.py`](prl_checks/) | E₁ exact, E₂ to 10⁻⁹; far zone exchange-invariant to 10⁻⁷⁴ |
| S7(b) | Exchange invariance of ν and of the far-zone factor | [`test1.py`](test1.py), [`test3.py`](test3.py) | ν to 3×10⁻⁵⁷, far zone to 3×10⁻⁵¹ |
| S7(d) | Eq. (13) for all m; q⁶ coefficients = ρ_\|m\|/64 | [`continuation_test/`](continuation_test/) | exact |
| S7(e) | [q⁵]𝒞 = 16b₂; 𝒞/q⁵ → 64/15 | [`conservative_continuation/`](conservative_continuation/) | 25 digits |
| S7(f) | Classical super-extremal a⁵ coefficients per m: continued closed forms vs direct MST | [`prl_checks/se_allm_analyze.py`](prl_checks/), [`continuation_test/d6_general.py`](continuation_test/) | 10⁻⁹ |
| S7(g), Eq. (S35) | ε⁷ operator for all m, both sectors | [`prl_checks/d7_phys.py`, `d7_fit2.py`, `e7_phys.py`, `e7_allm.py`](prl_checks/) | π and i(2γ_E + 2 ln 2 − 37/12 − ν₂) to 3×10⁻⁸; coefficient of q equal to 5m/9 to 2×10⁻⁸; no q² term |
| S7(h) | Classical a⁶ coefficients (m = 2) by extrapolating 64A₇/q⁶ to q → ∞ | [`prl_checks/se7_analyze.py`](prl_checks/) | 1.8964816 (exchange-even) and 1.8962963 (exchange-odd), vs 1.89648161 and 256/135 |
| S7(j) | Γ₊ − Γ₋ = −4ε⁶A₆⁻ (physical solution vs its image) | [`prl_checks/absorption_check.py`](prl_checks/) | ≤ 2.5×10⁻⁶ after linear extrapolation in ε |
| S7(k), S8 | Unit-flux normalization of the far-zone factor; flux balance; complex-ν branch; exchange-parity split of 𝒫_H | [`flux_test/`](flux_test/) | 10⁻⁵⁹ to 10⁻⁶¹; split to 2×10⁻⁶ |

Item S7(a) is covered by `exact_check.py`, `test3.py` and `figure_monodromy/`; the physical-region data
behind Eqs. (13), (15) and (19) and the comparison with arXiv:2610.40101 are produced in
[`physical_region_test/`](physical_region_test/).

## Repository layout

```
.
├── run_all_checks.py              runs all checks, writes check_report/report.pdf
├── mst.py                         core MST solver (see README)
├── proof_check.py                 symbolic proof of the exchange law (S1)
├── exact_check.py, test1.py, test3.py   exchange law and invariances vs full MST (S7 a, b)
├── classical.py                   classical-limit structure of the response (S2)
├── polar.py, polar_m.py, slopes.py      spheroidal eigenvalues; leading exchange-odd elements and R_2m slopes (S2, S7 i)
├── a7_exact.py, verify_a7.py      O(a^7) prediction of Eq. (22) and its ingredients
├── kinematics.py, fitcoeffs.py    general kinematics, exchange-odd rows of Table S1 (S5)
├── conservative_test/             polar classical amplitude, closed form (23)/(S12), higher orders (S3)
├── segment_test/                  segment representation (S16)–(S18) (S3)
├── physical_region_test/          physical-region double expansion in (eps, q), production data (S4)
├── continuation_test/             exchange-odd eps^6 function, Eq. (13), and its continuation (S7 d, f)
├── conservative_continuation/     exchange-even eps^6 function, C(q), Eqs. (15), (19) (S4, S7 e)
├── prl_checks/                    all-m functions, eps^7 operator, general kinematics, analytic control (S4–S6, S7)
├── flux_test/                     unit-flux normalization and flux balance (S8, S7 k)
└── figure_monodromy/              Fig. 1 (TikZ) and the direct check of its loop algebra
```

Each folder has its own README with script-by-script details.

## Running all checks

[`run_all_checks.py`](run_all_checks.py) runs every check script from its own folder, compares its output
with the numbers quoted in the paper and in these READMEs, and writes the results to `check_report/`:

| File | Content |
|---|---|
| `report.pdf` | summary table; for each check its purpose, the criteria (measured value, required value and tolerance) and the complete output of the script; Fig. 1 next to the check of its loops |
| `report.tex` | LaTeX source of the report, written even when `pdflatex` is not available |
| `summary.txt` | one line per check: status, running time, key result |
| `results.json` | the same in machine-readable form, with every criterion |
| `logs/` | the complete output of every script |

Each check ends with one of the statuses **PASS** (all criteria met), **FAIL** (a criterion not met),
**ERROR** (the script did not finish), **INFO** (exploratory script without criteria) or **SKIP** (a data
file or module is missing). The negative controls of the repository pass when the comparison fails by a
wide margin. The exit status is nonzero when a check fails or ends with an error.

```bash
python3 run_all_checks.py --list          # the 59 checks, with paper items and running times
python3 run_all_checks.py -j 4            # up to four scripts at a time (default: number of CPUs, at most 8)
python3 run_all_checks.py --quick         # skip slopes.py, verify_a7.py, d6_m.py and vs_production.py
python3 run_all_checks.py --only segment_test closedform_check   # checks whose path contains a pattern
python3 run_all_checks.py --production    # regenerate the data files first, then run all checks
```

With `--production`, the production drivers run first, one after the other (more than an hour on two
cores). Every shipped data file they overwrite is first copied to `check_report/backup/`. Without
`--production`, `physical_region_test/analyze_phys.py` is skipped because its input `phys_run.json` is
not shipped.

## Production runs and data files

The expensive MST evaluations are done once by *production drivers* (they use `multiprocessing` with
two workers and write JSON files). The analysis scripts read those files, so they run in seconds. All
data needed by the analysis scripts are included, except `physical_region_test/phys_run.json` (~5 MB),
which `physical_region_test/run_phys.py` regenerates.

| Driver | Output | Time on 2 cores |
|---|---|---|
| `conservative_test/run_polar2.py` | `polar_run2.json` (ℓ = 2…19, 28 nodes, both η, 110 digits) | ~15 min |
| `physical_region_test/run_phys.py` | `phys_run.json` (not included) | ~13 min |
| `physical_region_test/run_phys2.py 2e-5`, `run_phys2.py 1e-5` | `phys_l2_<radius>.json` (copies in `prl_checks/`) | ~4 min per radius |
| `prl_checks/run_allm.py` | `phys_allm.json`, `se_allm.json` | ~20 min |
| `conservative_continuation/e6_se.py`, `e6_se2.py`, `e6_se3.py` | `e6_se*.json` (copies in `prl_checks/`) | not timed |
| `continuation_test/d6_se.py` | `d6_se.json` | not timed |

Running a driver overwrites the corresponding shipped file. `python3 run_all_checks.py --production`
runs them all in order and keeps a backup. A few analysis scripts read data from a sibling folder via
relative paths (`../prl_checks/phys_l2_*.json`, `../conservative_test/polar_run2.json`,
`../conservative_test/taylor2.json`); keep the folder structure intact.

## Running times and precision

Measured on a two-core machine (Python 3.13, mpmath 1.3.0, no gmpy2):

- Most analysis scripts finish in under a minute, many in under a second.
- The longer ones: `slopes.py` ~1 min, `verify_a7.py` ~1.2 min, `continuation_test/d6_m.py` ~1.7 min
  (two worker processes), `segment_test/vs_production.py` ~2 min, `flux_test/decomp_check.py` ~45 s,
  `segment_test/segment_born.py` ~35 s, `polar_m.py` ~30 s.
- `run_all_checks.py` runs all 59 checks in about 8 minutes with two scripts at a time (`--quick`: ~4 min).

Each script sets its own working precision (`mp.mp.dps` = 30 to 110 digits). Installing
[gmpy2](https://pypi.org/project/gmpy2/) speeds up mpmath considerably. Small differences in the last
printed digits between machines are expected; the quoted accuracies are robust.

## Notes on the scripts

- **Equation numbers in the scripts.** Docstrings and comments of most scripts use the numbering of an
  earlier draft (v4); those in `segment_test/` and `conservative_test/higher_orders_eq17.py` use a v1
  arXiv numbering (where the polar series was Eq. (17) and the closed form Eq. (18); these are Eqs. (22)
  and (23) in the final paper). The correspondence with the final paper (LAPTH-050/26) is:

  | Draft v4 | Final paper | Content |
  |---|---|---|
  | (6), (7), (9) | (6), (7), (9) | exchange law, monodromy form, exchange-odd S-matrix |
  | (10) | (S34) | Boltzmann form of the exchange-odd amplitude |
  | (11) | (13) | Bose-factor form of the exchange-odd ε⁶ coefficient |
  | (12) | (S24) | ε⁶ function of the physical solution |
  | (13) | (15) | exchange-even ε⁶ function |
  | (14) | (19) | 𝒞(q) |
  | (15) | (S35) | ε⁷ operator |
  | (16), (17) | (S4), (S7) | H^odd_ℓm, R₂ₘ |
  | (18) | (22) | polar series in aω |
  | Supplement (S5) | (S12) | closed form of the polar amplitude |

- **Terminology in the scripts.** Earlier drafts called the horizon exchange an *irreducible-mass duality*.
  In the code, "self-dual" or "duality-even" means exchange-even, "anti-self-dual" or "duality-odd"
  means exchange-odd, "dual" or "inner" solution means the image, and the "swap law" is the exchange law.
  `D6`, `D7` denote exchange-odd ε⁶ and ε⁷ coefficients, `E6`, `E7` exchange-even ones. "BGKV" refers
  to Bautista, Guevara, Kavanagh and Vines, arXiv:2107.10179 and arXiv:2212.07965, and "BBCPSS" to
  Bjerrum-Bohr, Chen, Papadimos, Shah and Su, arXiv:2610.40101.

- **Negative controls.** Some scripts include deliberate negative tests: the opposite helicity in
  `kinematics.py` and `prl_checks/kin_even.py`, wrong loop orderings in
  `figure_monodromy/monodromy_loops_check.py`, and the shifted-argument hypothesis in
  `prl_checks/a6_hyp.py`. Their failure is the expected result.
