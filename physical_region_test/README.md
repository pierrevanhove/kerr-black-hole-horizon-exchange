# physical_region_test: double expansion in (ε, q) in the physical region

Physical region |a| < M, s = −2, m = 2 (polar projection). The Teukolsky S-matrix is expanded in
ε = 2GMω and q = a/GM by a two-dimensional Cauchy integral: complex ε on the circles |ε| = 2×10⁻⁵ and
10⁻⁵, complex q on |q| = 0.3–0.35, full MST at 90 digits, for the physical branch κ = +√(1 − q²) (outer
horizon) and the image κ = −√(1 − q²) (inner horizon). The coefficient of ε^k q^p corresponds to
ε̄^(k−p) |a|^p in the counting of Bjerrum-Bohr, Chen, Papadimos, Shah and Su, arXiv:2610.40101.

Paper: data behind Eqs. (10), (12) and (16); Supplemental S5.

## Scripts

| Script | Purpose | Time |
|---|---|---|
| `phys.py` | S-matrix in (ε, q) for either branch (module) | — |
| `run_phys.py` | Production: ℓ = 2…7, L = 2 polar projection, one radius → `phys_run.json` (not included, about 5 MB) | ~13 min |
| `run_phys2.py` | Production: ℓ = 2, both branches, radius given as argument (`python3 run_phys2.py 2e-5`, `python3 run_phys2.py 1e-5`) → `phys_l2_<radius>.json` (the two files are included in `../prl_checks/`) | ~4 min per radius |
| `analyze_phys.py` | Coefficients of the L = 2 projection, comparison with Eqs. (C.1), (4.6)–(4.7) of arXiv:2610.40101 and with the super-extremal data `../conservative_test/taylor2.json`. Needs `phys_run.json`: run `run_phys.py` first | seconds |
| `rich_phys.py` | Richardson-improved ℓ = 2 coefficients for both branches, from `../prl_checks/phys_l2_*.json` → `l2_rich.json` (included) | ~10 s |
| `pslq2.py` | Exploratory PSLQ identification of individual coefficients in a basis of zeta values; superseded by the closed forms (digamma function) of the paper | seconds |
| `mst.py`, `polar.py`, `csm.py` | Core modules (copies) | — |

`phys_run_coeffs.json` and `l2_rich.json` contain the extracted coefficients. `l2_rich.json` is also
used by `../continuation_test/d6_closed.py`.

## Expected output

`rich_phys.py` prints the Richardson changes, the q = 0 value of 1 − |S|² at order ε⁶, and the
exchange-odd and exchange-even coefficients for p ≤ 10. For the physical branch the q = 0 value is
0.0180219…, which is Page's 4/225 plus 1/4096 = (1/16)³: this S is the helicity-preserving
element, and at q = 0, |S|² = |S_P|²/(1 + ε²/16) (see `../flux_test/`). The exchange-odd q⁶
coefficient 0.00621387874520449… equals ρ₂/64 of arXiv:2610.40101.
