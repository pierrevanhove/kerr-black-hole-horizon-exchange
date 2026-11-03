# flux_test: unit-flux normalization and flux balance

Unit-flux normalization of the MST far-zone factor F_ℓm. s = −2, ℓ = 2, M = 1, ε = 2Mω, q = a/M;
60 digits, same MST code and spheroidal solver as all other checks.

Paper: Supplemental S10 and S9(k); main text Eq. (11).

## Scripts

| Script | Purpose | Time |
|---|---|---|
| `flux_check.py` | B_inc/ref[R_C^ν] and B_inc/ref[R_C^(−ν−1)] evaluated separately from the continued-fraction coefficients a_n; F_P = (−1)^(ℓ+1) 𝒞_P/(16ω⁴) B^ref[R_C^ν]/B^inc[R_C^ν] with 𝒞_P = D + 12iPMω; all m, q = 0…0.95, ε = 0.01…0.4, including superradiant cases | ~15 s |
| `highfreq_check.py` | The same flux balance at ε = 0.6 (Mω = 0.3), ν followed by continuation in ε | ~10 s |
| `complexnu_check.py` | The complex-ν branch (ν = 1/2 + iy, 3/2 + iy, −1/2 + iy), ε = 0.8…1.2 | ~20 s |
| `decomp_check.py` | Exchange-parity split of the unit-flux absorption probability, 1 − \|S₊\|² = [1 − \|S_even\|²] + Δ_odd − \|S_odd\|², with the image S₋ computed directly from MST (κ → −κ, path η; the exchange law is not used); ε⁶ coefficients by Richardson extrapolation from ε = 4×10⁻³, 2×10⁻³, 10⁻³ | ~45 s |
| `mst.py`, `phys.py`, `csm.py`, `polar.py` | Core modules (copies) | — |

The files `out_*.txt` contain the reference output of each script.

## Expected output

`flux_check.py`:

- MST relations B[R_C^(−ν−1)]/B[R_C^ν] = −i exp(−iπν)ρ (incoming) and i exp(iπν) (outgoing): 10⁻⁵⁹.
- Factorized F_P × (response fraction) equals the direct 𝒞_P/(16ω⁴) B^ref/B^inc of R_in: 10⁻⁶⁰.
- \|𝒞\|² as written in S10 equals D² + 144M²ω²: exact.
- Flux balance \|S_P\|² + 𝒫_H = 1, with 𝒫_H = (2Mr₊)⁵ k(k² + ϰ²)(k² + 4ϰ²)\|B^trans/B^inc\|²/ω⁵: 10⁻⁶⁰.
- Helicity-preserving, D-normalized S: \|S\|² + 𝒫_H − 1 = −\|S_flip\|² ≈ −ε²/16, and \|S/S_P\| = D/\|𝒞\| to 10⁻⁶¹.

`complexnu_check.py`: flux balance to 2×10⁻⁵⁹ on the complex-ν branch, superradiant amplification up to
10.7% (m = 2, q = 0.95, Mω = 0.6).

`decomp_check.py`:

- Δ_odd → (π/225) κ⁵ ϖ(1 + ϖ²)(4 + ϖ²) n̄_η, with n̄₊₁ = 1/(exp(2πϖ) − 1) and n̄₋₁ = n̄₊₁ + 1, to
  2×10⁻⁶ for (m, q) = (2, 0), (0, 0.6), (0, 0.3), (2, 0.3), (1, 0.5), (−2, 0.4), both η.
- q = 0: Δ_odd = 2ε⁶/225 = 𝒫_H/2 (Page: 𝒫_H = 4ε⁶/225).
- m = 0: Δ_odd/𝒫_H = κ/(1 + κ), with 𝒫_H = (2/225) κ⁴(1 + κ) ε⁶, to 5×10⁻⁷.
