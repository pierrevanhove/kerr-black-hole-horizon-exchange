# segment_test: segment representation of the exchange-odd partial waves

Paper: Supplemental S4, Eqs. (S15)–(S17); main text after Eq. (18).

s = −2, m = 2 (incidence along the spin), classical limit ε → 0 at fixed c = aω, units a = 1.

The classical amplitude is first order in M. At M = 0 the radial Teukolsky equation is the spheroidal
equation, the horizons sit at r = ±ia, and the solution regular at both is R₀(r) = Δ₀(r) Ŝ_ℓ(ir/a), with
Δ₀ = r² + a² and Ŝ_ℓ the spin-weighted spheroidal function ₋₂S_ℓ2(x; c) continued to complex x.
The continuation η imposes the ingoing condition at r = η i a. The Green identity between the two
horizons gives, for the exchange-odd part of the partial-wave element,

```math
e^{2i\delta_{\ell 2}}\Big|_{\rm odd} = \frac{\epsilon\, G_\ell}{8 i \omega^2 B_{\rm inc} B_{\rm ref}},\qquad
G_\ell = {\rm FP}\!\int_{-ia}^{ia}\frac{R_0\, T_1 R_0}{\Delta_0^2}\,dr + \hat S_\ell(1)\big[(2c-1)\hat S_\ell(1) - 4\hat S_\ell'(1)\big],
```

where T₁ is the derivative of the radial operator with respect to M at M = 0, Eq. (S15), and FP is the
finite part at the double pole r = −ia.

## Scripts

| Script | Purpose | Time |
|---|---|---|
| `segment_born.py` | Eq. (S16): finite-part integral along r = ia′, endpoint term, asymptotic amplitudes B_inc, B_ref by Taylor integration of the M = 0 radial equation. Run directly: ℓ = 2 at c = 0.05, 0.1, 0.2 against full MST at ε = 10⁻¹² | ~35 s |
| `fp_integral.py` | Finite part: Cauchy coefficients near the double pole, quadrature elsewhere (module) | — |
| `vs_production.py` | ℓ = 2…6 at c = 0.1997 and 0.1631 against the 110-digit MST data `../conservative_test/polar_run2.json` | ~2 min |
| `nonlinear_check.py` | Size of the terms nonlinear in 𝒦 in the MST form S/F − 1 = iE𝒦/(1 + β𝒦) | ~10 s |
| `mst.py`, `polar.py`, `csm.py` | Core modules (copies) | — |

## Expected output

`segment_born.py` prints the internal checks (residual of the radial equation 10⁻⁴⁰, cancellation of
the double pole against the boundary Wronskian, vanishing residue, S = 1 at M = 0, the identity for
the endpoint term) and

```
   S_odd/eps (segment formula) = (-6.07309597459634e-40 - 1.06837526731227e-5j)
   full MST (eps=1e-12):          (0.0 - 1.0683752673e-5j)    ratio-1 = (1.15e-11 - 5.68e-35j)
```

The agreement is limited by the 12 printed digits of the MST reference values.

`vs_production.py` gives relative differences between 3×10⁻¹⁶ and 4.4×10⁻¹⁶ for ℓ = 2…6, the size of
the finite-ε systematic of the MST data:

```
l=2 c=0.19968536: MST (... - 0.0004264947435263j)   segment (... - 0.0004264947435263j)   |rel diff| 3.15e-16
l=6 c=0.16311217: MST (... - 2.014942080859e-22j)   segment (... - 2.014942080859e-22j)   |rel diff| 4.43e-16
```

`nonlinear_check.py` shows that β𝒦 ≈ 0.43c³ for ℓ = 2, so the terms nonlinear in 𝒦 change the
exchange-odd ℓ = 2 element at relative order c³ (−0.85c³ at c = 0.025, tending to about −0.82c³). The
segment representation contains them automatically.
