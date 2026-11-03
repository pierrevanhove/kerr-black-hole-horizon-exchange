# figure_monodromy: Fig. 1 and the direct check of its loop algebra

Paper: Fig. 1; Supplemental S1 and S9(a).

| File | Content |
|---|---|
| `fig_monodromy.tex`, `fig_monodromy.pdf` | TikZ source (compile with `pdflatex fig_monodromy.tex`) and the figure included in the paper |
| `monodromy_loops_check.py` | Direct integration of the radial equation along the loops drawn in Fig. 1 (needs NumPy) |
| `mst.py`, `phys.py`, `csm.py`, `polar.py` | Core modules (copies) |

## What the check does

RK4 integration of the s = −2 radial Teukolsky equation (ℓ = m = 2, q = 0.6, ε = 0.3) along the loops of
Fig. 1, based at r₀ = 3 > r₊, all counterclockwise. Transport matrices act on (R, R′) at r₀, so "γ then
γ′" has monodromy T_γ′ T_γ. Runs in about 15 s.

## Expected output

1. T_∞ = T₊ T₋(above) = T₋(below) T₊ to ~10⁻⁴ (entries ~10⁴; the wrong orderings miss by ~7×10⁴).
2. Eigenvalues: T₊ → exp(±2πε₊), T₋ → exp(±2πε₋), T_∞ → exp(±2πiν) with ν = 1.9554220852 from the MST
   continued fraction (8 digits).
3. 𝒦₋/𝒦₊ in the Floquet basis (eigenvectors of T_∞): the lasso passing above r₊ reproduces the MST image
   with log(−κ) = log κ + iπ (η = +1), the lasso passing below gives η = −1 (8 digits):

```
(3) K_-/K_+ from monodromy: above r+ -> (0.99999989+0.00047113j)   below r+ -> (0.84690767+0.53173997j)
    MST image, eta=+1        -> (0.99999989+0.00047113j)   eta=-1 -> (0.84690768+0.53173996j)
```
