# Kerr horizon-exchange law: numerical checks
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23922017.svg)](https://doi.org/10.5281/zenodo.23922017)
[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/XXXX.XXXXX)
[![Author: Pierre Vanhove](https://img.shields.io/badge/author-Pierre_Vanhove-blue)](https://pierrevanhove.github.io)
[![Language: Python](https://img.shields.io/badge/Language-Python-yellow?logo=python)](https://www.python.org/downloads/)
[![mpmath](https://img.shields.io/badge/mpmath-1.3%2B-orange.svg)](https://mpmath.org)

Scripts and data for the numerical checks in

> **Horizon exchange and the universal all-spin response of Kerr**\
> Pierre Vanhove, LAPTH-050/26, 

The radial Teukolsky equation for Kerr is invariant under exchange of the outer and inner horizons,
although the physical boundary condition is not. The paper derives an exact monodromy relation between
the tidal responses selected at the two horizons (the *exchange law*, Eq. (6)) and uses it to split the
Kerr–Compton response into exchange-odd and exchange-even parts. This repository reproduces the
numerical statements of the paper and of its Supplemental Material.

Equation, section and item numbers refer to the final paper (LAPTH-050/26): main-text equations are
(1)–(23), Supplemental equations (S1)–(S48), and Supplemental sections S1–S8. See [CHECKS.md](CHECKS.md)
for the complete map of checks with paper references.

## Requirements

- Python 3 (tested with Python 3.13; `run_all_checks.py` needs Python 3.8 or later)
- [mpmath](https://mpmath.org) for all numerics (tested with 1.3.0)
- [SymPy](https://www.sympy.org) for the symbolic checks (tested with 1.14)
- [NumPy](https://numpy.org), only for `figure_monodromy/monodromy_loops_check.py`
- `pdflatex`, only for the PDF report of `run_all_checks.py` and (with TikZ) to rebuild Fig. 1

```bash
pip install -r requirements.txt
```

## Quick start

```bash
python3 run_all_checks.py         # ~8 min on two cores; writes check_report/report.pdf
python3 proof_check.py            # symbolic proof of the exchange law (S1), seconds
python3 exact_check.py            # exchange law vs full MST at 60 digits, ~10 s

cd conservative_test
python3 closedform_check.py       # closed form (23) vs full MST, order by order, seconds
python3 higher_orders_eq17.py     # a^8 term of Eq. (22) and coefficients of Eq. (S14), ~30 s

cd ../segment_test
python3 vs_production.py          # segment representation (S16)–(S18) vs MST for l=2...6, ~2 min
```

Each folder carries copies of the core modules it imports; no `PYTHONPATH` setup is needed.

## Conventions

| Quantity | Meaning in the scripts |
|---|---|
| `M = 1` | units; `eps` $= \varepsilon = 2M\omega$, `q` $= a/M$ |
| `kappa` | $\kappa = \sqrt{1-q^2}$; horizons $r_\pm = M(1\pm\kappa)$ |
| `tau` | $\tau = (\varepsilon - mq)/\kappa$ |
| `varpi`, `w` | $\varpi = mq/\kappa = m\Omega_H/(2\pi T_H)$ |
| `c` | $c = a\omega$; the *fixed-*$a\omega$*, small-*$M\omega$ *limit* is $\varepsilon\to 0$ at fixed $c$ |
| `X` | the response $\mathcal{K} = K_{-\nu-1}/K_\nu$ of Eq. (5), `MST.amplitudes()['X']` |
| `eta` | continuation label $\eta = \pm 1$ |

- **Exchange.** The image solution is the MST solution with $\kappa\to-\kappa$ and $\log(-\kappa) = \log\kappa + i\pi\eta$, computed directly without using the exchange law.
- **Super-extremal branches.** For $q>1$ the two continuations are $\kappa = \eta\,i\,q\sqrt{1-1/q^2}$, so that in the fixed-$a\omega$ limit the ingoing condition sits at $r = \eta\,i\,a$.
- **S-matrix normalization.** The helicity-preserving element used for the Compton amplitude is $S = (-1)^{\ell+1} D/(16\omega^4)\cdot B_\text{ref}/B_\text{inc}$, with $D$ the real Teukolsky–Starobinsky constant (`csm.D_TS`). The unit-flux element $S_P$ carries $\mathcal{C}_P = D + 12iPM\omega$ instead (Supplemental S8, [`flux_test/`](flux_test/)).
- **Exchange-odd and exchange-even parts** of any quantity $F$ are $(F^{(+)} - F^{(-)})/2$ and $(F^{(+)} + F^{(-)})/2$.

## Core modules

The same files are copied into every folder that needs them (identical copies).

- **`mst.py`**: minimal MST solver for the radial Teukolsky equation, following Sasaki and Tagoshi, Living Rev. Rel. 6 (2003) 6, Sec. 4. `MST(eps, q, m, s, lam, kappa, logkappa=None, N=50)` builds the three-term recurrence; `solve_nu` finds the renormalized angular momentum $\nu$ from the continued fractions; `amplitudes()` returns $\nu$, the coefficients $a_n$, $K_\nu$, $K_{-\nu-1}$, the response `X` $= \mathcal{K}$, near-zone factor `H`, asymptotic amplitudes $B_\text{inc}$/$B_\text{ref}$/$B_\text{tr}$, and $S_\text{far}$/$S = B_\text{ref}/B_\text{inc}$.
- **`polar.py`**: matrix of $\cos\theta$ in the spin-weighted spherical basis (`xmat`) and the spin-weighted spheroidal eigenvalue (`E_spheroidal`); run directly, it checks the leading exchange-odd $\ell=m=2$ element.
- **`csm.py`**: classical-limit S-matrix: spheroidal eigenvalue/eigenvector (`spheroidal`), Teukolsky–Starobinsky constant (`D_TS`), robust $\nu$ solvers, and the physical helicity-preserving `S_phys(l, m, s, c, eta, eps, x0)`.
- **`phys.py`**: physical-region S-matrix in $(\varepsilon, q)$ for either branch (`S_tilde`), used by the two-dimensional Cauchy-integral expansions.

For the full check tables, repository layout, production run details, running times, and the draft-to-final equation-number mapping, see **[CHECKS.md](CHECKS.md)**.

## Citation

If you use these scripts, please cite the paper:

```bibtex
@article{Vanhove:horizonexchange,
  author        = {Vanhove, Pierre},
  title         = {Horizon exchange and the universal all-spin response of {K}err},
  eprint        = {XXXX.XXXXX},
  archivePrefix = {arXiv},
  year          = {2026}
}
```

Archived release: Zenodo DOI to be added.


## Acknowledgments

This work has been supported by the Agence Nationale de la Recherche (ANR) under grant Observables
(ANR-24-CE31-7996). The Python code was written and run with the assistance of Claude, an AI model
developed by Anthropic; the author takes full responsibility for its content.
