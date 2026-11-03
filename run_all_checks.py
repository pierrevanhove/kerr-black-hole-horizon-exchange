#!/usr/bin/env python3
"""Run all numerical checks of the repository and collect the results in a PDF report.

Every check script is run from its own folder. Its output is compared with the numbers stated in the
paper and in the READMEs (each comparison is a *criterion* with an explicit tolerance), and everything is
collected in the output folder (default check_report/):

    report.pdf     the report: summary table, every criterion, the complete output of every script
    report.tex     its LaTeX source (written even when pdflatex is not available)
    summary.txt    one line per check
    results.json   machine-readable results, including every criterion
    logs/          the complete output of every script

Usage:
    python3 run_all_checks.py                  all checks (about 8 min on two cores)
    python3 run_all_checks.py --quick          skip the four longest scripts (about 4 min on two cores)
    python3 run_all_checks.py -j 4             run up to four scripts at a time
    python3 run_all_checks.py --only segment_test closedform_check
    python3 run_all_checks.py --list           list the checks and exit
    python3 run_all_checks.py --production     regenerate the data files with the production drivers
                                               first (about 2 h on two cores), then run all checks

The exit status is 0 when no check fails or ends with an error.
"""
import argparse
import concurrent.futures as cf
import datetime
import fnmatch
import json
import math
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ----------------------------------------------------------------------------------------------------
# Constants used as targets (all quoted in the paper)
# ----------------------------------------------------------------------------------------------------
EULER_GAMMA = 0.57721566490153286061
NU2 = -107 / 210                                   # nu = 2 + nu2 eps^2 + ..., S7(c)
IM_U0 = 2 * EULER_GAMMA + 2 * math.log(2) - 37 / 12 - NU2   # Im u0 of Eq. (S23)
RHO2_64 = 0.0062138787452045028789                 # rho_2/64 of arXiv:2610.40101
PAGE = 4 / 225                                     # Page: ell = 2 absorption at q = 0, order eps^6


def fmtg(x, d=3):
    """Compact %g formatting: 2.7e-9 rather than 2.7e-09."""
    if isinstance(x, complex):
        sign = '+' if x.imag >= 0 else '-'
        return f'({fmtg(x.real, d)} {sign} {fmtg(abs(x.imag), d)}i)'
    if x is None:
        return 'n/a'
    s = f'{x:.{d}g}'
    s = re.sub(r'e\+0*(\d)', r'e\1', s)
    return re.sub(r'e-0*(\d)', r'e-\1', s)


# ----------------------------------------------------------------------------------------------------
# Parsing of the script outputs
# ----------------------------------------------------------------------------------------------------
# Patterns are regular expressions in which {N} captures a real number and {C} a complex number printed
# by mpmath, NumPy or Python, for example (0.99999989+0.00047113j) or (1.0 - 3.1e-31j); {n} and {c}
# match the same without capturing. Patterns must not contain other capturing groups.
_NUM = r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?'
_NUMU = r'(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?'
_TOKENS = {
    '{N}': ('(' + _NUM + ')', 'N'),
    '{C}': (r'\(\s*(' + _NUM + r')\s*([-+])\s*(' + _NUMU + r')j\s*\)', 'C'),
    '{n}': ('(?:' + _NUM + ')', None),
    '{c}': (r'\(\s*' + _NUM + r'\s*[-+]\s*' + _NUMU + r'j\s*\)', None),
}
_COMPILED = {}


def _compile(pattern):
    if pattern not in _COMPILED:
        regex, kinds = '', []
        for part in re.split(r'(\{[NCnc]\})', pattern):
            if part in _TOKENS:
                rx, kind = _TOKENS[part]
                regex += rx
                if kind:
                    kinds.append(kind)
            else:
                regex += part
        _COMPILED[pattern] = (re.compile(regex, re.M), kinds)
    return _COMPILED[pattern]


def _values(match, kinds):
    vals, g = [], 1
    for kind in kinds:
        if kind == 'N':
            vals.append(float(match.group(g)))
            g += 1
        else:
            im = float(match.group(g + 2))
            vals.append(complex(float(match.group(g)), im if match.group(g + 1) == '+' else -im))
            g += 3
    return vals


def complexes(text):
    """All complex numbers printed in a piece of text, e.g. a printed list."""
    rx, kinds = _compile('{C}')
    return [_values(m, kinds)[0] for m in rx.finditer(text)]


@dataclass
class Criterion:
    desc: str
    measured: str
    required: str
    ok: bool


class Ev:
    """Evaluation of one script output: a list of criteria, each with a measured value and a tolerance."""

    def __init__(self, out):
        self.out = out
        self.crit = []
        self.headline = None

    # -- extraction
    def find(self, pattern, text=None):
        rx, kinds = _compile(pattern)
        return [_values(m, kinds) for m in rx.finditer(self.out if text is None else text)]

    def rows(self, pattern, desc, n=1, text=None):
        """All matches of `pattern`; records a failed criterion when fewer than n lines match."""
        r = self.find(pattern, text)
        if len(r) < n:
            self.add(desc, f'{len(r)} matching line(s) in the output', f'at least {n}', False)
            return None
        return r

    def first(self, pattern, desc, text=None):
        r = self.rows(pattern, desc, 1, text)
        return r[0] if r else None

    def section(self, start, desc):
        """The block of output that follows the line containing `start`, up to the next blank line."""
        i = self.out.find(start)
        if i < 0:
            self.add(desc, 'section not found in the output', start, False)
            return None
        rest = self.out[i:].split('\n', 1)
        rest = rest[1] if len(rest) > 1 else ''
        return re.split(r'\n\s*\n', rest, maxsplit=1)[0]

    # -- criteria
    def add(self, desc, measured, required, ok):
        self.crit.append(Criterion(desc, measured, required, bool(ok)))

    def has(self, desc, text):
        ok = text in self.out
        self.add(desc, 'printed' if ok else 'not found', f'"{text}"', ok)

    def le(self, desc, values, tol, d=3):
        """max |v| <= tol."""
        if not values:
            return
        mx = max(abs(v) for v in values)
        meas = fmtg(mx, d) if len(values) == 1 else f'max {fmtg(mx, d)} ({len(values)} values)'
        self.add(desc, meas, f'<= {fmtg(tol, 2)}', mx <= tol)
        return mx

    def ge(self, desc, values, tol, d=3):
        """min |v| >= tol (negative controls)."""
        if not values:
            return
        mn = min(abs(v) for v in values)
        meas = fmtg(mn, d) if len(values) == 1 else f'min {fmtg(mn, d)} ({len(values)} values)'
        self.add(desc, meas, f'>= {fmtg(tol, 2)}', mn >= tol)
        return mn

    def near(self, desc, pairs, tol, rel=False, d=3, show=True):
        """max |v - t| (or |v - t|/|t|) <= tol over the (value, target) pairs."""
        if not pairs:
            return
        diffs = [abs(v - t) / (abs(t) if rel and t != 0 else 1) for v, t in pairs]
        mx = max(diffs)
        word = 'rel. diff' if rel else 'diff'
        if len(pairs) == 1 and show:
            v = pairs[0][0]
            meas = f'{fmtg(v, 12) if not isinstance(v, complex) else fmtg(v, 10)} ({word} {fmtg(mx, 2)})'
        else:
            meas = f'max {word} {fmtg(mx, d)}' + (f' ({len(pairs)} values)' if len(pairs) > 1 else '')
        self.add(desc, meas, f'{word} <= {fmtg(tol, 2)}', mx <= tol)
        return mx

    def within(self, desc, values, lo, hi, d=3):
        if not values:
            return
        ok = all(lo <= v <= hi for v in values)
        meas = fmtg(values[0], d) if len(values) == 1 else \
            f'{fmtg(min(values), d)} ... {fmtg(max(values), d)} ({len(values)} values)'
        self.add(desc, meas, f'in [{fmtg(lo)}, {fmtg(hi)}]', ok)

    def count(self, desc, pattern, n):
        k = len(self.find(pattern))
        self.add(desc, f'{k} line(s)', f'{n} line(s)', k == n)


def _num_list(text):
    return [float(v.strip().strip("'\"")) for v in text.split(',') if v.strip()]


def _c5_closed(m):
    """Leading exchange-odd element, s = -2, l = 2, Eq. (S4): -(-1)^m (8/225) i (2+m)!(2-m)!."""
    return -(-1) ** m * 8 / 225 * 1j * math.factorial(2 + m) * math.factorial(2 - m)


# ----------------------------------------------------------------------------------------------------
# Evaluators: one function per script, called with an Ev on the script output
# ----------------------------------------------------------------------------------------------------
# -- Exchange law and monodromy
def ev_proof_check(e):
    e.has('Loop product M_- = M_inf M_+^(-1)', '= -lam^-1 (lam a - b)/(lam b - a)')
    e.has('Loop product M_- = M_+^(-1) M_inf', '= -lam^+1 (lam a - b)/(lam b - a)')
    r = e.first(r'sin pi\(nu \+ i tau\): {N} [-+] {N}\*I', 'numerical residual of the sine identity')
    if r:
        e.le('(lam a - b)/(lam b - a) - sin pi(nu - i tau)/sin pi(nu + i tau), numerical residual', r, 1e-100)
    e.headline = 'both loop orderings; sine identity holds'


def ev_exact_check(e):
    r = e.rows(r'\|X_dual/X / formula - 1\| = {N}', 'exchange-law lines', 10)
    if r:
        mx = e.le('|K(image)/K(physical) / exchange law - 1|: s = 0, -1, -2; l <= 3; both eta', [x[0] for x in r],
                  1e-55)
        e.headline = f'law holds to {fmtg(mx, 2)} (10 cases)'
    r = e.rows(r'\(\|X_dual/X - 1\| = {N}\)', 'non-triviality lines', 10)
    if r:
        e.ge('|K(image)/K(physical) - 1| (the image response differs from the physical one)', [x[0] for x in r], 1e-12)


def ev_test1(e):
    r = e.rows(r'nu_dual - nu_phys = {N}', 'nu lines', 9)
    if r:
        mx = e.le('|nu(image) - nu(physical)|, s = 0 and -2, 9 cases', [x[0] for x in r], 1e-50)
        e.headline = f'nu invariant to {fmtg(mx, 2)}'
    r = e.rows(r'\|Sfar_dual/Sfar - 1\| = {N}', 'far-zone lines', 9)
    if r:
        mx = e.le('|S_far(image)/S_far(physical) - 1|', [x[0] for x in r], 1e-45)
        e.headline += f'; far zone to {fmtg(mx, 2)}'


def ev_test3(e):
    r = e.rows(r'\|nu\(\+\)-nu\(-\)\|={N}', 'nu lines', 3)
    if r:
        e.le('super-extremal q = 1.5, 3, 2 (s = 0): |nu(+) - nu(-)|', [x[0] for x in r], 1e-50)
    r = e.rows(r'\|X\(-\)/X\(\+\)/Phi - 1\| = {N}', 'law lines', 3)
    if r:
        mx = e.le('super-extremal: |K(-)/K(+) / exchange law - 1|', [x[0] for x in r], 1e-45)
        e.headline = f'law to {fmtg(mx, 2)} at q > 1'
    r = e.rows(r'\|Sfar\(-\)/Sfar\(\+\)-1\| = {N}', 'far-zone lines', 3)
    if r:
        e.le('super-extremal: |S_far(-)/S_far(+) - 1|', [x[0] for x in r], 1e-45)
    r = e.rows(r'X_odd/X_even = {C}\s+delta/u = {N}', 'classical-like lines', 12)
    if r:
        e.near('a*omega fixed, eps -> 0: Re(X_odd/X_even) vs delta/u (12 cases)', [(x[0].real, x[1]) for x in r], 1e-3,
               rel=True)


def ev_monodromy(e):
    r = e.first(r'\|T_inf - T_\+ T_-\^above\| = {N}\s+\|T_inf - T_-\^below T_\+\| = {N}\s+'
                r'\(wrong orders:\s+{N}\s+{N}\s*\)', 'loop-product line')
    if r:
        e.le('|T_inf - T_+ T_-(above)| and |T_inf - T_-(below) T_+| (entries ~ 1e4, RK4)', r[:2], 1e-2)
        e.ge('wrong orderings T_-(above) T_+ and T_+ T_-(below) (negative control)', r[2:], 1e3)
    a = e.first(r'above r\+ -> {C}\s+below r\+ -> {C}', 'lasso line')
    b = e.first(r'eta=\+1\s+-> {C}\s+eta=-1 -> {C}', 'MST image line')
    if a and b:
        m1 = e.near('K_-/K_+ from the lasso above r_+ vs the MST image with eta = +1', [(a[0], b[0])], 1e-6)
        m2 = e.near('K_-/K_+ from the lasso below r_+ vs the MST image with eta = -1', [(a[1], b[1])], 1e-6)
        e.ge('|above - below| (the two lassos select different images)', [a[0] - a[1]], 0.1)
        e.headline = f'lassos reproduce eta = +1, -1 to {fmtg(max(m1, m2), 2)}'


# -- Classical limit, H_odd and R_2m
def ev_classical(e):
    r = e.rows(r'eps=0\.005: .*?X_odd/X_even={C}\s+delta/u={N}', 'eps = 0.005 lines', 2)
    if r:
        mx = e.near('Re(X_odd/X_even) vs delta/u at eps = 0.005 (m = 2, 1)', [(x[0].real, x[1]) for x in r], 1e-6,
                    rel=True)
        e.headline = ('Re X_odd/X_even = delta/u to all printed digits' if mx == 0 else
                      f'X_odd/X_even -> delta/u (rel. diff {fmtg(mx, 2)})')
    r = e.rows(r'eps={N}: .*?X_odd/X_even={C}', 'all lines', 8)
    if r:
        ratios = [abs(r[k + 2][1].imag / r[k + 3][1].imag) for k in (0, 4)]
        e.within('Im(X_odd/X_even) ~ eps^2: ratio of the values at eps = 0.01 and 0.005', ratios, 3.8, 4.2, d=5)


def ev_polar(e):
    r = e.first(r'C5 = {C}\s+expected {C}', 'C5 line')
    if r:
        mx = e.near('Im C5 vs -64/75, Eq. (S4) (l = m = 2)', [(r[0].imag, r[1].imag)], 1e-5, rel=True)
        e.le('|Re C5| (vanishes in the limit)', [r[0].real], 1e-4)
        e.headline = f'C5 = -64i/75 to {fmtg(mx, 2)}'
    r = e.first(r'kappa6 = C6/C5 = {C}', 'C6/C5 line')
    if r:
        e.near('Re C6/C5 vs 20/9, Eq. (S7) (m = 2)', [(r[0].real, 20 / 9)], 1e-2)


def ev_polar_m(e):
    r = e.rows(r's={N} l={N} m={N}: C5 = {C}\s+closed form {C} \| C6/C5 = {C}\s+expected {N}', 'table rows', 6)
    if r:
        mx = e.near('C5 vs closed form (S4): s = -2 (m = 1, 0, -1, -2), s = 0 (m = 2, -1)', [(x[3], x[4]) for x in r],
                    2e-4, rel=True)
        e.near('Re C6/C5 vs 10m/9 (s = -2) and 0 (s = 0), two-point fit in a*omega', [(x[5].real, x[6]) for x in r],
               0.05)
        e.headline = f'C5 to {fmtg(mx, 2)}; slopes to 0.05'


def ev_slopes(e):
    r = e.rows(r'^m={N}: C5={C}\s+C6/C5={C}\s+\(10m/9={N}\)', 'table rows', 4)
    if r:
        e.near('C5 vs -(-1)^m (8/225) i (2+m)!(2-m)!, Eq. (S4)', [(x[1], _c5_closed(int(x[0]))) for x in r], 1e-4,
               rel=True)
        mx = e.near('Re C6/C5 vs 10m/9, Eq. (S7), m = 1, 0, -1, -2 (quartic fit in a*omega)',
                    [(x[2].real, x[3]) for x in r], 0.02)
        e.headline = f'slopes 10m/9 to {fmtg(mx, 2)}'


def ev_verify_a7(e):
    r = e.rows(r'm={N} c={N}: ratio={C}\s+series={N}\s+\(ratio-series\)/c\^3={C}', 'l = 2 lines', 9)
    if r:
        e.le('|(H_odd/H_0 - R_2m)/c^3|, Eq. (S7), m = 2, 0, -2, c = 0.01, 0.02, 0.03', [x[4].real for x in r], 5)
        var = [abs(r[k + 2][4].real - r[k][4].real) / abs(r[k + 2][4].real) for k in (0, 6)]
        e.le('relative variation of (H_odd/H_0 - R_2m)/c^3 over c = 0.01...0.03, m = +-2 (remainder ~ c^3)', var,
             0.05)
    r = e.rows(r'm={N}: MST -> {C}\s+closed form {C}', 'l = 3 lines', 4)
    if r:
        mx = e.near('l = 3 leading coefficient vs closed form, m = 2, 0, -1, 3', [(x[1], x[2]) for x in r], 2.5e-4,
                    rel=True)
        e.headline = f'R_2m remainder ~ c^3; l = 3 to {fmtg(mx, 2)}'


# -- Polar amplitude
def ev_analyze3(e):
    r = e.first(r'worst relative \|G\| at root: {N};\s+max \|x\(eta=\+1\)-x\(eta=-1\)\| = {N}', 'root line')
    if r:
        e.le('worst relative |G| at the root nu (continued-fraction equation)', [r[0]], 1e-40)
        e.le('max |x(eta = +1) - x(eta = -1)| (nu is exchange-invariant)', [r[1]], 1e-40)
    r = e.first(r'max \|Re s1\| = {N}', 's1 line')
    if r:
        e.le('max |Re s_1|', r, 1e-13)
    r = e.first(r'max \|Re T_n\| \(n<=8\): {N}', 'T_n line')
    if r:
        e.le('max |Re T_n|, n <= 8 (purely imaginary Taylor coefficients)', r, 1e-13)
    e.headline = 'Taylor coefficients of the l = 2...19 MST data'


def ev_fit_contact2(e):
    r = e.rows(r'n=\s*{N}: eta=\+1 residual\s+{N}\s+eta=-1 contact size\s+{N}\s+ratio {N}', 'direct-test lines', 7)
    if r:
        mx = e.le('eta = +1 residual relative to the eta = -1 contact part, n = 5...11', [x[3] for x in r], 1e-14)
        e.headline = f'eta = +1 has no contact part ({fmtg(mx, 2)})'


def ev_closedform(e):
    m = e.rows(r'n=\s*{N} eta=-1\s+L<={n}:\s+max deviation\s+{N}\s+\(data scale {N}, rel {N}\)', 'eta = -1 lines', 9)
    p = e.rows(r'n=\s*{N} eta=\+1\s+L<={n}:\s+max deviation\s+{N}\s+\(data scale {N}, rel {N}\)', 'eta = +1 lines', 9)
    if m:
        mx = e.le('eta = -1: relative deviation of MST from the closed form (S12), O(a^5)...O(a^10)',
                  [x[3] for x in m if x[0] <= 10], 5e-16)
        e.le('eta = -1: relative deviation, O(a^11)...O(a^13)', [x[3] for x in m if x[0] > 10], 2e-15)
        e.headline = f'closed form to {fmtg(mx, 2)} through O(a^10)'
    if p:
        e.le('eta = +1 (pure exponential exp(-2cv)): relative deviation, O(a^5)...O(a^10)',
             [x[3] for x in p if x[0] <= 10], 5e-16)
        e.le('eta = +1: absolute deviation, O(a^5)...O(a^13)', [x[1] for x in p], 1e-14)


def ev_higher_orders(e):
    e.has('elementary form (S13) equals the series', 'elementary form - series: 0')
    e.has('odd part of f^(eta) equals f^odd', 'odd part of f^(eta) - fodd: 0')
    e.has('odd orders b_1, b_3, ..., b_9 have the factor 3 + cos(Theta)', 'odd b_k at x=-3: [0, 0, 0, 0, 0]')
    e.has('b_3 of Eq. (S14)', 'b_3(x) = (x + 3)*(3*x**2 + 14*x + 23)/42')
    e.has('[c^5] of Eq. (17)', '[c^5] f^odd/f0 = -64*v/15')
    e.has('[c^6] of Eq. (17)', '[c^6] f^odd/f0 = 256*v*(v - 2)/45')
    e.has('[c^7] of Eq. (17)', '[c^7] f^odd/f0 = -128*v*(11*v**2 - 40*v + 40)/315')
    e.has('[c^8] of Eq. (17)', '[c^8] f^odd/f0 = 256*v*(v - 2)*(3*v**2 - 10*v + 10)/315')
    e.headline = 'c^5...c^8 of Eq. (17); b_0...b_9 of (S14)'


def ev_a7_exact(e):
    e.has('R_2m through O(c^2) from the limiting MST recurrence, Eq. (S7)',
          'check vs 1 + (10m/9)c + 5(475m^2+261)c^2/3969 : 0')
    e.has('normalization of the spin-weighted basis', 'norm check: [1, 1, 1, 1, 1]')
    e.has('O(c) bracket of Eq. (17), b_1 = 2(x + 3)/3', 'O(c):   2*(x + 3)/3')
    e.has('O(c^2) bracket of Eq. (17), b_2', 'O(c^2): 11*x**2/42 + 29*x/21 + 13/6')
    e.headline = 'R_2m and the O(c), O(c^2) brackets exact'


def ev_a7_polar(e):
    e.count('c^5, c^6, c^7: closed form expanded minus prediction of Eq. (17) = 0', r'^c\^[567]:.*difference: 0\s*$', 3)
    e.headline = 'difference 0 at c^5, c^6, c^7'


# -- Segment representation
def ev_segment_born(e):
    r = e.rows(r'l=2 c={N}: ODE residual {N}; C2/\(A\^2 i/2a\) = {C}; C1 = {C}; S0\(M=0\) = {C}; '
               r'endpoint identity {N}', 'internal-check lines', 3)
    if r:
        e.le('residual of the M = 0 radial equation', [x[1] for x in r], 1e-30)
        e.le('double pole vs boundary Wronskian: |C2/(A^2 i/2a) - 1|', [x[2] - 1 for x in r], 1e-30)
        e.le('residue of the simple pole, |C1|', [x[3] for x in r], 1e-30)
        e.le('|S(M = 0) - 1|', [x[4] - 1 for x in r], 1e-25)
        e.le('endpoint identity', [x[5] for x in r], 1e-25)
    r = e.rows(r'ratio-1 = {C}', 'comparison lines', 3)
    if r:
        mx = e.le('|segment formula (S16) / MST - 1|, l = 2, c = 0.1, 0.2, 0.05 (MST printed to 12 digits)',
                  [x[0] for x in r], 1e-9)
        e.headline = f'(S16) vs MST to {fmtg(mx, 2)}'


def ev_vs_production(e):
    r = e.rows(r'\|rel diff\| {N}', 'table rows', 10)
    if r:
        mx = e.le('|segment (S16) / MST - 1|, l = 2...6 at c = 0.1997, 0.1631 (110-digit data)', [x[0] for x in r],
                  1e-15)
        e.headline = f'l = 2...6 to {fmtg(mx, 2)}'


def ev_nonlinear(e):
    r = e.rows(r'\( /c\^3 = {C} \)', 'nonlinear-term lines', 4)
    if r:
        e.within('relative change of the l = 2 exchange-odd element from the terms nonlinear in K, divided by c^3',
                 [x[0].real for x in r], -1.1, -0.7)
        e.near('same at c = 0.025 vs -0.85', [(r[0][0].real, -0.85)], 0.01)
        e.headline = f'nonlinear terms {fmtg(r[0][0].real, 2)} c^3 at c = 0.025'
    r = e.rows(r'beta K\(\+1\)/c\^3 = {C}', 'beta K lines', 4)
    if r:
        e.within('beta K/c^3 for l = 2', [x[0].real for x in r], 0.3, 0.7)
    r = e.rows(r'\|Sfar\+/Sfar- - 1\|={N}', 'far-zone lines', 4)
    if r:
        e.le('|S_far(+)/S_far(-) - 1|', [x[0] for x in r], 1e-20)


# -- Physical-region double expansion
def ev_rich_phys(e):
    gp = e.first(r'sign \+1: Gamma\(q=0\) at eps\^6 = {N}', 'sign +1 line')
    gm = e.first(r'sign -1: Gamma\(q=0\) at eps\^6 = {N}', 'sign -1 line')
    if gp and gm:
        e.near('q = 0, order eps^6: Gamma_+ - Gamma_- vs Page 4/225', [(gp[0] - gm[0], PAGE)], 1e-10, rel=True)
        e.near('Gamma_+ vs 4/225 + 1/4096 (helicity-preserving normalization)', [(gp[0], PAGE + 1 / 4096)], 1e-10,
               rel=True)
    sec = e.section('duality-odd (A+ - A-)/2 at eps^6', 'exchange-odd table')
    if sec is not None:
        r = e.first(r'p= 6:\s+{N}', 'p = 6 row', text=sec)
        if r:
            mx = e.near('exchange-odd q^6 coefficient vs rho_2/64 of arXiv:2610.40101', [(r[0], RHO2_64)], 1e-13,
                        rel=True)
            e.headline = f'[q^6] = rho_2/64 to {fmtg(mx, 2)}; Page 4/225'


def ev_b6(e):
    r = e.rows(r'm={N}: .*?naive two-radius {N}\s+joint bookkeeping {N}', 'table rows', 5)
    if r:
        mx = e.le('B6 vs i Pi_m/450 with joint log bookkeeping, m = 2...-2, p <= 9 (max rel. deviation)',
                  [x[2] for x in r], 5e-10)
        e.within('same with the naive two-radius fit (for comparison)', [x[1] for x in r], 1e-6, 1e-3)
        e.headline = f'B6 = i Pi_m/450 to {fmtg(mx, 2)}'


# -- Exchange-odd eps^6
def ev_d6_closed(e):
    r = e.rows(r'^\s*{N}\s+{N}\s+{N}\s+{N}\s*$', 'table rows', 11)
    if r:
        mx = e.le('Eq. (10) - exchange-odd MST Taylor coefficient, q^0...q^10', [x[3] for x in r], 1e-15)
        e.headline = f'Eq. (10) to {fmtg(mx, 2)}'
    r = e.first(r'rho_2/64: {N} {N}', 'q^6 line')
    if r:
        e.near('q^6 coefficient vs rho_2/64 of arXiv:2610.40101', [(r[0], r[1])], 1e-15, rel=True)


def ev_d6_m(e):
    r = e.rows(r'rel\.diff {N}\s+\(max \|D_k<6\| {N}\)', 'table rows', 8)
    if r:
        mx = e.le('D6 from MST vs Eq. (10), m = 2, 1, -1, -2 at q = 0.3, 0.7 (rel. diff)', [x[0] for x in r], 1e-12)
        e.le('lower orders |D_k|, k < 6', [x[1] for x in r], 1e-15)
        e.headline = f'Eq. (10) for m != 0 to {fmtg(mx, 2)}'


def ev_d6_m0(e):
    r = e.rows(r'D4={C}\s+D5={C}\s+D6={C}\s+D7={C}\s+\| -kappa\^5/225 = {N}', 'table rows', 2)
    if r:
        mx = e.near('m = 0: D6 vs -kappa^5/225, q = 0.3, 0.7', [(x[2], x[4]) for x in r], 1e-10, rel=True)
        e.le('|D4|, |D5|', [x[0] for x in r] + [x[1] for x in r], 1e-20)
        e.headline = f'D6 = -kappa^5/225 to {fmtg(mx, 2)}'


def ev_d6_general(e):
    r = e.rows(r'm={N}:\s+{N}\s+rho_\|m\|/64 = {N}', 'q^6 rows', 5)
    if r:
        e.near('q^6 Taylor coefficient of the coth part vs rho_|m|/64, m = 2...-2', [(x[1], x[2]) for x in r], 1e-15,
               rel=True)
    r = e.rows(r"m={N}: q=1e2,1e3,1e4 -> \['{C}', '{C}', '{C}'\]\s+target {C}", 'large-q rows', 5)
    if r:
        mx = e.near('32 D6/q^5 at q = 1e4 vs (-1)^m (8/225) i (2+m)!(2-m)!', [(x[3], x[4]) for x in r], 1e-7,
                    rel=True)
        e.headline = f'[q^6] = rho_|m|/64; large q to {fmtg(mx, 2)}'


def ev_paper_numbers(e):
    r = e.first(r's1\(c=0\) = {C}\s+i\(2gE\+2ln2-37/12\) = {C}', 'line (1)')
    if r:
        e.near('(1) eps-constant s_1(c = 0) vs i(2 gamma_E + 2 ln 2 - 37/12)', [(r[0], r[1])], 1e-15)
    r = e.first(r'u0 = pi \+ i\(2gE \+ 2ln2 - 37/12 - nu2\) = {C}', 'u0 line')
    if r:
        e.near('(1) u0 vs pi + i(2 gamma_E + 2 ln 2 - 37/12 + 107/210)', [(r[0], complex(math.pi, IM_U0))], 1e-12)
    e.count('(2) Eq. (10): Bose-factor form = Boltzmann form (S22)', r'^\(2\) .*: OK\s*$', 1)
    r = e.first(r'pi\(tau\+eps\) = {N}\s+\(w - m OmegaH\)/\(2TH\) = {N}', 'line (3)')
    if r:
        e.near('(3) pi(tau + eps) = (omega - m Omega_H)/(2 T_H)', [(r[0], r[1])], 1e-12, rel=True)
    r = e.first(r'160\[q\^5\]F = {N}\s+16 b2 = {N}', 'line (4)')
    if r:
        e.near('(4) 160 [q^5] F = 16 b_2 of arXiv:2610.40101', [(r[0], r[1])], 1e-15, rel=True)
    r = e.first(r'160 F/q\^5 at q = 1e8 \(continued\): {N}', 'line (4), large q')
    if r:
        e.near('(4) 160 F/q^5 at q = 1e8 vs 64/15', [(r[0], 64 / 15)], 1e-12, rel=True)
    r = e.rows(r'\(5\) m={N}: \[q\^6\]D6 =\s+{N}\s+rho/64 =\s+{N} \| 32 D6/q\^5 -> {C}\s+target {C}', 'lines (5)', 5)
    if r:
        e.near('(5) [q^6] D6 = rho_|m|/64, m = 2...-2', [(x[1], x[2]) for x in r], 1e-15, rel=True)
        e.near('(5) 32 D6/q^5 at large q vs the classical a^5 coefficient', [(x[3], x[4]) for x in r], 1e-9, rel=True)
    r = e.rows(r'\(6\) m={N}: pole_m = {C}\s+x75/i = {C}', 'lines (6)', 4)
    if r:
        tgt = {2: 64, 1: -16, -1: 16, -2: -64}
        e.near('(6) 75 pole_m/i = 64, -16, 16, -64 for m = 2, 1, -1, -2', [(x[2], tgt[int(x[0])]) for x in r], 1e-9)
    e.headline = 'eps-constant, u0, Eq. (10), [q^5]C = 16 b_2, rho_|m|/64'


def ev_se_allm(e):
    r = e.rows(r'^\s*{N}\s+{N}i\s+{N}i\s+{N}\s', 'table rows', 5)
    if r:
        mx = e.le('super-extremal exchange-even c^5 (classical MST) - [32 i a_5 + pole_m], m = 2...-2',
                  [x[3] for x in r], 2e-9)
        e.headline = f'continuation per m to {fmtg(mx, 2)}'
    r = e.first(r'm=\+2: eta=\+1\s+{N}', 'per-branch line')
    if r:
        e.near('m = 2, eta = +1 classical c^5 coefficient / 32 vs -2.3457997302e-5 (S5(d))',
               [(r[0] / 32, -2.3457997302e-5)], 1e-14)


# -- Exchange-even eps^6
def ev_allm(e):
    r = e.rows(r'450\*lam_log={N}\s+450\*lam_psi={N}', 'free-fit lines', 3)
    if r:
        mx = e.near('free fit: 450 x coefficients of Pi_m ln(kappa) and Pi_m Re psi(1 + i varpi), m = 1, -1, -2',
                    [(v, -1) for x in r for v in x], 1e-8)
        e.headline = f'both coefficients -1/450 to {fmtg(mx, 2)}'
    r = e.rows(r'^\s*{N}\s+{N}\s+{N}\s+{N}\s+{N}\s*$', 'closed-form tables', 16)
    if r:
        e.le('residual of the closed form, Taylor coefficients p = 6...15, m = 1, -1, -2', [x[3] for x in r if x[0] >= 6],
             1e-11)
    sec = e.section('should terminate: polynomial', 'm = 0 table')
    if sec is not None:
        r = e.rows(r'^\s*{N}\s+{C}\s*$', 'm = 0 rows', 16, text=sec)
        if r:
            e.le('m = 0: Taylor coefficients p >= 5 vanish (E6 is a polynomial)', [x[1] for x in r if x[0] >= 5], 1e-15)


def ev_unified(e):
    r = e.rows(r'max p=6\.\.{n}: {N}', 'table rows', 5)
    if r:
        mx = e.le('Taylor coefficients 6 <= p <= 20 after subtracting Eq. (S18), m = 2, 1, -1, -2, 0',
                  [x[0] for x in r], 2e-11)
        e.headline = f'remainder below {fmtg(mx, 2)}'


def ev_uniqueness(e):
    sec = e.section('rel. difference (continued closed form vs direct super-extremal MST)', 'continuation table')
    if sec is not None:
        r = e.rows(r'^\s*{N}\s+{N}\s*$', 'continuation rows', 17, text=sec)
        if r:
            e.le('polynomial fixed at |q| <= 0.3, continued to 50 <= q <= 5000, vs direct MST',
                 [x[1] for x in r if x[0] >= 50], 4e-10)
    sec = e.section('rel. difference with joint bookkeeping on both sides', 'joint-bookkeeping table')
    if sec is not None:
        r = e.rows(r'^\s*{N}\s+{N}\s*$', 'joint-bookkeeping rows', 8, text=sec)
        if r:
            mx = e.le('same for 1.5 <= q <= 30, joint log bookkeeping on both sides', [x[1] for x in r], 5e-12)
            e.headline = f'continuation to {fmtg(mx, 2)} (q <= 30)'
    r = e.first(r'\|c_6\| <~ {N}', 'c_6 bound')
    if r:
        e.le('bound on an additional term c_6 q^6', r, 1e-17)
    r = e.first(r'\|c_7\| <~ {N}', 'c_7 bound')
    if r:
        e.le('bound on an additional term c_7 q^7', r, 1e-20)
    r = e.first(r'^\s*5\s+{C}\s+{C}\s+{N}\s*$', 'q^5 row')
    if r:
        e.near('q^5 coefficient (joint bookkeeping), Im, vs -2.3457997302e-5 (S5(d))',
               [(r[0].imag, -2.3457997302e-5)], 1e-14)


def ev_phys_hi_p(e):
    r = e.first(r'p= 6: even A6 = {C}\s+odd A6 = {C}', 'p = 6 row')
    if r:
        mx = e.near('exchange-odd q^6 coefficient (one radius, no extrapolation) vs rho_2/64', [(r[1].real, RHO2_64)],
                    1e-8, rel=True)
        e.headline = f'coefficients to q^24; [q^6] odd to {fmtg(mx, 2)}'
    e.rows(r'p=\s*{N}: even A6 = {C}', 'rows p = 0...24', 25)


def ev_e6_fit(e):
    r = e.first(r'lambda \* 450 = {N}', 'lambda line')
    if r:
        mx = e.near('fitted coefficient of the psi + log term, 450 lambda, vs -1', [(r[0], -1)], 1e-8)
        e.headline = f'450 lambda = -1 to {fmtg(mx, 2)}'
    r = e.rows(r'^\s*{N}\s+{N}\s+{N}\s+{N}\s*$', 'table rows', 12)
    if r:
        e.le('fit residuals of the odd coefficients p = 7...23', [x[3] for x in r if x[0] >= 7], 1e-11)


def ev_e6_check_se2(e):
    r = e.rows(r'^\s*{N}\s+{C}\s+{C}\s+{N}\s*$', 'table rows', 17)
    if r:
        mx = e.le('continued closed form vs direct super-extremal MST, 1.5 <= q <= 5000 (naive log bookkeeping)',
                  [x[3] for x in r], 2e-6)
        e.headline = f'continuation to {fmtg(mx, 2)}'


def ev_e6_identities(e):
    a = re.search(r'160 \* \[q\^5\] F \(Taylor, physical\)\s+= (\S+)', e.out)
    b = re.search(r'16 b2 \(2610\.40101, eq\. 4\.2\)\s+= (\S+)', e.out)
    if a and b:
        e.add('160 [q^5] F = 16 b_2 of arXiv:2610.40101 (all printed digits)', a.group(1), b.group(1),
              a.group(1) == b.group(1))
    else:
        e.add('160 [q^5] F = 16 b_2 of arXiv:2610.40101', 'not found', '', False)
    r = e.first(r'q=1\.0e\+7: {N}', 'q = 1e7 line')
    if r:
        mx = e.near('160 F/q^5 at q = 1e7 vs 64/15', [(r[0], 64 / 15)], 1e-11, rel=True)
        e.headline = f'[q^5]C = 16 b_2; C/q^5 -> 64/15 ({fmtg(mx, 2)})'


def ev_se_taylor_l2(e):
    r = e.rows(r'^c\^{N}: eta-even {C}\s+eta-odd {C}', 'table rows', 9)
    if r:
        d = {int(x[0]): (x[1], x[2]) for x in r}
        e.le('exchange-odd coefficients of c^0...c^4', [d[k][1] for k in range(5)], 1e-25)
        e.near('exchange-odd c^5: Im vs -64/75', [(d[5][1].imag, -64 / 75)], 1e-14)
        e.near('exchange-odd c^6: Im vs -256/135', [(d[6][1].imag, -256 / 135)], 1e-14)
        e.headline = 'odd part starts at c^5; c^5, c^6 exact'


# -- eps^7 operator and a^6
def ev_d7_phys(e):
    r = e.first(r'poly_7\(q\):\s+residual {N}', 'fit line')
    if r:
        e.le('fit residual of D7 - (1/kappa) dD6/dtau0 = (u0 + u1 q + u2 q^2 + u3 q^3) D6 + poly_7(q)', r, 1e-12)
    u = [e.first(rf'^\s*u{k} = {{C}}', f'u{k} line') for k in range(4)]
    if all(u):
        u0, u1, u2, u3 = (x[0] for x in u)
        e.near('Re u0 vs pi (fit with u3 free)', [(u0.real, math.pi)], 1e-4)
        e.near('Im u0 vs 2 gamma_E + 2 ln 2 - 37/12 - nu2', [(u0.imag, IM_U0)], 1e-9)
        mx = e.near('u1 vs (5/9) m = 10/9', [(u1, 10 / 9)], 1e-8)
        e.le('|u2|, |u3| (no q^2, q^3 terms)', [u2, u3], 1e-9)
        e.headline = f'u1 = 10/9 to {fmtg(mx, 2)}; no q^2, q^3'
    m = re.search(r'check D6 closed form vs data \(p=0\.\.8\): \[(.*?)\]', e.out)
    if m:
        e.le('D6 closed form vs data, p = 0...8', complexes(m.group(1)), 1e-11)


def ev_d7_fit2(e):
    r = e.rows(r'm={N}: .*?u0 = {C}\s+u1 = {C}\s+u2 = {C}\s+resid {N}\s+\[\(5/9\)m = {N}\]', 'table rows', 5)
    if r:
        e.near('Re u0 vs pi, m = 2, 1, -1, -2, 0', [(x[1].real, math.pi) for x in r], 1e-8)
        e.near('Im u0 vs 2 gamma_E + 2 ln 2 - 37/12 - nu2', [(x[1].imag, IM_U0) for x in r], 1e-11)
        mx = e.near('u1 vs (5/9) m', [(x[2], x[5]) for x in r], 1e-9)
        e.le('|u2| (no q^2 term)', [x[3] for x in r], 1e-9)
        e.le('fit residual', [x[4] for x in r], 1e-11)
        e.headline = f'u0, u1 = 5m/9 ({fmtg(mx, 2)}) for all m'


def ev_e7_phys(e):
    r = e.first(r'B6 check \(p=1,3\):\s+{C}\s+{C}', 'B6 line')
    if r:
        e.near('B6 at p = 1, 3 vs 4i/225, 12i/225', [(r[0], 4j / 225), (r[1], 12j / 225)], 1e-10)
    r = e.first(r'fit - closed: {C}', 'u0 line')
    if r:
        mx = e.le('u0 (fit) - [pi + i(2 gamma_E + 2 ln 2 - 37/12 - nu2)]', r, 1e-7)
        e.headline = f'u0 to {fmtg(mx, 2)}'
    r = e.first(r'A7 fit:\s+u0 = {C}\s+u1 = {C}\s+u2 = {C}\s+residual {N}', 'A7 fit line')
    if r:
        e.near('u1 vs (5/9) m = 10/9', [(r[1], 10 / 9)], 1e-7)
        e.le('|u2| (no q^2 term)', [r[2]], 1e-9)
        e.le('fit residual', [r[3]], 1e-11)
    r = e.rows(r'p={N}: {C}\s+{C}\s+diff {C}', 'B7 rows', 10)
    if r:
        e.le("B7 data vs -(i/450)[P'/kappa + (u0 + 5mq/9) P], p = 0...9", [x[3] for x in r], 5e-6)


def ev_e7_allm(e):
    r = e.rows(r'm={N}: u0 = {C} \(closed {C}\)\s+u1 = {C} \(5m/9 = {N}\)\s+resid {N}', 'm != 0 rows', 3)
    if r:
        m1 = e.near('u0 vs closed form, m = 1, -1, -2', [(x[1], x[2]) for x in r], 5e-8)
        m2 = e.near('u1 vs 5m/9', [(x[3], x[4]) for x in r], 2e-8)
        e.le('fit residual', [x[5] for x in r], 1e-10)
        e.headline = f'u0 to {fmtg(m1, 2)}, u1 = 5m/9 to {fmtg(m2, 2)}'
    r = e.first(r'm=\+0: .*?resid {N};\s+max\|B7 - pred\| = {N}', 'm = 0 line')
    if r:
        e.le('m = 0: residual of A7 - (1/kappa) dT6/dtau0 = polynomial', [r[0]], 1e-10)
        e.le('m = 0: |B7 - prediction|, B7 = -(2i/225) kappa^4', [r[1]], 1e-9)


def ev_se7_analyze(e):
    t = e.first(r'targets \(classical c\^6, Im\): eta-even {N}\s+eta-odd {N}', 'targets line')
    r = e.first(r'extrapolation q -> infinity \(5 points, q >= {N}\):\s+eta-even {N}\s+eta-odd {N}',
                '5-point extrapolation line')
    if t and r:
        m1 = e.near('eta-even: 64 A7/q^6 extrapolated from 300 <= q <= 5000 vs classical c^6', [(r[1], t[0])], 1e-7)
        m2 = e.near('eta-odd: same vs 256/135', [(r[2], 256 / 135)], 1e-8)
        e.near('classical eta-odd c^6 = -256/135', [(t[1], -256 / 135)], 1e-9)
        e.headline = f'limits commute at a^6 ({fmtg(m1, 2)}, {fmtg(m2, 2)})'


def ev_se7_compare(e):
    r = e.rows(r'^\s*{N}\s+{C}\s+{C}\s+{N} \|\s+{C}\s+{C}\s+{N}\s*$', 'table rows', 17)
    if r:
        big = [x for x in r if x[0] >= 50]
        m1 = e.le('eta-even, 50 <= q <= 5000: continued closed form vs direct MST (rel.)', [x[3] for x in big], 1e-8)
        m2 = e.le('eta-odd, 50 <= q <= 5000', [x[6] for x in big], 1e-9)
        e.le('both sectors, q <= 30 (less precise direct data)', [x[k] for x in r if x[0] < 50 for k in (3, 6)], 1e-3)
        e.headline = f'{fmtg(m1, 2)} (even), {fmtg(m2, 2)} (odd) for q >= 50'


def ev_a6_hyp(e):
    r = e.first(r"m=\+2:\s+64 X7_even/q\^6 -> \['{C}', '{C}', '{C}'\]", 'm = 2 line')
    if r:
        e.ge('negative control: shifted singular function, |64 X7/q^6 (m = 2, large q) - 256/135|',
             [r[2].imag - 256 / 135], 0.1)
        e.headline = f'hypothesis rejected ({fmtg(r[2].imag, 6)} vs 256/135)'


# -- Analytic control of the exchange-odd response
def ev_analytic_odd(e):
    r = e.rows(r'c1 = {C}\s+\(5mq/9 = {N}\) \| E-coefficient = {C}\s+\(pi - i nu2 = {C}\) \| law residual {N}',
               'table rows', 4)
    if r:
        e.near('normalization ratio: O(eps) coefficient c1 vs 5mq/9 (finite eps)', [(x[0].real, x[1]) for x in r],
               5e-3)
        mx = e.near('E coefficient vs pi - i nu2', [(x[2], x[3]) for x in r], 1e-5)
        e.le('exchange-law residual', [x[4] for x in r], 1e-55)
        e.headline = f'E coefficient pi - i nu2 to {fmtg(mx, 2)}'


def ev_series_ratio(e):
    r = e.rows(r'c1 = {C}\s+\(5mq/9 = {N}\)\s+c2 = {C}.*?\(nu-2\)/eps\^2 at eps=1e-3: {C}\s+-> nu3 ~ {C} '
               r'\(mq/12 = {N}\)', 'table rows', 4)
    if r:
        mx = e.near('X-ratio coefficient c1 vs 5mq/9', [(x[0].real, x[1]) for x in r], 2e-6)
        e.near('(nu - 2)/eps^2 at eps = 1e-3 vs nu2 = -107/210, four (m, q)', [(x[3].real, NU2) for x in r], 2e-4)
        e.near('nu3 vs mq/12', [(x[4].real, x[5]) for x in r], 1e-3)
        e.headline = f'5mq/9 to {fmtg(mx, 2)}; nu2 = -107/210'


def ev_eps4(e):
    r = e.rows(r'm={N}: E1 = {N} \(-4m/3 = {N}\)\s+E2 = {N}\s+\(\(10/189\)\(9-m\^2\)-1 = {N}\)', 'E1, E2 rows', 4)
    if r:
        e.near('spheroidal coefficient E1 vs -4m/3', [(x[1], x[2]) for x in r], 1e-10)
        mx = e.near('spheroidal coefficient E2 vs (10/189)(9 - m^2) - 1', [(x[3], x[4]) for x in r], 5e-9)
        e.headline = f'E1 exact, E2 to {fmtg(mx, 2)}; pole term absent'
    lists = {}
    for k in (3, 4, 5):
        for sg in ('+1', '-1'):
            m = re.search(rf'eps\^{k} sign {re.escape(sg)}: q\^0\.\.q\^5 = \[(.*?)\]', e.out)
            lists[(k, sg)] = complexes(m.group(1)) if m else []
    if all(len(v) == 6 for v in lists.values()):
        e.le('eps^4 coefficient of S: q^4 and q^5 terms (cubic in q: the pole term is absent)',
             [lists[(4, s)][p] for s in ('+1', '-1') for p in (4, 5)], 1e-15)
        d = [abs(a - b) / (abs(a) + 1e-8) for k in (3, 4, 5) for a, b in zip(lists[(k, '+1')], lists[(k, '-1')])]
        e.le('exchange-odd part of S through eps^5: |S(+) - S(-)|/(|S| + 1e-8), q^0...q^5', [max(d)], 1e-6)
    else:
        e.add('eps^3...eps^5 coefficient lists', 'not found', '6 entries each', False)


def ev_h_scaling(e):
    r = e.rows(r'\|H_even\|/eps\^4 = {N}\s+\|H_even\|/eps\^5 = {N}\s+\|H_odd\|/eps\^6 = {N}\s+'
               r'\|Sfar_\+ - Sfar_-\| = {N}', 'table rows', 3)
    if r:
        def spread(v):
            return (max(v) - min(v)) / max(v)
        e.le('variation of |H_even|/eps^4 over eps = 1e-3, 5e-4, 2.5e-4', [spread([x[0] for x in r])], 1e-2)
        e.le('variation of |H_odd|/eps^6', [spread([x[2] for x in r])], 1e-2)
        mx = e.le('|S_far(+) - S_far(-)| (far zone exchange-invariant)', [x[3] for x in r], 1e-60)
        e.headline = f'H_even ~ eps^4, H_odd ~ eps^6; far zone {fmtg(mx, 2)}'


def ev_absorption(e):
    r = e.rows(r'linear extrapolation eps -> 0: {N}\s+\|difference\| = {N}', 'extrapolation lines', 4)
    if r:
        mx = e.le('(Gamma_+ - Gamma_-)/eps^6 extrapolated vs -4 A6^-, (m, q) = (2, 0), (0, 0.3), (2, 0.3), (-1, 0.5)',
                  [x[1] for x in r], 5e-6)
        e.headline = f'Gamma_+ - Gamma_- = -4 eps^6 A6^- to {fmtg(mx, 2)}'
    r = e.first(r'm=\+2 q=0 eps=1e-3: .*?-4 A6\^- = {N}', 'q = 0 line')
    if r:
        e.near('q = 0: -4 A6^- vs Page 4/225', [(r[0], PAGE)], 1e-9, rel=True)


# -- General kinematics
def _kin_block(e, order):
    sec = e.section(f'=== O(a^{order}) ===', f'O(a^{order}) block')
    if sec is None:
        return [], []
    p = re.search(r'helicity \+1: .*?\[(.*?)\]', sec)
    m = re.search(r'helicity -1: .*?\[(.*?)\]', sec)
    return (complexes(p.group(1)) if p else []), (complexes(m.group(1)) if m else [])


def ev_kinematics(e):
    r = e.first(r'addition theorem check: {N}', 'addition-theorem line')
    if r:
        e.near('addition theorem for the spin-weighted harmonics', [(r[0], 1)], 1e-12)
    p5, m5 = _kin_block(e, 5)
    p6, m6 = _kin_block(e, 6)
    if len(p5) == 6 and len(p6) == 6:
        m1 = e.near('helicity +1, a^5: BHPT / contact-term ansatz at 6 random kinematic points', [(z, 1) for z in p5],
                    1e-25)
        m2 = e.near('helicity +1, a^6: BHPT / ansatz = -1 (sign convention a -> -a)', [(z, -1) for z in p6], 1e-25)
        e.headline = f'ratios 1 to {fmtg(max(m1, m2), 2)} (helicity +1)'
    else:
        e.add('helicity +1 ratios', 'not found', '6 values at a^5 and at a^6', False)
    if len(m5) == 6 and len(m6) == 6:
        e.ge('helicity -1 (negative control): |ratio - 1| at a^5 and |ratio + 1| at a^6',
             [z - 1 for z in m5] + [z + 1 for z in m6], 1e-3)
    else:
        e.add('helicity -1 ratios', 'not found', '6 values at a^5 and at a^6', False)


def ev_fitcoeffs(e):
    r = e.rows(r'^\s*c\d+\^\(\d\)\s+fitted =\s+{N}\s+table = {N}', 'coefficient rows', 9)
    if r:
        mx = e.near('nine exchange-odd coefficients: least-squares fit vs Table S1', [(x[0], x[1]) for x in r], 1e-15)
        e.headline = ('nine coefficients of Table S1, all 20 digits' if mx == 0 else
                      f'nine coefficients of Table S1 (diff {fmtg(mx, 2)})')
    r = e.rows(r'O\(a\^[56]\): residual {N}', 'fit residuals', 2)
    if r:
        e.le('fit residuals at a^5 and a^6', [x[0] for x in r], 1e-20)
    r = e.first(r'perturbed c4\^\(1\) by 0\.01 -> BHPT/ansatz = {C}', 'perturbation line')
    if r:
        e.ge('negative control: c4^(1) perturbed by 0.01, |BHPT/ansatz - 1|', [r[0] - 1], 1e-3)


def ev_kin_even(e):
    i = e.out.find('helicity -1')
    plus, minus = (e.out[:i], e.out[i:]) if i >= 0 else (e.out, '')
    r = e.rows(r'BHPT/ansatz\(c4, alpha=1\) = {C}', 'helicity +1 ratios', 8, text=plus)
    if r:
        e.near('helicity +1: exchange-odd BHPT / ansatz (c4, alpha = 1) at 8 points', [(x[0], 1) for x in r], 1e-25)
    m = re.search(r'fitted c3 .*?: \[(.*?)\]\s+residual (' + _NUM + ')', plus)
    if m:
        c3 = complexes(m.group(1))
        e.le('fit residual of the exchange-even a^5 contact term', [float(m.group(2))], 1e-25)
        mx = e.near('fitted c3 vs -(64/15, 16/3, 16/15): the alpha-tagged psi-pole terms', list(zip(c3, [-64 / 15, -16 / 3,
                    -16 / 15])), 1e-12)
        e.headline = f'c3 = alpha(64/15, 16/3, 16/15) ({fmtg(mx, 2)})'
    else:
        e.add('fitted c3 (helicity +1)', 'not found', '', False)
    m = re.search(r'residual (' + _NUM + ')', minus)
    if m:
        e.ge('helicity -1 (negative control): fit residual', [float(m.group(1))], 1e-3)
    else:
        e.add('helicity -1 fit (negative control)', 'not found', '', False)


def ev_kin_tags(e):
    targets = [('eta-odd  full', (64 / 15, 16 / 5, 4 / 15), '(64/15, 16/5, 4/15)'),
               ('eta-odd  m != 0', (64 / 15, 16 / 5, 0), '(64/15, 16/5, 0)'),
               ('eta-odd  m = 0', (0, 0, 4 / 15), '(0, 0, 4/15)'),
               ('eta-even psi-pole', (64 / 15, 16 / 3, 16 / 15), '(64/15, 16/3, 16/15)')]
    for label, tgt, txt in targets:
        m = re.search(re.escape(label) + r".*?-> \(\[([^\]]*)\], '(" + _NUM + r")'\)", e.out)
        if not m:
            e.add(label, 'not found', '', False)
            continue
        vals = _num_list(m.group(1))
        name = ' '.join(label.split())
        e.near(f'{name} part: coefficients vs {txt}', list(zip(vals, tgt)), 1e-10)
        e.le(f'{name} part: fit residual', [float(m.group(2))], 1e-25)
    r = e.first(r'untagged c3\^\(2\)=4/15 as l=2 pattern .*?residual {N}', 'untagged line')
    if r:
        e.ge('negative control: untagged c3^(2) = 4/15 as a pure l = 2 pattern, fit residual', r, 0.1)
    e.headline = 'tagged and untagged terms of Table S1'


def ev_kin_a6(e):
    lines = re.findall(r'^\s*(eta-\S+ .*?):\s+\[(.*)\]\s*$', e.out, re.M)
    if len(lines) < 4:
        e.add('ratio lines', f'{len(lines)} found', '4', False)
        return
    mx = 0
    for label, lst in lines:
        z = complexes(lst)
        sign = -1 if label.startswith('eta-odd') else 1
        if len(z) != 6:
            e.add(' '.join(label.split()), f'{len(z)} values', '6 values', False)
            continue
        mx = max(mx, e.near(f"{' '.join(label.split())}: ratio = {sign:+d}", [(v, sign) for v in z], 1e-25))
    e.headline = f'a^6 pole and polynomial parts, ratios +-1 to {fmtg(mx, 2)}'


# -- Unit-flux normalization
def ev_flux_check(e):
    rx = (r'^{N}\s+{N}\s+{N}\s*\|\s*{N}\s+{N}\s+{N}\s*\|\s*{N}\s*\|\s*{N}\s+{N}(?:\s*\(superradiant\))?\s*\|'
          r'\s*{N}\s+{N}\s*$')
    r = e.rows(rx, 'table rows', 17)
    if r:
        e.le('MST continued-fraction residual g(nu)', [x[3] for x in r], 1e-55)
        e.le('B[R_C^(-nu-1)]/B[R_C^nu] relations (incoming and outgoing)', [x[4] for x in r], 1e-55)
        e.le('factorized vs direct S_P', [x[5] for x in r], 1e-55)
        e.le('|C|^2 of S10 vs D^2 + 144 M^2 omega^2', [x[6] for x in r], 1e-55)
        mx = e.le('flux balance |S_P|^2 + P_H - 1', [x[7] for x in r], 1e-55)
        e.le('|S/S_P| - D/|C| (D-normalized element)', [x[10] for x in r], 1e-55)
        n_sr = sum(1 for x in r if x[8] < 0)
        e.add('superradiant cases included (P_H < 0)', f'{n_sr} of {len(r)}', '>= 1', n_sr >= 1)
        e.headline = f'flux balance to {fmtg(mx, 2)}'


def ev_highfreq(e):
    r = e.rows(r'\|g\|={N}\s+R=\s*{N}\s+P_H=\s*{N}.*?R\+P_H-1=\s*{N}', 'table rows', 3)
    if r:
        mx = e.le('flux balance R + P_H - 1 at eps = 0.6', [x[3] for x in r], 1e-55)
        e.le('continued-fraction residual |g|', [x[0] for x in r], 1e-55)
        e.headline = f'flux balance to {fmtg(mx, 2)}'


def ev_complexnu(e):
    r = e.rows(r'nu=\s*{C}\s+\|g\|={N}\s+R=\s*{N}\s+P_H=\s*{N}.*?R\+P_H-1=\s*{N}', 'table rows', 6)
    if r:
        mx = e.le('flux balance R + P_H - 1 on the complex-nu branch', [x[4] for x in r], 1e-55)
        e.le('continued-fraction residual |g|', [x[1] for x in r], 1e-55)
        e.near('Re nu half-integer on the complex branch', [(x[0].real, math.floor(x[0].real) + 0.5) for x in r], 1e-12)
        amp = max(-x[3] for x in r)
        e.add('largest superradiant amplification -P_H', fmtg(amp, 3), '> 0', amp > 0)
        e.headline = f'flux balance to {fmtg(mx, 2)}; amplification {fmtg(100 * amp, 3)}%'


def ev_decomp(e):
    r = e.rows(r'^{N}\s+{N}\s+{N}\s*\|\s*{N}\s+{N}\s+\(\s*{N}\)', 'table rows', 12)
    if r:
        mx = e.le('Delta_odd/eps^6 vs (pi/225) kappa^5 varpi (1 + varpi^2)(4 + varpi^2) n_eta (rel.), both eta',
                  [x[5] for x in r], 5e-6)
        e.headline = f'Delta_odd to {fmtg(mx, 2)}; half of Page at q = 0'
    r = e.rows(r'{N}\s+\(k/\(1\+k\) = {N}\)', 'm = 0 and q = 0 rows', 6)
    if r:
        e.near('Delta_odd/P_H vs kappa/(1 + kappa) (m = 0; 1/2 at q = 0)', [(x[0], x[1]) for x in r], 1e-6)
    r = e.rows(r'\|\s*{N} id {N}\s*$', 'flux columns', 12)
    if r:
        e.le('flux balance and the decomposition identity', [v for x in r for v in x], 1e-55)


# ----------------------------------------------------------------------------------------------------
# Registry of the checks
# ----------------------------------------------------------------------------------------------------
@dataclass
class Group:
    key: str
    title: str          # LaTeX
    paper: str
    intro: str          # LaTeX


@dataclass
class Check:
    group: str
    folder: str         # '.' for the top-level folder
    script: str
    paper: str          # paper items, plain text
    title: str          # LaTeX
    purpose: str        # LaTeX
    evaluate: object = None
    args: tuple = ()
    est: float = 10.0   # running time in seconds, measured with two scripts on two cores
    cpus: int = 1       # worker processes used by the script
    long: bool = False  # skipped by --quick
    deps: tuple = ()    # checks that write files read by this one
    needs: tuple = ()   # Python modules besides mpmath and sympy
    files: tuple = ()   # data files that are not shipped or are produced by a driver
    info: bool = False  # informational: no pass/fail criterion

    @property
    def id(self):
        return self.script if self.folder == '.' else f'{self.folder}/{self.script}'


GROUPS = [
    Group('law', r'Exchange law and monodromy', 'Eq. (6), Fig. 1; S1, S9(a), S9(b)',
          r'The exchange law, Eq.~(6), relates the response $\mathcal K_-$ selected at the inner horizon (the '
          r'image solution, continuation label $\eta=\pm1$) to the physical response $\mathcal K_+$. In all '
          r'numerical checks the image is computed directly from MST with $\kappa\to-\kappa$ and '
          r'$\log(-\kappa)=\log\kappa+i\pi\eta$, without using the law.'),
    Group('classical', r'Classical limit: $H^{\rm odd}_{\ell m}$ and $R_{2m}$',
          'S2, S3; Eqs. (S3), (S4), (S7); S9(i)',
          r'The classical limit is $\epsilon\to0$ at fixed $c=a\omega$, on the super-extremal branches '
          r'$\kappa=\eta\,iq\sqrt{1-1/q^2}$. The exchange-odd element starts at $c^{2\ell+1}$ with the '
          r'closed-form coefficient of Eq.~(S4), and $R_{2m}=1+\frac{10m}{9}c+\frac{5(475m^2+261)}{3969}c^2'
          r'+O(c^3)$, Eq.~(S7).'),
    Group('polar', r'Polar classical amplitude: Eqs.~(17) and (18)', 'Eqs. (17), (18); S4, Eqs. (S12)-(S14)',
          r'For incidence along the spin the classical amplitude is known in closed form, Eq.~(18) '
          r'[Supplemental Eq.~(S12)]. These checks compare it with full MST (110 digits, $\ell=2\ldots19$, '
          r'both $\eta$) order by order in $a$, and expand it symbolically to obtain Eq.~(17), its $a^8$ term '
          r'and the coefficients of Eq.~(S14).'),
    Group('segment', r'Segment representation', 'S4, Eqs. (S15)-(S17)',
          r'In the classical limit the horizons sit at $r=\pm ia$, and the exchange-odd part of each '
          r'partial-wave element is a finite-part integral along the segment joining them plus an endpoint '
          r'term, Eq.~(S16).'),
    Group('physreg', r'Physical-region double expansion in $(\epsilon,q)$', 'S5; data behind Eqs. (10)-(12), (16)',
          r'Two-dimensional Cauchy integrals of the $\ell=2$, $m=2$ S-matrix element in $\epsilon$ (two radii) '
          r'and $q$, full MST at 90 digits, for the physical branch and its image. The data files are '
          r'produced by the production drivers (option \texttt{--production}).'),
    Group('odd6', r'Exchange-odd $\epsilon^6$ coefficient', 'Eq. (10); S9(d), S9(f)',
          r'The exchange-odd $\epsilon^6$ coefficient $D_6$, Eq.~(10), its Boltzmann form (S22), and its '
          r'continuation to super-extremal spin, where it gives the classical $a^5$ coefficients.'),
    Group('even6', r'Exchange-even $\epsilon^6$ function and $\mathcal C(q)$', 'Eqs. (12), (16); S5; S9(c), S9(e)',
          r'The exchange-even $\epsilon^6$ function, Eq.~(12): a polynomial of degree five plus '
          r'$-(i/450)P^{(m)}[\ln(\epsilon\kappa)+\mathrm{Re}\,\psi(1+i\tau_0)]$; the uniqueness of this form, '
          r'and the contact coefficient $\mathcal C(q)$ of Eq.~(16).'),
    Group('eps7', r'The $\epsilon^7$ operator and the $a^6$ coefficients', 'Eq. (S23); S9(g), S9(h)',
          r'At order $\epsilon^7$ both sectors follow from $\epsilon^6$ by the operator '
          r'$(1/\kappa)\partial_{\tau_0}+u_0+\frac{5m}{9}q$, with $u_0=\pi+i(2\gamma_E+2\ln2-\frac{37}{12}'
          r'-\nu_2)$ and $\nu_2=-107/210$, Eq.~(S23).'),
    Group('analytic', r'Analytic control of the exchange-odd response', 'S7; S9(j)',
          r'Full-MST checks of the analytic structure behind Eqs.~(10) and (11): the normalization of '
          r'$\delta^2\mathcal K_+$, the $\epsilon$ coefficient $\pi-i\nu_2$ with $\nu_2=-107/210$, the absence '
          r'of the near-zone pole term from $S$, and the absorption difference '
          r'$\Gamma_+-\Gamma_-=-4\epsilon^6A_6^-$.'),
    Group('kin', r'General kinematics', 'S6, Table S1',
          r'The exchange-odd BHPT amplitude in general kinematics against the contact-term ansatz of '
          r'Bautista, Guevara, Kavanagh and Vines, the nine exchange-odd coefficients of Table~S1, and the '
          r'tagged and untagged exchange-even terms. The helicity $-1$ comparisons are negative controls.'),
    Group('flux', r'Unit-flux normalization and the absorption split', 'Eq. (11); S9(k), S10',
          r'Unit-flux normalization of the MST far-zone factor, flux balance $|S_P|^2+\mathcal P_H=1$ '
          r'including superradiant cases and the complex-$\nu$ branch, and the exchange-parity split of the '
          r'absorption probability, Eq.~(11).'),
]

C = Check
CHECKS = [
    # -- exchange law
    C('law', '.', 'proof_check.py', 'Eq. (6), S1', r'Symbolic proof of the exchange law',
      r'SymPy: the monodromy algebra for both orderings of the loop product gives '
      r'$\mathcal K_-/\mathcal K_+=-\lambda^{\mp1}(\lambda a-b)/(\lambda b-a)$, and '
      r'$(\lambda a-b)/(\lambda b-a)=\sin\pi(\nu-i\tau)/\sin\pi(\nu+i\tau)$.', ev_proof_check, est=1),
    C('law', '.', 'exact_check.py', 'Eq. (6), S9(a)', r'Exchange law against full MST',
      r'Ratio of the image response to the physical one, both computed with MST at 60 digits, divided by '
      r'the exchange law, for $s=0,-1,-2$, $\ell\le3$, $q\le0.9$, $\epsilon\le0.6$ and both $\eta$.',
      ev_exact_check, est=11),
    C('law', '.', 'test1.py', 'S9(b)', r'Exchange invariance of $\nu$ and of the far-zone factor',
      r'The renormalized angular momentum $\nu$ and the far-zone factor are the same for the physical and '
      r'image solutions (Schwarzschild and Kerr, $s=0$ and $-2$).', ev_test1, est=17),
    C('law', '.', 'test3.py', 'Eq. (6), S9(a), S9(b)', r'Super-extremal spin and the classical-like regime',
      r'The exchange law and the invariances for $q=1.5$, $2$, $3$ ($s=0$, 50 digits); '
      r'$\mathrm{Re}\,\mathcal X_{\rm odd}/\mathcal X_{\rm even}\to\delta/u$ as $\epsilon\to0$ at fixed '
      r'$a\omega$.', ev_test3, est=17),
    C('law', 'figure_monodromy', 'monodromy_loops_check.py', 'Fig. 1, S1',
      r'The loops of Fig.~1 by direct integration',
      r'RK4 integration of the $s=-2$ radial Teukolsky equation ($\ell=m=2$, $q=0.6$, $\epsilon=0.3$) along '
      r'the loops of Fig.~1: $T_\infty=T_+T_-^{\rm above}=T_-^{\rm below}T_+$, and the lassos passing above '
      r'and below $r_+$ give the MST image with $\eta=+1$ and $\eta=-1$.', ev_monodromy, est=19,
      needs=('numpy',)),
    # -- classical limit
    C('classical', '.', 'classical.py', 'S2, Eq. (S3)', r'Structure of the response in the classical limit',
      r'For $s=0$, $c=0.3$ and $m=2,1$: $\mathcal X_{\rm odd}/\mathcal X_{\rm even}\to\delta/u$ as '
      r'$\epsilon\to0$, with an imaginary part vanishing as $\epsilon^2$.', ev_classical, est=13),
    C('classical', '.', 'polar.py', 'Eqs. (S4), (S7); S9(i)', r'Leading exchange-odd element for $\ell=m=2$',
      r'$H^{\rm odd}/(\epsilon c^5)$ from MST at $a\omega=0.004\ldots0.02$, extrapolated to '
      r'$C_5=-64i/75$; the slope $C_6/C_5=20/9$.', ev_polar, est=17),
    C('classical', '.', 'polar_m.py', 'Eqs. (S4), (S7); S9(i)', r'Leading elements and slopes for all $m$',
      r'$C_5$ against the closed form (S4) for $s=-2$ ($m=1,0,-1,-2$) and $s=0$ ($m=2,-1$); slopes '
      r'$C_6/C_5=10m/9$ for $s=-2$ and $0$ for $s=0$, from a two-point fit in $a\omega$ (accuracy about '
      r'$0.03$).', ev_polar_m, est=35),
    C('classical', '.', 'slopes.py', 'Eq. (S7), S9(i)', r'Slopes $10m/9$ of $R_{2m}$',
      r'Quartic fit in $a\omega$ of $H^{\rm odd}/(\epsilon c^5)$ at five values of $a\omega$ ($s=-2$, '
      r'$\ell=2$): $C_5$ and $C_6/C_5=10m/9$ for $m=1,0,-1,-2$.', ev_slopes, est=59, long=True),
    C('classical', '.', 'verify_a7.py', 'Eq. (S7), S9(i)', r'Remainder of $R_{2m}$ and the $\ell=3$ coefficients',
      r'$H^{\rm odd}/H_0$ against $R_{2m}$ through $O(c^2)$: the remainder is proportional to $c^3$. '
      r'Leading $\ell=3$ coefficients for $m=2,0,-1,3$ against the closed form.', ev_verify_a7, est=79,
      long=True),
    # -- polar amplitude
    C('polar', 'conservative_test', 'analyze3.py', 'Eq. (18), S4', r'Taylor coefficients of the polar MST data',
      r'Cauchy-integral Taylor coefficients in $c$ of the 110-digit MST data \texttt{polar\_run2.json} '
      r'($\ell=2\ldots19$, both $\eta$); writes \texttt{taylor2.json}, read by the next two scripts.',
      ev_analyze3, est=8),
    C('polar', 'conservative_test', 'fit_contact2.py', 'Eq. (18), S4', r'$\eta=+1$ is the pure exponential',
      r'Order by order, the $\eta=+1$ amplitude has no contact part: its Taylor coefficients are those of '
      r'the exponential $e^{-2cv}$.', ev_fit_contact2, est=3, deps=('conservative_test/analyze3.py',)),
    C('polar', 'conservative_test', 'closedform_check.py', 'Eq. (18), Eq. (S12)',
      r'Closed form (18) against MST through $O(a^{13})$',
      r'Partial-wave sums of the closed form against the MST Taylor coefficients, order by order from $a^5$ '
      r'to $a^{13}$, both $\eta$.', ev_closedform, est=6, deps=('conservative_test/analyze3.py',)),
    C('polar', 'conservative_test', 'higher_orders_eq17.py', 'Eq. (17), Eqs. (S13), (S14)',
      r'Expansion of the closed form: Eqs.~(17) and (S14)',
      r'SymPy: the elementary form (S13), the coefficients $b_0\ldots b_9$ of (S14), the factor '
      r'$3+\cos\Theta$ of the odd orders, and the $c^5\ldots c^8$ coefficients of $f^{\rm odd}/f_0$, '
      r'Eq.~(17).', ev_higher_orders, est=30),
    C('polar', '.', 'a7_exact.py', 'Eq. (17), Eq. (S7)', r'Exact $O(c^2)$ ingredients of Eq.~(17)',
      r'SymPy: $R_{2m}$ from the limiting MST recurrence, the spin-weighted spheroidal harmonic $S_{22}$ to '
      r'$O(c^2)$ and the $\ell=3$ term, assembled into the $O(c)$ and $O(c^2)$ brackets of Eq.~(17).',
      ev_a7_exact, est=4),
    C('polar', 'prl_checks', 'a7_polar_check.py', 'Eq. (17)', r'Closed form against Eq.~(17)',
      r'SymPy: the exchange-odd part of the closed form (S12) has the $c^5$, $c^6$, $c^7$ coefficients of '
      r'Eq.~(17).', ev_a7_polar, est=2),
    # -- segment representation
    C('segment', 'segment_test', 'segment_born.py', 'Eq. (S16)', r'Segment formula (S16) against MST',
      r'Internal checks of the construction (radial equation, double pole, residue, $S=1$ at $M=0$, '
      r'endpoint identity) and the $\ell=2$ element at $c=0.05$, $0.1$, $0.2$ against MST at '
      r'$\epsilon=10^{-12}$.', ev_segment_born, est=41),
    C('segment', 'segment_test', 'vs_production.py', 'Eq. (S16), S4', r'Segment formula for $\ell=2\ldots6$',
      r'$\ell=2\ldots6$ at $c=0.1997$ and $0.1631$ against the 110-digit MST data; the agreement is the '
      r'finite-$\epsilon$ systematic of the data.', ev_vs_production, est=146, long=True,
      files=('../conservative_test/polar_run2.json',)),
    C('segment', 'segment_test', 'nonlinear_check.py', 'S4', r'Terms nonlinear in $\mathcal K$',
      r'Size of the terms nonlinear in $\mathcal K$ in $S/F-1=iE\mathcal K/(1+\beta\mathcal K)$: they '
      r'change the $\ell=2$ exchange-odd element at relative order $c^3$ (about $-0.85c^3$ at $c=0.025$).',
      ev_nonlinear, est=13),
    # -- physical region
    C('physreg', 'physical_region_test', 'rich_phys.py', 'Eqs. (10), (11); S5', r'Richardson-improved coefficients',
      r"Richardson extrapolation between the two radii; Page's $4/225$ at $q=0$; the exchange-odd $q^6$ "
      r'coefficient equals $\rho_2/64$ of arXiv:2610.40101. Writes \texttt{l2\_rich.json}.', ev_rich_phys, est=9),
    C('physreg', 'physical_region_test', 'pslq2.py', 'S5', r'PSLQ identification (exploratory)',
      r'Exploratory PSLQ fits of individual coefficients in a basis of zeta values, superseded by the '
      r'closed forms of the paper (digamma function). Informational: no pass/fail criterion.', None, est=1,
      deps=('physical_region_test/rich_phys.py',), info=True),
    C('physreg', 'physical_region_test', 'analyze_phys.py', 'S5', r'Coefficients of the polar projection',
      r'Coefficients of the $L=2$ projection compared with arXiv:2610.40101 and with the super-extremal data. '
      r'Needs \texttt{phys\_run.json} (about 5\,MB, not shipped), produced by \texttt{run\_phys.py}. '
      r'Informational.', None, est=10, deps=('conservative_test/analyze3.py',), files=('phys_run.json',),
      info=True),
    C('physreg', 'prl_checks', 'b6_precision.py', 'S5, S9(c)', r'The logarithm $B_6=i\Pi_m/450$',
      r'Coefficient of $\epsilon^6\ln\epsilon$ for every $m$: the naive two-radius fit reaches '
      r'$10^{-5}$--$10^{-4}$, the joint log bookkeeping $10^{-10}$.', ev_b6, est=2),
    # -- exchange-odd eps^6
    C('odd6', 'continuation_test', 'd6_closed.py', 'Eq. (10), S9(d)', r'Eq.~(10) against the physical-region data',
      r'Taylor coefficients $q^0\ldots q^{10}$ of Eq.~(10) for $m=2$ against the exchange-odd MST '
      r'coefficients; $[q^6]=\rho_2/64$.', ev_d6_closed, est=1),
    C('odd6', 'continuation_test', 'd6_m.py', 'Eq. (10), S9(d)', r'Eq.~(10) for $m\neq0$',
      r'$D_6$ from MST (Cauchy integral in $\epsilon$) against Eq.~(10) for $m=2,1,-1,-2$ at $q=0.3$ and '
      r'$0.7$; the lower orders vanish. Uses two worker processes.', ev_d6_m, est=113, cpus=2, long=True),
    C('odd6', 'continuation_test', 'd6_m0.py', 'Eq. (10), S9(d)', r'Eq.~(10) for $m=0$',
      r'For $m=0$: $D_6=-\kappa^5/225$, and $D_4$, $D_5$ vanish.', ev_d6_m0, est=58),
    C('odd6', 'continuation_test', 'd6_general.py', 'Eq. (10), S9(d), S9(f)',
      r'$[q^6]D_6=\rho_{|m|}/64$ and the large-spin limit',
      r'The $q^6$ coefficients of the coth part equal $\rho_{|m|}/64$ of arXiv:2610.40101 for all $m$, and '
      r'$32D_6/q^5$ continued along the upper path tends to the classical $a^5$ coefficient '
      r'$(-1)^m\frac{8}{225}i(2+m)!(2-m)!$.', ev_d6_general, est=1),
    C('odd6', 'prl_checks', 'paper_numbers.py', 'Eqs. (10), (16), (S22), (S23)',
      r'Closed forms and identities of the main text',
      r'The $\epsilon$-constant $i(2\gamma_E+2\ln2-\frac{37}{12})$, $u_0$, Eq.~(10) in its Bose and '
      r'Boltzmann forms, $\pi(\tau+\epsilon)=(\omega-m\Omega_H)/2T_H$, $160[q^5]F=16b_2$, '
      r'$[q^6]D_6=\rho_{|m|}/64$ and the large-$q$ limits, and the pole terms ${\rm pole}_m$.',
      ev_paper_numbers, est=1),
    C('odd6', 'prl_checks', 'se_allm_analyze.py', 'S5(d), S9(f)', r'Classical $a^5$ coefficients for every $m$',
      r'Super-extremal exchange-even $c^5$ coefficient for every $m$ against $32ia_5^{(m)}+{\rm pole}_m$ '
      r'from the physical region; per-branch coefficients used in S5(d). Reads '
      r'\texttt{allm\_phys\_result.json}.', ev_se_allm, est=1, deps=('prl_checks/allm_analyze.py',)),
    # -- exchange-even eps^6
    C('even6', 'prl_checks', 'allm_analyze.py', 'Eq. (12), S5', r'Free fit of the transcendental coefficients',
      r'Exchange-even $\epsilon^6$ function for every $m$: a free fit gives $-1/450$ for both the '
      r'$\Pi_m\ln\kappa$ and the $\Pi_m\,\mathrm{Re}\,\psi$ coefficients; residuals of the closed form; for '
      r'$m=0$ a polynomial of degree four. Writes \texttt{allm\_phys\_result.json}.', ev_allm, est=4),
    C('even6', 'prl_checks', 'unified_check.py', 'Eq. (S18), S5', r'Remainder after subtracting Eq.~(S18)',
      r'Physical branch: after subtracting the closed form (S18), the Taylor coefficients $6\le p\le20$ '
      r'vanish for every $m$ (uniqueness (b)).', ev_unified, est=1),
    C('even6', 'prl_checks', 'uniqueness_check.py', 'S5, uniqueness (c), (d)',
      r'Continuation to super-extremal spin; bounds on extra terms',
      r'The $m=2$ polynomial fixed at $|q|\le0.3$ and continued to $1.5\le q\le5000$ reproduces direct MST; '
      r'bounds on additional terms $c_pq^p$; the $q^5$ coefficient $-2.3457997302\times10^{-5}i$ of S5(d).',
      ev_uniqueness, est=1),
    C('even6', 'conservative_continuation', 'phys_hi_p.py', 'Eq. (12), S5', r'Physical-region coefficients to $q^{24}$',
      r'Coefficients up to $q^{24}$ with the $\epsilon^6\ln\epsilon$ term separated, from '
      r'\texttt{../prl\_checks/phys\_l2\_*.json}; writes \texttt{phys\_e6\_hi.json}.', ev_phys_hi_p, est=1),
    C('even6', 'conservative_continuation', 'e6_fit.py', 'Eq. (12)', r'Fit of the $\psi+\log$ term',
      r'The odd coefficients $q^7\ldots q^{23}$ are reproduced by the $\psi+\log$ ansatz with '
      r'$\lambda=-1/450$. Writes \texttt{e6\_fit.json}.', ev_e6_fit, est=2,
      deps=('conservative_continuation/phys_hi_p.py',)),
    C('even6', 'conservative_continuation', 'e6_check_se2.py', 'Eq. (12), S5',
      r'Continued closed form against super-extremal MST',
      r'The closed form continued to $1.5\le q\le5000$ against direct MST at super-extremal spin, with the '
      r'naive log bookkeeping ($10^{-6}$ to $10^{-8}$; the refined test is \texttt{uniqueness\_check.py}).',
      ev_e6_check_se2, est=1, deps=('conservative_continuation/e6_fit.py',)),
    C('even6', 'conservative_continuation', 'e6_identities.py', 'Eq. (16), S9(e)',
      r'$[q^5]\mathcal C=16b_2$ and $\mathcal C/q^5\to64/15$',
      r'$160[q^5]F=16b_2$ of arXiv:2610.40101 to all printed digits, and $160F/q^5\to64/15$ at large $q$.',
      ev_e6_identities, est=1),
    C('even6', 'conservative_continuation', 'se_taylor_l2.py', 'Eq. (17), S9(h)',
      r'Classical $\ell=2$ Taylor coefficients',
      r'Classical ($\epsilon\to0$ at fixed $a\omega$) $\ell=2$ Taylor coefficients for both branches, from '
      r'\texttt{polar\_run2.json}: the exchange-odd part starts at $c^5$ with $-64i/75$, followed by '
      r'$-256i/135$. Writes \texttt{se\_taylor\_l2.json}.', ev_se_taylor_l2, est=1,
      files=('../conservative_test/polar_run2.json',)),
    # -- eps^7
    C('eps7', 'prl_checks', 'd7_phys.py', 'Eq. (S23), S9(g)', r'Exchange-odd $\epsilon^7$ for $m=2$',
      r'Fit of $D_7-(1/\kappa)\partial_{\tau_0}D_6$ to $(u_0+u_1q+u_2q^2+u_3q^3)D_6+{\rm poly}_7(q)$: '
      r'$u_1=10/9$, no $q^2$ and $q^3$ terms.', ev_d7_phys, est=1),
    C('eps7', 'prl_checks', 'd7_fit2.py', 'Eq. (S23), S9(g)', r'Exchange-odd $\epsilon^7$ for all $m$',
      r'$u_0$, $u_1=5m/9$ and $u_2=0$ for $m=2,1,-1,-2,0$.', ev_d7_fit2, est=2),
    C('eps7', 'prl_checks', 'e7_phys.py', 'Eq. (S23), S9(g)', r'Exchange-even $\epsilon^7$ for $m=2$',
      r'The same operator acting on $T_6$, plus a polynomial of degree six; the $\epsilon^7\ln\epsilon$ '
      r'coefficient $B_7$. Writes \texttt{e7\_fit.json} and \texttt{e7\_phys.json}.', ev_e7_phys, est=4),
    C('eps7', 'prl_checks', 'e7_allm.py', 'Eq. (S23), S9(g)', r'Exchange-even $\epsilon^7$ for all $m$',
      r'$u_0$ and $u_1=5m/9$ for $m=1,-1,-2$; for $m=0$, $B_7=-(2i/225)\kappa^4$.', ev_e7_allm, est=5),
    C('eps7', 'prl_checks', 'se7_analyze.py', 'S9(h)', r'Classical $a^6$ coefficients by extrapolation in $q$',
      r'$\epsilon^7$ coefficient at super-extremal $q=1.5\ldots5000$ and the extrapolation of $64A_7/q^6$ to '
      r'$q\to\infty$, against the classical values: the limits commute at $a^6$. Writes '
      r'\texttt{se7\_A7.json}.', ev_se7_analyze, est=1),
    C('eps7', 'prl_checks', 'se7_compare.py', 'S9(h)', r'$\epsilon^7$ closed forms continued to super-extremal $q$',
      r'The physical-region $\epsilon^7$ closed forms continued to $1.5\le q\le5000$ against the direct MST '
      r'coefficients.', ev_se7_compare, est=1, deps=('prl_checks/e7_phys.py', 'prl_checks/se7_analyze.py')),
    C('eps7', 'prl_checks', 'a6_hyp.py', 'S9(h)', r'Negative test at $a^6$',
      r'Shifting only the singular function (trigamma or csch$^2$) gives the wrong large-$q$ limit; the '
      r'check passes when this hypothesis is rejected.', ev_a6_hyp, est=1),
    # -- analytic control
    C('analytic', 'prl_checks', 'analytic_odd.py', 'S7(c)-(e), S9(j)', r'Normalization and the $\epsilon$ coefficient',
      r'$\delta^2\mathcal K_+/[(i/1800)\epsilon^7\kappa^5\tau(1+\tau^2)(4+\tau^2)]\to1+\frac{5mq}{9}\epsilon$, '
      r'and $(e^{i\pi\nu}+e^{-i\pi\nu}\rho)\epsilon/(2i\delta)\to1+(\pi-i\nu_2)\epsilon$.', ev_analytic_odd,
      est=6),
    C('analytic', 'prl_checks', 'series_ratio_check.py', 'S7, S9(j)', r'$\mathcal X$ ratio and $\nu$ to $O(\epsilon^3)$',
      r'Full MST at 60 digits: $\mathcal X$ ratio $1+\frac{5mq}{9}\epsilon+\ldots$ and '
      r'$\nu=2-\frac{107}{210}\epsilon^2+\frac{mq}{12}\epsilon^3$ for four $(m,q)$; $\nu_2$ depends on '
      r'neither $m$ nor $q$.', ev_series_ratio, est=5),
    C('analytic', 'prl_checks', 'eps4_check.py', 'S7(f), S9(j)', r'Pole term absent from $S$; spheroidal coefficients',
      r'The pole term $(i/900)\epsilon^4\Pi_m/\nu_2$ of the near-zone factor is absent from $S$ (its '
      r'$\epsilon^4$ coefficient is cubic in $q$); $E_1=-4m/3$ and $E_2=\frac{10}{189}(9-m^2)-1$.', ev_eps4,
      est=1),
    C('analytic', 'prl_checks', 'h_scaling.py', 'S7, S9(j)', r'Scaling of the near-zone factor',
      r'The exchange-even part of the near-zone factor scales as $\epsilon^4$ (the pole term, cancelled by '
      r'the far zone), the exchange-odd part as $\epsilon^6$; the far-zone factor is exchange-invariant.',
      ev_h_scaling, est=3),
    C('analytic', 'prl_checks', 'absorption_check.py', 'S9(j)', r'$\Gamma_+-\Gamma_-=-4\epsilon^6A_6^-$',
      r'Physical solution against its image, including $q=0$ (Page, $4/225$), after linear extrapolation in '
      r'$\epsilon$.', ev_absorption, est=6),
    # -- general kinematics
    C('kin', '.', 'kinematics.py', 'S6, Table S1', r'BHPT against the contact-term ansatz',
      r'Exchange-odd BHPT amplitude over the ansatz at random kinematic points: ratios $1$ ($a^5$) and $-1$ '
      r'($a^6$, sign convention $a\to-a$) for helicity $+1$; failure for helicity $-1$.', ev_kinematics, est=1),
    C('kin', '.', 'fitcoeffs.py', 'S6, Table S1', r'The nine exchange-odd coefficients by a fit',
      r'A least-squares fit of the ansatz structures to the BHPT amplitude at random points gives the nine '
      r'coefficients of Table~S1; perturbing one of them spoils the agreement.', ev_fitcoeffs, est=1),
    C('kin', 'prl_checks', 'kinematics.py', 'S6, Table S1', r'Copy of the top-level \texttt{kinematics.py}',
      r'Same as the top-level script.', ev_kinematics, est=1),
    C('kin', 'prl_checks', 'kin_even.py', 'S6, Table S1', r'Exchange-even $a^5$ contact term',
      r'The $\psi$-pole contributions are the $\alpha$-tagged terms $c_3=\alpha(64/15,16/3,16/15)$; '
      r'helicity $-1$ is a negative control.', ev_kin_even, est=1),
    C('kin', 'prl_checks', 'kin_tags.py', 'S6, Table S1', r'Tagged and untagged exchange-odd terms',
      r'The $m\neq0$ (coth-pole) part gives the $\alpha$-tagged $c_4=\eta\alpha(64/15,16/5,0)$, the $m=0$ part '
      r'the untagged $c_4=\eta(0,0,4/15)$; the exchange-even $\psi$-pole part gives $c_3$.', ev_kin_tags, est=1),
    C('kin', 'prl_checks', 'kin_a6.py', 'S6, Table S1', r'Tagged and untagged terms at $a^6$',
      r'At $a^6$ the $m\neq0$ (pole) part equals the $\alpha$-tagged terms and the $m=0$ part the untagged '
      r'$c_9^{(2)}=-8/45$; all ratios are $\pm1$.', ev_kin_a6, est=1),
    # -- unit flux
    C('flux', 'flux_test', 'flux_check.py', 'S10, S9(k)', r'Flux balance and factorization',
      r'Floquet amplitudes evaluated separately from the continued-fraction coefficients: MST relations, '
      r'factorized against direct $S_P$, $|\mathcal C|^2=D^2+144M^2\omega^2$, flux balance; all $m$, '
      r'$q=0\ldots0.95$, $\epsilon=0.01\ldots0.4$.', ev_flux_check, est=19),
    C('flux', 'flux_test', 'highfreq_check.py', 'S10', r'Flux balance at $\epsilon=0.6$',
      r'The same flux balance at $M\omega=0.3$, with $\nu$ followed by continuation in $\epsilon$.',
      ev_highfreq, est=18),
    C('flux', 'flux_test', 'complexnu_check.py', 'S10', r'The complex-$\nu$ branch',
      r'Flux balance for $\nu=\frac12+iy$, $\frac32+iy$, $-\frac12+iy$ at $\epsilon=0.8\ldots1.2$, with '
      r'superradiant amplification up to about 11\%.', ev_complexnu, est=28),
    C('flux', 'flux_test', 'decomp_check.py', 'Eq. (11), S10', r'Exchange-parity split of the absorption probability',
      r'$1-|S_+|^2=[1-|S_{\rm even}|^2]+\Delta_{\rm odd}-|S_{\rm odd}|^2$ with the image computed directly: '
      r'$\Delta_{\rm odd}\to\frac{\pi}{225}\kappa^5\varpi(1+\varpi^2)(4+\varpi^2)\bar n_\eta$, one half of '
      r"Page's value at $q=0$, and $\Delta_{\rm odd}/\mathcal P_H=\kappa/(1+\kappa)$ for $m=0$.",
      ev_decomp, est=53),
]
del C


@dataclass
class Production:
    folder: str
    script: str
    args: tuple = ()
    outputs: tuple = ()     # files written in `folder`
    copies: tuple = ()      # (file in `folder`, copy relative to the repository root)
    est: float = 0.0        # seconds on two cores; 0 = not timed

    @property
    def id(self):
        return f'{self.folder}/{self.script}' + ''.join(f' {a}' for a in self.args)


PRODUCTION = [
    Production('conservative_test', 'run_polar2.py', (), ('polar_run2.json',), est=900),
    Production('physical_region_test', 'run_phys.py', (), ('phys_run.json',), est=780),
    Production('physical_region_test', 'run_phys2.py', ('2e-5',), ('phys_l2_2e-5.json',),
               (('phys_l2_2e-5.json', 'prl_checks/phys_l2_2e-5.json'),), est=240),
    Production('physical_region_test', 'run_phys2.py', ('1e-5',), ('phys_l2_1e-5.json',),
               (('phys_l2_1e-5.json', 'prl_checks/phys_l2_1e-5.json'),), est=240),
    Production('prl_checks', 'run_allm.py', (), ('phys_allm.json', 'se_allm.json'), est=1200),
    Production('conservative_continuation', 'e6_se.py', (), ('e6_se.json',),
               (('e6_se.json', 'prl_checks/e6_se.json'),)),
    Production('conservative_continuation', 'e6_se2.py', (), ('e6_se2.json',),
               (('e6_se2.json', 'prl_checks/e6_se2.json'),)),
    Production('conservative_continuation', 'e6_se3.py', (), ('e6_se3.json',),
               (('e6_se3.json', 'prl_checks/e6_se3.json'),)),
    Production('continuation_test', 'd6_se.py', (), ('d6_se.json',)),
]


# ----------------------------------------------------------------------------------------------------
# Running the scripts
# ----------------------------------------------------------------------------------------------------
_PROCS = set()
_PROCS_LOCK = threading.Lock()
_INTERRUPTED = threading.Event()


def _popen(cmd, cwd, **kw):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONIOENCODING='utf-8', PYTHONUNBUFFERED='1')
    if os.name == 'posix':
        kw['start_new_session'] = True          # own process group: a timeout also stops worker processes
    elif hasattr(subprocess, 'CREATE_NEW_PROCESS_GROUP'):
        kw['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP
    proc = subprocess.Popen(cmd, cwd=str(cwd), env=env, text=True, encoding='utf-8', errors='replace', **kw)
    with _PROCS_LOCK:
        _PROCS.add(proc)
    return proc


def _kill(proc):
    try:
        if os.name == 'posix':
            os.killpg(proc.pid, signal.SIGKILL)
        else:
            proc.kill()
    except OSError:
        pass


def kill_all():
    with _PROCS_LOCK:
        procs = list(_PROCS)
    for p in procs:
        _kill(p)


def run_capture(folder, script, args, timeout):
    """Run a script from its folder; returns (stdout, stderr, returncode, timed_out, seconds)."""
    t0 = time.time()
    proc = _popen([sys.executable, script, *args], ROOT / folder, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    timed_out = False
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        _kill(proc)
        out, err = proc.communicate()
    finally:
        with _PROCS_LOCK:
            _PROCS.discard(proc)
    return out or '', err or '', proc.returncode, timed_out, time.time() - t0


def run_streaming(folder, script, args, logfile):
    """Run a production driver, echoing its progress; returns (output, returncode, seconds)."""
    t0 = time.time()
    proc = _popen([sys.executable, script, *args], ROOT / folder, stdout=subprocess.PIPE,
                  stderr=subprocess.STDOUT)
    lines = []
    try:
        with open(logfile, 'w', encoding='utf-8') as f:
            for line in proc.stdout:
                f.write(line)
                f.flush()
                lines.append(line)
                print('    | ' + line.rstrip(), flush=True)
        proc.wait()
    except KeyboardInterrupt:
        _kill(proc)
        proc.wait()
        raise
    finally:
        with _PROCS_LOCK:
            _PROCS.discard(proc)
    return ''.join(lines), proc.returncode, time.time() - t0


@dataclass
class Result:
    id: str
    status: str = 'SKIP'          # PASS, FAIL, ERROR, INFO or SKIP
    reason: str = ''
    time: float = 0.0
    returncode: object = None
    headline: str = ''
    criteria: list = field(default_factory=list)
    log: str = ''                 # log file, relative to the output folder
    output: str = ''


def _log_name(cid):
    return re.sub(r'[^A-Za-z0-9]+', '-', cid).strip('-') + '.txt'


def _missing_modules(names):
    import importlib.util
    return [n for n in names if importlib.util.find_spec(n) is None]


def run_check(c, outdir, timeout):
    res = Result(c.id)
    if _INTERRUPTED.is_set():
        res.reason = 'not run (interrupted)'
        return res
    missing = _missing_modules(c.needs)
    if missing:
        res.reason = f'needs {", ".join(missing)} (not installed)'
        return res
    absent = [f for f in c.files if not (ROOT / c.folder / f).exists()]
    if absent:
        res.reason = f'needs {", ".join(absent)}, produced by the production drivers (--production)'
        return res
    out, err, rc, timed_out, dt = run_capture(c.folder, c.script, c.args, timeout)
    res.time, res.returncode, res.output = dt, rc, out
    log = out
    if err.strip():
        log += ('\n' if out and not out.endswith('\n') else '') + '--- stderr ---\n' + err
    res.log = 'logs/' + _log_name(c.id)
    (outdir / res.log).write_text(log, encoding='utf-8')
    if _INTERRUPTED.is_set() and (timed_out or rc != 0):
        res.status, res.reason = 'SKIP', 'interrupted'
        return res
    if c.evaluate is not None and out.strip():
        e = Ev(out)
        try:
            c.evaluate(e)
        except Exception as ex:                      # a malformed output must not stop the run
            e.add('evaluation of the output', f'{type(ex).__name__}: {ex}', 'no error', False)
        res.criteria = e.crit
        res.headline = e.headline or (e.crit[0].measured if e.crit else '')
    if timed_out:
        res.status, res.reason = 'ERROR', f'timed out after {timeout:.0f} s'
    elif rc != 0:
        last = (err.strip().splitlines() or ['no message'])[-1]
        res.status, res.reason = 'ERROR', f'exit status {rc}: {last}'
    elif c.info:
        res.status, res.reason = 'INFO', 'informational, no pass/fail criterion'
    elif not out.strip():
        res.status, res.reason = 'FAIL', 'the script printed nothing'
    elif all(k.ok for k in res.criteria) and res.criteria:
        res.status = 'PASS'
    else:
        res.status = 'FAIL'
        res.reason = 'failed: ' + '; '.join(k.desc for k in res.criteria if not k.ok)
    return res


def run_checks(checks, jobs, outdir, timeout, on_done):
    """Run the checks, at most `jobs` worker processes at a time, respecting the file dependencies."""
    ids = {c.id for c in checks}
    pending, running, results, load = list(checks), {}, {}, 0
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        try:
            while pending or running:
                started = True
                while started and not _INTERRUPTED.is_set():
                    started = False
                    ready = [c for c in pending if all(d in results or d not in ids for d in c.deps)]
                    for c in sorted(ready, key=lambda c: (-min(c.cpus, jobs), -c.est)):   # longest first
                        need = min(c.cpus, jobs)
                        if load + need <= jobs or load == 0:
                            pending.remove(c)
                            load += need
                            running[ex.submit(run_check, c, outdir, timeout)] = c
                            started = True
                            break
                if not running:
                    break
                done, _ = cf.wait(list(running), timeout=0.5, return_when=cf.FIRST_COMPLETED)
                for f in done:
                    c = running.pop(f)
                    load -= min(c.cpus, jobs)
                    results[c.id] = f.result()
                    on_done(c, results[c.id])
        except KeyboardInterrupt:
            _INTERRUPTED.set()
            print('\ninterrupted: stopping the running scripts ...', flush=True)
            kill_all()
            for f, c in running.items():
                results[c.id] = f.result()
            running.clear()
    for c in pending:
        results[c.id] = Result(c.id, reason='not run (interrupted)')
    return results


def _mtime(path):
    try:
        return path.stat().st_mtime
    except OSError:
        return None


def run_production(steps, outdir):
    """Run the production drivers one after the other (each uses two worker processes)."""
    backup, results = outdir / 'backup', []
    total = sum(p.est for p in steps)
    untimed = sum(1 for p in steps if not p.est)
    print(f'Production: {len(steps)} driver(s), roughly {fmt_time(total)} on two cores'
          + (f' plus {untimed} untimed driver(s)' if untimed else '') +
          f'. Data files that are overwritten are first copied to {backup}.', flush=True)
    for p in steps:
        print(f'\n== production: {p.id}', flush=True)
        targets = [Path(p.folder) / o for o in p.outputs] + [Path(d) for _, d in p.copies]
        for t in targets:
            if (ROOT / t).exists() and not (backup / t).exists():
                (backup / t).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT / t, backup / t)
        before = {o: _mtime(ROOT / p.folder / o) for o in p.outputs}
        log = outdir / 'logs' / ('production-' + _log_name(p.id))
        try:
            _, rc, dt = run_streaming(p.folder, p.script, p.args, log)
        except KeyboardInterrupt:
            _INTERRUPTED.set()
            kill_all()
            rc, dt = None, 0.0
        written = all(_mtime(ROOT / p.folder / o) not in (None, before[o]) for o in p.outputs)
        if rc == 0 and written:
            for src, dst in p.copies:
                shutil.copy2(ROOT / p.folder / src, ROOT / dst)
            status, reason = 'DONE', 'wrote ' + ', '.join(p.outputs) + \
                ''.join(f'; copied to {d}' for _, d in p.copies)
        else:
            for t in targets:                              # put the shipped files back
                if (backup / t).exists():
                    shutil.copy2(backup / t, ROOT / t)
            status = 'ERROR'
            reason = 'interrupted' if _INTERRUPTED.is_set() else \
                (f'exit status {rc}' if rc != 0 else 'the output file was not written')
            reason += '; the shipped data files were restored'
        print(f'== {status} ({fmt_time(dt)}): {reason}', flush=True)
        results.append(dict(id=p.id, status=status, reason=reason, time=dt, log=f'logs/{log.name}'))
        if _INTERRUPTED.is_set():
            break
    return results


# ----------------------------------------------------------------------------------------------------
# Reports
# ----------------------------------------------------------------------------------------------------
def fmt_time(t):
    if not t or t <= 0:
        return '--'
    if t < 1:
        return '<1 s'
    if t < 59.5:
        return f'{t:.0f} s'
    m, s = divmod(int(round(t)), 60)
    if m < 60:
        return f'{m} min {s:02d} s'
    h, m = divmod(m, 60)
    return f'{h} h {m:02d} min'


def git_state():
    if not (ROOT / '.git').exists() or not shutil.which('git'):
        return 'not a git checkout'
    try:
        head = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, capture_output=True, text=True,
                              timeout=20).stdout.strip()
        dirty = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT,
                               capture_output=True, text=True, timeout=20).stdout.strip()
        return (head or 'unknown') + (' (with local modifications)' if dirty else '')
    except (OSError, subprocess.SubprocessError):
        return 'unknown'


def environment(jobs, argv):
    info = {'date': datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
            'python': f'{platform.python_implementation()} {platform.python_version()}',
            'platform': platform.platform(), 'cpus': os.cpu_count(), 'jobs': jobs}
    for mod in ('mpmath', 'sympy', 'numpy'):
        try:
            m = __import__(mod)
            v = getattr(m, '__version__', '?')
            if mod == 'mpmath':
                try:
                    from mpmath.libmp import BACKEND
                    v += f' (backend: {BACKEND})'
                except ImportError:
                    pass
            info[mod] = v
        except ImportError:
            info[mod] = 'not installed'
    info['git'] = git_state()
    info['command'] = ' '.join(['python3', Path(__file__).name] + list(argv))
    return info


def counts(results):
    c = {k: 0 for k in ('PASS', 'FAIL', 'ERROR', 'INFO', 'SKIP')}
    for r in results:
        c[r.status] += 1
    return c


def tex2txt(s):
    """Plain-text version of a short LaTeX string (titles)."""
    s = s.replace(r'\ldots', '...').replace(r'\,', ' ').replace('~', ' ').replace(r'\%', '%')
    s = re.sub(r'\\t?frac\{([^{}]*)\}\{([^{}]*)\}', r'\1/\2', s)
    s = s.replace(r'\epsilon', 'eps').replace(r'\ell', 'l').replace(r'\neq', '!=').replace(r'\le', '<=')
    s = re.sub(r'\\(?:mathrm|texttt|mathcal|rm|quad)\b\s*', '', s)
    s = re.sub(r'\\([A-Za-z]+)', r'\1', s)
    s = s.replace('$', '').replace('{', '').replace('}', '').replace('\\', '')
    return re.sub(r'\s+', ' ', s).strip()


_UNICODE = {'\u2212': '-', '\u2013': '-', '\u2014': '-', '\u2264': '<=', '\u2265': '>=', '\u00d7': 'x',
            '\u2248': '~', '\u2192': '->', '\u221e': 'inf', '\u00b1': '+-', '\u2026': '...'}


def ascii_text(s):
    s = ''.join(_UNICODE.get(ch, ch) for ch in s)
    return s.encode('ascii', 'replace').decode('ascii')


_TEX_SPECIAL = {'\\': r'\textbackslash{}', '{': r'\{', '}': r'\}', '$': r'\$', '&': r'\&', '#': r'\#',
                '_': r'\_', '%': r'\%', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}',
                '<': r'\textless{}', '>': r'\textgreater{}', '|': r'\textbar{}',
                '[': '{[}', ']': '{]}'}         # braces: a cell may not start with [ after \\


def tex_escape(s):
    s = ascii_text(str(s)).replace('<=', '\x00').replace('>=', '\x01')
    s = ''.join(_TEX_SPECIAL.get(ch, ch) for ch in s)
    return s.replace('\x00', r'$\le$ ').replace('\x01', r'$\ge$ ')


def bookmark(title):
    """PDF-bookmark text for a LaTeX title: plain ASCII without TeX special characters."""
    return re.sub(r'[\\{}$&#_%~^<>|]', '', ascii_text(tex2txt(title)))


def _breakable(tt):
    return tt.replace('/', r'/\allowbreak{}').replace(r'\_', r'\_\allowbreak{}')


def tex_path(s):
    """A file path or command in typewriter type, with line breaks allowed after / and _."""
    return r'\texttt{' + _breakable(tex_escape(s)) + '}'


def tex_text(latex):
    """LaTeX text written in this file: allow line breaks inside its \\texttt{...} paths."""
    return re.sub(r'\\texttt\{([^{}]*)\}', lambda m: r'\texttt{' + _breakable(m.group(1)) + '}', latex)


_STATUS_COLOR = {'PASS': 'pass', 'DONE': 'pass', 'FAIL': 'fail', 'ERROR': 'fail', 'INFO': 'info', 'SKIP': 'skip'}


def tex_status(s):
    return rf'\textcolor{{{_STATUS_COLOR[s]}}}{{\textbf{{{s}}}}}'


def log_listing(text, logname, max_lines=320, head=260, tail=40):
    lines = ascii_text(text).replace('\t', '        ').split('\n')
    while lines and not lines[-1].strip():
        lines.pop()
    if not lines:
        return r'{\small\emph{(no output)}}' + '\n'
    if len(lines) > max_lines:
        lines = lines[:head] + [f'[... {len(lines) - head - tail} lines omitted; complete output in '
                                f'{logname} ...]'] + lines[-tail:]
    width = max(len(x) for x in lines)
    size = r'\scriptsize' if width <= 124 else r'\tiny'
    body = '\n'.join(x.replace(r'\end{lstlisting}', r'\end {lstlisting}') for x in lines)
    return f'\\begin{{lstlisting}}[basicstyle=\\ttfamily{size}]\n{body}\n\\end{{lstlisting}}\n'


PREAMBLE = r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=2cm]{geometry}
\IfFileExists{lmodern.sty}{\usepackage[T1]{fontenc}\usepackage{lmodern}}{}
\usepackage{amsmath,amssymb}
\usepackage{array,longtable,booktabs}
\usepackage{xcolor}
\usepackage{listings}
\usepackage{graphicx}
\usepackage[hidelinks,bookmarksnumbered]{hyperref}
\definecolor{pass}{rgb}{0.0,0.50,0.0}
\definecolor{fail}{rgb}{0.80,0.0,0.0}
\definecolor{info}{rgb}{0.0,0.30,0.70}
\definecolor{skip}{rgb}{0.45,0.45,0.45}
\definecolor{lstrule}{gray}{0.75}
\definecolor{lstback}{gray}{0.97}
\definecolor{lstbreak}{gray}{0.55}
\lstset{basicstyle=\ttfamily\scriptsize,columns=fullflexible,keepspaces=true,breaklines=true,
  frame=single,framerule=0.3pt,rulecolor=\color{lstrule},backgroundcolor=\color{lstback},
  xleftmargin=3pt,xrightmargin=3pt,aboveskip=6pt,belowskip=6pt,showstringspaces=false,
  postbreak=\mbox{\textcolor{lstbreak}{$\hookrightarrow$}\space}}
\setlength{\parindent}{0pt}
\setlength{\parskip}{4pt}
\setcounter{tocdepth}{2}
\setlength{\emergencystretch}{3em}
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}
\newcommand{\ok}{\textcolor{pass}{\checkmark}}
\newcommand{\notok}{\textcolor{fail}{$\boldsymbol{\times}$}}
\hypersetup{pdftitle={Numerical checks: run report}}
'''


def write_tex(path, meta, ordered, production, figure):
    res = [r for _, r in ordered]
    n = counts(res)
    bad = n['FAIL'] + n['ERROR']
    W = []
    w = W.append
    w(PREAMBLE)
    w(r'\begin{document}' + '\n')
    w(r'\begin{center}{\LARGE\bfseries Numerical checks: run report}\\[6pt]' + '\n')
    w(r"{\large ``Conservative and absorptive Kerr scattering from a horizon-exchange law''}\\[4pt]" + '\n')
    w(r'{\normalsize Generated by \texttt{run\_all\_checks.py} on ' + tex_escape(meta['date']) + r'}\end{center}' + '\n')
    ran = len(res) - n['SKIP']
    if bad:
        verdict = (rf'\textbf{{{bad} of {ran} checks failed}}' if ran else r'\textbf{No check ran}')
        color = 'fail'
    else:
        verdict = rf'\textbf{{All {n["PASS"]} evaluated checks passed}}'
        color = 'pass'
    detail = (f'{n["PASS"]} passed, {n["FAIL"]} failed, {n["ERROR"]} ended with an error, '
              f'{n["INFO"]} informational, {n["SKIP"]} skipped; wall time {fmt_time(meta["wall"])}.')
    w(r'\begin{center}\fcolorbox{' + color + '}{' + color + r'!7}{\parbox{0.86\textwidth}{\centering\large '
      + verdict + r'\\[3pt]\normalsize ' + tex_escape(detail) + r'}}\end{center}' + '\n')

    w(r'\section*{Run}' + '\n' + r'{\small\begin{tabular}{@{}ll@{}}' + '\n')
    for key, label in (('date', 'Date'), ('command', 'Command'), ('python', 'Python'), ('mpmath', 'mpmath'),
                       ('sympy', 'SymPy'), ('numpy', 'NumPy'), ('platform', 'Platform'), ('cpus', 'CPUs'),
                       ('jobs', 'Scripts in parallel'), ('git', 'Repository state')):
        w(f'{label} & {tex_escape(meta[key])}\\\\\n')
    w(r'Wall time & ' + tex_escape(fmt_time(meta['wall'])) + r'\\' + '\n' + r'\end{tabular}}' + '\n')

    w(r'\section*{How to read this report}' + '\n')
    w(r'Each script is run from its own folder with the Python interpreter above, and its complete output '
      r'is reproduced below its check. The output is parsed, and every number that the paper or the '
      r'READMEs quote is compared with its expected value within the tolerance given in the column '
      r'\emph{Required}; a check passes when all its criteria are met. \textbf{ERROR} means that the script '
      r'did not run to completion, \textbf{INFO} marks exploratory scripts without a pass/fail criterion, '
      r'and \textbf{SKIP} marks checks that need data or modules that are not available. Negative controls '
      r'(the opposite helicity, wrong loop orderings, deliberately wrong hypotheses) pass when the '
      r'comparison fails by a wide margin. Equation and section numbers refer to the paper and its '
      r'Supplemental Material.' + '\n')

    w(r'\section*{Summary}' + '\n')
    w(r'{\small\setlength{\tabcolsep}{3pt}\begin{longtable}{@{}r L{5.3cm} L{2.6cm} L{5.3cm} l r@{}}' + '\n')
    w(r'\toprule \# & Check & Paper & Key result & Status & Time\\ \midrule\endhead' + '\n')
    w(r'\bottomrule\endfoot' + '\n')
    for k, (c, r) in enumerate(ordered, 1):
        key = r.headline if r.status in ('PASS', 'FAIL') and r.headline else (r.reason or r.headline)
        w(rf'\hyperlink{{chk{k}}}{{{k}}} & \hyperlink{{chk{k}}}{{{tex_path(c.id)}}}\newline'
          rf'{{\footnotesize {tex_text(c.title)}}} & {tex_escape(c.paper)} & {tex_escape(key)} & {tex_status(r.status)} & '
          rf'{tex_escape(fmt_time(r.time))}\\' + '\n')
    w(r'\end{longtable}}' + '\n')

    if production:
        w(r'\section*{Production runs}' + '\n')
        w(r'The production drivers regenerated the data files before the checks; the shipped files they '
          r'overwrote are kept in \texttt{backup/}. The last lines of each log follow.' + '\n\n')
        w(r'{\small\begin{longtable}{@{}L{6.2cm} l r L{6.2cm}@{}}\toprule Driver & Status & Time & Result\\'
          r'\midrule\endhead\bottomrule\endfoot' + '\n')
        for p in production:
            w(f'{tex_path(p["id"])} & {tex_status(p["status"])} & {tex_escape(fmt_time(p["time"]))} & '
              f'{tex_escape(p["reason"])}\\\\\n')
        w(r'\end{longtable}}' + '\n')
        for p in production:
            try:
                text = (path.parent / p['log']).read_text(encoding='utf-8', errors='replace')
            except OSError:
                text = ''
            tail_lines = text.rstrip('\n').split('\n')[-25:]
            w(r'\subsection*{' + tex_path(p['id']) + '}\n')
            w(log_listing('\n'.join(tail_lines), p['log']))

    w(r'\clearpage\tableofcontents' + '\n')
    groups = {g.key: g for g in GROUPS}
    index = {c.id: k for k, (c, _) in enumerate(ordered, 1)}
    for g in GROUPS:
        members = [(c, r) for c, r in ordered if c.group == g.key]
        if not members:
            continue
        w('\n' + r'\clearpage\phantomsection\section*{' + g.title + '}\n')
        w(r'\addcontentsline{toc}{section}{\texorpdfstring{' + g.title + '}{' + bookmark(g.title) + '}}\n')
        w(r'{\small Paper: ' + tex_escape(g.paper) + r'}\par' + '\n' + tex_text(g.intro) + '\n')
        for c, r in members:
            k = index[c.id]
            w('\n' + r'\phantomsection\hypertarget{chk' + str(k) + r'}{}\subsection*{' + f'{k}.\\quad ' +
              tex_path(c.id) + '}\n')
            w(r'\addcontentsline{toc}{subsection}{\texorpdfstring{' + f'{k}.\\ ' + tex_path(c.id) + '}{' +
              f'{k}. ' + tex_escape(c.id) + '}}\n')
            w(r'\textbf{' + tex_text(c.title) + r'.} ' + tex_text(c.purpose) + r'\par' + '\n')
            cmd = f'cd {c.folder} && python3 {c.script}' if c.folder != '.' else f'python3 {c.script}'
            cmd += ''.join(f' {a}' for a in c.args)
            w(r'{\small Paper: ' + tex_escape(c.paper) + r' \quad Status: ' + tex_status(r.status) +
              r' \quad Time: ' + tex_escape(fmt_time(r.time)) + r'\\ Command: ' + tex_path(cmd) + r'}\par' + '\n')
            if r.reason and r.status != 'PASS':
                w(r'{\small\emph{' + tex_escape(r.reason) + r'}}\par' + '\n')
            if r.criteria:
                w(r'{\small\begin{longtable}{@{}L{6.9cm}L{4.5cm}L{3.6cm}c@{}}' + '\n')
                w(r'\toprule Criterion & Measured & Required & \\ \midrule\endhead' + '\n')
                w(r'\bottomrule\endfoot' + '\n')
                for k2 in r.criteria:
                    mark = r'\ok' if k2.ok else r'\notok'
                    w(f'{tex_escape(k2.desc)} & {tex_escape(k2.measured)} & {tex_escape(k2.required)} & {mark}\\\\\n')
                w(r'\end{longtable}}' + '\n')
            if r.log:
                w(r'{\small Output (\texttt{' + tex_escape(r.log) + r'}):}' + '\n')
                try:
                    text = (path.parent / r.log).read_text(encoding='utf-8', errors='replace')
                except OSError:
                    text = r.output
                w(log_listing(text, r.log))
            if c.id == 'figure_monodromy/monodromy_loops_check.py' and figure:
                w(r'\begin{center}\includegraphics[width=0.62\textwidth]{' + figure + r'}\\[2pt]' + '\n' +
                  r'{\small Fig.~1 of the paper (\texttt{figure\_monodromy/fig\_monodromy.pdf}); the check above '
                  r'integrates the radial equation along these loops.}\end{center}' + '\n')
    w(r'\end{document}' + '\n')
    path.write_text(''.join(W), encoding='utf-8')


def write_summary(path, meta, ordered, production):
    res = [r for _, r in ordered]
    n = counts(res)
    L = ['Numerical checks of "Conservative and absorptive Kerr scattering from a horizon-exchange law"',
         f'{meta["date"]}   {meta["python"]}, mpmath {meta["mpmath"]}, sympy {meta["sympy"]}, '
         f'numpy {meta["numpy"]}',
         f'{n["PASS"]} PASS, {n["FAIL"]} FAIL, {n["ERROR"]} ERROR, {n["INFO"]} INFO, {n["SKIP"]} SKIP; '
         f'wall time {fmt_time(meta["wall"])}', '']
    for p in production:
        L.append(f'production  {p["status"]:5s} {fmt_time(p["time"]):>11s}  {p["id"]}: {p["reason"]}')
    if production:
        L.append('')
    L.append(f'{"#":>3s}  {"status":6s} {"time":>11s}  {"check":50s} {"paper":28s} key result')
    for k, (c, r) in enumerate(ordered, 1):
        key = r.headline if r.status in ('PASS', 'FAIL') and r.headline else (r.reason or r.headline)
        L.append(f'{k:3d}  {r.status:6s} {fmt_time(r.time):>11s}  {c.id:50s} {c.paper:28s} {ascii_text(key)}')
        if r.status == 'FAIL':
            for k2 in r.criteria:
                if not k2.ok:
                    L.append(f'{"":25s}failed: {k2.desc}: {k2.measured} (required {k2.required})')
    path.write_text('\n'.join(L) + '\n', encoding='utf-8')


def write_json(path, meta, ordered, production):
    data = {'meta': meta, 'production': production, 'checks': []}
    for c, r in ordered:
        data['checks'].append(dict(
            id=c.id, folder=c.folder, script=c.script, args=list(c.args), group=c.group, paper=c.paper,
            title=tex2txt(c.title), status=r.status, reason=r.reason, time=round(r.time, 2),
            returncode=r.returncode, headline=r.headline, log=r.log,
            criteria=[dict(desc=k.desc, measured=k.measured, required=k.required, ok=k.ok) for k in r.criteria]))
    path.write_text(json.dumps(data, indent=1), encoding='utf-8')


def build_pdf(outdir):
    exe = shutil.which('pdflatex')
    if not exe:
        return False, 'pdflatex was not found: compile report.tex with any LaTeX installation'
    for _ in range(2):                                 # twice, for the table of contents and longtables
        try:
            p = subprocess.run([exe, '-interaction=nonstopmode', '-halt-on-error', 'report.tex'], cwd=outdir,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace',
                               timeout=900)
        except subprocess.TimeoutExpired:
            return False, 'pdflatex timed out'
        if p.returncode != 0:
            return False, 'pdflatex failed (see report.log):\n' + '\n'.join(p.stdout.splitlines()[-20:])
    for ext in ('aux', 'out', 'toc', 'log'):
        (outdir / f'report.{ext}').unlink(missing_ok=True)
    return True, str(outdir / 'report.pdf')


# ----------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------
def select(checks, patterns):
    if not patterns:
        return list(checks)
    return [c for c in checks
            if any((any(ch in p for ch in '*?[') and fnmatch.fnmatch(c.id, p)) or p in c.id for p in patterns)]


def parse_args(argv):
    ap = argparse.ArgumentParser(
        description='Run the numerical checks of the paper and write a report (PDF via pdflatex).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='examples:\n  python3 run_all_checks.py\n  python3 run_all_checks.py --quick -j 4\n'
               '  python3 run_all_checks.py --only segment_test closedform_check\n'
               '  python3 run_all_checks.py --production')
    ap.add_argument('-j', '--jobs', type=int, default=None,
                    help='number of scripts run at the same time (default: number of CPUs, at most 8)')
    ap.add_argument('--quick', action='store_true',
                    help='skip the four longest scripts (about 4 min instead of 8 on two cores)')
    ap.add_argument('--only', nargs='+', metavar='PATTERN',
                    help='run only the checks whose path contains PATTERN (shell wildcards allowed)')
    ap.add_argument('--list', action='store_true', help='list the checks and exit')
    ap.add_argument('--timeout', type=float, default=1800, help='time limit per script in seconds (default 1800)')
    ap.add_argument('--outdir', default=None, help='output folder (default: check_report next to this script)')
    ap.add_argument('--no-pdf', action='store_true', help='write report.tex but do not run pdflatex')
    ap.add_argument('--production', nargs='*', metavar='STEP',
                    help='first regenerate the data files with the production drivers (all of them, or those '
                         'whose name contains STEP); the files they overwrite are backed up in OUTDIR/backup')
    return ap.parse_args(argv)


def _raise_interrupt(signum, frame):
    raise KeyboardInterrupt


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)
    # Ctrl-C, termination and hang-up all stop the running scripts (they run in their own process groups)
    signal.signal(signal.SIGINT, signal.default_int_handler)
    for name in ('SIGTERM', 'SIGHUP'):
        if hasattr(signal, name):
            signal.signal(getattr(signal, name), _raise_interrupt)
    jobs = max(1, args.jobs or min(os.cpu_count() or 1, 8))
    checks = select(CHECKS, args.only)

    if args.list:
        print(f'{"#":>3s}  {"check":50s} {"time":>11s}  {"paper":30s} notes')
        for k, c in enumerate(CHECKS, 1):
            notes = ', '.join(x for x in (
                'long (skipped by --quick)' if c.long else '', 'informational' if c.info else '',
                f'{c.cpus} processes' if c.cpus > 1 else '', ('needs ' + ', '.join(c.needs)) if c.needs else '',
                ('needs ' + ', '.join(c.files)) if c.files else '') if x)
            print(f'{k:3d}  {c.id:50s} {fmt_time(c.est):>11s}  {c.paper:30s} {notes}')
        print('\nproduction drivers (--production):')
        for p in PRODUCTION:
            print(f'     {p.id:50s} {fmt_time(p.est) if p.est else "not timed":>11s}  writes {", ".join(p.outputs)}')
        return 0
    if not checks:
        print('no check matches', ' '.join(args.only or []))
        return 2

    outdir = Path(args.outdir).resolve() if args.outdir else ROOT / 'check_report'
    (outdir / 'logs').mkdir(parents=True, exist_ok=True)
    for f in (outdir / 'logs').glob('*.txt'):
        f.unlink()
    for name in ('report.pdf', 'report.tex', 'summary.txt', 'results.json'):
        (outdir / name).unlink(missing_ok=True)
    t0 = time.time()
    meta = environment(jobs, argv)

    production = []
    if args.production is not None:
        steps = [p for p in PRODUCTION if not args.production or any(s in p.id for s in args.production)]
        production = run_production(steps, outdir)

    to_run = [c for c in checks if not (args.quick and c.long)]
    skipped = {c.id: Result(c.id, reason='skipped (--quick)') for c in checks if args.quick and c.long}
    est = sum(c.est * min(c.cpus, jobs) for c in to_run) / jobs
    print(f'\nRunning {len(to_run)} checks, up to {jobs} at a time (roughly {fmt_time(1.15 * est)}); '
          f'output in {outdir}', flush=True)
    color = sys.stdout.isatty() and os.name == 'posix'
    paint = {'PASS': '\033[32m', 'FAIL': '\033[31m', 'ERROR': '\033[31m', 'INFO': '\033[34m', 'SKIP': '\033[90m'}
    done = [0]

    def on_done(c, r):
        done[0] += 1
        st = f'{paint[r.status]}{r.status:5s}\033[0m' if color else f'{r.status:5s}'
        msg = r.headline if r.status in ('PASS', 'FAIL') and r.headline else (r.reason or r.headline)
        print(f'[{done[0]:2d}/{len(to_run)}] {st} {fmt_time(r.time):>11s}  {c.id:48s} {ascii_text(msg)}', flush=True)

    results = {}
    if not _INTERRUPTED.is_set():
        results = run_checks(to_run, jobs, outdir, args.timeout, on_done)
    else:
        results = {c.id: Result(c.id, reason='not run (interrupted)') for c in to_run}
    results.update(skipped)
    ordered = [(c, results.get(c.id, Result(c.id, reason='not run'))) for c in checks]
    meta['wall'] = time.time() - t0

    figure = None
    fig = ROOT / 'figure_monodromy' / 'fig_monodromy.pdf'
    if fig.exists() and any(c.folder == 'figure_monodromy' for c in checks):
        shutil.copy2(fig, outdir / 'fig_monodromy.pdf')
        figure = 'fig_monodromy.pdf'
    write_json(outdir / 'results.json', meta, ordered, production)
    write_summary(outdir / 'summary.txt', meta, ordered, production)
    write_tex(outdir / 'report.tex', meta, ordered, production, figure)
    n = counts([r for _, r in ordered])
    print(f'\n{n["PASS"]} passed, {n["FAIL"]} failed, {n["ERROR"]} errors, {n["INFO"]} informational, '
          f'{n["SKIP"]} skipped; wall time {fmt_time(meta["wall"])}')
    for c, r in ordered:
        if r.status in ('FAIL', 'ERROR'):
            print(f'  {r.status}: {c.id}: {ascii_text(r.reason)}')
    if not args.no_pdf:
        ok, msg = build_pdf(outdir)
        print(('Report: ' if ok else '') + msg)
    print(f'Summary: {outdir / "summary.txt"}')
    if _INTERRUPTED.is_set():
        return 130
    return 1 if n['FAIL'] + n['ERROR'] else 0


if __name__ == '__main__':
    sys.exit(main())
