"""Segment formula vs the 110-digit production MST data (conservative_test/polar_run2.json)."""
import json, sys
import mpmath as mp
mp.mp.dps = 40
import segment_born as sb
d = json.load(open('../conservative_test/polar_run2.json'))
want = [(l, k) for l in (2, 3, 4, 5, 6) for k in (0, 5)]
for r in sorted(d, key=lambda r: (r['l'], r['k'])):
    if (r['l'], r['k']) not in want: continue
    c = mp.mpf(r['c'])
    s1p = mp.mpc(*map(mp.mpf, r['s1_1'])); s1m = mp.mpc(*map(mp.mpf, r['s1_-1']))
    mst = (s1p - s1m)/2
    seg = sb.S_odd(r['l'], c, verbose=False)
    print(f"l={r['l']} c={mp.nstr(c,8)}: MST {mp.nstr(mst,13)}   segment {mp.nstr(seg,13)}   |rel diff| {mp.nstr(abs(seg/mst-1),3)}")
    sys.stdout.flush()
