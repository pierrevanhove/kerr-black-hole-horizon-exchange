import mpmath as mp
mp.mp.dps = 60
from polar_m import H_odd
s, l = -2, 2
for m in (1, 0, -1, -2):
    rows = []
    for aw in ('0.004', '0.008', '0.012', '0.016', '0.02'):
        e1, e2 = mp.mpf('2e-5'), mp.mpf('1e-5')
        h0 = 2*H_odd(s, l, m, aw, e2)/e2 - H_odd(s, l, m, aw, e1)/e1
        rows.append((mp.mpf(aw), h0/mp.mpf(aw)**5))
    A = mp.matrix([[x**k for k in range(5)] for x, _ in rows]); b = mp.matrix([r for _, r in rows])
    C = mp.lu_solve(A, b)
    print(f"m={m:+d}: C5={mp.nstr(C[0],10)}  C6/C5={mp.nstr(C[1]/C[0],8)}  (10m/9={mp.nstr(mp.mpf(10*m)/9,8)})", flush=True)
