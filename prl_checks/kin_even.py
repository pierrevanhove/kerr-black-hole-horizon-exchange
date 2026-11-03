"""General-kinematics test of the eta-even (conservative) a^5 contact term.
BHPT side: only l=2 carries a near-zone a^5 contact; per-m coefficient = pole_m (the psi pole picked up by the
continued self-dual function), i.e. w_even[m] = -sgn(m) w_odd[m], w_odd[m] = (-1)^m (2+m)!(2-m)!.
Ansatz side (BGKV): P5_even = sum_j xi^(j-1) (w.a)^(4-2j) ((w-k2).a)^j ((w+k3).a)^j c3^(j) (w.a)."""
import io, contextlib, mpmath as mp, random
with contextlib.redirect_stdout(io.StringIO()):
    import kinematics as K
mp.mp.dps = 30
I = mp.mpc(0, 1); f = mp.factorial
wodd = {m: (-1)**m*f(2+m)*f(2-m) for m in range(-2, 3)}
weven = {m: -mp.sign(m)*wodd[m] for m in range(-2, 3)}
def bhpt(gam, th, ph, w):
    Th = mp.acos(mp.sin(gam)*mp.sin(th)*mp.cos(ph) + mp.cos(gam)*mp.cos(th))
    D0 = sum(K.sY(-2, 2, m, gam)*K.sY(-2, 2, m, th)*mp.exp(I*m*ph) for m in range(-2, 3))
    N = sum(w[m]*K.sY(-2, 2, m, gam)*K.sY(-2, 2, m, th)*mp.exp(I*m*ph) for m in range(-2, 3))
    return -mp.mpf(8)/45*mp.sin(Th/2)**2*N/D0
def structures(gam, th, ph, hel):
    nin = [mp.sin(gam), 0, mp.cos(gam)]
    nout = [mp.sin(th)*mp.cos(ph), mp.sin(th)*mp.sin(ph), mp.cos(th)]
    cr = [nin[1]*nout[2]-nin[2]*nout[1], nin[2]*nout[0]-nin[0]*nout[2], nin[0]*nout[1]-nin[1]*nout[0]]
    cosT = sum(nin[i]*nout[i] for i in range(3))
    wa = (nin[2] + nout[2] + hel*I*cr[2])/(1 + cosT)
    k2a, mk3a = nin[2], nout[2]
    xi = -2/(1 - cosT)
    A, B = wa - k2a, wa - mk3a
    return [xi**(j-1)*wa**(4-2*j)*A**j*B**j for j in range(3)], wa
random.seed(11)
pts = [(mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0, 6.2))) for _ in range(8)]
c4 = [mp.mpf(64)/15, mp.mpf(16)/5, mp.mpf(4)/15]
for hel in (1, -1):
    print(f"helicity {hel:+d}")
    rows, rhs = [], []
    for (g, t, p) in pts:
        pref, wa = structures(g, t, p, hel)
        odd_ratio = bhpt(g, t, p, wodd)/sum(pref[j]*c4[j] for j in range(3))
        e = bhpt(g, t, p, weven)
        rows.append([pref[j]*wa for j in range(3)]); rhs.append(e)
        print(f"   eta-odd BHPT/ansatz(c4, alpha=1) = {mp.nstr(odd_ratio, 12)}")
    # least-squares fit of c3^(j) to the BHPT eta-even contact (complex system, 8 points, 3 unknowns), with the odd ratio
    Mr = mp.matrix(rows); b = mp.matrix(rhs)
    sol, res = mp.qr_solve(Mr, b)
    print(f"   fitted c3 (in units of the eta-odd normalization): {[mp.nstr(x/odd_ratio, 14) for x in sol]}   residual {mp.nstr(res, 3)}")
    print(f"   BGKV alpha=1: [64/15, 16/3, 4/3] = {[mp.nstr(x, 14) for x in (mp.mpf(64)/15, mp.mpf(16)/3, mp.mpf(4)/3)]};  alpha=0: [0, 0, 4/15]")
