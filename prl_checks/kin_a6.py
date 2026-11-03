"""a^6: do BGKV's alpha-tagged terms coincide with the single-branch pole contributions (m != 0 sector) in general kinematics?
eta-odd : BHPT weights w[m] (10m/9 diagonal + l=2-3 spheroidal mixing) split into m!=0 and m=0 parts.
eta-even: pole part = -sgn(m) x eta-odd weights (same structure)."""
import io, contextlib, mpmath as mp, random
with contextlib.redirect_stdout(io.StringIO()):
    import kinematics as K
mp.mp.dps = 30
I = mp.mpc(0, 1); f = mp.factorial; s = -2
def bhpt6(gam, th, ph, w):
    Th = mp.acos(mp.sin(gam)*mp.sin(th)*mp.cos(ph) + mp.cos(gam)*mp.cos(th))
    D0 = sum(K.sY(s,2,m,gam)*K.sY(s,2,m,th)*mp.exp(I*m*ph) for m in range(-2,3))
    N = sum(w[m]*(mp.mpf(10*m)/9*K.sY(s,2,m,gam)*K.sY(s,2,m,th)
                  + K.b[m]*(K.sY(s,3,m,gam)*K.sY(s,2,m,th) + K.sY(s,2,m,gam)*K.sY(s,3,m,th)))*mp.exp(I*m*ph) for m in range(-2,3))
    return -mp.mpf(8)/45*mp.sin(Th/2)**2*N/D0
def struct(gam, th, ph, hel=1):
    nin = [mp.sin(gam), 0, mp.cos(gam)]; nout = [mp.sin(th)*mp.cos(ph), mp.sin(th)*mp.sin(ph), mp.cos(th)]
    cr = [nin[1]*nout[2]-nin[2]*nout[1], nin[2]*nout[0]-nin[0]*nout[2], nin[0]*nout[1]-nin[1]*nout[0]]
    cosT = sum(nin[i]*nout[i] for i in range(3))
    wa = (nin[2] + nout[2] + hel*I*cr[2])/(1 + cosT); k2a, mk3a = nin[2], nout[2]
    xi = -2/(1 - cosT); A, B = wa - k2a, wa - mk3a
    pref = [xi**(j-1)*wa**(4-2*j)*A**j*B**j for j in range(3)]
    return dict(pref=pref, wa=wa, k2a=k2a, mk3a=mk3a, xi=xi, A=A, B=B)
def ans_odd(S, c9, c10):
    return sum(S['pref'][j]*(c9[j]*(S['k2a'] + S['mk3a']) + c10[j]*S['wa']) for j in range(3))
def ans_even(S, c6, c8, d1, f1):
    X = 2*S['wa'] - S['k2a'] - S['mk3a']
    return (sum(S['pref'][j]*(c6[j]*X*S['wa'] + c8[j]*S['wa']**2) for j in range(3))
            + d1*S['wa']**6/S['xi']**2 + f1*S['xi']**2*S['A']**3*S['B']**3)
F = mp.mpf
wodd = {m: (-1)**m*f(2+m)*f(2-m) for m in range(-2, 3)}
w_pole = {m: (wodd[m] if m else 0) for m in range(-2, 3)}
w_zero = {m: (wodd[m] if m == 0 else 0) for m in range(-2, 3)}
w_even = {m: -mp.sign(m)*wodd[m] for m in range(-2, 3)}
# BGKV a^6 (alpha=1, eta factored out): tagged / untagged pieces
c9_tag, c9_un = [-F(128)/45, -F(32)/15, 0], [0, 0, -F(8)/45]
c10_tag, c10_un = [-F(256)/45, -F(352)/45, -F(32)/15], [0, 0, 0]
c6_tag, c8_tag, f1_tag = [F(128)/45, F(32)/9, F(32)/45], [-F(512)/45, -F(160)/9, -F(304)/45], -F(16)/45
random.seed(3)
pts = [(F(random.uniform(0.2, 2.9)), F(random.uniform(0.2, 2.9)), F(random.uniform(0, 6.2))) for _ in range(6)]
print("ratios BHPT / BGKV-piece at 6 random kinematic points (helicity +1):")
for lab, w, fn in (("eta-odd full      / (tagged+untagged)", wodd, lambda S: ans_odd(S, [a+b for a, b in zip(c9_tag, c9_un)], [a+b for a, b in zip(c10_tag, c10_un)])),
                   ("eta-odd m!=0 part / alpha-tagged     ", w_pole, lambda S: ans_odd(S, c9_tag, c10_tag)),
                   ("eta-odd m=0 part  / untagged         ", w_zero, lambda S: ans_odd(S, c9_un, c10_un)),
                   ("eta-even pole part/ alpha-tagged     ", w_even, lambda S: ans_even(S, c6_tag, c8_tag, 0, f1_tag))):
    r = [bhpt6(g, t, p, w)/fn(struct(g, t, p)) for (g, t, p) in pts]
    print(f"  {lab}: {[mp.nstr(x, 8) for x in r]}")
