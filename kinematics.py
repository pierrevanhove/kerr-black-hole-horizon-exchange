"""
General-kinematics check of the eta-odd Compton coefficients (BGKV Table 1 / Bautista Table V).
BHPT side (from the swap law, s=-2): only l=2 contributes to the eta-odd part through O(a^6):
   H_odd(2,m) = -(8/225) i (-1)^m (2+m)!(2-m)! (aw)^5 eps [1 + (10m/9) aw]
   f_odd = (2pi/(i w)) sum_m H_odd(2,m) S_2m(gamma;c) S_2m(theta;c) e^{i m phi},   c = a w
   f_0   = M (4pi/5) D0 / sin^2(Theta/2),  D0 = sum_m Y_2m(gamma) Y_2m(theta) e^{i m phi}
Compton-ansatz side: P_odd = sum_m xi^{m-1} (w.a)^{4-2m} ((w-k2).a)^m ((w+k3).a)^m r_m  (eta-odd pieces of r_m).
"""
import mpmath as mp, itertools, random
mp.mp.dps = 30
I = mp.mpc(0, 1)
f = mp.factorial

def wigner_d(j, mp_, m, beta):
    tot = 0
    for k in range(0, 2*j+1):
        a1, a2, a3, a4 = j+m-k, k, j-k-mp_, k-m+mp_
        if min(a1, a2, a3, a4) < 0: continue
        tot += ((-1)**(k-m+mp_) * mp.sqrt(f(j+m)*f(j-m)*f(j+mp_)*f(j-mp_)) / (f(a1)*f(a2)*f(a3)*f(a4))
                * mp.cos(beta/2)**(2*j-2*k+m-mp_) * mp.sin(beta/2)**(2*k-m+mp_))
    return tot

def sY(s, l, m, th):          # spin-weighted harmonic without e^{i m phi}; convention-independent where used
    return (-1)**m * mp.sqrt((2*l+1)/(4*mp.pi)) * wigner_d(l, m, -s, th)

s = -2
# <3m|x|2m> by quadrature, and first-order spheroidal mixing b_m (S_2m = Y_2m + c b_m Y_3m)
def xel(m):
    return 2*mp.pi*mp.quad(lambda t: sY(s,3,m,t)*mp.cos(t)*sY(s,2,m,t)*mp.sin(t), [0, mp.pi])
b = {m: 2*s*xel(m)/(6-12) for m in range(-2, 3)}

def bhpt(gam, th, ph, order):
    """returns coefficient of (aw)^order in f_odd/f_0"""
    Th = mp.acos(mp.sin(gam)*mp.sin(th)*mp.cos(ph) + mp.cos(gam)*mp.cos(th))
    D0 = sum(sY(s,2,m,gam)*sY(s,2,m,th)*mp.exp(I*m*ph) for m in range(-2,3))
    w = {m: (-1)**m*f(2+m)*f(2-m) for m in range(-2,3)}
    if order == 5:
        N = sum(w[m]*sY(s,2,m,gam)*sY(s,2,m,th)*mp.exp(I*m*ph) for m in range(-2,3))
    else:
        N = sum(w[m]*(mp.mpf(10*m)/9*sY(s,2,m,gam)*sY(s,2,m,th)
                      + b[m]*(sY(s,3,m,gam)*sY(s,2,m,th) + sY(s,2,m,gam)*sY(s,3,m,th)))*mp.exp(I*m*ph) for m in range(-2,3))
    return -mp.mpf(8)/45*mp.sin(Th/2)**2*N/D0, Th, D0

def ansatz(gam, th, ph, order, hel, c):
    nin = mp.matrix([mp.sin(gam), 0, mp.cos(gam)])
    nout = mp.matrix([mp.sin(th)*mp.cos(ph), mp.sin(th)*mp.sin(ph), mp.cos(th)])
    cr = mp.matrix([nin[1]*nout[2]-nin[2]*nout[1], nin[2]*nout[0]-nin[0]*nout[2], nin[0]*nout[1]-nin[1]*nout[0]])
    cosT = sum(nin[i]*nout[i] for i in range(3))
    v = [(nin[i]+nout[i] + hel*I*cr[i])/(1+cosT) for i in range(3)]
    az = mp.matrix([0, 0, 1])
    k2a, mk3a, wa = nin[2], nout[2], v[2]            # k2.a = n_in.a, -k3.a = n_out.a, w.a = v.a   (a = z-hat, w = 1)
    xi = -2/(1-cosT)
    A, B = wa - k2a, wa - mk3a                      # (w-k2).a , (w+k3).a
    pref = [xi**(mm-1)*wa**(4-2*mm)*A**mm*B**mm for mm in range(3)]
    if order == 5:
        return sum(pref[mm]*c['c4'][mm] for mm in range(3))
    return sum(pref[mm]*(c['c9'][mm]*(k2a+mk3a) + c['c10'][mm]*wa) for mm in range(3))

# Table 1 of 2212.07965 / Table V of 2304.04287 with alpha = 1, eta factored out
table = {'c4': [mp.mpf(64)/15, mp.mpf(16)/5, mp.mpf(4)/15],
         'c9': [-mp.mpf(128)/45, -mp.mpf(32)/15, -mp.mpf(8)/45],
         'c10': [-mp.mpf(256)/45, -mp.mpf(352)/45, -mp.mpf(32)/15]}

random.seed(7)
pts = [(mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0.2, 2.9)), mp.mpf(random.uniform(0, 6.2))) for _ in range(6)]

# sanity: addition theorem |D0| = (5/4pi) cos^4(Theta/2)  (helicity-preserving combination)
g, t, p = pts[0]; _, Th, D0 = bhpt(g, t, p, 5)
print("addition theorem check:", mp.nstr(abs(D0)/(5/(4*mp.pi)*mp.cos(Th/2)**4), 15))

for order, keys in ((5, ['c4']), (6, ['c9', 'c10'])):
    print(f"\n=== O(a^{order}) ===")
    for hel in (+1, -1):
        ratios = []
        for (g, t, p) in pts:
            lhs = bhpt(g, t, p, order)[0]
            rhs = ansatz(g, t, p, order, hel, table)
            ratios.append(lhs/rhs)
        print(f"helicity {hel:+d}: BHPT/ansatz at 6 random kinematic points:", [mp.nstr(r, 10) for r in ratios])
