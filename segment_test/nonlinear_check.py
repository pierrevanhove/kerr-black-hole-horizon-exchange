"""Size of the terms nonlinear in K in the classical limit (eps -> 0 at fixed c = a omega), l = 2, s = -2.
S/F - 1 = H = i E K/(1 + beta K),  beta = -i e^{-i pi nu} rho,  E = e^{i pi nu} + e^{-i pi nu} rho.
The eta-odd part is (i/2) E (K+ - K-) / ((1 + beta K+)(1 + beta K-))."""
import sys
import mpmath as mp
mp.mp.dps = 60
from csm import S_phys, spheroidal
I = mp.mpc(0, 1)
l, s = 2, -2
eps = mp.mpf('1e-12')
for m in (2,):
    for c in ('0.025', '0.05', '0.1', '0.2'):
        c = mp.mpf(c)
        E0, _ = spheroidal(s, l, m, c)
        res = {}
        for eta in (1, -1):
            S, x, A = S_phys(l, m, s, c, eta, eps, mp.mpf('-0.5'), E=E0)
            nu, rho, X, H, Sfar = A['nu'], A['rho'], A['X'], A['H'], A['Sfar']
            beta = -I*mp.exp(-I*mp.pi*nu)*rho
            E = mp.exp(I*mp.pi*nu) + mp.exp(-I*mp.pi*nu)*rho
            res[eta] = dict(S=S, x=x, nu=nu, X=X, H=H, Sfar=Sfar, beta=beta, E=E)
        p, q_ = res[1], res[-1]
        Hodd = (p['H'] - q_['H'])/2
        lin = I/2*p['E']*(p['X'] - q_['X'])
        Sodd = (p['S'] - q_['S'])/2
        print(f"m={m} c={mp.nstr(c,3)}: x={mp.nstr(p['x'],8)}  |Sfar+/Sfar- - 1|={mp.nstr(abs(p['Sfar']/q_['Sfar']-1),3)}")
        print(f"   beta K(+1) = {mp.nstr(p['beta']*p['X'],8)}   beta K(-1) = {mp.nstr(q_['beta']*q_['X'],8)}")
        print(f"   beta K(+1)/c^3 = {mp.nstr(p['beta']*p['X']/c**3,8)}   (K odd part)/c^5 = {mp.nstr((p['X']-q_['X'])/2/c**5,8)}")
        print(f"   H_odd exact/linear - 1 = {mp.nstr(Hodd/lin - 1, 8)}   ( /c^3 = {mp.nstr((Hodd/lin-1)/c**3, 8)} )")
        print(f"   S_odd/eps = {mp.nstr(Sodd/eps, 12)}")
        sys.stdout.flush()
