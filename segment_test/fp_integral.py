"""Robust finite-part integral of F along r = i a', a' in (-a, a): Taylor expansion of the regular part
near the double pole at -ia (Cauchy coefficients), ordinary quadrature away from it."""
import mpmath as mp
I = mp.mpc(0, 1)

def fp_segment(F, a, rho=mp.mpf('0.1'), Rc=mp.mpf('0.6'), K=70, nC=256):
    rm = -I*a
    # Cauchy coefficients of h(y) = y^2 F(rm + y) on |y| = Rc (trapezoid rule is spectrally accurate)
    pts = [Rc*mp.expj(2*mp.pi*j/nC) for j in range(nC)]
    vals = [y**2*F(rm + y) for y in pts]
    h = [mp.fsum(v*mp.expj(-2*mp.pi*j*k/nC) for j, v in enumerate(vals))/nC/Rc**k for k in range(K+2)]
    C2, C1 = h[0], h[1]
    # regular part g(y) = sum_{k>=0} h[k+2] y^k ; integrate from y = 0 to y = i rho
    yr = I*rho
    near = mp.fsum(h[k+2]*yr**(k+1)/(k+1) for k in range(K))
    # far part: y from i rho to 2ia, i.e. a' from -a + rho to a
    far = mp.quad(lambda ap: F(I*ap)*I, [-a + rho, 0, a])
    # FP of the singular part from 0 to i rho: C2*FP int dy/y^2 = -C2/(i rho); C1 log term
    sing = -C2/yr
    return near + far + sing, C2, C1
