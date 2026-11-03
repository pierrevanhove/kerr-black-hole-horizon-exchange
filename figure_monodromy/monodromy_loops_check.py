"""Numerical check of the loop algebra drawn in Fig. 1 for the actual radial Teukolsky equation (s=-2, l=2):
transport matrices T_gamma of (R, R') along explicit loops based at r0 > r+, all counterclockwise, so that the
monodromy of 'gamma then gamma'' is T_gamma' T_gamma.  Checks
  (1) T_inf = T_+ T_-^(above) = T_-^(below) T_+          (the relation M_inf = M_+ M_-^(+) = M_-^(-) M_+),
  (2) eigenvalues: T_+ -> exp(+-2 pi eps+), T_-^(.) -> exp(+-2 pi eps-), T_inf -> exp(+-2 i pi nu) with nu from MST,
  (3) the image defined by the 'above' lasso has K_-/K_+ equal to the MST image with log(-kappa) = log kappa + i pi,
      i.e. eta = +1 is the class passing above r+ (panel b)."""
import numpy as np, mpmath as mp
from mst import MST
from phys import spheroidal_c
from csm import secant
mp.mp.dps = 30

s, l, m, q, eps = -2, 2, 2, 0.6, 0.3
om = eps/2; kap = np.sqrt(1 - q**2); rp, rm = 1 + kap, 1 - kap
c = q*eps/2
E, _ = spheroidal_c(s, l, m, mp.mpf(c))
lam = complex(E - 2*m*c + c**2 - s*(s + 1))

def rhs(r, Y):
    D = r*r - 2*r + q*q; Dp = 2*r - 2
    K = (r*r + q*q)*om - q*m
    V = (K*K - 2j*s*(r - 1)*K)/D + 4j*s*om*r - lam
    return np.array([Y[1], -((s + 1)*Dp*Y[1] + V*Y[0])/D])

def transport(path, n=20000):
    """path: list of callables r(t), t in [0,1] (with derivative by finite difference); RK4 on the 2x2 matrix."""
    T = np.eye(2, dtype=complex)
    for seg in path:
        h = 1.0/n
        for k in range(n):
            t = k*h
            def f(tt, Y):
                r = seg(tt); dr = (seg(tt + 1e-7) - seg(tt - 1e-7))/2e-7
                return dr*np.column_stack([rhs(r, Y[:, j]) for j in range(2)])
            k1 = f(t, T); k2 = f(t + h/2, T + h/2*k1); k3 = f(t + h/2, T + h/2*k2); k4 = f(t + h, T + h*k3)
            T = T + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return T

line = lambda a, b: (lambda t: a + (b - a)*t)
circ = lambda c0, rad, th0: (lambda t: c0 + rad*np.exp(1j*(th0 + 2*np.pi*t)))
r0, rho, H = 3.0, 0.4, 1.2
def lasso(tail, center, th0):
    back = [line(p.__call__(1), p.__call__(0)) for p in tail[::-1]]
    return tail + [circ(center, rho, th0)] + back
loop_p = lasso([line(r0, rp + rho)], rp, 0.0)
up = [line(r0, r0 + 1j*H), line(r0 + 1j*H, rm + 1j*H), line(rm + 1j*H, rm + 1j*rho)]
dn = [line(r0, r0 - 1j*H), line(r0 - 1j*H, rm - 1j*H), line(rm - 1j*H, rm - 1j*rho)]
loop_m_above = lasso(up, rm, np.pi/2)
loop_m_below = lasso(dn, rm, -np.pi/2)
loop_inf = [circ(1.0, r0 - 1.0, 0.0)]

Tp, Ta, Tb, Ti = (transport(L) for L in (loop_p, loop_m_above, loop_m_below, loop_inf))
print("(1) |T_inf - T_+ T_-^above| =", np.abs(Ti - Tp @ Ta).max(), "  |T_inf - T_-^below T_+| =", np.abs(Ti - Tb @ Tp).max(),
      "  (wrong orders: ", np.abs(Ti - Ta @ Tp).max(), np.abs(Ti - Tp @ Tb).max(), ")")
tau = (eps - m*q)/kap; ep, em = (eps + tau)/2, (eps - tau)/2
P = MST(eps, q, m, s, lam, kap, N=40)
x = secant(lambda x: P.g(l + x*eps**2)/eps**2, -0.5, mp.mpf(10)**-20); nu = complex(l + x*eps**2)
ev = lambda T: sorted(np.linalg.eigvals(T), key=lambda z: abs(z))
print("(2) eig T_+ :", np.round(ev(Tp), 8), " expected", np.round(sorted([np.exp(2*np.pi*ep), np.exp(-2*np.pi*ep)], key=abs), 8))
print("    eig T_- :", np.round(ev(Ta), 8), np.round(ev(Tb), 8), " expected", np.round(sorted([np.exp(2*np.pi*em), np.exp(-2*np.pi*em)], key=abs), 8))
lamb = np.exp(2j*np.pi*nu)
print("    eig T_inf:", np.round(ev(Ti), 8), " expected", np.round(sorted([lamb, 1/lamb], key=abs), 8), " (nu =", np.round(nu, 10), ")")
# (3) responses in the Floquet basis: eigenvectors of T_inf for lambda and 1/lambda
w, V = np.linalg.eig(Ti); iL = np.argmin(abs(w - lamb)); iI = np.argmin(abs(w - 1/lamb))
Fb = np.column_stack([V[:, iL], V[:, iI]])
def K(T, mu):
    w, V = np.linalg.eig(T); v = V[:, np.argmin(abs(w - mu))]
    cf = np.linalg.solve(Fb, v); return cf[1]/cf[0]
Kp = K(Tp, np.exp(2*np.pi*ep)); Kabove = K(Ta, np.exp(2*np.pi*em)); Kbelow = K(Tb, np.exp(2*np.pi*em))
A = P.amplitudes() if setattr(P, 'nu', mp.mpc(nu)) is None else None
Pm = {}
for eta in (1, -1):
    Q = MST(eps, q, m, s, lam, -kap, logkappa=mp.log(kap) + 1j*mp.pi*eta, N=40); Q.nu = mp.mpc(nu)
    Pm[eta] = complex(Q.amplitudes()['X']/A['X'])
print("(3) K_-/K_+ from monodromy: above r+ ->", np.round(Kabove/Kp, 8), "  below r+ ->", np.round(Kbelow/Kp, 8))
print("    MST image, eta=+1        ->", np.round(Pm[1], 8), "  eta=-1 ->", np.round(Pm[-1], 8))
