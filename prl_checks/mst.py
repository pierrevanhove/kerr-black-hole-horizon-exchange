"""
Minimal MST (Mano-Suzuki-Takasugi) solver for the Teukolsky radial equation,
following Sasaki & Tagoshi, Living Rev. Rel. 6 (2003) 6, Sec. 4.
Units: M = 1, so eps = 2 M omega = 2 omega, q = a/M.
The horizon branch enters only through kappa (r_pm = 1 pm kappa) and log(kappa);
r+ <-> r- is kappa -> -kappa (with an explicit choice of log branch).
"""
import mpmath as mp

I = mp.mpc(0, 1)


class MST:
    def __init__(self, eps, q, m, s, lam, kappa, logkappa=None, N=50):
        self.eps = mp.mpmathify(eps)
        self.q = mp.mpmathify(q)
        self.m = m
        self.s = s
        self.lam = mp.mpmathify(lam)
        self.kappa = mp.mpmathify(kappa)
        self.logk = mp.log(self.kappa) if logkappa is None else mp.mpmathify(logkappa)
        self.tau = (self.eps - m * self.q) / self.kappa
        self.epp = (self.eps + self.tau) / 2
        self.N = N

    # ---- three-term recurrence -------------------------------------------
    def al(self, n, nu):
        e, k, t, s = self.eps, self.kappa, self.tau, self.s
        return I*e*k*(n+nu+1+s+I*e)*(n+nu+1+s-I*e)*(n+nu+1+I*t)/((n+nu+1)*(2*n+2*nu+3))

    def be(self, n, nu):
        e, s, m, q = self.eps, self.s, self.m, self.q
        return (-self.lam - s*(s+1) + (n+nu)*(n+nu+1) + e**2 + e*(e-m*q)
                + e*(e-m*q)*(s**2+e**2)/((n+nu)*(n+nu+1)))

    def ga(self, n, nu):
        e, k, t, s = self.eps, self.kappa, self.tau, self.s
        return -I*e*k*(n+nu-s+I*e)*(n+nu-s-I*e)*(n+nu-I*t)/((n+nu)*(2*n+2*nu-1))

    def Rn(self, nu):
        """R_n = a_n/a_{n-1}, n = 1..N (downward continued fraction)."""
        R = {self.N+1: mp.mpf(0)}
        for n in range(self.N, 0, -1):
            R[n] = -self.ga(n, nu)/(self.be(n, nu) + self.al(n, nu)*R[n+1])
        return R

    def Ln(self, nu):
        """L_n = a_n/a_{n+1}, n = -1..-N (upward continued fraction)."""
        L = {-self.N-1: mp.mpf(0)}
        for n in range(-self.N, 0):
            L[n] = -self.al(n, nu)/(self.be(n, nu) + self.ga(n, nu)*L[n-1])
        return L

    def g(self, nu):
        return self.be(0, nu) + self.al(0, nu)*self.Rn(nu)[1] + self.ga(0, nu)*self.Ln(nu)[-1]

    def solve_nu(self, guess):
        self.nu = mp.findroot(self.g, mp.mpmathify(guess), tol=mp.mpf(10)**(-mp.mp.dps+10))
        return self.nu

    def coeffs(self, nu=None):
        nu = self.nu if nu is None else nu
        R, L = self.Rn(nu), self.Ln(nu)
        a = {0: mp.mpf(1)}
        for n in range(1, self.N+1):
            a[n] = R[n]*a[n-1]
        for n in range(-1, -self.N-1, -1):
            a[n] = L[n]*a[n+1]
        return a

    # ---- matching coefficient K_nu (r = 0 form) -------------------------
    def K(self, nu, a):
        e, k, t, epp, s = self.eps, self.kappa, self.tau, self.epp, self.s
        pref = (mp.exp(I*e*k) * mp.exp((s-nu)*(mp.log(2*e)+self.logk)) * mp.power(2, -s)
                * mp.gamma(1-s-2*I*epp) * mp.gamma(2*nu+2)
                / (mp.gamma(nu+1-s+I*e)*mp.gamma(nu+1+I*t)*mp.gamma(nu+1+s+I*e)))
        S1 = mp.fsum((-1)**n * mp.gamma(n+2*nu+1)/mp.factorial(n)
                     * mp.gamma(n+nu+1+s+I*e)/mp.gamma(n+nu+1-s-I*e)
                     * mp.gamma(n+nu+1+I*t)/mp.gamma(n+nu+1-I*t) * a[n]
                     for n in range(0, self.N+1))
        S2 = mp.fsum((-1)**n/(mp.factorial(-n)*mp.rf(2*nu+2, n))
                     * mp.rf(nu+1+s-I*e, n)/mp.rf(nu+1-s+I*e, n) * a[n]
                     for n in range(-self.N, 1))
        return pref*S1/S2

    def amplitudes(self):
        """Return dict with nu, a_n, K_nu, K_{-nu-1}, A_pm, B_inc, B_ref, B_trans."""
        nu, s, e, k = self.nu, self.s, self.eps, self.kappa
        a = self.coeffs(nu)
        aneg = {n: a[-n] for n in a}            # a_n^{-nu-1} = a_{-n}^{nu}
        Kp = self.K(nu, a)
        Km = self.K(-nu-1, aneg)
        sum_a = mp.fsum(a[n] for n in a)
        sum_b = mp.fsum((-1)**n*mp.rf(nu+1+s-I*e, n)/mp.rf(nu+1-s+I*e, n)*a[n] for n in a)
        Ap = (mp.exp(-mp.pi*e/2)*mp.exp(I*mp.pi*(nu+1-s)/2)*mp.power(2, -1+s-I*e)
              * mp.gamma(nu+1-s+I*e)/mp.gamma(nu+1+s-I*e)*sum_a)
        Am = mp.power(2, -1-s+I*e)*mp.exp(-I*mp.pi*(nu+1+s)/2)*mp.exp(-mp.pi*e/2)*sum_b
        om = e/2
        ph = e*mp.log(e) - (1-k)*e/2
        rho = mp.sin(mp.pi*(nu-s+I*e))/mp.sin(mp.pi*(nu+s-I*e))
        Binc = (Kp - I*mp.exp(-I*mp.pi*nu)*rho*Km)*Ap*mp.exp(-I*ph)/om
        Bref = (Kp + I*mp.exp(I*mp.pi*nu)*Km)*Am*mp.exp(I*ph)/om**(1+2*s)
        Btr = ((e*k/om)**(2*s)*mp.exp(I*k*self.epp*(1+2*self.logk/(1+k)))*sum_a
               if abs(1+k) > 0 else mp.nan)   # horizon at r=0 (dual of Schwarzschild): not needed
        X = Km/Kp
        H = (1 + I*mp.exp(I*mp.pi*nu)*X)/(1 - I*mp.exp(-I*mp.pi*nu)*rho*X) - 1
        Sfar = (Am/Ap)*mp.exp(2*I*ph)/om**(2*s)
        return dict(nu=nu, a=a, Kp=Kp, Km=Km, X=X, H=H, Ap=Ap, Am=Am,
                    Binc=Binc, Bref=Bref, Btr=Btr, Sfar=Sfar, S=Bref/Binc, rho=rho)


def X_leading(nu, eps, kappa, logk, tau, s):
    """Leading-order (a_n = delta_n0) K_{-nu-1}/K_nu."""
    e = mp.mpmathify(eps)
    return (mp.exp((2*nu+1)*(mp.log(2*e)+logk)) * mp.gamma(-2*nu)/mp.gamma(2*nu+2)
            * mp.gamma(nu+1-s+I*e)*mp.gamma(nu+1+I*tau)*mp.gamma(nu+1+s+I*e)
            / (mp.gamma(-nu-s+I*e)*mp.gamma(-nu+I*tau)*mp.gamma(-nu+s+I*e)))
