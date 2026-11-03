"""Exchange-parity decomposition of the physical (unit-flux) absorption probability, s=-2, l=2:
   1 - |S_+|^2 = [1 - |S_even|^2] + Delta_odd + |S_odd|^2,   Delta_odd = -2 Re[S_even^* S_odd],
   S_even/odd = (S_+ +- S_-^(eta))/2, both with the unit-flux far-zone factor F_P (P = +1).
The image S_- is computed directly from MST with kappa -> -kappa, log(-kappa) = log kappa + i pi eta (the law is
not used).  Coefficients of eps^6 are extracted by quadratic Richardson extrapolation from eps = 4e-3, 2e-3, 1e-3.
Predictions:  Delta_odd -> (pi/225) kappa^5 w(1+w^2)(4+w^2) n_eta(w) eps^6,  w = m q/kappa,
              n_{+1} = 1/(e^{2 pi w}-1),  n_{-1} = e^{2 pi w}/(e^{2 pi w}-1)   [= (2/225) kappa^5 at w = 0];
              m = 0:  P_H -> (2/225) kappa^4 (1+kappa) eps^6  (Starobinsky-Page), so Delta_odd/P_H = kappa/(1+kappa);
              q = 0:  P_H -> (4/225) eps^6 (Page).
"""
import mpmath as mp
from mst import I
from flux_check import analyse, branch
from csm import D_TS
mp.mp.dps = 60


def parts(l, m, q, eps, eta):
    o = analyse(l, m, q, eps)
    d, A = o['d'], o['A']
    om = mp.mpf(eps)/2
    CP = o['D'] + 12*I*om
    Sp = (-1)**(l+1)*CP/(16*om**4)*A['Bref']/A['Binc']          # unit-flux physical S (P = +1)
    _, Am = branch(d, -1, eta)
    Sm = (-1)**(l+1)*CP/(16*om**4)*Am['Bref']/Am['Binc']        # image, same far zone and normalization
    Se, So = (Sp + Sm)/2, (Sp - Sm)/2
    return dict(PH_flux=o['PH'], PH_S=1 - abs(Sp)**2, Ev=1 - abs(Se)**2, Dodd=-2*mp.re(mp.conj(Se)*So),
                So2=abs(So)**2, flux=o['flux_res'])


def richardson(vals, hs):
    # fit v(h) = c0 + c1 h + c2 h^2 through three points
    M = mp.matrix([[1, h, h**2] for h in hs]); return mp.lu_solve(M, mp.matrix(vals))[0]


if __name__ == '__main__':
    l = 2
    eps_list = [mp.mpf('4e-3'), mp.mpf('2e-3'), mp.mpf('1e-3')]
    cases = [(2, 0), (0, '0.6'), (0, '0.3'), (2, '0.3'), (1, '0.5'), (-2, '0.4')]
    print(" m    q   eta |  Delta_odd/eps^6      prediction           |  P_H/eps^6 (eps->0)   [1-|S_even|^2]/eps^6"
          "   Delta_odd/P_H   | flux res.")
    for m, q in cases:
        q = mp.mpf(q); k = mp.sqrt(1 - q**2); w = m*q/k
        for eta in (1, -1):
            rows = [parts(l, m, q, e, eta) for e in eps_list]
            D6 = richardson([r['Dodd']/e**6 for r, e in zip(rows, eps_list)], eps_list)
            E6 = richardson([r['Ev']/e**6 for r, e in zip(rows, eps_list)], eps_list)
            if m == 0 or q == 0:
                pred = mp.mpf(2)/225*k**5
                PH6 = richardson([r['PH_flux']/e**6 for r, e in zip(rows, eps_list)], eps_list)
                ph = f"{mp.nstr(PH6, 10):>14}"
                ratio = f"{mp.nstr(D6/PH6, 8):>10} (k/(1+k) = {mp.nstr(k/(1+k), 8)})"
            else:
                n = 1/(mp.exp(2*mp.pi*w) - 1) if eta == 1 else mp.exp(2*mp.pi*w)/(mp.exp(2*mp.pi*w) - 1)
                pred = mp.pi/225*k**5*w*(1 + w**2)*(4 + w**2)*n
                ph = f"{'(eps^5 term)':>14}"
                ratio = ''
            fl = max(r['flux'] for r in rows)
            idt = max(abs(r["PH_S"] - r["Ev"] - r["Dodd"] + r["So2"]) for r in rows)
            print(f"{m:+d} {mp.nstr(q,2):>4} {eta:+d} | {mp.nstr(D6, 12):>16} {mp.nstr(pred, 12):>16}"
                  f" ({mp.nstr(D6/pred - 1, 2):>8}) | {ph} {mp.nstr(E6, 10):>16}   {ratio} | {mp.nstr(fl, 2)} id {mp.nstr(idt, 2)}")
