import json, mpmath as mp
from fractions import Fraction
mp.mp.dps = 60
R = {tuple(int(t) for t in k.split('_')): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in json.load(open('l2_rich.json')).items()}
odd = lambda p: mp.re((R[(1, 6, p)] - R[(-1, 6, p)])/2)
evn = lambda p: mp.im((R[(1, 6, p)] + R[(-1, 6, p)])/2)
z = {n: mp.zeta(n) for n in range(2, 9)}
def ident(name, x, basis, labels, tol=mp.mpf('1e-15'), maxc=10**7):
    mp.mp.dps = 25
    rel = mp.pslq([x] + basis, tol=tol, maxcoeff=maxc, maxsteps=10**7)
    mp.mp.dps = 60
    if rel is None:
        print(f"  {name}: no relation with {labels}"); return None
    out = {lab: Fraction(-r, rel[0]) for r, lab in zip(rel[1:], labels)}
    resid = x - mp.fsum(mp.mpf(f.numerator)/f.denominator*b for f, b in zip(out.values(), basis))
    print(f"  {name} = " + " + ".join(f"({f}){lab}" for lab, f in out.items() if f != 0) + f"    [residual {mp.nstr(resid, 2)}]")
    return out
print("duality-odd (absorption) at eps^6 q^p:")
ident("p=0", odd(0), [mp.mpf(1), mp.pi], ['1', 'pi'])
ident("p=1", odd(1), [mp.mpf(1), mp.pi], ['1', 'pi'])
ident("p=2", odd(2), [mp.mpf(1), mp.pi, mp.pi**2], ['1', 'pi', 'pi^2'])
ident("p=3", odd(3), [mp.mpf(1), mp.pi, mp.pi**3], ['1', 'pi', 'pi^3'])
ident("p=4", odd(4), [mp.mpf(1), mp.pi**2, mp.pi**4], ['1', 'pi^2', 'pi^4'])
ident("p=6", odd(6), [mp.mpf(1), z[2], z[4], z[6]], ['1', 'z2', 'z4', 'z6'])
print("duality-even, imaginary (phase) at eps^6 q^p:")
x5 = evn(5)
ident("p=5", x5, [mp.mpf(1), z[3], z[5]], ['1', 'z3', 'z5'])
ident("p=5 (test with pi^2 added)", x5, [mp.mpf(1), z[3], z[5], mp.pi**2], ['1', 'z3', 'z5', 'pi^2'], tol=mp.mpf('1e-16'))
ident("p=7", evn(7), [mp.mpf(1), z[3], z[5], z[7]], ['1', 'z3', 'z5', 'z7'], tol=mp.mpf('1e-16'))
# l=2 partial wave alone in BBCPSS units: -160 i A_{6,5}
print("\nl=2 partial wave alone, BBCPSS units (-160 i A_65):", mp.nstr(160*x5, 16),
      "  minus (2048/45)(z3-z5):", mp.nstr(160*x5 - mp.mpf(2048)/45*(z[3]-z[5]), 16))
