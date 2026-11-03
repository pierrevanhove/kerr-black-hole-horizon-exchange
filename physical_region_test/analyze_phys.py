import json, sys, mpmath as mp
mp.mp.dps = 90
fname = sys.argv[1] if len(sys.argv) > 1 else 'phys_run.json'
RE, NE, RQ, NQ = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('2e-5'), 12, mp.mpf('0.35'), 28
data = json.load(open(fname))
D = {(r['l'], r['sign'], r['j']): r['pts'] for r in data}
C = lambda z: mp.mpc(mp.mpf(z[0]), mp.mpf(z[1]))
eps = [RE*mp.expjpi(2*(j + mp.mpf(1)/2)/NE) for j in range(NE)]
qs = [RQ*mp.expjpi(2*(k + mp.mpf(1)/2)/NQ) for k in range(NQ)]
Y0 = lambda L: mp.sqrt(mp.mpf(2*L+1)/(4*mp.pi))

def coeffs(F, KMAX=9, PMAX=9):
    return {(k, p): mp.fsum(F[j][kk]*eps[j]**(-k)*qs[kk]**(-p) for j in range(NE) for kk in range(NQ))/(NE*NQ)
            for k in range(KMAX) for p in range(PMAX)}

def St(l, sg):
    return [[C(D[(l, sg, j)][k]['St']) for k in range(NQ)] for j in range(NE)]

def P2(sg2):
    F = [[mp.mpc(0)]*NQ for _ in range(NE)]
    for l in range(2, 8):
        sg = sg2 if l == 2 else 1
        for j in range(NE):
            for k in range(NQ):
                pt = D[(l, sg, j)][k]
                b = {int(L): C(v) for L, v in pt['b'].items()}
                S0 = mp.fsum(v*Y0(L) for L, v in b.items())
                F[j][k] += (C(pt['St']) - 1)*S0*b[2]      # drop the forward delta: sum_l S_l(0) b_2 = Y_2(0) exactly
    return F

A2p, A2m = coeffs(St(2, 1)), coeffs(St(2, -1))
Pp, Pm = coeffs(P2(1)), coeffs(P2(-1))
z3, z5 = mp.zeta(3), mp.zeta(5)
print("=== L=2 polar projection, physical branch: f^(2,2) coefficient of eps_bar w^(p-1)|a|^p = (2^(p+1) pi/i) A_{p+1,p} / sqrt(pi/5)")
for p in range(0, 7):
    v = (2**(p+1)*mp.pi/I if False else 2**(p+1)*mp.pi/mp.j)*Pp[(p+1, p)]/mp.sqrt(mp.pi/5)
    print(f"  p={p}: {mp.nstr(v, 16)}")
bb = (mp.mpf(2048)/45)*(z3 - z5) - mp.mpf(4706)/945
print(f"  BBCPSS (C.1) at gamma=0, p=5:  (2048/45)(z3-z5) - 4706/945 = {mp.nstr(bb, 16)}")
TJ = json.load(open('../conservative_test/taylor2.json'))
for eta in (1, -1):
    T2 = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in TJ[f"{eta}_2"]]
    print(f"  super-extremal eta={eta:+d}: L=2 Taylor coeffs vs physical 2^p A_(p+1,p):")
    for p in range(1, 7):
        phys = 2**p*Pp[(p+1, p)]
        print(f"     p={p}:  SE {mp.nstr(T2[p], 14):>32}   phys {mp.nstr(phys, 14):>32}   diff {mp.nstr(phys - T2[p], 6)}")
print("\n=== eps_bar^0 |a|^6 term:  (128 pi/i) A_{6,6} / sqrt(pi/5)")
v66 = 128*mp.pi/mp.j*Pp[(6, 6)]/mp.sqrt(mp.pi/5)
z2, z4, z6 = mp.zeta(2), mp.zeta(4), mp.zeta(6)
c1 = mp.j/15*(10*z2 - 92*z4 + 80*z6 - 1); c2 = mp.j/120*(120*z2 - 1152*z4 + 1024*z6 - 1); c3 = mp.j/90*(60*z2 - 600*z4 + 544*z6 - 5)
print(f"  ours: {mp.nstr(v66, 16)}   BBCPSS (4.6) at gamma=0: (27c1+92c2+27c3)/6 = {mp.nstr((27*c1 + 92*c2 + 27*c3)/6, 16)}")
print("\n=== all nonzero coefficients A_{k,p} of the L=2 projection (physical), |A|>1e-12, k<=8, p<=8:")
for k in range(9):
    row = [(p, Pp[(k, p)]) for p in range(9) if abs(Pp[(k, p)]) > mp.mpf('1e-12')]
    print(f"  k={k}: " + "  ".join(f"p{p}:{mp.nstr(v, 10)}" for p, v in row))
print("\n=== duality: A(+) - A(-) for the l=2 partial wave (|diff|>1e-14):")
for k in range(9):
    row = [(p, A2p[(k, p)] - A2m[(k, p)]) for p in range(9) if abs(A2p[(k, p)] - A2m[(k, p)]) > mp.mpf('1e-14')]
    if row: print(f"  k={k}: " + "  ".join(f"p{p}:{mp.nstr(v, 10)}" for p, v in row))
json.dump({f"{s}_{k}_{p}": [str(mp.re(v)), str(mp.im(v))] for s, A in (('P+', Pp), ('P-', Pm), ('l2+', A2p), ('l2-', A2m)) for (k, p), v in A.items()},
          open(fname.replace('.json', '_coeffs.json'), 'w'))
