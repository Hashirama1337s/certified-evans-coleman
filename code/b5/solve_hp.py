"""High-precision NON-RIGOROUS Newton for (sig0, zs): centre Chebyshev-Picard guesses (no validation) + sonic
parameterization head; matching at z = zm.  Used only to centre the rigorous box.
usage: solve_hp.py k_num k_den sig0 zs zm [prec]"""
import sys
import json
from flint import arb, arb_poly, ctx
ctx.prec = int(sys.argv[6]) if len(sys.argv) > 6 else 200
import systems as SY
import sonic as SO
from rigode import guess

kk = arb(int(sys.argv[1])) / int(sys.argv[2])
sig0 = arb(sys.argv[3])
zs = arb(sys.argv[4])
zm = arb(sys.argv[5])
M = 30
CS = SY.centre_system(False)
RS = SY.regular_system(False)
z1 = arb('0.1')


ZS_APPROX = float(sys.argv[4])


def mesh(z1, zm):
    """uniform 0.02 steps, then geometric refinement toward the sonic point zs ~ ZS_APPROX"""
    pts = [float(z1.mid())]
    zmf = float(zm.mid())
    while True:
        z = pts[-1]
        h = min(0.02, 0.25 * (ZS_APPROX - z))
        if z + h >= zmf - 1e-12:
            break
        pts.append(z + h)
    return [arb(repr(p)) for p in pts] + [zm]


def centre_at(s0):
    U0 = -2 / (3 * (1 + kk))
    ub, dif = guess(CS, kk, arb(0), z1, [s0, U0], M, iters=150, tol=1e-70)
    u = [z1**2 * ub[0].right(), z1 * ub[1].right()]
    pts = mesh(z1, zm)
    for a, b in zip(pts[:-1], pts[1:]):
        ub, dif = guess(RS, kk, a, b, u, M, iters=150, tol=1e-70)
        u = [ub[0].right(), ub[1].right()]
    return [arb(x.mid()) for x in u]


WHICH = sys.argv[7] if len(sys.argv) > 7 else 'unstable'


def sonic_at(z_s, N=40):
    B = SO.branch(kk, z_s, N=N, which=WHICH, scale=arb('0.2'))
    H = [arb_poly(B['p'][i]) for i in range(3)]
    # Newton for H_z(th) = zm (sign: inner side S>0 -> z < zs)
    th = (zm - z_s) / B['p'][0][1]
    dH = H[0].derivative()
    for it in range(60):
        th = th - (H[0](th) - zm) / dH(th)
        th = arb(th.mid())
    return [arb(H[1](th).mid()), arb(H[2](th).mid())], th, B


def G(s0, z_s):
    c = centre_at(s0)
    s, th, B = sonic_at(z_s)
    return [c[0] - s[0], c[1] - s[1]]


x = [sig0, zs]
for it in range(12):
    g = G(*x)
    eps = arb(10) ** (-(ctx.prec // 8))
    J = []
    g1 = G(x[0] + eps, x[1])
    g2 = G(x[0], x[1] + eps)
    J = [[(g1[0] - g[0]) / eps, (g2[0] - g[0]) / eps], [(g1[1] - g[1]) / eps, (g2[1] - g[1]) / eps]]
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    dx0 = (J[1][1] * g[0] - J[0][1] * g[1]) / det
    dx1 = (-J[1][0] * g[0] + J[0][0] * g[1]) / det
    x = [arb((x[0] - dx0).mid()), arb((x[1] - dx1).mid())]
    print('it', it, 'resid', float(abs(g[0]).mid()), float(abs(g[1]).mid()), 'step', float(abs(dx0).mid()), float(abs(dx1).mid()), flush=True)
    if float(abs(dx0).mid()) + float(abs(dx1).mid()) < 1e-50:
        break
print('sig0', x[0].str(60))
print('zs  ', x[1].str(60))
s, th, B = sonic_at(x[1])
print('theta*', th.str(20), 'lam', B['lam'].str(15), 'class', B['class']['type'])
print('J', [[J[i][j].str(15) for j in range(2)] for i in range(2)])
json.dump({'k': sys.argv[1] + '/' + sys.argv[2], 'sig0': x[0].str(70), 'zs': x[1].str(70), 'zm': sys.argv[5],
           'J': [[J[i][j].str(20) for j in range(2)] for i in range(2)]},
          open('<LAB>/b5/hp_%s_%s%s.json' % (sys.argv[1], sys.argv[2], '' if WHICH == 'unstable' else '_' + WHICH), 'w'), indent=1)
