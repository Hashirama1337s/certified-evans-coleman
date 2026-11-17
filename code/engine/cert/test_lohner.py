"""Lohner integrator tests.
usage: test_lohner.py a b K q [tmax]
 - thin start from 50-digit mpmath shooting at x=a (engine/bg.py), certified Lohner to x=b,
 - compare with mpmath at b (must be contained), report width, steps, time,
 - also carry the 3 unit tangent vectors (flow derivative) and compare with mp finite differences."""
import sys, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, ctx
from mpmath import mp, mpf
import bg as BG
from pipeline import Problem
import lohner as LO

a, b, K, q = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
tmax = float(sys.argv[5]) if len(sys.argv) > 5 else None
k = sp.Rational(1, 3)
mp.dps = 50
B = BG.Background(k, dps=50, K=50)
V0 = mpf('0.112439401388098336223763027603135130845')
N0, om0 = B.sonic_point(V0, om_guess=mpf('0.3256888894230684'))
co, info = B.sonic_series(V0, om0, 0, K=50)
xa, ya = B.integrate(mpf(-0.05), B.evalser(co, mpf(-0.05)), mpf(a), frac=mpf('0.15'))
xb, yb = B.integrate(mpf(a), ya, mpf(b), frac=mpf('0.15'))
Pb = Problem(k, prec=200)
S = LO.LohnerSet([arb(mp.nstr(v, 55)) for v in ya],
                 vecs=[[arb(1) if i == j else arb(0) for i in range(3)] for j in range(3)])
t0 = time.time()
info = LO.integrate(Pb.F, S, a, b, K=K, q=q, log=print, tmax=tmax)
Y = S.hull()
print('done=%s steps=%d time=%.1fs reached x=%s' % (info['done'], info['steps'], info['time'], info['x'].str(6)))
if info['done']:
    ok = all(Y[i].contains(arb(mp.nstr(yb[i], 55))) for i in range(3))
    print('width at b:', max(v.rad() for v in Y).str(3), ' contains mp50:', ok)
    for i in range(3):
        print('  y%d cert %s  mp %s' % (i, Y[i].str(20), mp.nstr(yb[i], 20)))
    # flow derivative check (mp central differences, column 0)
    eps = mpf(10)**-20
    yp = list(ya); yp[0] += eps; ym = list(ya); ym[0] -= eps
    _, bp = B.integrate(mpf(a), yp, mpf(b), frac=mpf('0.15')); _, bm = B.integrate(mpf(a), ym, mpf(b), frac=mpf('0.15'))
    fd = [(bp[i] - bm[i])/(2*eps) for i in range(3)]
    V = S.vecs()[0]
    print('  dPhi/dN(a) cert', [v.str(10) for v in V], ' mp-fd', [mp.nstr(x, 10) for x in fd])
    print('  tangent width', max(x.rad() for x in V).str(3))
