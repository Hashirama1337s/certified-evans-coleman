"""Lohner control for lohner2.lohner_step2(mode='qrw') (the step every E3/E4 run uses) at the SEALED x-leg
parameters (REFEREE-E item 6; test_lohner.py only exercises lohner.py).  New file.
  k = 1/3, x-leg x0 = -0.05 -> xbar (= the run's matching abscissa), K_t = 36, q = 0.25, h0 = 0.002,
  driver2.integrate_rec (the run's integrator), weighted QR installed by importing driver_v3.
Reference: engine/bg.py mpmath Taylor integration at 60 digits (independent code; enclosure widths are
>= 1e-22, so 60 digits resolve them by > 30 orders).
usage: test_lohner2.py thin   thin start: containment at xbar + flow-derivative columns vs mp central differences
       test_lohner2.py box    start cube of half-width 1e-15: all 8 corner trajectories (mp) inside the end hull"""
import sys, os, time
os.environ['LB_NGRID'] = '24'
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from mpmath import mp, mpf
from flint import arb, ctx
ctx.prec = 200
import driver_v3 as V3                     # installs lohner2 qrw as LO.lohner_step
import driver2 as D2
import driver_v2 as V2
import lohner as LO
import param as PM
from core import sym
from param import val
import bg as BG
assert LO.lohner_step is V3._step and V3.MODE == 'qrw'

mode = sys.argv[1]
DPS = 60
t0 = time.time()
k = sp.Rational(1, 3)
Pb = V2.Problem2(k, prec=200, K_sonic=16, K_centre=60, K_ty=36, log=lambda *a: None)
# xbar exactly as the run: mid of ln(thm) + chi(thm) + xc at x~ (thm = 0.9 nu_c)
thm = arb('0.9')*arb('0.03')
cm = Pb.centre_manifold(arb('8.90132327775564351911721605677568608689'), arb('0.03'), tails=False)
xcs = -arb('0.808668956318442829476786547218464827352').log()
xbar = arb(val(thm.log() + PM.horner([val(a) for a in cm['co']['X']], thm) + xcs).mid())
x0 = arb('-0.05')
B = BG.Background(k, dps=DPS, K=60)
mp.dps = DPS
V0 = mpf('0.112439401388098336223763027603135130845')
N0, om0 = B.sonic_point(V0, om_guess=mpf('0.3256888894230684'))
co, _ = B.sonic_series(V0, om0, 0, K=60)
co2, _ = B.sonic_series(V0, om0, 0, K=90)
XA = mpf('-0.05')      # decimal -0.05, the same point as arb('-0.05'); attempt 1 used the binary double -0.05 (2.8e-18 off)
ya = B.evalser(co2, XA)
print('mp start series check |K60 - K90| at x=-0.05: %s' % mp.nstr(max(abs(a - b) for a, b in zip(B.evalser(co, XA), ya)), 3))
xb = mpf(xbar.mid().str(70, radius=False))
print('arb x0 contains mp x0: %s; arb xbar mid == mp xbar: %s' % (x0.contains(arb(mp.nstr(XA, 70))), arb(mp.nstr(xb, 70)).contains(xbar.mid())))
def mpflow(y):
    _, yb = B.integrate(XA, y, xb, frac=mpf('0.15'), maxsteps=5000)
    return yb
def arbs(y): return [arb(mp.nstr(v, DPS + 5)) for v in y]
print('xbar = %s (k=1/3, K_t=36, q=0.25, h0=0.002, mode qrw)' % xbar.str(15))

if mode == 'thin':
    S = LO.LohnerSet(arbs(ya), vecs=[[arb(1) if i == j else arb(0) for i in range(3)] for j in range(3)])
    info = D2.integrate_rec(Pb.F, S, x0, xbar, K=36, q='0.25', h0=0.002, log=print)
    Y = S.hull(); V = S.vecs()
    print('cert: steps=%d time=%.1fs width at xbar %s' % (info['steps'], info['time'], max(v.rad() for v in Y).str(3)))
    yb = mpflow([mpf(mp.nstr(v, DPS + 5)) for v in ya])
    ok = all(Y[i].contains(arb(mp.nstr(yb[i], DPS + 5))) for i in range(3))
    for i in range(3):
        d_ = abs(mpf(Y[i].mid().str(70, radius=False)) - yb[i])
        print('  y%d cert %s  mp %s  |mid - mp| = %s  rad = %s' % (i, Y[i].str(22), mp.nstr(yb[i], 22), mp.nstr(d_, 3), Y[i].rad().str(3)))
    print('THIN containment at xbar:', ok)
    h = mpf(10)**-25; okd = True
    for j in range(3):
        yp = [mpf(mp.nstr(v, DPS + 5)) for v in ya]; ym = list(yp)
        yp[j] += h; ym[j] -= h
        bp, bm = mpflow(yp), mpflow(ym)
        fd = [(bp[i] - bm[i])/(2*h) for i in range(3)]
        cj = all((V[j][i] + sym(arb('1e-30'))).contains(arb(mp.nstr(fd[i], DPS))) for i in range(3))
        okd &= cj
        print('    |mid - fd| = %s  rad = %s' % ([mp.nstr(abs(mpf(V[j][i].mid().str(70, radius=False)) - fd[i]), 3) for i in range(3)], [V[j][i].rad().str(3) for i in range(3)]))
        print('  dPhi/dy%d cert %s  mp-fd %s  contained(+-1e-30 fd error): %s' % (j, [v.str(12) for v in V[j]], [mp.nstr(x, 12) for x in fd], cj))
    print('LOHNER2 THIN TEST:', 'PASS' if (ok and okd) else 'FAIL', ' total %.1fs' % (time.time() - t0))
else:
    w = arb('1e-15')
    c = arbs(ya)
    S = LO.LohnerSet([ci + sym(w) for ci in c])
    info = D2.integrate_rec(Pb.F, S, x0, xbar, K=36, q='0.25', h0=0.002, log=print)
    Y = S.hull()
    print('cert box: steps=%d time=%.1fs width start %s end %s' % (info['steps'], info['time'], w.str(3), max(v.rad() for v in Y).str(3)))
    allok = True; n = 0
    for s0 in (-1, 1):
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                yc = [mpf(mp.nstr(v, DPS + 5)) for v in ya]
                for i, s in enumerate((s0, s1, s2)):
                    yc[i] += s*mpf('1e-15')
                yb = mpflow(yc)
                ok = all(Y[i].contains(arb(mp.nstr(yb[i], DPS + 5))) for i in range(3))
                allok &= ok; n += 1
                print('  corner %s: contained %s' % ((s0, s1, s2), ok))
    print('LOHNER2 BOX TEST (8 corners):', 'PASS' if allok else 'FAIL', ' total %.1fs' % (time.time() - t0))
