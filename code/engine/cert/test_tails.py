"""Tail-bound checks (BAR-LB-E controls): certified manifold enclosures vs an
independent 50-digit mpmath evaluation (engine/bg.py x-series at the sonic point,
engine/bgc.py u-series at the centre - different charts, different recursions).
usage: test_tails.py [K] [nu_s] [nu_c]"""
import sys, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from mpmath import mp, mpf, exp as mexp, log as mlog
from flint import arb, ctx
from pipeline import Problem
from core import sym
from param import val, horner
import bg as BG, bgc as BGC

K = int(sys.argv[1]) if len(sys.argv) > 1 else 30
nus = arb(sys.argv[2]) if len(sys.argv) > 2 else arb('0.25')
nuc = arb(sys.argv[3]) if len(sys.argv) > 3 else arb('0.2')
t0 = time.time()
def M(a): return mpf(a.mid().str(60, radius=False))
Kc = int(sys.argv[4]) if len(sys.argv) > 4 else K
Pb = Problem(sp.Rational(1, 3), prec=200, K_sonic=K, K_centre=Kc)
V0s = '0.112439401388098336223763027603135130845'
W0s = '8.90132327775564351911721605677568608689'
sm = Pb.sonic_manifold(arb(V0s), '0.3256888894230684', nus)
print('sonic: lam=%s others=%s eps=%s epsX=%s  [%.1fs]' % (val(sm['lam']).str(10), sm['others'], sm['eps'].str(5), sm['epsX'].str(5), time.time() - t0))
print('   tail info: delta=%s beta=%s L=%s low_defect=%s' % (sm['tinfo']['delta'].str(4), sm['tinfo']['beta'].str(4), sm['tinfo']['L'].str(4), sm['tinfo']['low_defect'].str(3)))
lam_other = sm['others']
cm = Pb.centre_manifold(arb(W0s), nuc)
print('centre: eps=%s epsX=%s  [%.1fs]' % (cm['eps'].str(5), cm['epsX'].str(5), time.time() - t0))
print('   tail info: delta=%s beta=%s L=%s' % (cm['tinfo']['delta'].str(4), cm['tinfo']['beta'].str(4), cm['tinfo']['L'].str(4)))

# independent mpmath reference, 50 digits
mp.dps = 50
B = BG.Background(sp.Rational(1, 3), dps=50, K=80)
V0 = mpf(V0s)
N0, om0 = B.sonic_point(V0, om_guess=mpf('0.3256888894230684'))
co_mp, info = B.sonic_series(V0, om0, 0, K=80)
fails = 0
for fr in ['-0.95', '-0.5', '0.3', '0.95']:
    tha = arb(fr)*nus; th = tha.mid().str(60, radius=False)
    Y = [horner([val(a) for a in sm['co']['Y'][i]], tha) for i in range(3)]
    X = horner([val(a) for a in sm['co']['X']], tha)
    xt = M(X)
    ymp = B.evalser(co_mp, xt)
    # admissible: |Y_true - Y_K| <= eps ; x uncertain by epsX + rad(X) -> |f| * that
    fmax = max(abs(v) for v in [mpf(1)])  # placeholder scale, refined below
    dx = float(sm['epsX'].mid()) + float(X.rad())
    ok = True; worst = 0
    for i in range(3):
        diff = abs(M(Y[i]) - ymp[i])
        allow = float(sm['eps'].mid()) + float(Y[i].rad()) + 10*dx + mpf(10)**-45
        worst = max(worst, float(diff))
        ok &= diff <= allow
    print('  sonic theta=%5s x=%s  max|cert - mp50| = %.3e  (tail eps %.2e)  %s' % (th, mp.nstr(xt, 8), worst, float(sm['eps'].mid()), 'OK' if ok else 'VIOLATION'))
    fails += (not ok)
# centre
C = BGC.Centre(B, K=110)
Zmp_series = C.series(mpf(W0s))
for fr in ['0.1', '0.5', '0.95']:
    tha = arb(fr)*nuc; th = tha.mid().str(60, radius=False)
    Z = [horner([val(a) for a in cm['co']['Y'][i]], tha) for i in range(3)]
    chi = horner([val(a) for a in cm['co']['X']], tha)
    vt = M(tha)*mexp(M(chi))
    zmp = B.evalser(Zmp_series, vt)
    worst = 0; ok = True
    for i in range(3):
        diff = abs(M(Z[i]) - zmp[i])
        allow = float(cm['eps'].mid()) + float(Z[i].rad()) + 10*float(cm['epsX'].mid()) + mpf(10)**-45
        worst = max(worst, float(diff)); ok &= diff <= allow
    print('  centre theta=%5s  max|cert - mp50| = %.3e (tail eps %.2e) %s' % (th, worst, float(cm['eps'].mid()), 'OK' if ok else 'VIOLATION'))
    fails += (not ok)
print('TAIL CHECK:', 'PASS' if fails == 0 else 'FAIL (%d)' % fails, ' total %.1fs' % (time.time() - t0))
