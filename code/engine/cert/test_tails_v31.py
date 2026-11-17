"""Tail control at the SEALED E3/E4 parameters (REFEREE-E item 6; new file, test_tails.py untouched).
Sonic manifold exactly as the runs build it: Problem2, K_s = 16, nu_s = 0.004, defect to N2 = 12 K,
LB_NGRID = 24 (driver_v2.sonic_manifold_v2); centre: K_c = 60, nu_c = 0.03 (Problem.centre_manifold).
Reference: independent mpmath series at 90 digits (engine/bg.py x-series at the sonic point, engine/bgc.py
v-series at the centre), two truncation orders; their difference bounds the reference error.
Tolerance (per component):  eps + rad(cert) + |dy/dx|*(epsX + rad(X)) + 10*referr + 1e-85,
so a violation of order eps ~ 4e-59 is detectable.  Detection power is shown by two can-fail controls:
 INJ: the certified value shifted by 3*allow must be flagged VIOLATION at every point;
 TRUNC: the K-1 truncation (last coefficient dropped) checked against the K tail bound (reported).
usage: test_tails_v31.py [e3|v31]   (e3 = tail_bound exactly as sealed in E3; v31 = with repair R7)"""
import sys, os, time
os.environ['LB_NGRID'] = '24'
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from mpmath import mp, mpf, exp as mexp
from flint import arb, ctx
ctx.prec = 200
import driver_v2 as V2
import param as PM
CODE = sys.argv[1] if len(sys.argv) > 1 else 'v31'
if CODE == 'v31':
    import driver_v31                       # noqa: F401  installs R7 (exact series inputs in tail_bound)
    assert PM.arb_series is driver_v31.arb_series_exact
print('tail_bound code:', 'v3.1 (R7 installed)' if CODE == 'v31' else 'E3 as sealed (param.py unpatched)')
from param import val, horner
import bg as BG, bgc as BGC

DPS = 90
t0 = time.time()
k = sp.Rational(1, 3)
V0s = '0.112439401388098336223763027603135130845'
W0s = '8.90132327775564351911721605677568608689'
om_t = '0.3256888894230684'
nus, nuc = arb('0.004'), arb('0.03')
Pb = V2.Problem2(k, prec=200, K_sonic=16, K_centre=60, K_ty=36, log=lambda *a: None)
sm = V2.sonic_manifold_v2(Pb, arb(V0s), om_t, nus, 16)
cm = Pb.centre_manifold(arb(W0s), nuc)
print('sealed sonic: K=16 nu=0.004 N2=192 ngrid=%s  eps=%s epsX=%s' % (os.environ['LB_NGRID'], sm['eps'].str(5), sm['epsX'].str(5)))
print('sealed centre: K=60 nu=0.03  eps=%s epsX=%s  [%.1fs]' % (cm['eps'].str(5), cm['epsX'].str(5), time.time() - t0))

B = BG.Background(k, dps=DPS, K=120)          # sets mp.dps = 90
mp.dps = DPS
def M(a): return mpf(a.mid().str(110, radius=False))
V0 = mpf(V0s)
N0, om0 = B.sonic_point(V0, om_guess=mpf(om_t))
co1, _ = B.sonic_series(V0, om0, 0, K=120)
co2, _ = B.sonic_series(V0, om0, 0, K=160)
dco = [[n*c[n] for n in range(1, len(c))] for c in co2]
eps, epsX = float(sm['eps'].mid()), float(sm['epsX'].mid())
fails = 0; inj_missed = 0; trunc_flag = []
Yc = [[val(a) for a in sm['co']['Y'][i]] for i in range(3)]
Xc = [val(a) for a in sm['co']['X']]
for fr in ['-0.95', '-0.9', '-0.5', '0.3', '0.9', '0.95']:
    tha = arb(fr)*nus
    Y = [horner(Yc[i], tha) for i in range(3)]
    Yt = [horner(Yc[i][:-1], tha) for i in range(3)]
    X = horner(Xc, tha)
    xt = M(X)
    y1 = B.evalser(co1, xt); y2 = B.evalser(co2, xt); sl = B.evalser(dco, xt)
    referr = max(abs(a - b) for a, b in zip(y1, y2))
    dx = mpf(sm['epsX'].upper().str(30, radius=False)) + M(X.rad()) + abs(M(X) - xt)
    ok = True; worst = mpf(0); ratio = mpf(0); tr = False
    for i in range(3):
        allow = mpf(sm['eps'].upper().str(30, radius=False)) + M(Y[i].rad()) + 2*abs(sl[i])*dx + 10*referr + mpf(10)**-85
        diff = abs(M(Y[i]) - y2[i])
        worst = max(worst, diff); ratio = max(ratio, diff/eps)
        ok &= diff <= allow
        inj_missed += int(abs(M(Y[i]) + 3*allow - y2[i]) <= allow)
        tr |= abs(M(Yt[i]) - y2[i]) > allow
    trunc_flag.append(tr)
    print('  sonic theta=%5s*nu x=%s  max|cert - mp90| = %s  (= %s eps)  referr %s  %s  | TRUNC(K-1) flagged: %s'
          % (fr, mp.nstr(xt, 8), mp.nstr(worst, 4), mp.nstr(ratio, 3), mp.nstr(referr, 3), 'OK' if ok else 'VIOLATION', tr))
    fails += (not ok)

C1 = BGC.Centre(B, K=140).series(mpf(W0s))
C2 = BGC.Centre(B, K=180).series(mpf(W0s))
dC = [[n*c[n] for n in range(1, len(c))] for c in C2]
Zc = [[val(a) for a in cm['co']['Y'][i]] for i in range(3)]
Xcc = [val(a) for a in cm['co']['X']]
for fr in ['0.1', '0.5', '0.9', '0.95']:
    tha = arb(fr)*nuc
    Z = [horner(Zc[i], tha) for i in range(3)]
    Zt = [horner(Zc[i][:-1], tha) for i in range(3)]
    chi = horner(Xcc, tha)
    vt = M(tha)*mexp(M(chi))
    z1 = B.evalser(C1, vt); z2 = B.evalser(C2, vt); sl = B.evalser(dC, vt)
    referr = max(abs(a - b) for a, b in zip(z1, z2))
    dv = abs(vt)*(2*mpf(cm['epsX'].upper().str(30, radius=False)) + M(chi.rad()))
    ok = True; worst = mpf(0); ratio = mpf(0); tr = False
    for i in range(3):
        allow = mpf(cm['eps'].upper().str(30, radius=False)) + M(Z[i].rad()) + 2*abs(sl[i])*dv + 10*referr + mpf(10)**-85
        diff = abs(M(Z[i]) - z2[i])
        worst = max(worst, diff); ratio = max(ratio, diff/float(cm['eps'].mid()))
        ok &= diff <= allow
        inj_missed += int(abs(M(Z[i]) + 3*allow - z2[i]) <= allow)
        tr |= abs(M(Zt[i]) - z2[i]) > allow
    trunc_flag.append(tr)
    print('  centre theta=%4s*nu  max|cert - mp90| = %s  (= %s eps)  referr %s  %s  | TRUNC(K-1) flagged: %s'
          % (fr, mp.nstr(worst, 4), mp.nstr(ratio, 3), mp.nstr(referr, 3), 'OK' if ok else 'VIOLATION', tr))
    fails += (not ok)
print('INJ control (shift 3*allow must be flagged at every component): missed %d of %d -> %s'
      % (inj_missed, 3*10, 'PASS' if inj_missed == 0 else 'FAIL'))
print('TRUNC control (K-1 truncation flagged at %d of %d points)' % (sum(trunc_flag), len(trunc_flag)))
print('TAIL CHECK (sealed parameters, mp 90 digits):', 'PASS' if (fails == 0 and inj_missed == 0) else 'FAIL (%d)' % fails,
      ' total %.1fs' % (time.time() - t0))
