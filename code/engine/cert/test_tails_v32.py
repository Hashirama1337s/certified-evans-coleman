"""Tail control for repair R8 (new file; test_tails_v31.py untouched).   usage: test_tails_v32.py [13|12] [r]
REFEREE-R8 fixes applied (seats/REFEREE-R8.md sec. 5): INJ shifts by 3*eps_v32 (not 3*allow) and max allow/eps is
printed; check (C) added (Dual derivative vs central difference + frozen-eigenpair null that must FAIL, logic of the
referee's ref_r8_probe3b.py); mp DPS 200 at k=1/2 (the 90-digit 'resonance' was cancellation, sec. 5d).
The R8 bound eps_v32(box) must bound the TRUE tail  y(V0; theta) - P_K(V0; theta)  at EVERY V0 in the box,
where P_K(V0) has the exact coefficients at V0.  Checked at sample points V0 = V0~ + f r (f in -1..1):
 (A) mp reference (engine/bg.py x-series at the sonic point of that V0, two truncation orders) against the POINT
     polynomial P_K(V0s):  |P_K(V0s; theta) - y_mp(x(theta))| <= allow = eps_v32 + rad + 2|y'| dx + 10 referr + 1e-85.
     NOTE (referee 5b): the true tail is point-sized, so (A) is a point-tail control; it cannot see the r-term.
     INJ: the certified value shifted by 3*eps_v32 must be flagged (meaningful only where allow < 2 eps; allow/eps printed).
 (B) mean-value head at each sample V0s: every c_n(V0s), n > K, meets [c_n(V0~) +- r|c_n'(box)|], and
     |c_n(V0s)| <= h_n = min(|box c_n|, |c_n(V0~)| + r|c_n'(box)|)  (the bound the driver feeds to tail_bound).
 (C) sensitive derivative control (the part (A)/(B) cannot detect): at the thin V0~ the Dual derivative of every
     head coefficient must overlap the central difference (c_n(V0~+h) - c_n(V0~-h))/2h, h = 1e-25 (truncation error
     O(h^2) ~ 1e-50 relative, below the ball radius ~1e-35); the same test with the eigenpair derivative FROZEN
     (lam, v treated as V0-independent) must FAIL on at least one coefficient.
 k=13: sealed target parameters (V0~ = 40-digit data, nu_s = 0.004).  k=12: float64 warm-up seed, nu_s = 0.001."""
import sys, os, time
os.environ['LB_NGRID'] = '24'
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from mpmath import mp, mpf
from flint import arb, ctx
ctx.prec = 200
import driver_v32 as V32
import driver_v31 as V31
import driver_v2 as V2
import param as PM
from param import val, horner
from core import Dual, sym
import bg as BG
assert PM.arb_series is V31.arb_series_exact

KK = sys.argv[1] if len(sys.argv) > 1 else '13'
RAD = sys.argv[2] if len(sys.argv) > 2 else '1e-12'
t0 = time.time()
if KK == '13':
    k = sp.Rational(1, 3); V0t = arb('0.112439401388098336223763027603135130845'); om_t = '0.3256888894230684'; nus = arb('0.004'); DPS = 90
else:
    k = sp.Rational(1, 2); V0t = arb(float((-2*arb(2).sqrt()/5).mid())); om_t = '0.7058823529411765'; nus = arb('0.001'); DPS = 200
K = 16; N2 = 12*K
r = arb(RAD)
Pb = V2.Problem2(k, prec=200, K_sonic=16, K_centre=60, K_ty=36, log=lambda *a: None)
smb = V32.sonic_manifold_v32(Pb, Dual(V0t + sym(r), [arb(1)]), om_t, nus, K)
eps_b, epsX_b = smb['eps'], smb['epsX']
sm31 = V2.sonic_manifold_v2(Pb, Dual(V0t + sym(r), [arb(1)]), om_t, nus, K)
print('k=%s nu_s=%s K=%d box radius %s: eps_v32 %s epsX_v32 %s | v3.1 eps %s epsX %s  [%.1fs]' % (
    k, nus, K, RAD, eps_b.str(4), epsX_b.str(4), sm31['eps'].str(4), sm31['epsX'].str(4), time.time() - t0), flush=True)
# head enclosures exactly as the driver builds them (point at V0~, derivative over the box, box series)
pt_ins = [V32._ser([val(a) for a in smb['smc']['co']['Y'][i]]) for i in range(3)]
du_ins = [V32._dual_series(smb['co']['Y'][i]) for i in range(3)]
box_ins = [V32._ser([val(a) for a in smb['co']['Y'][i]]) for i in range(3)]
cp_ = Pb.Gs.series_eval(pt_ins, N2)
cd_ = V32.series_eval_dual(Pb.Gs, du_ins, N2)
cb_ = Pb.Gs.series_eval(box_ins, N2)

# ---------------- (C) derivative vs central difference at the thin V0~, and the frozen-eigenpair null
V0c = arb(V0t.mid()); hh = arb(arb('1e-25').mid())
def cser(V0):
    sm = V2.sonic_manifold_v2(Pb, Dual(V0, [arb(1)]), om_t, nus, K, tails=False)
    return Pb.Gs.series_eval([V32._ser([val(a) for a in sm['co']['Y'][i]]) for i in range(3)], N2)
def dser():
    sm = V2.sonic_manifold_v2(Pb, Dual(V0c, [arb(1)]), om_t, nus, K, tails=False)
    return V32.series_eval_dual(Pb.Gs, [V32._dual_series(sm['co']['Y'][i]) for i in range(3)], N2)
cplus, cminus = cser(V0c + hh), cser(V0c - hh)
built = dser()
orig = Pb.eigen
def eig_frozen(y0, branch='fast'):
    lam, v, DG, others = orig(y0, branch)
    fz = lambda x: Dual(val(x), [arb(0)]) if isinstance(x, Dual) else x
    return fz(lam), [fz(x) for x in v], DG, others
Pb.eigen = eig_frozen
frozen = dser()
Pb.eigen = orig
cres = {}
for tag, du in [('built', built), ('frozen-eigenpair null', frozen)]:
    nb = 0; ov = 0; worst = 0.0
    for o in range(3):
        a = V32._coeffs(cplus[o], N2); b = V32._coeffs(cminus[o], N2); d = V32._coeffs(du[o].d[0], N2)
        for n in range(K + 1, N2):
            fd = (a[n] - b[n])/(2*hh)
            if float(fd.mid()) == 0:
                continue
            nb += 1; ov += int(d[n].overlaps(fd))
            worst = max(worst, abs(float((d[n] - fd).mid()))/abs(float(fd.mid())))
    cres[tag] = (ov, nb)
    print('(C) derivative "%s" vs central difference (h=1e-25): %d of %d head coefficients overlap; max rel err %.3g' % (tag, ov, nb, worst), flush=True)
C_ok = cres['built'][0] == cres['built'][1] and cres['frozen-eigenpair null'][0] < cres['frozen-eigenpair null'][1]
print('(C) built derivative agrees everywhere: %s; frozen null FAILS (>=1 non-overlap): %s -> %s' % (
    cres['built'][0] == cres['built'][1], cres['frozen-eigenpair null'][0] < cres['frozen-eigenpair null'][1], 'PASS' if C_ok else 'FAIL'), flush=True)

B = BG.Background(k, dps=DPS, K=120)
mp.dps = DPS
DIG = DPS + 20
def M(a): return mpf(a.mid().str(DIG, radius=False))
EPS = mpf(eps_b.upper().str(30, radius=False))
fails = 0; inj_missed = 0; ntests = 0; tm_fail = 0; h_fail = 0; ratio_max = mpf(0); allow_eps = mpf(0)
branch = None
for f in ['-1', '-0.5', '0', '0.7', '1']:
    V0s = V0t + arb(f)*r                                  # exact sample point
    sp_ = V2.sonic_manifold_v2(Pb, Dual(V0s, [arb(1)]), om_t, nus, K, tails=False)
    Yc = [[val(a) for a in sp_['co']['Y'][i]] for i in range(3)]
    Xc = [val(a) for a in sp_['co']['X']]
    # (B) mean-value head containment and |c_n(V0s)| <= h_n
    cs = Pb.Gs.series_eval([V32._ser(Yc[i]) for i in range(3)], N2)
    for o in range(3):
        a = V32._coeffs(cs[o], N2); b = V32._coeffs(cp_[o], N2); d = V32._coeffs(cd_[o].d[0], N2); bx = V32._coeffs(cb_[o], N2)
        for n in range(K + 1, N2):
            if not (b[n] + sym(r*d[n].abs_upper())).overlaps(a[n]):
                tm_fail += 1
            hn = bx[n].abs_upper().min(b[n].abs_upper() + r*d[n].abs_upper())
            if not (a[n].abs_lower() <= hn):
                h_fail += 1
    # (A) mp reference at this V0
    V0m = M(V0s)
    N0, om0 = B.sonic_point(V0m, om_guess=mpf(om_t))
    if branch is None:
        best = None
        for br in (0, 1):
            co, _ = B.sonic_series(V0m, om0, br, K=40)
            x1 = M(horner(Xc, arb('0.5')*nus)); y1 = B.evalser(co, x1)
            dif = max(abs(M(horner(Yc[i], arb('0.5')*nus)) - y1[i]) for i in range(3))
            if best is None or dif < best[0]: best = (dif, br)
        branch = best[1]
        print('  mp branch matched to the certified manifold: %d (diff %s); mp dps %d' % (branch, mp.nstr(best[0], 3), DPS), flush=True)
    for K1, K2 in [(120, 160), (60, 80), (40, 50)]:
        try:
            co1, _ = B.sonic_series(V0m, om0, branch, K=K1)
            co2, _ = B.sonic_series(V0m, om0, branch, K=K2)
            break
        except ZeroDivisionError:
            print('  (mp reference: zero slope at K=%d/%d, trying lower orders)' % (K1, K2), flush=True)
    dco = [[n*c[n] for n in range(1, len(c))] for c in co2]
    worst = mpf(0)
    for fr in ['-0.95', '-0.9', '-0.5', '0.5', '0.95']:
        tha = arb(fr)*nus
        Y = [horner(Yc[i], tha) for i in range(3)]
        X = horner(Xc, tha)
        xt = M(X)
        y1 = B.evalser(co1, xt); y2 = B.evalser(co2, xt); sl = B.evalser(dco, xt)
        referr = max(abs(a - b) for a, b in zip(y1, y2))
        dx = mpf(epsX_b.upper().str(30, radius=False)) + M(X.rad()) + abs(M(X) - xt)
        ok = True
        for i in range(3):
            allow = EPS + M(Y[i].rad()) + 2*abs(sl[i])*dx + 10*referr + mpf(10)**-85
            allow_eps = max(allow_eps, allow/EPS)
            diff = abs(M(Y[i]) - y2[i])
            worst = max(worst, diff); ratio_max = max(ratio_max, diff/EPS)
            ok &= diff <= allow
            inj_missed += int(abs(M(Y[i]) + 3*EPS - y2[i]) <= allow)
            ntests += 1
        fails += (not ok)
    print('  V0 = V0~ %+s r: max|P_K(V0) - mp| over 5 theta = %s  (%s x eps_v32); reference orders %d/%d  %s' % (
        f, mp.nstr(worst, 4), mp.nstr(worst/EPS, 3), K1, K2, 'OK' if fails == 0 else 'VIOLATION'), flush=True)
nh = 5*3*(N2 - K - 1)
print('(B) mean-value head: enclosure failures %d of %d; |c_n(V0s)| <= h_n failures %d of %d' % (tm_fail, nh, h_fail, nh))
print('(A) max |error| / eps_v32 = %s; max allow/eps_v32 = %s; INJ (3*eps_v32 shift) missed %d of %d -> %s' % (
    mp.nstr(ratio_max, 3), mp.nstr(allow_eps, 4), inj_missed, ntests, 'PASS' if inj_missed == 0 else 'FAIL'))
allok = fails == 0 and inj_missed == 0 and tm_fail == 0 and h_fail == 0 and C_ok
print('TAIL CHECK R8 (k=%s, box %s, mp %d digits):' % (k, RAD, DPS), 'PASS' if allok else 'FAIL', ' total %.1fs' % (time.time() - t0))
