"""BAR-LB-E5 runs (v3.2 = v3.1 + R8 mean-value sonic tail head, driver_v32.MatcherV32).   usage: run_v32.py target | nullV0 | warmup | null2
  target : k=1/3: item-2 exact identities; radius ladder (PE1); saddle; PE2 v3.1 (det sign, x monotone,
           mono-run direction, inflated boxes); CONTENT NULLS along W0 and xc (radius 1e-15, centre 1e-13
           off x~, Krawczyk must say K not inside X; reuses the cached sonic legs); outward leg + m/r (PE3).
  nullV0 : k=1/3 CONTENT NULL along V0 (own sonic legs; ~2 residual runs).
  warmup : k=1/2 positive control (PE4), seeded from the float64-rounded Friedmann data
           (AMEND-LB-E-1: 'seeded ... from a float guess'); E3 seeded at the 200-bit data (disclosed).
  null2  : SMOKE CHECK ONLY (not a control): singular centre branch, P = 0 invariant -> R[0] != 0 by construction.
Sealed parameters as BAR-LB-E3 (PREC 200, K_s 16, K_c 60, K_t 36, nu_c 0.03, nu_s ladder 0.004/0.002/0.001
chosen at box radius 1e-9, LB_NGRID 24, x_switch -0.05 / -0.1, X1 = 12, radius ladder 1e-15..1e-9).
No environment flag may be set (checked below)."""
import sys, os, time
FLAGS = ['LB_SMOKE', 'LB_EXPANDED', 'LB_MVSTART', 'LB_LADDER', 'LB_V3_LADDER', 'LB_NULADDER', 'LB_NU_AT',
         'LB_X1', 'LB_LOHNER_MODE', 'LB_V3ALIGN', 'LB_PREC', 'LB_KS', 'LB_KC', 'LB_KT', 'LB_NUC']
_set = [f for f in FLAGS if os.environ.get(f) is not None]
if _set:
    raise SystemExit('run_v32: environment flags must be unset: %s' % _set)
if os.environ.get('LB_NGRID', '24') != '24':
    raise SystemExit('run_v32: LB_NGRID must be 24 (sealed)')
os.environ['LB_NGRID'] = '24'
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, ctx
ctx.prec = 200
import driver_v31 as V31
import driver_v32 as V32
import driver_v3 as V3
import driver_v2 as V2
import driver2 as D2
import driver as DR
import lohner as LO
import param as PM
from param import val
from core import sym
assert LO.lohner_step is V3._step and V3.MODE == 'qrw', 'weighted-QR Lohner not installed'
assert PM.arb_series is V31.arb_series_exact, 'R7 (exact series inputs in tail_bound) not installed'
assert V32.MatcherV32.start_parts is not V31.MatcherV31.start_parts, 'R8 (mean-value sonic tail head) not installed'

mode = sys.argv[1]
PREC, KS, KC, KT, NUC = 200, 16, 60, 36, '0.03'
NU_LADDER = ['0.004', '0.002', '0.001']
LADDER = ['1e-15', '1e-14', '1e-13', '1e-12', '1e-11', '1e-10', '1e-9']
NU_AT = arb('1e-9')
NULL_OFF, NULL_RAD = '1e-13', '1e-15'
T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
log('[run_v32] Lohner step = lohner2.lohner_step2(mode=%s); matcher = MatcherV32 (R8); ladder %s' % (V3.MODE, ','.join(LADDER)))

def setup(k, V0, om, xsw):
    Pb = V31.Problem2(k, prec=PREC, K_sonic=KS, K_centre=KC, K_ty=KT)
    nu, _ = V2.choose_nu(Pb, V0 + sym(NU_AT), om, NU_LADDER, KS, log=log)
    log('sealed rule: nu_s = %s (first in ladder contracting at box radius %s)' % (nu, NU_AT.str(3)))
    return Pb, V32.MatcherV32(Pb, nu, NUC, Ks=KS, x_switch=xsw, Kt=KT)

def show(Kx, Xb):
    for i, nm in enumerate(['V0', 'W0', 'xc']):
        log('  %s: K = %s   X = %s' % (nm, Kx[i].str(25), Xb[i].str(25)))

def show_repairs(d, tag):
    log('%s repairs: leak4a(enclosed)=%s  e_t=%s  e^{L dmax}(switch)=%s  e0=%s  e_end=%s  e^{L dmx}(end)=%s  Erad=%s'
        % (tag, d['leak4a'].str(3), d['et_max'].str(3), d['g0'].str(6), d['e0'].str(3), d['e_end'].str(3), d['gt'].str(6), d['Erad'].str(3)))
    tw = max(t.rad() for t in d['tan']); dw = max((a - b).abs_upper() for a, b in zip(d['tan'], d['tan_sur']))
    log('%s 4c switch tangent: rigorous width %s, |rigorous - v3 surrogate| <= %s' % (tag, tw.str(3), dw.str(3)))

def null_report(res, name):
    if res['void']:
        log('CONTENT NULL %s: VOID (exception, not a pass): %s' % (name, res['error'])); return False
    log('CONTENT NULL %s: Krawczyk contraction=%s (must be False); K inside X per component %s; |mid K - mid X|/r = %.3g; rad K/rad X = %.3g; Ebox %s; cached sonic legs point/box %s/%s'
        % (name, res['contracts'], res['inside'], res['miss_over_r'], res['excess'], res['Ebox'].str(3), res['cached_point'], res['cached_box']))
    log('CONTENT NULL %s PASS: %s' % (name, res['passed']))
    return res['passed']

k13 = sp.Rational(1, 3)
V0s = arb('0.112439401388098336223763027603135130845'); W0s = arb('8.90132327775564351911721605677568608689')
xcs = -arb('0.808668956318442829476786547218464827352').log(); om_t = '0.3256888894230684'

if mode == 'target':
    log('item 2 exact identities (k=1/3):', V31.exact_chart_checks(k13))
    Pb, M = setup(k13, V0s, om_t, '-0.05')
    Xt = [V0s, W0s, xcs]
    L = D2.krawczyk_ladder(M, Xt, om_t, ladder=LADDER, log=log)
    log('PE1 ladder:', L['tried'])
    show_repairs(L['point'], 'point')
    if L['ok']:
        log('PE1: first contracting radius r =', L['r']); show(L['Kx'], L['Xb'])
        log('V0 certified width:', (2*L['Kx'][0].rad()).str(3))
        show_repairs(L['box'], 'box')
        sc = DR.saddle_check(L['box']['sm'])
        log('PE2 saddle:', sc['saddle'], ' lam=%s lam_b=%s' % (sc['lam'].str(12), sc['lam_b'].str(12)))
        pe = V31.pe2_v31(Pb, M, L['box'])
        log('PE2 v3.1:', {k: (v.str(4) if hasattr(v, 'str') else v) for k, v in pe.items()})
        n_ok = []
        for ax, nm in [(1, 'W0'), (2, 'xc')]:
            n_ok.append(null_report(V31.content_null(M, Xt, ax, om_t, off=NULL_OFF, rad=NULL_RAD), nm))
        log('CONTENT NULLS (W0, xc) ALL PASS:', all(n_ok))
        V0box = L['Kx'][0].intersection(L['Xb'][0])
        ow = V2.outward_leg_v2(Pb, M, V0box, om_t, X1=12.0, log=log)
        log('outward: steps=%d D_sonic<0: %s V>0: %s lemma: %s' % (ow['steps'], ow['ok_D'], ow['ok_V'], ow.get('lemma')))
        if ow.get('lemma'):
            mr = ow['m_over_r']
            log('PE3 m/r(t=0) =', mr.str(10), ' inside 0.0596 +- 0.0005:', bool(mr > arb('0.0591')) and bool(mr < arb('0.0601')))
    else:
        log('PE1: no contraction on the ladder -> FAILED')

elif mode == 'nullV0':
    Pb, M = setup(k13, V0s, om_t, '-0.05')
    null_report(V31.content_null(M, [V0s, W0s, xcs], 0, om_t, off=NULL_OFF, rad=NULL_RAD, log=log), 'V0')

elif mode == 'warmup':
    k = sp.Rational(1, 2)
    log('item 2 exact identities (k=1/2):', V31.exact_chart_checks(k))
    r2 = arb(2).sqrt()
    V0x = -2*r2/5; W0x = arb(8)/27; px = (9*r2/10)*(arb(7)/5)**(arb(-9)/10); xcx = -px.log()
    seed = [arb(float(V0x.mid())), arb(float(W0x.mid())), arb(float(xcx.mid()))]     # float64 guess
    log('warm-up seed (float64-rounded exact data):', [s.str(20) for s in seed])
    om_g = '0.7058823529411765'
    Pb, M = setup(k, seed[0], om_g, '-0.1')
    L = D2.krawczyk_ladder(M, seed, om_g, ladder=LADDER, log=log)
    log('PE4 ladder:', L['tried'])
    if L['ok']:
        show(L['Kx'], L['Xb']); show_repairs(L['box'], 'box')
        exact = [V0x, W0x, xcx]
        log('PE4 (K(X) contains the exact Friedmann data):', all(L['Kx'][i].contains(exact[i]) for i in range(3)))
    else:
        log('PE4: no contraction on the ladder')

elif mode == 'null2':
    log('null2 is a SMOKE CHECK, not a control (P = 0 invariant: R[0] != 0 by construction; REFEREE-E item 6)')
    Pb, M = setup(k13, V0s, om_t, '-0.05')
    lam_c = Pb.lam_c
    def centre_U(W0, nu, tails=True):
        z = 0*W0
        co = PM.coefficients(Pb.Gc, Pb.Dc, [z, W0, Pb.Ustar + z], -3*lam_c + z, [z, z, 1 + z], Pb.Kc)
        Yc = [[val(a) for a in co['Y'][i]] for i in range(3)]
        eps, epsX, info = PM.tail_bound(Pb.Gc, Pb.Dc, Yc, [val(a) for a in co['X']], val(-3*lam_c), co['DG0'], Pb.Kc, nu)
        return dict(co=co, eps=eps, epsX=epsX)
    Pb.centre_manifold = centre_U
    R, E, d = M.residual(V0s, W0s, xcs, om_t, log=log)
    log('null2 smoke: R =', [val(r).str(6) for r in R], ' certified nonzero:', any(not val(r).contains(0) for r in R))
log('total time %.1fs' % (time.time() - T0))
