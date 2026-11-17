"""B5G certified run.  usage:
  memcap 8 py -3.11 -u b5g/cert_g2.py A      (main configuration: items 1-4 + controls)
  memcap 8 py -3.11 -u b5g/cert_g2.py B      (convergence configuration: item 1 + FRW control + gauge checks)
Writes b5g/cert_<cfg>.json and prints a log.
"""
import json
import sys
import time
sys.path.insert(0, '<LAB>/b5g')
import sympy as sp
from flint import arb, ctx
import pipeline_g2 as PG
import lohner_g2 as LG
import frw_exact_g2

ctx.cap = 80
which = sys.argv[1]
T0 = time.time()
SIG0 = '[8.9013232777556435191172160567756856 +/- 5.31e-35]'      # b5/result_ec.json "sig0"
ZS = '[0.49718669611129068678126112853988237484 +/- 2.61e-39]'     # b5/result_ec.json "zs"
CF = {
    'A': dict(Nbg=30, Np=28, h0=0.02, frac=0.15, frac0=0.2, z1='0.05', zm_ec='0.4947', zm_frw='0.8635', prec=200),
    'B': dict(Nbg=36, Np=34, h0=0.016, frac=0.12, frac0=0.15, z1='0.04', zm_ec='0.4935', zm_frw='0.8620', prec=256),
}[which]


def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)


def s(x, n=25):
    return x.str(n, radius=True)


def mk_cfg(kind):
    ctx.prec = CF['prec']
    c = dict(Nbg=CF['Nbg'], Np=CF['Np'], h0=CF['h0'], frac=CF['frac'], frac0=CF['frac0'], z1=CF['z1'], prec=CF['prec'])
    if kind == 'ec':
        c.update(sig0=SIG0, zs=ZS, zm=arb(CF['zm_ec']))
    else:
        c.update(sig0=arb(3) / 8, zs=arb(3).sqrt() / 2, zm=arb(CF['zm_frw']), sonic_exact_bg=frw_exact_g2.frw_sonic_bg)
    return c


# momentum-constraint gradient (control: kappa dmu = dMo, never used by the pipeline)
mu_, Om_, V_, N_ = PG.mu_, PG.Om_, PG.V_, PG.N_
Mo = PG.S('Mo')
f_gMo = PG.mk([sp.diff(Mo, v) for v in (mu_, Om_, V_, N_)], (mu_, Om_, V_, N_), 'f_gMo')


def bg_at_zm(steps, cfg):
    st = steps[-1]
    tm = cfg['zm'] - st['z0']
    Om, V = st['y'][0].ev(tm), st['y'][1].ev(tm)
    dOm, dV = st['u'][0].ev(tm), st['u'][1].ev(tm)
    z = cfg['zm']
    N = 1 / z
    mu = PG.f_MU(Om, V, N)[0]
    SL = PG.f_SL(mu, Om)[0]
    zx = -z * SL
    Hh = PG.f_Mo(mu, Om, V, N)[1]
    return dict(Om=Om, V=V, N=N, mu=mu, SL=SL, zx=zx, dOm=dOm, dV=dV, H=Hh)


def gauge_check(P, steps, cfg):
    """kappa=1: centre-regular solution at zm must be parallel to the closed-form time-translation mode
    (mu_x, Om_x, V_x, (lnN)_x + 1)."""
    Fv, dF, lam, allc, allk = P.F(arb(1), want_vectors=True)
    b = bg_at_zm(steps, cfg)
    g = [b['H'], b['dOm'] * b['zx'], b['dV'] * b['zx'], b['SL'] + 1]
    Y = allc[0]
    ok = True
    ratios = []
    for i in range(4):
        ratios.append(Y[i] / g[i])
    for i in range(4):
        for j in range(i + 1, 4):
            cr = Y[i] * g[j] - Y[j] * g[i]
            if not cr.contains(0):
                ok = False
    ok = ok and not Y[1].contains(0)
    return ok, [s(r, 15) for r in ratios], s(Fv, 12)


def constraint_check(P, steps, cfg, kv):
    Fv, dF, lam, allc, allk = P.F(arb(kv), want_vectors=True)
    b = bg_at_zm(steps, cfg)
    gM = f_gMo(b['mu'], b['Om'], b['V'], b['N'])
    Y = allc[0]
    res = arb(kv) * Y[0] - (gM[0] * Y[0] + gM[1] * Y[1] + gM[2] * Y[2] + gM[3] * b['N'] * Y[3])
    scale = abs(arb(kv) * Y[0]) + abs(gM[1] * Y[1])
    return res.contains(0), s(res, 8), s(scale, 8)


def sgn(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def certify_sign(P, a, b, kg=1.0, w0=0.02, wg=2e-4, maxn=4000):
    """certified sign of G(kappa) = F(kappa)/(kappa - kg) on [a,b]; F(kg) = 0 exactly (gauge argument).
    On K not containing kg: F(K) in F(c) + F'(K)(K-c).  On K containing kg: G(K) in F'(K)."""
    pts = []
    x = a
    while x < b - 1e-15:
        if x < kg - wg / 2 < x + w0:
            pts.append((x, kg - wg / 2))
            pts.append((kg - wg / 2, kg + wg / 2))
            x = kg + wg / 2
            continue
        y = min(b, x + w0)
        pts.append((x, y))
        x = y
    todo = list(reversed(pts))
    done = []
    nev = 0
    signs = set()
    while todo:
        lo, hi = todo.pop()
        c = (lo + hi) / 2
        K = arb(c, (hi - lo) / 2)
        Fk, dFk, lamK = P.F(K)
        nev += 1
        if lo < kg < hi:
            sg = sgn(dFk)
            how = 'F(%g)=0 & F\'(K)' % kg
        else:
            Fc = P.F(arb(c))[0]
            nev += 1
            enc = Fc + dFk * (K - c)
            sF = sgn(enc)
            sg = sF * (1 if lo > kg else -1) if sF != 0 else 0
            how = 'mean value'
        if sg == 0:
            if nev > maxn:
                raise RuntimeError('too many evaluations')
            m = (lo + hi) / 2
            todo.append((m, hi))
            todo.append((lo, m))
            continue
        if not (lamK < 0):
            raise RuntimeError('sonic exponent not certified negative on K')
        signs.add(sg)
        done.append((lo, hi, sg, how))
    return signs, done, nev


def bisect_root(P, a, b, tol=4e-7):
    fa = P.F(arb(a))[0]
    fb = P.F(arb(b))[0]
    assert sgn(fa) * sgn(fb) == -1, (s(fa), s(fb))
    while b - a > tol:
        m = (a + b) / 2
        fm = P.F(arb(m))[0]
        if sgn(fm) == sgn(fa):
            a, fa = m, fm
        elif sgn(fm) == sgn(fb):
            b, fb = m, fm
        else:
            break
    K = arb((a + b) / 2, (b - a) / 2)
    Fa, Fb = P.F(arb(a))[0], P.F(arb(b))[0]
    FK, dFK, lamK = P.F(K)
    cert = sgn(Fa) * sgn(Fb) == -1 and (dFK > 0 or dFK < 0) and lamK < 0
    return dict(lo=a, hi=b, F_lo=s(Fa, 12), F_hi=s(Fb, 12), dF_K=s(dFK, 10), lam_K=s(lamK, 10), certified=bool(cert), K=K)


out = {'config': which, 'params': CF}
# ------------------------------------------------------------------ EC background
cfg = mk_cfg('ec')
steps, info = LG.build_background2(cfg, log=log)
log('EC background:', json.dumps(info))
out['ec_background'] = info
P = LG.Pert2(steps, cfg)
# gauge mode
ok, ratios, F1 = gauge_check(P, steps, cfg)
log('EC kappa=1: F(1)=%s  Y/g ratios %s  parallel=%s' % (F1, ratios, ok))
out['ec_gauge_eigenfunction_parallel'] = ok
out['ec_F_at_1'] = F1
out['ec_gauge_ratios'] = ratios
for kv in ('0.5', '2.0'):
    okc, res, sc = constraint_check(P, steps, cfg, kv)
    log('EC momentum-constraint residual at kappa=%s: %s (scale %s) contains0=%s' % (kv, res, sc, okc))
    out['ec_constraint_kappa_' + kv] = {'residual': res, 'scale': sc, 'contains0': okc}
# gauge zero enclosure
eps = 1e-4
Fm, Fp = P.F(arb(1 - eps))[0], P.F(arb(1 + eps))[0]
FK, dFK, lamK = P.F(arb(1, eps))
gz = sgn(Fm) * sgn(Fp) == -1 and (dFK < 0 or dFK > 0)
out['ec_gauge_zero'] = {'enclosure': '[%r, %r]' % (1 - eps, 1 + eps), 'F_lo': s(Fm, 12), 'F_hi': s(Fp, 12), 'dF_K': s(dFK, 10),
                        'certified_unique_simple_zero': bool(gz), 'contains_exact_1': True}
log('EC gauge zero:', out['ec_gauge_zero'])
# item 1
r = bisect_root(P, 2.8, 3.0)
K = r.pop('K')
ks = K
out['kappa_star'] = {'enclosure': s(K, 12), 'lo': repr(r['lo']), 'hi': repr(r['hi']), 'width': repr(r['hi'] - r['lo']), **r}
out['gamma'] = s(1 / K, 12)
log('kappa* :', out['kappa_star'], ' gamma =', out['gamma'])
if which == 'A':
    # item 3/4: certified sign of G = F/(kappa-1) on [0.40, 2.70]
    rng = (0.95, 1.10) if (len(sys.argv) > 2 and sys.argv[2] == 'probe') else (0.40, 2.70)
    signs, done, nev = certify_sign(P, *rng)
    out['deflated_sign_0p40_2p70'] = {'signs': sorted(signs), 'n_intervals': len(done), 'n_evaluations': nev,
                                     'certified_constant_sign': len(signs) == 1,
                                     'sign': {-1: 'negative', 1: 'positive'}[list(signs)[0]] if len(signs) == 1 else 'mixed'}
    json.dump([[lo, hi, sg, how] for lo, hi, sg, how in done], open('<LAB>/b5g/sign_intervals_A.json', 'w'))
    log('deflated sign on [0.40,2.70]:', out['deflated_sign_0p40_2p70'])
    # coarse scan record (how kappa* was located)
    scan = []
    for i in range(0, 21):
        kv = 0.2 * i
        Fv = P.F(arb(kv))[0]
        scan.append([round(kv, 3), s(Fv, 10)])
    out['coarse_scan_0_4'] = scan
    log('coarse scan', scan)
# ------------------------------------------------------------------ FRW control (same perturbation pipeline)
cfgF = mk_cfg('frw')
stepsF, infoF = LG.build_background2(cfgF, log=log)
log('FRW background:', json.dumps(infoF))
out['frw_background'] = infoF
PF = LG.Pert2(stepsF, cfgF)
okF, ratiosF, F1F = gauge_check(PF, stepsF, cfgF)
out['frw_gauge_eigenfunction_parallel'] = okF
out['frw_gauge_ratios'] = ratiosF
log('FRW kappa=1: F(1)=%s ratios %s parallel=%s' % (F1F, ratiosF, okF))
eps = 1e-5
Fm, Fp = PF.F(arb(1 - eps))[0], PF.F(arb(1 + eps))[0]
FK, dFK, lamK = PF.F(arb(1, eps))
fz = sgn(Fm) * sgn(Fp) == -1      # continuity on K: the interval evaluation PF.F(K) succeeded (all denominators nonzero on K)
fu = bool(dFK < 0 or dFK > 0)
out['frw_control'] = {'closed_form_rate': '1', 'enclosure': '[%r, %r]' % (1 - eps, 1 + eps), 'F_lo': s(Fm, 12), 'F_hi': s(Fp, 12),
                      'dF_K': s(dFK, 10), 'certified_sign_change': bool(fz), 'unique_zero_in_enclosure': fu, 'contains_closed_form': bool(fz and (1 - eps) < 1 < (1 + eps))}
log('FRW control:', out['frw_control'])
okc, res, sc = constraint_check(PF, stepsF, cfgF, '0.3')
out['frw_constraint_kappa_0.3'] = {'residual': res, 'scale': sc, 'contains0': okc}
log('FRW constraint residual kappa=0.3:', res, sc, okc)
near0 = []
for kv in ('-0.01', '-0.001', '0.001', '0.01'):
    Fv, dF, lam = PF.F(arb(kv))
    near0.append([kv, s(Fv, 10), s(lam, 8)])
out['frw_near_kappa0'] = near0
log('FRW near kappa=0 (sonic resonance lambda=1):', near0)
out['time_s'] = round(time.time() - T0, 1)
json.dump(out, open('<LAB>/b5g/cert_%s%s.json' % (which, '_probe' if len(sys.argv) > 2 else ''), 'w'), indent=1)
log('DONE cert', which)
