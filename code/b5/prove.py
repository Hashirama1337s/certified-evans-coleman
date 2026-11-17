"""B5 blind seat: computer-assisted existence proof for the CSS radiation-fluid solution (and the FRW control).

Unknowns: sig0 = lim Om/z^2 at the centre, zs = sonic-point value of z (Maison's B).  Matching at z = zm:
  Gamma(sig0) = (Om, V)(zm) of the centre-regular solution          [rigode: Volterra (K3) + Chebyshev-Picard]
  c(zs)       = (Om, V) of the analytic sonic branch where z = zm    [sonic: parameterization + weighted contraction]
Existence: Poincare-Miranda for H = C (Gamma(sig0) - c(zs)) on the box (C = approximate inverse Jacobian).
usage: prove.py <config.json>
"""
import sys
import json
import time
from flint import arb, arb_poly, ctx
import systems as SY
import sonic as SO
from rigode import guess, validate, integrate, step_enclosure
from cheb import aup, hull

cfg = json.load(open(sys.argv[1]))
OUT = cfg['out']
LOGF = open(OUT.replace('.json', '.log'), 'w')
T0 = time.time()


def log(*a):
    s = ' '.join(str(x) for x in a)
    print('[%7.1fs] %s' % (time.time() - T0, s), flush=True)
    LOGF.write('[%7.1fs] %s\n' % (time.time() - T0, s))
    LOGF.flush()


ODE_PREC = cfg.get('ode_prec', 256)
SON_PREC = cfg.get('son_prec', 600)
ctx.prec = SON_PREC
kn, kd = [int(x) for x in cfg['k'].split('/')]


def K():
    return arb(kn) / kd


hp = json.load(open(cfg['hp']))
sig_c = arb(hp['sig0'])
zs_c = arb(hp['zs'])
sig_c = arb(sig_c.mid())
zs_c = arb(zs_c.mid())
wS = arb(cfg['w_sig'])
wZ = arb(cfg['w_zs'])
zm = arb(cfg['zm'])
z1 = arb(cfg['z1'])
which = cfg['which']
NSUB = cfg['nsub']
TH1 = arb(cfg['th1'])
NHEAD = cfg.get('nhead', 30)
SCALE = arb(cfg.get('scale', '0.2'))
result = {'config': cfg, 'resume_tested': False}

# ------------------------------------------------------------------ sonic side
ctx.prec = SON_PREC
kk = K()
Bc = SO.branch(kk, zs_c, N=NHEAD, which=which, scale=SCALE)
HEAD, NORM = Bc['head'], Bc['norm']
log('sonic head at centre zs computed; class', Bc['class']['type'], 'lam', Bc['lam'].str(12))


def sonic_c(zsb, zt, side):
    """c(zs) for zs in ball zsb: (Om, V) at z = zt on the given side (+1: S>0 inner, -1: outer)."""
    hd = SO.branch(kk, arb(zsb.mid()), N=NHEAD, which=which, scale=SCALE, norm=NORM)['head']
    B = SO.branch(kk, zsb, N=NHEAD, which=which, scale=SCALE, head=hd, norm=NORM)
    Vt = SO.validate_tail(B, TH1, pieces=16, kk=kk)
    v_z = B['p'][0][1]
    th0 = (zt - zsb) / v_z
    th0 = float(th0.mid())
    th, ok = SO.solve_theta(B, Vt, kk, zt, th0, abs(th0) * 0.3)
    if not ok:
        raise ValueError('theta not certified')
    if abs(th).upper() > TH1:
        raise ValueError('theta outside validated disk')
    Pv = SO.P_eval(B, Vt, th)
    return Pv, th, B, Vt


# centre and endpoints
t1 = time.time()
zlo, zhi = zs_c - wZ, zs_c + wZ
c_lo = sonic_c(arb(zlo.mid()), zm, 1)[0]
c_hi = sonic_c(arb(zhi.mid()), zm, 1)[0]
c_mid, th_mid, B_mid, Vt_mid = sonic_c(zs_c, zm, 1)
log('sonic endpoints done in %.1fs; theta* = %s, tail r = %s' % (time.time() - t1, th_mid.str(15),
    [x.str(3) for x in Vt_mid['r']]))
# sub-boxes covering [zlo, zhi]
subs = []
for j in range(NSUB):
    a = zlo + 2 * wZ * j / NSUB
    b = zlo + 2 * wZ * (j + 1) / NSUB
    zb = hull(a, b)
    Pv, th, B, Vt = sonic_c(zb, zm, 1)
    subs.append((zb, Pv))
log('sonic sub-boxes done (%d)' % NSUB)

# ------------------------------------------------------------------ centre side
ctx.prec = ODE_PREC
kk = K()
U0 = -2 / (3 * (1 + kk))
CSv = SY.centre_system(True)
CS0 = SY.centre_system(False)
RSv = SY.regular_system(True)
RS0 = SY.regular_system(False)
M = cfg.get('M', 28)
GT = cfg.get('gtol', 1e-62)
CT = cfg.get('ctol', 1e-64)


def run_centre(sig, var, record=None):
    cs = CSv if var else CS0
    bases = [sig, U0] + ([arb(1), arb(0)] if var else [])
    ub, dif = guess(cs, kk, arb(0), z1, bases, M, iters=300, tol=1e-90)
    res = validate(cs, kk, arb(0), z1, bases, ub)
    if record is not None:
        record.append(('centre', res))
    e = res['end']
    u0 = [z1**2 * e[0], z1 * e[1]] + ([z1**2 * e[2], z1 * e[3]] if var else [])
    rs = RSv if var else RS0
    rec = [] if record is not None else None
    u = integrate(rs, kk, z1, zm, u0, M=M, h0=cfg.get('h0', 0.02), record=rec, log=None, gtol=GT, ctol=CT)
    if record is not None:
        record.extend(('reg', r) for r in rec)
    return u


def inner_checks(REC):
    kk = K()
    Vzeros = 0
    checks = {'mu_lt_1': True, 'P_neg': True, 'V2_lt_1': True}
    muf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.MU)
    Pf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.P)
    FVf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.FV)
    sign_prev = None
    for tag, res in REC:
        a, b = res['a'], res['b']
        zr = hull(a, b)
        if tag == 'centre':
            Ur = step_enclosure(res, 1)
            sr = step_enclosure(res, 0)
            assert Ur < 0, 'U not negative on centre step'
            Vr = zr * Ur
            Omr = zr**2 * sr
            sign_prev = -1
            # mu on centre step: Om(1+kV^2+(1+k)V/z) -> use U form: V/z = U
            mur = 2 * Omr * (1 + kk * Vr**2 + (1 + kk) * Ur) / (1 - Vr**2)
            if not mur < 1:
                checks['mu_lt_1'] = False
            continue
        Omr = step_enclosure(res, 0)
        Vr = step_enclosure(res, 1)
        if not (Vr**2 < 1):
            checks['V2_lt_1'] = False
        if not arb(muf(kk, zr, Omr, Vr)) < 1:
            checks['mu_lt_1'] = False
        if not arb(Pf(kk, zr, Omr, Vr)) < 0:
            checks['P_neg'] = False
        if Vr > 0:
            if sign_prev == -1:
                raise ValueError('sign change of V not resolved inside a step')
            sign_prev = 1
        elif Vr < 0:
            if sign_prev == 1:
                raise ValueError('sign change of V not resolved inside a step')
            sign_prev = -1
        else:
            # V encloses 0 on this step: need monotone V and opposite endpoint signs
            dV = arb(FVf(kk, zr, Omr, Vr))
            va = res['ubar'][1].left() + arb(0, res['r'][1])
            vb = res['end'][1]
            if dV.contains(0) or not ((va < 0 and vb > 0) or (va > 0 and vb < 0)):
                if (va < 0 and vb < 0) or (va > 0 and vb > 0):
                    if dV.contains(0):
                        raise ValueError('V near zero, not monotone, step [%s,%s]' % (a, b))
                    # monotone and same sign at both ends: no zero
                    continue
                raise ValueError('cannot certify V zero count on step [%s,%s]' % (a.str(8), b.str(8)))
            Vzeros += 1
            sign_prev = 1 if vb > 0 else -1
    return {'Vzeros': Vzeros, 'sign_prev': sign_prev, 'checks': checks}


ctx.prec = ODE_PREC
sig_lo, sig_hi = arb((sig_c - wS).mid()), arb((sig_c + wS).mid())
sig_box = sig_c + arb(0, aup(wS))
CKEY = '|'.join(str(x) for x in (cfg['k'], sig_c.str(80), cfg['w_sig'], cfg['zm'], cfg['z1'], M, ODE_PREC, GT, CT))
CFILE = OUT.replace('.json', '_centre_cache.json')
CACHE = None
try:
    cj = json.load(open(CFILE))
    if cj.get('key') == CKEY:
        CACHE = cj
        log('centre cache found: reusing enclosures (arb strings)')
except FileNotFoundError:
    pass
if CACHE is None:
    t1 = time.time()
    G_lo = run_centre(sig_lo, False)
    G_hi = run_centre(sig_hi, False)
    G_mid = run_centre(sig_c, False)
    log('thin centre runs done in %.1fs; G_mid = %s' % (time.time() - t1, [x.str(20) for x in G_mid]))
    REC = []
    G_box = run_centre(sig_box, True, record=REC)
    log('box centre run (variational) done: dG/dsig = %s ; steps %d' % ([x.str(10) for x in G_box[2:]], len(REC)))
    inner = inner_checks(REC)
    CACHE = {'key': CKEY, 'G_lo': [x.str(radius=True) for x in G_lo], 'G_hi': [x.str(radius=True) for x in G_hi],
             'G_mid': [x.str(radius=True) for x in G_mid], 'G_box': [x.str(radius=True) for x in G_box],
             'inner': inner, 'nsteps': len(REC)}
    json.dump(CACHE, open(CFILE, 'w'), indent=1)
    # resume test: reload and compare
    back = json.load(open(CFILE))
    same = all(arb(a).str(radius=True) == b for nm_ in ('G_lo', 'G_hi', 'G_mid', 'G_box')
               for a, b in zip(back[nm_], CACHE[nm_]))
    contains = all(arb(a).contains(x) for a, x in zip(back['G_mid'], G_mid))
    CACHE['reload_contains_original'] = bool(contains)
    CACHE['reload_string_identical'] = bool(same)
    json.dump(CACHE, open(CFILE, 'w'), indent=1)
    log('cache write+reload test: strings identical', same, '; rebuilt balls contain originals', contains)
result['resume_tested'] = bool(CACHE.get('reload_contains_original', False)) or 'unknown (cache written by earlier run; see its log)'
G_lo = [arb(x) for x in CACHE['G_lo']]
G_hi = [arb(x) for x in CACHE['G_hi']]
G_mid = [arb(x) for x in CACHE['G_mid']]
G_box = [arb(x) for x in CACHE['G_box']]
log('G radii: lo %s mid %s box %s' % ([float(x.rad()) for x in G_lo], [float(x.rad()) for x in G_mid], [float(x.rad()) for x in G_box]))

# ------------------------------------------------------------------ Miranda
ctx.prec = ODE_PREC
Jn = [[float(G_box[2].mid()), None], [float(G_box[3].mid()), None]]
# numerical dc/dzs from thin endpoints (preconditioner only)
dcz = [float(((c_hi[1] - c_lo[1]) / (2 * wZ)).mid()), float(((c_hi[2] - c_lo[2]) / (2 * wZ)).mid())]
Jn[0][1], Jn[1][1] = -dcz[0], -dcz[1]
det = Jn[0][0] * Jn[1][1] - Jn[0][1] * Jn[1][0]
C = [[arb(Jn[1][1] / det), arb(-Jn[0][1] / det)], [arb(-Jn[1][0] / det), arb(Jn[0][0] / det)]]
log('Jacobian (num) =', Jn, 'cond-ish', abs(det))


def H(G, c):
    d0, d1 = G[0] - c[1], G[1] - c[2]
    return [C[0][0] * d0 + C[0][1] * d1, C[1][0] * d0 + C[1][1] * d1]


# faces sig = lo/hi: H_1 over all zs sub-boxes
f1_lo = [H(G_lo, Pv)[0] for (zb, Pv) in subs]
f1_hi = [H(G_hi, Pv)[0] for (zb, Pv) in subs]
# faces zs = lo/hi: H_2 with G over the sig box (mean value form)
# H_2(sig, zs+-) = C2 (G_mid - c(zs+-)) + (C2 . G'(xi)) (sig - sig_c), xi in box  (mean value; scalar first)
C2Gp = C[1][0] * G_box[2] + C[1][1] * G_box[3]
f2_lo = H(G_mid, c_lo)[1] + C2Gp * arb(0, aup(wS))
f2_hi = H(G_mid, c_hi)[1] + C2Gp * arb(0, aup(wS))
log('C2.dG/dsig =', C2Gp.str(5))


def sgn(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


s1lo = set(sgn(x) for x in f1_lo)
s1hi = set(sgn(x) for x in f1_hi)
s2lo, s2hi = sgn(f2_lo), sgn(f2_hi)
mir = (len(s1lo) == 1 and len(s1hi) == 1 and 0 not in s1lo and 0 not in s1hi and s1lo != s1hi
       and s2lo != 0 and s2hi != 0 and s2lo != s2hi)
log('Miranda: H1 signs on sig_lo face', s1lo, 'sig_hi face', s1hi, '| H2 sign zs_lo', s2lo, 'zs_hi', s2hi, '=>', mir)
log('  H1(lo) worst', max(float(x.upper().mid()) if sgn(x) < 0 else -float(x.lower().mid()) for x in f1_lo),
    'H2 lo/hi', f2_lo.str(5), f2_hi.str(5))
result['miranda'] = mir
result['H2_lo'] = f2_lo.str(10)
result['H2_hi'] = f2_hi.str(10)
result['sig0_box'] = [sig_lo.str(70), sig_hi.str(70)]
result['zs_box'] = [arb(zlo.mid()).str(70), arb(zhi.mid()).str(70)]
json.dump(result, open(OUT, 'w'), indent=1)
if not mir:
    log('MIRANDA FAILED')
    sys.exit(2)

# ------------------------------------------------------------------ certified consequences
ctx.prec = SON_PREC
kk = K()
q = kk.sqrt()
ZS = hull(zlo, zhi)
SIG = hull(sig_c - wS, sig_c + wS)
cl, e, Jz, t = SO.classify(kk, ZS)
tr, c2, disc = cl['tr'], cl['c2'], cl['disc']
sq = disc.sqrt()
lp, lm = (tr + sq) / 2, (tr - sq) / 2
Vs = (q - ZS) / (1 - q * ZS)
Oms = e[1]
mus = arb(SO.L((SY.k, SY.z, SY.Om, SY.V), SY.MU)(kk, ZS, Oms, Vs))
result.update({'sig0': SIG.str(40, radius=True), 'zs': ZS.str(40, radius=True), 'V_s': Vs.str(40, radius=True),
               'Om_s': Oms.str(40, radius=True), 'mu_s': mus.str(40, radius=True),
               'sonic_type': cl['type'], 'sonic_eig_pos': lp.str(30, radius=True), 'sonic_eig_neg_or_second': lm.str(30, radius=True),
               'sonic_c2': c2.str(20, radius=True), 'sonic_disc': disc.str(20, radius=True)})
log('certified: sig0 in', SIG.str(30), ' zs in', ZS.str(30), ' V_s in', Vs.str(30))
log('sonic point type', cl['type'], 'eigenvalues', lp.str(20), lm.str(20), '(third = 0 exactly)')
json.dump(result, open(OUT, 'w'), indent=1)

muf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.MU)
Pf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.P)
FVf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.FV)
IC = CACHE['inner']
Vzeros, sign_prev, checks = IC['Vzeros'], IC['sign_prev'], IC['checks']
log('inner part (centre -> zm): V zeros =', Vzeros, 'final sign', sign_prev, 'checks', checks)
# sonic piece theta in [theta*, 0]: V, S/theta, P, mu
ctx.prec = SON_PREC
kk = K()
Sf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.S)
son_ok = True
sonV = None
for (zb, Pv) in subs[:1] + subs[-1:]:
    B = SO.branch(kk, ZS, N=NHEAD, which=which, scale=SCALE, head=HEAD, norm=NORM)
    Vt = SO.validate_tail(B, TH1, pieces=16, kk=kk)
    th, ok = SO.solve_theta(B, Vt, kk, zm, float(th_mid.mid()), abs(float(th_mid.mid())) * 0.3)
    npc = 64
    lo = th.lower() if th < 0 else arb(0)
    hi = arb(0) if th < 0 else th.upper()
    H = [arb_poly(B['p'][i]) for i in range(3)]
    for pc in range(npc):
        tb = hull(lo + (hi - lo) * pc / npc, lo + (hi - lo) * (pc + 1) / npc)
        Pv2 = SO.P_eval(B, Vt, tb)
        Vb = Pv2[2]
        sonV = Vb if sonV is None else hull(sonV, Vb)
        if not (arb(muf(kk, *Pv2)) < 1 and arb(Pf(kk, *Pv2)) < 0):
            son_ok = False
    # S(P(th))/th nonzero on [lo,hi]: S(H)/th polynomial + |grad S| |rho|/|th|
    DEG = 4 * NHEAD + 4
    Ht = [SO.TP(H[i], DEG) for i in range(3)]
    SH = Sf(kk, *Ht).p
    cs_ = SH.coeffs()
    assert cs_[0].contains(0)
    SH1 = arb_poly([arb(0)] + cs_[1:]) if False else arb_poly(cs_[1:])
    tb = hull(lo, hi)
    Pbox = SO.P_eval(B, Vt, tb)
    import sympy as sp
    gradS = [SO.L((SY.k, SY.z, SY.Om, SY.V), sp.diff(SY.S, x)) for x in (SY.z, SY.Om, SY.V)]
    gS = [abs(arb(g(kk, *Pbox))) for g in gradS]
    rb = sum((gS[j] * Vt['rho_coef'][j] for j in range(3)), arb(0)) * abs(tb)
    s_over = SH1(tb) + arb(0, aup(rb))
    if s_over.contains(0):
        son_ok = False
    result['S_over_theta_on_sonic_piece'] = s_over.str(10)
log('sonic piece: V range', sonV.str(10), 'checks ok', son_ok)
if sonV.contains(0):
    log('WARNING: V may vanish on sonic piece')
inner_zeros = Vzeros + (0 if (sonV > 0 or sonV < 0) else None)
result.update({'V_zeros_inner': inner_zeros, 'inner_checks': checks, 'sonic_piece_ok': son_ok,
               'V_on_sonic_piece': sonV.str(10)})
json.dump(result, open(OUT, 'w'), indent=1)

# ---- outer part: from the other side of the sonic point to t = 0 (N = 1/z = 0)
if cfg.get('outer', True):
    ctx.prec = SON_PREC
    kk = K()
    zo = arb(cfg['zo'])
    outs = None
    for (zb, _) in subs:
        Pv, th, B, Vt = sonic_c(zb, zo, -1)
        outs = Pv if outs is None else [hull(outs[i], Pv[i]) for i in range(3)]
    log('outer start at z=%s: Om %s V %s (theta %s)' % (zo.str(8), outs[1].str(25), outs[2].str(25), th.str(8)))
    ctx.prec = ODE_PREC
    kk = K()
    rec = []
    Zsw = arb(cfg.get('zswitch', '3'))
    u = integrate(RS0, kk, zo, Zsw, [outs[1], outs[2]], M=M, h0=0.005, record=rec, gtol=GT, ctol=CT)
    oc = {'mu_lt_1': True, 'P_neg': True, 'V_pos': True}
    for res in rec:
        zr = hull(res['a'], res['b'])
        Omr, Vr = step_enclosure(res, 0), step_enclosure(res, 1)
        oc['mu_lt_1'] &= bool(arb(muf(kk, zr, Omr, Vr)) < 1)
        oc['P_neg'] &= bool(arb(Pf(kk, zr, Omr, Vr)) < 0)
        oc['V_pos'] &= bool(Vr > 0)
    log('outer z-part to z=%s: Om %s V %s checks %s' % (Zsw.str(5), u[0].str(25), u[1].str(25), oc))
    FS = SY.far_system()
    rec2 = []
    Nsw = 1 / Zsw
    # integrate in N from Nsw down to 0: use variable s = Nsw - N >= 0 -> reuse integrator with sign flip
    import sympy as sp
    sN = sp.symbols('sN')
    from rigode import System
    # variable s = -N runs from -1/zswitch up to 0 (t = 0); dOm/ds = -dOm/dN at N = -s (integer coefficients)
    G = [(-g).subs(SY.N, -sN) for g in [SY.E['dOm_dN'], SY.E['dV_dN']]]
    G = [sp.together(g) for g in G]
    FSs = System(sN, [SY.Om, SY.V], SY.k, G, ['int', 'int'])
    u2 = integrate(FSs, kk, -1 / Zsw, arb(0), [u[0], u[1]], M=M, h0=0.02, record=rec2, gtol=GT, ctol=CT)
    for res in rec2:
        Omr, Vr = step_enclosure(res, 0), step_enclosure(res, 1)
        sr = hull(res['a'], res['b'])
        Nr = -sr
        mu_r = 2 * Omr * (1 + kk * Vr**2 + (1 + kk) * Vr * Nr) / (1 - Vr**2)
        oc['mu_lt_1'] &= bool(mu_r < 1)
        oc['V_pos'] &= bool(Vr > 0)
    Om0, V0 = u2
    mu0 = 2 * Om0 * (1 + kk * V0**2) / (1 - V0**2)
    result.update({'t0_Om': Om0.str(30, radius=True), 't0_V': V0.str(30, radius=True),
                   't0_2m_over_r': mu0.str(30, radius=True), 't0_m_over_r': (mu0 / 2).str(30, radius=True),
                   't0_a': (1 / (1 - mu0).sqrt()).str(30, radius=True), 'outer_checks': oc})
    log('t=0: m/r =', (mu0 / 2).str(25), ' Om =', Om0.str(25), ' V =', V0.str(25), ' a =', (1 / (1 - mu0).sqrt()).str(20), oc)
json.dump(result, open(OUT, 'w'), indent=1)
log('DONE')
