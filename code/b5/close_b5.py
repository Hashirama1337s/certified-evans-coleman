"""B5 closure checks, sealed as BAR-LB-B5C.md (seal 71904980c4), after seats/REFEREE-B5.md.
usage:  memcap 8 py -u b5/close_b5.py {c1|c2|c3|c4|c5|neg|all}

Acts on the cached certified objects of the B5 proof (result_ec.json, result_ec_centre_cache.json, cfg_ec.json,
hp_1_3_refined.json). The proof itself is not re-run. Old files are only read. Results go to close_b5_results.json,
under a cache layer whose key carries the sha256 of every producing code/data file (C4).

C1  Outer sonic stretch theta in [0, theta_o] (theta_o: P_z(theta_o) = zo = 0.4995, the outer-leg start of the proof):
    S/theta != 0, P < 0, mu < 1, V > 0.  Pieces are bisected up to depth 12; ambiguous pieces left are counted.
C2  det C != 0 in arb, for the Miranda preconditioner C rebuilt exactly as prove.py builds it. Jn is asserted
    bit-identical to the value logged in run_ec.out.
C3  P < 0 on z in (0, 0.05] (z monotone on the centre step). The centre step is not stored in the cache (only its
    end values are), so the step is re-validated here over the whole certified sig0 box. P/z is a polynomial in
    (z, sig, U) after Om = z^2 sig, V = z U, and it is evaluated on the validated enclosure.
C4  Cache integrity: every cached string fed to Miranda (G_lo, G_hi, G_mid, G_box) is re-parsed and re-printed, and the
    rebuilt ball must contain the ball it came from. Also records the sha256 of the producing files and checks that the
    cache file is newer than every producing file.
C5  FRW k = 1/2 outer leg (positive control).
    CLOSED FORM (flat FRW, p = rho/2, beta = 2/(3(1+k)) = 4/9, c = beta(1-beta) = 20/81, y = r/T < 0):
        V = beta y,  Om = (3/2) V^2,  2m/r = V^2,  z = -y/(1 + c y^2).
      - The sonic point is at y_s = -9 sqrt2/10, z_s = 9 sqrt2/14.
      - z is NOT monotone on the outer side: z has a maximum z_max = 1/(2 sqrt c) = 9/(4 sqrt5) at y = -1/sqrt c, and
        P = 0 there (z_x = 0).
      - The polar chart ends at the apparent horizon V = -1 (y = -9/4, 2m/r = 1, z = 1).
      - So FRW has NO t = 0 (N = 1/z = 0) endpoint in these coordinates. The sealed C5 wording ("reaches t = 0") cannot
        hold for FRW, and this script reports it as NOT APPLICABLE.
      - What is checked instead: the validated outer leg (the outer half of the certified sonic curve, then the
        z-integrator) from the exact FRW sonic point to z = 1 must contain the exact FRW values there.
      - z = 1 on the outer branch before z_max: c y^2 - |y| + 1 = 0 gives |y| = 9/5, so V = -4/5, Om = 24/25,
        2m/r = 16/25.
NEG C1's sign tests on a shifted interval past the next zero of S/theta or V of the certified solution curve.
    - S/theta and V have no zero on the validated sonic disk |theta| <= 0.05.
    - The literal theta shift is attempted first: locate the zero of the head polynomial, then try to validate the tail
      on a disk that reaches it. If that tail validation fails, the shift is reported as not certifiable.
    - The certified analogue is then run: the same V > 0 sign test on the certified inner solution curve (a thin
      validated run at the certified sig0 centre), on the stretch that straddles its V zero. It must report a certified
      sign change.
Doctrine #115 note: only arb_poly is used, never arb_series, so ctx.cap truncation does not apply.
"""
import sys
import os
import json
import time
import hashlib
from flint import arb, arb_poly, ctx
import sympy as sp

B5 = '<LAB>/b5/'
sys.path.insert(0, B5)
import systems as SY
import sonic as SO
from rigode import guess, validate, integrate, step_enclosure
from cheb import aup, hull

T0 = time.time()
MODE = sys.argv[1] if len(sys.argv) > 1 else 'all'
OUTF = B5 + 'close_b5_results.json'
PRODUCERS = ['prove.py', 'rigode.py', 'cheb.py', 'sonic.py', 'systems.py', 'series.py', 'system2.json',
             'sonic_sym.json', 'cfg_ec.json', 'hp_1_3_refined.json', 'close_b5.py']


def log(*a):
    print('[%7.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a)), flush=True)


def sha(p):
    return hashlib.sha256(open(B5 + p, 'rb').read()).hexdigest()


CODE_HASH = {p: sha(p) for p in PRODUCERS}
KEY = hashlib.sha256(json.dumps(CODE_HASH, sort_keys=True).encode()).hexdigest()

cfg = json.load(open(B5 + 'cfg_ec.json'))
hp = json.load(open(cfg['hp']))
RES = json.load(open(B5 + 'result_ec.json'))
SON_PREC, ODE_PREC = cfg['son_prec'], cfg['ode_prec']
NHEAD, which = cfg['nhead'], cfg['which']
DEPTH = 12


def K13():
    return arb(1) / 3


# ------------------------------------------------------------ shared sonic objects (as prove.py, SON_PREC)
def sonic_setup():
    ctx.prec = SON_PREC
    kk = K13()
    zs_c = arb(arb(hp['zs']).mid())
    wZ = arb(cfg['w_zs'])
    SCALE = arb(cfg['scale'])
    TH1 = arb(cfg['th1'])
    Bc = SO.branch(kk, zs_c, N=NHEAD, which=which, scale=SCALE)
    return dict(kk=kk, zs_c=zs_c, wZ=wZ, SCALE=SCALE, TH1=TH1, NORM=Bc['norm'])


def branch_at(S, zsb):
    kk = S['kk']
    hd = SO.branch(kk, arb(zsb.mid()), N=NHEAD, which=which, scale=S['SCALE'], norm=S['NORM'])['head']
    B = SO.branch(kk, zsb, N=NHEAD, which=which, scale=S['SCALE'], head=hd, norm=S['NORM'])
    Vt = SO.validate_tail(B, S['TH1'], pieces=16, kk=kk)
    return B, Vt


def theta_for(S, B, Vt, zt):
    th0 = float(((zt - B['e'][0]) / B['p'][0][1]).mid())
    th, ok = SO.solve_theta(B, Vt, S['kk'], zt, th0, abs(th0) * 0.3)
    assert ok, 'theta not certified'
    assert abs(th).upper() <= S['TH1'], 'theta outside validated disk'
    return th


muf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.MU)
Pf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.P)
Sf = SO.L((SY.k, SY.z, SY.Om, SY.V), SY.S)
gradS = [SO.L((SY.k, SY.z, SY.Om, SY.V), sp.diff(SY.S, x)) for x in (SY.z, SY.Om, SY.V)]


def piece_status(S, B, Vt, SH1, tb):
    """status of each condition on theta-ball tb: +1 holds, -1 certified violated, 0 ambiguous"""
    kk = S['kk']
    Pv = SO.P_eval(B, Vt, tb)
    out = {}
    mu = arb(muf(kk, *Pv))
    out['mu<1'] = 1 if mu < 1 else (-1 if mu >= 1 else 0)
    P = arb(Pf(kk, *Pv))
    out['P<0'] = 1 if P < 0 else (-1 if P >= 0 else 0)
    Vv = Pv[2]
    out['V>0'] = 1 if Vv > 0 else (-1 if Vv <= 0 else 0)
    gS = [abs(arb(g(kk, *Pv))) for g in gradS]
    rb = sum((gS[j] * Vt['rho_coef'][j] for j in range(3)), arb(0)) * abs(tb)
    s_over = SH1(tb) + arb(0, aup(rb))
    out['S/theta!=0'] = 0 if s_over.contains(0) else 1
    out['_S_over_theta'] = s_over
    out['_V'] = Vv
    return out


def bisect_check(S, B, Vt, SH1, lo, hi, npieces=8):
    """check conditions on [lo, hi] (arb endpoints); bisect ambiguous pieces up to DEPTH"""
    conds = ['mu<1', 'P<0', 'V>0', 'S/theta!=0']
    amb = {c: 0 for c in conds}
    viol = {c: 0 for c in conds}
    sign_seen = {'V': set(), 'S/theta': set()}
    stack = [(lo + (hi - lo) * i / npieces, lo + (hi - lo) * (i + 1) / npieces, 0) for i in range(npieces)]
    npc = 0
    while stack:
        a, b, d = stack.pop()
        st = piece_status(S, B, Vt, SH1, hull(a, b))
        npc += 1
        for key, nm in (('_V', 'V'), ('_S_over_theta', 'S/theta')):
            x = st[key]
            if x > 0:
                sign_seen[nm].add(1)
            elif x < 0:
                sign_seen[nm].add(-1)
        bad = [c for c in conds if st[c] == 0]
        for c in conds:
            if st[c] == -1:
                viol[c] += 1
        if bad:
            if d < DEPTH:
                m = (a + b) / 2
                stack.append((a, m, d + 1))
                stack.append((m, b, d + 1))
            else:
                for c in bad:
                    amb[c] += 1
    return {'ambiguous': amb, 'violated': viol, 'pieces_evaluated': npc,
            'sign_change_V': len(sign_seen['V']) == 2, 'sign_change_S_over_theta': len(sign_seen['S/theta']) == 2}


def sh1_poly(S, B):
    """S(H(theta))/theta as polynomial (S(e) = 0 exactly: sonic_sym.py)"""
    H = [arb_poly(B['p'][i]) for i in range(3)]
    DEG = 4 * NHEAD + 4
    Ht = [SO.TP(H[i], DEG) for i in range(3)]
    cs = Sf(S['kk'], *Ht).p.coeffs()
    assert cs[0].contains(0)
    return arb_poly(cs[1:])


# ------------------------------------------------------------ C1
def c1():
    S = sonic_setup()
    zlo, zhi = S['zs_c'] - S['wZ'], S['zs_c'] + S['wZ']
    ZS = hull(zlo, zhi)
    assert arb(RES['zs']).contains(ZS), 'stored certified zs must contain the box ZS'
    B, Vt = branch_at(S, ZS)
    zo = arb(cfg['zo'])
    tho = theta_for(S, B, Vt, zo)
    assert tho > 0, 'outer side expected at theta > 0'
    SH1 = sh1_poly(S, B)
    r = bisect_check(S, B, Vt, SH1, arb(0), tho.upper())
    hold = all(v == 0 for v in r['ambiguous'].values()) and all(v == 0 for v in r['violated'].values())
    r.update({'theta_o': tho.str(20, radius=True), 'interval': '[0, theta_o.upper]', 'HOLD': hold,
              'note': 'theta in [theta*, 0] (inner sonic piece) was already checked in prove.py'})
    log('C1', json.dumps({k: v for k, v in r.items()}, default=str))
    return r


# ------------------------------------------------------------ C2
def c2():
    S = sonic_setup()
    kk = S['kk']
    zm = arb(cfg['zm'])
    zlo, zhi = S['zs_c'] - S['wZ'], S['zs_c'] + S['wZ']

    def sonic_c(zsb):
        B, Vt = branch_at(S, zsb)
        th = theta_for(S, B, Vt, zm)
        return SO.P_eval(B, Vt, th)
    c_lo = sonic_c(arb(zlo.mid()))
    c_hi = sonic_c(arb(zhi.mid()))
    ctx.prec = ODE_PREC
    cache = json.load(open(B5 + 'result_ec_centre_cache.json'))
    G_box = [arb(x) for x in cache['G_box']]
    wZ = S['wZ']
    Jn = [[float(G_box[2].mid()), None], [float(G_box[3].mid()), None]]
    dcz = [float(((c_hi[1] - c_lo[1]) / (2 * wZ)).mid()), float(((c_hi[2] - c_lo[2]) / (2 * wZ)).mid())]
    Jn[0][1], Jn[1][1] = -dcz[0], -dcz[1]
    logged = None
    for line in open(B5 + 'run_ec.out'):
        if 'Jacobian (num) =' in line:
            logged = line.split('Jacobian (num) =')[1].split(']] ')[0] + ']]'
    identical = (logged is not None and repr(Jn) == logged.strip())
    det = Jn[0][0] * Jn[1][1] - Jn[0][1] * Jn[1][0]
    C = [[arb(Jn[1][1] / det), arb(-Jn[0][1] / det)], [arb(-Jn[1][0] / det), arb(Jn[0][0] / det)]]
    detC = C[0][0] * C[1][1] - C[0][1] * C[1][0]
    r = {'Jn': repr(Jn), 'Jn_logged_in_run_ec_out': logged, 'Jn_bit_identical_to_proof_run': identical,
         'C': [[x.str(20, radius=True) for x in row] for row in C], 'detC': detC.str(20, radius=True),
         'HOLD': bool(identical and not detC.contains(0))}
    log('C2', json.dumps(r))
    return r


# ------------------------------------------------------------ C3
def c3():
    ctx.prec = ODE_PREC
    kk = K13()
    U0 = -2 / (3 * (1 + kk))
    sig_c = arb(arb(hp['sig0']).mid())
    wS = arb(cfg['w_sig'])
    sig_box = sig_c + arb(0, aup(wS))
    assert sig_box.contains(arb(RES['sig0'])) or arb(RES['sig0']).overlaps(sig_box)
    z1 = arb(cfg['z1'])
    M = cfg['M']
    cs = SY.centre_system(False)
    ub, dif = guess(cs, kk, arb(0), z1, [sig_c, U0], M, iters=300, tol=1e-90)
    res = validate(cs, kk, arb(0), z1, [sig_box, U0], ub)
    # P/z with Om = z^2 sig, V = z U
    Pz_expr = sp.cancel(SY.P.subs({SY.Om: SY.z**2 * SY.sig, SY.V: SY.z * SY.U}) / SY.z)
    assert sp.denom(sp.together(Pz_expr)) == 1, 'P/z must be a polynomial'
    Pz = SO.L((SY.k, SY.z, SY.sig, SY.U), sp.expand(Pz_expr))
    pi = arb.pi()
    pieces = 64
    h = z1 / 2
    worst = None
    ok = True
    for p in range(pieces):
        th = hull(pi * p / pieces, pi * (p + 1) / pieces)
        zr = h + h * th.cos()
        rr = []
        for i in range(2):
            v = res['ubar'][i].c[0]
            for kk_ in range(1, len(res['ubar'][i].c)):
                v += res['ubar'][i].c[kk_] * (kk_ * th).cos()
            rr.append(v + arb(0, res['r'][i]))
        val = arb(Pz(kk, zr, rr[0], rr[1]))
        worst = val if worst is None else hull(worst, val)
        if not val < 0:
            ok = False
    r = {'P_over_z_range_on_[0,0.05]': worst.str(10, radius=True), 'centre_step_r': [x.str(5) for x in res['r']],
         'HOLD': ok, 'note': 'P = z * (P/z) with P/z < 0 on [0, 0.05], so P < 0 on (0, 0.05]; step re-validated over the sig0 box'}
    log('C3', json.dumps(r))
    return r


# ------------------------------------------------------------ C4
def c4():
    out = {'objects': {}, 'producer_sha256': CODE_HASH}
    allok = True
    for cf in ('result_ec_centre_cache.json', 'result_frw_centre_cache.json'):
        cache = json.load(open(B5 + cf))
        rows = {}
        for prec in (ODE_PREC, SON_PREC):
            ctx.prec = prec
            for nm in ('G_lo', 'G_hi', 'G_mid', 'G_box'):
                for i, s in enumerate(cache[nm]):
                    x = arb(s)
                    y = arb(x.str(radius=True))
                    c = bool(y.contains(x))
                    rows['%s[%d]@%d' % (nm, i, prec)] = c
                    allok &= c
        t_cache = os.path.getmtime(B5 + cf)
        newer = {p: t_cache > os.path.getmtime(B5 + p) for p in PRODUCERS if p not in ('close_b5.py',)
                 and not p.startswith('cfg') and not p.startswith('hp')}
        out['objects'][cf] = {'reload_contains': all(rows.values()), 'n_objects': len(rows),
                              'cache_newer_than_producers': newer,
                              'in_run_reload_test': {k: cache.get(k) for k in ('reload_contains_original', 'reload_string_identical')},
                              'cache_sha256': sha(cf)}
    ec_newer = all(out['objects']['result_ec_centre_cache.json']['cache_newer_than_producers'].values())
    out['HOLD'] = bool(allok and ec_newer)
    out['note'] = ('Miranda consumed arb(s) of these strings, so the objects fed to Miranda are the reloaded copies; '
                   'in-run original-vs-reload containment exists for G_mid only (EC run). The FRW cache predates the '
                   'last prove.py edit (see cache_newer_than_producers).')
    log('C4', json.dumps(out, default=str)[:1500])
    return out


# ------------------------------------------------------------ C5
def c5():
    ctx.prec = SON_PREC
    kk = arb(1) / 2
    zsF = 9 * arb(2).sqrt() / 14
    SCALE = arb('0.2')
    TH1 = arb('0.008')
    Bc = SO.branch(kk, zsF, N=50, which='strong', scale=SCALE)
    B = SO.branch(kk, zsF, N=50, which='strong', scale=SCALE, head=Bc['head'], norm=Bc['norm'])
    Vt = SO.validate_tail(B, TH1, pieces=16, kk=kk)
    zo = zsF + arb('0.0005')
    S = dict(kk=kk, TH1=TH1)
    th = theta_for(S, B, Vt, zo)
    Pv = SO.P_eval(B, Vt, th)
    log('C5 outer start: theta', th.str(10), 'Om', Pv[1].str(20), 'V', Pv[2].str(20))
    ctx.prec = ODE_PREC
    kk = arb(1) / 2
    rs = SY.regular_system(False)
    rec = []
    zend = arb(1)
    u = integrate(rs, kk, zo, zend, [Pv[1], Pv[2]], M=cfg['M'], h0=0.002, record=rec, gtol=cfg['gtol'], ctol=cfg['ctol'],
                  hmin=1e-9)
    exOm, exV = arb(24) / 25, -arb(4) / 5
    checks = {'mu_lt_1': True, 'P_neg': True}
    for res in rec:
        zr = hull(res['a'], res['b'])
        Omr, Vr = step_enclosure(res, 0), step_enclosure(res, 1)
        checks['mu_lt_1'] &= bool(arb(muf(kk, zr, Omr, Vr)) < 1)
        checks['P_neg'] &= bool(arb(Pf(kk, zr, Omr, Vr)) < 0)
    contain = bool(u[0].contains(exOm) and u[1].contains(exV))
    r = {'t0_endpoint': 'NOT APPLICABLE: flat FRW has no N = 0 (t = 0) point in the polar chart (z_max = 9/(4 sqrt5), '
                        'apparent horizon at z = 1, V = -1); see header',
         'checked_instead': 'outer leg from exact sonic point to z = 1, exact (Om, V) = (24/25, -4/5)',
         'Om_at_z1': u[0].str(25, radius=True), 'V_at_z1': u[1].str(25, radius=True), 'contains_exact': contain,
         'steps': len(rec), 'checks': checks, 'HOLD_as_sealed': False, 'HOLD_substitute': bool(contain and all(checks.values()))}
    log('C5', json.dumps(r))
    return r


# ------------------------------------------------------------ NEG
def neg():
    S = sonic_setup()
    zlo, zhi = S['zs_c'] - S['wZ'], S['zs_c'] + S['wZ']
    ZS = hull(zlo, zhi)
    B, Vt = branch_at(S, ZS)
    kk = S['kk']
    H = [arb_poly(B['p'][i]) for i in range(3)]
    SH1 = sh1_poly(S, B)
    # (1) literal: first zero of V(H) or S(H)/theta for |theta| <= 1.5 (head polynomial, non-rigorous locator)
    first = None
    for sgn in (-1, 1):
        prev = None
        for i in range(0, 301):
            t = arb(sgn * 1.5 * i / 300)
            vals = (float(H[2](t).mid()), float(SH1(t).mid()))
            if prev is not None and (vals[0] * prev[0] < 0 or vals[1] * prev[1] < 0):
                cand = abs(float(t.mid()))
                if first is None or cand < first[0]:
                    first = (cand, sgn, 'V' if vals[0] * prev[0] < 0 else 'S/theta')
                break
            prev = vals
    literal = {'head_zero_locator': first}
    if first is not None:
        try:
            SO.validate_tail(B, arb(first[0] * 1.05), pieces=16, kk=kk)
            literal['tail_validated_to_zero'] = True
        except Exception as ex:
            literal['tail_validated_to_zero'] = False
            literal['reason'] = str(ex)[:120]
    else:
        literal['tail_validated_to_zero'] = False
        literal['reason'] = 'no zero of V or S/theta of the head polynomial for |theta| <= 1.5'
    if literal.get('tail_validated_to_zero'):
        s = first[1]
        lo, hi = (arb(0.9 * first[0]), arb(1.1 * first[0])) if s > 0 else (arb(-1.1 * first[0]), arb(-0.9 * first[0]))
        literal['C1_on_shift'] = bisect_check(S, B, Vt, SH1, lo, hi)
    log('NEG literal', json.dumps(literal, default=str))
    # (2) certified analogue: thin validated inner solution at the certified sig0 centre, V-sign test across its zero
    ctx.prec = ODE_PREC
    kk = K13()
    U0 = -2 / (3 * (1 + kk))
    sig_c = arb(arb(hp['sig0']).mid())
    z1 = arb(cfg['z1'])
    M = cfg['M']
    cs = SY.centre_system(False)
    ub, dif = guess(cs, kk, arb(0), z1, [sig_c, U0], M, iters=300, tol=1e-90)
    res0 = validate(cs, kk, arb(0), z1, [sig_c, U0], ub)
    e = res0['end']
    u0 = [z1**2 * e[0], z1 * e[1]]
    rs = SY.regular_system(False)
    rec = []
    integrate(rs, kk, z1, arb(cfg['zm']), u0, M=M, h0=0.02, record=rec, gtol=cfg['gtol'], ctol=cfg['ctol'])
    signs = []
    for res in rec:
        Vr = step_enclosure(res, 1)
        signs.append((float(res['a'].mid()), float(res['b'].mid()), 1 if Vr > 0 else (-1 if Vr < 0 else 0),
                      1 if res['end'][1] > 0 else (-1 if res['end'][1] < 0 else 0)))
    iz = [i for i, s in enumerate(signs) if s[2] == 0]
    if iz:
        a0 = signs[max(0, iz[0] - 2)][0]
        b0 = signs[min(len(signs) - 1, iz[-1] + 2)][1]
    else:
        a0 = b0 = None
    window = [s for s in signs if a0 is not None and s[0] >= a0 and s[1] <= b0]
    neg_seen = any(s[2] == -1 for s in window)
    pos_seen = any(s[2] == 1 for s in window)
    analogue = {'window_z': [a0, b0], 'step_signs_in_window(V>0 test)': window,
                'certified_sign_change': bool(neg_seen and pos_seen),
                'V>0_test_passes_on_window': bool(all(s[2] == 1 for s in window))}
    r = {'literal_theta_shift': literal, 'certified_analogue_on_solution_curve': analogue,
         'NEG_reports_sign_change': bool(analogue['certified_sign_change'] and not analogue['V>0_test_passes_on_window'])}
    log('NEG', json.dumps(r, default=str))
    return r


MODES = {'c1': c1, 'c2': c2, 'c3': c3, 'c4': c4, 'c5': c5, 'neg': neg}
if __name__ == '__main__':
    old = {}
    if os.path.exists(OUTF):
        old = json.load(open(OUTF))
        if old.get('key') != KEY:
            old = {}
    todo = list(MODES) if MODE == 'all' else [MODE]
    out = {'key': KEY, 'code_sha256': CODE_HASH, 'results': old.get('results', {})}
    for m in todo:
        t1 = time.time()
        out['results'][m] = MODES[m]()
        out['results'][m]['seconds'] = round(time.time() - t1, 1)
        json.dump(out, open(OUTF, 'w'), indent=1, default=str)
    flags = ('HOLD', 'HOLD_as_sealed', 'HOLD_substitute', 'NEG_reports_sign_change')
    log('DONE', {m: {f: out['results'][m][f] for f in flags if f in out['results'][m]} for m in out['results']})
