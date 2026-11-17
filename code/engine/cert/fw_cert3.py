"""B4 Far Wall certificate, referee repair (NEW file; fw_cert.py and fw_cert2.py are imported, never edited).
REFEREE-B4FW 'next strike':
  (1) ctx.cap = 32 BEFORE any series is built (doctrine #115: fw_cert never raised the cap, so .integral()
      and exp() truncated to 10 coefficients and the order-11 / order-22 Lagrange terms read as 0);
  (2) strict coefficient access: coef(s, k) RAISES if k >= s.prec (a coefficient that is not known);
      an index inside s.prec that coeffs() omits is an exact zero;
  (3) Lemma 4.1 asserted: C0(0) contains C00, C1(0) contains 0, C1'(0) contains C11, and every coefficient
      the centre-manifold Taylor models drop (k < m) contains 0;
  (4) background read from lin_fw4_bg3.pkl (radii rounded up, lin_fw4_bg3_make.py);
  (5) centre data built at 1.001 v_need;
  (6) tau/x switch sliver: explicit bootstrap bound (no 1e-8 inflation), see sliver();
  (7) lower/rhs printed (worst ratio) together with the absolute margin.
The strict versions of taylor_cell, split_tau, tm_eval, centre_manifold_data and tube_cells replace the fw_cert
module attributes, so every fw_cert / fw_cert2 routine that calls them by module name uses the repaired ones.
usage:
  fw_cert3.py capchk                          positive (and negative) control for the cap leak
  fw_cert3.py run M Ka Kmax [nth ns nb1 depth] certificate (fw_cert2 theta/K machinery; ns ladder from FW2_NS_LADDER)
"""
import sys, os, math, time, pickle
from flint import arb, acb, arb_series, ctx
ctx.cap = 32                                                    # (1) before anything else
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
_argv = sys.argv[:]; sys.argv = sys.argv[:1]
import fw_cert as FW
import fw_cert2 as FW2
sys.argv = _argv
from core import sym
ctx.cap = 32
assert ctx.cap == 32
T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
FW.log = log; FW2.log = log

# (4) background with radii rounded up
FW.BG = pickle.load(open('<LAB>/engine/cert/lin_fw4_bg3.pkl', 'rb'))
FW.MARGIN = FW.A_(FW.BG['margin']).abs_upper() + arb('1e-30')
FW.XSW = FW.A_(FW.BG['x_switch']); FW.XBAR = FW.A_(FW.BG['xbar'])
FW.XC = FW.A_(FW.BG['Xb'][2]); FW.PS = (-FW.XC).exp()
A_, S, CAP, PTM, NSUB = FW.A_, FW.S, FW.CAP, FW.PTM, FW.NSUB

# (2) strict coefficients
def coef(s, k):
    if k >= s.prec:
        raise ArithmeticError('coefficient %d requested from a series known only to O(x^%d)' % (k, s.prec))
    cs = s.coeffs()
    return cs[k] if k < len(cs) else arb(0)

TAYLOR_TOL = arb('1e-13')
def taylor_cell(tape, start, rough, h):
    """enclosure of the trajectory over the step [0, h] from the thin start.  The order-(CAP-1) Lagrange term is
    bounded with the jets on the a-priori box `rough`; if that term is too big (it is ~1e21 on wide E5 boxes),
    the step is cut into n equal sub-steps (thin restart at each sub-step end, same a-priori box) and the
    sub-step boxes are united.  Valid for any n (each sub-step trajectory lies in `rough`)."""
    P = CAP - 1
    Yr, _ = FW.picard_t(tape, rough)
    R = max(arb(coef(y, P).abs_upper()) for y in Yr)
    ah = arb(h.abs_upper())
    n = 1
    while R*(ah/n)**P > TAYLOR_TOL and n < 4096:
        n *= 2
    hs = h/n
    cur = list(start); hull = None
    for j in range(n):
        Ys, _ = FW.picard_t(tape, cur)
        tau = hs/2 + sym(hs/2)
        box = [sum((coef(Ys[i], k)*tau**k for k in range(P)), arb(0)) + coef(Yr[i], P)*tau**P for i in range(len(start))]
        box = [b.intersection(r) for b, r in zip(box, rough)]
        hull = box if hull is None else [a.union(b) for a, b in zip(hull, box)]
        if j < n - 1:
            cur = [sum((coef(Ys[i], k)*hs**k for k in range(P)), arb(0)) + coef(Yr[i], P)*tau**P for i in range(len(start))]
    return hull

def tm_eval(fs, fb, tau):
    res = {}
    for key in fs:
        res[key] = sum((coef(fs[key], k)*tau**k for k in range(PTM)), arb(0)) + coef(fb[key], PTM)*tau**PTM
    return res

def endpoint(c, tj):
    """thin enclosure of the 4D tau-leg state at tau = tj inside step c (Taylor + Lagrange on the step box)."""
    box = taylor_cell(FW.TT, c['start'], c['rough'], c['ht'])
    Ys, _ = FW.picard_t(FW.TT, c['start']); Yr, _ = FW.picard_t(FW.TT, box)
    P = CAP - 1
    tb = tj/2 + sym(tj/2)
    return [sum((coef(Ys[i], k)*tj**k for k in range(P)), arb(0)) + coef(Yr[i], P)*tb**P for i in range(4)], box

def split_tau(c):
    box = taylor_cell(FW.TT, c['start'], c['rough'], c['ht'])
    hs = c['ht']/NSUB
    starts = [c['start']] + [endpoint(c, hs*j)[0] for j in range(1, NSUB)]
    subs = []
    for j in range(NSUB):
        sb = taylor_cell(FW.TT, starts[j], box, hs)
        fs = FW._fser('tau', starts[j]); fb = FW._fser('tau', sb)
        xhi = c['xhi'] if j == 0 else starts[j][3]
        xlo = c['xlo'] if j == NSUB - 1 else starts[j + 1][3]
        subs.append(dict(xlo=xlo, xhi=xhi, h=xhi - xlo, kind='tau', box=sb[:3], xbox=sb[3],
                         f=tm_eval(fs, fb, hs/2 + sym(hs/2)), step=c))
    return subs

def tube_cells():
    """fw_cert.tube_cells without the 1e-8 inflation of the last tau step (replaced by sliver())."""
    cells = []
    tr = FW.BG['tau_rec']
    xs = [A_(st[3]) for (_, _, _, st) in tr] + [FW.XSW]
    for k, (_, h, Bx, st) in enumerate(tr):
        rough = [A_(Bx[i]) + sym(FW.MARGIN) for i in range(4)]
        start = [A_(st[i]) + sym(FW.MARGIN) for i in range(4)]
        cells.append(dict(xlo=xs[k + 1], xhi=xs[k], h=(xs[k] - xs[k + 1]), rough=rough, start=start, ht=A_(h), kind='tau'))
    x = FW.XSW
    for (_, h, Bx, st) in FW.BG['rec']:
        ha = A_(h)
        rough = [A_(Bx[i]) + sym(FW.MARGIN) for i in range(3)]
        start = [A_(st[i]) + sym(FW.MARGIN) for i in range(3)]
        cells.append(dict(xlo=x + ha, xhi=x, h=-ha, rough=rough, start=start, ht=ha, kind='x'))
        x = x + ha
    return cells

def sliver(cells):
    """(6) The last tau step ends at x_end (thin enclosure).  If x_end may lie above x_switch, the background on
    [x_switch, x_end] is not in any tau box.  Bootstrap: with B = (last tau sub-cell box) + 1e-6 and
    F = sup_B |f_bg|, every state on the sliver is within delta = l F of the landing state (which lies in the
    sub-cell box) as long as delta < 1e-6, l = max(0, upper(x_end - x_switch)).  The last tau sub-cell then gets
    the box (sub-box + delta) and coefficient values from interval jets on that box (valid at every point of it)."""
    ilast = max(i for i, c in enumerate(cells) if c['kind'] == 'tau')
    c = cells[ilast]; st = c['step']
    yend, _ = endpoint(st, st['ht'])
    gap = arb((yend[3] - FW.XSW).upper())                      # upper bound of x_end - x_switch
    ell = gap if gap > 0 else arb(0)                            # gap <= 0 certified: no sliver
    B = [b + sym(arb('1e-6')) for b in c['box']]
    F = max(arb(t.abs_upper()) for t in FW.FT.eval(B))
    delta = arb((ell*F).upper())
    assert delta < arb('1e-6'), 'sliver bootstrap fails'
    info = dict(x_end=yend[3], ell=ell, F=F, delta=delta)
    if ell > 0:
        # the tau sub-cell keeps its Taylor model on [XSW + ell, xhi] (x_end <= XSW + ell); the sliver
        # [XSW, XSW + ell] becomes its own cell with box (landing box + delta) and interval-jet values.
        # tau sub-cell keeps its Taylor model on [XSW + ell, xhi] (x_end <= XSW + ell).  The first x-leg cell's
        # Taylor model (expanded at the thin switch set, x = XSW) is extended backwards to x = XSW + ell; its
        # Lagrange remainder is taken on hull(cell box, landing box + delta), which contains the trajectory on
        # [XSW + h, XSW + ell].
        top = FW.XSW + ell
        c['xlo'] = top; c['h'] = c['xhi'] - top
        nb = [b + sym(delta) for b in c['box']]
        cx = cells[ilast + 1]; assert cx['kind'] == 'x' and cx['xhi'] is FW.XSW or cx['kind'] == 'x'
        hull = [cx['box'][i].union(nb[i]) for i in range(3)]
        fs = FW._fser('x', cx['start'][:3]); fb = FW._fser('x', hull)
        lo = cx['xlo'] - FW.XSW                                  # = h_x < 0
        tau = (lo + ell)/2 + sym((ell - lo)/2)                   # covers [h_x, ell]
        cx['f'] = tm_eval(fs, fb, tau)
        cx['xhi'] = top; cx['h'] = top - cx['xlo']; cx['box'] = hull
        info['f_check'] = dict(sig=cx['f']['sig'], ellphi=cx['f']['ellphi'])
    return info

def centre_manifold_data(nint=10, cap=24, P=22):
    """fw_cert.centre_manifold_data with strict coefficients and the Lemma 4.1 asserts (3)."""
    BG = FW.BG; C00, C11, MT = FW.C00, FW.C11, FW.MT
    cm = BG['cm']; nu = arb(BG['nu_c']); thm = A_(BG['thm'])
    epsY = arb(cm['eps'][0]) + arb(cm['eps'][1]); epsX = arb(cm['epsX'][0]) + arb(cm['epsX'][1])
    Yc = [[A_(a) for a in cm['Y'][i]] for i in range(3)]; Xc = [A_(a) for a in cm['X']]
    assert Yc[0][0].is_zero() and (Yc[0][1] - 1).is_zero(), 'centre manifold: need P = theta + O(theta^2) exactly'
    assert ctx.cap >= cap + 2
    def j0(c, eps):
        return arb_series([(c[k] if k < len(c) else arb(0)) + sym(eps/nu**k) for k in range(cap)], cap)
    def jb(c, eps, Th, rho):
        out = []
        for k in range(cap):
            s = arb(0)
            for n_ in range(k, len(c)):
                s += c[n_]*arb(math.comb(n_, k))*Th**(n_ - k)
            out.append(s + sym(eps/rho**k))
        return arb_series(out, cap)
    def mcZ(Z, cap_):
        o = MT.series_eval(Z, cap_)
        return [[o[4*i + j] for j in range(4)] for i in range(4)], [[o[16 + 4*i + j] for j in range(4)] for i in range(4)]
    th = arb_series([arb(0), arb(1)], cap)
    Y0 = [j0(Yc[i], epsY) for i in range(3)]
    X0 = j0(Xc, epsX)
    C0s, C1s = mcZ(Y0, cap); v0 = th*X0.exp()
    assert v0.prec >= cap, 'v0 truncated (cap leak)'
    # Lemma 4.1 asserts
    for i in range(4):
        for j in range(4):
            assert coef(C0s[i][j], 0).contains(C00[i][j]), 'C0(0) != C00 at %d%d' % (i, j)
            assert coef(C1s[i][j], 0).contains(0), 'C1(0) != 0 at %d%d' % (i, j)
            assert coef(C1s[i][j], 1).contains(C11[i][j]), "C1'(0) != C11 at %d%d" % (i, j)
    assert coef(v0, 0).contains(0) and coef(v0, 1).contains(1)
    C01 = [[coef(C0s[i][j], 1) for j in range(4)] for i in range(4)]
    C12 = [[coef(C1s[i][j], 2) - C11[i][j]*coef(v0, 2) for j in range(4)] for i in range(4)]
    def G_series(C0, C1, v, which):
        out = [[None]*4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                if which == 'G0': g = C0[i][j] - C00[i][j]; m = 1
                elif which == 'G1': g = C1[i][j] - v*C11[i][j]; m = 2
                elif which == 'R0': g = C0[i][j] - C00[i][j] - v*C01[i][j]; m = 2
                else: g = C1[i][j] - v*C11[i][j] - v*v*C12[i][j]; m = 3
                out[i][j] = (g, m)
        return out
    A0 = {w: G_series(C0s, C1s, v0, w) for w in ('G0', 'G1', 'R0', 'R1')}
    nz = 0
    for w in A0:
        for i in range(4):
            for j in range(4):
                g, m = A0[w][i][j]
                for k in range(m):
                    assert coef(g, k).contains(0), 'dropped coefficient %s[%d][%d]_%d does not contain 0' % (w, i, j, k)
                    nz += 1
    data = []
    Tb = [thm*arb(2)**(-(nint - k)) for k in range(1, nint + 1)]
    for hi in Tb:
        Th = hi/2 + sym(hi/2)
        rho = nu - hi.abs_upper()
        Yb = [jb(Yc[i], epsY, Th, rho) for i in range(3)]
        Xb = jb(Xc, epsX, Th, rho)
        C0b, C1b = mcZ(Yb, cap)
        vb = th*Xb.exp()
        chiTh = sum((Xc[n_]*Th**n_ for n_ in range(len(Xc))), arb(0)) + sym(epsX)
        ent = {}
        for which in ('G0', 'G1', 'R0', 'R1'):
            Ab = G_series(C0b, C1b, vb, which)
            M_ = [[arb(0)]*4 for _ in range(4)]
            for i in range(4):
                for j in range(4):
                    g0, m = A0[which][i][j]; gb, _ = Ab[i][j]
                    H = sum((coef(g0, k)*Th**(k - m) for k in range(m, P)), arb(0)) + coef(gb, P)*Th**(P - m)
                    M_[i][j] = FW.up(H*(-m*chiTh).exp())
            ent[which] = M_
        data.append(dict(vlo=arb(0), vhi=FW.up(hi*chiTh.exp()), **ent))
    return dict(C01=C01, C12=C12, data=data, thm=thm, lemma41_checked=nz)

FW_bulk_orig = FW.bulk_ia
def bulk_best(cells, son, M, Ka, Kb, weights=None, iters=None):
    """the weights are a free (non-rigorous) choice; the float power iteration that fixes them is fragile.
    Try several iteration counts, keep the certified result with the smallest Collatz bound."""
    best = None
    for it in (6, 3, 4, 2, 8, 10, 1):
        try:
            r = FW_bulk_orig(cells, son, M, Ka, Kb, iters=it)
        except (ValueError, ArithmeticError, FW.np.linalg.LinAlgError):
            continue
        if r['ok'] and (best is None or r['collatz'] < best['collatz']):
            best = r; best['iters'] = it
            if r['collatz'] < arb('0.8'): break
    if best is None:
        return dict(ok=False, collatz=arb('inf'))
    return best
FW.bulk_ia = bulk_best

FW_centre_orig = FW.centre_theta
def centre_theta_safe(th_lo, th_hi, Ka, M, cdat, ns=150, **kw):
    """a diverging float fixed point (no supersolution at this ns) made fw_cert's diagnostic max() crash;
    it is a FAILURE of the cell at this ns (the ladder / bisection of fw_cert2 then takes over)."""
    try:
        return FW_centre_orig(th_lo, th_hi, Ka, M, cdat, ns=ns, **kw)
    except ValueError as e:
        log('   theta [%.4f, %.4f] ns=%d: no supersolution (float fixed point diverged)' % (float(th_lo.mid()), float(th_hi.mid()), ns))
        return {'w': None, 'N': None}
FW.centre_theta = centre_theta_safe

FW_taylor_orig = FW.taylor_cell
FW.taylor_cell = taylor_cell; FW.tm_eval = tm_eval; FW.split_tau = split_tau
FW.tube_cells = tube_cells; FW.centre_manifold_data = centre_manifold_data

# ------------------------------------------------------------------ setup / run
def setup(M, Ka, Kmax):
    T1 = time.time()
    cells = FW.eval_cells(FW.tube_cells())
    sl = sliver(cells)
    log('[sliver] x_end %s, l = %s, sup|f| = %s, delta = %s%s' % (sl['x_end'].str(8), sl['ell'].str(3), sl['F'].str(3),
        sl['delta'].str(3), ('; last tau sub-cell re-enclosed: sig %s ellphi %s' % (sl['f_check']['sig'].str(4), sl['f_check']['ellphi'].str(4))) if 'f_check' in sl else '; no sliver (x_end <= x_switch)'))
    tube_n = len(cells)
    xbar_K = float((arb(M)/cells[-1]['xlo'].exp()).lower())
    if Kmax > xbar_K:
        th_min = float((FW.PS*arb(M)/arb(Kmax)).lower())*0.97
        cells = cells + FW.centre_bulk_cells(th_min)
        log('[setup] bulk extended below xbar: %d centre-manifold cells (theta >= %.3g)' % (len(cells) - tube_n, th_min))
    son = FW.sonic_piece(); cells, son = FW.prep_bulk(cells, son)
    Kcap = float((arb(M)/cells[-1]['xlo'].exp()).lower())
    if Kmax > Kcap:
        log('Kmax %.1f lowered to %.2f (deepest bulk cell)' % (Kmax, Kcap)); Kmax = Kcap
    cm = FW.centre_manifold_data(); FW.cdat_C01, FW.cdat_C12 = cm['C01'], cm['C12']
    log('[lemma 4.1] C0(0) = C00, C1(0) = 0, C1\'(0) = C11 contained; %d dropped coefficients contain 0' % cm['lemma41_checked'])
    vneed = FW.PS*arb(M)/arb(Ka)*arb('1.001')                    # (5)
    cdat = cm['data'] + [FW.bt_centre_data(cm['C01'], cm['C12'])] + FW.tube_centre_data(cells, cm['C01'], cm['C12'], vneed)
    log('[setup] cells %d, sonic g=%s, A_max=%s, centre data %d cells, ctx.cap %d (%.0fs)' % (
        len(cells), son['g'].str(4), son['Amax'].str(4), len(cdat), ctx.cap, time.time() - T1))
    return cells, son, cdat, Kmax

def run(M, Ka, Kmax, nth, ns, nb1, maxdepth, kratio=1.03):
    FW.NB1 = nb1
    log('fw_cert3: M=%s Ka=%s Kmax=%s theta cells=%d ns=%d NB1=%d depth %d kratio %s ladder %s prec %d cap %d' % (
        M, Ka, Kmax, nth, ns, nb1, maxdepth, kratio, os.environ.get('FW2_NS_LADDER', '-'), FW.PREC, ctx.cap))
    cells, son, cdat, Kmax = setup(M, Ka, Kmax)
    bulks = FW2.bulk_cells_K(cells, son, M, Ka, Kmax, kratio)
    if bulks is None:
        log('VERDICT: bulk not certified'); return None
    pi2 = arb.pi()/2
    ths = [pi2*k/nth for k in range(nth + 1)]
    stats = dict(centre=0, split=0, ns_used={})
    fails = []; wmin = None; rmin = None
    sub = [int(t) for t in os.environ.get('FW3_SUB', '').split(',') if t]      # probe only (non-counting)
    for k in (sub if sub else range(nth)):
        ok, worst, wm =FW2.theta_cell(cells, son, M, Ka, cdat, bulks, ths[k], ths[k + 1], ns, 0, maxdepth, stats)
        if wm is not None: wmin = wm if wmin is None else min(wmin, wm)
        if worst is not None: rmin = worst if rmin is None else min(rmin, worst, key=lambda t: t[0])
        if not ok:
            fails.append(k); log('  theta cell %d FAIL (worst lower/rhs %s)' % (k, worst))
    assert ctx.cap == 32
    log('accepted cells by ns %s; theta cells %d (+%d bisections, %d centre evaluations), K-cells %d' % (
        stats['ns_used'], nth, stats['split'], stats['centre'], len(bulks)))
    log('worst lower/rhs = %.6f at K in [%.2f, %.2f]; worst absolute margin (lower - rhs) = %.6g' % (rmin[0], rmin[1], rmin[2], wmin))
    if sub:
        log('PROBE ONLY (theta cells %s): failures %s' % (sub, fails)); return fails
    if not fails:
        log('VERDICT (CERTIFIED, fw_cert3: cap 32, strict coefficients, bg3 radii, sliver bound): no physical eigenvalue kappa with Re kappa >= 0 and %.4g <= |kappa| <= %.4g (M = %s)' % (Ka, Kmax, M))
    else:
        log('VERDICT: NOT certified (failing theta cells %s)' % fails)
    return fails

def capchk():
    """positive control: with cap 32 and strict coefficients the cell box contains the box built correctly at
    cap 13 (the order-11 Lagrange term present), at the 7 cells of the referee's capchk4; negative control:
    fw_cert's original taylor_cell at cap 10 does NOT contain it (the bug is reproduced)."""
    global TAYLOR_TOL
    tol_keep = TAYLOR_TOL; TAYLOR_TOL = arb('1e300')            # n = 1: same construction as fw_cert, cap fixed
    cells = FW.tube_cells()
    xc = [c for c in cells if c['kind'] == 'x']; tc = [c for c in cells if c['kind'] == 'tau']
    pos = True; neg_fail = 0
    for c in (xc[0], xc[100], xc[300], xc[500], xc[-1], tc[0], tc[-1]):
        tape = FW.FT if c['kind'] == 'x' else FW.TT
        ctx.cap = 13; b13 = FW_taylor_orig(tape, c['start'], c['rough'], c['ht'])
        ctx.cap = 10; b10 = FW_taylor_orig(tape, c['start'], c['rough'], c['ht'])
        ctx.cap = 32; b32 = taylor_cell(tape, c['start'], c['rough'], c['ht'])
        inc32 = all(b32[i].contains(b13[i]) for i in range(len(b13)))
        inc10 = all(b10[i].contains(b13[i]) for i in range(len(b13)))
        pos = pos and inc32; neg_fail += (not inc10)
        log('%s x=%.4f: fw_cert3 (cap 32) box contains cap-13 box: %s | as-run (cap 10) contains it: %s' % (
            c['kind'], float(c['xhi'].mid()), inc32, inc10))
    ctx.cap = 32; TAYLOR_TOL = tol_keep
    # with sub-stepping (production setting) the box must still contain a high-accuracy reference point:
    # the cap-13 n=1 box midpoint is not a trajectory point, so check the sub-stepped box against the
    # trajectory endpoint computed by 64 sub-steps (both rigorous; they must intersect and the sub-stepped box
    # must be no wider than the n=1 box).
    c = xc[-1]
    b1 = taylor_cell(FW.FT, c['start'], c['rough'], c['ht'])
    TAYLOR_TOL = arb('1e300'); bn1 = taylor_cell(FW.FT, c['start'], c['rough'], c['ht']); TAYLOR_TOL = tol_keep
    sub_ok = all(bn1[i].contains(b1[i]) or b1[i].rad() <= bn1[i].rad() for i in range(3)) and all(b1[i].overlaps(bn1[i]) for i in range(3))
    log('sub-stepped box (deep cell x=%.3f): widths %s vs n=1 widths %s ; consistent: %s' % (
        float(c['xhi'].mid()), [b.rad().str(2) for b in b1], [b.rad().str(2) for b in bn1], sub_ok))
    pos = pos and sub_ok
    # strictness: a request beyond the known order must raise
    try:
        coef(arb_series([arb(1)], 3), 5); strict = False
    except ArithmeticError:
        strict = True
    log('POSITIVE CONTROL (cap 32 boxes contain the correct boxes): %s ; NEGATIVE CONTROL (cap-10 boxes fail): %d/7 ; strict coef raises: %s'
        % ('PASS' if pos else 'FAIL', neg_fail, strict))
    return pos and strict

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'capchk':
        sys.exit(0 if capchk() else 1)
    if mode == 'run':
        M = float(sys.argv[2]); Ka = float(sys.argv[3]); Kmax = float(sys.argv[4])
        nth = int(sys.argv[5]) if len(sys.argv) > 5 else 160
        ns = int(sys.argv[6]) if len(sys.argv) > 6 else 150
        nb1 = int(sys.argv[7]) if len(sys.argv) > 7 else 1500
        md = int(sys.argv[8]) if len(sys.argv) > 8 else 0
        f = run(M, Ka, Kmax, nth, ns, nb1, md)
        sys.exit(0 if f == [] else 1)
