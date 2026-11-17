"""B4 Far Wall TAIL, referee-repaired version (REFEREE-B4TAIL).  NEW file; fw_tail.py is kept as the historical run.
Changes against fw_tail.py: ctl mode with the referee's three controls (envelope-vs-cell overlap incl. power check,
Mono/cell integrand ratios, measured outer monotonicity at K = 4000 and 1e9); asserts m = 0 for A, bp, bsum; strict
theta_env assert; deep A_max from the envelope (no hard-coded constant); dead line removed; ctx.cap and sha256 of
fw_tail2.py, fw_cert3.py and lin_fw4_bg3.pkl logged; the old 'positive control' is dropped (its K-cell is
[1942.69, 2000], used here only as the source of outer WEIGHTS, a free choice).
Original description follows.
B4 Far Wall TAIL: no physical eigenvalue for Re kappa >= 0, |kappa| >= K_T (all K, uniformly).
PROOF-B4-v4-AMEND2-TAIL.md is the text.  Imports fw_cert3 (cap 32, strict coefficients, bg3, sliver) and reuses
fw_cert's integrands()/functions() unchanged, fed with new number types:
  * Sc   : theta-Laurent quantities  q = theta^m s(theta), s regular (arb_series); exact blow-up of the bulk
           coefficient functions on the E5 centre manifold (N = 1/P has a simple pole).
  * Mono : nonnegative sums  sum c xi^a y^-b  (y = K xi, xi = e^x); the deep-region bookkeeping.
Structure (u = 1/K in (0, u_T], u_T = 1/K_T):
  [O] outer cells [x_T, 0] (fw_cert3 cells down to x_T = x_c(K_T) - one cell): with weights rz = rz_T (u/u_T),
      nu = nu_T every integrand is a polynomial in u with nonnegative coefficients and every c2 term carries a
      factor u, so all outer quantities are <= their values at K_T (PROVED by inspection, sec. 2 of the text).
  [D] deep region [x_c(K), x_T]: coefficient envelopes q = xi^m qhat, |qhat| <= Qhat on theta in [0, theta_T]
      (Taylor models at theta = 0 + Lagrange on the theta-box); integrals in closed form in y (sec. 3).
  [T] centre: fw_cert.centre_theta at Ka = K_T (uniform in K >= K_T); contradiction test with u in [0, u_T].
usage: fw_tail2.py env | bulk [KT] | ctl [KT] | run [KT] [nth]
"""
import sys, os, math, time
from fractions import Fraction as Fr
ARGV = sys.argv[:]; sys.argv = ['fw_cert3.py']
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import fw_cert3 as F3
sys.argv = ARGV
import numpy as np
from flint import arb, acb, arb_series, ctx
from core import sym, sym_const
FW = F3.FW; coef = F3.coef
T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
M = 3.5
CAPS = 16          # length of the theta series in the blow-up pipeline
PT = 8             # Taylor-model degree of the envelopes

# ================================================================== Sc: theta-Laurent numbers
class Sc:
    """q = TH^m * s, TH = arb_series([theta*, 1]) (expansion point theta* = 0 or a ball), s an arb_series."""
    TH = None
    __slots__ = ('m', 's')
    def __init__(self, m, s): self.m = m; self.s = s
    @staticmethod
    def c(o):
        if isinstance(o, Sc): return o
        if isinstance(o, arb_series): return Sc(0, o)
        return Sc(0, arb_series([arb(o) if not isinstance(o, arb) else o], CAPS))
    def __add__(a, o):
        b = Sc.c(o)
        if a.m > b.m: a, b = b, a
        d = b.m - a.m
        return Sc(a.m, a.s + (Sc.TH**d)*b.s if d else a.s + b.s)
    __radd__ = __add__
    def __neg__(a): return Sc(a.m, -a.s)
    def __sub__(a, o): return a + (-Sc.c(o))
    def __rsub__(a, o): return Sc.c(o) + (-a)
    def __mul__(a, o): b = Sc.c(o); return Sc(a.m + b.m, a.s*b.s)
    __rmul__ = __mul__
    def __truediv__(a, o): b = Sc.c(o); return Sc(a.m - b.m, a.s/b.s)
    def __rtruediv__(a, o): return Sc.c(o)/a
    def __pow__(a, n):
        r = Sc(0, arb_series([arb(1)], CAPS))
        for _ in range(n): r = r*a
        return r
    def derivative(a):          # d/dtheta (TH^m s) = TH^{m-1} (m s + TH s')
        return Sc(a.m - 1, a.m*a.s + Sc.TH*a.s.derivative()) if a.m else Sc(-1, Sc.TH*a.s.derivative())

_GT_TAPE = FW.GT
class _GTSc:
    """fw_cert.GT evaluated on Sc inputs (constants as Sc)."""
    def series_eval(self, ins, cap):
        t = _GT_TAPE; vals = [None]*len(t.nodes)
        for i, (op, a, b) in enumerate(t.nodes):
            if op == 'var': vals[i] = Sc.c(ins[a])
            elif op == 'const': vals[i] = Sc.c(sym_const(a))
            elif op == 'add': vals[i] = vals[a] + vals[b]
            elif op == 'mul': vals[i] = vals[a]*vals[b]
            elif op == 'div': vals[i] = vals[a]/vals[b]
        return [vals[j] for j in t.outs]

def deep_functions(mode, Theta=None):
    """fw_cert.functions on the centre manifold in blown-up form.  mode 'thin': expansion at theta = 0;
    'box': coefficients enclose the Taylor coefficients at every theta* in Theta (Cauchy tails)."""
    BG = FW.BG; cm = BG['cm']; nu = arb(BG['nu_c'])
    epsY = arb(cm['eps'][0]) + arb(cm['eps'][1]); epsX = arb(cm['epsX'][0]) + arb(cm['epsX'][1])
    Yc = [[FW.A_(a) for a in cm['Y'][i]] for i in range(3)]; Xc = [FW.A_(a) for a in cm['X']]
    assert Yc[0][0].is_zero() and (Yc[0][1] - 1).is_zero()
    Pt = Yc[0][1:]                                           # P = theta * Ptilde exactly (P_0 = 0)
    if mode == 'thin':
        Sc.TH = arb_series([arb(0), arb(1)], CAPS)
        jet = lambda c, eps, rr: arb_series([(c[k] if k < len(c) else arb(0)) + sym(eps/rr**k) for k in range(CAPS)], CAPS)
        rr = nu
    else:
        Sc.TH = arb_series([Theta, arb(1)], CAPS)
        rr = nu - Theta.abs_upper()
        def jet(c, eps, rr):
            out = []
            for k in range(CAPS):
                s = arb(0)
                for n_ in range(k, len(c)):
                    s += c[n_]*arb(math.comb(n_, k))*Theta**(n_ - k)
                out.append(s + sym(eps/rr**k))
            return arb_series(out, CAPS)
    Ptil = jet(Pt, epsY/nu, rr)                              # |R_P/theta| <= eps/nu on |theta| <= nu
    W = jet(Yc[1], epsY, rr); U = jet(Yc[2], epsY, rr); chi = jet(Xc, epsX, rr)
    N = Sc(-1, 1/Ptil); om = Sc(2, W*Ptil*Ptil); V = Sc(1, U*Ptil)
    Xd = Sc(-1, 1 + Sc.TH*chi.derivative())                  # dx/dtheta, x = log theta + chi + xc
    Y = [N, om, V]
    Fx = [y.derivative()/Xd for y in Y]
    keepS, keepGT = FW.S, FW.GT
    FW.S = lambda c_: Sc.c(c_); FW.GT = _GTSc()
    try:
        f = FW.functions(Y, Fx, Xd)
    finally:
        FW.S, FW.GT = keepS, keepGT
    f['phi'] = f['ellphi']/f['ell']; f['bpb'] = f['bp'] + f['b11']/f['ell']; f['bm'] = f['b11']/f['ell']
    f['tau'] = f['b01']*f['invlp']; f['aeAalm'] = f['eAalm']/f['ell']; f['bsum'] = f['bpb']
    return f, chi

def envelopes(theta_T):
    """q = xi^m qhat on theta in [0, theta_T], qhat = (p e^{-chi})^m s: ball enclosure of qhat (Taylor model)."""
    Theta = theta_T/2 + sym(theta_T/2)
    ft, chit = deep_functions('thin'); fb, chib = deep_functions('box', Theta)
    p = FW.PS
    env = {}
    for k in ft:
        m = ft[k].m; assert fb[k].m == m
        qt = (p**m)*(-m*chit).exp()*ft[k].s if m else ft[k].s
        qb = (p**m)*(-m*chib).exp()*fb[k].s if m else fb[k].s
        val = sum((coef(qt, j)*Theta**j for j in range(PT)), arb(0)) + coef(qb, PT)*Theta**PT
        env[k] = (m, val)
    return env

# ================================================================== Mono: sum c xi^a y^-b (c >= 0)
class Mono:
    __slots__ = ('t',)
    def __init__(self, t=None): self.t = dict(t or {})
    @staticmethod
    def c(o):
        if isinstance(o, Mono): return o
        a = arb(o) if not isinstance(o, arb) else o
        return Mono({(Fr(0), Fr(0)): arb(a.abs_upper())}) if not a.is_zero() else Mono()
    def __add__(a, o):
        b = Mono.c(o); t = dict(a.t)
        for k, v in b.t.items(): t[k] = t[k] + v if k in t else v
        return Mono(t)
    __radd__ = __add__
    def __mul__(a, o):
        b = Mono.c(o); t = {}
        for (a1, b1), v1 in a.t.items():
            for (a2, b2), v2 in b.t.items():
                k = (a1 + a2, b1 + b2); t[k] = t[k] + v1*v2 if k in t else v1*v2
        return Mono(t)
    __rmul__ = __mul__
    def __pow__(a, n):
        r = Mono({(Fr(0), Fr(0)): arb(1)})
        for _ in range(n): r = r*a
        return r
    def __repr__(self): return ' + '.join('%s xi^%s y^-%s' % (v.str(3), a, b) for (a, b), v in sorted(self.t.items()))

def mono(m, val):
    return Mono({(Fr(m), Fr(0)): arb(val.abs_upper())})

def integ(I, e, xiT):
    """bound of int_t^{x_T} (y_s/y_t)^e I(s) ds for every t in the deep region, as Mono in (xi_t, y_t).
    xi_s = xi_t y_s/y_t exactly (same K).  term c xi^a y^-b:  cc = a - b + e;
      cc < 0 : c xi_t^a y_t^-b / (-cc)          (integral to infinity)
      cc > 0 : c xi_T^cc xi_t^(b - e) y_t^-b / cc  (integral to y_T = y_t xi_T/xi_t)."""
    out = {}
    for (a, b), c in I.t.items():
        ee = e
        cc = a - b + ee
        if abs(cc) < CCMIN:             # raise the exponent: (y_s/y_t)^e <= (y_s/y_t)^(e + d), y_s >= y_t
            ee = e + (CCMIN - cc); cc = a - b + ee
        e_ = ee
        if cc < 0:
            k = (a, b); v = c*arb(cc.denominator)/arb(-cc.numerator)
        else:
            k = (b - e_, b); v = c*_pw(xiT, cc)*arb(cc.denominator)/arb(cc.numerator)
        out[k] = out[k] + v if k in out else v
    return Mono(out)

def _pw(x, q):
    """upper bound of x^q, x > 0 a ball (its upper end is used when q >= 0, lower end when q < 0), q a Fraction."""
    q = Fr(q)
    if q == 0: return arb(1)
    xb = arb(x.upper()) if q > 0 else arb(x.lower())
    assert xb > 0
    return arb((xb.log()*arb(q.numerator)/arb(q.denominator)).exp().upper())

def ratio_sup(C, W, xiT, Ym):
    """sup over xi in (0, xiT], y >= Ym of C/W (C, W Monos, W > 0):  sum over terms of C of
    min over terms of W of (c/w) xi_T^(a-a') Ym^(b'-b), finite only if a >= a', b >= b'."""
    tot = arb(0)
    for (a, b), c in C.t.items():
        if c.is_zero(): continue
        best = None
        for (aw, bw), w in W.t.items():
            if a >= aw and b >= bw:
                val = c/w*_pw(xiT, a - aw)*_pw(1/Ym, b - bw)
                best = val if best is None or val < best else best
        if best is None:
            raise ArithmeticError('term %s xi^%s y^-%s not dominated by the weight %s' % (c.str(3), a, b, W))
        tot += best
    return arb(tot.upper())

def eval_at(C, xiT, Ym):
    """sup over xi in (0, xiT] of C at y = Ym (needs a >= 0)."""
    tot = arb(0)
    for (a, b), c in C.t.items():
        if c.is_zero(): continue
        assert a >= 0, 'output term with negative xi power'
        tot += c*_pw(xiT, a)*_pw(1/Ym, b)
    return arb(tot.upper())

def rat(x, digits=6, upward=True):
    """a Fraction upper (or lower) bound of the arb x."""
    import decimal
    decimal.getcontext().prec = 50
    v = x.upper() if upward else x.lower()
    d = decimal.Decimal(v.mid().str(30, radius=False).replace('e', 'E'))
    q = decimal.Decimal(1).scaleb(d.adjusted() - digits + 1)
    d2 = d.quantize(q, rounding=decimal.ROUND_CEILING if upward else decimal.ROUND_FLOOR) + (q if upward else -q)
    return Fr(str(d2))

# ================================================================== [O] outer cells at K_T with u-monotone weights
UKEYS = dict(invlp='invlp', lp='lp', sig='sig', dsig='dsig', chi='chi', dchi='dchi', phi='phi', bpb='bpb', bm='bm',
             tau='tau', eAp='eAp', eNp='eNp', eAm='eAm', eNm='eNm', deAm='deAm', deNm='deNm', alp='al_p', alm='al_m',
             cN='cN', gap='ga_p', gam='ga_m', A='A', dA='dA', H3='H3', H3lp='H3lp', dH3lp='dH3lp', bpp='b00', bpm='b01',
             dAalp='dAalp', dAalm='dAalm', Aalmlp='Aalmlp', aqd='q', dqd='dq', eAcN='eAcN', deAcN='deAcN', aeAalm='aeAalm')

def build_cells(KT):
    cells = FW.eval_cells(FW.tube_cells()); sl = F3.sliver(cells)
    th_min = float((FW.PS*arb(M)/arb(KT)).lower())*0.97
    cells = cells + FW.centre_bulk_cells(th_min)
    son = FW.sonic_piece(); cells, son = FW.prep_bulk(cells, son)
    xcT = (arb(M)/arb(KT)).log()
    outer = [c for c in cells if c['xlo'].lower() >= xcT.upper()]         # cells entirely above x_c(K_T)
    assert outer == cells[:len(outer)]
    return outer, son, xcT, cells

def outer_part(outer, son, KT, weights):
    """fw_cert.bulk_ia's recursion at K = K_T on the outer cells with FIXED weights (rz_T, nu_T); returns the
    outer rows, the junction values at x_T and the sonic constant.  All outer quantities are nondecreasing in
    u = 1/K when rz = rz_T u/u_T and nu = nu_T (every operation is + or x of nonnegative numbers; iK, iK2, iKA,
    rz increase with u; Ku rz = K_T rz_T is constant), and every c2/I2 term carries a factor u."""
    Ka = arb(KT)
    kp = dict(iK=1/Ka, iK2=1/Ka**2, iKA=1/(Ka - son['Amax']), Ku=Ka)
    rz, nu, rzS, nuS = weights
    s0, g, eSP = son['s0'], son['g'], son['eSP']
    J = FW.integrands(son['u'], kp, rzS, nuS)
    G2 = FW.vsc(eSP, FW.vadd(FW.vsc(1/g, J['S2']), FW.vsc(s0, J['B2'])))
    G3 = FW.vsc(eSP*s0, J['I3']); Wb = eSP; Fw0lp0 = J['Fwlp']
    R1 = FW.vsc(s0, J['r1'])
    c2 = [FW.vadd(G2, J['p2'])]; c3 = [FW.vadd(G3, J['p3'], FW.vsc(Fw0lp0*eSP, FW.E1))]; wts = [(rzS, nuS)]
    for jj, c in enumerate(outer):
        u = c['u']; h = c['hU']
        Ji = FW.integrands(u, kp, rz[jj], nu[jj])
        EG = (c['h']*c['bG']).upper().exp(); OG = max(arb(0), (c['h']*c['bG']).upper()).exp()
        EP = (c['h']*c['bP']).upper().exp(); OP = max(arb(0), (c['h']*c['bP']).upper()).exp()
        c2.append(FW.vadd(FW.vsc(OG, FW.vadd(FW.vsc(h, Ji['I2']), G2)), Ji['p2']))
        c3.append(FW.vadd(FW.vsc(OP, FW.vadd(FW.vsc(h, Ji['I3']), G3)), Ji['p3'], FW.vsc(Fw0lp0*OP*Wb, FW.E1)))
        G2 = FW.vadd(FW.vsc(OG*h, Ji['I2']), FW.vsc(EG, G2))
        G3 = FW.vadd(FW.vsc(OP*h, Ji['I3']), FW.vsc(EP, G3))
        Wb = EP*Wb
        R1 = FW.vadd(R1, FW.vsc(h, Ji['r1']))
        wts.append((rz[jj], nu[jj]))
    A2 = [max(FW.up(c2[j][k]/wts[j][0]) for j in range(len(c2))) for k in range(3)]
    A3 = [max(FW.up(c3[j][k]/wts[j][1]) for j in range(len(c3))) for k in range(3)]
    return dict(R1=[FW.up(t) for t in R1], A2=A2, A3=A3, G2=[FW.up(t) for t in G2], G3=[FW.up(t) for t in G3],
                Wb=FW.up(Wb), Fw0lp0=FW.up(Fw0lp0))

# ================================================================== [D] deep region
def theta_env_bound(xT):
    """theta_env >= theta(x) for every x <= x_T on the centre manifold (x = log theta + chi + xc)."""
    cm = FW.BG['cm']; Xc = [FW.A_(a) for a in cm['X']]; epsX = arb(cm['epsX'][0]) + arb(cm['epsX'][1])
    th = arb((xT - FW.XC).exp().upper())*arb('1.1')
    for _ in range(3):
        chib = sum((arb(Xc[n_].abs_upper())*th**n_ for n_ in range(len(Xc))), arb(0)) + epsX
        th = arb(((xT - FW.XC).exp()*chib.exp()).upper())*arb('1.01')
    chib = sum((arb(Xc[n_].abs_upper())*th**n_ for n_ in range(len(Xc))), arb(0)) + epsX
    F = arb(((xT - FW.XC).exp()*chib.exp()).upper())
    assert th > F, 'theta_env not strictly self-consistent'          # strict: x(theta_env) > x_T (AMEND2 sec. 1.2)
    assert th < arb(FW.BG['nu_c']), 'theta_env outside the centre-manifold disc'
    return th

CCMIN = Fr(1, 20)

def deep_part(env, KT, xiT, junc, rho1, rho2, nu1, alpha=None, beta=None):
    uT = 1/arb(KT)
    ud = {k: mono(*env[v]) for k, v in UKEYS.items()}
    assert env['A'][0] == 0 and env['bp'][0] == 0 and env['bsum'][0] == 0, 'A, bp, bsum must have m = 0'
    Amax = arb(env['A'][1].abs_upper())                                   # deep region: local sup of A
    aT = 1/(1 - Amax*uT)
    BG_ = rat(arb(max(arb(0), env['bsum'][1].upper()).upper())) + Fr(1, 997)    # Gamma(s,t) <= (y_s/y_t)^BG_
    eP = rat(arb(env['bp'][1].upper())) + Fr(1, 991)                             # P_+(s,t) <= (y_s/y_t)^eP (eP < 0)
    assert eP < 0
    an = -eP if alpha is None else min(Fr(alpha), -eP)
    br = (1 - BG_) if beta is None else min(Fr(beta), 1 - BG_)
    kp =dict(iK=Mono({(Fr(1), Fr(1)): arb(1)}), iK2=Mono({(Fr(2), Fr(2)): arb(1)}), iKA=Mono({(Fr(1), Fr(1)): arb(aT.upper())}),
              Ku=Mono({(Fr(-1), Fr(-1)): arb(1)}))
    rz = Mono({(br, Fr(1)): rho1, (Fr(0), Fr(2)): rho2})
    nu = Mono({(an, Fr(0)): nu1})
    J = FW.integrands(ud, kp, rz, nu)
    c2 = []; c3 = []; r1 = []
    for k in range(3):
        j2 = Mono({(1 - BG_, Fr(1)): junc['G2'][k]*arb(KT)*_pw(xiT, BG_)})
        c2.append(integ(J['I2'][k], BG_, xiT) + J['p2'][k] + j2)
        j3v = junc['G3'][k] + (junc['Fw0lp0']*junc['Wb'] if k == 0 else arb(0))
        j3 = Mono({(an, Fr(0)): j3v*_pw(xiT, eP)})
        c3.append(integ(J['I3'][k], eP, xiT) + J['p3'][k] + j3)
        r1.append(eval_at(integ(J['r1'][k], Fr(1, 1000), xiT), xiT, arb(M)))
    return dict(c2=c2, c3=c3, r1=r1, rz=rz, nu=nu, BG=BG_, eP=eP, aT=aT)

def assemble(outer_res, deep, xiT):
    A = [[outer_res['R1'][k] + deep['r1'][k] for k in range(3)],
         [max(outer_res['A2'][k], ratio_sup(deep['c2'][k], deep['rz'], xiT, arb(M))) for k in range(3)],
         [max(outer_res['A3'][k], ratio_sup(deep['c3'][k], deep['nu'], xiT, arb(M))) for k in range(3)]]
    Af = np.array([[float(a.mid()) for a in row] for row in A])
    w_, V_ = np.linalg.eig(Af); pv = np.abs(V_[:, np.argmax(abs(w_))].real) + 1e-6
    y = [arb(float(t)) for t in pv]
    coll = max(FW.up(sum((A[i][j]*y[j] for j in range(3)), arb(0))/y[i]) for i in range(3))
    res = dict(A=A, collatz=coll, ok=bool(coll < 1))
    if not res['ok']: return res
    I_A = [[(1 if i == j else 0) - A[i][j] for j in range(3)] for i in range(3)]
    det = (I_A[0][0]*(I_A[1][1]*I_A[2][2] - I_A[1][2]*I_A[2][1]) - I_A[0][1]*(I_A[1][0]*I_A[2][2] - I_A[1][2]*I_A[2][0])
           + I_A[0][2]*(I_A[1][0]*I_A[2][1] - I_A[1][1]*I_A[2][0]))
    v = (FW.up((I_A[1][1]*I_A[2][2] - I_A[1][2]*I_A[2][1])/det), FW.up(-(I_A[1][0]*I_A[2][2] - I_A[1][2]*I_A[2][0])/det),
         FW.up((I_A[1][0]*I_A[2][1] - I_A[1][1]*I_A[2][0])/det))
    eb = eval_at(sum((deep['c2'][k]*v[k] for k in range(3)), Mono()), xiT, arb(M))
    nbar = eval_at(sum((deep['c3'][k]*v[k] for k in range(3)), Mono()), xiT, arb(M))
    res.update(v=v, eb=eb, nbar=nbar)
    return res

# ================================================================== [T] contradiction test, u = 1/K in [0, u_T]
def tail_test(bk, env, cen, KT, xiT, th_lo, th_hi):
    thb = (th_lo + th_hi)/2 + sym((th_hi - th_lo)/2)
    eps = acb(0, thb).exp()
    uT = 1/arb(KT); u = uT/2 + sym(uT/2)
    Mb = arb(M)
    mq = env['q'][0]; assert mq + 2 >= 0
    Xq = _pw(xiT, Fr(mq + 2))/2 + sym(_pw(xiT, Fr(mq + 2))/2)
    cb = acb(env['sig'][1])/(Mb*eps) + acb(env['chi'][1])/(Mb**2*eps**2) - acb(env['q'][1]*Xq)/(Mb**2*eps*(eps - acb(env['A'][1]*u)))
    assert env['sig'][0] == -1 and env['chi'][0] == -2
    def Yball(m):
        Yf, B1, E2 = cen[m]
        return [Yf[a] + acb(u)*B1[a] + acb(sym(E2[a]*uT**2), sym(E2[a]*uT**2)) for a in range(4)]
    Yw, YN = Yball('w'), Yball('N')
    L = [env[k][1] for k in ('L00', 'L01', 'L10', 'L11')]
    assert all(env[k][0] == 0 for k in ('L00', 'L01', 'L10', 'L11'))
    Aw = [acb(L[0])*Yw[2] + acb(L[1])*Yw[3], acb(L[2])*Yw[2] + acb(L[3])*Yw[3]]
    Bw = [acb(L[0])*YN[2] + acb(L[1])*YN[3], acb(L[2])*YN[2] + acb(L[3])*YN[3]]
    a, bN = Yw[1], YN[1]
    nK = bk['nbar']*uT
    den = arb(bN.abs_lower()) - nK*arb(Bw[0].abs_upper())
    if not den > 0: return False, None
    mu0 = -a/bN
    r = nK*(arb(Aw[0].abs_upper()) + arb(mu0.abs_upper())*arb(Bw[0].abs_upper()))/den
    g0 = Aw[1] - cb*Aw[0]; g1 = Bw[1] - cb*Bw[0]
    lower = arb((g0 + mu0*g1).abs_lower()) - r*arb(g1.abs_upper())
    wp = arb((Aw[0] + mu0*Bw[0]).abs_upper()) + r*arb(Bw[0].abs_upper())
    rhs = bk['eb']*wp
    return bool(lower > rhs), (lower, rhs)

def bulk_tail(KT, verbose=True):
    outer, son, xcT, allc = build_cells(KT)
    xT = outer[-1]['xlo']; xiT = arb(xT.exp().upper())
    thT = theta_env_bound(arb(xT.upper()))
    if verbose: log('[tail] K_T=%g: %d outer cells, x_T=%s (x_c(K_T)=%s), xi_T<=%s, theta_env=%s' % (
        KT, len(outer), xT.str(6), xcT.str(6), xiT.str(4), thT.str(4)))
    # weights at K_T from the fw_cert3 bulk (bulk_best) on exactly these cells
    rT = F3.bulk_best(allc, son, M, arb(KT)/arb('1.0295'), arb(KT))     # K in [1942.69, 2000]: source of WEIGHTS only
    assert rT['ok']
    rz, nu, rzS, nuS = rT['weights']
    out = outer_part(outer, son, KT, (rz[:len(outer)], nu[:len(outer)], rzS, nuS))
    env = envelopes(thT)
    if verbose: log('[env] xi sig -> %s, xi^2 chi -> %s, beta_+ in %s, beta_+ + beta_- in %s, A in %s' % (
        env['sig'][1].str(4), env['chi'][1].str(4), env['bp'][1].str(4), env['bsum'][1].str(3), env['A'][1].str(4)))
    # deep weights (any positive weights are admissible): deterministic coordinate descent on the rigorous
    # Collatz bound over (rho1, rho2, nu1) on a log grid and the exponents (alpha_nu, beta_rz)
    def trial(w):
        try:
            dp = deep_part(env, KT, xiT, out, arb(w[0]), arb(w[1]), arb(w[2]), w[3], w[4])
            res = assemble(out, dp, xiT)
        except ArithmeticError as e:
            return None
        res['w'] = w; res['deep'] = dp
        return res
    w = [100.0, 0.3, 1000.0, Fr(9, 10), Fr(9, 10)]
    best = trial(w)
    steps = [('m', 0, 2.0), ('m', 1, 2.0), ('m', 2, 2.0), ('e', 3, Fr(1, 20)), ('e', 4, Fr(1, 20))]
    for sweep in range(8):
        improved = False
        for kind, i, s in steps:
            for sgn in (1, -1):
                w2 = list(w)
                if kind == 'm': w2[i] = w[i]*(s if sgn > 0 else 1/s)
                else: w2[i] = min(Fr(99, 100), max(Fr(1, 2), w[i] + sgn*s))
                r = trial(w2)
                if r is not None and (best is None or r['collatz'] < best['collatz']):
                    best, w, improved = r, w2, True
        if verbose: log('   sweep %d: collatz %s at weights %s' % (sweep, best['collatz'].str(5) if best else None, w))
        if not improved: break
    best.update(outer=out, env=env, xiT=xiT, thT=thT, rT=rT, son=son)
    return best

# ================================================================== main
def provenance():
    import hashlib
    def h(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
    d = '<LAB>/engine/cert/'
    log('[provenance] ctx.cap=%d prec=%d sha256 fw_tail2.py=%s fw_cert3.py=%s lin_fw4_bg3.pkl=%s' % (
        ctx.cap, ctx.prec, h(d + 'fw_tail2.py'), h(d + 'fw_cert3.py'), h(d + 'lin_fw4_bg3.pkl')))

def controls(KT):
    """REFEREE-B4TAIL controls (reft_chk1/2 logic), each with a PASS criterion:
    C1 envelope-vs-cell overlap: on every fw_cert3 cell below x_T, xi^m * envelope must overlap the cell enclosure
       for every key (0 non-overlaps); power check: on the 3 cells just above x_T at least one key must NOT overlap.
    C2 Mono/cell integrand ratios: Mono integrands at the sup over each deep cell divided by the cell integrands at
       K_T with the smallest weight on the cell, for r1, I2, p2, I3, p3: must lie in [0.9, 1.25] (a factor-K or
       exponent error gives ratios off by >= y or xi^{-1}).
    C3 measured outer monotonicity: outer_part at K = 4000 and 1e9 (rz = rz_T K_T/K, nu = nu_T) must give
       R1, A2, A3, G3, Wb, Fw0lp0 <= their K_T values and G2 <= (u/u_T) G2(K_T)."""
    outer, son, xcT, allc = build_cells(KT)
    n = len(outer)
    xT = outer[-1]['xlo']; thT = theta_env_bound(arb(xT.upper())); env = envelopes(thT)
    def extra(f):
        g = dict(f)
        g['phi'] = f['ellphi']/f['ell']; g['bpb'] = f['bp'] + f['b11']/f['ell']; g['bm'] = f['b11']/f['ell']
        g['tau'] = f['b01']*f['invlp']; g['aeAalm'] = f['eAalm']/f['ell']; g['bsum'] = g['bpb']
        return g
    deep = allc[n:]
    bad = 0; tested = 0; power = 0
    for c in deep + outer[-3:]:
        xi = c['xlo'].union(c['xhi']).exp()
        f = extra(c['f'])
        inside = bool(c['xhi'].upper() <= xT.upper() + arb('1e-12'))
        for k in sorted(env):
            if k not in f: continue
            m, val = env[k]
            pred = xi**m*val if m >= 0 else val/xi**(-m)
            ov = pred.overlaps(f[k])
            if inside:
                tested += 1; bad += (not ov)
            else:
                power += (not ov)
    c1 = (bad == 0 and tested > 0 and power > 0)
    log('[ctl C1] deep cells %d: %d key comparisons, non-overlaps %d; power check above x_T: %d non-overlaps -> %s' % (
        len(deep), tested, bad, power, 'PASS' if c1 else 'FAIL'))
    # C2
    rT = F3.bulk_best(allc, son, M, arb(KT)/arb('1.0295'), arb(KT))
    rz, nu, rzS, nuS = rT['weights']
    out = outer_part(outer, son, KT, (rz[:n], nu[:n], rzS, nuS))
    b = bulk_tail_from(outer, son, xcT, allc, KT, rT, out, verbose=False)
    w = b['w']; dp = b['deep']
    ud = {k: mono(*env[v]) for k, v in UKEYS.items()}
    uT = 1/arb(KT); aT = 1/(1 - arb(env['A'][1].abs_upper())*uT)
    an = min(Fr(w[3]), -dp['eP']); br = min(Fr(w[4]), 1 - dp['BG'])
    kpM = dict(iK=Mono({(Fr(1), Fr(1)): arb(1)}), iK2=Mono({(Fr(2), Fr(2)): arb(1)}), iKA=Mono({(Fr(1), Fr(1)): arb(aT.upper())}),
               Ku=Mono({(Fr(-1), Fr(-1)): arb(1)}))
    rzM = Mono({(br, Fr(1)): arb(w[0]), (Fr(0), Fr(2)): arb(w[1])}); nuM = Mono({(an, Fr(0)): arb(w[2])})
    JM = FW.integrands(ud, kpM, rzM, nuM)
    Ka = arb(KT)
    def pw(x, q): return (x.log()*arb(q.numerator)/arb(q.denominator)).exp()
    def msup(M_, xlo, xhi):
        tot = arb(0)
        for (a, bb), v in M_.t.items():
            xi = arb(xhi.upper()).exp() if a - bb >= 0 else arb(xlo.lower()).exp()
            tot += v*pw(xi, a)*pw(Ka*xi, -bb)
        return tot
    def minf(M_, xlo, xhi):
        tot = arb(0)
        for (a, bb), v in M_.t.items():
            xi = arb(xlo.lower()).exp() if a - bb >= 0 else arb(xhi.upper()).exp()
            tot += v*pw(xi, a)*pw(Ka*xi, -bb)
        return tot
    kpC = dict(iK=1/Ka, iK2=1/Ka**2, iKA=1/(Ka - son['Amax']), Ku=Ka)
    lo, hi = 9.0, 0.0
    for c in deep:
        JC = FW.integrands(c['u'], kpC, minf(rzM, c['xlo'], c['xhi']), minf(nuM, c['xlo'], c['xhi']))
        for key in ('r1', 'I2', 'p2', 'I3', 'p3'):
            for i in range(3):
                if JC[key][i] > 0:
                    r = float((msup(JM[key][i], c['xlo'], c['xhi'])/JC[key][i]).mid())
                    lo, hi = min(lo, r), max(hi, r)
    c2 = (0.9 <= lo and hi <= 1.25)
    log('[ctl C2] Mono(sup)/cell integrand ratios on %d deep cells: min %.4f max %.4f (band [0.9, 1.25]) -> %s' % (
        len(deep), lo, hi, 'PASS' if c2 else 'FAIL'))
    # C3
    c3 = True
    for K in (4000.0, 1e9):
        s = arb(KT)/arb(K)
        o = outer_part(outer, son, K, ([r*s for r in rz[:n]], nu[:n], rzS*s, nuS))
        okk = all(bool(o[key][i] <= out[key][i]) for key in ('R1', 'A2', 'A3', 'G3') for i in range(3))
        okk = okk and all(bool(o['G2'][i] <= out['G2'][i]*s) for i in range(3))
        okk = okk and bool(o['Wb'] <= out['Wb']) and bool(o['Fw0lp0'] <= out['Fw0lp0'])
        g2r = max(float((o['G2'][i]/(out['G2'][i]*s)).mid()) for i in range(3))
        log('[ctl C3] K=%g: outer quantities <= K_T values: %s (max G2/((u/u_T) G2_T) = %.4f, A2/A2_T = %s)' % (
            K, okk, g2r, [round(float((o['A2'][i]/out['A2'][i]).mid()), 4) for i in range(3)]))
        c3 = c3 and okk
    log('[ctl C3] measured outer monotonicity -> %s' % ('PASS' if c3 else 'FAIL'))
    allp = c1 and c2 and c3
    log('CONTROLS: C1 %s, C2 %s, C3 %s -> %s' % (c1, c2, c3, 'ALL PASS' if allp else 'FAIL'))
    return allp

def bulk_tail_from(outer, son, xcT, allc, KT, rT, out, verbose=True):
    xT = outer[-1]['xlo']; xiT = arb(xT.exp().upper())
    thT = theta_env_bound(arb(xT.upper()))
    env = envelopes(thT)
    def trial(w):
        try:
            dp = deep_part(env, KT, xiT, out, arb(w[0]), arb(w[1]), arb(w[2]), w[3], w[4])
            res = assemble(out, dp, xiT)
        except ArithmeticError:
            return None
        res['w'] = w; res['deep'] = dp
        return res
    w = [100.0, 0.3, 1000.0, Fr(9, 10), Fr(9, 10)]
    best = trial(w)
    steps = [('m', 0, 2.0), ('m', 1, 2.0), ('m', 2, 2.0), ('e', 3, Fr(1, 20)), ('e', 4, Fr(1, 20))]
    for sweep in range(8):
        improved = False
        for kind, i, s in steps:
            for sgn in (1, -1):
                w2 = list(w)
                if kind == 'm': w2[i] = w[i]*(s if sgn > 0 else 1/s)
                else: w2[i] = min(Fr(99, 100), max(Fr(1, 2), w[i] + sgn*s))
                r = trial(w2)
                if r is not None and (best is None or r['collatz'] < best['collatz']):
                    best, w, improved = r, w2, True
        if not improved: break
    best.update(outer=out, env=env, xiT=xiT, thT=thT, rT=rT, son=son)
    return best

if __name__ == '__main__':
    mode = sys.argv[1]
    KT = float(sys.argv[2]) if len(sys.argv) > 2 else 2000.0
    provenance()
    if mode == 'ctl':
        sys.exit(0 if controls(KT) else 1)
    if mode == 'env':
        thT = arb('0.0014')
        env = envelopes(thT)
        for k in sorted(env):
            log('  %-8s m=%2d  qhat in %s' % (k, env[k][0], env[k][1].str(5)))
    elif mode in ('bulk', 'run'):
        nth = int(sys.argv[3]) if len(sys.argv) > 3 else 160
        b = bulk_tail(KT)
        A = b['A']
        log('[bulk tail] weights %s, A = %s' % (b.get('w'), [[a.str(4) for a in row] for row in A]))
        log('[bulk tail] Collatz %s ok=%s' % (b['collatz'].str(5), b['ok']))
        if not b['ok']:
            log('VERDICT: tail bulk NOT certified'); sys.exit(1)
        log('[bulk tail] UNIFORM for all K >= %g: v=%s e_b <= %s, nbar <= %s' % (KT, [t.str(4) for t in b['v']], b['eb'].str(5), b['nbar'].str(5)))
        if mode == 'bulk': sys.exit(0)
        cdat_cm = FW.centre_manifold_data(); FW.cdat_C01, FW.cdat_C12 = cdat_cm['C01'], cdat_cm['C12']
        cdat = cdat_cm['data'] + [FW.bt_centre_data(cdat_cm['C01'], cdat_cm['C12'])]
        pi2 = arb.pi()/2; ths = [pi2*k/nth for k in range(nth + 1)]
        sub = [int(t) for t in os.environ.get('FWT_SUB', '').split(',') if t]
        fails = []; worst = None
        for k in (sub if sub else range(nth)):
            cen = FW.centre_theta(ths[k], ths[k + 1], KT, M, cdat, ns=150)
            if cen['w'] is None or cen['N'] is None:
                fails.append(k); log('  theta cell %d: centre supersolution FAIL' % k); continue
            ok, info = tail_test(b, b['env'], cen, KT, b['xiT'], ths[k], ths[k + 1])
            if info is not None:
                rt = float((info[0]/info[1]).lower()) if info[1] > 0 else float('inf')
                worst = rt if worst is None else min(worst, rt)
            if not ok:
                fails.append(k); log('  theta cell %d FAIL %s' % (k, info))
        log('theta cells %s, failures %s, worst lower/rhs %s' % (len(sub) if sub else nth, fails, worst))
        if sub:
            log('PROBE ONLY')
        elif not fails:
            log('VERDICT (CERTIFIED, fw_tail2): no physical eigenvalue kappa with Re kappa >= 0 and |kappa| >= %g (M = %s)' % (KT, M))
        else:
            log('VERDICT: tail NOT certified (failing theta cells %s)' % fails)
