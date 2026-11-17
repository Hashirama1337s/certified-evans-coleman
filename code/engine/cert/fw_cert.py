"""B4 Far Wall CERTIFICATE (proof seat v4, new file).  PROOF-B4-v4.md is the text; this script evaluates
every constant of the v4 chain in ball arithmetic (python-flint arb/acb) on the BAR-LB-E5 background cover
(lin_fw4_bg.pkl, produced by lin_fw4_bg.py with the sealed E5 instrument) and prints the rigorous R0.

Structure (labels as in PROOF-B4-v4 sec. 6):
  [B] bulk  : cells = E5 tau-leg + x-leg Lohner records (+ PE2 margin), sonic piece theta in [th1, 0] by the
              certified sonic-manifold series with Cauchy tails.  Interval jets (Picard on the x-field with
              arb_series) give every coefficient function and its x-derivatives as enclosures over each cell.
              Bootstrap matrix, weights, Collatz bound rho(A) < 1, (I - A)^{-1} e1, outputs e_b, nbar, c_b,
              L(x_c) -- all upper bounds, per K-cell [Ka, Kb].
  [C] centre: flat closed-form kernel on theta-cells (acb), Delta majorant uniform in K >= Ka on the real
              v-segment (tube cells + centre-manifold cells), Volterra supersolution, K-uniform signed term
              (1/K) B1^(1)(theta); valid for all K >= Ka at once.
  [T] test  : closed-form mu-disc contradiction in ball arithmetic per (K-cell, theta-cell).
usage: fw_cert.py M Klo Khi [dK] [ntheta]
"""
import sys, os, math, time, pickle
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import numpy as np
import sympy as sp
from flint import arb, acb, arb_series, ctx
from core import Tape, sym, sym_const
PREC = 128
ctx.prec = PREC
T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)

# ------------------------------------------------------------------ symbolic model (sealed derivation)
import lin_float as LF                       # HKA pencil (checked in lin_sym.py) + x-field f_bg (system.py)
N_, om_, V_ = LF.N, LF.om, LF.V
cs_ = sp.Symbol('cs')                        # cs = sqrt(1/3), passed as a ball (Tape needs rational constants)
lam_p_s = 1 + N_*(V_ + cs_)/(1 + cs_*V_)
lam_m_s = 1 + N_*(V_ - cs_)/(1 - cs_*V_)
FT = Tape(LF.f_bg, [N_, om_, V_])
G_OUT = list(LF.Ms_s) + list(LF.Mx_s) + list(LF.P0_s) + list(LF.I0_s) + [LF.A_s, lam_p_s, lam_m_s]
GT = Tape(G_OUT, [N_, om_, V_, LF.omx_s, LF.Vx_s, cs_])
CS = (arb(1)/3).sqrt()
CAP = 12
PTM = 5

def S(c):                                    # constant series
    return arb_series([c], CAP)

def ser_inv2(m):
    a, b, c, d = m
    det = a*d - b*c
    return [d/det, -b/det, -c/det, a/det]
def mm(A, B):                                 # 2x2 (row-major lists) product
    return [A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3], A[2]*B[0] + A[3]*B[2], A[2]*B[1] + A[3]*B[3]]

import system as S0
_Ms, _bs = S0.Msym, S0.bsym
_det = _Ms[0, 0]*_Ms[1, 1] - _Ms[0, 1]*_Ms[1, 0]
TT = Tape([_det*S0.gN, _Ms[1, 1]*_bs[0] - _Ms[0, 1]*_bs[1], _Ms[0, 0]*_bs[1] - _Ms[1, 0]*_bs[0], _det],
          [S0.N, S0.om, S0.V, sp.Symbol('X')], {S0.k: sp.Rational(1, 3)})       # = driver_v2.tau_tape(1/3)

def picard_t(tape, box):
    """interval jets (orders 0..CAP-1) of the autonomous ODE y' = tape(y) through every point of box."""
    n = len(box)
    Y = [S(b) for b in box]
    for _ in range(CAP + 1):
        Fs = tape.series_eval(Y, CAP)
        Y = [S(box[i]) + Fs[i].integral() for i in range(n)]
    return Y, tape.series_eval(Y, CAP)

def picard(box3):
    return picard_t(FT, box3)

def taylor_cell(tape, start, rough, h):
    """enclosure of the trajectory over the step: Taylor polynomial at the (thin) start hull of degree
    CAP-2 + Lagrange remainder with the order CAP-1 jet evaluated on the a-priori (rough) box."""
    Ys, _ = picard_t(tape, start)
    Yr, _ = picard_t(tape, rough)
    tau = h/2 + sym(h/2)
    P = CAP - 1
    out = []
    for i in range(len(start)):
        cs = Ys[i].coeffs(); cr = Yr[i].coeffs()
        c = lambda L, k: L[k] if k < len(L) else arb(0)
        s = sum((c(cs, k)*tau**k for k in range(P)), arb(0)) + c(cr, P)*tau**P
        out.append(s.intersection(rough[i]))
    return out

def functions(Y, Fx, Xd):
    """all coefficient functions as series in the local parameter; x-derivative D(G) = G'/Xd.
    Regular forms only (no division by l_- = |lambda_-|), so this also works on the sonic piece."""
    D = lambda G: G.derivative()/Xd
    g = GT.series_eval(Y + [Fx[1], Fx[2], S(CS)], CAP)
    Ms, Mx, P0, I0 = g[0:4], g[4:8], g[8:24], g[24:28]
    A, lp, lm = g[28], g[29], g[30]
    ell = -lm                                   # l_- = |lambda_-|  (lambda_- < 0 on x < 0)
    Msi = ser_inv2(Ms)
    L0 = mm(Msi, Mx)
    sp_ = (lp - L0[0])/L0[1]; sm_ = (lm - L0[0])/L0[1]
    R = [S(arb(1)), S(arb(1)), sp_, sm_]
    Lm = ser_inv2(R)
    Eu = [P0[10], P0[11], P0[14], P0[15]]; Em = [P0[8], P0[9], P0[12], P0[13]]
    B0 = mm(mm(Lm, mm(Msi, Eu)), R)
    em = mm(Lm, mm(Msi, Em))
    dL = [D(t) for t in Lm]
    dLR = mm(dL, R)
    b = [B0[0] + lp*dLR[0], B0[1] + lp*dLR[1], B0[2] + lm*dLR[2], B0[3] + lm*dLR[3]]
    lpl = lp + ell
    bp = b[0]/lp
    sig = b[2]*lp/lpl                           # sigma-hat = b_-+/(l_- T)
    dsig = D(sig)
    invT = lp*ell/lpl                           # 1/T
    chi = (dsig + bp*sig)*invT + b[3]*sig*lp/lpl
    dchi = D(chi)
    ellphi = ell*dchi + ell*bp*chi + b[3]*chi   # l_- phi
    cN, cw, cV = I0[1], I0[2], I0[3]
    al_p = cw + cV*sp_; al_m = cw + cV*sm_
    ga_p = P0[2] + P0[3]*sp_; ga_m = P0[2] + P0[3]*sm_
    H3 = -Y[1]*sp.Rational(2, 3) if False else Y[1]*arb(-2)/3
    eAm, eNm = em[2], em[3]
    q = eAm*al_p*lp/lpl                         # e_-A al_+/(l_- T)
    out = dict(lp=lp, ell=ell, b00=b[0], b01=b[1], b10=b[2], b11=b[3], bp=bp, sig=sig, dsig=dsig, chi=chi, dchi=dchi,
               ellphi=ellphi, eAp=em[0], eNp=em[1], eAm=eAm, eNm=eNm, deAm=D(eAm), deNm=D(eNm),
               al_p=al_p, al_m=al_m, cN=cN, ga_p=ga_p, ga_m=ga_m, A=A, dA=D(A), H3=H3, H3lp=H3*lp, dH3lp=D(H3*lp),
               dAalp=D(A*al_p*lp), dAalm=D(A*al_m*lp), Aalmlp=A*al_m*lp, q=q, dq=D(q), eAcN=eAm*cN, deAcN=D(eAm*cN),
               eAalm=eAm*al_m, L00=Lm[0], L01=Lm[1], L10=Lm[2], L11=Lm[3], dlm=D(lm), invlp=1/lp, lpl=lpl)
    return out

def tm_eval(fs, fb, tau):
    """cell value enclosure: sum_{k<PTM} c_k(start) tau^k + c_PTM(cell box) tau^PTM (Lagrange form)."""
    res = {}
    for key in fs:
        a = fs[key].coeffs(); bb = fb[key].coeffs()
        c = lambda L, k: L[k] if k < len(L) else arb(0)
        res[key] = sum((c(a, k)*tau**k for k in range(PTM)), arb(0)) + c(bb, PTM)*tau**PTM
    return res

# ------------------------------------------------------------------ background cover
BG = pickle.load(open('<LAB>/engine/cert/lin_fw4_bg.pkl', 'rb'))
def A_(p): return arb(p[0]) + sym(arb(p[1]))
MARGIN = A_(BG['margin']).abs_upper() + arb('1e-30')
XSW = A_(BG['x_switch']); XBAR = A_(BG['xbar'])

def tube_cells():
    """cells ordered from the sonic side: dict(xlo, xhi, h, box).  tau leg: x interval between consecutive
    start points (x monotone on the tau leg: PE2 det sign certified); x leg: [x_k + h_k, x_k]."""
    cells = []
    tr = BG['tau_rec']
    xs = [A_(st[3]) for (_, _, _, st) in tr] + [XSW]
    for k, (_, h, Bx, st) in enumerate(tr):
        infl = MARGIN + (arb('1e-8') if k == len(tr) - 1 else 0)   # last: landing gap at x_switch (L dmax <= 1e-10)
        rough = [A_(Bx[i]) + sym(infl) for i in range(4)]
        start = [A_(st[i]) + sym(infl) for i in range(4)]
        xhi, xlo = xs[k], xs[k + 1]
        cells.append(dict(xlo=xlo, xhi=xhi, h=(xhi - xlo), rough=rough, start=start, ht=A_(h), kind='tau'))
    x = XSW
    for (_, h, Bx, st) in BG['rec']:
        ha = A_(h)
        rough = [A_(Bx[i]) + sym(MARGIN) for i in range(3)]
        start = [A_(st[i]) + sym(MARGIN) for i in range(3)]
        cells.append(dict(xlo=x + ha, xhi=x, h=-ha, rough=rough, start=start, ht=ha, kind='x'))
        x = x + ha
    return cells

def _fser(kind, box):
    if kind == 'x':
        Y, Fx = picard(box)
        return functions(Y, Fx, S(arb(1)))
    Y4, _ = picard_t(TT, box)
    Xd = Y4[3].derivative()
    return functions(Y4[:3], [Y4[i].derivative()/Xd for i in range(3)], Xd)

NSUB = 8
def split_tau(c):
    """tau step -> NSUB sub-cells in tau-time, each with a thin start (Taylor + Lagrange on the step box) and
    tau-parametrised function Taylor models (regular field; x-jets are useless next to the sonic point)."""
    box = taylor_cell(TT, c['start'], c['rough'], c['ht'])
    Ys, _ = picard_t(TT, c['start']); Yr, _ = picard_t(TT, box)
    P = CAP - 1
    hs = c['ht']/NSUB
    def at(tj):
        out = []
        for i in range(4):
            cs = Ys[i].coeffs(); cr = Yr[i].coeffs()
            g = lambda L, k: L[k] if k < len(L) else arb(0)
            tb = tj/2 + sym(tj/2)
            out.append(sum((g(cs, k)*tj**k for k in range(P)), arb(0)) + g(cr, P)*tb**P)
        return out
    starts = [c['start']] + [at(hs*j) for j in range(1, NSUB)]
    subs = []
    for j in range(NSUB):
        sb = taylor_cell(TT, starts[j], box, hs)
        fs = _fser('tau', starts[j]); fb = _fser('tau', sb)
        xhi = c['xhi'] if j == 0 else starts[j][3]
        xlo = c['xlo'] if j == NSUB - 1 else starts[j + 1][3]
        subs.append(dict(xlo=xlo, xhi=xhi, h=xhi - xlo, kind='tau', box=sb[:3], xbox=sb[3],
                         f=tm_eval(fs, fb, hs/2 + sym(hs/2))))
    return subs

def eval_cells(cells, log_every=0):
    """x-leg cells in place; tau cells are replaced by their NSUB sub-cells (returns the new list)."""
    out = []
    for j, c in enumerate(cells):
        if c['kind'] == 'tau':
            out += split_tau(c); continue
        out.append(c)
        tape = FT if c['kind'] == 'x' else TT
        box = taylor_cell(tape, c['start'], c['rough'], c['ht'])
        # function Taylor models in x on [xlo, xhi], expanded at the thin start point (x = xhi):
        fs = _fser('x', c['start'][:3])
        fb = _fser('x', box[:3])
        hx = c['xlo'] - c['xhi']
        c['f'] = tm_eval(fs, fb, hx/2 + sym(hx/2))
        c['box'] = box[:3]
        if log_every and j % log_every == 0:
            log('  cell %d/%d x=%.4f' % (j, len(cells), float(c['xhi'].mid())))
    return out

# ------------------------------------------------------------------ sonic piece theta in [th1, 0]
def sonic_piece():
    sm = BG['sm']; nu = arb(BG['nu_s']); th1 = A_(BG['th1'])
    Th = th1/2 + sym(th1/2)                       # the theta box [th1, 0]
    rho = nu - th1.abs_upper()                    # Cauchy radius for the tails
    def jets(co, eps):
        c = [A_(a) for a in co]
        out = []
        for k in range(CAP):
            s = arb(0)
            for n in range(k, len(c)):
                s += c[n]*arb(math.comb(n, k))*Th**(n - k)
            out.append(s + sym(arb(eps)/rho**k))      # |R^(k)/k!| <= eps/rho^k  (Cauchy on |theta| <= nu)
        return arb_series(out, CAP)
    def jets0(co, eps):                           # expansion at theta = 0 (thin), Cauchy tails radius nu
        c = [A_(a) for a in co]
        return arb_series([(c[k] if k < len(c) else arb(0)) + sym(arb(eps)/nu**k) for k in range(CAP)], CAP)
    def fset(jf):
        Y = [jf(sm['Y'][i], arb(sm['eps'][0])) for i in range(3)]
        X = jf(sm['X'], arb(sm['epsX'][0]))
        Xd = X.derivative()
        return functions(Y, [Yi.derivative()/Xd for Yi in Y], Xd), X
    fs, _ = fset(jets0)
    fb, X = fset(jets)
    f = tm_eval(fs, fb, Th)
    s0 = (X.coeffs()[0]).abs_upper()              # |x| <= s0 on the piece (X(0) = 0)
    # exact: x(th1) is the tau-leg start
    s0 = max(s0, A_(BG['tau_rec'][0][3][3]).abs_upper())
    return dict(f=f, s0=s0, Th=Th)

# ================================================================== [B] bulk bootstrap in ball arithmetic
def up(x):
    return arb(x.abs_upper())

def udict(f, sonic=None):
    """nonnegative upper bounds of every bulk coefficient on one cell.  sonic = (s0, lam'_lo): singular
    entries are returned as sup |s| f (|s| <= s0, l_- >= lam'_lo |s|)."""
    if sonic is None:
        ell_lo = f['ell'].lower()
        assert ell_lo > 0, 'l_- not bounded away from 0 on a tube cell'
        ie = 1/arb(ell_lo)
        sing = lambda R: up(R)*ie
        bpb = up(f['bp'] + f['b11']/f['ell'])
    else:
        s0, lam_lo = sonic
        sing = lambda R: up(R)/lam_lo
        bpb = s0*up(f['bp']) + sing(f['b11'])
    return dict(invlp=up(f['invlp']), lp=up(f['lp']), sig=up(f['sig']), dsig=up(f['dsig']), chi=up(f['chi']),
                dchi=up(f['dchi']), phi=sing(f['ellphi']), bpb=bpb, bm=sing(f['b11']), tau=up(f['b01']*f['invlp']),
                eAp=up(f['eAp']), eNp=up(f['eNp']), eAm=up(f['eAm']), eNm=up(f['eNm']), deAm=up(f['deAm']), deNm=up(f['deNm']),
                alp=up(f['al_p']), alm=up(f['al_m']), cN=up(f['cN']), gap=up(f['ga_p']), gam=up(f['ga_m']), A=up(f['A']),
                dA=up(f['dA']), H3=up(f['H3']), H3lp=up(f['H3lp']), dH3lp=up(f['dH3lp']), bpp=up(f['b00']), bpm=up(f['b01']),
                dAalp=up(f['dAalp']), dAalm=up(f['dAalm']), Aalmlp=up(f['Aalmlp']), aqd=up(f['q']), dqd=up(f['dq']),
                eAcN=up(f['eAcN']), deAcN=up(f['deAcN']), aeAalm=sing(f['eAalm']))

Z3 = (arb(0), arb(0), arb(0))
def vadd(*vs): return tuple(sum((v[i] for v in vs), arb(0)) for i in range(3))
def vsc(a, v): return tuple(a*v[i] for i in range(3))
def vdot(c, v): return sum((c[i]*v[i] for i in range(3)), arb(0))
E1, E2, E3 = (arb(1), arb(0), arb(0)), (arb(0), arb(1), arb(0)), (arb(0), arb(0), arb(1))

def integrands(u, kp, rz, nu, w=arb(1)):
    """row 1 (pint), row 2 Gamma-integrand split (S: singular part, B: bounded part), row 2 pointwise,
    row 3 P+-integrand, row 3 pointwise -- all as (U, Z, N)-coefficient vectors per unit |E_+(s, x_c)|."""
    iK, iK2, iKA, Ku = kp['iK'], kp['iK2'], kp['iKA'], kp['Ku']
    D0 = u['aqd']*iK*iKA
    D1 = (u['dqd']*iKA + u['aqd']*u['dA']*iKA**2)*iK
    P1 = u['sig']*iK + u['chi']*iK2 + D0
    dP = u['dsig']*iK + u['dchi']*iK2 + D1
    EmN = u['eNm'] + u['eAcN']*iKA
    dEmN = u['deNm'] + u['deAcN']*iKA + u['eAcN']*u['dA']*iKA**2
    wm = (P1, rz, arb(0))
    Ab = ((u['alp'] + u['alm']*P1)*iKA, u['alm']*rz*iKA, u['cN']*nu*iK*iKA)
    Nb = (arb(0), arb(0), nu*iK)
    Nbp = vadd(vsc(u['A'], Ab), vsc(u['H3'], vadd(E1, wm)))
    gp = vadd(vsc(u['eAp'], Ab), vsc(u['eNp'], Nb))
    r1 = vadd(vsc(u['tau'], wm), vsc(u['invlp'], gp))
    S2 = vadd(vsc(u['phi']*iK2 + u['bpb']*D0, E1), vsc(u['aeAalm']*iKA, wm), vsc(iK*u['bm']*EmN, Nb))
    B2 = vadd(vsc(D1, E1), vsc(P1, vadd(vsc(u['tau'], wm), vsc(u['invlp'], gp))), vsc(iK*dEmN, Nb), vsc(iK*EmN, Nbp))
    p2 = vsc(EmN*iK, Nb)
    Fw = u['H3']*(1 + P1) + u['A']*(u['alp'] + u['alm']*P1)*iKA
    dFw = u['dH3lp']*(1 + P1) + u['H3lp']*dP + (u['dAalp'] + u['dAalm']*P1 + u['Aalmlp']*dP)*iKA \
        + u['A']*(u['alp'] + u['alm']*P1)*u['lp']*u['dA']*iKA**2
    Fz = u['H3'] + u['A']*u['alm']*iKA
    Fn = u['A']*u['cN']*iKA
    I3 = vadd(vsc(Fw*u['bpp'] + dFw, E1), vsc(Fw*u['bpm'], wm), vsc(Fw, gp), vsc(Ku*Fz*rz, E2), vsc(Fn*nu, E3))
    p3 = vsc(Fw*u['lp'], E1)
    return dict(r1=r1, S2=S2, B2=B2, I2=vadd(S2, B2), p2=p2, I3=I3, p3=p3, Fwlp=Fw*u['lp'])

def bulk_ia(cells, son, M, Ka, Kb, weights=None, iters=6):
    """rigorous bulk outputs for every K in [Ka, Kb] (K = |kappa|, Re kappa >= 0)."""
    s0, g, eSP = son['s0'], son['g'], son['eSP']
    xc_lo = (arb(M)/Kb).log(); xc_hi = (arb(M)/Ka).log()
    rng = [i for i, c in enumerate(cells) if c['xhi'].upper() >= xc_lo.lower()]
    n = len(rng)
    out_idx = [i for i in rng if cells[i]['xlo'].lower() <= xc_hi.upper()]
    kp = dict(iK=1/Ka, iK2=1/Ka**2, iKA=1/(Ka - son['Amax']), Ku=Kb)
    if weights is None:
        Kf = float(Ka.mid())
        rz = [arb(1/(Kf*math.exp(float(cells[i]['xhi'].mid())))**2 + 1/Kf**2) for i in rng]
        nu = [arb(float(cells[i]['u']['H3lp'].mid()) + 1/Kf) for i in rng]
        rzS, nuS = arb(1/Kf**2), arb(float(son['u']['H3lp'].mid()) + 1/Kf)
    else:
        rz, nu, rzS, nuS = weights
    for it in range(iters + 1):
        # sonic piece (integrals over [-s0, 0] and t inside the piece)
        uS = son['u']
        J = integrands(uS, kp, rzS, nuS)
        G2 = vsc(eSP, vadd(vsc(1/g, J['S2']), vsc(s0, J['B2'])))
        G3 = vsc(eSP*s0, J['I3'])
        Wb = eSP
        Fw0lp0 = J['Fwlp']
        R1 = vsc(s0, J['r1'])
        c2 = [vadd(G2, J['p2'])]; c3 = [vadd(G3, J['p3'], vsc(Fw0lp0*eSP, E1))]
        wts = [(rzS, nuS)]
        for jj, i in enumerate(rng):
            c = cells[i]; u = c['u']; h = c['hU']
            Ji = integrands(u, kp, rz[jj], nu[jj])
            EG = (c['h']*c['bG']).upper().exp(); OG = max(arb(0), (c['h']*c['bG']).upper()).exp()
            EP = (c['h']*c['bP']).upper().exp(); OP = max(arb(0), (c['h']*c['bP']).upper()).exp()
            c2.append(vadd(vsc(OG, vadd(vsc(h, Ji['I2']), G2)), Ji['p2']))
            c3.append(vadd(vsc(OP, vadd(vsc(h, Ji['I3']), G3)), Ji['p3'], vsc(Fw0lp0*OP*Wb, E1)))
            G2 = vadd(vsc(OG*h, Ji['I2']), vsc(EG, G2))
            G3 = vadd(vsc(OP*h, Ji['I3']), vsc(EP, G3))
            Wb = EP*Wb
            R1 = vadd(R1, vsc(h, Ji['r1']))
            wts.append((rz[jj], nu[jj]))
        A = [[up(R1[k]) for k in range(3)],
             [max(up(c2[j][k]/wts[j][0]) for j in range(len(c2))) for k in range(3)],
             [max(up(c3[j][k]/wts[j][1]) for j in range(len(c3))) for k in range(3)]]
        Af = np.array([[float(a.mid()) for a in row] for row in A])
        ev = np.linalg.eigvals(Af); rho = max(abs(ev))
        if rho < 1:
            vf = np.linalg.solve(np.eye(3) - Af, np.array([1.0, 0, 0]))
        else:
            vf = np.array([1.0, 1.0, 1.0])
        if it < iters:            # weights := current profiles (any fixed positive weights are admissible)
            rzS = arb(max(float(vdot(c2[0], vf).mid()), 1e-300)); nuS = arb(max(float(vdot(c3[0], vf).mid()), 1e-300))
            rz = [arb(max(float(vdot(c2[j + 1], vf).mid()), 1e-300)) for j in range(n)]
            nu = [arb(max(float(vdot(c3[j + 1], vf).mid()), 1e-300)) for j in range(n)]
    # ---- rigorous closure: Collatz bound with the float Perron vector, then (I - A)^{-1} e1
    w_, V_ = np.linalg.eig(Af); pv = np.abs(V_[:, np.argmax(abs(w_))].real) + 1e-6
    y = [arb(float(t)) for t in pv]
    coll = max(up(sum((A[i][j]*y[j] for j in range(3)), arb(0))/y[i]) for i in range(3))
    res = dict(A=A, rho_float=rho, collatz=coll, ok=bool(coll < 1), weights=(rz, nu, rzS, nuS), n=n)
    if not res['ok']:
        return res
    I_A = [[(1 if i == j else 0) - A[i][j] for j in range(3)] for i in range(3)]
    det = (I_A[0][0]*(I_A[1][1]*I_A[2][2] - I_A[1][2]*I_A[2][1]) - I_A[0][1]*(I_A[1][0]*I_A[2][2] - I_A[1][2]*I_A[2][0])
           + I_A[0][2]*(I_A[1][0]*I_A[2][1] - I_A[1][1]*I_A[2][0]))
    v = (up((I_A[1][1]*I_A[2][2] - I_A[1][2]*I_A[2][1])/det), up(-(I_A[1][0]*I_A[2][2] - I_A[1][2]*I_A[2][0])/det),
         up((I_A[1][0]*I_A[2][1] - I_A[1][1]*I_A[2][0])/det))
    pos = {i: j + 1 for j, i in enumerate(rng)}
    eb = max(up(vdot(c2[pos[i]], v)) for i in out_idx)
    nbar = max(up(vdot(c3[pos[i]], v)) for i in out_idx)
    def hull(key):
        h_ = None
        for i in out_idx:
            a = cells[i]['f'][key]; h_ = a if h_ is None else h_.union(a)
        return h_
    res.update(v=v, eb=eb, nbar=nbar, sig_c=hull('sig'), chi_c=hull('chi'), q_c=hull('q'), A_c=hull('A'),
               L_c=[hull('L00'), hull('L01'), hull('L10'), hull('L11')], out_idx=out_idx)
    return res

def prep_bulk(cells, son):
    for c in cells:
        c['u'] = udict(c['f'])
        c['hU'] = up(c['h'])
        c['bP'] = c['f']['bp']; c['bG'] = c['f']['bp'] + c['f']['b11']/c['f']['ell']
    f = son['f']
    lam_lo = arb(f['dlm'].lower()); lam_hi = arb(f['dlm'].upper())
    assert lam_lo > 0 and f['b11'].upper() < 0, 'sonic piece: need dlambda_-/dx > 0 and b_-- < 0'
    son['g'] = arb((-f['b11']).lower())/lam_hi
    son['u'] = udict(f, sonic=(son['s0'], lam_lo))
    son['eSP'] = (son['s0']*max(arb(0), arb(f['bp'].upper()))).exp()
    son['Amax'] = max([up(c['f']['A']) for c in cells] + [up(f['A'])])
    return cells, son

# ================================================================== [C] centre zone (uniform in K >= Ka)
_fl = pickle.load(open('<LAB>/engine/cert/lin_fields.pkl', 'rb'))
_Mc = sp.Matrix(4, 4, [sp.sympify(s) for s in _fl['Mc']]); _kap = sp.Symbol('kappa')
MT = Tape(list(_Mc.subs(_kap, 0)) + list(_Mc.diff(_kap)), [S0.P, S0.W, S0.U], {S0.k: sp.Rational(1, 3)})
assert all(e.free_symbols <= {S0.P, S0.W, S0.U} for e in _Mc.diff(_kap)), 'Mc must be affine in kappa'
C00 = [[-1, 0, 0, 0], [1, 0, 0, 0], [-3, 0, 0, 0], [0, 0, 0, -2]]
C11 = [[0]*4 for _ in range(4)]; C11[2][3] = arb(-4); C11[3][2] = arb(-3)/4
XC = A_(BG['Xb'][2]); PS = (-XC).exp()                      # p = e^{-xc}  (v = p e^x = theta e^{chi})
R3 = arb(3).sqrt()

def mc_series(Y3, cap):
    """C0, C1 (4x4 lists of series) from background series (N, om, V) -> (P, W, U)."""
    N, om, V = Y3
    P = 1/N; W = om*N*N; U = V*N
    o = MT.series_eval([P, W, U], cap)
    return [[o[4*i + j] for j in range(4)] for i in range(4)], [[o[16 + 4*i + j] for j in range(4)] for i in range(4)]

def centre_manifold_data(nint=10, cap=24, P=22):
    """C01, C12 and the running-sup data on theta in [0, thm] (Taylor models from theta = 0)."""
    cm = BG['cm']; nu = arb(BG['nu_c']); thm = A_(BG['thm'])
    epsY = arb(cm['eps'][0]); epsX = arb(cm['epsX'][0])
    Yc = [[A_(a) for a in cm['Y'][i]] for i in range(3)]; Xc = [A_(a) for a in cm['X']]
    assert Yc[0][0].is_zero() and (Yc[0][1] - 1).is_zero(), 'centre manifold: need P = theta + O(theta^2) exactly'
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
    def gser(Y, X):
        C0, C1 = mc_series(Y, cap)
        th = arb_series([arb(0), arb(1)], cap)
        v = th*X.exp()
        return C0, C1, v
    # thin series at theta = 0: Taylor coefficients
    Y0 = [j0(Yc[i], epsY) for i in range(3)]
    # (P, W, U) are the manifold variables themselves here: Y = Z; feed MT directly
    def mcZ(Z, cap_):
        o = MT.series_eval(Z, cap_)
        return [[o[4*i + j] for j in range(4)] for i in range(4)], [[o[16 + 4*i + j] for j in range(4)] for i in range(4)]
    th = arb_series([arb(0), arb(1)], cap)
    X0 = j0(Xc, epsX)
    C0s, C1s = mcZ(Y0, cap); v0 = th*X0.exp()
    vc = v0.coeffs()
    cf = lambda s_, k: s_.coeffs()[k] if k < len(s_.coeffs()) else arb(0)
    C01 = [[cf(C0s[i][j], 1) for j in range(4)] for i in range(4)]
    C12 = [[cf(C1s[i][j], 2) - C11[i][j]*cf(v0, 2) for j in range(4)] for i in range(4)]
    # exact facts (PROVED, lin_sym / Mc algebra): C0(0) = C00, C1(0) = 0, C1'(0) = C11  ->  zero coefficients
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
    data = []
    Tb = [thm*arb(2)**(-(nint - k)) for k in range(1, nint + 1)]
    lo = arb(0)
    for hi in Tb:
        Th = hi/2 + sym(hi/2)                       # [0, hi] (running sup from 0)
        rho = nu - hi.abs_upper()
        Yb = [jb(Yc[i], epsY, Th, rho) for i in range(3)]
        Xb = jb(Xc, epsX, Th, rho)
        C0b, C1b = mcZ(Yb, cap)
        vb = th*Xb.exp()
        chiTh = sum((Xc[n_]*Th**n_ for n_ in range(len(Xc))), arb(0)) + sym(epsX)
        ent = {}
        for which in ('G0', 'G1', 'R0', 'R1'):
            A0 = G_series(C0s, C1s, v0, which); Ab = G_series(C0b, C1b, vb, which)
            M_ = [[arb(0)]*4 for _ in range(4)]
            for i in range(4):
                for j in range(4):
                    g0, m = A0[i][j]; gb, _ = Ab[i][j]
                    a = g0.coeffs(); b = gb.coeffs()
                    c = lambda L, k: L[k] if k < len(L) else arb(0)
                    # H = G/theta^m: coefficients below m are exactly zero (PROVED); Lagrange remainder at order P
                    H = sum((c(a, k)*Th**(k - m) for k in range(m, P)), arb(0)) + c(b, P)*Th**(P - m)
                    M_[i][j] = up(H*(-m*chiTh).exp())
            ent[which] = M_
        vhi = up(hi*chiTh.exp())
        data.append(dict(vlo=arb(0), vhi=vhi, **ent))
    return dict(C01=C01, C12=C12, data=data, thm=thm)

def tube_centre_data(cells, C01, C12, v_need):
    """same four ratios on tube cells with v_lo <= v_need (x-Taylor models, thin start + cell box)."""
    out = []
    for c in cells:
        if c['kind'] != 'x': continue
        vlo = PS*c['xlo'].exp(); vhi = PS*c['xhi'].exp()
        if vlo.lower() > v_need.upper(): continue
        box = c['box']
        res = {}
        for which in ('G0', 'G1', 'R0', 'R1'):
            vals = []
            for Y3 in (c['start'], box):
                Y, _ = picard(Y3)
                C0, C1 = mc_series(Y, CAP)
                v = (PS*c['xhi'].exp())*arb_series([arb(0), arb(1)], CAP).exp()
                M_ = [[None]*4 for _ in range(4)]
                for i in range(4):
                    for j in range(4):
                        if which == 'G0': g = C0[i][j] - C00[i][j]
                        elif which == 'G1': g = C1[i][j] - v*C11[i][j]
                        elif which == 'R0': g = C0[i][j] - C00[i][j] - v*C01[i][j]
                        else: g = C1[i][j] - v*C11[i][j] - v*v*C12[i][j]
                        M_[i][j] = g
                vals.append(M_)
            m = {'G0': 1, 'G1': 2, 'R0': 2, 'R1': 3}[which]
            hx = c['xlo'] - c['xhi']; tau = hx/2 + sym(hx/2)
            Mo = [[up(tm_eval({'g': vals[0][i][j]}, {'g': vals[1][i][j]}, tau)['g'])/arb(vlo.lower())**m
                   for j in range(4)] for i in range(4)]
            res[which] = Mo
        out.append(dict(vlo=vlo, vhi=vhi, **res))
    return out

def bt_centre_data(C01, C12):
    """the transport gap [x_m, xbar]: box Bt (N, om, V), v between v(x_m) and v(xbar)."""
    Bt = [A_(b) for b in BG['Bt']]
    N, om, V = Bt
    o = MT.eval([1/N, om*N*N, V*N])
    C0 = [[o[4*i + j] for j in range(4)] for i in range(4)]; C1 = [[o[16 + 4*i + j] for j in range(4)] for i in range(4)]
    cm = BG['cm']; thm = A_(BG['thm'])
    chi_m = sum((A_(a)*thm**n_ for n_, a in enumerate(cm['X'])), arb(0)) + sym(arb(cm['epsX'][0]))
    v = (thm*chi_m.exp()).union(PS*XBAR.exp())
    res = {}
    for which, m in (('G0', 1), ('G1', 2), ('R0', 2), ('R1', 3)):
        M_ = [[None]*4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                if which == 'G0': g = C0[i][j] - C00[i][j]
                elif which == 'G1': g = C1[i][j] - v*C11[i][j]
                elif which == 'R0': g = C0[i][j] - C00[i][j] - v*C01[i][j]
                else: g = C1[i][j] - v*C11[i][j] - v*v*C12[i][j]
                M_[i][j] = up(g)/arb(v.lower())**m
        res[which] = M_
    return dict(vlo=v, vhi=v, **res)

def running_sup(data, vmax):
    """entrywise sup of each ratio over all cells with u <= vmax."""
    out = {w: [[arb(0)]*4 for _ in range(4)] for w in ('G0', 'G1', 'R0', 'R1')}
    for d in data:
        lo = d['vlo'] if not hasattr(d['vlo'], 'lower') else arb(d['vlo'].lower())
        if lo > vmax: continue
        for w in out:
            for i in range(4):
                for j in range(4):
                    out[w][i][j] = max(out[w][i][j], d[w][i][j])
    return out

# ---- flat model, closed form (PROOF-B4-v4 sec. 4.2): entire functions of d, no cancellation
def _ser_fn(d, kind, nterm=40):
    if arb(d.abs_lower()) >= 1:                  # direct closed forms (no cancellation for |d| >= 1)
        if kind == 'sh1': return d.sinh()/d
        if kind == 'ch': return d.cosh()
        if kind == 'g3': return (d*d.cosh() - d.sinh())/d**3
        if kind == 'c2': return (d.cosh() - 1)/d**2
    w = d*d
    s = acb(0); p = acb(1)
    fac = {'sh1': lambda k: 1/arb.fac_ui(2*k + 1), 'ch': lambda k: 1/arb.fac_ui(2*k),
           'g3': lambda k: arb(2*(k + 1))/arb.fac_ui(2*k + 3), 'c2': lambda k: 1/arb.fac_ui(2*k + 2)}[kind]
    for k in range(nterm):
        s += fac(k)*p; p *= w
    aw = w.abs_upper()
    tail = arb(aw)**nterm/arb.fac_ui(2*nterm)*arb(2)          # |sum_{k>=n} w^k c_k| <= |w|^n/(2n)! * 2 for |w| <= n
    assert aw < nterm
    return s + acb(sym(tail), sym(tail))

def Kf_ball(zx, rho):
    """4x4 acb ball of K_f(x, y) for zeta_x in zx, rho = t_y/t_x in rho (real ball in [0, 1])."""
    d = (1 - rho)*zx
    sh, ch, g3, c2 = _ser_fn(d, 'sh1'), _ser_fn(d, 'ch'), _ser_fn(d, 'g3'), _ser_fn(d, 'c2')
    o = acb(0); K = [[o]*4 for _ in range(4)]
    K[0][0] = acb(rho); K[1][0] = acb(1 - rho); K[1][1] = acb(1)
    K[2][0] = -3*(1 - rho)*sh
    K[2][2] = (1 - rho)*sh + rho*ch
    K[2][3] = -(4/R3)*rho*d*sh
    K[3][0] = (3*R3/4)*zx*(rho*(1 - rho)**2*c2 + (1 - rho)**3*g3)
    K[3][2] = -(R3/4)*zx*((1 - rho)**3*g3 + rho*(1 - rho)*sh)
    K[3][3] = rho**2*ch + rho*(1 - rho)**3*zx*zx*g3
    return K

def _rser(x, kind, nterm=60):
    """real positive series value at x >= 0 (arb), with tail bound: sh1, ch, g3, c2 at |d| = x."""
    w = x*x; s = arb(0); p = arb(1)
    fac = {'sh1': lambda k: 1/arb.fac_ui(2*k + 1), 'ch': lambda k: 1/arb.fac_ui(2*k),
           'g3': lambda k: arb(2*(k + 1))/arb.fac_ui(2*k + 3), 'c2': lambda k: 1/arb.fac_ui(2*k + 2)}[kind]
    for k in range(nterm):
        s += fac(k)*p; p *= w
    assert w.upper() < nterm
    return s + sym(arb(w.upper())**nterm/arb.fac_ui(2*nterm)*2)

def Kabs_bound(tx, rho, cth):
    """entrywise upper bounds of |K_f(x, y)|, t_x in tx, rho = t_y/t_x in rho (in [0,1]), cos(theta) in cth.
    |sinh d|, |cosh d| <= cosh(Re d); |sh1(d)| <= sh1(|d|); |g3(d)| <= g3(|d|); |c2(d)| <= c2(|d|) (positive series)."""
    om = 1 - rho
    ad = arb((om*tx).upper())                    # |d| upper
    adl = arb(max(arb(0), om*tx).lower()) if (om*tx).lower() > 0 else arb(0)
    a = arb((ad*arb(cth.abs_upper())).upper())      # |Re d| <= |d| sup|cos theta|  (cosh even)
    cha = a.cosh()
    S1 = _rser(ad, 'sh1')
    if adl > 0: S1 = min(arb(S1.upper()), arb((cha/adl).upper()))
    G3 = _rser(ad, 'g3')
    if adl > 0: G3 = min(arb(G3.upper()), arb(((ad + 1)*cha/adl**3).upper()))
    C2 = _rser(ad, 'c2')
    if adl > 0: C2 = min(arb(C2.upper()), arb(((cha + 1)/adl**2).upper()))
    CH = cha
    r = arb(rho.upper()); o = arb(om.upper()); t = arb(tx.upper())
    K = [[arb(0)]*4 for _ in range(4)]
    K[0][0] = r; K[1][0] = o; K[1][1] = arb(1)
    K[2][0] = 3*o*S1; K[2][2] = o*S1 + r*CH; K[2][3] = (4/R3)*r*ad*S1
    K[3][0] = (3*R3/4)*t*(r*o**2*C2 + o**3*G3)
    K[3][2] = (R3/4)*t*(o**3*G3 + r*o*S1)
    K[3][3] = r**2*CH + r*o**3*t**2*G3
    return [[arb(x.upper()) for x in row] for row in K]

def Yf_ball(z, mode):
    if mode == 'N':
        return [acb(0), acb(1), acb(0), acb(0)]
    return [acb(0), acb(0), _ser_fn(z, 'sh1'), -(R3/4)*z*_ser_fn(z, 'g3')]

def centre_theta(th_lo, th_hi, Ka, M, cdat, ns=150, tmin=arb('1e-3')):
    """uniform in K >= Ka: rho_m(x_c) in Yf(Z e^{i th}) + u B1 +- u^2 E2 (u = 1/K), m = 'w', 'N'.
    Returns dict m -> (Yf, B1 (acb 4-vectors), E2 (arb 4-vector)) for the theta-cell [th_lo, th_hi]."""
    thb = (th_lo + th_hi)/2 + sym((th_hi - th_lo)/2)
    eps = acb(0, thb).exp()
    Zt = R3*PS*M
    s0_ = tmin.log(); s1 = Zt.log()
    ds = (s1 - s0_)/ns
    tg = [tmin] + [(s0_ + ds*k).exp() for k in range(1, ns)] + [Zt]
    cells_t = [(tg[k], tg[k + 1]) for k in range(ns)]
    Ka_ = arb(Ka)
    # Delta-hat per t-cell (entrywise), running sups at u <= t_hi/(sqrt3 Ka)
    def dhat(thi, rem=False):
        rs = running_sup(cdat, thi/(R3*Ka_))
        if not rem:
            return [[thi/R3*rs['G0'][i][j] + thi**2/3*rs['G1'][i][j] for j in range(4)] for i in range(4)]
        return [[thi**2/3*rs['R0'][i][j] + thi**3/(3*R3)*rs['R1'][i][j] for j in range(4)] for i in range(4)]
    DH = [dhat(arb(c[1].upper())) for c in cells_t]
    DR = [dhat(arb(c[1].upper()), rem=True) for c in cells_t]
    rsT = running_sup(cdat, tmin/(R3*Ka_))
    dT = [[tmin/R3*rsT['G0'][i][j] + tmin**2/6*rsT['G1'][i][j] for j in range(4)] for i in range(4)]       # int_0^tmin DH ds
    dTR = [[tmin**2/6*rsT['R0'][i][j] + tmin**3/(9*R3)*rsT['R1'][i][j] for j in range(4)] for i in range(4)]
    tb = [(a + b)/2 + sym((b - a)/2) for (a, b) in cells_t]
    tT = tmin/2 + sym(tmin/2)
    # kernel sups: KA[i][j] (j <= i), KT[i] (y in tail), KTT (x, y in tail)
    cth = thb.cos()
    def kab(tx, rho):
        return Kabs_bound(tx, rho, cth)
    one01 = arb('0.5') + sym(arb('0.5'))
    KA = []
    for i in range(ns):
        row = []
        for j in range(i + 1):
            rho = (tb[j]/tb[i]).intersection(one01)
            row.append(kab(tb[i], rho))
        KA.append(row)
    KT = [kab(tb[i], (tT/tb[i]).intersection(one01)) for i in range(ns)]
    KTT = kab(tT, one01)
    def mv(Mx, v): return [sum((Mx[a][b]*v[b] for b in range(4)), arb(0)) for a in range(4)]
    def mm4(A, B): return [[sum((A[a][c]*B[c][b] for c in range(4)), arb(0)) for b in range(4)] for a in range(4)]
    out = {}
    for mode in ('w', 'N'):
        Yfc = [Yf_ball(tb[k]*eps, mode) for k in range(ns)]
        aY = [[arb(y.abs_upper()) for y in Yk] for Yk in Yfc]
        aYT = [arb(y.abs_upper()) for y in Yf_ball(tT*eps, mode)]
        # majorant fixed point eps (supersolution), all /Ka
        def rhs(e, eT):
            fT = mv(dT, [aYT[b] + eT[b] for b in range(4)])
            eTn = [up(t)/Ka_ for t in mv(KTT, fT)]
            f = [[t*ds for t in mv(DH[j], [aY[j][b] + e[j][b] for b in range(4)])] for j in range(ns)]
            new = []
            for i in range(ns):
                acc = mv(KT[i], fT)
                for j in range(i + 1):
                    Kij = KA[i][j]; fj = f[j]
                    acc = [acc[a] + Kij[a][0]*fj[0] + Kij[a][1]*fj[1] + Kij[a][2]*fj[2] + Kij[a][3]*fj[3] for a in range(4)]
                new.append([up(t)/Ka_ for t in acc])
            return new, eTn
        # non-rigorous float iteration to the fixed point, then ONE rigorous supersolution check in arb
        fl = lambda x: float(x.upper())
        KAf = [[np.array([[fl(t) for t in r_] for r_ in KA[i][j]]) for j in range(i + 1)] for i in range(ns)]
        KTf = [np.array([[fl(t) for t in r_] for r_ in KT[i]]) for i in range(ns)]
        KTTf = np.array([[fl(t) for t in r_] for r_ in KTT]); dTf = np.array([[fl(t) for t in r_] for r_ in dT])
        DHf = [np.array([[fl(t) for t in r_] for r_ in DH[j]]) for j in range(ns)]
        aYf = np.array([[fl(t) for t in r_] for r_ in aY]); aYTf = np.array([fl(t) for t in aYT])
        dsf = fl(ds); Kaf = float(Ka)
        ef = np.zeros((ns, 4)); eTf = np.zeros(4)
        for it in range(200):
            fT = dTf @ (aYTf + eTf)
            eTn = KTTf @ fT/Kaf
            fj = np.array([DHf[j] @ (aYf[j] + ef[j]) for j in range(ns)])*dsf
            en = np.array([KTf[i] @ fT + sum(KAf[i][j] @ fj[j] for j in range(i + 1)) for i in range(ns)])/Kaf
            dif = np.max(np.abs(en - ef)) + np.max(np.abs(eTn - eTf))
            ef, eTf = en, eTn
            if dif < 1e-15*(1 + np.max(ef)): break
        fl0 = 1e-12*(1 + float(np.max(ef)))
        es = [[arb(float(t)*1.002 + fl0) for t in row] for row in ef]; eTs = [arb(float(t)*1.002 + fl0) for t in eTf]
        chk, chkT = rhs(es, eTs)
        sup_ok = all(chk[i][a] <= es[i][a] for i in range(ns) for a in range(4)) and all(chkT[a] <= eTs[a] for a in range(4))
        if not sup_ok:
            worst = max(float((chk[i][a]/es[i][a]).upper()) for i in range(ns) for a in range(4) if es[i][a] > 0)
            log('   supersolution check failed (%s): max rhs/eps = %.6f, iterations %d' % (mode, worst, it))
            out[mode] = None; continue
        # E2 at x_c (last cell): K_R |Yf| + Ka K|Delta| eps   (both with the 1/K^2 scaling made explicit)
        last = ns - 1
        acc = mv(KT[last], [mv(dTR, aYT)[a] + Ka_*mv(dT, eTs)[a] for a in range(4)])
        for j in range(ns):
            f = [mv(DR[j], aY[j])[a] + Ka_*mv(DH[j], es[j])[a] for a in range(4)]
            acc = [acc[a] + mv(KA[last][j], f)[a]*ds for a in range(4)]
        E2 = [up(t) for t in acc]
        if os.environ.get('FW_DEBUG'):
            accR = mv(KT[last], mv(dTR, aYT))
            for j in range(ns):
                accR = [accR[a] + mv(KA[last][j], mv(DR[j], aY[j]))[a]*ds for a in range(4)]
            log('   debug %s: eps(x_c)*Ka=%s  K_R|Yf|=%s  |Yf(Z)|=%s' % (mode, [float(t.mid())*float(Ka) for t in es[-1]],
                [float(t.mid()) for t in accR], [float(t.mid()) for t in aY[-1]]))
            fj = [[float(t.mid()) for t in mv(DH[j], [aY[j][b] + es[j][b] for b in range(4)])] for j in range(ns)]
            contrib = np.zeros((4, 4))
            for j in range(ns):
                for a_ in range(4):
                    for b_ in range(4):
                        contrib[a_, b_] += float(KA[last][j][a_][b_].mid())*fj[j][b_]*float(ds.mid())
            log('   debug %s: eps*Ka contributions (row a <- source component b):\n%s' % (mode, np.round(contrib, 2)))
            log('   debug %s: es*Ka profile comp: %s' % (mode, [[round(float(es[j][b].mid())*float(Ka), 2) for b in range(4)] for j in range(0, ns, 25)]))
        # B1^(1)(theta) = int K_f(Z, y) L(t) Yf(y) ds,  L = C01 t/sqrt3 + eps C12 t^2/3 ; tail cell with ds = dt/t
        Lm = lambda t: [[acb(cdat_C01[a][b]*t/R3) + eps*acb(cdat_C12[a][b]*t*t/3) for b in range(4)] for a in range(4)]
        B1 = [acb(0)]*4
        zX = Zt*eps
        nb = NB1
        dsb = (s1 - s0_)/nb
        for j in range(nb):
            ta = (s0_ + dsb*j).exp(); tc = (s0_ + dsb*(j + 1)).exp() if j < nb - 1 else Zt
            tj = (ta + tc)/2 + sym((tc - ta)/2)
            K = Kf_ball(zX, (tj/Zt).intersection(one01))
            Lj = Lm(tj); Yj = Yf_ball(tj*eps, mode)
            LY = [sum((Lj[a][b]*Yj[b] for b in range(4)), acb(0)) for a in range(4)]
            KLY = [sum((K[a][b]*LY[b] for b in range(4)), acb(0)) for a in range(4)]
            B1 = [B1[a] + KLY[a]*(tc.log() - ta.log()) for a in range(4)]
        # tail [0, tmin]: |L|/t bounded -> |int| <= tmin * sup|K| sup|L/t| sup|Yf|
        LT = [[arb(cdat_C01[a][b].abs_upper())/R3 + arb(cdat_C12[a][b].abs_upper())*tmin/3 for b in range(4)] for a in range(4)]
        KTZ = kab(Zt, (tT/Zt).intersection(one01))
        tl = mv(KTZ, mv(LT, aYT))
        B1 = [B1[a] + acb(sym(tl[a]*tmin), sym(tl[a]*tmin)) for a in range(4)]
        out[mode] = (Yf_ball(zX, mode), B1, E2)
    return out

cdat_C01 = None; cdat_C12 = None
NB1 = 1500

# ================================================================== [T] contradiction test (mu-disc)
def centre_test_ia(bk, cen, Ka, Kb, th_lo, th_hi):
    thb = (th_lo + th_hi)/2 + sym((th_hi - th_lo)/2)
    eps = acb(0, thb).exp()
    Kb_ = arb(Ka)/2 + arb(Kb)/2 + sym((arb(Kb) - arb(Ka))/2)
    kap = Kb_*eps
    cb = acb(bk['sig_c'])/kap + acb(bk['chi_c'])/kap**2 - acb(bk['q_c'])/(kap*(kap - acb(bk['A_c'])))
    u = 1/Kb_
    def Yball(m):
        Yf, B1, E2 = cen[m]
        return [Yf[a] + u*B1[a] + acb(sym(E2[a]*u.upper()**2), sym(E2[a]*u.upper()**2)) for a in range(4)]
    Yw, YN = Yball('w'), Yball('N')
    L = bk['L_c']
    Aw = [acb(L[0])*Yw[2] + acb(L[1])*Yw[3], acb(L[2])*Yw[2] + acb(L[3])*Yw[3]]
    Bw = [acb(L[0])*YN[2] + acb(L[1])*YN[3], acb(L[2])*YN[2] + acb(L[3])*YN[3]]
    a, bN = Yw[1], YN[1]
    nK = bk['nbar']/arb(Ka)
    den = arb(bN.abs_lower()) - nK*arb(Bw[0].abs_upper())
    if not den > 0:
        return False, 'lapse-unbounded'
    mu0 = -a/bN
    r = nK*(arb(Aw[0].abs_upper()) + arb(mu0.abs_upper())*arb(Bw[0].abs_upper()))/den
    g0 = Aw[1] - cb*Aw[0]; g1 = Bw[1] - cb*Bw[0]
    lower = arb((g0 + mu0*g1).abs_lower()) - r*arb(g1.abs_upper())
    wp = arb((Aw[0] + mu0*Bw[0]).abs_upper()) + r*arb(Bw[0].abs_upper())
    ok = bool(lower > bk['eb']*wp)
    return ok, dict(lower=lower, rhs=bk['eb']*wp, r=r, cb=cb)

def centre_bulk_cells(theta_min, ratio=1.02):
    """bulk cells below xbar from the E5 centre manifold: theta from thm down to theta_min (geometric),
    x = log theta + chi(theta) + xc, (N, om, V) = (1/P, W P^2, U P); theta-Taylor models (thin start at the
    upper end + Lagrange remainder on the theta-box, Cauchy tails radius nu_c - theta).  Plus the transport gap
    [x_m, xbar] as one cell with the box Bt."""
    cm = BG['cm']; nu = arb(BG['nu_c']); thm = A_(BG['thm'])
    epsY = arb(cm['eps'][0]); epsX = arb(cm['epsX'][0])
    Yc = [[A_(a) for a in cm['Y'][i]] for i in range(3)]; Xc = [A_(a) for a in cm['X']]
    def jets(c, eps, Th, rho):
        out = []
        for k in range(CAP):
            s = arb(0)
            for n_ in range(k, len(c)):
                s += c[n_]*arb(math.comb(n_, k))*Th**(n_ - k)
            out.append(s + sym(eps/rho**k))
        return arb_series(out, CAP)
    def fset(Th, rho):
        Z = [jets(Yc[i], epsY, Th, rho) for i in range(3)]
        chi = jets(Xc, epsX, Th, rho)
        xs = arb_series([Th, arb(1)], CAP).log() + chi + XC
        Xd = xs.derivative()
        P, W, U = Z
        Y = [1/P, W*P*P, U*P]
        return functions(Y, [Yi.derivative()/Xd for Yi in Y], Xd), xs
    def xat(th):
        return th.log() + sum((Xc[n_]*th**n_ for n_ in range(len(Xc))), arb(0)) + sym(epsX) + XC
    cells = []
    # transport gap cell
    Bt = [A_(b) for b in BG['Bt']]
    Y, Fx = picard(Bt)
    fg = {k_: (v_.coeffs()[0] if v_.coeffs() else arb(0)) for k_, v_ in functions(Y, Fx, S(arb(1))).items()}
    xm = xat(thm)
    gl, gh = (xm, XBAR) if xm.mid() < XBAR.mid() else (XBAR, xm)
    cells.append(dict(xlo=arb(gl.lower()), xhi=arb(gh.upper()), h=arb(gh.upper()) - arb(gl.lower()), kind='gap', f=fg, box=Bt))
    th = thm
    while th.upper() > theta_min:
        tl = th/ratio
        rho = nu - th.abs_upper()
        fs, _ = fset(th, rho)
        fb, _ = fset((th + tl)/2 + sym((th - tl)/2), rho)
        hs = tl - th
        f = tm_eval(fs, fb, hs/2 + sym(hs/2))
        cells.append(dict(xlo=xat(tl), xhi=xat(th), h=xat(th) - xat(tl), kind='cm', f=f, box=None))
        th = tl
    return cells

def certify(M, Ka, Kmax, nth, kratio=1.03, ns=150, th_sub=None, log_each=False):
    """main certificate: for all kappa with Re kappa >= 0, Ka <= |kappa| <= Kmax: no physical eigenvalue.
    theta in [0, pi/2] by nth equal cells (conjugate symmetry covers [-pi/2, 0])."""
    global cdat_C01, cdat_C12
    T1 = time.time()
    cells = eval_cells(tube_cells())
    tube_n = len(cells)
    xbar_K = float((arb(M)/cells[-1]['xlo'].exp()).lower())
    if Kmax > xbar_K:                                  # extend the bulk below xbar with centre-manifold cells
        th_min = float((PS*arb(M)/arb(Kmax)).lower())*0.97
        cells = cells + centre_bulk_cells(th_min)
        log('[setup] bulk extended below xbar: %d centre-manifold cells (theta >= %.3g)' % (len(cells) - tube_n, th_min))
    son = sonic_piece(); cells, son = prep_bulk(cells, son)
    xmin = cells[-1]['xlo']
    Kcap = float((arb(M)/xmin.exp()).lower())
    if Kmax > Kcap:
        log('Kmax %.1f lowered to %.2f (deepest bulk cell)' % (Kmax, Kcap))
        Kmax = Kcap
    cm = centre_manifold_data(); cdat_C01, cdat_C12 = cm['C01'], cm['C12']
    vneed = PS*arb(M)/arb(Ka)
    cdat = cm['data'] + [bt_centre_data(cm['C01'], cm['C12'])] + tube_centre_data(cells, cm['C01'], cm['C12'], vneed)
    log('[setup] cells %d, sonic g=%s, A_max=%s, centre data %d cells (%.0fs)' % (len(cells), son['g'].str(4), son['Amax'].str(4), len(cdat), time.time() - T1))
    # K cells
    Ks = [Ka]
    while Ks[-1] < Kmax:
        Ks.append(min(Kmax, Ks[-1]*kratio))
    Kcells = list(zip(Ks[:-1], Ks[1:]))
    bulks = []
    for (a, b) in Kcells:
        r = bulk_ia(cells, son, M, arb(a), arb(b))
        if not r['ok']:
            log('BULK FAILS on K in [%.3f, %.3f]: collatz %s' % (a, b, r['collatz'].str(4))); return None
        bulks.append(r)
    log('[bulk] %d K-cells in [%.2f, %.2f] all contractive; worst collatz %s, max e_b %s, max nbar %s (%.0fs)' % (
        len(Kcells), Ka, Kmax, max(r['collatz'] for r in bulks).str(4), max(r['eb'] for r in bulks).str(4),
        max(r['nbar'] for r in bulks).str(4), time.time() - T1))
    # theta cells
    pi2 = arb.pi()/2
    ths = [pi2*k/nth for k in range(nth + 1)]
    th_idx = range(nth) if th_sub is None else th_sub
    fails = []; worst_margin = None
    for k in th_idx:
        cen = centre_theta(ths[k], ths[k + 1], Ka, M, cdat, ns=ns)
        if cen['w'] is None or cen['N'] is None:
            log('CENTRE supersolution fails on theta cell %d' % k); fails.append((k, 'super')); continue
        for (a, b), bk in zip(Kcells, bulks):
            ok, info = centre_test_ia(bk, cen, a, b, ths[k], ths[k + 1])
            if not ok:
                fails.append((k, a, b)); break
            m = float((info['lower'] - info['rhs']).lower())
            worst_margin = m if worst_margin is None else min(worst_margin, m)
        if log_each:
            log('  theta cell %d/%d [%.4f, %.4f]: %s (%.0fs)' % (k, nth, float(ths[k].mid()), float(ths[k + 1].mid()),
                'ok' if not fails or fails[-1][0] != k else 'FAIL K in [%.2f,%.2f]' % fails[-1][1:3] if len(fails[-1]) == 3 else 'FAIL', time.time() - T1))
    return dict(fails=fails, Kcells=Kcells, Kmax=Kmax, worst_margin=worst_margin, bulks=bulks)

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'run':
    # usage: run M Ka Kmax nth      -- the certificate
    M = float(sys.argv[2]); Ka = float(sys.argv[3]); Kmax = float(sys.argv[4]); nth = int(sys.argv[5])
    log('fw_cert run: M=%s Ka=%s Kmax=%s theta cells=%d, prec=%d bits, CAP=%d PTM=%d NSUB=%d NB1=%d' % (M, Ka, Kmax, nth, PREC, CAP, PTM, NSUB, NB1))
    res = certify(M, Ka, Kmax, nth, log_each=False)
    if res is None:
        log('VERDICT: bulk not certified'); sys.exit(1)
    log('theta cells: %d, K-cells: %d, failures: %s, worst margin (lower - e_b|w+|) %s' % (nth, len(res['Kcells']), res['fails'], res['worst_margin']))
    if not res['fails']:
        log('VERDICT (CERTIFIED, interval arithmetic on the E5 enclosure): no physical eigenvalue kappa with Re kappa >= 0 and %.4g <= |kappa| <= %.4g' % (Ka, res['Kmax']))
        log('RIGOROUS R0 (for the range above) = %.4g ; |kappa| > %.4g : bulk cells below xbar not built (OPEN, see PROOF-B4-v4 sec. 7)' % (Ka, res['Kmax']))
    else:
        log('VERDICT: NOT certified at Ka=%s (failing cells above)' % Ka)

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'bulkrange':
    # usage: bulkrange M K1 K2 ratio   -- bulk-only certificate on [K1, K2] (incl. centre-manifold bulk cells)
    M = float(sys.argv[2]); K1 = float(sys.argv[3]); K2 = float(sys.argv[4]); rt = float(sys.argv[5])
    cells = eval_cells(tube_cells()); tube_n = len(cells)
    th_min = float((PS*arb(M)/arb(K2)).lower())*0.97
    cells = cells + centre_bulk_cells(th_min)
    son = sonic_piece(); cells, son = prep_bulk(cells, son)
    log('bulk cells %d (centre-manifold %d), deepest x %s' % (len(cells), len(cells) - tube_n, cells[-1]['xlo'].str(5)))
    Ks = [K1]
    while Ks[-1] < K2: Ks.append(min(K2, Ks[-1]*rt))
    worst = None
    for a, b in zip(Ks[:-1], Ks[1:]):
        r = bulk_ia(cells, son, M, arb(a), arb(b))
        if not r['ok']:
            log('FAIL [%.2f, %.2f] collatz %s' % (a, b, r['collatz'].str(4))); break
        if worst is None or r['collatz'] > worst: worst = r['collatz']
        if a == K1 or b == K2 or len(Ks) < 12:
            log('K in [%.2f, %.2f]: collatz %s e_b %s nbar %s' % (a, b, r['collatz'].str(4), r['eb'].str(4), r['nbar'].str(4)))
    log('bulkrange done: %d K-cells, worst collatz %s' % (len(Ks) - 1, worst.str(4) if worst else None))

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'probe':
    # usage: probe M Ka Kmax nth step   (theta cells 0, step, 2 step, ... and the last)
    M = float(sys.argv[2]); Ka = float(sys.argv[3]); Kmax = float(sys.argv[4]); nth = int(sys.argv[5]); st = int(sys.argv[6])
    sub = sorted(set(list(range(0, nth, st)) + [nth - 1]))
    res = certify(M, Ka, Kmax, nth, th_sub=sub, log_each=True)
    if res: log('PROBE M=%s Ka=%s: %d theta cells tested, failures %s, worst margin %s' % (M, Ka, len(sub), res['fails'], res['worst_margin']))

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'centretest':
    M = float(sys.argv[2]); Ka = float(sys.argv[3])
    cells = eval_cells(tube_cells())
    cm = centre_manifold_data()
    cdat_C01, cdat_C12 = cm['C01'], cm['C12']
    log('C01 =', [[t.str(4) for t in row] for row in cm['C01']])
    log('C12 =', [[t.str(4) for t in row] for row in cm['C12']])
    vneed = PS*arb(M)/arb(Ka)
    cdat = cm['data'] + [bt_centre_data(cm['C01'], cm['C12'])] + tube_centre_data(cells, cm['C01'], cm['C12'], vneed)
    rs = running_sup(cdat, vneed)
    log('running sups at v_a=%s: max G0 %s G1 %s R0 %s R1 %s' % (vneed.str(4), max(max(r) for r in rs['G0']).str(3),
        max(max(r) for r in rs['G1']).str(3), max(max(r) for r in rs['R0']).str(3), max(max(r) for r in rs['R1']).str(3)))
    for w in ('G0', 'G1', 'R0', 'R1'):
        log('  %s:' % w, [[float(t.upper()) for t in row] for row in rs[w]])
    for th in (0.0, 1.0, math.pi/2 - 0.01):
        t0 = time.time()
        cen = centre_theta(arb(th), arb(th + 0.01), Ka, M, cdat, ns=int(sys.argv[4]) if len(sys.argv) > 4 else 120)
        for m in ('w', 'N'):
            if cen[m] is None: log('theta %.2f %s: supersolution FAILED' % (th, m)); continue
            Yf, B1, E2 = cen[m]
            log('theta %.2f %s: Yf=%s B1=%s E2=%s (%.0fs)' % (th, m, [y.str(4) for y in Yf], [b.str(4) for b in B1],
                [e.str(3) for e in E2], time.time() - t0))

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'bulktest':
    cells = eval_cells(tube_cells())
    son = sonic_piece()
    cells, son = prep_bulk(cells, son)
    log('cells %d, g=%s, Amax=%s' % (len(cells), son['g'].str(5), son['Amax'].str(5)))
    M = float(sys.argv[2])
    for K in sys.argv[3].split(','):
        Ka = arb(K); Kb = arb(K)
        r = bulk_ia(cells, son, M, Ka, Kb)
        log('K=%s M=%s: rho_float=%.3f collatz=%s ok=%s' % (K, M, r['rho_float'], r['collatz'].str(4), r['ok']))
        if r['ok']:
            log('   v=%s e_b=%s nbar=%s sig_c=%s chi_c=%s' % ([t.str(4) for t in r['v']], r['eb'].str(4), r['nbar'].str(4),
                r['sig_c'].str(4), r['chi_c'].str(4)))
            log('   A=%s' % [[a.str(3) for a in row] for row in r['A']])

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'selftest':
    cells = tube_cells()
    log('cells: %d (tau %d, x %d); x range [%s, %s]' % (len(cells), sum(c['kind'] == 'tau' for c in cells),
        sum(c['kind'] == 'x' for c in cells), cells[-1]['xlo'].str(6), cells[0]['xhi'].str(6)))
    # contiguity
    gaps = [(cells[i]['xlo'] - cells[i + 1]['xhi']) for i in range(len(cells) - 1)]
    log('max |gap| between consecutive cells:', max(g.abs_upper() for g in gaps).str(3))
    sel = eval_cells(cells[:3] + cells[64:66] + cells[300:302] + cells[-2:])
    sel = sel[:3] + sel[8:10] + sel[23:26] + sel[-4:]
    for c in sel:
        f = c['f']
        log('   box widths %s' % [b.rad().str(2) for b in c['box']])
        log('x in [%.5f, %.5f]: lp %s ell %s sig %s chi %s ellphi %s eAm %s q %s' % (float(c['xlo'].mid()), float(c['xhi'].mid()),
            f['lp'].str(5), f['ell'].str(5), f['sig'].str(5), f['chi'].str(5), f['ellphi'].str(5), f['eAm'].str(5), f['q'].str(5)))
    sp_ = sonic_piece()
    f = sp_['f']
    log('sonic piece s0=%s: lp %s ell %s dlm %s b11 %s sig %s chi %s dchi %s ellphi %s q %s' % (sp_['s0'].str(5), f['lp'].str(5),
        f['ell'].str(5), f['dlm'].str(5), f['b11'].str(5), f['sig'].str(5), f['chi'].str(5), f['dchi'].str(5), f['ellphi'].str(5), f['q'].str(5)))
