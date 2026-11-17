"""B5G pipeline: certified background re-integration + linear perturbations + characteristic function F(kappa).

Formulation (own): polar-areal slicing, central-lapse gauge alpha(t, r=0) = 1, self-similar (tau, x), independent
variable z = 1/N (N = alpha(-t)/(a r)); perturbation Y = (dmu, dOm, dV, eta=dN/N) ~ e^{kappa tau}, tau = -ln(-t).
Algorithm (own): Taylor models with interval coefficients, causal bootstrap for the tail (no contraction needed),
analytic (exponent-0) Frobenius solutions at the centre z=0 and at the sonic point z=zs.
F(kappa) = det[ Y_centre(zm), Y_sonic^(1..3)(zm) ];  F'(kappa) from the kappa-derivative system.
"""
import json
import sys
import time
import sympy as sp
import flint
from flint import arb, arb_mat, ctx
sys.path.insert(0, '<LAB>/b5g')
from tm import TM, R, codegen, horner

ctx.cap = 80

EX = json.load(open('<LAB>/b5g/exprs_g2.json'))
mu_, Om_, V_, N_, Omx_, Vx_, z_, sig_, U_, kap_ = sp.symbols('mu Om V N Om_x V_x z sig U kappa')
LOC = dict(mu=mu_, Om=Om_, V=V_, N=N_, Om_x=Omx_, V_x=Vx_, z=z_, sig=sig_, U=U_, kappa=kap_)


def S(name):
    return sp.sympify(EX[name], locals=LOC)


def SQRT(x):
    return x.sqrt()


NS = {'R': R, 'SQRT': SQRT}


def mk(exprs, args, name):
    src = codegen(exprs, args, name)
    exec(src, NS)
    raw = NS[name]

    def f(*xs):
        ref = next((x for x in xs if isinstance(x, TM)), None)
        out = raw(*xs)
        res = []
        for o in out:
            if ref is not None and not isinstance(o, TM):
                o = TM.const(o, ref.n, ref.T)
            elif ref is None and not isinstance(o, arb):
                o = arb(o)
            res.append(o)
        return res
    return f


f_MU = mk([S('MU')], (Om_, V_, N_), 'f_MU')
f_SL = mk([S('SL')], (mu_, Om_), 'f_SL')
f_GOV = mk([S('GO'), S('GV')], (mu_, Om_, V_, N_), 'f_GOV')
f_Sz = mk([S('Sz')], (V_, z_), 'f_Sz')
f_det = mk([S('detfac')], (mu_, Om_, V_), 'f_det')
FLN = ['Ax00', 'Ax01', 'Ax10', 'Ax11', 'At00', 'At01', 'At10', 'At11', 'Bu00', 'Bu01', 'Bu10', 'Bu11',
       'Bmu0', 'Bmu1', 'BNN0', 'BNN1']
f_fl = mk([S(n) for n in FLN], (mu_, Om_, V_, N_, Omx_, Vx_), 'f_fl')
f_grad = mk([S(n) for n in ('gH0', 'gH1', 'gH2', 'gS0', 'gS1', 'gS2')], (mu_, Om_, V_), 'f_grad')
f_Mo = mk([S('Mo'), S('H')], (mu_, Om_, V_, N_), 'f_Mo')
C = EX['centre']
f_cbg = mk([sp.sympify(C['dsig'], locals=LOC), sp.sympify(C['gU'], locals=LOC)], (z_, sig_, U_), 'f_cbg')
Mc = [[sp.sympify(C['M'][i][j], locals=LOC) for j in range(4)] for i in range(4)]
Mca, Mcb = [], []
for i in range(4):
    for j in range(4):
        e = Mc[i][j]
        assert sp.simplify(sp.diff(e, kap_, 2)) == 0
        Mca.append(e.subs(kap_, 0))
        Mcb.append(sp.diff(e, kap_))
f_cMa = mk(Mca, (z_, sig_, U_), 'f_cMa')
f_cMb = mk(Mcb, (z_, sig_, U_), 'f_cMb')
def SQ3f():
    return arb(3).sqrt()      # at the current working precision


def Vs_of(zs):
    SQ3 = SQ3f()
    return (1 - SQ3 * zs) / (SQ3 - zs)


def Oms_of(zs):
    SQ3 = SQ3f()
    return zs * (3 * zs - 2 * SQ3) / (2 * (-zs**2 + 4 * SQ3 * zs - 6))


def ball(lo, hi):
    lo, hi = arb(lo), arb(hi)
    return lo.union(hi)


def inside(a, b):
    """a strictly inside b (componentwise lists ok)"""
    return all(bb.contains_interior(aa) for aa, bb in zip(a, b))


def widen(x, rel, ab):
    out = []
    for xi in x:
        r = abs(xi).upper() * rel + ab
        out.append(xi + arb(0, r))
    return out


# --------------------------------------------------------------------------------------------- generic nonlinear solver
def nl_solve(Psi, y0, T, N, d, mode, u0=None, maxit=30, Lmat=None):
    """Unknown u = y' on t in T, y = y0 + t D, D = A[u].  Psi(D) -> list of d TMs (order >= len(D)-1).
    mode 'explicit' (u_n = Psi_n), 'affine' (Psi_n affine in D_n), 'sonic' (like affine, u0 supplied).
    Returns (yTM list, uTM list)."""
    u = [[None] * (N + 1) for _ in range(d)]

    def Dlow(n, extra=None):
        Ds = []
        for i in range(d):
            cs = [u[i][j] / (j + 1) for j in range(n)]
            cs.append(arb(0) if extra is None else extra[i])
            cs.append(arb(0))
            Ds.append(TM(cs, T))
        return Ds

    # L = d Psi_0 / d D_0 at the solution: the coefficient of D_n in Psi_n for every n (computed ONCE; computing it
    # by differences at order n would inherit the radius of c_n and blow up super-exponentially)
    L = Lmat
    if mode == 'affine':
        base0 = Psi(Dlow(0))
        L = arb_mat(d, d)
        for jj in range(d):
            e = [arb(0)] * d
            e[jj] = arb(1)
            col = Psi(Dlow(0, e))
            for i in range(d):
                L[i, jj] = col[i].c[0] - base0[i].c[0]
    for n in range(N + 1):
        if mode == 'sonic' and n == 0:
            for i in range(d):
                u[i][0] = u0[i]
            continue
        base = Psi(Dlow(n))
        c = [base[i].c[n] for i in range(d)]
        if mode == 'explicit':
            for i in range(d):
                u[i][n] = c[i]
            continue
        Mm = arb_mat(d, d)
        for i in range(d):
            for jj in range(d):
                Mm[i, jj] = (1 if i == jj else 0) - L[i, jj] / (n + 1)
        sol = Mm.solve(arb_mat([[ci] for ci in c]))
        for i in range(d):
            u[i][n] = sol[i, 0]
    uN = [u[i][N] for i in range(d)]
    W = widen(uN, 0.5, arb(10) ** (-ctx.dps + 5))
    for it in range(maxit):
        Ds = []
        for i in range(d):
            Ds.append(TM([u[i][j] / (j + 1) for j in range(N)] + [W[i] / (N + 1)], T))
        Uu = [p.lower(N) for p in Psi(Ds)]
        Rr = [p.c[N] for p in Uu]
        ok = inside(Rr, W) and inside(uN, W)
        if ok:
            break
        W = [Rr[i].union(uN[i]) for i in range(d)]
        W = widen(W, 0.25 * (1.5 ** it), arb(10) ** (-ctx.dps + 5))
    else:
        raise RuntimeError('tail bootstrap failed')
    uT = [TM([u[i][j] for j in range(N)] + [Rr[i]], T) for i in range(d)]
    yT = [(uT[i].average().mul_t()) + y0[i] for i in range(d)]
    return yT, uT


# --------------------------------------------------------------------------------------------- background
def psi_regular(z0, y0):
    def Psi(D):
        T = D[0].T
        t = TM([arb(0), arb(1)] + [arb(0)] * (D[0].n), T)
        Om = D[0].mul_t() + y0[0]
        V = D[1].mul_t() + y0[1]
        z = t + z0
        return bg_rhs(Om, V, z)
    return Psi


def bg_fields(Om, V, z):
    N = 1 / z
    mu = f_MU(Om, V, N)[0]
    zx = -z * f_SL(mu, Om)[0]
    GO, GV = f_GOV(mu, Om, V, N)
    Sz = f_Sz(V, z)[0]
    return N, mu, zx, GO, GV, Sz


def bg_rhs(Om, V, z):
    N, mu, zx, GO, GV, Sz = bg_fields(Om, V, z)
    z2 = z * z
    return [GO * z2 / Sz / zx, GV * z2 / Sz / zx]


def psi_centre(sig0):
    U0 = arb(-1) / 2

    def Psi(D):
        T = D[0].T
        t = TM([arb(0), arb(1)] + [arb(0)] * (D[0].n), T)
        sig = D[0].mul_t() + sig0
        U = D[1].mul_t() + U0
        ds, gU = f_cbg(t, sig, U)
        return [ds, gU.div_t_zero()]      # exact: gU(0, sig, -1/2) = 0 (centre_sym.py)
    return Psi


def psi_sonic(zs, b0, parts=False):
    def Psi(D):
        T = D[0].T
        t = TM([arb(0), arb(1)] + [arb(0)] * (D[0].n), T)
        Om = D[0].mul_t() + b0[0]
        V = D[1].mul_t() + b0[1]
        z = t + zs
        N, mu, zx, GO, GV, Sz = bg_fields(Om, V, z)
        z2 = z * z
        Szt = Sz.div_t_zero()             # exact: Sz(e(zs)) = 0
        GOt = GO.div_t_zero()             # exact: GO(e(zs)) = 0 (sonic_sym_g2.py)
        GVt = GV.div_t_zero()
        num = [GOt * z2 / zx, GVt * z2 / zx]
        if parts:
            return num, Szt
        return [num[0] / Szt, num[1] / Szt]
    return Psi


def sonic_slopes(zs, b0, approx=None):
    """roots of u*Shat0(u) = Ghat0(u); returns list of float approximations and a verifier."""
    T = arb(0)
    Pp = psi_sonic(zs, b0, parts=True)

    def ev(Dv):
        D = [TM([Dv[0], arb(0)], T), TM([Dv[1], arb(0)], T)]
        num, den = Pp(D)
        return [num[0].c[0], num[1].c[0]], den.c[0]
    n0, d0 = ev([arb(0), arb(0)])
    n1, d1 = ev([arb(1), arb(0)])
    n2, d2 = ev([arb(0), arb(1)])
    a = d0
    b = [d1 - d0, d2 - d0]
    c = n0
    G = [[n1[0] - n0[0], n2[0] - n0[0]], [n1[1] - n0[1], n2[1] - n0[1]]]
    return a, b, c, G


def krawczyk_slope(a, b, c, G, ubar, rad):
    """verify a unique root of F(u) = u (a + b.u) - c - G u in the box ubar +/- rad."""
    def Fv(u):
        s = a + b[0] * u[0] + b[1] * u[1]
        return [u[0] * s - c[0] - G[0][0] * u[0] - G[0][1] * u[1],
                u[1] * s - c[1] - G[1][0] * u[0] - G[1][1] * u[1]]

    def J(u):
        s = a + b[0] * u[0] + b[1] * u[1]
        return [[s + u[0] * b[0] - G[0][0], u[0] * b[1] - G[0][1]],
                [u[1] * b[0] - G[1][0], s + u[1] * b[1] - G[1][1]]]
    ub = [arb(x) for x in ubar]
    for _ in range(8):     # Newton polish on midpoints
        Jm = J(ub)
        Fm = Fv(ub)
        M = arb_mat([[Jm[0][0].mid(), Jm[0][1].mid()], [Jm[1][0].mid(), Jm[1][1].mid()]])
        dx = M.solve(arb_mat([[Fm[0].mid()], [Fm[1].mid()]]))
        ub = [arb((ub[0] - dx[0, 0]).mid()), arb((ub[1] - dx[1, 0]).mid())]
    X = [ub[0] + arb(0, rad), ub[1] + arb(0, rad)]
    Jm = J(ub)
    Y = arb_mat([[Jm[0][0].mid(), Jm[0][1].mid()], [Jm[1][0].mid(), Jm[1][1].mid()]]).inv()
    Y = arb_mat([[Y[0, 0].mid(), Y[0, 1].mid()], [Y[1, 0].mid(), Y[1, 1].mid()]])
    Fm = Fv(ub)
    JX = J(X)
    K = []
    for i in range(2):
        v = ub[i] - (Y[i, 0] * Fm[0] + Y[i, 1] * Fm[1])
        for j in range(2):
            Mij = (1 if i == j else 0) - (Y[i, 0] * JX[0][j] + Y[i, 1] * JX[1][j])
            v += Mij * (X[j] - ub[j])
        K.append(v)
    ok = inside(K, X)
    return ok, K


def build_background(cfg, log=print):
    ctx.prec = cfg['prec']
    sig0 = arb(cfg['sig0'])
    zs = arb(cfg['zs'])
    z1 = arb(cfg['z1'])
    zm = arb(cfg['zm'])
    Nb = cfg['Nbg']
    steps = []
    t0 = time.time()
    # centre step on [0, z1]
    T = ball(0, z1)
    yC, uC = nl_solve(psi_centre(sig0), [sig0, arb(-1) / 2], T, Nb, 2, 'affine')
    steps.append({'kind': 'centre', 'z0': arb(0), 'T': T, 'y': yC, 'u': uC, 'h': z1})
    sig1, U1 = yC[0].ev(z1), yC[1].ev(z1)
    y = [z1 * z1 * sig1, z1 * U1]
    log('centre step done  Om(z1)=%s V(z1)=%s  (%.1fs)' % (y[0].str(8, radius=True), y[1].str(8, radius=True), time.time() - t0))
    zc = z1
    zsm = arb(zs.mid())
    while True:
        dist = float((zsm - zc).mid())
        h = min(cfg['h0'], cfg['frac'] * dist, cfg['frac0'] * float(zc.mid()))
        hz = arb(float('%.5g' % h))     # exact dyadic step
        if (zc + hz).upper() >= zm.lower() or float((zm - zc - hz).mid()) < 0.3 * h:
            hz = zm - zc
            last = True
        else:
            last = False
        T = ball(0, hz)
        yT, uT = nl_solve(psi_regular(zc, y), y, T, Nb, 2, 'explicit')
        steps.append({'kind': 'reg', 'z0': zc, 'T': T, 'y': yT, 'u': uT, 'h': hz})
        if cfg.get('verbose'):
            log('  reg step z0=%s h=%s V(end)=%s' % (zc.str(8), hz.str(8), yT[1].ev(hz).str(10, radius=True)))
        y = [yT[0].ev(hz), yT[1].ev(hz)]
        zc = zc + hz
        if last:
            break
    log('regular steps: %d, at zm: Om=%s V=%s (%.1fs)' % (len(steps) - 1, y[0].str(12, radius=True), y[1].str(12, radius=True), time.time() - t0))
    # sonic step on [zm - zs, 0]
    b0 = [Oms_of(zs), Vs_of(zs)]
    a, b, c, G = sonic_slopes(zs, b0)
    # candidate roots via sigma-cubic (float)
    import mpmath as mp
    mp.mp.dps = 40
    A_ = mp.mpf(a.mid().str(40)); B_ = [mp.mpf(x.mid().str(40)) for x in b]; C_ = [mp.mpf(x.mid().str(40)) for x in c]
    G_ = [[mp.mpf(G[i][j].mid().str(40)) for j in range(2)] for i in range(2)]

    def ufromsig(s):
        Mm = mp.matrix([[s - G_[0][0], -G_[0][1]], [-G_[1][0], s - G_[1][1]]])
        return mp.lu_solve(Mm, mp.matrix(C_))

    def cub(s):
        u = ufromsig(s)
        return s - A_ - B_[0] * u[0] - B_[1] * u[1]
    roots = []
    for s0 in [mp.mpf(x) / 10 for x in range(-40, 41)]:
        try:
            r = mp.findroot(cub, s0)
            if all(abs(r - q) > 1e-20 for q in roots):
                roots.append(r)
        except Exception:
            pass
    cands = []
    for s in roots:
        u = ufromsig(s)
        cands.append((float(u[0]), float(u[1]), s))
    log('sonic slope candidates (dOm/dz, dV/dz, Shat0): %s' % [(round(x, 6), round(y_, 6), float(s)) for x, y_, s in cands])
    Ts = ball(zm - zs, arb(0))
    T_lo = ball((zm - zs).lower(), 0)
    chosen = None
    for (u0a, u1a, s) in cands:
        ok, K = krawczyk_slope(a, b, c, G, [u0a, u1a], 1e-25)
        if not ok:
            log('  candidate %s: Krawczyk failed' % ((u0a, u1a),))
            continue
        try:
            ySo, uSo = nl_solve(psi_sonic(zs, b0), b0, T_lo, Nb, 2, 'sonic', u0=K)
        except Exception as ex:
            log('  candidate %s: sonic TM failed (%s)' % ((u0a, u1a), ex))
            continue
        ym = [ySo[0].ev(zm - zs), ySo[1].ev(zm - zs)]
        ov = all((ym[i] - y[i]).contains(0) for i in range(2))
        log('  candidate slope %s: sonic-side Om(zm)=%s V(zm)=%s overlap=%s' % ((u0a, u1a), ym[0].str(12, radius=True), ym[1].str(12, radius=True), ov))
        if ov:
            chosen = (K, ySo, uSo, ym)
            break
    if chosen is None:
        raise RuntimeError('no sonic branch overlaps the centre integration')
    K, ySo, uSo, ym = chosen
    steps.append({'kind': 'sonic', 'z0': zs, 'T': T_lo, 'y': ySo, 'u': uSo, 'h': zm - zs, 'b0': b0, 'slope': K})
    info = {'steps': len(steps), 'Om_zm_centre': y[0].str(20, radius=True), 'V_zm_centre': y[1].str(20, radius=True),
            'Om_zm_sonic': ym[0].str(20, radius=True), 'V_zm_sonic': ym[1].str(20, radius=True),
            'slope': [K[0].str(20, radius=True), K[1].str(20, radius=True)], 'overlap': True, 'time_s': time.time() - t0}
    return steps, info


# --------------------------------------------------------------------------------------------- perturbation matrices
def tm_t(T, n):
    return TM([arb(0), arb(1)] + [arb(0)] * (n - 1), T)


def pert_tms(step):
    """Q0,Q1 (regular, dY/dz = (Q0 + kappa Q1) Y) or M0,M1 (centre/sonic: t dY/dz = M Y) as 4x4 lists of TMs."""
    kind = step['kind']
    T = step['T']
    y, u = step['y'], step['u']
    n = u[0].n
    t = tm_t(T, n + 1)
    if kind == 'centre':
        sig, U = y[0], y[1]
        Ma = f_cMa(t, sig, U)
        Mb = f_cMb(t, sig, U)
        return [[Ma[4 * i + j] for j in range(4)] for i in range(4)], [[Mb[4 * i + j] for j in range(4)] for i in range(4)], None
    Om, V = y[0], y[1]
    z = t + step['z0']
    N, mu, zx, GO, GV, Sz = bg_fields(Om, V, z)
    Omx = u[0] * zx
    Vx = u[1] * zx
    fl = dict(zip(FLN, f_fl(mu, Om, V, N, Omx, Vx)))
    g = f_grad(mu, Om, V)
    dfac = f_det(mu, Om, V)[0]
    z2 = z * z
    if kind == 'sonic':
        den = dfac * Sz.div_t_zero() / z2          # det Ax / t  (exact zero constant term of Sz)
    else:
        den = dfac * Sz / z2                        # det Ax
    adj = [[fl['Ax11'], -fl['Ax01']], [-fl['Ax10'], fl['Ax00']]]
    Cm = [[fl['Bmu0'], fl['Bu00'], fl['Bu01'], fl['BNN0']], [fl['Bmu1'], fl['Bu10'], fl['Bu11'], fl['BNN1']]]
    Ct = [[None, fl['At00'], fl['At01'], None], [None, fl['At10'], fl['At11'], None]]
    Z = TM.const(0, n, T)
    P0 = [[Z] * 4 for _ in range(4)]
    P1 = [[Z] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(4):
            P0[1 + i][j] = -(adj[i][0] * Cm[0][j] + adj[i][1] * Cm[1][j]) / den / zx
            if Ct[0][j] is not None:
                P1[1 + i][j] = -(adj[i][0] * Ct[0][j] + adj[i][1] * Ct[1][j]) / den / zx
    for j in range(3):
        if kind == 'sonic':
            P0[0][j] = g[j].mul_t() / zx
            P0[3][j] = g[3 + j].mul_t() / zx
        else:
            P0[0][j] = g[j] / zx
            P0[3][j] = g[3 + j] / zx
    extra = None
    if kind == 'sonic':
        # rank-one data at t = 0 for the kernel basis (exact structure: adj(Ax(0)) has proportional rows)
        A11, A01 = fl['Ax11'].c[0], fl['Ax01'].c[0]
        wa = [A11 * Cm[0][j].c[0] - A01 * Cm[1][j].c[0] for j in range(4)]
        wb = [arb(0), A11 * Ct[0][1].c[0] - A01 * Ct[1][1].c[0], A11 * Ct[0][2].c[0] - A01 * Ct[1][2].c[0], arb(0)]
        extra = {'wa': wa, 'wb': wb, 'A11': A11}
    # constraint (momentum) data for checks: mu_tau = Mo
    step['bgf'] = {'mu': mu, 'Om': Om, 'V': V, 'N': N, 'zx': zx}
    return P0, P1, extra


def lower_all(Mx, n):
    return [[e.lower(n) for e in row] for row in Mx]


def coef_mats(Mx, n):
    """list over j<=n of arb_mat with coefficient j"""
    out = []
    for j in range(n + 1):
        out.append(arb_mat([[Mx[i][k].c[j] for k in range(len(Mx[0]))] for i in range(len(Mx))]))
    return out


def matvec_tm(Mx, v):
    out = []
    for i in range(len(Mx)):
        s = None
        for k in range(len(v)):
            term = Mx[i][k] * v[k]
            s = term if s is None else s + term
        out.append(s)
    return out


def comb(A, B, kap):
    return [[A[i][j] + B[i][j] * kap for j in range(4)] for i in range(4)]


def lin_regular(Q, Q1, y0, yk0, T, N, h):
    """dY/dz = Q Y, dYk/dz = Q Yk + Q1 Y on t in T; returns values at t=h."""
    Q = lower_all(Q, N + 1)
    Q1 = lower_all(Q1, N + 1)
    Qc = coef_mats(Q, N)
    Q1c = coef_mats(Q1, N)
    Y = [arb_mat([[x] for x in y0])]
    Yk = [arb_mat([[x] for x in yk0])]
    u, uk = [], []
    for n in range(N + 1):
        s = Qc[0] * Y[n]
        sk = Qc[0] * Yk[n] + Q1c[0] * Y[n]
        for j in range(1, n + 1):
            s += Qc[j] * Y[n - j]
            sk += Qc[j] * Yk[n - j] + Q1c[j] * Y[n - j]
        u.append(s)
        uk.append(sk)
        Y.append(s / (n + 1))
        Yk.append(sk / (n + 1))
    uN = [u[N][i, 0] for i in range(4)] + [uk[N][i, 0] for i in range(4)]
    W = widen(uN, 0.5, arb(10) ** (-ctx.dps + 5))
    for it in range(30):
        D = [TM([u[j][i, 0] / (j + 1) for j in range(N)] + [W[i] / (N + 1)], T) for i in range(4)]
        Dk = [TM([uk[j][i, 0] / (j + 1) for j in range(N)] + [W[4 + i] / (N + 1)], T) for i in range(4)]
        Yt = [D[i].mul_t() + y0[i] for i in range(4)]
        Ykt = [Dk[i].mul_t() + yk0[i] for i in range(4)]
        ut = matvec_tm(Q, Yt)
        ukt = [a + b for a, b in zip(matvec_tm(Q, Ykt), matvec_tm(Q1, Yt))]
        Rr = [p.lower(N).c[N] for p in ut + ukt]
        if inside(Rr, W) and inside(uN, W):
            break
        W = widen([Rr[i].union(uN[i]) for i in range(8)], 0.25 * 1.5 ** it, arb(10) ** (-ctx.dps + 5))
    else:
        raise RuntimeError('linear tail bootstrap failed (regular)')
    yh, ykh = [], []
    for i in range(4):
        D = TM([u[j][i, 0] / (j + 1) for j in range(N)] + [Rr[i] / (N + 1)], T)
        Dk = TM([uk[j][i, 0] / (j + 1) for j in range(N)] + [Rr[4 + i] / (N + 1)], T)
        yh.append(y0[i] + h * D.ev(h))
        ykh.append(yk0[i] + h * Dk.ev(h))
    return yh, ykh


def lin_singular(M, M1, y0, yk0, T, N, tev):
    """t Y' = M Y with M(0) y0 = 0 (exact), t Yk' = M Yk + M1 Y with M(0)yk0 + M1(0)y0 = 0 (exact).
    Analytic (exponent-0) solution; returns values at t = tev.  Needs M of order >= N+2."""
    M = lower_all(M, N + 2)
    M1 = lower_all(M1, N + 2)
    Mh = [[e.drop0() for e in row] for row in M]       # (M - M(0))/t, exact identity, order N+1
    M1h = [[e.drop0() for e in row] for row in M1]
    Mc = coef_mats(M, N)
    M1c = coef_mats(M1, N)
    Mhc = coef_mats(Mh, N)
    M1hc = coef_mats(M1h, N)
    y0m = arb_mat([[x] for x in y0])
    yk0m = arb_mat([[x] for x in yk0])
    u, uk = [], []
    I4 = arb_mat(4, 4)
    for i in range(4):
        I4[i, i] = 1
    for n in range(N + 1):
        rhs = Mhc[n] * y0m
        rhsk = Mhc[n] * yk0m + M1hc[n] * y0m
        for j in range(n):
            rhs += Mc[n - j] * u[j] / (j + 1)
            rhsk += Mc[n - j] * uk[j] / (j + 1) + M1c[n - j] * u[j] / (j + 1)
        Ls = I4 - Mc[0] / (n + 1)
        un = Ls.solve(rhs)
        rhsk += M1c[0] * un / (n + 1)
        ukn = Ls.solve(rhsk)
        u.append(un)
        uk.append(ukn)
    uN = [u[N][i, 0] for i in range(4)] + [uk[N][i, 0] for i in range(4)]
    W = widen(uN, 0.5, arb(10) ** (-ctx.dps + 5))
    for it in range(40):
        D = [TM([u[j][i, 0] / (j + 1) for j in range(N)] + [W[i] / (N + 1)], T) for i in range(4)]
        Dk = [TM([uk[j][i, 0] / (j + 1) for j in range(N)] + [W[4 + i] / (N + 1)], T) for i in range(4)]
        y0t = [TM.const(x, N, T) for x in y0]
        yk0t = [TM.const(x, N, T) for x in yk0]
        ut = [a + b for a, b in zip(matvec_tm(Mh, y0t), matvec_tm(M, D))]
        ukt = [a + b + c + d for a, b, c, d in zip(matvec_tm(Mh, yk0t), matvec_tm(M1h, y0t), matvec_tm(M, Dk), matvec_tm(M1, D))]
        Rr = [p.lower(N).c[N] for p in ut + ukt]
        if inside(Rr, W) and inside(uN, W):
            break
        W = widen([Rr[i].union(uN[i]) for i in range(8)], 0.25 * 1.5 ** it, arb(10) ** (-ctx.dps + 5))
    else:
        raise RuntimeError('linear tail bootstrap failed (singular)')
    yv, ykv = [], []
    for i in range(4):
        D = TM([u[j][i, 0] / (j + 1) for j in range(N)] + [Rr[i] / (N + 1)], T)
        Dk = TM([uk[j][i, 0] / (j + 1) for j in range(N)] + [Rr[4 + i] / (N + 1)], T)
        yv.append(y0[i] + tev * D.ev(tev))
        ykv.append(yk0[i] + tev * Dk.ev(tev))
    return yv, ykv, [u[0][i, 0] for i in range(4)]


class Pert:
    def __init__(self, steps, cfg):
        self.steps = steps
        self.cfg = cfg
        self.mats = []
        t0 = time.time()
        for st in steps:
            r = pert_tms(st)
            self.mats.append(r)
        self.tprep = time.time() - t0

    def sonic_data(self, kap):
        P0, P1, ex = self.mats[-1]
        w = [ex['wa'][j] + kap * ex['wb'][j] for j in range(4)]
        M00 = [[(P0[i][j] + P1[i][j] * kap).c[0] for j in range(4)] for i in range(4)]
        lam = M00[1][1] + M00[2][2]
        return w, ex['wb'], lam

    def F(self, kap, want_vectors=False):
        cfg = self.cfg
        Np = cfg['Np']
        # centre
        st = self.steps[0]
        Ma, Mb, _ = self.mats[0]
        sig0 = st['y'][0].c[0]
        M = comb(Ma, Mb, kap)
        y0 = [arb(2) / 3, arb(1), -kap / (4 * sig0), arb(0)]
        yk0 = [arb(0), arb(0), -1 / (4 * sig0), arb(0)]
        z1 = st['h']
        ys, yks, _ = lin_singular(M, Mb, y0, yk0, st['T'], Np, z1)
        sc = [z1 * z1, z1 * z1, z1, z1 * z1]
        y = [ys[i] * sc[i] for i in range(4)]
        yk = [yks[i] * sc[i] for i in range(4)]
        for idx in range(1, len(self.steps) - 1):
            st = self.steps[idx]
            P0, P1, _ = self.mats[idx]
            Q = comb(P0, P1, kap)
            y, yk = lin_regular(Q, P1, y, yk, st['T'], Np, st['h'])
        yc, ykc = y, yk
        # sonic basis
        st = self.steps[-1]
        P0, P1, ex = self.mats[-1]
        M = comb(P0, P1, kap)
        w, wb, lam = self.sonic_data(kap)
        piv = 1 if abs(w[1]).lower() > abs(w[2]).lower() else 2
        if w[piv].contains(0):
            raise RuntimeError('sonic pivot contains 0')
        cols, colsk = [], []
        for a in range(4):
            if a == piv:
                continue
            y0 = [arb(0)] * 4
            yk0 = [arb(0)] * 4
            y0[a] = w[piv]
            y0[piv] = -w[a]
            yk0[a] = wb[piv]
            yk0[piv] = -wb[a]
            yv, ykv, _ = lin_singular(M, P1, y0, yk0, st['T'], Np, st['h'])
            cols.append(yv)
            colsk.append(ykv)
        allc = [yc] + cols
        allk = [ykc] + colsk
        Mt = arb_mat([[allc[j][i] for j in range(4)] for i in range(4)])
        Fv = Mt.det()
        dF = arb(0)
        for jj in range(4):
            Mj = arb_mat([[(allk[j][i] if j == jj else allc[j][i]) for j in range(4)] for i in range(4)])
            dF += Mj.det()
        if want_vectors:
            return Fv, dF, lam, allc
        return Fv, dF, lam
