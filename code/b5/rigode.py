"""Validated Chebyshev-Picard integration (sup-norm contraction) for rational ODE systems.

Component i on a step [a,b] (x = mid + h t):
  op 'int':  u_i(x) = base_i + int_a^x G_i(s,u(s)) ds
  op 'K3' :  u_i(x) = base_i + x^{-3} int_0^x s^2 G_i(s,u(s)) ds      (only on [0,b]; centre)
G_i = num_i/den_i (polynomials, evaluated on Chebyshev models).  Proof: for the numerical guess ubar,
  Y_i >= sup|Phi(ubar)_i - ubar_i|  (rigorous: exact polynomial arithmetic + division residual bound),
  L_ij >= sup over the ball of |dG_i/du_j| (interval evaluation on sub-pieces),
  if Y_i + opnorm_i * sum_j L_ij r_j < r_i for all i then Phi maps the sup-ball B(ubar, r) into itself and is a
  contraction in the r-weighted max norm, so the unique solution lies in ubar +- r on the whole step.
opnorm: 'int' -> (b-a), 'K3' -> 1/3.
"""
from flint import arb, arb_poly, ctx
from cheb import CM, cheb_nodes, interp, aup, hull
import sympy as sp


def lamb(args, expr):
    # guard: no non-integer rationals may be baked in (they would become python floats)
    for at in sp.sympify(expr).atoms(sp.Rational):
        assert at.q == 1, ('non-integer rational in expression', at)
    return sp.lambdify(args, expr, modules=[{}])


class System:
    """symbolic rational system with parameter k; ops, bases supplied at run time"""

    def __init__(self, zsym, usyms, ksym, G, ops):
        self.zs, self.us, self.ks = zsym, usyms, ksym
        self.n = len(usyms)
        self.ops = ops
        args = (ksym, zsym) + tuple(usyms)
        self.num, self.den, self.full, self.jac = [], [], [], []
        for g in G:
            nn, dd = sp.fraction(sp.together(g))
            self.num.append(lamb(args, sp.expand(nn)))
            self.den.append(lamb(args, sp.expand(dd)))
            self.full.append(lamb(args, g))
            row = []
            for u in usyms:
                row.append(lamb(args, sp.diff(g, u)))
            self.jac.append(row)


def K3_poly(cm, z1):
    """apply K3 f = x^-3 int_0^x s^2 f(s) ds to the polynomial part of a model on [0,z1] (t <-> x = z1(1+t)/2)."""
    c = cm.c
    s = arb_poly([arb(-1), 2 / z1])          # t as polynomial in x
    T0, T1 = arb_poly([arb(1)]), s
    p = arb_poly([c[0]])
    if len(c) > 1:
        p += c[1] * T1
    for kk in range(2, len(c)):
        T0, T1 = T1, 2 * s * T1 - T0
        p += c[kk] * T1
    m = p.coeffs()
    m = [m[j] / (j + 3) for j in range(len(m))]
    # back to t: x = z1 (1+t)/2 ; Horner in Chebyshev basis
    res = [arb(0)]

    def times_t(v):
        out = [arb(0)] * (len(v) + 1)
        for j, a in enumerate(v):
            if j == 0:
                out[1] += a
            else:
                out[j + 1] += a / 2
                out[j - 1] += a / 2
        return out
    half = z1 / 2
    for a in reversed(m):
        # res = x*res + a = half*(1+t)*res + a
        tr = times_t(res)
        res = [half * ((res[j] if j < len(res) else 0) + tr[j]) for j in range(len(tr))]
        res[0] += a
    return CM(res, arb(0), cm.M).trunc()


def eval_models(fn, k, zcm, ucms):
    return fn(k, zcm, *ucms)


def guess(sys_, k, a, b, bases, M, u_init=None, iters=80, tol=None):
    """non-rigorous Chebyshev-Picard iteration (Clenshaw-Curtis); returns list of coefficient lists (exact arbs)"""
    h = (b - a) / 2
    mid = (a + b) / 2
    nodes = cheb_nodes(M)
    xs = [mid + h * t for t in nodes]
    n = sys_.n
    if u_init is None:
        cur = [CM([arb(bases[i].mid())] + [arb(0)] * M, arb(0), M) for i in range(n)]
    else:
        cur = u_init
    for it in range(iters):
        vals = [[cur[i].eval_t(t) for t in nodes] for i in range(n)]
        new = []
        for i in range(n):
            gv = [arb(sys_.full[i](k, xs[j], *[vals[q][j] for q in range(n)])) for j in range(M + 1)]
            gv = [arb(v.mid()) for v in gv]
            qc = CM(interp(gv, M), arb(0), M)
            if sys_.ops[i] == 'int':
                nm = qc.integral_from_left(h) + arb(bases[i].mid())
            else:
                nm = K3_poly(qc, b) + arb(bases[i].mid())
            nm = CM([arb(x.mid()) for x in nm.c], arb(0), M)
            new.append(nm)
        diff = max(abs(float((new[i] - cur[i]).l1().mid())) for i in range(n))
        cur = new
        if tol is not None and diff < tol:
            break
    return cur, diff


def validate(sys_, k, a, b, bases, ubar, Mx=None, pieces=16, rguess=None, verbose=False, kappa_max=None):
    """rigorous step.  bases: list of arb (may be balls).  returns dict with r, Y, L, endpoint enclosures."""
    n = sys_.n
    M = ubar[0].M
    Mx = Mx or 2 * M + 8
    h = (b - a) / 2
    mid = (a + b) / 2
    zcm = CM([mid, h] + [arb(0)] * (Mx - 1), arb(0), Mx)
    ub = [CM(u.c + [arb(0)] * (Mx - M), arb(0), Mx) for u in ubar]
    nodes = cheb_nodes(Mx)
    xs = [mid + h * t for t in nodes]
    vals = [[ubar[i].eval_t(t) for t in nodes] for i in range(n)]
    Y = []
    for i in range(n):
        num = eval_models(sys_.num[i], k, zcm, ub)
        den = eval_models(sys_.den[i], k, zcm, ub)
        gv = [arb(arb(sys_.full[i](k, xs[j], *[vals[q][j] for q in range(n)])).mid()) for j in range(Mx + 1)]
        q = CM(interp(gv, Mx), arb(0), Mx)
        R = num - den * q
        drange = den.range(pieces)
        if drange.contains(0):
            raise ValueError('denominator range contains 0 on step [%s,%s] comp %d' % (a, b, i))
        dmin = drange.abs_lower()
        gerr = aup(R.supbound() / dmin)
        qq = CM(q.c, gerr, Mx)
        if sys_.ops[i] == 'int':
            phi = qq.integral_from_left(h) + bases[i]
        else:
            ph = K3_poly(CM(q.c, arb(0), Mx), b)
            phi = CM(ph.c, aup(ph.err + gerr / 3), Mx) + bases[i]
        Y.append(aup((phi - ub[i]).supbound()))
    opn = [aup(b - a) if sys_.ops[i] == 'int' else arb(1) / 3 for i in range(n)]

    def lip(r):
        # sup |dG_i/du_j| over pieces: x-range of piece, u-range of piece +- r
        L = [[arb(0)] * n for _ in range(n)]
        pi = arb.pi()
        for p in range(pieces):
            th = hull(pi * p / pieces, pi * (p + 1) / pieces)
            tt = th.cos()
            xr = mid + h * tt
            ur = []
            for i in range(n):
                v = ubar[i].c[0]
                for kk in range(1, len(ubar[i].c)):
                    v += ubar[i].c[kk] * (kk * th).cos()
                ur.append(v + arb(0, r[i]))
            for i in range(n):
                for j in range(n):
                    d = abs(arb(sys_.jac[i][j](k, xr, *ur)))
                    L[i][j] = hull(L[i][j], arb(aup(d))) if p else arb(aup(d))
        return [[aup(L[i][j]) for j in range(n)] for i in range(n)]

    # choose r: smallest self-consistent radius (start from Y, inflate slowly)
    tiny = arb(10) ** (-(ctx.prec * 3 // 10))
    r = rguess or [aup(Y[i] * arb('1.0001') + tiny) for i in range(n)]
    ok = False
    for attempt in range(40):
        L = lip(r)
        lhs = [aup(Y[i] + opn[i] * sum((L[i][j] * r[j] for j in range(n)), arb(0))) for i in range(n)]
        if all(lhs[i] < r[i] for i in range(n)):
            ok = True
            break
        r = [aup(lhs[i] * arb('1.05') + tiny) for i in range(n)]
    if not ok:
        raise ValueError('contraction failed on [%s,%s]: Y=%s' % (a, b, [float(y.mid()) for y in Y]))
    kap = max(float((opn[i] * sum((L[i][j] * r[j] for j in range(n)), arb(0)) / r[i]).mid()) for i in range(n))
    if kappa_max is not None and kap > kappa_max:
        raise ValueError('contraction factor %.3f above %.3f' % (kap, kappa_max))
    endR = [ubar[i].right() + arb(0, r[i]) for i in range(n)]
    return {'r': r, 'Y': Y, 'L': L, 'end': endR, 'ubar': ubar, 'a': a, 'b': b, 'kappa': kap}


def step_enclosure(res, i, pieces=16):
    """enclosure of component i over the whole step"""
    return res['ubar'][i].range(pieces) + arb(0, res['r'][i])


def integrate(sys_, k, a, bfinal, u0, M=24, h0=None, hmin=None, log=None, record=None, maxstep=None, gtol=1e-40, ctol=1e-45, kappa_max=0.5):
    """adaptive validated integration from a to bfinal (a < bfinal).  record: list to append step results."""
    a = arb(a)
    bfinal = arb(bfinal)
    h = arb(h0) if h0 is not None else (bfinal - a) / 8
    ub = None
    while True:
        b = a + h
        last = False
        if b.mid() >= bfinal.mid() or abs(float((b - bfinal).mid())) < 1e-40:
            b = bfinal
            last = True
        try:
            ub, dif = guess(sys_, k, a, b, u0, M, iters=120, tol=gtol * 1e-10)
            if dif > gtol:
                raise ValueError('guess not converged %g' % dif)
            tail = max(float(abs(u.c[-1]).mid()) + float(abs(u.c[-2]).mid()) for u in ub)
            rad0 = max(float(x.rad()) for x in u0)
            if tail > max(ctol, 1e-4 * rad0):
                raise ValueError('coefficients not decayed %g' % tail)
            res = validate(sys_, k, a, b, u0, ub, kappa_max=kappa_max)
        except ValueError as ex:
            h = h / 2
            if log:
                log('  halve at a=%.6g: %s' % (float(a.mid()), str(ex)[:80]))
            if hmin is not None and float(h.mid()) < hmin:
                raise
            continue
        if record is not None:
            record.append(res)
        u0 = res['end']
        if log:
            log('  step [%.6g, %.6g] r=%s' % (float(a.mid()), float(b.mid()), ['%.2e' % float(x.mid()) for x in res['r']]))
        if last:
            return u0
        a = b
        h = h * 1.5
        if maxstep is not None and float(h.mid()) > maxstep:
            h = arb(maxstep)
