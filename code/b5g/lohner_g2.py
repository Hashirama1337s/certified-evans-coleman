"""Wrapping-free propagation (Lohner parallelepiped + QR) for the regular region, on top of pipeline_g2.

Background: mean-value form  y(h; y0) in y(h; c) + J_box (y0 - c), J from the variational equation integrated
with the box initial data; set kept as c + B r (B orthonormal from QR of mid(J B)).
Perturbation (linear): per step fundamental matrix Phi8 = [[Phi,0],[Phi',Phi]] (kappa-derivative block), set
propagated the same way.
"""
import json
import sys
import time
import numpy as np
import sympy as sp
from flint import arb, arb_mat, ctx
sys.path.insert(0, '<LAB>/b5g')
import pipeline_g2 as PG
from pipeline_g2 import TM, nl_solve, ball, inside, widen, lower_all, coef_mats, matvec_tm, comb, lin_singular

ctx.cap = 80
BJ = json.load(open('<LAB>/b5g/bgjac_g2.json'))
_Om, _V, _z = sp.symbols('Om V z')
_loc = dict(Om=_Om, V=_V, z=_z)
f_bgJ = PG.mk([sp.sympify(BJ['f0'], locals=_loc), sp.sympify(BJ['f1'], locals=_loc)] +
              [sp.sympify(BJ['J'][i][j], locals=_loc) for i in range(2) for j in range(2)], (_Om, _V, _z), 'f_bgJ')


def psi_aug(z0, y0):
    """unknowns (Om, V, Phi00, Phi01, Phi10, Phi11); Phi' = Jf Phi, Phi(0) = I"""
    def Psi(D):
        T = D[0].T
        t = TM([arb(0), arb(1)] + [arb(0)] * (D[0].n), T)
        y = [D[i].mul_t() + y0[i] for i in range(6)]
        z = t + z0
        f0, f1, J00, J01, J10, J11 = f_bgJ(y[0], y[1], z)
        P = [[y[2], y[3]], [y[4], y[5]]]
        Jf = [[J00, J01], [J10, J11]]
        out = [f0, f1]
        for i in range(2):
            for j in range(2):
                out.append(Jf[i][0] * P[0][j] + Jf[i][1] * P[1][j])
        return out
    return Psi


def psi_pt(z0, y0):
    def Psi(D):
        T = D[0].T
        t = TM([arb(0), arb(1)] + [arb(0)] * (D[0].n), T)
        y = [D[i].mul_t() + y0[i] for i in range(2)]
        z = t + z0
        r = f_bgJ(y[0], y[1], z)
        return [r[0], r[1]]
    return Psi


def mid_arb(x):
    return arb(x.mid())


def qr_basis(A):
    """orthonormal point basis from QR of mid(A) (float), returned as exact arb_mat."""
    n = A.nrows()
    M = np.array([[float(A[i, j].mid()) for j in range(n)] for i in range(n)])
    Q, Rr = np.linalg.qr(M)
    return arb_mat([[arb(float(Q[i, j])) for j in range(n)] for i in range(n)])


def lohner_update(phi_c, A, r):
    """phi_c: interval image of the centre point (list); A = Phi*B (interval); r: interval list.
    returns new (c, B, r)."""
    n = len(phi_c)
    c = [mid_arb(x) for x in phi_c]
    B = qr_basis(A)
    Binv = B.inv()
    d = arb_mat([[phi_c[i] - c[i]] for i in range(n)])
    rv = arb_mat([[x] for x in r])
    rn = Binv * d + (Binv * A) * rv
    return c, B, [rn[i, 0] for i in range(n)]


def set_eval(c, B, r):
    n = len(c)
    v = B * arb_mat([[x] for x in r])
    return [c[i] + v[i, 0] for i in range(n)]


def build_background2(cfg, log=print):
    ctx.prec = cfg['prec']
    sig0 = arb(cfg['sig0'])
    zs = arb(cfg['zs'])
    z1 = arb(cfg['z1'])
    zm = arb(cfg['zm'])
    Nb = cfg['Nbg']
    steps = []
    t0 = time.time()
    T = ball(0, z1)
    yC, uC = nl_solve(PG.psi_centre(sig0), [sig0, arb(-1) / 2], T, Nb, 2, 'affine')
    steps.append({'kind': 'centre', 'z0': arb(0), 'T': T, 'y': yC, 'u': uC, 'h': z1})
    sig1, U1 = yC[0].ev(z1), yC[1].ev(z1)
    y = [z1 * z1 * sig1, z1 * U1]
    log('centre step: Om(z1)=%s V(z1)=%s' % (y[0].str(25, radius=True), y[1].str(25, radius=True)))
    c = [mid_arb(x) for x in y]
    B = arb_mat([[1, 0], [0, 1]])
    r = [y[i] - c[i] for i in range(2)]
    zc = z1
    zsm = arb(zs.mid())
    while True:
        dist = float((zsm - zc).mid())
        h = min(cfg['h0'], cfg['frac'] * dist, cfg['frac0'] * float(zc.mid()))
        hz = arb(float('%.5g' % h))
        if float((zm - zc - hz).mid()) < 0.3 * h:
            hz = zm - zc
            last = True
        else:
            last = False
        T = ball(0, hz)
        box = set_eval(c, B, r)
        yA, uA = nl_solve(psi_aug(zc, box + [arb(1), arb(0), arb(0), arb(1)]), box + [arb(1), arb(0), arb(0), arb(1)],
                          T, Nb, 6, 'explicit')
        yP, uP = nl_solve(psi_pt(zc, c), c, T, Nb, 2, 'explicit')
        phi_c = [yP[i].ev(hz) for i in range(2)]
        J = arb_mat([[yA[2].ev(hz), yA[3].ev(hz)], [yA[4].ev(hz), yA[5].ev(hz)]])
        c, B, r = lohner_update(phi_c, J * B, r)
        steps.append({'kind': 'reg', 'z0': zc, 'T': T, 'y': yA[:2], 'u': uA[:2], 'h': hz})
        zc = zc + hz
        if cfg.get('verbose'):
            yy = set_eval(c, B, r)
            log('  reg z=%s  Om=%s V=%s' % (zc.str(8), yy[0].str(25, radius=True), yy[1].str(25, radius=True)))
        if last:
            break
    y = set_eval(c, B, r)
    log('regular steps: %d, at zm: Om=%s V=%s (%.1fs)' % (len(steps) - 1, y[0].str(30, radius=True), y[1].str(30, radius=True), time.time() - t0))
    son = sonic_background(cfg, zs, zm, y, Nb, log)
    steps.append(son['step'])
    info = {'steps': len(steps), 'Om_zm_centre': y[0].str(30, radius=True), 'V_zm_centre': y[1].str(30, radius=True),
            'Om_zm_sonic': son['ym'][0].str(30, radius=True), 'V_zm_sonic': son['ym'][1].str(30, radius=True),
            'slope': [son['K'][0].str(25, radius=True), son['K'][1].str(25, radius=True)], 'overlap': True,
            'sonic_slope_candidates': son['cands'], 'time_s': round(time.time() - t0, 1)}
    return steps, info


def sonic_background(cfg, zs, zm, y, Nb, log):
    b0 = [PG.Oms_of(zs), PG.Vs_of(zs)]
    a, b, c, G = PG.sonic_slopes(zs, b0)
    import mpmath as mp
    mp.mp.dps = 50
    A_ = mp.mpf(a.mid().str(50, radius=False)); B_ = [mp.mpf(x.mid().str(50, radius=False)) for x in b]; C_ = [mp.mpf(x.mid().str(50, radius=False)) for x in c]
    G_ = [[mp.mpf(G[i][j].mid().str(50, radius=False)) for j in range(2)] for i in range(2)]
    # u (A + B.u) = C + G u : eliminate via s = A + B.u ; u = (sI - G)^{-1} C ; cubic in s (polynomial form)
    s = sp.symbols('s')
    Gm = sp.Matrix([[sp.Float(str(G_[0][0]), 50), sp.Float(str(G_[0][1]), 50)], [sp.Float(str(G_[1][0]), 50), sp.Float(str(G_[1][1]), 50)]])
    Cv = sp.Matrix([sp.Float(str(C_[0]), 50), sp.Float(str(C_[1]), 50)])
    Mm = s * sp.eye(2) - Gm
    poly = sp.expand((s - sp.Float(str(A_), 50)) * Mm.det() - (sp.Matrix([[sp.Float(str(B_[0]), 50), sp.Float(str(B_[1]), 50)]]) * Mm.adjugate() * Cv)[0])
    rts = [complex(x) for x in np.roots([float(cc) for cc in sp.Poly(poly, s).all_coeffs()])]
    cands = []
    for sr in rts:
        if abs(sr.imag) > 1e-6 * (1 + abs(sr)):
            continue
        uu = Mm.subs(s, sr.real).inv() * Cv
        cands.append((float(uu[0]), float(uu[1])))
    log('sonic slope candidates (dOm/dz, dV/dz): %s' % [(round(x, 8), round(y_, 8)) for x, y_ in cands])
    T_lo = ball((zm - zs).lower(), 0)
    if cfg.get('sonic_exact_bg') is not None:
        ySo, uSo = cfg['sonic_exact_bg'](T_lo, zs, Nb + 1)
        ym = [ySo[0].ev(zm - zs), ySo[1].ev(zm - zs)]
        ov = all((ym[i] - y[i]).contains(0) for i in range(2))
        log('  exact sonic background: Om(zm)=%s V(zm)=%s overlap=%s' % (ym[0].str(30, radius=True), ym[1].str(30, radius=True), ov))
        if not ov:
            raise RuntimeError('exact sonic background does not overlap')
        K = [uSo[0].c[0], uSo[1].c[0]]
        step = {'kind': 'sonic', 'z0': zs, 'T': T_lo, 'y': ySo, 'u': uSo, 'h': zm - zs, 'b0': b0, 'slope': K}
        return {'step': step, 'ym': ym, 'K': K, 'cands': cands}
    if cfg.get('sonic_slope_exact') is not None:
        cands = [('exact', None)]
    for (u0a, u1a) in cands:
        if u0a == 'exact':
            ok, K = True, cfg['sonic_slope_exact']()
        else:
            ok, K = PG.krawczyk_slope(a, b, c, G, [u0a, u1a], 1e-28)
        if not ok:
            log('  candidate %s: Krawczyk failed' % ((u0a, u1a),))
            continue
        try:
            Sh0 = a + b[0] * K[0] + b[1] * K[1]
            Lm = arb_mat([[(G[i][j] - K[i] * b[j]) / Sh0 for j in range(2)] for i in range(2)])
            ySo, uSo = nl_solve(PG.psi_sonic(zs, b0), b0, T_lo, Nb, 2, 'sonic', u0=K, Lmat=Lm)
        except Exception as ex:
            log('  candidate %s: sonic TM failed (%s)' % ((u0a, u1a), ex))
            continue
        ym = [ySo[0].ev(zm - zs), ySo[1].ev(zm - zs)]
        ov = all((ym[i] - y[i]).contains(0) for i in range(2))
        log('  candidate slope %s: sonic-side Om(zm)=%s V(zm)=%s overlap=%s' % ((u0a, u1a), ym[0].str(30, radius=True), ym[1].str(30, radius=True), ov))
        if ov:
            step = {'kind': 'sonic', 'z0': zs, 'T': T_lo, 'y': ySo, 'u': uSo, 'h': zm - zs, 'b0': b0, 'slope': K}
            return {'step': step, 'ym': ym, 'K': K, 'cands': cands}
    raise RuntimeError('no sonic branch overlaps the centre integration')


# ------------------------------------------------------------------------------------------------- perturbation
def lin_regular_mat(Q, Q1, T, N, h):
    """fundamental matrices: Phi' = Q Phi, Phi(0)=I ; Psi' = Q Psi + Q1 Phi, Psi(0)=0. Returns (Phi(h), Psi(h))."""
    Q = lower_all(Q, N + 1)
    Q1 = lower_all(Q1, N + 1)
    Qc = coef_mats(Q, N)
    Q1c = coef_mats(Q1, N)
    I4 = arb_mat(4, 4)
    for i in range(4):
        I4[i, i] = 1
    Y = [I4]
    Yk = [arb_mat(4, 4)]
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
    # tail bootstrap on all 32 entries
    uN = [u[N][i, j] for j in range(4) for i in range(4)] + [uk[N][i, j] for j in range(4) for i in range(4)]
    W = widen(uN, 0.5, arb(10) ** (-ctx.dps + 5))
    for it in range(30):
        Rr = []
        Rk = []
        for j in range(4):
            D = [TM([u[m][i, j] / (m + 1) for m in range(N)] + [W[4 * j + i] / (N + 1)], T) for i in range(4)]
            Dk = [TM([uk[m][i, j] / (m + 1) for m in range(N)] + [W[16 + 4 * j + i] / (N + 1)], T) for i in range(4)]
            Yt = [D[i].mul_t() + (1 if i == j else 0) for i in range(4)]
            Ykt = [Dk[i].mul_t() for i in range(4)]
            ut = matvec_tm(Q, Yt)
            ukt = [a + b for a, b in zip(matvec_tm(Q, Ykt), matvec_tm(Q1, Yt))]
            Rr += [p.lower(N).c[N] for p in ut]
            Rk += [p.lower(N).c[N] for p in ukt]
        Rall = Rr + Rk
        if inside(Rall, W) and inside(uN, W):
            break
        W = widen([Rall[i].union(uN[i]) for i in range(32)], 0.25 * 1.5 ** it, arb(10) ** (-ctx.dps + 5))
    else:
        raise RuntimeError('matrix tail bootstrap failed')
    Phi = arb_mat(4, 4)
    Psi = arb_mat(4, 4)
    for j in range(4):
        for i in range(4):
            D = TM([u[m][i, j] / (m + 1) for m in range(N)] + [Rall[4 * j + i] / (N + 1)], T)
            Dk = TM([uk[m][i, j] / (m + 1) for m in range(N)] + [Rall[16 + 4 * j + i] / (N + 1)], T)
            Phi[i, j] = (1 if i == j else 0) + h * D.ev(h)
            Psi[i, j] = h * Dk.ev(h)
    return Phi, Psi


def det4(m):
    from itertools import permutations
    tot = arb(0)
    for perm in permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if perm[i] > perm[j]:
                    sg = -sg
        term = m[0][perm[0]] * m[1][perm[1]] * m[2][perm[2]] * m[3][perm[3]]
        tot = tot + term if sg > 0 else tot - term
    return tot


class Pert2(PG.Pert):
    def F(self, kap, want_vectors=False):
        cfg = self.cfg
        Np = cfg['Np']
        st = self.steps[0]
        Ma, Mb, _ = self.mats[0]
        sig0 = st['y'][0].c[0]
        M = comb(Ma, Mb, kap)
        y0 = [arb(2) / 3, arb(1), -kap / (4 * sig0), arb(0)]
        yk0 = [arb(0), arb(0), -1 / (4 * sig0), arb(0)]
        z1 = st['h']
        ys, yks, _ = lin_singular(M, Mb, y0, yk0, st['T'], Np, z1)
        sc = [z1 * z1, z1 * z1, z1, z1 * z1]
        v = [ys[i] * sc[i] for i in range(4)] + [yks[i] * sc[i] for i in range(4)]
        c = [arb(x.mid()) for x in v]
        B = arb_mat(8, 8)
        for i in range(8):
            B[i, i] = 1
        r = [v[i] - c[i] for i in range(8)]
        for idx in range(1, len(self.steps) - 1):
            st = self.steps[idx]
            P0, P1, _ = self.mats[idx]
            Q = comb(P0, P1, kap)
            Phi, Psi = lin_regular_mat(Q, P1, st['T'], Np, st['h'])
            Phi8 = arb_mat(8, 8)
            for i in range(4):
                for j in range(4):
                    Phi8[i, j] = Phi[i, j]
                    Phi8[4 + i, 4 + j] = Phi[i, j]
                    Phi8[4 + i, j] = Psi[i, j]
            phic = Phi8 * arb_mat([[x] for x in c])
            c, B, r = lohner_update([phic[i, 0] for i in range(8)], Phi8 * B, r)
        v = set_eval(c, B, r)
        yc, ykc = v[:4], v[4:]
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
        # explicit Leibniz expansion (arb_mat.det falls back to a crude bound when an LU pivot contains 0)
        Fv = det4([[allc[j][i] for j in range(4)] for i in range(4)])
        dF = arb(0)
        for jj in range(4):
            dF += det4([[(allk[j][i] if j == jj else allc[j][i]) for j in range(4)] for i in range(4)])
        # pivot-independent normalisation: det[x, kernel basis(p)] = c_p w_p^2 (w.x), c_1 = -1, c_2 = +1
        cp = -1 if piv == 1 else 1
        wp, wbp = w[piv], wb[piv]
        Fn = Fv / (cp * wp * wp)
        dFn = dF / (cp * wp * wp) - 2 * Fv * wbp / (cp * wp * wp * wp)
        self.last_piv = piv
        if want_vectors:
            return Fn, dFn, lam, allc, allk
        return Fn, dFn, lam
