"""Rigorous analytic branch through a sonic equilibrium (parameterization method + weighted sup-norm Volterra
contraction).  tau-field F = (z P S, Om nO, (V^2-1) nV/(k+1)) (system2.py); equilibrium e(zs) on the sonic line;
J = DF(e) has eigenvalues {0 (line), lam_a, lam_b}.  For the chosen eigenvalue lam (>0 here):
    P(th) = e + th v + sum_{n>=2} p_n th^n,     F(P(th)) = lam th P'(th).
Head p_2..p_N: (J - n lam) p_n = -[F(P_<n)]_n  (exact via eigen-basis T).  Remainder rho = P - H satisfies
    lam th rho' - J rho = d + Nl(rho),  d = F(H) - lam th H',  Nl = F(H+rho) - F(H) - J rho,
    q = T^-1 rho :  q_i = K_i[(T^-1(d + Nl))_i],  K_i f(th) = (1/lam) int_0^1 s^(-kap_i/lam - 1) f(s th) ds,
    |K_i f(th)| <= ||f||_m |th|^m / (m lam - kap_i)    (m = N+1 > kap_i/lam),  ||f||_m = sup |f(th)|/|th|^m.
Contraction in the weighted norm on |th| <= th1 gives  |rho_j(th)| <= sum_i |T_ji| r_i |th|^m.
"""
import json
import sympy as sp
from flint import arb, arb_poly, ctx
from cheb import aup, hull
import systems as SY

_sym = json.load(open('<LAB>/b5/sonic_sym.json'))
q_, zs_ = sp.symbols('q zs', positive=True)
k, z, Om, V = SY.k, SY.z, SY.Om, SY.V
loc = {'q': q_, 'zs': zs_, 'k': k, 'z': z, 'Om': Om, 'V': V}


def L(args, s):
    e = sp.sympify(s, locals=loc) if isinstance(s, str) else s
    for at in e.atoms(sp.Rational):
        assert at.q == 1, at
    return sp.lambdify(args, e, modules=[{}])


Vs_f = L((q_, zs_), _sym['Vs'])
Oms_f = L((q_, zs_), _sym['Oms'])
tang_f = [L((q_, zs_), t) for t in _sym['tang']]
F_f = [L((k, z, Om, V), s) for s in _sym['F']]
J_f = [[L((k, z, Om, V), _sym['J'][i][j]) for j in range(3)] for i in range(3)]


class TP:
    """truncated power series via arb_poly (mullow)"""
    __slots__ = ('p', 'N')

    def __init__(self, p, N):
        self.p = p if isinstance(p, arb_poly) else arb_poly(p)
        self.N = N

    def _c(self, o):
        return o if isinstance(o, TP) else TP(arb_poly([arb(o)]), self.N)

    def __add__(self, o):
        return TP(self.p + self._c(o).p, self.N)
    __radd__ = __add__

    def __sub__(self, o):
        return TP(self.p - self._c(o).p, self.N)

    def __rsub__(self, o):
        return TP(self._c(o).p - self.p, self.N)

    def __neg__(self):
        return TP(-self.p, self.N)

    def __mul__(self, o):
        if isinstance(o, TP):
            return TP((self.p * o.p).truncate(self.N + 1), self.N)
        return TP(self.p * arb(o), self.N)
    __rmul__ = __mul__

    def __truediv__(self, o):
        return TP(self.p * (1 / arb(o)), self.N)

    def __pow__(self, e):
        r = TP(arb_poly([arb(1)]), self.N)
        b = self
        while e:
            if e & 1:
                r = r * b
            b = b * b
            e >>= 1
        return r

    def coef(self, n):
        c = self.p.coeffs()
        return c[n] if n < len(c) else arb(0)


def mat_inv3(A):
    a, b, c = A[0]
    d, e, f = A[1]
    g, h, i = A[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    adj = [[e * i - f * h, c * h - b * i, b * f - c * e],
           [f * g - d * i, a * i - c * g, c * d - a * f],
           [d * h - e * g, b * g - a * h, a * e - b * d]]
    return [[adj[r][s] / det for s in range(3)] for r in range(3)], det


def matvec(A, v):
    return [sum((A[i][j] * v[j] for j in range(3)), arb(0)) for i in range(3)]


def equilibrium(kk, zs):
    q = kk.sqrt()
    e = [zs, Oms_f(q, zs), Vs_f(q, zs)]
    J = [[arb(J_f[i][j](kk, *e)) for j in range(3)] for i in range(3)]
    t = [arb(tf(q, zs)) for tf in tang_f]
    return e, J, t


def eigen(J):
    """J has an exact zero eigenvalue (sonic_sym.py).  char poly = -x (x^2 - tr x + c2)."""
    tr = J[0][0] + J[1][1] + J[2][2]
    c2 = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) + (J[0][0] * J[2][2] - J[0][2] * J[2][0]) + \
         (J[1][1] * J[2][2] - J[1][2] * J[2][1])
    disc = tr * tr - 4 * c2
    return tr, c2, disc


def eigvec(J, lam):
    M = [[J[i][j] - (lam if i == j else 0) for j in range(3)] for i in range(3)]
    best = None
    for (r1, r2) in ((0, 1), (0, 2), (1, 2)):
        a, b = M[r1], M[r2]
        v = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        nrm = sum(float(abs(x).mid()) for x in v)
        if best is None or nrm > best[0]:
            best = (nrm, v)
    return best[1]


def classify(kk, zs):
    e, J, t = equilibrium(kk, zs)
    tr, c2, disc = eigen(J)
    out = {'tr': tr, 'c2': c2, 'disc': disc}
    if c2 < 0:
        out['type'] = 'saddle'
    elif c2 > 0 and disc > 0:
        out['type'] = 'node'
    elif disc < 0:
        out['type'] = 'focus'
    else:
        out['type'] = 'undecided'
    return out, e, J, t


def _select(out, which):
    tr, c2, disc = out['tr'], out['c2'], out['disc']
    sq = disc.sqrt()
    lp, lm = (tr + sq) / 2, (tr - sq) / 2
    if which == 'unstable':
        assert out['type'] == 'saddle'
        return lp, lm
    if which == 'strong':
        assert out['type'] == 'node' and lm > 0
        return lp, lm
    if which == 'weak':
        assert out['type'] == 'node' and lm > 0
        return lm, lp
    raise ValueError(which)


def branch(kk, zs, N=60, scale=None, which='unstable', head=None, norm=None):
    """Rigorous leading objects at the sonic equilibrium for zs (ball):
         e(zs), J, eigenvalues (closed form, J has exact eigenvalue 0), branch eigenvector v (cross product),
         eigen-basis T = [tangent, v, v_other].
       Head coefficients p_2..p_N: if head is None they are computed by the recursion (J - n lam) p_n = -[F]_n in
       ball arithmetic and then REPLACED BY THEIR MIDPOINTS (any polynomial is admissible: validate_tail bounds the
       true curve around it using the exact defect).  norm fixes the eigenvector normalisation (mid value)."""
    out, e, J, t = classify(kk, zs)
    lam, kap_other = _select(out, which)
    if which == 'weak':
        for n in range(2, N + 1):
            assert not (kap_other - n * lam).contains(0), 'resonance'
    v = eigvec(J, lam)
    vo = eigvec(J, kap_other)
    if norm is None:
        nv = max(abs(x) for x in v)
        norm = arb(aup(nv).mid())
    s = scale if scale is not None else arb(1)
    v = [x / norm * s for x in v]
    nvo = max(abs(x) for x in vo)
    vo = [x / arb(aup(nvo).mid()) for x in vo]
    T = [[t[i], v[i], vo[i]] for i in range(3)]
    kap = [arb(0), lam, kap_other]
    Ti, detT = mat_inv3(T)
    assert not detT.contains(0)
    if head is None:
        p = [[e[i], v[i]] for i in range(3)]
        for n in range(2, N + 1):
            P = [TP(arb_poly(p[i]), n) for i in range(3)]
            Fv = [f(kk, *P) for f in F_f]
            R = [Fv[i].coef(n) for i in range(3)]
            w = matvec(Ti, R)
            w = [-w[i] / (kap[i] - n * lam) for i in range(3)]
            pn = matvec(T, w)
            for i in range(3):
                p[i].append(pn[i])
        head = [[arb(x.mid()) for x in p[i][2:]] for i in range(3)]
    p = [[e[i], v[i]] + list(head[i]) for i in range(3)]
    return {'p': p, 'lam': lam, 'kap': kap, 'T': T, 'Ti': Ti, 'e': e, 'J': J, 'class': out, 'N': N,
            'head': head, 'norm': norm, 'which': which}


def validate_tail(B, th1, pieces=16, kk=None):
    """weighted sup-norm contraction for rho on |th| <= th1.  returns radii r (eigen-coords) and m."""
    p, lam, kap, T, Ti, J, N = B['p'], B['lam'], B['kap'], B['T'], B['Ti'], B['J'], B['N']
    m = 2      # rho = O(th^2): p_0 = e(zs), p_1 = v exact (balls); d_0 = F(e) = 0, d_1 = (J - lam) v = 0 exactly
    for i in range(3):
        assert (m * lam - kap[i]) > 0
    H = [arb_poly(p[i]) for i in range(3)]
    # exact defect polynomial d = F(H) - lam th H'
    DEG = 12 * N + 4
    Ht = [TP(H[i], DEG) for i in range(3)]
    Fv = [f(kk, *Ht).p for f in F_f]
    th = arb_poly([arb(0), arb(1)])
    d = [Fv[i] - lam * th * H[i].derivative() for i in range(3)]
    # d_0 = F(e(zs)) = 0 and d_1 = (J - lam I) v = 0 hold exactly (sonic_sym.py identities; v is a cross product
    # of two rows of J - lam I, which has rank 2 because lam is a simple nonzero eigenvalue): set them to 0.
    dc = [d[i].coeffs() for i in range(3)]
    for i in range(3):
        for n in range(min(2, len(dc[i]))):
            assert dc[i][n].contains(0), ('identity violated', i, n)
            dc[i][n] = arb(0)
    Td = [[sum((Ti[a][b] * (dc[b][n] if n < len(dc[b]) else arb(0)) for b in range(3)), arb(0))
           for n in range(max(len(x) for x in dc))] for a in range(3)]
    D = []
    for a in range(3):
        s = arb(0)
        for n in range(m, len(Td[a])):
            s += abs(Td[a][n]) * th1 ** (n - m)
        D.append(aup(s))
    Kn = [1 / (m * lam - kap[i]) for i in range(3)]
    absT = [[abs(T[i][j]) for j in range(3)] for i in range(3)]

    def lip(r):
        # A(th) = DF(H(th) + rho(th)) - J over |th| <= th1 ; |rho_j| <= sum_i |T_ji| r_i th1^m
        rb = [aup(sum((absT[j][i] * r[i] for i in range(3)), arb(0)) * th1 ** m) for j in range(3)]
        Lm = [[arb(0)] * 3 for _ in range(3)]
        first = True
        for sgn in (-1, 1):
            for pc in range(pieces):
                tb = hull(th1 * pc / pieces, th1 * (pc + 1) / pieces) * sgn
                X = [H[j](tb) + arb(0, rb[j]) for j in range(3)]
                A = [[arb(J_f[i][j](kk, *X)) - J[i][j] for j in range(3)] for i in range(3)]
                # eigen coords: Ti A T
                AT = [[sum((A[i][l] * T[l][j] for l in range(3)), arb(0)) for j in range(3)] for i in range(3)]
                TAT = [[sum((Ti[i][l] * AT[l][j] for l in range(3)), arb(0)) for j in range(3)] for i in range(3)]
                for i in range(3):
                    for j in range(3):
                        val = arb(aup(abs(TAT[i][j])))
                        Lm[i][j] = val if first else hull(Lm[i][j], val)
                first = False
        return [[aup(Lm[i][j]) for j in range(3)] for i in range(3)]

    tiny = arb(10) ** (-(ctx.prec * 3 // 10))
    r = [aup(Kn[i] * D[i] * arb('1.0001') + tiny) for i in range(3)]
    for it in range(40):
        Lm = lip(r)
        lhs = [aup(Kn[i] * (D[i] + sum((Lm[i][j] * r[j] for j in range(3)), arb(0)))) for i in range(3)]
        if all(lhs[i] < r[i] for i in range(3)):
            rho = [aup(sum((absT[j][i] * r[i] for i in range(3)), arb(0))) for j in range(3)]
            FH1 = []
            for i in range(3):
                cc = Fv[i].coeffs()
                assert cc[0].contains(0)          # F(e) = 0 exactly
                FH1.append(arb_poly(cc[1:]))       # F(H(th))/th
            return {'r': r, 'rho_coef': rho, 'm': m, 'D': D, 'L': Lm, 'th1': th1, 'FH1': FH1}
        r = [aup(lhs[i] * arb('1.05') + tiny) for i in range(3)]
    raise ValueError('tail contraction failed: D=%s KL=%s' % ([float(x.mid()) for x in D], [[float((Kn[i] * Lm[i][j]).mid()) for j in range(3)] for i in range(3)]))


def P_eval(B, V_, thv):
    """enclosure of P(th) for a ball th with |th| <= th1"""
    H = [arb_poly(B['p'][i]) for i in range(3)]
    m = V_['m']
    w = abs(thv) ** m
    return [H[i](thv) + arb(0, aup(V_['rho_coef'][i] * aup(w))) for i in range(3)]


def dP_eval(B, V_, thv, kk):
    """P'(th) = F(P(th))/(lam th) = [F(H)/th](th)/lam + (F(H+rho)-F(H))/(lam th),
    |F(H+rho)-F(H)| <= |DF(P box)| |rho| <= |DF| rho_coef th^2   (th ball, may contain 0)."""
    Pv = P_eval(B, V_, thv)
    DFb = [[abs(arb(J_f[i][j](kk, *Pv))) for j in range(3)] for i in range(3)]
    at = abs(thv)
    out = []
    for i in range(3):
        bnd = sum((DFb[i][j] * V_['rho_coef'][j] for j in range(3)), arb(0)) * at
        out.append((V_['FH1'][i](thv) + arb(0, aup(bnd))) / B['lam'])
    return out


def solve_theta(B, V_, kk, zt, th_guess, width):
    """interval Newton for P_z(th) = zt in th_guess +- width.  proven=True once N(X) lies strictly inside X
    (then X contains exactly one root); afterwards X is replaced by N(X) (still contains the root)."""
    X = arb(th_guess, width)
    proven = False
    for it in range(40):
        c = arb(X.mid())
        fc = P_eval(B, V_, c)[0] - zt
        dX = dP_eval(B, V_, X, kk)[0]
        if dX.contains(0):
            raise ValueError('derivative of P_z contains 0')
        Nx = c - fc / dX
        if Nx.lower() > X.lower() and Nx.upper() < X.upper():
            proven = True
        if not Nx.overlaps(X):
            raise ValueError('no root in box')
        Xn = Nx.intersection(X)
        if float(Xn.rad()) >= 0.999 * float(X.rad()) and proven:
            return Xn, proven
        X = Xn
    return X, proven
