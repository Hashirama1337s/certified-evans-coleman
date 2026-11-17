"""Chebyshev models with python-flint arb: polynomial in Chebyshev basis on an interval + sup-norm remainder.

A model  f(x) = sum_k c_k T_k(t) + E,  t = (2x - a - b)/(b - a) in [-1,1],  |E(x)| <= err for all x in [a,b].
All arithmetic is rigorous: c_k are arb balls, products are exact (arb_poly) and truncation moves the dropped
coefficients' l1 mass into err.  sup|sum c_k T_k| <= sum |c_k| (l1 bound) and a subdivided Clenshaw range is
available for tighter ranges.
"""
from flint import arb, arb_poly, ctx


def hull(a, b):
    """smallest-ish ball containing balls a and b (arb.union in this build over-widens)"""
    lo = a.lower() if a.lower() < b.lower() else b.lower()
    hi = a.upper() if a.upper() > b.upper() else b.upper()
    mid = (lo + hi) / 2
    rad = (hi - lo) / 2
    return mid + arb(0, aup(rad + mid.rad()))


def l1(cs):
    s = arb(0)
    for c in cs:
        s += abs(c)
    return s.upper() if hasattr(s, 'upper') else s


def aup(x):
    """upper bound of an arb as arb (exact)"""
    return arb(x.mid() + x.rad())


class CM:
    __slots__ = ('c', 'err', 'M')

    def __init__(self, c, err=None, M=None):
        self.c = list(c)
        self.err = arb(0) if err is None else err
        self.M = M if M is not None else len(self.c) - 1

    @staticmethod
    def const(v, M):
        return CM([arb(v)] + [arb(0)] * M, arb(0), M)

    def copy(self):
        return CM(list(self.c), self.err, self.M)

    def l1(self):
        s = arb(0)
        for c in self.c:
            s += abs(c)
        return aup(s)

    def supbound(self):
        return aup(self.l1() + self.err)

    def trunc(self, M=None):
        M = self.M if M is None else M
        if len(self.c) <= M + 1:
            return CM(self.c + [arb(0)] * (M + 1 - len(self.c)), self.err, M)
        extra = arb(0)
        for c in self.c[M + 1:]:
            extra += abs(c)
        return CM(self.c[:M + 1], aup(self.err + extra), M)

    def __add__(self, o):
        if not isinstance(o, CM):
            c = list(self.c)
            c[0] = c[0] + o
            return CM(c, self.err, self.M)
        n = max(len(self.c), len(o.c))
        c = [(self.c[i] if i < len(self.c) else 0) + (o.c[i] if i < len(o.c) else 0) for i in range(n)]
        c = [arb(x) for x in c]
        return CM(c, aup(self.err + o.err), max(self.M, o.M)).trunc()
    __radd__ = __add__

    def __neg__(self):
        return CM([-x for x in self.c], self.err, self.M)

    def __sub__(self, o):
        return self + (-o if isinstance(o, CM) else -arb(o))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, CM):
            o = arb(o)
            return CM([x * o for x in self.c], aup(self.err * abs(o)), self.M)
        M = max(self.M, o.M)
        n1, n2 = len(self.c) - 1, len(o.c) - 1
        # Laurent trick: f(cos th) = sum_{j=-n}^{n} d_j e^{ij th}, d_0 = c_0, d_{+-j} = c_j/2
        p1 = [self.c[abs(j)] * (1 if j == 0 else arb(0.5)) for j in range(-n1, n1 + 1)]
        p2 = [o.c[abs(j)] * (1 if j == 0 else arb(0.5)) for j in range(-n2, n2 + 1)]
        pr = arb_poly(p1) * arb_poly(p2)
        coeffs = pr.coeffs()
        N = n1 + n2
        while len(coeffs) < 2 * N + 1:
            coeffs.append(arb(0))
        c = [coeffs[N]] + [2 * coeffs[N + j] for j in range(1, N + 1)]
        err = aup(self.l1() * o.err + o.l1() * self.err + self.err * o.err)
        return CM(c, err, M).trunc()
    __rmul__ = __mul__

    def __pow__(self, e):
        assert isinstance(e, int) and e >= 0
        r = CM.const(1, self.M)
        b = self
        while e:
            if e & 1:
                r = r * b
            b = b * b
            e >>= 1
        return r

    def __truediv__(self, o):
        assert not isinstance(o, CM)
        return self * (1 / arb(o))

    # ---- calculus on [-1,1] in t; caller supplies the half-length h=(b-a)/2 ----
    def integral_from_left(self, h):
        """G(t) = h * int_{-1}^t f(s) ds  (as a model); err -> err * 2h"""
        c = self.c + [arb(0), arb(0)]
        n = len(self.c) - 1
        g = [arb(0)] * (n + 2)
        # standard formulas: int T0 = T1 ; int T1 = T2/4 (+const) ; int Tk = T_{k+1}/(2(k+1)) - T_{k-1}/(2(k-1))
        g[1] += c[0]
        if n >= 1:
            g[2] += c[1] / 4
        for kk in range(2, n + 1):
            g[kk + 1] += c[kk] / (2 * (kk + 1))
            g[kk - 1] -= c[kk] / (2 * (kk - 1))
        # fix constant so that G(-1) = 0 : T_j(-1) = (-1)^j
        val = arb(0)
        for j in range(1, n + 2):
            val += g[j] * (1 if j % 2 == 0 else -1)
        g[0] = -val
        g = [x * h for x in g]
        return CM(g, aup(self.err * 2 * abs(h)), self.M).trunc()

    def eval_t(self, t):
        """Clenshaw at t (arb ball allowed); remainder added"""
        b1 = arb(0)
        b2 = arb(0)
        for c in reversed(self.c[1:]):
            b1, b2 = 2 * t * b1 - b2 + c, b1
        v = t * b1 - b2 + self.c[0]
        return v + arb(0, self.err) if self.err != 0 else v

    def right(self):
        s = arb(0)
        for c in self.c:
            s += c
        return s + arb(0, self.err) if self.err != 0 else s

    def left(self):
        s = arb(0)
        for j, c in enumerate(self.c):
            s += c if j % 2 == 0 else -c
        return s + arb(0, self.err) if self.err != 0 else s

    def range(self, pieces=8):
        """enclosure of the range over t in [-1,1]: t = cos(th), th in [0,pi] split in pieces,
        T_k(t) = cos(k th) enclosed by arb cos on the th-interval (no Clenshaw blow-up)."""
        out = None
        pi = arb.pi()
        for i in range(pieces):
            th = hull(pi * i / pieces, pi * (i + 1) / pieces)
            v = self.c[0]
            for kk in range(1, len(self.c)):
                v += self.c[kk] * (kk * th).cos()
            out = v if out is None else hull(out, v)
        if self.err != 0:
            out = out + arb(0, self.err)
        return out


def cheb_nodes(M):
    """Chebyshev points of the first kind t_j = cos(pi (j+1/2)/(M+1)), j=0..M (as arb)"""
    from flint import arb as A
    pi = A.pi()
    return [((pi * (j + A(0.5))) / (M + 1)).cos() for j in range(M + 1)]


def interp(values, M):
    """Chebyshev coefficients interpolating values at first-kind nodes (non-rigorous use only: mids)"""
    from flint import arb as A
    pi = A.pi()
    n = M + 1
    c = []
    for kk in range(n):
        s = A(0)
        for j in range(n):
            s += values[j] * ((pi * kk * (j + A(0.5))) / n).cos()
        s = s * 2 / n
        if kk == 0:
            s = s / 2
        c.append(A(s.mid()))
    return c
