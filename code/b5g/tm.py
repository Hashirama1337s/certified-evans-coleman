"""Taylor models with interval coefficients (own implementation, python-flint arb).

A TM of order n on a domain T (an arb ball containing 0 and every t used) is a list c[0..n] of arb:
    f(t) = sum_{i<n} c[i] t^i + t^n r(t),   r(t) in c[n] for every t in T (pointwise).
Arithmetic keeps the low coefficients as the exact (ball-enclosed) Taylor coefficients of the true functions and
folds everything else into the top coefficient pointwise (Horner over T), so the remainder of a result is an
enclosure relative to the TRUE low coefficients (needed by the causal bootstrap in solver code).
"""
import flint
from flint import arb, arb_poly, arb_series, ctx

ctx.cap = 80          # above every series length used (orders <= 40); default 10 would silently truncate


def horner(cs, t):
    v = arb(0)
    for c in reversed(cs):
        v = v * t + c
    return v


class TM:
    __slots__ = ('c', 'T')

    def __init__(self, c, T):
        self.c = list(c)
        self.T = T

    @property
    def n(self):
        return len(self.c) - 1

    @staticmethod
    def const(v, n, T):
        return TM([arb(v)] + [arb(0)] * n, T)

    def lower(self, n):
        if n >= self.n:
            return self
        return TM(self.c[:n] + [horner(self.c[n:], self.T)], self.T)

    def __neg__(self):
        return TM([-x for x in self.c], self.T)

    def _co(self, o):
        if isinstance(o, TM):
            return o
        return None

    def __add__(self, o):
        if not isinstance(o, TM):
            c = list(self.c)
            c[0] = c[0] + o
            return TM(c, self.T)
        n = min(self.n, o.n)
        a, b = self.lower(n), o.lower(n)
        return TM([x + y for x, y in zip(a.c, b.c)], self.T)

    __radd__ = __add__

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, TM):
            return TM([x * o for x in self.c], self.T)
        n = min(self.n, o.n)
        a, b = self.lower(n), o.lower(n)
        pc = (arb_poly(a.c) * arb_poly(b.c)).coeffs()
        pc = pc + [arb(0)] * (2 * n + 1 - len(pc))
        return TM(pc[:n] + [horner(pc[n:], self.T)], self.T)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, TM):
            return TM([x / o for x in self.c], self.T)
        n = min(self.n, o.n)
        a, b = self.lower(n), o.lower(n)
        if n == 0:
            br = b.c[0]
            if br.contains(0):
                raise ZeroDivisionError('TM division: denominator range contains 0')
            return TM([a.c[0] / br], self.T)
        q = (arb_series(a.c[:n], prec=n) / arb_series(b.c[:n], prec=n)).coeffs()
        q = q + [arb(0)] * (n - len(q))
        qt = TM(q + [arb(0)], self.T)
        E = a - qt * b
        br = horner(b.c, self.T)
        if br.contains(0):
            raise ZeroDivisionError('TM division: denominator range contains 0')
        return TM(q + [E.c[n] / br], self.T)

    def __rtruediv__(self, o):
        return TM.const(o, self.n, self.T) / self

    def __pow__(self, e):
        assert isinstance(e, int) and e >= 0
        r = None
        b = self
        while e:
            if e & 1:
                r = b if r is None else r * b
            e >>= 1
            if e:
                b = b * b
        return r if r is not None else TM.const(1, self.n, self.T)

    def ev(self, t):
        return horner(self.c, t)

    def rng(self):
        return horner(self.c, self.T)

    def mul_t(self):
        """t*f: exact, order+1"""
        return TM([arb(0)] + self.c, self.T)

    def div_t_zero(self):
        """(f - f(0))/t, with the constant term DECLARED exactly zero by an exact identity (caller's duty)."""
        return TM(self.c[1:], self.T)

    def drop0(self):
        """(f - c0)/t (exact identity, no assumption): order-1"""
        return TM(self.c[1:], self.T)

    def average(self):
        """A[f](t) = (1/t) int_0^t f ; top r -> hull(r)/(n+1) (valid pointwise, causal)."""
        n = self.n
        return TM([self.c[i] / (i + 1) for i in range(n)] + [self.c[n] / (n + 1)], self.T)


def R(p, q):
    return arb(p) / q


def codegen(exprs, args, name='f'):
    """sympy -> python source; Rationals printed as R(p,q) (exact arb), cse for speed."""
    import sympy as sp
    from sympy.printing.str import StrPrinter

    class P(StrPrinter):
        def _print_Rational(self, e):
            return 'R(%d,%d)' % (e.p, e.q)

        def _print_Integer(self, e):
            return str(int(e))

        def _print_Pow(self, e, rational=False):
            b, x = e.as_base_exp()
            if x.is_Integer:
                xi = int(x)
                if xi >= 0:
                    return '(%s)**%d' % (self._print(b), xi)
                return '(1/((%s)**%d))' % (self._print(b), -xi)
            if x == sp.Rational(1, 2):
                return 'SQRT(%s)' % self._print(b)
            raise ValueError('non-integer power %s' % e)

        def _print_Float(self, e):
            raise ValueError('float in expression')

    pr = P()
    reps, red = sp.cse(exprs)
    lines = ['def %s(%s):' % (name, ', '.join(str(a) for a in args))]
    for s, e in reps:
        lines.append('    %s = %s' % (s, pr.doprint(e)))
    lines.append('    return [%s]' % ', '.join(pr.doprint(e) for e in red))
    return '\n'.join(lines)
