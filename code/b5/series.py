"""Truncated power-series arithmetic (coefficient lists) usable with mpmath mpf, python-flint arb/acb.
Also loads the polynomial right-hand sides from system2.json as python functions."""
import json
import sympy as sp


class TS:
    __slots__ = ('c', 'N')

    def __init__(self, c, N):
        self.c = list(c)[:N + 1]
        self.N = N
        while len(self.c) < N + 1:
            self.c.append(self.c[0] * 0 if self.c else 0)

    def _co(self, o):
        if isinstance(o, TS):
            return o
        z = self.c[0] * 0
        return TS([z + o] + [z] * self.N, self.N)

    def __add__(self, o):
        o = self._co(o)
        return TS([a + b for a, b in zip(self.c, o.c)], self.N)
    __radd__ = __add__

    def __neg__(self):
        return TS([-a for a in self.c], self.N)

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) + (-self)

    def __mul__(self, o):
        if not isinstance(o, TS):
            return TS([a * o for a in self.c], self.N)
        N = self.N
        out = []
        for n in range(N + 1):
            s = self.c[0] * 0
            for i in range(n + 1):
                s += self.c[i] * o.c[n - i]
            out.append(s)
        return TS(out, N)
    __rmul__ = __mul__

    def __truediv__(self, o):
        assert not isinstance(o, TS)
        return TS([a / o for a in self.c], self.N)

    def __pow__(self, e):
        assert isinstance(e, int) and e >= 0
        r = self._co(1)
        b = self
        while e:
            if e & 1:
                r = r * b
            b = b * b
            e >>= 1
        return r


def load_funcs():
    d = json.load(open('<LAB>/b5/system2.json'))
    k, z, Om, V, sig, U, N = sp.symbols('k z Om V sig U N')
    loc = {'k': k, 'z': z, 'Om': Om, 'V': V, 'sig': sig, 'U': U, 'N': N}
    out = {}
    for nm in ('nO', 'nV', 'P', 'S', 'mu'):
        e = sp.sympify(d[nm], locals=loc)
        out[nm] = sp.lambdify((k, z, Om, V), e, modules=[{}])
    for nm in ('dsig', 'dU'):
        e = sp.sympify(d[nm], locals=loc)
        out[nm] = sp.lambdify((k, z, sig, U), e, modules=[{}])
        nn, dd = sp.fraction(sp.factor(e))
        if nm == 'dU':
            dd = sp.cancel(dd / z)   # dU = nn/(z*dd)
        out[nm + '_num'] = sp.lambdify((k, z, sig, U), sp.expand(nn), modules=[{}])
        out[nm + '_den'] = sp.lambdify((k, z, sig, U), sp.expand(dd), modules=[{}])
        out[nm + '_numexpr'] = sp.expand(nn)
        out[nm + '_denexpr'] = sp.expand(dd)
    for nm in ('dOm_dN', 'dV_dN'):
        e = sp.sympify(d[nm], locals=loc)
        out[nm] = sp.lambdify((k, N, Om, V), e, modules=[{}])
    return out
