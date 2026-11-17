"""Exact flat-FRW (k=1/3) background as Taylor models (used ONLY at the FRW sonic point, where the
background slope equation has a double root and the Taylor recursion is resonant, so the analytic solution through
e(zs) is not isolated and cannot be shot; the closed form V = -(1 - sqrt(1-z^2))/z, Om = 3 V^2 / 2 is used there).
Checked against the integrated background at the matching point (overlap test in the pipeline)."""
from flint import arb, arb_series, ctx
from tm import TM, horner

ctx.cap = 80


def tm_sqrt(a):
    n = a.n
    q = arb_series(a.c[:n], prec=n).sqrt().coeffs()
    q = q + [arb(0)] * (n - len(q))
    qt = TM(q + [arb(0)], a.T)
    E = a - qt * qt
    ra = a.rng()
    if not ra > 0:
        raise ValueError('sqrt of non-positive range')
    den = ra.sqrt() + horner(q, a.T)
    if den.contains(0):
        raise ZeroDivisionError
    return TM(q + [E.c[n] / den], a.T)


def frw_sonic_bg(T, zs, N):
    t = TM([arb(0), arb(1)] + [arb(0)] * (N), T)
    z = t + zs
    w = tm_sqrt(1 - z * z)
    V = (w - 1) / z
    Vp = (1 - w - z * z / w) / (z * z)
    Om = V * V * 3 / 2
    Omp = V * Vp * 3
    return [Om, V], [Omp.lower(N), Vp.lower(N)]
