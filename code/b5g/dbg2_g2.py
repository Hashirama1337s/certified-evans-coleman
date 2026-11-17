import sys
sys.path.insert(0, '<LAB>/b5g')
from flint import arb, ctx, arb_series, arb_poly
from tm import TM
ctx.prec = 200
n = 26
b = [arb(-3), arb(0), arb('18.55')] + [arb(0)] * (n - 3)
q = (arb_series([arb(1)] + [arb(0)] * (n - 1), prec=n) / arb_series(b, prec=n)).coeffs()
print('div series', ['%.1e/%.1e' % (float(abs(x).mid()), float(x.rad())) for x in q[::4]])
# repeated products
T = arb(0.025, 0.025)
a = TM([arb(1), arb(-1) / 2] + [arb(1) / 7] * (n - 1), T)
p = a
for k in range(6):
    p = p * a
print('prod', ['%.1e/%.1e' % (float(abs(x).mid()), float(x.rad())) for x in p.c[::4]])
d = TM([arb(-3), arb(0), arb('18.55')] + [arb(0)] * (n - 2), T)
r = a / d
print('tmdiv', ['%.1e/%.1e' % (float(abs(x).mid()), float(x.rad())) for x in r.c[::4]])
