"""REFEREE-B4FW item B2 (NEW file): rewrite lin_fw4_bg.pkl as lin_fw4_bg3.pkl with radii rounded UP.
lin_fw4_bg.py stored (mid to 70 significant digits, rad to 10 significant digits, round-to-nearest).
  stored rad  >= true rad * (1 - 5e-10)        -> new rad = stored rad * (1 + 1e-9)
  stored mid differs from true mid by <= |mid| * 5e-70  -> add |mid| * 1e-69
Every (mid, rad) pair becomes (mid, rad3) with rad3 an upward-rounded decimal (arb upper bound, 20 digits).
usage: lin_fw4_bg3_make.py"""
import pickle, sys, decimal
decimal.getcontext().prec = 120
from flint import arb, ctx
ctx.prec = 256
src = '<LAB>/engine/cert/lin_fw4_bg.pkl'
dst = '<LAB>/engine/cert/lin_fw4_bg3.pkl'
cnt = [0]
def fix(p):
    m, r = p
    if r == '0' and decimal.Decimal(m) == decimal.Decimal(m).to_integral_value():
        cnt[0] += 1
        return (m, '0')                                   # exact integer ball (e.g. centre-manifold P_0 = 0, P_1 = 1)
    rr =arb(r)*(1 + arb('1e-9')) + abs(arb(m))*arb('1e-69') + arb('1e-300')
    up = rr.upper()                                      # exact binary point >= rr
    d = decimal.Decimal(up.mid().str(80, radius=False).replace('e', 'E'))
    q = decimal.Decimal(1).scaleb(d.adjusted() - 19)
    s = str(d.quantize(q, rounding=decimal.ROUND_CEILING) + q)   # decimal strictly above up
    cnt[0] += 1
    return (m, s)
def walk(o):
    if isinstance(o, tuple) and len(o) == 2 and all(isinstance(t, str) for t in o):
        return fix(o)
    if isinstance(o, tuple) and len(o) == 4 and isinstance(o[0], float):      # rec entry (x, h, Bx, st)
        return (o[0], walk(o[1]), walk(o[2]), walk(o[3]))
    if isinstance(o, list):
        return [walk(t) for t in o]
    if isinstance(o, dict):
        return {k: walk(v) for k, v in o.items()}
    return o
d = pickle.load(open(src, 'rb'))
d3 = walk(d)
d3['note'] = d['note'] + ' | bg3: radii rounded up (x(1+1e-9) + |mid| 1e-69), lin_fw4_bg3_make.py'
pickle.dump(d3, open(dst, 'wb'))
# check: every new ball contains the old ball
def chk(a, b):
    if isinstance(a, tuple) and len(a) == 2 and all(isinstance(t, str) for t in a):
        A = arb(a[0]) + arb(0, arb(a[1])); B = arb(b[0]) + arb(0, arb(b[1]))
        assert B.contains(A), (a, b)
        return 1
    if isinstance(a, (list, tuple)):
        return sum(chk(x, y) for x, y in zip(a, b))
    if isinstance(a, dict):
        return sum(chk(a[k], b[k]) for k in a)
    return 0
print('pairs rewritten', cnt[0], ' containment checks passed', chk(d, d3))
print('wrote', dst)
