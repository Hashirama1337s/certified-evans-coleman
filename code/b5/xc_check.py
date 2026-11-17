import json
from flint import arb, ctx
ctx.prec = 400
x = json.load(open('crosscheck50_ec.json')); r = json.load(open('result_ec.json'))
for a, b in (('t0_m_over_r', 't0_m_over_r'), ('t0_Om', 't0_Om'), ('t0_V', 't0_V')):
    print(a, 'cross-check value inside certified enclosure:', arb(r[b]).contains(arb(x[a])), ' |diff to mid| =', abs(arb(x[a]) - arb(arb(r[b]).mid())).str(3))
