import json
from flint import arb, ctx
ctx.prec = 300
sig0 = arb(8) / 27
zs = 9 * arb(2).sqrt() / 14
json.dump({'k': '1/2', 'sig0': sig0.str(80), 'zs': zs.str(80), 'note': 'exact flat FRW (k=1/2): sig0 = (3/2) beta^2 = 8/27, zs = 9 sqrt(2)/14'},
          open('<LAB>/b5/hp_frw_1_2.json', 'w'), indent=1)
print(sig0.str(30), zs.str(30))
