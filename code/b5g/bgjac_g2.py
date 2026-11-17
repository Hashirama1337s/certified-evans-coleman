"""Export background RHS f(Om,V,z) = (dOm/dz, dV/dz) and its Jacobian (for the Lohner mean-value propagation)."""
import json
import sympy as sp
from system_g2 import FOmx, FVx, MUalg, SL, mu, Om, V, N

z = sp.symbols('z')
zx = (-z * SL).subs(mu, MUalg).subs(N, 1 / z)
f0 = sp.factor(sp.cancel((FOmx.subs(mu, MUalg).subs(N, 1 / z)) / zx))
f1 = sp.factor(sp.cancel((FVx.subs(mu, MUalg).subs(N, 1 / z)) / zx))
J = [[sp.factor(sp.cancel(sp.diff(f, s))) for s in (Om, V)] for f in (f0, f1)]
print('f0 =', f0)
print('f1 =', f1)
json.dump({'f0': str(f0), 'f1': str(f1), 'J': [[str(e) for e in r] for r in J]},
          open('<LAB>/b5g/bgjac_g2.json', 'w'), indent=1)
print('written bgjac_g2.json')
