"""Symbolic facts at the sonic line (exact sympy; identities proved by cancelling to a zero numerator):
  sonic branch S=0 with (V+z) = q(1+Vz), q = sqrt(k)  ->  V_s(z_s) = (q - z_s)/(1 - q z_s)
  nO = 0 is linear in Om  ->  Om_s(z_s)
  identities: S(e)=0, nV(e(z_s)) == 0, det J(e(z_s)) == 0, J t == 0 (t = de/dz_s): line of equilibria.
Writes sonic_sym.json first (phase 1), then checks identities (phase 2)."""
import json, sys
import sympy as sp
import systems as SY

k, z, Om, V = SY.k, SY.z, SY.Om, SY.V
q, zs = sp.symbols('q zs', positive=True)
Vs = (q - zs) / (1 - q * zs)
nO = sp.expand(SY.nO)
c1 = nO.coeff(Om, 1)
c0 = nO.coeff(Om, 0)
assert sp.expand(nO - c1 * Om - c0) == 0
Oms = sp.cancel((-c0 / c1).subs({V: Vs, z: zs, k: q**2}))
F = SY.tau_field()
X = [z, Om, V]
J = sp.Matrix([[sp.diff(f, x) for x in X] for f in F])
tang = [sp.Integer(1), sp.cancel(sp.diff(Oms, zs)), sp.cancel(sp.diff(Vs, zs))]
json.dump({'Vs': str(Vs), 'Oms': str(Oms), 'tang': [str(t) for t in tang],
           'J': [[str(J[i, j]) for j in range(3)] for i in range(3)],
           'F': [str(f) for f in F]}, open('<LAB>/b5/sonic_sym.json', 'w'), indent=1)
print('phase 1 written; Om_s =', sp.factor(Oms), flush=True)


def is_zero(e):
    n = sp.numer(sp.together(e))
    return sp.expand(n) == 0


sub = {V: Vs, z: zs, Om: Oms, k: q**2}
print('S(e)=0:', is_zero(SY.S.subs(sub)), flush=True)
print('nV(e)=0:', is_zero(SY.nV.subs(sub)), flush=True)
Je = J.subs(sub)
print('J t = 0:', all(is_zero(x) for x in (Je * sp.Matrix(tang))), flush=True)
# det J = 0 follows from J t = 0 with t != 0 (t[0] = 1).
print('det J = 0: implied by J t = 0, t != 0', flush=True)
