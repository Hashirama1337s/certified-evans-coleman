"""Symbolic checks + expression export for the arb pipeline.
1. sonic line e(zs): Sz(V_s, zs) = 0, Om_s from the V-numerator, and the identity that both background numerators
   (GO, GV) vanish on e(zs) for every zs (exact, sympy) -> justifies zeroing constant Taylor terms at the sonic point.
2. FRW (k=1/3) exact solution V = -(1-sqrt(1-z^2))/z, Om = 3V^2/2 satisfies the background system (exact).
3. FRW sonic point zs = sqrt(3)/2 is on the sonic line with Om_s = 1/2.
4. export all expression strings needed by pipeline_g2.py to exprs_g2.json.
"""
import json
import sympy as sp
from system_g2 import FOmx, FVx, MUalg, SL, H, fl, Ax, At, mu, Om, V, N, Om_x, V_x, Om_t, V_t, kap

z = sp.symbols('z', positive=True)
zs = sp.symbols('zs', positive=True)
Sx = 3 * N**2 * V**2 - N**2 + 4 * N * V - V**2 + 3
Sz = sp.expand(Sx * z**2).subs(N, 1 / z)
Sz = sp.expand(sp.simplify(Sz))
print('Sz =', Sz)
GO = sp.factor(sp.cancel(FOmx * Sx))
GV = sp.factor(sp.cancel(FVx * Sx))
print('GO =', GO)
print('GV =', GV)
# sonic line
Vs = (1 - sp.sqrt(3) * zs) / (sp.sqrt(3) - zs)
print('Sz(Vs,zs) =', sp.simplify(Sz.subs({V: Vs, z: zs})))
GVz = sp.together(GV.subs(mu, MUalg).subs(N, 1 / z))
numV = sp.factor(sp.fraction(GVz)[0])
print('GV numerator (mu algebraic) =', numV)
# linear in Om ?
lin = sp.Poly(sp.expand(sp.fraction(sp.cancel(GVz / (V**2 - 1)))[0]), Om)
print('degree of GV/(V^2-1) numerator in Om:', lin.degree())
sol = sp.solve(sp.fraction(sp.cancel(GVz / ((V - 1) * (V + 1))))[0], Om)
assert len(sol) == 1
Om_s_expr = sp.simplify(sol[0])
print('Om_s(V, z) =', Om_s_expr)
Oms = sp.simplify(Om_s_expr.subs({V: Vs, z: zs}))
print('Om_s(zs) =', Oms)
GOz = sp.together(GO.subs(mu, MUalg).subs(N, 1 / z))
chkO = sp.simplify(GOz.subs({Om: Oms, V: Vs, z: zs}))
chkV = sp.simplify(GVz.subs({Om: Oms, V: Vs, z: zs}))
print('GO(e(zs)) =', chkO, '  GV(e(zs)) =', chkV)
assert chkO == 0 and chkV == 0
# FRW
Vf = -(1 - sp.sqrt(1 - z**2)) / z
Omf = sp.Rational(3, 2) * Vf**2
zx = (-z * SL).subs(mu, MUalg).subs(N, 1 / z)
dOm = (FOmx.subs(mu, MUalg).subs(N, 1 / z) / zx)
dV = (FVx.subs(mu, MUalg).subs(N, 1 / z) / zx)
r1 = sp.simplify(sp.diff(Omf, z) - dOm.subs({Om: Omf, V: Vf}))
r2 = sp.simplify(sp.diff(Vf, z) - dV.subs({Om: Omf, V: Vf}))
r3 = sp.simplify(MUalg.subs({Om: Omf, V: Vf, N: 1 / z}) - Vf**2)
print('FRW residuals dOm/dz, dV/dz, mu-V^2:', r1, r2, r3)
assert r1 == 0 and r2 == 0 and r3 == 0
zsF = sp.sqrt(3) / 2
print('FRW at zs=sqrt3/2: V =', sp.simplify(Vf.subs(z, zsF)), ' Vs =', sp.simplify(Vs.subs(zs, zsF)),
      ' Om =', sp.simplify(Omf.subs(z, zsF)), ' Om_s =', sp.simplify(Oms.subs(zs, zsF)))
# exports
pd = {}
args6 = (mu, Om, V, N, Om_x, V_x)
for i in range(2):
    for j, s in enumerate((Om_x, V_x)):
        pd['Ax%d%d' % (i, j)] = sp.expand(sp.diff(fl[i], s))
    for j, s in enumerate((Om_t, V_t)):
        pd['At%d%d' % (i, j)] = sp.expand(sp.diff(fl[i], s))
    for j, s in enumerate((Om, V)):
        pd['Bu%d%d' % (i, j)] = sp.expand(sp.diff(fl[i], s).subs({Om_t: 0, V_t: 0}))
    pd['Bmu%d' % i] = sp.expand(sp.diff(fl[i], mu).subs({Om_t: 0, V_t: 0}))
    pd['BNN%d' % i] = sp.expand((sp.diff(fl[i], N) * N).subs({Om_t: 0, V_t: 0}))
for j, s in enumerate((mu, Om, V)):
    pd['gH%d' % j] = sp.factor(sp.diff(H, s))
    pd['gS%d' % j] = sp.factor(sp.diff(SL, s))
pd['MU'] = sp.factor(MUalg)
pd['SL'] = sp.factor(SL)
pd['GO'] = GO
pd['GV'] = GV
pd['Sz'] = Sz
pd['detfac'] = 36 * Om * (V**2 - 1)**2 * (mu - 1)**2
# check det Ax == detfac * Sx (with any Om_x, V_x)
dchk = sp.expand(Ax.det() - pd['detfac'] * Sx)
print('det Ax - detfac*Sx =', dchk)
assert dchk == 0
pd['Om_s'] = Om_s_expr
pd['H'] = H
pd['Mo'] = sp.sympify(json.load(open('<LAB>/b5g/eqs_g2.json'))['Mo'],
                      locals=dict(mu=mu, Om=Om, V=V, N=N, k=sp.Symbol('k', positive=True))).subs(sp.Symbol('k', positive=True), sp.Rational(1, 3))
c = json.load(open('<LAB>/b5g/centre_sys.json'))
out = {kk: str(vv) for kk, vv in pd.items()}
out['centre'] = c
json.dump(out, open('<LAB>/b5g/exprs_g2.json', 'w'), indent=1)
print('written exprs_g2.json;  SONIC/FRW SYMBOLIC CHECKS OK')
