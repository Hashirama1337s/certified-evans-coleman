"""B5: reduced polynomial forms.
From system.py:  Om_x = Om*nO/(S*Q), V_x = (V^2-1)*nV/((k+1)*S*Q), z_x = z*P/Q,
with S = k(1+Vz)^2 - (V+z)^2 (sonic factor), Q = -z(1-V^2)(1-mu).
Desingularised field (dx/dtau = S*Q):
   z' = z*P*S,  Om' = Om*nO,  V' = (V^2-1)*nV/(k+1)
Centre form: Om = z^2 sig, V = z U.   Far form: N = 1/z.
"""
import sympy as sp
import json

exec(open('<LAB>/b5/system.py').read().split('# (e) KHA')[0])
S = k * (1 + V * z)**2 - (V + z)**2
Q = 2 * Om * V**2 * k * z + 2 * Om * V * k + 2 * Om * V + 2 * Om * z + V**2 * z - z
assert sp.simplify(Q - (-z * (1 - V**2) * (1 - C))) == 0
Dc = sp.factor(D.subs(mu, C))
assert sp.simplify(Dc - (-4 * Om * (k + 1) * S * Q**2)) == 0
NOf = sp.factor(sp.together(Omx_sol.subs(mu, C) * Dc))
NVf = sp.factor(sp.together(Vx_sol.subs(mu, C) * Dc))
nO = sp.cancel(NOf / (-4 * Om**2 * (k + 1) * Q))
nV = sp.cancel(NVf / (-4 * Om * (V**2 - 1) * Q))
P = sp.cancel(sp.together(zx_sol.subs(mu, C) / z) * Q)
for nm_, e in (('nO', nO), ('nV', nV), ('P', P)):
    assert sp.denom(sp.together(e)) == 1, nm_
nO, nV, P = [sp.expand(e) for e in (nO, nV, P)]
# verify
assert sp.simplify(Omx_sol.subs(mu, C) - Om * nO / (S * Q)) == 0
assert sp.simplify(Vx_sol.subs(mu, C) - (V**2 - 1) * nV / ((k + 1) * S * Q)) == 0
assert sp.simplify(zx_sol.subs(mu, C) - z * P / Q) == 0
print('nO =', nO)
print('nV =', nV)
print('P  =', P)

# sonic line: S=0 and numerators proportional?  check: on S=0 (branch V+z = sqrt(k)(1+Vz)), nO and nV share the condition
# centre form
sig, U = sp.symbols('sig U')
cs = {Om: z**2 * sig, V: z * U}
dsig = sp.cancel((Om * nO / (S * z * P)).subs(cs) / z**2 - 2 * sig / z)
dU = sp.cancel(((V**2 - 1) * nV / ((k + 1) * S * z * P)).subs(cs) / z - U / z)
print('dsig/dz =', sp.factor(dsig))
print('dU/dz   =', sp.factor(dU))
# far form, N = 1/z
Nn = sp.symbols('N')
dOm_dN = sp.cancel((-(Om * nO / (S * z * P)) * z**2).subs(z, 1 / Nn))
dV_dN = sp.cancel((-((V**2 - 1) * nV / ((k + 1) * S * z * P)) * z**2).subs(z, 1 / Nn))
print('dOm/dN =', sp.factor(dOm_dN))
print('dV/dN  =', sp.factor(dV_dN))
json.dump({'nO': str(nO), 'nV': str(nV), 'P': str(P), 'S': str(sp.expand(S)), 'mu': str(C),
           'dsig': str(dsig), 'dU': str(dU), 'dOm_dN': str(dOm_dN), 'dV_dN': str(dV_dN)},
          open('<LAB>/b5/system2.json', 'w'), indent=1)
