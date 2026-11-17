"""Centre analysis (symbolic). Background: Om = z^2 sig, V = z U, N = 1/z, mu algebraic.
Perturbation scaled: dmu = z^2 m, dOm = z^2 p, dV = z q, eta = z^2 e.  z Y' = M(z) Y.
Prints M(0) and its eigenvalues as functions of kappa; writes centre_sys.json (cancelled rational functions).
"""
import json
import sympy as sp
from system_g2 import FOmx, FVx, MUalg, SL, pert_matrices, mu, Om, V, N, kap, H

z, sig, U = sp.symbols('z sig U')
sub = {Om: z**2 * sig, V: z * U, N: 1 / z}
MUc = sp.factor(MUalg.subs(sub))
print('mu =', MUc)
zx = sp.factor((-z * SL).subs(mu, MUalg).subs(sub))
print('z_x =', zx)
Omx = FOmx.subs(mu, MUalg).subs(sub)
Vx = FVx.subs(mu, MUalg).subs(sub)
dOm_dz = Omx / zx
dV_dz = Vx / zx
dsig = sp.factor(sp.cancel((dOm_dz - 2 * z * sig) / z**2))
gU = sp.factor(sp.cancel(z * (dV_dz - U) / z))   # z dU/dz = z*(V' - U)/z
print('dsig/dz =', dsig)
print('z dU/dz =', gU)
print('g(0) =', sp.factor(gU.subs(z, 0)))
P0, P1, det, adj, Cu, Cmu, CN, Ct = pert_matrices()
Pk = P0 + kap * P1
bgs = {mu: MUalg}
scale = sp.diag(z**2, z**2, z, z**2)  # Y = scale * Ys
# z Ys' = z scale^-1 (P/zx) scale Ys - z scale^-1 scale' Ys
dscale = sp.diag(2 * z, 2 * z, 1, 2 * z)
M = sp.zeros(4, 4)
for i in range(4):
    for j in range(4):
        e = Pk[i, j].subs(bgs).subs(sub)
        e = z * e / zx * scale[j, j] / scale[i, i]
        if i == j:
            e -= z * dscale[i, i] / scale[i, i]
        M[i, j] = sp.factor(sp.cancel(sp.together(e)))
M0 = M.subs(z, 0).subs(U, sp.Rational(-1, 2))
M0 = sp.simplify(M0)
print('M(0) =', M0)
ev = M0.eigenvals()
print('eigenvalues of M(0):', ev)
json.dump({'dsig': str(dsig), 'gU': str(gU), 'mu': str(MUc), 'zx': str(zx),
           'M': [[str(M[i, j]) for j in range(4)] for i in range(4)]},
          open('<LAB>/b5g/centre_sys.json', 'w'), indent=1)
print('written centre_sys.json')
