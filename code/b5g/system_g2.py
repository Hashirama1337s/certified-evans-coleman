"""Build background + linear-perturbation systems from eqs_g2.json (own derivation, derive_g2.py).
Perturbation unknowns Y = (dmu, dOm, dV, eta), eta = dN/N, all ~ exp(kappa tau) f(x).
x-system:  Y_x = P(x; kappa) Y  with  P = P0 + kappa P1   (fluid block solved with the 2x2 matrix Ax).
Everything is returned as sympy expressions in background symbols (mu, Om, V, N) and kappa.
"""
import json
import sympy as sp

d = json.load(open('<LAB>/b5g/eqs_g2.json'))
k = sp.Rational(1, 3)
mu, Om, V, N = sp.symbols('mu Om V N')
Om_t, V_t, Om_x, V_x = sp.symbols('Om_t V_t Om_x V_x')
kap = sp.symbols('kappa')
loc = dict(mu=mu, Om=Om, V=V, N=N, Om_t=Om_t, V_t=V_t, Om_x=Om_x, V_x=V_x, k=sp.Symbol('k', positive=True))
ks = loc['k']
H = sp.sympify(d['H'], locals=loc).subs(ks, k)
Mo = sp.sympify(d['Mo'], locals=loc).subs(ks, k)
SL = sp.sympify(d['SL'], locals=loc).subs(ks, k)
t_ = sp.Symbol('t')
fl = []
for nm in ('fl0', 'fl1'):
    e = sp.sympify(d[nm], locals={**loc, 't': t_}).subs(ks, k)
    e = sp.factor(e)
    # strip factors that never vanish in the domain (t, V-1, V+1, constants)
    e = sp.Mul(*[f for f in sp.Mul.make_args(e) if not (f.free_symbols <= {t_} or f in (V - 1, V + 1, -(V - 1)) or f.is_number)])
    fl.append(sp.expand(e))

# stationary background: solve fluid for Om_x, V_x
bg = sp.solve([f.subs({Om_t: 0, V_t: 0}) for f in fl], [Om_x, V_x], dict=True)
assert len(bg) == 1
FOmx = sp.factor(bg[0][Om_x])
FVx = sp.factor(bg[0][V_x])
MUalg = sp.solve(sp.Eq(Mo, 0), mu)[0]          # stationary momentum constraint -> algebraic mu
zx_over_z = -SL                                  # z = 1/N
Ax = sp.Matrix(2, 2, lambda i, j: sp.diff(fl[i], [Om_x, V_x][j]))
At = sp.Matrix(2, 2, lambda i, j: sp.diff(fl[i], [Om_t, V_t][j]))
detAx = sp.factor(Ax.det())


def pert_matrices():
    """P0, P1 (4x4) in terms of (mu, Om, V, N) with Om_x, V_x replaced by the background values."""
    bgsub = {Om_t: 0, V_t: 0}
    # fluid linearisation: Ax dU_x + kappa At dU + Bu dU + Bmu dmu + BN N eta = 0, dU = (dOm, dV)
    Bu = sp.Matrix(2, 2, lambda i, j: sp.diff(fl[i], [Om, V][j]))
    Bmu = sp.Matrix(2, 1, lambda i, j: sp.diff(fl[i], mu))
    BN = sp.Matrix(2, 1, lambda i, j: sp.diff(fl[i], N) * N)
    rep = {Om_x: FOmx, V_x: FVx, Om_t: 0, V_t: 0}
    Axb = Ax.subs(rep)
    adj = Axb.adjugate()
    det = sp.factor(Axb.det())
    P0 = sp.zeros(4, 4)
    P1 = sp.zeros(4, 4)
    # rows: 0 dmu, 1 dOm, 2 dV, 3 eta;  dmu_x = dH ; eta_x = dSL
    gradH = [sp.diff(H, s) for s in (mu, Om, V)]
    gradS = [sp.diff(SL, s) for s in (mu, Om, V)]
    for j in range(3):
        P0[0, j] = gradH[j]
        P0[3, j] = gradS[j]
    # fluid
    Cu = Bu.subs(rep)
    Cmu = Bmu.subs(rep)
    CN = BN.subs(rep)
    Ct = At.subs(rep)
    for i in range(2):
        for j in range(2):
            P0[1 + i, 1 + j] = -(adj * Cu)[i, j] / det
            P1[1 + i, 1 + j] = -(adj * Ct)[i, j] / det
        P0[1 + i, 0] = -(adj * Cmu)[i, 0] / det
        P0[1 + i, 3] = -(adj * CN)[i, 0] / det
    return P0, P1, det, adj, Cu, Cmu, CN, Ct


if __name__ == '__main__':
    print('Om_x =', FOmx)
    print('V_x =', FVx)
    print('mu (stationary) =', sp.factor(MUalg))
    print('det Ax =', detAx)
    print('det Ax at bg =', sp.factor(Ax.subs({Om_x: FOmx, V_x: FVx}).det()))
    print('At =', At)
