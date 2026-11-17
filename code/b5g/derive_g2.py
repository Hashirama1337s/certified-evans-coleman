"""B5G derivation (own, from scratch): spherical Einstein-Euler in polar-areal slicing, rewritten in
self-similar variables (tau, x), tau = -ln(-t), x = ln(r/(-t)), t<0 = central proper time (alpha(t,0)=1 gauge).

Metric ds^2 = -al^2 dt^2 + A^2 dr^2 + r^2 dOmega^2, fluid p = k rho, u = (W/al, W V/A, 0, 0), W = 1/sqrt(1-V^2).
Unknowns (all functions of tau, x):  mu = 1 - 1/A^2 (= 2m/r),  Om = 4 pi r^2 rho,  V,  N = al (-t)/(A r).
Chain rule: (-t) d/dt = d_tau + d_x,  r d/dr = d_x.
Outputs b5g/eqs_g2.json with strings for:
  H  : mu_x = H(mu, Om, V)                     (G_tt)
  Mo : mu_tau = Mo(mu, Om, V, N)               (G_tr, with mu_x from H)
  SL : (ln N)_x = SL(mu, Om, V)                (G_rr, with mu_x from H)
  fluid: two equations linear in (Om_tau, V_tau, Om_x, V_x), with mu_tau, mu_x, N_x eliminated.
usage: memcap 8 py -3.11 -u b5g/derive_g2.py
"""
import json
import sympy as sp

t, r, th = sp.symbols('t r theta')
k = sp.symbols('k', positive=True)
al = sp.Function('al')(t, r)
A = sp.Function('A')(t, r)
rho = sp.Function('rho')(t, r)
Vf = sp.Function('V')(t, r)
X = [t, r, th, sp.Symbol('ph')]
g = sp.diag(-al**2, A**2, r**2, r**2 * sp.sin(th)**2)
gi = g.inv()
n = 4


def christ():
    G = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                G[a][b][c] = sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b])
                                                       - sp.diff(g[b, c], X[d])) for d in range(n)) / 2)
    return G


Gam = christ()


def ricci(b, c):
    e = 0
    for a in range(n):
        e += sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c])
        for d in range(n):
            e += Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a]
    return sp.simplify(e)


Ric = sp.Matrix(n, n, lambda b, c: ricci(b, c))
Rs = sp.simplify(sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n)))
Gt = sp.simplify(Ric - Rs * g / 2)
W = 1 / sp.sqrt(1 - Vf**2)
p = k * rho
uU = [W / al, W * Vf / A, 0, 0]
TU = sp.Matrix(n, n, lambda a, b: (rho + p) * uU[a] * uU[b] + p * gi[a, b])
TL = sp.simplify(g * TU * g)
E_tt = sp.simplify(Gt[0, 0] - 8 * sp.pi * TL[0, 0])
E_tr = sp.simplify(Gt[0, 1] - 8 * sp.pi * TL[0, 1])
E_rr = sp.simplify(Gt[1, 1] - 8 * sp.pi * TL[1, 1])
# fluid: nabla_mu T^{mu nu}
divT = []
for nu in range(2):
    e = 0
    for mu_ in range(n):
        e += sp.diff(TU[mu_, nu], X[mu_])
        for l in range(n):
            e += Gam[mu_][mu_][l] * TU[l, nu] + Gam[nu][mu_][l] * TU[mu_, l]
    divT.append(sp.simplify(e))

# ---- self-similar substitution ----
tau_s = sp.symbols('tau')
names = ['mu', 'Om', 'V', 'N']
F = {nm: sp.Function(nm + 'f')(t, r) for nm in names}
sym = {nm: sp.symbols(nm) for nm in names}
st = {nm: sp.symbols(nm + '_t') for nm in names}   # d_tau
sx = {nm: sp.symbols(nm + '_x') for nm in names}   # d_x
Aexpr = (1 - F['mu'])**sp.Rational(-1, 2)
sub_fields = {A: Aexpr, al: F['N'] * Aexpr * r / (-t), rho: F['Om'] / (4 * sp.pi * r**2), Vf: F['V']}


def to_ss(e):
    e = e.subs(sub_fields).doit()
    rep = {}
    for nm in names:
        rep[sp.Derivative(F[nm], t)] = (st[nm] + sx[nm]) / (-t)
        rep[sp.Derivative(F[nm], r)] = sx[nm] / r
    e = e.subs(rep)
    e = e.subs({F[nm]: sym[nm] for nm in names})
    return sp.simplify(e)


Ett = to_ss(E_tt * (-t)**2 / (sym['N']**2) )
Etr = to_ss(E_tr)
Err = to_ss(E_rr)
d0 = to_ss(divT[0])
d1 = to_ss(divT[1])
mu, Om, V, N = [sym[nm] for nm in names]
# H from G_tt
solH = sp.solve(Ett, sx['mu'])
assert len(solH) == 1
H = sp.simplify(solH[0])
print('H  mu_x =', H)
Etr2 = sp.simplify(Etr.subs(sx['mu'], H))
solM = sp.solve(Etr2, st['mu'])
assert len(solM) == 1
Mo = sp.simplify(solM[0])
print('Mo mu_tau =', Mo)
Err2 = sp.simplify(Err.subs(sx['mu'], H))
solS = sp.solve(Err2, sx['N'])
assert len(solS) == 1
SL = sp.simplify(solS[0] / N)
print('SL (lnN)_x =', SL)
fl = []
for d in (d0, d1):
    e = d.subs({sx['mu']: H, st['mu']: Mo, sx['N']: SL * N})
    e = sp.together(sp.simplify(e))
    num, den = sp.fraction(e)
    num = sp.factor(sp.expand(num))
    fl.append(num)
    print('fluid numerator (factored):', num)
    print('   denominator:', den)
out = {'H': str(H), 'Mo': str(Mo), 'SL': str(SL), 'fl0': str(fl[0]), 'fl1': str(fl[1])}
json.dump(out, open('<LAB>/b5g/eqs_g2.json', 'w'), indent=1)
print('written eqs_g2.json')
