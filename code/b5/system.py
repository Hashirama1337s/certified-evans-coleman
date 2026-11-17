"""B5: build the reduced self-similar system from derive.py output and check it against the papers.

Writes system.json: polynomial right-hand sides (as sympy strings) for
  (i)  x-system:  mu_x, z_x, Om_x, V_x  (rational)
  (ii) desingularised field in (z, Om, V)  [dz/dtau, dOm/dtau, dV/dtau]  polynomial
  (iii) centre-regular form  Om = z^2 S, V = z U  in independent variable z
  (iv) far form  N = 1/z.
"""
import sympy as sp
import pickle
import json

raw = pickle.load(open('<LAB>/b5/derived_raw.pkl', 'rb'))
E = {k_: eval(v, dict(sp.__dict__)) for k_, v in raw.items()}

mu, z, Om, V, mux, zx, Omx, Vx = sp.symbols('mu z Om V mu_x z_x Om_x V_x', real=True)
k = sp.symbols('k', positive=True)
x = sp.symbols('x', real=True)
num = {kk: sp.numer(sp.together(v)) for kk, v in E.items()}
num = {kk: sp.expand(v.subs(sp.exp(x), 1).subs(sp.exp(-x), 1).subs(sp.exp(2 * x), 1)) for kk, v in num.items()}
for kk, v in num.items():
    assert not v.has(x), kk

# (a) Hamiltonian & momentum constraints
mux_H = sp.solve(num['tt'], mux)[0]
mux_M = sp.solve(num['tr'], mux)[0]
mu_C = sp.simplify(sp.solve(sp.Eq(mux_H, mux_M), mu)[0])
print('mu constraint (cf. Maison eq.(8)):', sp.factor(mu_C))
maison8 = 2 * Om * (1 + (1 + k) * V / z + k * V**2) / (1 - V**2)
assert sp.simplify(mu_C - maison8) == 0
print('  == Maison (8) with A^2=1-mu, B=z, rho~=Om : OK')

# (b) z_x from rr
zx_sol = sp.solve(num['rr'], zx)[0]
zx_sol = sp.simplify(zx_sol.subs(mux, mux_H))
print('z_x/z =', sp.factor(sp.simplify(zx_sol / z)))
Aa = 1 / (1 - mu)
maison4 = Aa * (2 * (1 - mu) - 1 + (1 - k) * Om)  # B'/B
assert sp.simplify(zx_sol / z - maison4) == 0
print('  == Maison (4) B\'/B : OK ; KHA (18b) with N=1/z, omega=Om*A, k=1/3:',
      sp.simplify((zx_sol / z).subs(k, sp.Rational(1, 3)) - (2 - Aa + sp.Rational(2, 3) * Om * Aa)) == 0)

# (c) fluid equations -> linear system in Om_x, V_x
f1 = sp.expand(num['divT_t'].subs({mux: mux_H, zx: zx_sol}))
f2 = sp.expand(num['divT_r'].subs({mux: mux_H, zx: zx_sol}))
Mm = sp.Matrix([[sp.diff(f, Omx), sp.diff(f, Vx)] for f in (f1, f2)])
bb = -sp.Matrix([f.subs({Omx: 0, Vx: 0}) for f in (f1, f2)])
sol = Mm.LUsolve(bb)
Omx_sol = sp.simplify(sol[0])
Vx_sol = sp.simplify(sol[1])
D = sp.factor(Mm.det())
print('det M =', D)

# (d) consistency: constraint is preserved: d/dx(mu - C) = 0 on mu = C
C = maison8
dC = sp.diff(C, Om) * Omx_sol + sp.diff(C, V) * Vx_sol + sp.diff(C, z) * zx_sol
chk = sp.simplify((dC - mux_H).subs(mu, C))
print('constraint propagation residual:', chk)
assert chk == 0

# (e) KHA eq. (18) check (self-similar: drop s-derivatives), k=1/3, N=1/z, A=1/(1-mu), omega = Om*A
kk = sp.Rational(1, 3)
A_ = 1 / (1 - mu)
N_ = 1 / z
om = Om * A_
subsk = {k: kk}
Ax = sp.diff(A_, mu) * mux_H
Nx = -zx_sol / z**2
omx = sp.diff(om, Om) * Omx_sol + sp.diff(om, mu) * mux_H
kha1 = (1 + N_ * V) * omx / om + 4 * (N_ + V) * Vx_sol / (3 * (1 - V**2)) - N_ * V * Ax / (3 * A_) + 4 * V * Nx / 3 \
    + 2 * N_ * (1 + 4 * om / (9 * (1 - V**2)))
kha2 = (4 * V + N_ + 3 * N_ * V**2) * omx / om + 4 * (1 + V**2 + 2 * N_ * V) * Vx_sol / (1 - V**2) \
    + N_ * (1 - V**2) * Ax / A_ + 4 * (1 + V**2) * Nx + 2 * N_ * (1 + 3 * V**2)
for nm_, e in (('KHA18c', kha1), ('KHA18d', kha2)):
    val = sp.simplify(e.subs(subsk).subs(mu, C.subs(subsk)))
    print(nm_, 'residual on our system:', val)

# (f) desingularised polynomial field in (z, Om, V): divide out positive factors
# Om_x = NO/D, V_x = NV/D, z_x = z*Zr.  Use mu = C.
Zr = sp.together((zx_sol / z).subs(mu, C))
NO = sp.together(Omx_sol.subs(mu, C) * D.subs(mu, C))
NV = sp.together(Vx_sol.subs(mu, C) * D.subs(mu, C))
Dc = sp.together(D.subs(mu, C))
print('Zr =', sp.factor(Zr))
print('D  =', sp.factor(Dc))
print('NO =', sp.factor(NO))
print('NV =', sp.factor(NV))
out = {'mu_C': str(C), 'Zr': str(sp.factor(Zr)), 'D': str(sp.factor(Dc)), 'NO': str(sp.factor(NO)),
       'NV': str(sp.factor(NV)), 'Omx': str(Omx_sol), 'Vx': str(Vx_sol), 'zx': str(zx_sol), 'muxH': str(mux_H)}
json.dump(out, open('<LAB>/b5/system.json', 'w'), indent=1)
