"""B5 blind seat: derive the self-similar Einstein-Euler ODEs from scratch (sympy).

Metric (Evans-Coleman gr-qc/9402041 eq. (1); KHA gr-qc/9503007 eq. (16)):
    ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2
Fluid T^{mu nu} = (rho+p) u^mu u^nu + p g^{mu nu}, p = k rho,
u^t = W/alpha, u^r = W V/a, W = (1-V^2)^(-1/2)  (V = Eulerian 3-velocity, KHA sec. 3).
Equations: G_{mu nu} = 8 pi T_{mu nu}, nabla_mu T^{mu nu} = 0.

Self-similar ansatz, t<0:  x = ln r - ln(-t),  xi = e^x,
  a = a(x), alpha = alpha(x), V = V(x), rho = Om(x)/(4 pi r^2).
Reduced variables: mu = 2m/r = 1 - 1/a^2 ;  z = a xi/alpha.
"""
import sympy as sp
import pickle

t = sp.symbols('t', negative=True)
r = sp.symbols('r', positive=True)
th, ph = sp.symbols('theta phi', real=True)
k = sp.symbols('k', positive=True)
x = sp.symbols('x', real=True)
X = [t, r, th, ph]
xe = sp.log(r) - sp.log(-t)
aS, alS, VS, OmS = [sp.Function(nm) for nm in ('aS', 'alS', 'VS', 'OmS')]
a = aS(xe)
al = alS(xe)
V = VS(xe)
rho = OmS(xe) / (4 * sp.pi * r**2)

g = sp.diag(-al**2, a**2, r**2, r**2 * sp.sin(th)**2)
gi = sp.diag(-1 / al**2, 1 / a**2, 1 / r**2, 1 / (r**2 * sp.sin(th)**2))
n = 4
Gam = [[[sum(gi[i, l] * (sp.diff(g[l, j], X[kk]) + sp.diff(g[l, kk], X[j]) - sp.diff(g[j, kk], X[l]))
             for l in range(n)) / 2 for kk in range(n)] for j in range(n)] for i in range(n)]


def ricci(i, j):
    res = 0
    for l in range(n):
        res += sp.diff(Gam[l][i][j], X[l]) - sp.diff(Gam[l][i][l], X[j])
        for m in range(n):
            res += Gam[l][l][m] * Gam[m][i][j] - Gam[l][j][m] * Gam[m][i][l]
    return res


Ric = sp.Matrix(n, n, lambda i, j: ricci(i, j))
R = sum(gi[i, i] * Ric[i, i] for i in range(n))
W = 1 / sp.sqrt(1 - V**2)
uU = [W / al, W * V / a, 0, 0]
p = k * rho
TU = sp.Matrix(n, n, lambda i, j: (rho + p) * uU[i] * uU[j] + p * gi[i, j])
TD = g * TU * g


def EIN(i, j):
    return Ric[i, j] - R * g[i, j] / 2 - 8 * sp.pi * TD[i, j]


divT = []
for nu in range(2):
    e = 0
    for m_ in range(n):
        e += sp.diff(TU[m_, nu], X[m_])
        for l in range(n):
            e += Gam[m_][m_][l] * TU[l, nu] + Gam[nu][m_][l] * TU[m_, l]
    divT.append(e)


def fix(e):
    e = e.doit()
    e = e.replace(lambda q: isinstance(q, sp.Subs), lambda q: q.expr.subs(q.variables[0], x))
    e = e.subs(xe, x)
    e = e.subs(t, -r * sp.exp(-x))
    return e


raw = {'tt': EIN(0, 0), 'tr': EIN(0, 1), 'rr': EIN(1, 1), 'thth': EIN(2, 2),
       'divT_t': divT[0], 'divT_r': divT[1]}
mu, z, Om, Vv, mup, zp, Omp, Vp = sp.symbols('mu z Om V mu_x z_x Om_x V_x', real=True)
aexp = 1 / sp.sqrt(1 - mu)
alexp = aexp * sp.exp(x) / z
# d/dx of a and alpha in terms of mu_x, z_x
ax = sp.diff(aexp, mu) * mup
alx = sp.diff(alexp, mu) * mup + sp.diff(alexp, z) * zp + alexp
axx_bad = sp.Symbol('SECOND_DERIV')
fin = {}
for nm, e in raw.items():
    e = fix(e)
    A2, L2 = sp.symbols('A_xx AL_xx')
    e = e.subs({sp.Derivative(aS(x), (x, 2)): A2, sp.Derivative(alS(x), (x, 2)): L2})
    if nm == 'thth':
        print(nm, 'second-order equation (implied by the others via Bianchi) -> not used')
        continue
    e = e.subs({sp.Derivative(aS(x), x): ax, sp.Derivative(alS(x), x): alx,
                sp.Derivative(VS(x), x): Vp, sp.Derivative(OmS(x), x): Omp})
    e = e.subs({aS(x): aexp, alS(x): alexp, VS(x): Vv, OmS(x): Om})
    e = sp.simplify(e)
    for pw in range(0, 6):
        e2 = sp.simplify(e * r**pw)
        if not e2.has(r):
            break
    e2 = sp.simplify(e2)
    assert not e2.has(r), (nm, e2)
    assert not e2.has(A2) and not e2.has(L2), ('second derivatives survive', nm)
    fin[nm] = e2
    print(nm, ':', e2)

with open('<LAB>/b5/derived_raw.pkl', 'wb') as fh:
    pickle.dump({kk: sp.srepr(v) for kk, v in fin.items()}, fh)
