"""Null controls (must give a CERTIFIED NONZERO matching residual):
  N1 off-solution data: sig0 -> sig0*(1+1e-10) at the certified zs (k=1/3)
  N2 wrong physics (perturbed constant): k = 17/50 instead of 1/3, at the k=1/3 solution's (sig0, zs)
  N3 FRW exact data (k=1/2 values) pushed through the k = 49/100 equations (k=51/100 makes that sonic point a focus, attempt 1)
Residual G = Gamma(sig0) - c(zs) at z = zm with rigorous enclosures; certified nonzero iff some component excludes 0.
usage: null.py"""
import json
import time
from flint import arb, ctx
import systems as SY
import sonic as SO
from rigode import guess, validate, integrate

T0 = time.time()
out = {}


def centre_G(kn, kd, sig, z1, zm, M=40):
    ctx.prec = 320
    kk = arb(kn) / kd
    U0 = -2 / (3 * (1 + kk))
    cs = SY.centre_system(False)
    rs = SY.regular_system(False)
    ub, dif = guess(cs, kk, arb(0), z1, [sig, U0], M, iters=300, tol=1e-90)
    res = validate(cs, kk, arb(0), z1, [sig, U0], ub)
    e = res['end']
    u0 = [z1**2 * e[0], z1 * e[1]]
    return integrate(rs, kk, z1, zm, u0, M=M, h0=0.02, gtol=1e-82, ctol=1e-84)


def sonic_G(kn, kd, zs, zm, which, th1, N=50):
    ctx.prec = 600
    kk = arb(kn) / kd
    B = SO.branch(kk, zs, N=N, which=which, scale=arb('0.2'))
    Vt = SO.validate_tail(B, arb(th1), pieces=16, kk=kk)
    th0 = float(((zm - zs) / B['p'][0][1]).mid())
    th, ok = SO.solve_theta(B, Vt, kk, zm, th0, abs(th0) * 0.3)
    assert ok
    return SO.P_eval(B, Vt, th)


def report(name, G, c):
    r = [G[0] - c[1], G[1] - c[2]]
    nz = any(not x.contains(0) for x in r)
    out[name] = {'residual': [x.str(8, radius=True) for x in r], 'certified_nonzero': nz}
    json.dump(out, open('<LAB>/b5/result_null.json', 'w'), indent=1)
    print('[%6.0fs] %s residual %s certified nonzero: %s' % (time.time() - T0, name, [x.str(6) for x in r], nz), flush=True)


ec = json.load(open('<LAB>/b5/hp_1_3.json'))
ctx.prec = 320
sig = arb(arb(ec['sig0']).mid())
zs = arb(arb(ec['zs']).mid())
zm = arb('0.4945')
# N1
c = sonic_G(1, 3, zs, zm, 'unstable', '0.05')
ctx.prec = 320
G = centre_G(1, 3, sig * (1 + arb('1e-10')), arb('0.05'), zm)
report('N1_offsolution_sig0_times_1+1e-10', G, c)
# N2
c = sonic_G(17, 50, zs, zm, 'unstable', '0.05')
ctx.prec = 320
G = centre_G(17, 50, sig, arb('0.05'), zm)
report('N2_k_17/50_at_k1/3_solution', G, c)
# N3 / P3: FRW closed form (k=1/2: V = beta y, Om = 3V^2/2, z = -y/(1+beta(1-beta)y^2)) at z = 0.6 compared with the
# validated centre-regular solution at sig0 = 8/27:  P3 (k = 1/2 equations) must CONTAIN it, N3 (k = 49/100) must EXCLUDE it.
ctx.prec = 320
sigF = arb(8) / 27
zt = arb('0.6')
be = 2 / (3 * (1 + arb(1) / 2))
cc = be * (1 - be)
y = (-1 + (1 - 4 * cc * zt**2).sqrt()) / (2 * cc * zt)
VF = be * y
OmF = 3 * VF**2 / 2
for name, kn, kd in (('P3_FRW_closed_form_vs_k1/2_solution_(must_contain)', 1, 2), ('N3_FRW_closed_form_vs_k49/100_solution', 49, 100)):
    G = centre_G(kn, kd, sigF, arb('0.05'), zt)
    report(name, G, [zt, OmF, VF])
json.dump(out, open('<LAB>/b5/result_null.json', 'w'), indent=1)
print('DONE')
