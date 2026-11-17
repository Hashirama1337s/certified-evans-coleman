"""Newton re-centring of the Miranda box from the certified attempt-2 data (G_mid from the centre cache, c(zs) thin).
x* ~ x_c - C (G_mid - c(z_c)),  C = numerical inverse Jacobian (non-rigorous; only moves the box centre)."""
import json
from flint import arb, ctx
ctx.prec = 600
import sonic as SO
cfg = json.load(open('cfg_ec.json'))
hp = json.load(open(cfg['hp']))
cache = json.load(open('result_ec_centre_cache.json'))
kk = arb(1) / 3
sig_c, zs_c = arb(arb(hp['sig0']).mid()), arb(arb(hp['zs']).mid())
zm = arb(cfg['zm'])
G_mid = [arb(x) for x in cache['G_mid']]
G_box = [arb(x) for x in cache['G_box']]


def c_at(z):
    B0 = SO.branch(kk, z, N=50, which='unstable', scale=arb('0.2'))
    Vt = SO.validate_tail(B0, arb('0.05'), pieces=16, kk=kk)
    th0 = float(((zm - z) / B0['p'][0][1]).mid())
    th, ok = SO.solve_theta(B0, Vt, kk, zm, th0, abs(th0) * 0.3)
    return SO.P_eval(B0, Vt, th)


h = arb('1e-30')
c0 = c_at(zs_c)
c1 = c_at(zs_c + h)
dc = [(c1[1] - c0[1]) / h, (c1[2] - c0[2]) / h]
J = [[G_box[2], -dc[0]], [G_box[3], -dc[1]]]
det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
C = [[J[1][1] / det, -J[0][1] / det], [-J[1][0] / det, J[0][0] / det]]
g = [G_mid[0] - c0[1], G_mid[1] - c0[2]]
H = [C[0][0] * g[0] + C[0][1] * g[1], C[1][0] * g[0] + C[1][1] * g[1]]
print('H at centre', [x.str(8) for x in H])
new = {'k': '1/3', 'sig0': arb((sig_c - H[0]).mid()).str(70), 'zs': arb((zs_c - H[1]).mid()).str(70),
       'note': 'Newton re-centred from attempt 2 (refine.py)'}
json.dump(new, open('hp_1_3_refined.json', 'w'), indent=1)
print(new)
