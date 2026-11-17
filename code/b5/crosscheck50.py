"""50-digit NON-RIGOROUS cross-check with an independent code path (mpmath, Taylor-series stepping in z,
finite-difference eigenvectors + power-series parameterization at the sonic point, Frobenius series at the centre).
Shares only the derived polynomial system (system2.json) with the rigorous pipeline.
  1. residual: integrate the sonic branch (at the certified zs) inward to z_c and compare with the centre series at the
     certified sig0  -> should vanish to ~50 digits;
  2. t = 0 values: outer branch -> z = 3 (Taylor in z) -> N = 1/z from 1/3 to 0 (mpmath odefun) -> m/r.
usage: crosscheck50.py k_num k_den result.json"""
import sys
import json
import mpmath as mp
import numerics as nm
from series import load_funcs

mp.mp.dps = 60
kk = mp.mpf(int(sys.argv[1])) / int(sys.argv[2])
R = json.load(open(sys.argv[3]))
F = load_funcs()


def mid_of(s):
    s = s.strip('[]').split('+/-')[0]
    return mp.mpf(s)


sig0 = mid_of(R['sig0'])
zs = mid_of(R['zs'])
which = 'pos'
p, lam, evr = nm.param_branch(kk, zs, which, 80, mp.mpf('0.25'))
print('eigs', [mp.nstr(e, 12) for e in evr], 'lam', mp.nstr(lam, 15))


def ev(th):
    return [nm.evalc(p[j], th) for j in range(3)]


def th_for(zt, th0):
    return mp.findroot(lambda t: ev(t)[0] - zt, th0)


zx = zs - mp.mpf('0.01')
thx = th_for(zx, (zx - zs) / p[0][1])
Px = ev(thx)
print('inner sonic point used: theta', mp.nstr(thx, 10), 'z', mp.nstr(Px[0], 20))
zc = mp.mpf('0.02')
Om_c, V_c = nm.integrate(kk, zx, Px[1], Px[2], zc, N=60, frac=mp.mpf('0.2'))
s, u = nm.centre_series(kk, sig0, 60)
Om_s = zc**2 * nm.evalc(s, zc)
V_s = zc * nm.evalc(u, zc)
print('centre-series tail check |last terms|:', mp.nstr(abs(s[-1]) * zc**60, 3), mp.nstr(abs(u[-1]) * zc**60, 3))
res = [Om_c - Om_s, V_c - V_s]
print('residual at z=0.02 (sonic->centre Taylor vs centre series):', [mp.nstr(x, 5) for x in res])
out = {'residual_Om': mp.nstr(res[0], 10), 'residual_V': mp.nstr(res[1], 10)}
if len(sys.argv) > 4 and sys.argv[4] == 'outer':
    zo = zs + mp.mpf('0.01')
    tho = th_for(zo, (zo - zs) / p[0][1])
    Po = ev(tho)
    Om3, V3 = nm.integrate(kk, zo, Po[1], Po[2], mp.mpf(3), N=60, frac=mp.mpf('0.2'))
    # s = -N runs forward from -1/3 to 0 (mpmath odefun integrates forward only)
    f = mp.odefun(lambda s_, y: [-F['dOm_dN'](kk, -s_, y[0], y[1]), -F['dV_dN'](kk, -s_, y[0], y[1])],
                  -mp.mpf(1) / 3, [Om3, V3])
    Om0, V0 = f(mp.mpf(0))
    mu0 = 2 * Om0 * (1 + kk * V0**2) / (1 - V0**2)
    print('t=0: Om', mp.nstr(Om0, 50), '\n     V ', mp.nstr(V0, 50), '\n     m/r', mp.nstr(mu0 / 2, 50))
    out.update({'t0_Om': mp.nstr(Om0, 50), 't0_V': mp.nstr(V0, 50), 't0_m_over_r': mp.nstr(mu0 / 2, 50)})
json.dump(out, open(sys.argv[3].replace('result_', 'crosscheck50_'), 'w'), indent=1)
