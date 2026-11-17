"""Development checks: codegen vs sympy; FRW background vs closed form; F(kappa) samples."""
import sys
import time
sys.path.insert(0, '<LAB>/b5g')
import sympy as sp
from flint import arb, ctx
import pipeline_g2 as PG
import lohner_g2 as LG

ctx.prec = 200
# 1. codegen positive control
vals = {PG.mu_: sp.Rational(1, 7), PG.Om_: sp.Rational(2, 9), PG.V_: sp.Rational(-1, 5), PG.N_: sp.Rational(5, 2),
        PG.Omx_: sp.Rational(3, 11), PG.Vx_: sp.Rational(-2, 13)}
got = PG.f_fl(*[arb(sp.Rational(v).p) / sp.Rational(v).q for v in vals.values()])
for nm, g in zip(PG.FLN, got):
    ex = PG.S(nm).subs(vals)
    assert g.contains(arb(ex.p) / ex.q), nm
print('codegen fl OK')
zz, ss, uu = sp.Rational(1, 10), sp.Rational(9, 1), sp.Rational(-1, 2)
ga = PG.f_cMa(arb(1) / 10, arb(9), arb(-1) / 2)
for idx, e in enumerate(PG.Mca):
    ev = sp.nsimplify(e.subs({PG.z_: zz, PG.sig_: ss, PG.U_: uu}))
    assert ga[idx].contains(arb(ev.p) / ev.q), idx
print('codegen centre M OK')

which = sys.argv[1] if len(sys.argv) > 1 else 'frw'
if which == 'frw':
    cfg = dict(sig0=arb(3) / 8, zs=arb(3).sqrt() / 2, z1='0.05', zm=None, Nbg=26, Np=24, h0=0.02, frac=0.15, frac0=0.2, prec=200, verbose=False)
    cfg['zm'] = arb('0.8635')
    import frw_exact_g2
    cfg['sonic_exact_bg'] = frw_exact_g2.frw_sonic_bg
else:
    cfg = dict(sig0='[8.9013232777556435191172160567756856 +/- 5.31e-35]', zs='[0.49718669611129068678126112853988237484 +/- 2.61e-39]',
               z1='0.05', zm=arb('0.4947'), Nbg=26, Np=24, h0=0.02, frac=0.15, frac0=0.2, prec=200, verbose=False)
t0 = time.time()
steps, info = LG.build_background2(cfg)
print(info)
if which == 'frw':
    for st in steps[1:-1:10]:
        z = st['z0']
        Vf = -(1 - (1 - z * z).sqrt()) / z
        print('z=%s  V=%s exact=%s ' % (z.str(6), st['y'][1].ev(arb(0)).str(15, radius=True), Vf.str(15, radius=True)))
P = LG.Pert2(steps, cfg)
print('prep %.1fs' % P.tprep)
for kv in sys.argv[2:] if len(sys.argv) > 2 else ['0.5', '1', '1.5']:
    t1 = time.time()
    Fv, dF, lam = P.F(arb(kv))
    print('kappa=%s  F=%s  dF=%s  lambda_sonic=%s  (%.1fs)' % (kv, Fv.str(12, radius=True), dF.str(8, radius=True), lam.str(10, radius=True), time.time() - t1))
