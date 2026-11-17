"""diagnostics: widths of centre/sonic vectors; bisection of sign changes (float-only exploration)."""
import sys
import time
sys.path.insert(0, '<LAB>/b5g')
from flint import arb, ctx
import pipeline_g2 as PG
import lohner_g2 as LG
import frw_exact_g2

which = sys.argv[1]
mode = sys.argv[2]
Nbg = int(sys.argv[3]) if len(sys.argv) > 3 else 26
frac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.15
if which == 'frw':
    ctx.prec = 200
    cfg = dict(sig0=arb(3) / 8, zs=arb(3).sqrt() / 2, z1='0.05', zm=arb('0.8635'), Nbg=Nbg, Np=Nbg - 2, h0=0.02, frac=frac, frac0=0.2,
               prec=200, sonic_exact_bg=frw_exact_g2.frw_sonic_bg)
else:
    cfg = dict(sig0='[8.9013232777556435191172160567756856 +/- 5.31e-35]', zs='[0.49718669611129068678126112853988237484 +/- 2.61e-39]',
               z1='0.05', zm=arb('0.4947'), Nbg=Nbg, Np=Nbg - 2, h0=0.02, frac=frac, frac0=0.2, prec=200)
steps, info = LG.build_background2(cfg, log=lambda s: None)
print(info)
P = LG.Pert2(steps, cfg)
if mode == 'widths':
    for kv in sys.argv[5:]:
        Fv, dF, lam, allc, allk = P.F(arb(kv), want_vectors=True)
        print('kappa', kv, 'F', Fv.str(15, radius=True), 'lam', lam.str(8))
        for v in allc:
            print('   ', ['%.2e' % float(x.rad()) for x in v])
elif mode == 'scan':
    a, b, n = float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])
    prev = None
    for i in range(n + 1):
        kv = a + (b - a) * i / n
        Fv, dF, lam = P.F(arb(kv))
        print('%.5f  F=%s  dF=%s  lam=%s piv=%d' % (kv, Fv.str(10, radius=True), dF.str(6, radius=True), lam.str(6), P.last_piv))
elif mode == 'bisect':
    a, b = float(sys.argv[5]), float(sys.argv[6])
    fa = P.F(arb(a))[0]
    fb = P.F(arb(b))[0]
    assert (fa > 0 and fb < 0) or (fa < 0 and fb > 0), (fa, fb)
    for it in range(int(sys.argv[7]) if len(sys.argv) > 7 else 30):
        m = (a + b) / 2
        fm = P.F(arb(m))[0]
        if (fm > 0) == (fa > 0) and not fm.contains(0):
            a, fa = m, fm
        elif not fm.contains(0):
            b, fb = m, fm
        else:
            break
    print('root in [%.15f, %.15f]' % (a, b))
