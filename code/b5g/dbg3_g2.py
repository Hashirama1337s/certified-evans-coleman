import sys, time
sys.path.insert(0, '<LAB>/b5g')
from flint import arb, ctx
import lohner_g2 as LG
cfg = dict(sig0='[8.9013232777556435191172160567756856 +/- 5.31e-35]', zs='[0.49718669611129068678126112853988237484 +/- 2.61e-39]',
           z1='0.05', zm=arb('0.4947'), Nbg=30, Np=28, h0=0.02, frac=0.15, frac0=0.2, prec=200)
steps, info = LG.build_background2(cfg, log=lambda s: None)
P = LG.Pert2(steps, cfg)
for (m, r) in [(0.5, 0.025), (1.0, 0.05), (2.0, 0.05), (2.6, 0.02), (1.0, 0.0)]:
    t0 = time.time()
    K = arb(m, r)
    Fk, dFk, lam = P.F(K)
    Fc, dFc, _ = P.F(arb(m))
    print(m, r, 'F(K)', Fk.str(5, radius=True), 'dF(K)', dFk.str(5, radius=True), 'F(c)', Fc.str(8, radius=True), 'lam(K)', lam.str(5, radius=True), '%.1fs' % (time.time() - t0))
