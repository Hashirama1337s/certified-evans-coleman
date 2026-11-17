"""UNTESTED curiosity check: sonic exponent of the perturbation lambda(kappa) vs background sonic eigenvalue ratio."""
import sys
sys.path.insert(0, '<LAB>/b5g')
from flint import arb, ctx
import lohner_g2 as LG
ctx.prec = 200
cfg = dict(sig0='[8.9013232777556435191172160567756856 +/- 5.31e-35]', zs='[0.49718669611129068678126112853988237484 +/- 2.61e-39]',
           Nbg=30, Np=28, h0=0.02, frac=0.15, frac0=0.2, z1='0.05', zm=arb('0.4947'), prec=200)
steps, info = LG.build_background2(cfg, log=lambda s: None)
P = LG.Pert2(steps, cfg)
l0 = P.sonic_data(arb(0))[2]
l1 = P.sonic_data(arb(1))[2]
print('lambda(0) =', l0.str(20, radius=True), ' lambda(1) =', l1.str(20, radius=True), ' slope =', (l1 - l0).str(20, radius=True))
lp = arb('[0.164293492237524824696953885384 +/- 4.15e-32]')   # b5/result_ec.json sonic_eig_pos
lm = arb('[-0.0148038245662013370930339955438 +/- 3.29e-32]')  # b5/result_ec.json sonic_eig_neg_or_second
print('lambda_-/lambda_+ (b5) =', (lm / lp).str(20, radius=True))
