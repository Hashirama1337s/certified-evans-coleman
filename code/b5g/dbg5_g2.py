import sys
sys.path.insert(0, '<LAB>/b5g')
from flint import arb, ctx
import lohner_g2 as LG, frw_exact_g2
ctx.prec = 200
cfg = dict(sig0=arb(3) / 8, zs=arb(3).sqrt() / 2, z1='0.05', zm=arb('0.8635'), Nbg=30, Np=28, h0=0.02, frac=0.15, frac0=0.2, prec=200, sonic_exact_bg=frw_exact_g2.frw_sonic_bg)
steps, info = LG.build_background2(cfg, log=lambda s: None)
P = LG.Pert2(steps, cfg)
for r in (1e-4, 1e-5, 1e-6):
    print(r, P.F(arb(1, r))[1].str(8, radius=True))
