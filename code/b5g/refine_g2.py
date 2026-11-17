"""Convergence evidence for kappa*: certified bisection to width ~1e-13 in configuration A or B.
usage: memcap 8 py -3.11 -u b5g/refine_g2.py A|B  -> refine_<cfg>.json"""
import json
import sys
import time
sys.path.insert(0, '<LAB>/b5g')
from flint import arb, ctx
import lohner_g2 as LG

which = sys.argv[1]
CF = {
    'A': dict(Nbg=30, Np=28, h0=0.02, frac=0.15, frac0=0.2, z1='0.05', zm='0.4947', prec=200),
    'B': dict(Nbg=36, Np=34, h0=0.016, frac=0.12, frac0=0.15, z1='0.04', zm='0.4935', prec=256),
}[which]
ctx.prec = CF['prec']
cfg = dict(sig0='[8.9013232777556435191172160567756856 +/- 5.31e-35]', zs='[0.49718669611129068678126112853988237484 +/- 2.61e-39]',
           Nbg=CF['Nbg'], Np=CF['Np'], h0=CF['h0'], frac=CF['frac'], frac0=CF['frac0'], z1=CF['z1'], zm=arb(CF['zm']), prec=CF['prec'])
t0 = time.time()
steps, info = LG.build_background2(cfg, log=lambda s: None)
P = LG.Pert2(steps, cfg)
a, b = arb('2.8105522155761715'), arb('2.8105525970458984')     # certified bracket from cert_A/cert_B
fa, fb = P.F(a)[0], P.F(b)[0]
assert fa < 0 and fb > 0
for it in range(80):
    if float((b - a).mid()) < 1e-13:
        break
    m = (a + b) / 2
    fm = P.F(m)[0]
    if fm < 0:
        a, fa = m, fm
    elif fm > 0:
        b, fb = m, fm
    else:
        print('sign of F(m) not resolved at width', float((b - a).mid()))
        break
K = a.union(b)
FK, dFK, lamK = P.F(K)
out = {'config': which, 'lo': a.str(20, radius=False), 'hi': b.str(20, radius=False), 'F_lo': fa.str(8, radius=True), 'F_hi': fb.str(8, radius=True),
       'dF_K': dFK.str(10, radius=True), 'certified': bool(fa < 0 and fb > 0 and dFK > 0), 'kappa_star': K.str(16, radius=True),
       'gamma': (1 / K).str(16, radius=True), 'time_s': round(time.time() - t0, 1)}
json.dump(out, open('<LAB>/b5g/refine_%s.json' % which, 'w'), indent=1)
print(json.dumps(out))
