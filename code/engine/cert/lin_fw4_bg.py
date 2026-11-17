"""B4 proof seat v4 (new file).  Certified BACKGROUND COVER for the Far Wall certificate, produced with the
BAR-LB-E5 instrument unchanged (run_v32 setup, MatcherV32.residual on a box that contains the certified
E5 K-box).  Stores the same objects the certified PE2 step (pe2_v31) consumes:
  * sonic-manifold piece  theta in [th1, 0]: Y(theta), X(theta) series (val) + eps, epsX, nu_s;
  * tau-leg Lohner records  (x_mid, h, Bx = (N, om, V, x) a-priori box, start hull);
  * x-leg  Lohner records   (x_start, h, Bx = (N, om, V) a-priori box over [x_start + h, x_start]);
  * transport box Bt on [xbar, x_m]; centre manifold (Z = (P, W, U), chi) + eps, epsX, nu_c, thm;
  * margin = max(e_t, e0, e_end)  (PE2 inflation).
Nothing is changed in the sealed files; this only re-runs the residual at one box and pickles the records.
usage: lin_fw4_bg.py  -> lin_fw4_bg.pkl
"""
import sys, os, time, pickle
for f in ['LB_SMOKE', 'LB_EXPANDED', 'LB_MVSTART', 'LB_LADDER', 'LB_V3_LADDER', 'LB_NULADDER', 'LB_NU_AT',
          'LB_X1', 'LB_LOHNER_MODE', 'LB_V3ALIGN', 'LB_PREC', 'LB_KS', 'LB_KC', 'LB_KT', 'LB_NUC']:
    assert os.environ.get(f) is None, f
os.environ['LB_NGRID'] = '24'
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, ctx
ctx.prec = 200
import driver_v31 as V31
import driver_v32 as V32
import driver_v3 as V3
import driver_v2 as V2
import lohner as LO
import param as PM
from param import val
from core import sym
assert LO.lohner_step is V3._step and V3.MODE == 'qrw'
assert PM.arb_series is V31.arb_series_exact
T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)

PREC, KS, KC, KT, NUC = 200, 16, 60, 36, '0.03'
k13 = sp.Rational(1, 3)
om_t = '0.3256888894230684'
# certified E5 K-box (runs6/target.log, PE1 rung 1e-15)
KV0 = arb('0.1124394013880983362') + sym(arb('3.31e-20') + arb('1e-19'))
KW0 = arb('8.90132327775564352') + sym(arb('7.92e-18') + arb('1e-17'))
Kxc = arb('0.212365646758945385') + sym(arb('7.95e-19') + arb('1e-18'))
R = arb('1e-15')
Xb = [arb(KV0.mid()) + sym(R), arb(KW0.mid()) + sym(R), arb(Kxc.mid()) + sym(R)]
assert all(Xb[i].contains(K_) for i, K_ in enumerate([KV0, KW0, Kxc])), 'box must contain the E5 K-box'

Pb = V31.Problem2(k13, prec=PREC, K_sonic=KS, K_centre=KC, K_ty=KT)
nu, _ = V2.choose_nu(Pb, arb('0.112439401388098336223763027603135130845') + sym(arb('1e-9')), om_t,
                     ['0.004', '0.002', '0.001'], KS, log=log)
M = V32.MatcherV32(Pb, nu, NUC, Ks=KS, x_switch='-0.05', Kt=KT)
log('nu_s =', nu)
Rr, E, d = M.residual(Xb[0], Xb[1], Xb[2], om_t, log=log)
margin = d['et_max'].max(d['e0']).max(d['e_end'])
log('residual done: margin =', margin.str(3), ' Erad =', d['Erad'].str(3))

def ser(a):
    a = val(a)
    return (a.mid().str(70, radius=False), a.rad().str(10, radius=False))
def ser_rec(rec):
    return [(float(x), ser(h), [ser(b) for b in Bx], [ser(s) for s in st]) for (x, h, Bx, st) in rec]
sm, cm = d['sm'], d['cm']
out = dict(
    note='E5 instrument background cover (lin_fw4_bg.py); boxes valid for every background with (V0,W0,xc) in Xb',
    Xb=[ser(t) for t in Xb], margin=ser(margin), nu_s=nu, nu_c=NUC, th1=ser(M.th1), thm=ser(M.thm),
    x_switch=ser(M.x_switch), xbar=ser(M.xbar),
    sm=dict(Y=[[ser(a) for a in sm['co']['Y'][i]] for i in range(3)], X=[ser(a) for a in sm['co']['X']],
            eps=ser(sm['eps']), epsX=ser(sm['epsX']), lam=ser(sm['lam'])),
    cm=dict(Y=[[ser(a) for a in cm['co']['Y'][i]] for i in range(3)], X=[ser(a) for a in cm['co']['X']],
            eps=ser(cm['eps']), epsX=ser(cm['epsX']), lam=ser(cm['lam'])),
    tau_rec=ser_rec(d['tau_rec']), rec=ser_rec(d['rec']), Bt=[ser(b) for b in d['Bt']],
    e_t=[ser(t) for t in d['e_t']], e0=ser(d['e0']), e_end=ser(d['e_end']))
pickle.dump(out, open('<LAB>/engine/cert/lin_fw4_bg.pkl', 'wb'))
log('saved lin_fw4_bg.pkl: tau steps %d, x steps %d' % (len(out['tau_rec']), len(out['rec'])))
