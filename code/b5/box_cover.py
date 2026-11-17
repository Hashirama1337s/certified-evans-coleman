import json
from flint import arb, ctx
from cheb import aup
for cf in ('cfg_ec.json', 'cfg_frw.json'):
    cfg = json.load(open(cf)); hp = json.load(open(cfg['hp']))
    for prec in (320, 600):
        ctx.prec = prec
        sig_c = arb(arb(hp['sig0']).mid()); wS = arb(cfg['w_sig'])
        sig_lo, sig_hi = arb((sig_c - wS).mid()), arb((sig_c + wS).mid())
        sig_box = sig_c + arb(0, aup(wS)); mv = arb(0, aup(wS))
        print(cf, prec, 'box covers faces:', sig_box.lower() <= sig_lo and sig_box.upper() >= sig_hi,
              ' mean-value term covers |sig-sig_c|:', mv.upper() >= (sig_hi - sig_c).upper() and mv.upper() >= (sig_c - sig_lo).upper())
