"""Collect the certified numbers (copied by code from the result JSONs) into results_summary.json and print them.
RESULT.md quotes this file verbatim."""
import json
import os
from flint import arb, ctx

ctx.prec = 400
B = '<LAB>/b5/'


def ld(n):
    p = B + n
    return json.load(open(p)) if os.path.exists(p) else None


ec, frw, nul = ld('result_ec.json'), ld('result_frw.json'), ld('result_null.json')
xc = ld('crosscheck50_ec.json')
S = {}


def lohi(s):
    x = arb(s)
    return {'lo': x.lower().str(45), 'hi': x.upper().str(45), 'width': arb(2 * x.rad()).str(3)}


if ec:
    S['EC'] = {k: ec.get(k) for k in ('miranda', 'sonic_type', 'V_zeros_inner', 'inner_checks', 'sonic_piece_ok',
                                       'outer_checks', 'resume_tested')}
    for k in ('sig0', 'zs', 'V_s', 'Om_s', 'mu_s', 'sonic_eig_pos', 'sonic_eig_neg_or_second', 't0_m_over_r',
              't0_2m_over_r', 't0_Om', 't0_V', 't0_a'):
        if k in ec:
            S['EC'][k] = lohi(ec[k])
if frw:
    S['FRW_k1/2'] = {k: frw.get(k) for k in ('miranda', 'sonic_type', 'V_zeros_inner', 'inner_checks', 'resume_tested')}
    for k in ('sig0', 'zs', 'V_s', 'Om_s', 'sonic_eig_pos', 'sonic_eig_neg_or_second'):
        if k in frw:
            S['FRW_k1/2'][k] = lohi(frw[k])
    ctx.prec = 400
    ex_sig = arb(8) / 27
    ex_zs = 9 * arb(2).sqrt() / 14
    ex_Vs = -arb(2).sqrt() * 2 / 5   # FRW: V_s = beta*y_s = -(4/9)(9 sqrt2/10)
    S['FRW_k1/2']['exact'] = {'sig0': ex_sig.str(40), 'zs': ex_zs.str(40), 'V_s': ex_Vs.str(40)}
    if 'sig0' in frw:
        S['FRW_k1/2']['exact_inside_box'] = {'sig0': bool(arb(frw['sig0']).contains(ex_sig)),
                                             'zs': bool(arb(frw['zs']).contains(ex_zs)),
                                             'V_s': bool(arb(frw['V_s']).contains(ex_Vs))}
if nul:
    S['nulls'] = nul
if xc:
    S['crosscheck50'] = xc
json.dump(S, open(B + 'results_summary.json', 'w'), indent=1)
print(json.dumps(S, indent=1))

# KHA-variable conversions at the sonic point (interval arithmetic): N = 1/z, A = 1/(1-mu), omega = Om/(1-mu)
if ec:
    ctx.prec = 400
    zs, mus, Oms = arb(ec['zs']), arb(ec['mu_s']), arb(ec['Om_s'])
    conv = {'N_s_KHA': 1 / zs, 'A_s_KHA': 1 / (1 - mus), 'omega_s_KHA': Oms / (1 - mus), 'a_s': 1 / (1 - mus).sqrt(),
            'eig_ratio_neg_over_pos': arb(ec['sonic_eig_neg_or_second']) / arb(ec['sonic_eig_pos'])}
    S['EC']['KHA_conversions'] = {k: lohi(v.str(radius=True)) for k, v in conv.items()}
    m = arb(ec['t0_m_over_r'])
    S['EC']['EC_comparison'] = {'EC_quoted_m_over_r': '0.0596 (gr-qc/9402041 p.4)', 'ours_minus_EC': (m - arb('0.0596')).str(6),
                                'EC_quoted_a': '1.07', 'EC_quoted_Omega': '9.56e-3'}
    json.dump(S, open(B + 'results_summary.json', 'w'), indent=1)
    print(json.dumps(S['EC']['KHA_conversions'], indent=1), S['EC']['EC_comparison'])
