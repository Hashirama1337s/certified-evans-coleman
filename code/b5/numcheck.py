"""check every [lo, hi] row of the RESULT.md tables against results_summary.json: printed lo <= certified lo, printed hi >= certified hi"""
import json, re
from flint import arb, ctx
ctx.prec = 400
S = json.load(open('results_summary.json'))['EC']
rows = {'z_s (Maison': S['zs'], "KHA's V at": S['V_s'], 'Ω_s = 4πr²ρ': S['Om_s'], '2m/r at the sonic': S['mu_s'],
        "KHA's N_s": S['KHA_conversions']['N_s_KHA'], "KHA's A_s": S['KHA_conversions']['A_s_KHA'],
        "KHA's ω_s": S['KHA_conversions']['omega_s_KHA'], 'σ0 = lim': S['sig0'],
        '| m/r |': S['t0_m_over_r'], '| Ω |': S['t0_Om'], '| V |': S['t0_V'], '| a |': S['t0_a']}
txt = open('RESULT.md', encoding='utf-8').read()
bad = 0
for line in txt.split('\n'):
    for key, ent in rows.items():
        if line.startswith('|') and key in line:
            cells = [c.strip().strip('*') for c in line.split('|')]
            lo, hi = arb(cells[2]), arb(cells[3])
            clo, chi = arb(ent['lo']), arb(ent['hi'])
            ok = (lo.upper() <= clo.lower()) and (hi.lower() >= chi.upper())
            print(key, 'OK' if ok else 'WRONG-DIRECTION')
            bad += (not ok)
# eigen lines
for key, ent in (('λ₊ ∈', S['sonic_eig_pos']), ('λ₋ ∈', S['sonic_eig_neg_or_second']), ('ratio λ₋/λ₊ ∈', S['KHA_conversions']['eig_ratio_neg_over_pos'])):
    m = re.search(re.escape(key) + r' \[([-0-9.e]+), ([-0-9.e]+)\]', txt)
    lo, hi = arb(m.group(1)), arb(m.group(2))
    ok = (lo.upper() <= arb(ent['lo']).lower()) and (hi.lower() >= arb(ent['hi']).upper())
    print(key, 'OK' if ok else 'WRONG-DIRECTION'); bad += (not ok)
print('bad', bad)
