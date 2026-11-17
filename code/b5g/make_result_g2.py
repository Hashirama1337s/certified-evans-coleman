"""Assemble b5g/result_g2.json from cert_A.json (main) and cert_B.json (convergence)."""
import json
from flint import arb, ctx

ctx.prec = 200
A = json.load(open('<LAB>/b5g/cert_A.json'))
B = json.load(open('<LAB>/b5g/cert_B.json'))
kA, kB = A['kappa_star'], B['kappa_star']
lo = max(float(kA['lo']), float(kB['lo']))
hi = min(float(kA['hi']), float(kB['hi']))
overlap = lo < hi
gA = arb(A['gamma'])
gB = arb(B['gamma'])
res = {
    'seat': 'b5g (blind G2 second implementation)',
    'formulation': 'polar-areal slicing, central-lapse gauge alpha(t,0)=1; variables (dmu, dOm, dV, eta=dN/N) ~ exp(kappa tau), '
                   'tau=-ln(-t), kappa>0 growing; independent variable z=1/N; interval Taylor models + causal tail bootstrap + '
                   'Lohner QR; analytic Frobenius solutions at centre and sonic point; F = normalised 4x4 determinant',
    'item1_kappa_star': {'enclosure': kA['enclosure'], 'lo': kA['lo'], 'hi': kA['hi'], 'width': kA['width'],
                         'F_lo': kA['F_lo'], 'F_hi': kA['F_hi'], 'dF_on_interval': kA['dF_K'],
                         'sonic_exponent_on_interval': kA['lam_K'], 'certified': kA['certified']},
    'item1_enclosure_endpoints': '[%s, %s]' % (kA['lo'], kA['hi']),
    'item1_gamma': A['gamma'],
    'item1_gamma_endpoints': '[%s, %s]' % ((1 / arb(kA['hi'])).str(12, radius=False), (1 / arb(kA['lo'])).str(12, radius=False)),
    'item1_gamma_note': 'gamma endpoints = 1/hi, 1/lo printed to 12 digits (round-to-nearest); the rigorous ball is item1_gamma',
    'item1_convergence_B': {'enclosure': kB['enclosure'], 'lo': kB['lo'], 'hi': kB['hi'], 'gamma': B['gamma'],
                            'certified': kB['certified'], 'overlap_with_A': overlap,
                            'common_interval': '[%r, %r]' % (lo, hi) if overlap else None},
    'item1_certified': bool(kA['certified'] and float(kA['width']) <= 1e-6),
    'item2_frw_control': {'closed_form_rate': '1 (time-translation mode t -> t + const; derived by hand)',
                          'A': A['frw_control'], 'B': B['frw_control'],
                          'eigenfunction_matches_closed_form_A': A['frw_gauge_eigenfunction_parallel'],
                          'eigenfunction_matches_closed_form_B': B['frw_gauge_eigenfunction_parallel'],
                          'near_kappa0_pole_not_zero': A['frw_near_kappa0']},
    'item2_certified': bool(A['frw_control']['contains_closed_form'] and B['frw_control']['contains_closed_form']),
    'item3_zeros_undeflated_0p40_2p70': [
        {'enclosure': A['ec_gauge_zero']['enclosure'], 'exact_value': '1', 'F_lo': A['ec_gauge_zero']['F_lo'],
         'F_hi': A['ec_gauge_zero']['F_hi'], 'dF_on_interval': A['ec_gauge_zero']['dF_K'],
         'identification': 'GAUGE: time translation t -> t + const (the only residual coordinate freedom of polar-areal '
                           'slicing with alpha(t,0)=1); eigenfunction e^tau (mu_x, Om_x, V_x, (lnN)_x + 1)',
         'eigenfunction_check': A['ec_gauge_eigenfunction_parallel']}],
    'item3_certified': bool(A['ec_gauge_zero']['certified_unique_simple_zero'] and
                            A['deflated_sign_0p40_2p70']['certified_constant_sign']),
    'item4_deflated': {'function': 'G(kappa) = F(kappa)/(kappa - 1)', **A['deflated_sign_0p40_2p70'],
                       'kappa_star_in_range': False},
    'item4_certified': bool(A['deflated_sign_0p40_2p70']['certified_constant_sign']),
    'controls': {'ec_background_overlap_at_zm_A': A['ec_background']['overlap'],
                 'ec_background_overlap_at_zm_B': B['ec_background']['overlap'],
                 'ec_momentum_constraint_A': {k: v for k, v in A.items() if k.startswith('ec_constraint')},
                 'frw_momentum_constraint_A': A['frw_constraint_kappa_0.3'],
                 'ec_gauge_eigenfunction_A': A['ec_gauge_eigenfunction_parallel'],
                 'ec_gauge_eigenfunction_B': B['ec_gauge_eigenfunction_parallel']},
    'inputs_reused_from_b5': {'sig0': '[8.9013232777556435191172160567756856 +/- 5.31e-35]',
                              'zs': '[0.49718669611129068678126112853988237484 +/- 2.61e-39]'},
    'item1_refined_convergence': {'A': json.load(open('<LAB>/b5g/refine_A.json')), 'B': json.load(open('<LAB>/b5g/refine_B.json'))},
    'resume_tested': False,
    'resume_note': 'no checkpoint/resume: every run is < 30 min and recomputes from scratch',
    'float_only_items': [],
}
json.dump(res, open('<LAB>/b5g/result_g2.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
