# B5G run book (exact commands, in order)

All runs from any directory; every Python call under memcap. Scripts live in `<LAB>/b5g/`.

```
# 1. own symbolic derivation (Einstein tensor + div T in polar-areal slicing, self-similar variables)
py -3.11 <MEMCAP> 8 py -3.11 -u <LAB>/b5g/derive_g2.py        # -> eqs_g2.json
# 2. centre analysis (Fuchs exponents of the perturbation system at z=0)
cd <LAB>/b5g
py -3.11 <MEMCAP> 8 py -3.11 -u centre_sym.py                                # -> centre_sys.json
# 3. sonic-line identities, FRW exactness check, expression export
py -3.11 <MEMCAP> 8 py -3.11 -u sonic_sym_g2.py                              # -> exprs_g2.json
py -3.11 <MEMCAP> 8 py -3.11 -u bgjac_g2.py                                  # -> bgjac_g2.json
# 4. certified runs (each < 30 min; A = main, B = convergence configuration)
py -3.11 <MEMCAP> 8 py -3.11 -u cert_g2.py A > cert_A.log 2>&1              # -> cert_A.json, sign_intervals_A.json
py -3.11 <MEMCAP> 8 py -3.11 -u cert_g2.py B > cert_B.log 2>&1              # -> cert_B.json
# 5. assemble deliverable
py -3.11 <MEMCAP> 8 py -3.11 -u make_result_g2.py                           # -> result_g2.json
```

Probe (doctrine #32): `cert_g2.py A probe` (sign range shrunk to [0.95, 1.10]) ran in 205 s and passed every check
except the FRW interval-derivative at eps=1e-4 (fixed: det via Leibniz expansion, eps=1e-5); log `probe_A.log`.

Module map: `tm.py` (interval Taylor models), `pipeline_g2.py` (codegen, nonlinear/linear singular solvers, sonic
slope + Krawczyk, perturbation coefficient matrices), `lohner_g2.py` (Lohner QR propagation, background build,
F(kappa) and F'(kappa)), `frw_exact_g2.py` (closed-form FRW background, used only at the FRW sonic point).
Development scripts (not part of the certified chain): `test_g2.py`, `diag_g2.py`, `dbg*_g2.py`.

Convergence refinement of kappa* (both configurations; certified bisection to width ~1e-13):
```
py -3.11 <MEMCAP> 8 py -3.11 -u refine_g2.py A > refine_A.log 2>&1         # -> refine_A.json
py -3.11 <MEMCAP> 8 py -3.11 -u refine_g2.py B > refine_B.log 2>&1         # -> refine_B.json
```
(run refine before make_result_g2.py, which reads refine_A/B.json)
