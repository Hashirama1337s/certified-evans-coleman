# AMEND-LB-M-1 (amends BAR-LB-M, seal b1b0ee1cde): speedup citation, ceiling, fixed parameters

Written 2026-11-16 by the lead, after sealing BAR-LB-M and before any sealed run. No data has been seen since the seal.

- **Speedup citation (doctrine #33):** SPEEDUP-LB-M.md. The projection went from 21.5 h to 1.7 h. The 3.4 h figure in BAR-LB-M is a kill cap (2x the projection), not a projection.
- **ceiling:** none - this is an integer count, so the best possible outcome is an exact N. The remaining slack is the open tail |kappa| > 2000, which is analytic work.
- **k:** fixed - the object is the k = 1/3 radiation-fluid solution. The k = 1/2 tube enters only through control C3, which is already run.
- **R1:** fixed - R1 = 40 is set by the lowest refereed Far Wall band. A second value (R1 = 3, the positive control) is run under this seal.
- **Gauge null:** kappa_g is a genuine zero of the undeflated function, coming from a coordinate freedom. C1b tests that deflation removes exactly that zero. It is not a symmetry of D_new.

## Seal fields
- formula_ref: BAR-LB-M.md:1-80
- calibrated_on: none (no data has been seen since the seal)
- depends_on: b1b0ee1cde
- null_not_symmetry: as above
- target_form: unchanged
- parameters: unchanged
- tolerances: unchanged
- positive control: unchanged
