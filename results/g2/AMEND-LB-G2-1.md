# AMEND-LB-G2-1 (amends BAR-LB-G2, seal 8a1c27c7fb): gauge-independent positive control

Written 2026-11-17 by the lead, after sealing and BEFORE the blind seat launched. No data exists.

**Why:** a gauge (coordinate) mode's growth rate depends on the slicing and variables used. kappa_g = 0.3556992037 belongs to our first engine's formulation. A blind seat in a different gauge may honestly have a different gauge eigenvalue, or none. PG2 as sealed could therefore fail for no error. The physical kappa* (PG1) is gauge-invariant, so it is unaffected.

**Change:**
- **PG2 is replaced by PG2'**, a known-answer control on the seat's own machinery. The blind seat runs its eigenvalue pipeline, unchanged, on the exact flat Friedmann solution with the same k = 1/3 (in its self-similar form). It must enclose at least one perturbation growth rate that it ALSO derives in closed form by hand. The enclosure must contain the closed-form value.
- **The old PG2 is kept as DISCLOSED only.** The seat reports every real zero of its undeflated function on [0.40, 2.70], with its gauge or physical identification and the argument for it. Overlap with kappa_g is reported, not scored.
- **PG3 reads:** after the seat's own gauge zeros are removed, there is a certified constant sign on [0.40, 2.70], apart from the zero it identifies as kappa* if kappa* lies in that interval.
- PG1, the kill conditions and the verdict rule are otherwise unchanged. The HOLD requires PG1, PG2' and PG3.

## Seal fields
- formula_ref: BAR-LB-G2.md:16-19
- calibrated_on: none (no data)
- depends_on: 8a1c27c7fb
- null_not_symmetry: as in BAR-LB-G2
- target_form: closed-form Friedmann rate contained in the blind enclosure (exact containment in arb)
- parameters: unchanged
- tolerances: unchanged
- positive control: PG2' (the exact Friedmann rate)
- toy_model: none - PG2' is the known-answer case
