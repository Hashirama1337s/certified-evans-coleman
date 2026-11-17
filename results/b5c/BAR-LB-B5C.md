# BAR-LB-B5C: closure checks for the blind second implementation (B5), after REFEREE-B5

Written 2026-11-16 by the lead, before the closure script exists.
- Source of the items: seats/REFEREE-B5.md (sections 3, 4, 5, 7, 9).
- The B5 proof run (b5/RESULT.md) is NOT re-run. These checks act on its cached certified objects, plus one new FRW outer-leg control.

instruments: b5/close_b5.py (new; pinned before run) plus the existing b5/ modules it imports

## Predictions (each a certified interval statement; all must hold)
- **C1. Outer sonic stretch.** On theta in [0, theta_o] of the certified sonic curve: S/theta != 0, P < 0, mu < 1 and V > 0. Bisect up to depth 12 if ambiguous; the number of ambiguous pieces left must be 0.
- **C2. Miranda preconditioner.** det C != 0 is certified by arb.
- **C3. Centre step.** P < 0 on z in [0, 0.05], from the cached certified centre step.
- **C4. Cache integrity.** Every cached object fed to Miranda (G_lo, G_hi, G_box, G_mid) contains its reloaded copy. Each cache key carries a hash of the code that produced it.
- **C5. FRW outer leg (positive control).** The k = 1/2 FRW run with the outer leg enabled reaches t = 0, and its enclosure contains the exact FRW values at t = 0.

## Controls
- **Positive.** C5 is the positive control.
- **Negative.** C1 run on a deliberately shifted theta interval past the sonic curve's next zero of S/theta (or of V, if that comes first) must report a sign change (not pass). This shows the check can fail.

## Kill conditions
- Any check is still ambiguous after depth 12: report it as GAP, not HOLD.
- A run exceeds 30 min.

## Verdict
- HOLD if C1-C5 all hold and the negative control fails as required.
- Wording: "B5 independently proves a k = 1/3 solution (centre-regular, analytic saddle sonic crossing, one V zero, reaches t = 0); uniqueness not claimed; agreement with the first implementation is POST HOC".

## Seal fields
- formula_ref: seats/REFEREE-B5.md:1-120
- calibrated_on: none - the predictions are qualitative sign and containment statements derived from the referee report; no numbers are tuned
- depends_on: 7928715efe
- null_not_symmetry: the negative control moves the theta interval off the claimed stretch (a data change)
- target_form: booleans per check, plus ambiguous-piece counts
- parameters: k = 1/3 (C1-C4); k = 1/2 (C5); bisection depth 12
- tolerances: interval arithmetic
- positive control: C5 FRW outer leg
- files: BAR-LB-B5C.md
- command: memcap 8 py -u b5/close_b5.py all
