# AMEND-LB-B5C-1 (amends BAR-LB-B5C, seal 71904980c4): C5 and neg made well-posed; C3 wording

Written 2026-11-16 by the lead, AFTER the build seat's non-counting smokes of every mode and BEFORE the sealed run. POST HOC with respect to those smokes, disclosed.

## Why
- **C5 cannot hold as sealed, and this is the lead's error.** Flat FRW at k = 1/2 never reaches t = 0 (N = 0) in this chart. z peaks at 9/(4 sqrt 5), and the chart ends at the apparent horizon V = -1 (z = 1, 2m/r = 1). This was already known (LOG 11-15, AMEND-LB-E-1: "chart edge (apparent horizon) ... outward leg not in PE4").
- **neg as sealed is not certifiable.** S/theta and V have no zero inside the validated sonic disk |theta| <= 0.05.

## Changes
- **C5'.** The FRW k = 1/2 outer leg runs from the exact sonic point to z = 1. Its enclosure must contain the exact values Omega = 24/25 and V = -4/5, with mu < 1 and P < 0 certified on every step.
  - The t = 0 leg (the far N-system) has NO exact-solution positive control. This is DISCLOSED as a remaining gap. Its only check is the 50-digit mpmath cross-check, which is not a control.
- **neg'.** The same V > 0 test used by C1 is run on the certified inner solution across its V zero (z ~ 0.4528-0.4541). It must report a certified sign change. This shows that the sign test can fail.
- **C3.** "From the cached centre step" becomes "re-validated over the whole sigma0 box on z in [0, 0.05]", because the cache stores only end values.
- **C4.** The HOLD rests on the EC cache. The FRW cache predates the current prove.py and is disclosed as such.
- **Verdict wording.** Unchanged, plus: "the far t = 0 leg lacks an exact positive control".

## Seal fields
- formula_ref: BAR-LB-B5C.md:1-46
- calibrated_on: build-seat smokes of every mode (disclosed above); no numbers tuned
- depends_on: 71904980c4
- null_not_symmetry: neg' tests a real V zero of the certified solution (data, not relabelling)
- target_form: booleans and ambiguous counts as in the BAR
- parameters: unchanged
- tolerances: interval arithmetic
- positive control: C5' (FRW to the horizon z = 1)
