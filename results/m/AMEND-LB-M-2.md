# AMEND-LB-M-2 (amends BAR-LB-M, seal b1b0ee1cde): conditional full half-plane wording

Written 2026-11-16 by the lead, before any sealed run of BAR-LB-M. No count data has been seen.

**Why:** the Far Wall tail has since been certified, NON-COUNTING with respect to this BAR.
- File: engine/cert/fw_tail.py, log runs6/fw_tail_KT2000.log.
- Result: no physical eigenvalue for Re kappa >= 0 and |kappa| >= 2000, with Collatz 0.5407 and worst lower/rhs 14.77.
- The argument is in PROOF-B4-v4-AMEND2-TAIL.md.
- A hostile referee (seats/REFEREE-B4TAIL.md) is in progress.

**Change, to the wording only. Predictions, controls and kill conditions are unchanged.**
- If PM1-PM3 hold AND REFEREE-B4TAIL returns no WRONG and no unresolved GAP, then the verdict wording becomes: "certified by interval arithmetic: the Evans-Coleman critical solution has exactly one growing non-gauge mode (Re kappa >= 0), kappa* in [...], gamma in [...]".
- Otherwise the BAR-LB-M wording, which is scoped to |kappa| <= 2000, stands.
- Any fix the tail referee requires must be rerun, and its log must exist before the wider wording is used.

instruments: engine/cert/fw_tail.py, runs6/fw_tail_KT2000.log

## Seal fields
- formula_ref: BAR-LB-M.md:1-80
- calibrated_on: none (no count data); the tail certificate is disclosed above
- depends_on: b1b0ee1cde
- null_not_symmetry: unchanged
- target_form: unchanged
- parameters: unchanged
- tolerances: unchanged
- positive control: unchanged (the tail's own control: uniform bounds no better than the fw_cert3 values at K = 2000)
