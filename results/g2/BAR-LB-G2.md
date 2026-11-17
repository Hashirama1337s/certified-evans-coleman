# BAR-LB-G2: blind second enclosure of the relevant eigenvalue kappa* (and so gamma)

Written 2026-11-17 by the lead, before the blind seat launches. The blind seat has seen no data.

**Why:** doctrine #31. Existence has two independent proofs (E5 and B5C), but kappa* and gamma (BAR-LB-M, HOLD b1b0ee1cde) come from one engine. Matching the published 0.3558019 is consistency, not independence. The where-else pass (WHERE-ELSE.md, sec. 5) asked for this. This test does not change any count claim.

**Object.** Spherical linear perturbations of the k = 1/3 Evans-Coleman solution. The seat may use the blind B5 background in b5/, which was already an independent formulation. The seat must:
- derive its own perturbation equations from Einstein-Euler, with its own variables and gauge;
- impose regularity at the centre and analyticity at the sonic point;
- certify a real growth rate kappa*, where perturbations go as exp(kappa tau) and tau = -ln(-t), together with its gauge-mode bookkeeping.

The brief gives no target number, no predicted form and no paths to our first engine.

instruments: b5g/ (written by the blind seat; every file pinned by sha256 with lab.py pin before its first certified run)

## Predictions
- **PG1.** The blind seat returns an interval enclosure of a real non-gauge growth rate, with width <= 1e-6. That enclosure overlaps kappa* = [2.810552549 +- 3.83e-10] (BAR-LB-M). Equivalently, gamma = 1/kappa* overlaps [0.3558019224 +- 6.97e-11].
- **PG2 (positive control).** The blind seat also encloses a real eigenvalue that it identifies as a gauge (coordinate) mode, and that enclosure overlaps kappa_g = [0.3556992037 +- 7.41e-11] (C1c, BAR-LB-M controls).
- **PG3 (null).** After removing gauge modes, the blind function has a certified constant sign on the real interval [0.40, 2.70]: no real growing mode other than kappa* below 2.70. This is predicted by the M count (N = 1 in the half-disc).

## Kill conditions
- An enclosure disjoint from ours, on PG1 or PG2, gives FAILED. Both engines then get a hostile referee, and no gamma is claimed until that is resolved.
- If the seat cannot certify (it returns floats only), the test is INCONCLUSIVE and is reported as such. It is not a HOLD.
- Compute must be under 30 min per run (memcap 8).

## Verdict
HOLD if PG1, PG2 and PG3 hold. The scorer reads only the seat's result JSON and logs. Overlap is checked in arb.

## Seal fields
- formula_ref: BAR-LB-M.md:35-43
- calibrated_on: none (the blind seat has not run)
- depends_on: b1b0ee1cde, 71904980c4
- null_not_symmetry: PG3 is a data-level sign test on a real interval, not a relabelling. Gauge removal deletes one known real zero (PG2); it is not a symmetry of the blind function
- target_form: interval for kappa*; interval for kappa_g; a sign certificate on [0.40, 2.70]
- parameters: k: fixed - the object is the k = 1/3 solution of the release; width cap 1e-6; real interval [0.40, 2.70]
- tolerances: interval overlap in arb, because both sides are rigorous enclosures, so any honest pair must overlap with no headroom needed. The width cap 1e-6 is 2500x looser than ours (7.7e-10), to give a blind seat headroom
- ceiling: none - this is an overlap test of two rigorous enclosures; the best possible outcome is overlap
- positive control: PG2 (the gauge eigenvalue)
- toy_model: none - the gauge eigenvalue (PG2) is the known-answer control on the same object
- files: BAR-LB-G2.md
- command: the blind seat's own (recorded in b5g/RUN.md), under memcap 8
