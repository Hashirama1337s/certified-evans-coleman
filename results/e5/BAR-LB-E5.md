# BAR-LB-E5: existence of the Evans-Coleman critical solution, attempt 5 (v3.2 = v3.1 + R8)

Written 2026-11-16 by the lead, before any sealed v3.2 run.

**Why attempt 5.** E4 (seal 7928715efe) was settled FAILED by kill rule PE4: the k=1/2 warm-up never contracted.
- **Diagnosis** (lead, non-counting, runs4/diag_e4w_1e-15.log; verified by the build seat): the sonic tail bound booked the width of box-evaluated head coefficients, which is first order in r, as error.
- **Repair R8.** A mean-value head bound, |c_n(V0~)| + r·|c_n'(box)|, applied to eps and epsX (PROPOSAL-V32.md, sealed with this bar).
- **Referee round** required by E4: seats/REFEREE-R8.md. Driver fit to seal unchanged; only the control test_tails_v32.py was fixed.
- Parameters are unchanged from E3/E4.

**Disclosed prior knowledge:**
- the E3 and E4 numbers;
- the v3.2 NON-COUNTING smokes: k=1/2 contracts at r=1e-15 with ratio 0.0247 on a smoke box, and k=1/3 is unchanged;
- the blind B5 proof (BAR-LB-B5C HOLD, seal 71904980c4).

The claim bands are unchanged and come from published sources and exact algebra. No prediction is made about which rung contracts beyond "some rung of the sealed ladder".

instruments: engine/cert/run_v32.py, engine/cert/driver_v32.py, engine/cert/test_tails_v32.py, plus everything pinned under BAR-LB-E4

## Object
Unchanged (the Evans-Coleman solution, p = rho/3). Existence only; local uniqueness is NOT claimed.

## Predictions
Identical to BAR-LB-E4:
- **PE1:** radius ladder 1e-15..1e-9, first contraction.
- **PE2:** pe2_v31 in full.
- **PE3:** m/r(t=0) in 0.0596 ± 0.0005.
- **PE4:** k=1/2 warm-up from float64 seeds; the contracted box contains the exact Friedmann sonic data.

## Controls (any failure = VOID)
- **Content nulls** for W0, xc and V0, as in E4. Each must refuse via the Krawczyk test itself; an exception is VOID.
- **Instrument controls at the sealed parameters,** all PASS under this seal:
  - test_tails_v32.py 13 1e-12 and 12 1e-12. These include the injected 3·eps test that must flag, derivative-vs-central-difference check (C), and the frozen-eigenpair null, which must FAIL.
  - test_tails_v31.py v31.
  - test_lohner2.py thin and box.
  - check_struct.py.
- null2 is a smoke check only. No LB_* flag may be set.

**Disclosed, not fixed:**
- outward_leg_v2 uses the v3.1 head (rigorous, loose; k=1/3 only).
- choose_nu uses the v3.1 head; the selection is identical.
- The k=1/2 warm-up may have only one rung of margin.

## Kill conditions
- PE4 fails: STOP, no k=1/3 claim. The warm-up runs FIRST.
- Any run exceeds 2× its projection: warm-up 21 min, target 25, nullV0 14.
- A sixth attempt needs a new referee round and <OWNER>'s word.

## Verdict
HOLD if PE1-PE4 hold and all controls pass.

Wording: "existence certified by interval arithmetic (first implementation), independently proved by a refereed blind second implementation (B5); HKA 1996 asserted existence via a shooting argument; local uniqueness not claimed".

## Seal fields
- formula_ref: PROPOSAL-V32.md:1-78
- calibrated_on: E3 (VOID) and E4 (FAILED) runs, disclosed; v3.2 smokes runs5/smoke/ NON-COUNTING, disclosed; claim bands from published sources and exact algebra
- depends_on: 7928715efe, 8aec628fd2, 921909428f, 71904980c4
- null_not_symmetry: content nulls shift the box off the solution; PE4 uses exact Friedmann data
- target_form: first contracting radius + V0 interval; PE2 booleans/counts; PE3 interval; PE4 containment; null booleans; control PASS flags
- parameters: k = 1/3 (target), k = 1/2 (warm-up); radius ladder 1e-15 .. 1e-9; null offset 1e-13
- tolerances: interval arithmetic; radius cap within BAR-LB-E width <= 1e-6
- positive control: PE4 Friedmann warm-up
- toy_model: flat Friedmann
- files: BAR-LB-E5.md, PROPOSAL-V32.md
- command: memcap 8 py -u engine/cert/run_v32.py warmup|target|nullV0|null2
