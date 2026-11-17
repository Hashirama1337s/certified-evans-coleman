# BAR-LB-M: certified count of growing modes of the Evans-Coleman critical solution, |kappa| <= 2000

Written 2026-11-16 by the lead, before any sealed count run.

**Object.** The linear perturbations of the certified k = 1/3 solution (BAR-LB-E5, HOLD d0af75803f). We count zeros of the gauge-deflated characteristic function D_new (zeros = physical non-gauge modes by Thm 1.6 of PROOF-B4.md, given c_V(x_m) != 0 and A-bar_ss'(0) != 0, both certified at the start of every run).

**Method:**
- Argument principle on a half-disc of radius R1 = 40, upper path only, using conjugate symmetry. The 8-chord arc has vertices at Rp = R1/cos(dphi/2)(1+1e-12); every chord distance is >= 40 and every vertex is <= 55, both asserted in exact arithmetic.
- kappa-Taylor-model transition matrices in the characteristic frame A1 (lin_tm_v3.py, LinearTM3).
- Point-value argument sharpening.
- Fixed refinement: at most 3 halvings; any failure means INCONCLUSIVE.
- A fresh ledger with a settings header.

**Combined with:** the Far Wall certificates, which are pinned as files, not rerun:
- runs6/fw_cert3_M2.5_Ka40_K55.log: no physical eigenvalue for Re kappa >= 0, 40 <= |kappa| <= 55;
- runs6/fw_cert3_M3.5_Ka55_K2000.log: the same for 55 <= |kappa| <= 2000.

These are refereed in REFEREE-B4FW.md and REFEREE-B4FW2.md, with the PROOF-B4-v4 and AMEND1 mathematics.

**Disclosed prior knowledge (all NON-COUNTING):**
- build and speedup smokes: b4build2/, b4speed/, including the R1 = 3 box and half-disc dry runs (N = 1 each);
- the controls C1-C3 and N1-N3;
- kappa* = 2.8105525488 +- 8.6e-11 from a smoke bisection;
- the published literature (Koike-Hara-Adachi 1995 and Hara-Koike-Adachi 1996), which reports exactly one growing mode with gamma ~ 0.3558 and a gauge mode.

The prediction below is the literature's claim. No parameter was tuned on the count itself.

**Cost driver** (doctrine #33):
- The cost is the number of discs on the imaginary axis.
- The speedup pass (PROPOSAL-B4-SPEEDUP.md) brought the projection from 21.5 h to 1.5-1.7 h with two levers: the characteristic frame, and the Far Wall lowered from 55 to 40.
- The projection for this run is 1.7 h on one process (99 discs), under the 2 h line.

instruments: engine/cert/count_b4v3.py, engine/cert/count_b4v3_seal.py, engine/cert/lin_tm_v3.py, b4speed/cacheL.py, b4speed/L13_nsub2.pkl, engine/cert/count_b4.py, engine/cert/lin_tm_v2.py, engine/cert/lin_tm.py, engine/cert/lin_stages_v2.py, engine/cert/lin_stages.py, engine/cert/core.py, engine/cert/lin_mc.py, engine/cert/lin_fields2.pkl, engine/cert/lin_tube_v2.py, engine/cert/lin_tube2_k1_3.pkl, engine/cert/lin_tube2_k1_2.pkl, engine/cert/lin_b4_smoke3.py, engine/cert/fw_cert3.py, engine/cert/lin_fw4_bg3.pkl, engine/cert/lin_fw4_bg3_make.py, runs6/fw_cert3_M2.5_Ka40_K55.log, runs6/fw_cert3_M3.5_Ka55_K2000.log, plus everything pinned under BAR-LB-E5

## Predictions
- **PM1 (the count).** N(R1 = 40) = 1. In words: D_new has exactly one zero with Re kappa >= 0 and |kappa| < 40, with Delta/pi within 0.25 of 1 and the parity check holding.
  - Combined with the Far Wall, this gives exactly one physical non-gauge growing mode with Re kappa >= 0 and |kappa| <= 2000.
- **PM2 (the exponent).** `count_b4v3.py kstar` encloses kappa* in 2.8106 +- 0.0005, so gamma = 1/kappa* lies in 0.3558 +- 0.0001. This band is set by the published value, not by our smoke.
- **PM3 (preconditions, printed at start).** All must hold:
  - c_V(x_m) != 0;
  - A-bar_ss'(0) != 0;
  - the detA root is < 0;
  - iw = 2 at every plan point.

## Controls (run FIRST, under this seal; any failure = VOID)
- **Positive:** the half-disc count at R1 = 3 on a fresh ledger gives N = 1 with parity True.
- **Gauge null and positive:**
  - `count_b4v3.py smoke 1/3 C1`: the non-deflated D_HKA winds once around kappa_g (C1a);
  - the deflated D_new winds zero times there (C1b; the null);
  - the bisection contains kappa_g (C1c).
- Every guard in count_b4v3_seal is active. LB_R_FIXED, LB_MAXDISCS and LB_SHAPE are unset for the main count.

## Kill conditions
- Any INCONCLUSIVE disc or parity mismatch gives FAILED. No count is claimed.
- The main count is killed if it runs past 3.4 h (2x projection); controls are killed past 2x their own times (C1 about 10 min, R1 = 3 about 30 min).
- A second attempt needs a new referee round.

## Verdict
HOLD if PM1-PM3 hold and all controls pass.

Wording: "certified by interval arithmetic: the Evans-Coleman critical solution has exactly one growing non-gauge mode with |kappa| <= 2000 (Re kappa >= 0); the tail |kappa| > 2000 is open (analytic work in progress); gamma enclosed in [...]". Never "exactly one growing mode" without the scope.

## Seal fields
- formula_ref: PROPOSAL-B4-SPEEDUP.md:1-267
- calibrated_on: NON-COUNTING smokes in b4build2/ and b4speed/ (disclosed); the claim bands come from the published literature
- depends_on: d0af75803f
- null_not_symmetry: C1b tests a real zero of the undeflated function (data), not a relabelling
- target_form: integer N plus Delta/pi interval; kappa* interval; precondition balls; control booleans
- parameters: k = 1/3; R1 = 40; 8-chord arc; radius rule 2:0.25,5:0.4,90:0.5, right side 1.0; frame A1; prec 128; cap 32; REFINE 3
- tolerances: interval arithmetic; the 0.25-of-an-integer rule
- positive control: half-disc count at R1 = 3 (N = 1); C1a
- toy_model: none (controls are the known gauge zero and the known kappa*)
- files: BAR-LB-M.md, PROPOSAL-B4-SPEEDUP.md
- command: memcap 8 py -u engine/cert/count_b4v3.py count 40 runs7/count_R40.jsonl ; kstar ; controls as above
