# PROPOSAL-V32: repair R8 (input to BAR-LB-E5)
Build seat, 2026-11-16. Proposal only. Every number below is a **NON-COUNTING** smoke (runs5/smoke/). No sealed file was edited. All new code lives in engine/cert/driver_v32.py, run_v32.py, test_tails_v32.py and smoke_v32_*.py. The sealed r=1e-15 chain has **NOT** been run.

## 1. Root cause, verified (runs5/smoke/tails_r8.log, dbox.log)
- **Mechanism.** param.tail_bound sums the defect head Σ_{K<n<N2} |[G(P_K)]_n| ν^n with abs_upper of coefficients evaluated over the V0 box (param.py:102-109). The width of those coefficients is first order in r, and it went into eps as error.
- **k=1/2 (ν_s 0.001) sonic eps.**
  - Point: 2.67e-42.
  - Box, v3.1: 4.32e-22 at r=1e-15, 5.26e-19 at 1e-12, 6.38e-13 at 1e-9, i.e. linear in r.
  - This reproduces the lead's chain number (diag_e4w: 4.3e-22).
- **k=1/3 (ν_s 0.004) sonic eps.** Point 1.86e-32. Box, v3.1: 2.91e-32 at 1e-15, 1.36e-29 at 1e-12. Also linear in r, but small next to the point value.
- **Why k=1/2 is so much worse: wrapping in the box Dual recursion, not true V0 dependence.**
  - At r=1e-15 the derivative enclosure of the order-15 coefficient, dY_{0,15}/dV0, has radius 2.8e26 at k=1/2, against 3.4e7 at k=1/3 (dbox.log).
  - The start derivative d_box has relative width 2.5e-7 at k=1/2, against 2e-11 at k=1/3.
  - (Refuted by REFEREE-R8 §5d: the k=1/2 bg.sonic_series ZeroDivisionError at high order was mp cancellation at 90 digits, not a resonance; at 200 digits orders 120/160 run.)
- **Centre manifold over the W0 box: no defect.** eps stays flat at 3.29e-56 (k=1/3) and 2.82e-57 (k=1/2) from r=0 to 1e-9.
- **epsX** had the same defect as eps and is fixed by the same proxy.

## 2. The change (R8): mean-value head. Clever-first, one function.
- **Bound.** For every V0 in [V0~ ± r]: |c_n(V0)| ≤ |c_n(V0~)| + r·|c_n'(box)|.
  - c_n(V0~) comes from the point coefficients (V0~ = mid, exact).
  - c_n'(box) comes from forward-mode Dual series arithmetic over the box. It is a rigorous enclosure of dc_n/dV0, valid by the mean value theorem in real V0.
  - Per coefficient, the code takes h_n = min(box ball, mean-value bound). Both are rigorous, so the head is ≤ the v3.1 head (the final eps comes from the doubling iteration, param.py:179-191, not proven monotone; measured eps_v32 ≤ eps_v31 in every smoke).
  - The derivative enclosure has width O(r), so the r-term is r·|true derivative| + O(r²). No separate Taylor-model remainder is needed. The start set already uses this same form for P_K itself (driver_v3.py:105-123).
- **Not changed.** The low-order check, the Cauchy tail beyond N2, the eigen-coordinate operator bound and the Lipschitz constant are all the pinned code, still evaluated over the box. Only the coefficients n > K are replaced, through a proxy tape. ctx.cap is raised before any series is built (doctrine #115).
- **Code (driver_v32.py):** R8 header :1-23; _ser (no truncation) :35; series_eval_dual :39-53; TMHeadTape proxy :61-94 (min at :85); _dual_series :97; sonic_manifold_v32 :104-128 (r=0: pinned tail_bound unchanged); MatcherV32.start_parts :131-147 (4a-4d, item 2, 5i-iii, R7, cache, nulls inherited); kr_split (diagnostic) :150.
- **run_v32.py** is run_v31.py with three changes: MatcherV32 (:53), an R8 assert (:37), and its banner. The sealed parameters, ν-ladder rule (pinned choose_nu, unchanged result: ν_s 0.004 at k=1/3, 0.001 at k=1/2), radius ladder, controls and nulls are unchanged.
- **4c switch tangent:** a derivative, so its width goes into J, never E. Linear in r at k=1/2 (rad 7.3e-5 at 1e-15, 0.072 at 1e-12); it feeds |I − C·J(X)|, O(r²) in K. No change.
- **epsX** is fixed by the same proxy, applied to the companion S. **No other first-order error terms found:** e0's dmax²·half_acc and eend's dmx²·half_acc2 are true second-order remainders; the centre manifold is flat in r.
- **Seed.** No Newton refinement is needed. At k=1/2 the Newton step C·R(x~) is (5.7e-18, -1.6e-17, -2e-18), ≪ r = 1e-15. The step lands within 4.7e-19 of the exact data (point_12.log). Using the identical float64-seed pipeline for both k is fine.

## 3. Smoke table (NON-COUNTING, runs5/smoke/)
| smoke | v3.1 | v3.2 (R8) |
|---|---|---|
| (a) k=1/2 Ebox, r=1e-15 | 7.17e-12 (≈7r) | **1.76e-22** (1.8e-7·r) |
| (a) k=1/2 Krawczyk, r=1e-15 | no contraction, ratio 2.66e5 | **contracts, ratio 0.0247.** E part 6.5e-6, J part 0.024. K contains the exact Friedmann data on all 3 components. V0 = -0.565685424949238020 ± 5.7e-19 (box_12_1e-15.log, 581 s) |
| (a) k=1/2 Ebox, r=1e-12 | 8.74e-9 | **1.88e-16** (1.9e-4·r). No contraction: ratio 24.4, from the J part 24.4; the E part is 0.007 |
| (b) k=1/3 Ebox, r=1e-12 | 2.58e-19, ratio 0.0114 | **2.57e-19, ratio 0.0114** (box_13_1e-12.log). Not worse |
| (b) k=1/3 sonic eps, r=1e-12 | 1.36e-29 | 1.86e-32 (= point eps) |
| point runs (unchanged path) | k=1/3: \|R\| 7.99e-18, Epoint 1.67e-24 | identical; k=1/2: \|R\| 1.24e-12, Epoint 8.66e-27 |
| (c) test_tails_v32.py 13 (sealed K, ν; box 1e-12; mp 90; REFEREE-R8 fixes) | | **PASS** (runs5/smoke/test_tails_v32_13_ref.log): (C) Dual derivative vs central difference 525/525 overlap, max rel err 7.1e-23; frozen-eigenpair null 0/525 overlap = FAILS as required. (B) enclosure 0/2625 and |c_n|≤h_n 0/2625 failures. (A) error ≤ 0.00197·eps; max allow/eps 1.16; INJ 3·eps shift flagged 75/75 |
| (c) test_tails_v32.py 12 (mp 200 digits, orders 120/160) | | **PASS** (test_tails_v32_12_ref.log): (C) 525/525 overlap (max rel err 14.5 on wide-ball coefficients; overlap is the criterion); frozen null 197/525 non-overlap = FAILS as required. (B) 0/2625, 0/2625. (A) error 9e-18·eps (point-tail control only, REFEREE-R8 §5b); max allow/eps 2.11; INJ flagged 75/75 |
| (c) test_tails_v31.py v31 (sealed control, rerun) | | **PASS** |

- **Reading the k=1/2 result.** The warm-up contracts at r = 1e-15, the first rung. The derivative-width term grows linearly with r (0.024 at 1e-15, 24.4 at 1e-12; about 0.24 at 1e-14 by extrapolation, POST HOC, UNTESTED: no prediction about 1e-14 is sealed). That is not a soundness problem: the ladder stops at the first contraction.

## 4. Sealed-run commands (after BAR-LB-E5 is sealed; all LB_* unset; run from <LAB>)
```
py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/run_v32.py warmup > runs6/warmup.log
py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/run_v32.py target > runs6/target.log
py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/run_v32.py nullV0 > runs6/nullV0.log
py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/run_v32.py null2  > runs6/null2_smoke.log   # smoke only
```
- **Chains.** A = warmup then target, so the PE4 kill rule fires first. B = nullV0 then null2. Never more than 2 jobs at once.
- **Controls under the seal:** test_tails_v32.py 13 1e-12; test_tails_v32.py 12 1e-12; test_tails_v31.py v31; test_lohner2.py thin and box; check_struct.py.

## 5. Runtime projections (kill at 2x)
| run | projection | basis |
|---|---|---|
| warmup | ~21 min, kill at 42 | point 641 s + one box rung 420-580 s; R8 adds ~3 s per residual |
| target | ~25 min, kill at 50 | unchanged from E4: 1504 s; the R8 tail adds ~2 s |
| nullV0 | ~14 min | |
| null2 | ~7.5 min | |
| controls | ~3 min total | |

- **2-minute probe:** the k=1/2 tau-step-50 line appears at about 50 s with width ≈ 1e-22 (point).

## 6. Files to pin
- **New:** engine/cert/driver_v32.py, run_v32.py, test_tails_v32.py (post-REFEREE-R8 version; hash at seal time).
- **Carried over, as in BAR-LB-E4:** driver_v31.py, test_tails_v31.py, test_lohner2.py, plus the E3/E4 list (PROPOSAL-V31.md §6).
- **Exclude:** smoke_v32_box.py, smoke_v32_tails.py, smoke_v32_dbox.py, runs5/smoke/*.

## 7. Open risks
1. **The k=1/2 warm-up has a single rung of margin.** It contracts at 1e-15 with ratio 0.0247, which is mostly the J width. If the sealed run differs from the smoke (it should be bit-identical apart from xbar, which the ladder recomputes from the box centre), PE4 fails at every rung.
   - The real fix would cut the wrapping in the box Dual recursion: a sharper eigen Krawczyk box (driver_v2.Problem2.eigen doubles the radius, driver_v2.py:384-397), or C1-Lohner on d_box. That is not proposed here.
2. **Mixed xbar.** The smoke took xbar from the point run in a separate process (exact man/exp transfer). The sealed ladder computes xbar in one process, which is the E4 logic, unchanged.
3. **outward_leg_v2 is not repaired.** It still calls sonic_manifold_v2 with the v3.1 head on the certified V0 box (radius ~1e-20). At k=1/3 the effect is harmless (eps 1.5e-26 or less, flat), but it is not repaired.
4. **(c) part (A) is a point-tail control**: it cannot see the R8 r-term (REFEREE-R8 §5b). The r-term is tested by (B) and (C). mp reference at k=1/2 now runs at 200 digits, orders 120/160.
5. **The v3.2 smoke timings ran under load from a parallel job.** k=1/2 box rungs took 581-597 s, against 420 s in E4.
