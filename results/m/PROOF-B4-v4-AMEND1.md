# PROOF-B4-v4-AMEND1: referee repairs (REFEREE-B4FW) and the [40, 55] band

Proof seat, 2026-11-16. New file. PROOF-B4-v4.md is unchanged; this amendment supersedes its sec. 0 result line, sec. 6.5/6.6 numbers and the "Lean" row. Labels as in v4.

## 0. Amended result

| band (Re κ >= 0) | M | instrument | verdict | worst lower/rhs | worst absolute margin | log |
|---|---|---|---|---|---|---|
| 40 <= abs(κ) <= 55 | 2.5 | fw_cert3 (fw_cert2 θ/K machinery, ns ladder 150/300/600, bisection depth 2) | **CERTIFIED** | 1.010825 at K in [40, 41.2] | 7.48e-4 | runs6/fw_cert3_M2.5_Ka40_K55.log |
| 55 <= abs(κ) <= 2000 | 3.5 | fw_cert3 (ns 150, no bisection) | **CERTIFIED** | 1.314537 at K in [55, 56.65] | 6.41e-3 | runs6/fw_cert3_M3.5_Ka55_K2000.log |

- Union: **no physical eigenvalue with Re κ >= 0 and 40 <= abs(κ) <= 2000** (CERTIFIED, fw_cert3). Runs: 858 s and 922 s; 122 K-cells all contractive on [55, 2000] (worst Collatz 0.5682), 11 on [40, 55] (worst 0.7566).
- Different M is used per band. This is legitimate: the lemma holds for each fixed M on each band separately.
- **Superseded:** the earlier logs runs6/fw_cert_Ka55*.log, fw_cert_Ka60.log and fw_cert2_Ka40_K55.log. They were run with the cap leak (sec. 1) and are **void as certificates**. Their conclusions are reproduced by the fw_cert3 runs above.
- Still OPEN: the tail abs(κ) > 2000; a referee pass on this amendment. **The half-plane count is still not rigorous.**

## 1. The cap leak (doctrine #115) and its repair (`engine/cert/fw_cert3.py`, NEW; fw_cert/fw_cert2 not edited)

**Defect (REFEREE-B4FW B1).** fw_cert.py never raised python-flint's series cap (default 10). As a result:
- `.integral()` in the Picard jets truncated at order 10;
- `exp()` on the centre-manifold series truncated at order 10;
- the default `c(L, k) = 0` therefore silently dropped the order-11 cell remainders and the order-22 centre-manifold remainders (and the terms 10..21).

**Repair.** fw_cert3 does the following.
1. It sets `ctx.cap = 32` before anything is built, and asserts it at the end of the run.
2. It replaces every coefficient access with `coef(s, k)`. This **raises** if k >= s.prec (an unknown coefficient). An index below s.prec that is absent from coeffs() is an exact zero.
3. It asserts Lemma 4.1 in ball arithmetic: C0(0) ∋ C00, C1(0) ∋ 0, C1'(0) ∋ C11, v'(0) ∋ 1. It also asserts that every coefficient dropped by the centre-manifold Taylor models (k < m) contains 0: 128 checks, all pass.
4. It reads the background from `lin_fw4_bg3.pkl`, written by `lin_fw4_bg3_make.py` (NEW):
   - every stored radius r becomes ceil20(r(1 + 1e-9) + abs(mid)·1e-69), rounded up as a decimal;
   - exact integer balls with radius 0 stay exact;
   - all 4826 new balls were checked to contain the old ones.
5. It builds the centre data at 1.001·v_need.
6. It bounds the switch sliver (sec. 2).
7. It prints the worst lower/rhs alongside the absolute margin.

**Positive control (REFEREE capchk4 test), PASS** (runs6/fw_cert3_capchk.log). On the referee's 7 cells (5 x-leg, 2 tau-leg):
- the fw_cert3 box (cap 32) contains the correctly built cap-13 box, 7/7;
- the as-run cap-10 box does not contain it, 7/7 (the bug is reproduced);
- a strict request beyond the known order raises.

**Consequence found by the repair: the old deep cells were tight only because of the bug.**
- With the honest order-11 remainder, the jets on the 30%-wide E5 a-priori boxes reach ~1e21 (cell x ≈ -3.0, h = 0.016).
- The one-step Taylor box then collapses to the a-priori box. The bulk iteration first produced absurd values there, and then a non-finite matrix.
- **Fix (PROVED valid for any n):** each step is cut into n = 2^j equal sub-steps, with n the smallest value making coef_11(rough)·(h/n)^11 <= 1e-13.
  - Each sub-step restarts from a thin Taylor endpoint and uses the same a-priori box, because the trajectory over every sub-step lies in it.
  - The sub-step boxes are united.
- With this, the bulk reproduces the old Collatz values exactly: 0.7566 on [40, 41.2] at M = 2.5.

**Weights.**
- The bootstrap weights are a free choice; any fixed positive weights are admissible (v4 sec. 3).
- The float power iteration that picks them is fragile. fw_cert3 tries iteration counts 6, 3, 4, 2, 8, 10, 1 and keeps the result with the smallest rigorous Collatz bound.

**Centre robustness.** If the float fixed point diverges, there is no supersolution at that ns. fw_cert's diagnostic line then crashed on an empty max(). fw_cert3 treats this as a **failure** of that cell at that ns, and fw_cert2's ladder/bisection takes over unchanged.

## 2. The tau/x switch sliver (REFEREE-B4FW A, GAP) — PROVED, evaluation CERTIFIED

The 1e-8 inflation of the last tau box is **removed**. It is replaced by an explicit bootstrap.

- Let x_end be the x-coordinate at the end of the last tau step. It is enclosed thinly by a Taylor model with a Lagrange remainder on the step box: x_end ∈ -0.050000000 ± 7e-14.
- Let ℓ := max(0, upper(x_end - x_switch)). In the run, ℓ <= 5.41e-14.
- Let B be the last tau sub-cell box + 1e-6, and F := sup_B abs(f_bg) = 0.727.
- Every state on the sliver [x_switch, x_switch + ℓ] is within δ = ℓF <= 3.93e-14 of the landing state, which lies in the sub-cell box. The usual a-priori argument applies: δ < 1e-6, so the trajectory never leaves B.
- **Assignment of cells.**
  - The last tau sub-cell keeps its Taylor model on [x_switch + ℓ, x_hi]; this is valid because x_end <= x_switch + ℓ.
  - The first x-leg cell's Taylor model is extended backwards over the sliver. It is still expanded at the thin switch set (x = x_switch), with τ ∈ [h_x, ℓ], and its Lagrange remainder is taken on hull(cell box, landing box + δ).
  - The cells stay contiguous, and every point is covered by an enclosure valid there.

## 3. Text fix: bulk_ia range (REFEREE-B4FW A, rows)

- PROOF-B4-v4 sec. 3 defines U, Z, N on [x_c, 0] with x_c = ln(M/K). For a K-cell [K_a, K_b], bulk_ia instead integrates and takes sups over every cell meeting [x_c(K_b), 0] ⊇ [x_c(K), 0]. This is valid for two reasons.
  - (i) All integrands are nonnegative coefficient bounds, and only s >= t >= x_c(K) enter the representation. A cell sup over a superset is still an upper bound.
  - (ii) The normalisation is e(s) = abs(E_+(s, x_c(K))). Moving the base point only multiplies every bound by the same factor, since abs(E_+(t,s)) e(s)/e(t) = 1. So the ratios entering 𝔄 do not depend on where in [x_c(K_b), x_c(K_a)] the true x_c lies.
- The outputs e_b and n̄ take the sup over cells meeting [x_c(K_b), x_c(K_a)], which contains x_c(K) for every K in the cell.

## 4. Lean (narrow)

Regenerated with fw_cert3: `engine/cert/fw3_collatz.lean`. It states the Collatz rows A_up·y < y and y > 0 for **one** K-cell, [55, 56.65] at M = 3.5. The IA Collatz bound there is 0.568212; the entries are decimal upper bounds, strictly above the interval uppers (asserted). Result: PASS, kernel-checked, standard axioms only, lab ledger **ae7bc39344** (log runs6/fw3_lean.log). The old extract 17158af6e3 came from the leaky run and is superseded.

Only these rational inequalities are kernel-checked. The implication "Collatz ⇒ ρ < 1", the other K-cells and the [40, 55] band are **not** formalised. **The Far Wall must not be described as Lean-checked.**

## 5. Commands (one process at a time, memcap 8)
```
py -3.11 .../memcap.py 8 py -3.11 engine/cert/lin_fw4_bg3_make.py
py -3.11 .../memcap.py 8 py -3.11 -u engine/cert/fw_cert3.py capchk
set FW2_NS_LADDER=150,300,600
py -3.11 .../memcap.py 8 py -3.11 -u engine/cert/fw_cert3.py run 2.5 40 55 160 150 1500 2      (858 s)
py -3.11 .../memcap.py 8 py -3.11 -u engine/cert/fw_cert3.py run 3.5 55 2000 160 150 1500 0
py -3.11 .../memcap.py 8 py -3.11 -u engine/cert/lin_fw4_lean_extract3.py 3.5 55 56.65 engine/cert/fw3_collatz.lean
py -3.11 .../memcap.py 12 py -3.11 <TOOLS>/lab/tools/lab.py lean engine/cert/fw3_collatz.lean --project <LEANPROJ>/k3lean
```
Lean needs 12 GB: the lab tool refuses below that, because Mathlib maps about 9 GB.

## 6. Status changes against PROOF-B4-v4 sec. 7

| step | new status |
|---|---|
| cell / centre-manifold enclosures | CERTIFIED with cap 32 + strict coefficients + sub-stepping (was: leaky) |
| background radii | rounded up (bg3) |
| Lemma 4.1 in code | asserted (was: assumed) |
| switch sliver | PROVED bootstrap + CERTIFIED (was: GAP) |
| Far Wall on [40, 55] (M = 2.5) and [55, 2000] (M = 3.5) | **CERTIFIED** (fw_cert3 runs, sec. 0) |
| tail abs(κ) > 2000 | OPEN |
| referee of this amendment | OPEN |
