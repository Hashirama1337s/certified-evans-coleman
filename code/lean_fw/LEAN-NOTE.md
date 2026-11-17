# LEAN-NOTE: what of the Far Wall Collatz step is kernel-checked (2026-11-16)

Lean 4 (toolchain v4.36.0-rc2) + Mathlib, project <LEANPROJ>/k3lean. Every check ran with `lab.py lean` under memcap 12. A PASS means: exit 0, no error, no sorry, no forbidden word, and only the axioms propext, Classical.choice and Quot.sound.

**Wording to use:** the Collatz step of the Far Wall is kernel-checked on every bulk cell and on the tail. The matrix entries are ball-arithmetic enclosures; they are not Lean-checked. The Far Wall as a whole is **not** Lean-checked.

## 1. Kernel-checked

| file (lean_fw/) | sha256 (first 16) | ledger | content |
|---|---|---|---|
| fw_collatz_all.lean | df11e3c378bffa9c | **16cb6d873a** PASS (lead) | for each of the 133 bulk K-cells of [40, 55] (M = 2.5) and [55, 2000] (M = 3.5): 3 rational rows `a_i0*y0 + a_i1*y1 + a_i2*y2 < y_i` and `y > 0` |
| fw_collatz_tail.lean | 4f1dca963105ed0a | **6a2a9e1d2f** PASS | the same 4 statements for the uniform tail matrix (every K >= 2000, M = 3.5) |
| fw_collatz_logic.lean | bc6d8d78c7df3440 | **750fcfda8b** PASS | the implication used by Corollary 3.1 (8 theorems, listed below) |
| fw_collatz_logic_control.lean | b0bafd4762e4e3f4 | 3078ef4ec1 FAIL, as intended | teeth control: the core lemma with `A y <= y` instead of `A y < y`. This version is false (take A = I), and the checker rejects it |

### 1a. Tail rows (fw_collatz_tail.lean, as written over ℚ)
- A_up = [[0.219429, 0.0982062, 0.000364644], [0.962082, 0.109831, 0.243731], [0.866933, 0.0363838, 0.000330783]].
- y = (0.3063, 1, 0.5587).
- Rational Collatz ratio max_i (A_up y)_i / y_i = 0.540738. The arb Collatz bound is [0.540706 +/- 3.36e-8].
- Each A_up entry is a 6-significant-digit decimal rounded up. lean_tail.py asserts `arb(A_up[i][j]) > A[i][j].upper()` (strictly above) and A_up >= 0.
- Source: `fw_tail2.bulk_tail(2000.0)`. fw_tail2.py was imported unchanged (sha256 a2ecd2ae…, the same as in runs6/fw_tail2_KT2000.log). The same weights were reproduced: [3.125, 0.6, 250.0, 19/20, 11/20].
- Logs: lean_fw/lean_tail.log (extract, 45 s, peak 0.84 GB) and lean_fw/lean_tail_check.log.

### 1b. Logic (fw_collatz_logic.lean, namespace FarWallCollatz)
All entries are real, on `Fin n`; `*ᵥ` is the matrix-vector product, and every inequality is componentwise.
1. `collatz_sub_nonneg`: if A >= 0, y > 0, A y < y and A d <= d, then d >= 0.
2. `isUnit_one_sub`: under the same hypotheses on A and y, I - A is invertible.
3. `neumann_bound`: if v <= A v + b, then v <= (I - A)⁻¹ b. This holds for every b, with no sign condition on v.
4. `inv_mulVec_nonneg`: if b >= 0, then (I - A)⁻¹ b >= 0. `inv_nonneg`: (I - A)⁻¹ >= 0 entrywise.
5. `eigen_norm_lt_one`: every complex eigenpair (μ, x) of A with x ≠ 0 has ‖μ‖ < 1. This is ρ(A) < 1 in eigenvalue form; Mathlib's `spectralRadius` itself is not used.
6. `collatz_mono` (the bridge): if A <= A_up entrywise and A_up y < y with y > 0, then A y < y.
7. `far_wall_corollary_3_1` (Fin 3, the form the certificate uses). The inputs are:
   - Atrue, M real, with 0 <= Atrue <= M <= A_up entrywise;
   - A_up and y rational, with y > 0;
   - the rows written literally as `Aup i 0 * y 0 + Aup i 1 * y 1 + Aup i 2 * y 2 < y i` over ℚ. This is the shape of the cell files.
   - v >= 0 with v <= Atrue v + e₁.

   The conclusion is: I - M is invertible, (I - M)⁻¹ >= 0, every complex eigenvalue of M has modulus < 1, and **v <= (I - M)⁻¹ e₁**.

**Match to the paper.** PROOF-B4-v4.md:69 states: "If a positive y has 𝔄y < y (Collatz), then ρ(𝔄) < 1 and v <= (I - 𝔄)^{-1}e_1."

Downstream, `assemble` (engine/cert/fw_tail2.py) and fw_cert.py:361-365 (bulk cells) compute v = (I - 𝔄)⁻¹e₁ in arb, by cofactors from the ball matrix 𝔄, and take `up(·)`. The code then uses v only in the nonnegative combinations e_b = Σ c2_k v_k and n̄ = Σ c3_k v_k.

To apply theorem 7, take M = the matrix of arb ball uppers of 𝔄. M lies in the balls, so the arb cofactor value encloses (I - M)⁻¹e₁. Also M <= A_up holds by the asserted strict rounding. Theorem 7 then gives exactly the componentwise bound that the code uses.

Inputs the theorem takes as hypotheses (not proved in Lean):
- 0 <= Atrue <= M, meaning the true bootstrap coefficients are bounded by the ball uppers;
- v = (U, Z, N) is finite and nonnegative;
- the bootstrap inequality v <= Atrue v + e₁.

## 2. Outside Lean (not kernel-checked)
- **Entries.** The ball-arithmetic (arb/python-flint) computation of the matrix entries: fw_cert3 cells (bulk) and fw_tail2 (outer part, blown-up envelopes, closed-form deep integrals, ratio sups). These are CERTIFIED enclosures, not proofs in Lean.
- **The analytic estimates that produce the matrix.** These are PROOF-B4-v4 sec. 2-3 (rows 1-3, slaving, IBP, discretisation, sonic piece), AMEND1 and AMEND2-TAIL (monotonicity in u, blow-up, junctions). They include:
  - the bootstrap inequality v <= Atrue v + e₁;
  - finiteness of the sups;
  - Atrue <= 𝔄 (monotonicity of the entries in the coefficient bounds).
- **The (I - 𝔄)⁻¹e₁ evaluation and what follows.** This covers the arb evaluation of (I - 𝔄)⁻¹e₁ and the bounds e_b and n̄ computed from it.
- **The centre test.** This covers the centre zone (sec. 4: Volterra supersolution, Taylor models, flat kernel) and the μ-disc contradiction test on the 160 θ-cells, both in the bands and in the tail.
- **The [40, 55] band.** Its Collatz rows are among the 133 cells of fw_collatz_all.lean. Everything else about that band is likewise outside Lean.
- **Statement match.** Per lean-check step 5, a second seat still has to compare these Lean statements with Corollary 3.1, hypothesis by hypothesis, before any report calls them PROVED.

## 3. Reproduce (from <LAB>)
```
py -3.11 <MEMCAP> 8 py -3.11 -u lean_fw/lean_tail.py lean_fw/fw_collatz_tail.lean
py -3.11 <MEMCAP> 12 py -3.11 <TOOLS>/lab/tools/lab.py lean lean_fw/fw_collatz_tail.lean --project <LEANPROJ>/k3lean
py -3.11 <MEMCAP> 12 py -3.11 <TOOLS>/lab/tools/lab.py lean lean_fw/fw_collatz_logic.lean --project <LEANPROJ>/k3lean
py -3.11 <MEMCAP> 12 py -3.11 <TOOLS>/lab/tools/lab.py lean lean_fw/fw_collatz_logic_control.lean --project <LEANPROJ>/k3lean   (must FAIL)
```
