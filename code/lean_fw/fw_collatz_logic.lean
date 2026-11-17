import Mathlib

/-! # Far Wall: the logic step of the Collatz test (PROOF-B4-v4, Corollary 3.1)

Source: <LAB>/PROOF-B4-v4.md line 69:
  "Corollary 3.1. If a positive y has 𝔄y < y (Collatz), then ρ(𝔄) < 1 and v <= (I - 𝔄)^{-1}e_1."
Downstream use (engine/cert/fw_tail2.py `assemble`, same in fw_cert2/fw_cert3): v := (I - 𝔄)^{-1} e_1 is
computed in arb from the ball matrix 𝔄 and then enters e_b, n̄ only through nonnegative combinations
(sum_k c_k v_k with c_k >= 0). So the needed statement is the componentwise bound v_i <= ((I - M)^{-1} e_1)_i
for a real matrix M that lies in the arb balls (e.g. the matrix of the ball uppers).

Here v = (U, Z, N) are the bootstrap sups (finite, nonnegative), satisfying the bootstrap inequality
v <= A_true v + e_1 with 0 <= A_true <= M entrywise (all entries of 𝔄 are monotone upper bounds).

Contents:
* `collatz_sub_nonneg`  - core: A >= 0, y > 0, A y < y, A d <= d  ==>  d >= 0.
* `isUnit_one_sub`      - I - A is invertible.
* `neumann_bound`       - v <= A v + b  ==>  v <= (I - A)^{-1} b   (any n, any b; no sign condition on v).
* `inv_nonneg`, `inv_mulVec_nonneg` - (I - A)^{-1} >= 0 entrywise; b >= 0 ==> (I - A)^{-1} b >= 0.
* `eigen_norm_lt_one`   - every complex eigenvalue μ of A has |μ| < 1  (i.e. ρ(A) < 1).
* `collatz_mono`        - bridge: 0 <= A <= A_up and A_up y < y  ==>  A y < y.
* `far_wall_corollary_3_1` - Corollary 3.1 on Fin 3 in the exact form used, with the rational A_up and y of the
  cell files (rows written literally as  a0*y0 + a1*y1 + a2*y2 < y_i  over ℚ, as in fw_collatz_all.lean /
  fw_collatz_tail.lean).
-/

open Matrix

namespace FarWallCollatz

variable {n : ℕ}

/-- Core lemma (Collatz / sub-solution argument). -/
theorem collatz_sub_nonneg (A : Matrix (Fin n) (Fin n) ℝ) (y d : Fin n → ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (A *ᵥ y) i < y i)
    (hd : ∀ i, (A *ᵥ d) i ≤ d i) : ∀ i, 0 ≤ d i := by
  by_contra h
  push_neg at h
  obtain ⟨i0, hi0⟩ := h
  obtain ⟨k, -, hk⟩ := Finset.exists_max_image (Finset.univ : Finset (Fin n))
    (fun i => -d i / y i) ⟨i0, Finset.mem_univ _⟩
  have ht0 : 0 < -d k / y k := by
    have h1 : 0 < -d i0 / y i0 := div_pos (by linarith) (hy i0)
    exact lt_of_lt_of_le h1 (hk i0 (Finset.mem_univ _))
  have hbound : ∀ j, -(-d k / y k * y j) ≤ d j := by
    intro j
    have h2 := hk j (Finset.mem_univ _)
    rw [div_le_iff₀ (hy j)] at h2
    linarith
  have hdk : -d k / y k * y k = -d k := div_mul_cancel₀ _ (hy k).ne'
  have key : 0 ≤ (A *ᵥ d) k + (-d k / y k) * (A *ᵥ y) k := by
    have e : (A *ᵥ d) k + (-d k / y k) * (A *ᵥ y) k
        = ∑ j, A k j * (d j + (-d k / y k) * y j) := by
      simp only [Matrix.mulVec, dotProduct]
      rw [Finset.mul_sum, ← Finset.sum_add_distrib]
      apply Finset.sum_congr rfl
      intro j _
      ring
    rw [e]
    exact Finset.sum_nonneg (fun j _ => mul_nonneg (hA k j) (by linarith [hbound j]))
  have hlt := mul_lt_mul_of_pos_left (hc k) ht0
  linarith [hd k]

/-- `I - A` is invertible. -/
theorem isUnit_one_sub (A : Matrix (Fin n) (Fin n) ℝ) (y : Fin n → ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (A *ᵥ y) i < y i) :
    IsUnit (1 - A) := by
  rw [← Matrix.mulVec_injective_iff_isUnit]
  intro x z hxz
  have h0 : (1 - A) *ᵥ (x - z) = 0 := by
    rw [Matrix.mulVec_sub]
    exact sub_eq_zero.mpr hxz
  have hfix : A *ᵥ (x - z) = x - z := by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec] at h0
    exact (sub_eq_zero.mp h0).symm
  have p := collatz_sub_nonneg A y (x - z) hA hy hc (fun i => le_of_eq (by rw [hfix]))
  have q := collatz_sub_nonneg A y (-(x - z)) hA hy hc
    (fun i => le_of_eq (by rw [Matrix.mulVec_neg, hfix]))
  funext i
  have p' := p i
  have q' := q i
  simp only [Pi.sub_apply, Pi.neg_apply] at p' q'
  linarith

/-- Neumann-type bound: if `v ≤ A v + b` componentwise then `v ≤ (I - A)⁻¹ b` componentwise. -/
theorem neumann_bound (A : Matrix (Fin n) (Fin n) ℝ) (y v b : Fin n → ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (A *ᵥ y) i < y i)
    (hv : ∀ i, v i ≤ (A *ᵥ v) i + b i) : ∀ i, v i ≤ ((1 - A)⁻¹ *ᵥ b) i := by
  have hU : IsUnit (1 - A).det :=
    (Matrix.isUnit_iff_isUnit_det _).mp (isUnit_one_sub A y hA hy hc)
  obtain ⟨x, hxd⟩ : ∃ x : Fin n → ℝ, x = (1 - A)⁻¹ *ᵥ b := ⟨_, rfl⟩
  have hx : x - A *ᵥ x = b := by
    have h1 : (1 - A) *ᵥ ((1 - A)⁻¹ *ᵥ b) = b := by
      rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hU, Matrix.one_mulVec]
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, ← hxd] at h1
    exact h1
  have hd : ∀ i, (A *ᵥ (x - v)) i ≤ (x - v) i := by
    intro i
    have e1 := congrFun hx i
    have e2 := hv i
    rw [Matrix.mulVec_sub]
    simp only [Pi.sub_apply] at e1 ⊢
    linarith
  intro i
  rw [← hxd]
  have h3 := collatz_sub_nonneg A y (x - v) hA hy hc hd i
  simp only [Pi.sub_apply] at h3
  linarith

/-- `b ≥ 0` implies `(I - A)⁻¹ b ≥ 0`. -/
theorem inv_mulVec_nonneg (A : Matrix (Fin n) (Fin n) ℝ) (y b : Fin n → ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (A *ᵥ y) i < y i)
    (hb : ∀ i, 0 ≤ b i) : ∀ i, 0 ≤ ((1 - A)⁻¹ *ᵥ b) i := by
  have h := neumann_bound A y 0 b hA hy hc (fun i => by simp [Matrix.mulVec_zero, hb i])
  intro i
  simpa using h i

/-- `(I - A)⁻¹ ≥ 0` entrywise. -/
theorem inv_nonneg (A : Matrix (Fin n) (Fin n) ℝ) (y : Fin n → ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (A *ᵥ y) i < y i) :
    ∀ i j, 0 ≤ (1 - A)⁻¹ i j := by
  intro i j
  have h := inv_mulVec_nonneg A y (Pi.single j 1) hA hy hc
    (fun k => by by_cases hk : k = j <;> simp [hk])
  have h' := h i
  rw [Matrix.mulVec_single_one] at h'
  simpa using h'

/-- Spectral radius `< 1`: every complex eigenvalue of `A` has modulus `< 1`. -/
theorem eigen_norm_lt_one (A : Matrix (Fin n) (Fin n) ℝ) (y : Fin n → ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (A *ᵥ y) i < y i)
    (μ : ℂ) (x : Fin n → ℂ) (hx : x ≠ 0)
    (heig : (A.map (fun r : ℝ => (r : ℂ))) *ᵥ x = μ • x) : ‖μ‖ < 1 := by
  by_contra hμ
  push_neg at hμ
  have hd : ∀ i, (A *ᵥ (fun j => -‖x j‖)) i ≤ -‖x i‖ := by
    intro i
    have e := congrFun heig i
    simp only [Matrix.mulVec, dotProduct, Matrix.map_apply, Pi.smul_apply, smul_eq_mul] at e
    have h1 : ‖μ * x i‖ ≤ ∑ j, A i j * ‖x j‖ := by
      rw [← e]
      refine (norm_sum_le _ _).trans (le_of_eq ?_)
      apply Finset.sum_congr rfl
      intro j _
      rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg (hA i j)]
    rw [norm_mul] at h1
    have h2 : ‖x i‖ ≤ ‖μ‖ * ‖x i‖ := le_mul_of_one_le_left (norm_nonneg _) hμ
    simp only [Matrix.mulVec, dotProduct, mul_neg, Finset.sum_neg_distrib]
    linarith
  have h3 := collatz_sub_nonneg A y (fun j => -‖x j‖) hA hy hc hd
  apply hx
  funext i
  have h4 : 0 ≤ -‖x i‖ := h3 i
  exact norm_eq_zero.mp (le_antisymm (by linarith) (norm_nonneg _))

/-- Monotonicity bridge: a Collatz vector for an entrywise upper matrix is a Collatz vector for `A`. -/
theorem collatz_mono (A Aup : Matrix (Fin n) (Fin n) ℝ) (y : Fin n → ℝ)
    (hle : ∀ i j, A i j ≤ Aup i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (Aup *ᵥ y) i < y i) :
    ∀ i, (A *ᵥ y) i < y i := by
  intro i
  have h : (A *ᵥ y) i ≤ (Aup *ᵥ y) i := by
    simp only [Matrix.mulVec, dotProduct]
    exact Finset.sum_le_sum (fun j _ => mul_le_mul_of_nonneg_right (hle i j) (hy j).le)
  linarith [hc i]

/-- **Corollary 3.1 in the form used by the Far Wall certificate (Fin 3).**
`Atrue`: the true bootstrap matrix (0 ≤ Atrue ≤ M); `M`: any real matrix below the rational upper bounds
`Aup` of the cell file (e.g. the matrix of arb ball uppers, so M lies in the arb balls); `y`: the rational
Collatz vector of the cell file; `v = (U, Z, N)` nonnegative with `v ≤ Atrue v + e₁`.
Conclusion: `I - M` invertible, `(I - M)⁻¹ ≥ 0`, `ρ(M) < 1` (eigenvalue form) and `v ≤ (I - M)⁻¹ e₁`. -/
theorem far_wall_corollary_3_1 (Atrue M : Matrix (Fin 3) (Fin 3) ℝ)
    (Aup : Matrix (Fin 3) (Fin 3) ℚ) (y : Fin 3 → ℚ) (v : Fin 3 → ℝ)
    (hA0 : ∀ i j, 0 ≤ Atrue i j) (hAM : ∀ i j, Atrue i j ≤ M i j)
    (hMup : ∀ i j, M i j ≤ (Aup i j : ℝ))
    (hy : ∀ i, 0 < y i)
    (hrow : ∀ i, Aup i 0 * y 0 + Aup i 1 * y 1 + Aup i 2 * y 2 < y i)
    (hv0 : ∀ i, 0 ≤ v i)
    (hv : ∀ i, v i ≤ (Atrue *ᵥ v) i + (Pi.single 0 1 : Fin 3 → ℝ) i) :
    IsUnit (1 - M) ∧ (∀ i j, 0 ≤ (1 - M)⁻¹ i j) ∧
    (∀ (μ : ℂ) (x : Fin 3 → ℂ), x ≠ 0 → (M.map (fun r : ℝ => (r : ℂ))) *ᵥ x = μ • x → ‖μ‖ < 1) ∧
    (∀ i, v i ≤ ((1 - M)⁻¹ *ᵥ (Pi.single 0 1 : Fin 3 → ℝ)) i) := by
  have hyR : ∀ i, 0 < (fun k => (y k : ℝ)) i := fun i => by
    simp only
    exact_mod_cast hy i
  have hcR : ∀ i, ((Aup.map (fun q : ℚ => (q : ℝ))) *ᵥ (fun k => (y k : ℝ))) i
      < (fun k => (y k : ℝ)) i := by
    intro i
    have h' : ((Aup i 0 * y 0 + Aup i 1 * y 1 + Aup i 2 * y 2 : ℚ) : ℝ) < (y i : ℝ) := by
      exact_mod_cast hrow i
    push_cast at h'
    simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_three, Matrix.map_apply]
    linarith
  have hM0 : ∀ i j, 0 ≤ M i j := fun i j => le_trans (hA0 i j) (hAM i j)
  have hcM := collatz_mono M (Aup.map (fun q : ℚ => (q : ℝ))) (fun k => (y k : ℝ))
    (fun i j => by simpa using hMup i j) hyR hcR
  have hvM : ∀ i, v i ≤ (M *ᵥ v) i + (Pi.single 0 1 : Fin 3 → ℝ) i := by
    intro i
    have h1 : (Atrue *ᵥ v) i ≤ (M *ᵥ v) i := by
      simp only [Matrix.mulVec, dotProduct]
      exact Finset.sum_le_sum (fun j _ => mul_le_mul_of_nonneg_right (hAM i j) (hv0 j))
    linarith [hv i]
  exact ⟨isUnit_one_sub M _ hM0 hyR hcM, inv_nonneg M _ hM0 hyR hcM,
    fun μ x hx heig => eigen_norm_lt_one M _ hM0 hyR hcM μ x hx heig,
    neumann_bound M _ v _ hM0 hyR hcM hvM⟩

end FarWallCollatz
