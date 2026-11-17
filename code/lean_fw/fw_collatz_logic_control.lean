import Mathlib

/-! TEETH CONTROL (expected to FAIL): fw_collatz_logic.lean's core lemma with the strict Collatz hypothesis
    A y < y weakened to A y <= y. The weakened statement is FALSE (A = identity, y = 1, d = -1), so no proof
    can pass; the same proof script is used. -/

open Matrix

theorem control_collatz_weak (n : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) (y d : Fin n → ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hy : ∀ i, 0 < y i) (hc : ∀ i, (A *ᵥ y) i ≤ y i)
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
  have hlt := mul_le_mul_of_nonneg_left (hc k) ht0.le
  linarith [hd k]
