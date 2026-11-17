import Mathlib

/-! Far Wall tail (fw_tail2, K_T = 2000, M = 3.5): Collatz rows A_up * y < y and y > 0 for the UNIFORM
    bulk matrix valid for every K >= 2000 (PROOF-B4-v4-AMEND2-TAIL sec. 2). Entries are decimal upper bounds
    strictly above the ball-arithmetic (arb) uppers of fw_tail2.bulk_tail(2000) (asserted in lean_fw/lean_tail.py).
    Only these rational inequalities are kernel-checked here. -/

theorem fw_tail_row0 : (219429/1000000 : ℚ) * (3063/10000 : ℚ) + (491031/5000000 : ℚ) * (1/1 : ℚ) + (91161/250000000 : ℚ) * (5587/10000 : ℚ) < (3063/10000 : ℚ) := by norm_num
theorem fw_tail_row1 : (481041/500000 : ℚ) * (3063/10000 : ℚ) + (109831/1000000 : ℚ) * (1/1 : ℚ) + (243731/1000000 : ℚ) * (5587/10000 : ℚ) < (1/1 : ℚ) := by norm_num
theorem fw_tail_row2 : (866933/1000000 : ℚ) * (3063/10000 : ℚ) + (181919/5000000 : ℚ) * (1/1 : ℚ) + (330783/1000000000 : ℚ) * (5587/10000 : ℚ) < (5587/10000 : ℚ) := by norm_num
theorem fw_tail_ypos : (0 : ℚ) < (3063/10000 : ℚ) ∧ (0 : ℚ) < (1/1 : ℚ) ∧ (0 : ℚ) < (5587/10000 : ℚ) := by norm_num
