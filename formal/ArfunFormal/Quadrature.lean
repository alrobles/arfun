module

public import ArfunFormal.Likelihood

@[expose] public section

namespace ArfunFormal.Quadrature

open ArfunFormal.Likelihood

theorem log_relative_error_bound {z zhat ε : ℝ}
    (hz : 0 < z) (hzhat : 0 < zhat) (hε : 0 ≤ ε) (hε₁ : ε < 1)
    (herr : |zhat / z - 1| ≤ ε) :
    |Real.log zhat - Real.log z| ≤ -Real.log (1 - ε) := by
  obtain ⟨hlo, hhi⟩ := abs_le.mp herr
  have hminus : 0 < 1 - ε := by linarith
  have hplus : 0 < 1 + ε := by linarith
  have hr : 0 < zhat / z := div_pos hzhat hz
  have hlower : Real.log (1 - ε) ≤ Real.log (zhat / z) :=
    Real.log_le_log hminus (by linarith)
  have hupper : Real.log (zhat / z) ≤ Real.log (1 + ε) :=
    Real.log_le_log hr (by linarith)
  have hprod : (1 + ε) * (1 - ε) ≤ 1 := by nlinarith [sq_nonneg ε]
  have hlogs : Real.log (1 + ε) + Real.log (1 - ε) ≤ 0 := by
    rw [← Real.log_mul (ne_of_gt hplus) (ne_of_gt hminus)]
    exact Real.log_nonpos (le_of_lt (mul_pos hplus hminus)) hprod
  rw [← Real.log_div (ne_of_gt hzhat) (ne_of_gt hz)]
  apply abs_le.mpr
  constructor <;> linarith

theorem logScore_error_bound (n : ℕ) (a : ℝ) {z zhat ε : ℝ}
    (hz : 0 < z) (hzhat : 0 < zhat) (hε : 0 ≤ ε) (hε₁ : ε < 1)
    (herr : |zhat / z - 1| ≤ ε) :
    |logScore n a zhat - logScore n a z| ≤
      (n : ℝ) * (-Real.log (1 - ε)) := by
  have hd : logScore n a zhat - logScore n a z =
      -(n : ℝ) * (Real.log zhat - Real.log z) := by
    unfold logScore
    ring
  rw [hd, abs_mul, abs_neg, abs_of_nonneg (Nat.cast_nonneg n)]
  exact mul_le_mul_of_nonneg_left
    (log_relative_error_bound hz hzhat hε hε₁ herr) (Nat.cast_nonneg n)

theorem bound_nonneg (n : ℕ) {ε : ℝ} (hε : 0 ≤ ε) (hε₁ : ε < 1) :
    0 ≤ (n : ℝ) * (-Real.log (1 - ε)) := by
  apply mul_nonneg (Nat.cast_nonneg n)
  exact neg_nonneg.mpr (Real.log_nonpos (by linarith) (by linarith))

end ArfunFormal.Quadrature
