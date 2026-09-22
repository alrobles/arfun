module

public import Mathlib.Analysis.SpecialFunctions.Artanh
public import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
public import Mathlib.LinearAlgebra.Matrix.PosDef
public import Mathlib.Tactic

@[expose] public section

namespace ArfunFormal.Geometry

noncomputable def covariance (s₁ s₂ ρ : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![s₁ ^ 2, ρ * s₁ * s₂; ρ * s₁ * s₂, s₂ ^ 2]

noncomputable def cholesky (s₁ s₂ ρ : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![s₁, 0; ρ * s₂, s₂ * Real.sqrt (1 - ρ ^ 2)]

noncomputable def quadratic (s₁ s₂ ρ x y : ℝ) : ℝ :=
  s₁ ^ 2 * x ^ 2 + 2 * ρ * s₁ * s₂ * x * y + s₂ ^ 2 * y ^ 2

theorem correlation_margin_pos {ρ : ℝ} (hρ : |ρ| < 1) : 0 < 1 - ρ ^ 2 := by
  obtain ⟨hl, hu⟩ := abs_lt.mp hρ
  have h : 0 < (1 - ρ) * (1 + ρ) := mul_pos (by linarith) (by linarith)
  nlinarith

theorem covariance_det (s₁ s₂ ρ : ℝ) :
    (covariance s₁ s₂ ρ).det = s₁ ^ 2 * s₂ ^ 2 * (1 - ρ ^ 2) := by
  simp [covariance, Matrix.det_fin_two]
  ring

theorem covariance_det_pos {s₁ s₂ ρ : ℝ}
    (h₁ : 0 < s₁) (h₂ : 0 < s₂) (hρ : |ρ| < 1) :
    0 < (covariance s₁ s₂ ρ).det := by
  rw [covariance_det]
  exact mul_pos (mul_pos (sq_pos_of_pos h₁) (sq_pos_of_pos h₂))
    (correlation_margin_pos hρ)

theorem quadratic_completed_square (s₁ s₂ ρ x y : ℝ) :
    quadratic s₁ s₂ ρ x y =
      (s₁ * x + ρ * s₂ * y) ^ 2 + (1 - ρ ^ 2) * (s₂ * y) ^ 2 := by
  unfold quadratic
  ring

theorem quadratic_pos {s₁ s₂ ρ x y : ℝ}
    (h₁ : 0 < s₁) (h₂ : 0 < s₂) (hρ : |ρ| < 1) (hxy : x ≠ 0 ∨ y ≠ 0) :
    0 < quadratic s₁ s₂ ρ x y := by
  rw [quadratic_completed_square]
  by_cases hy : y = 0
  · have hx : x ≠ 0 := hxy.resolve_right (not_not.mpr hy)
    simpa [hy] using sq_pos_of_ne_zero (mul_ne_zero (ne_of_gt h₁) hx)
  · exact add_pos_of_nonneg_of_pos (sq_nonneg _)
      (mul_pos (correlation_margin_pos hρ)
        (sq_pos_of_ne_zero (mul_ne_zero (ne_of_gt h₂) hy)))

theorem covariance_posDef {s₁ s₂ ρ : ℝ}
    (h₁ : 0 < s₁) (h₂ : 0 < s₂) (hρ : |ρ| < 1) :
    (covariance s₁ s₂ ρ).PosDef := by
  apply Matrix.posDef_iff_dotProduct_mulVec.mpr
  constructor
  · ext i j
    fin_cases i <;> fin_cases j <;> simp [covariance, Matrix.conjTranspose_apply]
  · intro x hx
    have hxy : x 0 ≠ 0 ∨ x 1 ≠ 0 := by
      by_contra h
      have hzero := not_or.mp h
      apply hx
      funext i
      fin_cases i <;> simp_all
    have heq : dotProduct (star x) ((covariance s₁ s₂ ρ).mulVec x) =
        quadratic s₁ s₂ ρ (x 0) (x 1) := by
      simp [covariance, quadratic, dotProduct, Matrix.mulVec, Fin.sum_univ_two]
      ring
    rw [heq]
    exact quadratic_pos h₁ h₂ hρ hxy

theorem cholesky_factorization {s₁ s₂ ρ : ℝ} (hρ : |ρ| < 1) :
    cholesky s₁ s₂ ρ * (cholesky s₁ s₂ ρ).transpose = covariance s₁ s₂ ρ := by
  have hsqrt := Real.sq_sqrt (le_of_lt (correlation_margin_pos hρ))
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [cholesky, covariance, Matrix.mul_apply, Fin.sum_univ_two]
  all_goals ring_nf
  nlinarith [congrArg (fun x : ℝ => s₂ ^ 2 * x) hsqrt]

theorem transformed_covariance_posDef (η₁ η₂ ζ : ℝ) :
    (covariance (Real.exp η₁) (Real.exp η₂) (Real.tanh ζ)).PosDef := by
  apply covariance_posDef (Real.exp_pos _) (Real.exp_pos _)
  exact abs_lt.mpr ⟨Real.neg_one_lt_tanh ζ, Real.tanh_lt_one ζ⟩

end ArfunFormal.Geometry
