module

public import ArfunFormal.Geometry
public import Mathlib.Algebra.Order.Star.Real

@[expose] public section

namespace ArfunFormal.Evolution

noncomputable def branchCovariance {S B : Type*} [Fintype B] [DecidableEq B]
    (paths : Matrix S B ℝ) (lengths : B → ℝ) : Matrix S S ℝ :=
  paths * Matrix.diagonal lengths * paths.transpose

noncomputable def ouCovariance (α q ti tj shared : ℝ) : ℝ :=
  q / (2 * α) * Real.exp (-α * (ti + tj - 2 * shared)) *
    (1 - Real.exp (-2 * α * shared))

theorem branchCovariance_posSemidef {S B : Type*} [Fintype S] [Fintype B] [DecidableEq B]
    (paths : Matrix S B ℝ) (lengths : B → ℝ) (hl : ∀ b, 0 ≤ lengths b) :
    (branchCovariance paths lengths).PosSemidef := by
  classical
  have hd : (Matrix.diagonal lengths).PosSemidef := Matrix.posSemidef_diagonal_iff.mpr hl
  simpa [branchCovariance] using hd.mul_mul_conjTranspose_same paths

theorem noncentered_covariance {S : Type*} [Fintype S]
    (L : Matrix S S ℝ) {v : ℝ} (hv : 0 ≤ v) :
    (Real.sqrt v • L) * (Real.sqrt v • L).transpose = v • (L * L.transpose) := by
  rw [Matrix.transpose_smul, smul_mul_assoc, mul_smul_comm, smul_smul,
    Real.mul_self_sqrt hv]

theorem bm_contrast_variance (v ti tj shared : ℝ) :
    v * ti + v * tj - 2 * (v * shared) = v * (ti + tj - 2 * shared) := by
  ring

theorem bm_time_rescaling {c : ℝ} (hc : 0 < c) (v t : ℝ) :
    (v / c) * (c * t) = v * t := by
  field_simp [ne_of_gt hc]

theorem ouCovariance_symm (α q ti tj shared : ℝ) :
    ouCovariance α q ti tj shared = ouCovariance α q tj ti shared := by
  unfold ouCovariance
  rw [add_comm ti tj]

theorem ouCovariance_diagonal (α q t : ℝ) :
    ouCovariance α q t t t = q / (2 * α) * (1 - Real.exp (-2 * α * t)) := by
  unfold ouCovariance
  have ht : t + t - 2 * t = 0 := by ring
  rw [ht]
  simp

theorem ouCovariance_nonneg {α q shared : ℝ} (ti tj : ℝ)
    (hα : 0 < α) (hq : 0 ≤ q) (ht : 0 ≤ shared) :
    0 ≤ ouCovariance α q ti tj shared := by
  have he : Real.exp (-2 * α * shared) ≤ 1 :=
    Real.exp_le_one_iff.mpr (by nlinarith [mul_nonneg (le_of_lt hα) ht])
  exact mul_nonneg
    (mul_nonneg (div_nonneg hq (by linarith)) (le_of_lt (Real.exp_pos _)))
    (sub_nonneg.mpr he)

theorem stationary_root_identity (α q ti tj shared : ℝ) :
    ouCovariance α q ti tj shared + q / (2 * α) * Real.exp (-α * (ti + tj)) =
      q / (2 * α) * Real.exp (-α * (ti + tj - 2 * shared)) := by
  have he : Real.exp (-α * (ti + tj - 2 * shared)) * Real.exp (-2 * α * shared) =
      Real.exp (-α * (ti + tj)) := by
    rw [← Real.exp_add]
    congr 1
    ring
  unfold ouCovariance
  rw [← he]
  ring

end ArfunFormal.Evolution
