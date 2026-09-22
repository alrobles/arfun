module

public import Mathlib.Analysis.SpecialFunctions.Log.Basic
public import Mathlib.Tactic

@[expose] public section

namespace ArfunFormal.Likelihood

noncomputable def logScore (n : ℕ) (logNumerator mass : ℝ) : ℝ :=
  logNumerator - (n : ℝ) * Real.log mass

theorem normalizer_cancels (n : ℕ) (a : ℝ) {c z : ℝ}
    (hc : 0 < c) (hz : 0 < z) :
    logScore n (a + n * Real.log c) (c * z) = logScore n a z := by
  simp only [logScore, Real.log_mul (ne_of_gt hc) (ne_of_gt hz)]
  ring

theorem finite_normalizer_cancels {ι : Type*} [Fintype ι]
    (f : ι → ℝ) {c z : ℝ} (hf : ∀ i, 0 < f i) (hc : 0 < c) (hz : 0 < z) :
    logScore (Fintype.card ι) (∑ i, Real.log (c * f i)) (c * z) =
      logScore (Fintype.card ι) (∑ i, Real.log (f i)) z := by
  have hs : (∑ i, Real.log (c * f i)) =
      (∑ i, Real.log (f i)) + (Fintype.card ι : ℝ) * Real.log c := by
    calc
      _ = ∑ i, (Real.log c + Real.log (f i)) :=
        Finset.sum_congr rfl (fun i _ => Real.log_mul (ne_of_gt hc) (ne_of_gt (hf i)))
      _ = _ := by simp [Finset.sum_add_distrib, add_comm]
  rw [hs]
  exact normalizer_cancels _ _ hc hz

theorem background_scale {ι : Type*} [Fintype ι]
    (w f : ι → ℝ) (c : ℝ) :
    (∑ i, w i * (c * f i)) = c * (∑ i, w i * f i) := by
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem background_pos {ι : Type*} [Fintype ι] [Nonempty ι]
    (w f : ι → ℝ) (hw : ∀ i, 0 < w i) (hf : ∀ i, 0 < f i) :
    0 < ∑ i, w i * f i := by
  exact Finset.sum_pos (fun i _ => mul_pos (hw i) (hf i)) Finset.univ_nonempty

theorem quadrature_constant (n : ℕ) (a : ℝ) {area m z : ℝ}
    (ha : 0 < area) (hm : 0 < m) (hz : 0 < z) :
    logScore n a ((area / m) * z) =
      logScore n a z - (n : ℝ) * Real.log (area / m) := by
  simp only [logScore, Real.log_mul (ne_of_gt (div_pos ha hm)) (ne_of_gt hz)]
  ring

theorem common_constant_cancels (n : ℕ) (a b : ℝ) {c z₁ z₂ : ℝ}
    (hc : 0 < c) (h₁ : 0 < z₁) (h₂ : 0 < z₂) :
    logScore n a (c * z₁) - logScore n b (c * z₂) =
      logScore n a z₁ - logScore n b z₂ := by
  simp only [logScore, Real.log_mul (ne_of_gt hc) (ne_of_gt h₁),
    Real.log_mul (ne_of_gt hc) (ne_of_gt h₂)]
  ring

theorem log_sum_exp_shift {ι : Type*} [Fintype ι] [Nonempty ι]
    (a : ι → ℝ) (c : ℝ) :
    Real.log (∑ i, Real.exp (c + a i)) = c + Real.log (∑ i, Real.exp (a i)) := by
  have hz : 0 < ∑ i, Real.exp (a i) :=
    Finset.sum_pos (fun i _ => Real.exp_pos _) Finset.univ_nonempty
  simp_rw [Real.exp_add]
  rw [← Finset.mul_sum, Real.log_mul (Real.exp_ne_zero c) (ne_of_gt hz), Real.log_exp]

end ArfunFormal.Likelihood
