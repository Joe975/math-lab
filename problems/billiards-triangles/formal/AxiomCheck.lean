/-
Axiom audit for the BilliardsFormal certificates.  Run with
`lake env lean AxiomCheck.lean` from `problems/billiards-triangles/formal/`;
every theorem must report exactly [propext, Classical.choice, Quot.sound]
(and never `sorryAx`).  Not part of the library build.
-/
import BilliardsFormal

-- 009: Lemma L1
#print axioms BilliardsFormal.cot_lemma_L1
#print axioms BilliardsFormal.strictAntiOn_id_mul_cot

-- 011, Tier A
#print axioms BilliardsFormal.Laurent.pair_collapse_A
#print axioms BilliardsFormal.Laurent.pair_collapse_B
#print axioms BilliardsFormal.Laurent.half_word_closed_form
#print axioms BilliardsFormal.Laurent.iterate_collapse_A
#print axioms BilliardsFormal.Laurent.iterate_collapse_B
#print axioms BilliardsFormal.Laurent.half_word_letters

-- 011, Tier B
#print axioms BilliardsFormal.Laurent.mu_eq_delta_sq
#print axioms BilliardsFormal.Laurent.tau_parallel_axis
#print axioms BilliardsFormal.Laurent.mu_mul_conj_mu
#print axioms BilliardsFormal.Laurent.half_word_sq_is_translation
#print axioms BilliardsFormal.Laurent.glide_action
#print axioms BilliardsFormal.Laurent.identity_I1
#print axioms BilliardsFormal.Laurent.identity_I2
#print axioms BilliardsFormal.Laurent.identity_I3
#print axioms BilliardsFormal.Laurent.identity_I4
#print axioms BilliardsFormal.Laurent.identity_D1
#print axioms BilliardsFormal.Laurent.identity_D2

-- 011, Tier C
#print axioms BilliardsFormal.Laurent.half_word_closed_form_spec
#print axioms BilliardsFormal.Laurent.identity_I1_spec
#print axioms BilliardsFormal.Laurent.identity_I2_spec
#print axioms BilliardsFormal.Laurent.identity_I3_spec
#print axioms BilliardsFormal.Laurent.identity_I4_spec
#print axioms BilliardsFormal.Laurent.glide_facts_spec
