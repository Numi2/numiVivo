# What the fixed-architecture campaign learned

All numbers below are **exposed development**, not independent biological
validation. The prior three-condition EIF2AK3/ATF6 result is unchanged and was not
used in this campaign. No new independent evaluation was reserved.

The registered schedule extension executed 13 arms through 5,760 updates, keeping
all five checkpoints and source × condition × target diagnostics. The original
three-target qualification remains historical evidence; it was not reclassified.

* Individually, the two K562 studies fit distinct responses: mean-only final
  equal-group response RMSE 0.0150 (GSE90063) and 0.0158 (GSE90546), with every
  training target distinguished. Distributional training fits worse, 0.0573 and
  0.0361, without losing target identity on those studies.
* GSE92872 alone remains substantially harder: mean-only RMSE 0.2671 versus
  distributional 0.3270. Its 64 source/condition/target groups account for most
  cohort difficulty. This failure exists within a study, so cross-study response
  conflict cannot be its sole explanation.
* Across 81 groups, source balancing at learning rate .001 worsens mean-only fit
  (0.3173 versus group-balanced 0.2803). A smaller .0003 learning rate improves
  source-balanced fit to 0.2648. Restored distributional fits are worse.
* Supported condition metadata improves fit: source-balanced mean-only RMSE
  0.3173 versus RNA-only 0.3440 and shuffled metadata 0.3341. In matched condition
  probes, descriptor swaps move predicted responses much more than control-RNA
  swaps (RMS 0.2054 versus 0.0050). The RNA-only response head is effectively
  invariant in this probe. Condition is present in the control profile, but its
  presence is not evidence the learned response uses it effectively.
* A preregistered three-arm follow-up keeps the architecture fixed. Target-ID
  conditioning reaches RMSE 0.1174 with 69/81 nearest-target identities correct
  and 89.7% direction agreement outside the fixed 0.1 deadband. This is known-
  target fit, **not unseen-target transfer**. Group-mean descriptor fitting helps
  modestly (0.2423); a larger .003 learning rate collapses toward an average
  response (0.3748). Thus bag noise contributes, but neither noise removal nor a
  larger learning rate resolves descriptor-conditioned underfitting.
* The sole validation capture favors early, flatter checkpoints while training
  fit continues improving. The preregistered selector chooses the .0003 arm at
  step 240, validation RMSE 0.70575. It is a stimulated technical capture from
  GSE92872, not validation for either K562 study or independent biological units.
  The target-ID model's better fit worsens that validation result. Selecting a
  later checkpoint for its response variation would hide this failure.

These comparisons discriminate optimization/conditioning limitations from a
universal inability to learn intervention effects. They do not establish useful
transfer, an information-gain policy, or a validated intervention ranking. The
next independent evaluation should await a frozen model that resolves these
limitations on development data; no success threshold was lowered here.

Native historical inference remains bit-exact. New schedule admission tests
reject unregistered, unordered, duplicate, empty, excessive and over-budget plans.
Batch-aligned chunked inference matches a full 648-row native prediction exactly.
An interrupted spatial-scoring attempt and a concurrent source-directory loss
are retained; source restoration matched all original cohort/input/runtime hashes.

## Separately retrained spatial development

Three matched spatial configurations used the corrected Adam optimizer at fixed
architecture, with checkpoints 240/1440/4320. Validation chose step 240 in each.
Exposed specimen equal-target population RMSE is 0.76276 with neighborhood context,
0.76398 without it, 0.76141 with shuffled neighborhoods, and 0.73608 for no change.
The small neighborhood-versus-no-neighborhood difference does not survive the
shuffle control or beat no change. Spatial benefit remains unearned.

The separately retrained receiving-population configurations score 0.92838 for
receiver-conditioned pretraining, 0.92991 for spatial-only receiver training,
0.78407 for the anchor-reference comparison and 0.91833 for no change (mean group
RMSE across the same 26 exposed groups). This does not support the hypothesis
that receiver-reference substitution fixes the measured failure. Verified section
edges and receiving-population preservation controls remain unavailable.

The exact corrected spatial artifact is `spatial-response-v04-corrected`, runtime
SHA-256 `bca3896e2145deee72db2ad83981eddf81e0c369a68aff2073628437e60375de`.
It remains MODEL INFERENCE. Its three native weight artifacts and plan are bound
on each new shared comparison card; v0.3 remains separately selectable.
