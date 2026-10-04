# v0.4: target learning and the shared Codex laboratory

The bounded training failure is diagnosed and corrected without changing the
width-64 architecture. The existing NumiLab plugin now operates editable cards
beside the tissue canvas. These are experimental software results. Biological
promotion remains false; all earlier failures remain unchanged.

## Fixed-architecture diagnosis

Three stimulated GSE92872 source-target groups met the outcome-independent
threshold of 100 treated and 100 explicit control cells: NFKB1, NFATC1 and EGR1.
The attempted eight-group preparation failed that admission check before training
and is retained. All diagnostic expression was previously exposed development
data. Eight arms used seed 271828 and steps 240, 720, 1440. No architecture search.

| Final diagnostic arm | Training response RMSE | Mean between-target SD | Identity accuracy |
|---|---:|---:|---:|
| Legacy Gaussian / SGD | 0.233298 | 0.000482 | 1/3 |
| Mean-only / SGD / no decay | 0.220877 | 0.000571 | 1/3 |
| Mean-only / AdamW / no decay | 0.000044 | 0.194381 | 3/3 |
| Mean-only / AdamW / training scaling | 0.000044 | 0.194381 | 3/3 |
| Mean-only / AdamW / target ID | 0.000044 | 0.194381 | 3/3 |
| Restored Gaussian / AdamW / decay | 0.000315 | 0.194539 | 3/3 |
| Measured diagnostic responses | — | 0.194381 | — |

Removing variance estimation or regularization alone did not fix collapse.
Changing the optimizer under matched mean-only loss did. Scaling was unnecessary.
This establishes a sufficient correction under this budget, not that AdamW is
uniquely optimal. The decoder previously had a feature-count-amplified SGD step
while encoders received the much smaller base step. Instrumentation confirms
small but nonzero legacy encoder updates, distinct target descriptors, the
intended target mask, a 16×515 measured-element denominator, masked-outcome
invariance, finite-difference decoder gradients, actual parameter updates, and
an unchanged frozen prior. Target-input exchanges reproduce exchanged predictions
exactly. Replaying the retained native Gaussian prediction is bit-exact.

The descriptor ridge comparator fits the same small training groups. Native ID
and biological descriptors both fit; this test does not show descriptor superiority.
Exposed capture-six transfer still loses to no change: restored Gaussian RMSE
0.61475 versus 0.55219. Technical captures are not independent biological units.

## Broader campaign and separately reserved result

The unchanged three-study training cohort was rerun with the diagnosed optimizer,
original checkpoint budget and validation selection. Selected descriptor-model
between-group SD is **0.07532**, up from **0.00422**, versus **0.40325** observed
across 81 source/context/target groups. Remaining broad-cohort underfit is explicit.
The target-ID and shuffled-descriptor variants give 0.03906 and 0.01581.

Only after the bounded diagnosis passed, GSM2406677 was reserved. This is the
separately named epistasis preparation in the already exposed GSE90546 study,
not a new independent study. Original barcode, guide, read and UMI identities
restore capture membership; ambiguous identities remain excluded. The author's
notebook maps captures 1/2/3 to tunicamycin/thapsigargin/DMSO and maps PERK/IRE1 to
EIF2AK3/ERN1. Conditions use their own explicit controls. Captures do not become
three biological replicates. No pretrained expression model was used.

The proposed HSPA5/XBP1 objective was rejected because XBP1 is outside the fixed
feature axis. **UPR-HSPA5-v2** is an explicitly new objective, registered before
treated RNA access. It does not replace the failed immediate-early objective.
HSPA5 is detected in 99.77%, 99.91% and 90.74% of 1,317/1,065/983 controls.
All eligible single targets were retained; combination perturbations are outside
this fixed single-target test. Source-bound biological descriptors represent the
unseen ATF6, EIF2AK3 and ERN1 interventions; target-ID lookup receives unknown IDs.

| Condition | Selected | Observed utility | Gain over random | Regret |
|---|---|---:|---:|---:|
| Tunicamycin | EIF2AK3 | 0.68243 | 0.22152 | 0.47239 |
| Thapsigargin | EIF2AK3 | 0.71107 | 0.20633 | 0.49813 |
| DMSO | EIF2AK3 | 0.05348 | -0.12137 | 0.37484 |

Utility is reduction in HSPA5 mean log1p(CPM), not treatment benefit. The training
mean chose no intervention. ATF6 was best observed in every condition. Predictions
and the exact objective were sealed before treated expression was summarized.
There is no independent-unit confidence interval and no reliable-winner claim.
The next discriminating development experiment is condition-specific fitting and
checkpoint diagnostics on admitted development units, followed by a newly reserved
independent preparation; do not tune against and reuse these exposed results.

## Reproduction and evidence

`evidence/v04-recovery/` retains preregistrations, ablations, coverage, results and
local artifact identities. `retained-paths.json` locates native weights and runtime.

```sh
python Tools/VirtualWetLab/diagnose_target_learning.py prepare --inputs INPUTS --output NEW_DIAGNOSIS
python Tools/VirtualWetLab/diagnose_target_learning.py fit --output NEW_DIAGNOSIS --binary NATIVE
python Tools/VirtualWetLab/train_intervention_design.py --inputs ADAM_INPUTS --binary NATIVE --output NEW_CAMPAIGN
python Tools/VirtualWetLab/reserve_target_recovery.py prepare --root RESERVATION --inputs ADAM_INPUTS
python Tools/VirtualWetLab/reserve_target_recovery.py predict --root RESERVATION --campaign NEW_CAMPAIGN --binary NATIVE
# Explicit observation-opening stage, after authorization:
python Tools/VirtualWetLab/reserve_target_recovery.py reveal --root RESERVATION
```

Use retained source files and hashes; the reservation is now exposed, so reproducing
it is replay, never a new independent validation. Legacy plans retain their original
SGD semantics. Old v0.4 weights replay bit-exactly under the updated native code.

The NumiLab plugin guide is `references/virtual-wet-lab.md`: use `numi wet-lab
context`, typed proposals, revision-checked edits, seal, explicit authorized reveal,
and replay. The shared spatial demonstration deliberately retains the original
v0.3 model and its failed biological status; it does not silently promote that
model using this cellular training diagnosis. State PerturbMean remains distinct
from State Transition. Spatial preservation remains unavailable without verified
sections and receiving controls. Safari remains blocked by disabled remote automation.
