# Fixed-architecture cohort recovery (development)

The native spatial learner now accepts an explicitly registered `trainingBudget`
and strictly increasing `steps`, ending at that budget (1–100,000 updates, at most
32 checkpoints). Without a budget, only the historical 240/720/1440 schedule is
accepted. Absent `sampling` retains the exact group-sampling random-number path;
`sampling: source` first samples a source, then a group within that source. Source
identity must be supplied for every row, and a group cannot cross sources.

`train_cohort_recovery.py register` freezes the input/runtime/owner identities and
13-arm budget before `run`. It retains every checkpoint and training/validation
source × condition × target diagnostic. Validation exists only for the GSE92872
stimulated technical capture; missing study validation is explicitly unavailable.
No GSM2406677 outcomes or new independent evaluation enter this campaign.

The initial cohort result selects an early checkpoint while later training fit
improves. `diagnose_cohort_fit.py` preregisters three additional development-only
comparisons at the same architecture. `probe_cohort_condition.py` retains matched
control-profile and metadata swaps. `retrain_spatial_development.py` separately
requalifies the exposed spatial/receiver configurations with Adam. These are
engineering and exposed-development results, not biological promotion.

Local evidence: `/Users/home/virtual-wet-lab-v04-cohort-20261004`.
