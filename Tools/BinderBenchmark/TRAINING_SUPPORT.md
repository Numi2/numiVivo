# Opt-in training-support gate

`VivoBinderTrainingSupport` evaluates an explicitly supplied development policy
before calling the existing fitted benchmark. It shares the legacy target split,
exact-group purge and complete-case admission. No second fitter is introduced.

```text
numivivo binder-evaluate-supported IMPORT_BUNDLE PLAN.json POLICY.json NEW_RESULT
numivivo binder-verify NEW_RESULT
```

Exit 0 means an experimentally eligible fit completed. Exit 2 means the requested
support policy was not met: the complete source-bound bundle is still published,
with counts and reasons, but without any fitted model. Exit 65 means invalid input
or numerical failure, not a successful fallback. When support is declined, retain
the separately produced fixed ranking (`binder-rank`); this command neither refills
nor modifies that selection. Legacy `binder-evaluate` semantics remain unchanged.

`report.json` contains the support report and, only when admitted, `fittedReport`.
`support-policy.json` is retained verbatim alongside original source/import/plan
files. Verification reconstructs cohort admission, group counts, policy decision
and any fitted result from the original source. Rehashing a modified policy does
not make an unchanged dependent report valid.

## Policy meaning

All requirements must be supplied: an identifier/rationale, minimum exclusively
positive and negative training groups, minimum represented training targets,
minimum targets containing both group classes, and maximum feature count.
There is no silent default or universal statistically sufficient sample size.
`plans/support-coverage-example.json` is an illustrative coverage policy, NOT a
calibrated prediction-accuracy rule. It asks for two positive and two negative
training groups, two represented targets with both group classes, and at most
three features. Choose and justify the policy without tuning it to test scores.

Repeated rows in one declared group do not add group support. A group containing
both outcomes is recorded as mixed and contributes to neither exclusive group
class. Missing predictors, nonbinary training outcomes and groups shared with ANY
held-out candidate are excluded before counting. Only training-side counts and
exclusions affect admission; held-out outcomes do not authorize a fit.

Groups remain caller-declared. The existing exact-sequence grouping does not
establish design-family, near-homology or target independence. Admission alone
also does not establish coefficient stability, prediction calibration, sufficient
power, or reliable generalization. Those remain separate required experiments.

## Portable and published-table checks

The native test set covers accepted/declined fits, single-class data, repeated
and mixed groups, target coverage, missing features, group leakage, test-label
isolation, unchanged accepted legacy reports and policy-tamper replay rejection.

After running the existing six-fold published campaign with the same executable:

```sh
python3 Tools/BinderBenchmark/check_support.py CAMPAIGN_DIR NATIVE_CLI \
  Tools/BinderBenchmark/plans/support-coverage-example.json NEW_CHECK_DIR
```

This independently reconstructs training candidates/group counts in Python,
checks the support decisions, verifies every resulting native bundle, and requires
accepted fitted reports to equal the corresponding original legacy reports.
The six source bundles are retained in the check directory. It is an example-policy
execution check on reused development data, not a new model or biological uplift.
