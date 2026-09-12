# Protected integration design admission

Ridge and MNN accept optional `protectedSampleGroups`, a complete observed-sample
map. Covariate levels must connect through shared joint condition/group strata.
The guard runs before correction. Missing/extra samples, invalid labels, excessive
mapping bytes and disconnected designs reject. Omission preserves prior behavior.
See [the interface](../INTEGRATION.md#optional-protected-sample-groups).

## Evidence and boundaries

- Actual native owners run on all 13,863 retained Hagai PCA rows (20 components).
  Old/new unprotected results are byte-identical for both methods. With a constant
  meaningful protected group, all result fields except the added option match.
  This reruns integration, not source H5AD import or PCA fitting.
- Source-bound HIRISA metadata represent all 131 deposited libraries. Donor design
  admission passes; batch correction rejects preparation confounding in both
  owners. These are one-row-per-library zero-coordinate rejection sentinels, not
  a new 1.6-million-cell integration or preservation measurement.
- Missing, extra, empty, unreported, oversized and unknown maps/options reject.
  A crossed design whose marginal graphs connect but joint strata do not rejects.
- Fresh synthetic H5AD/PCA CLI bundles publish and replay for both methods with
  and without the option. Confounded maps reject without a published destination.
  These are software lifecycle checks, not experimental biology.
- Scoped actual owner and single-cell router compilation passed. The full app and
  package were not built. The compiler inputs and binary hashes are retained.

This guard does not protect an undeclared variable or mixed cell-level labels
within a sample. It does not constrain the numerical correction, establish
preservation, support simultaneous nuisance covariates, or predict unseen donors.
The historical HIRISA preservation failures remain failures.

## Reproduction

`evidence.tar.gz` retains inputs, compressed numerical outputs, source provenance,
the Swift driver, checks, fixture bundles and failures. Its manifest hashes every
member. `verify.py` independently checks retained numerical equivalence and hash
bindings without rerunning a model. `run_checks.py` reruns native owners; its
`--reuse-retained` switch is only for inspecting retained evidence, not fresh
execution qualification. Build the driver against the documented old/new scoped
modules; local binary paths in the scripts must be adjusted on another machine.

The first run stopped with ENOSPC, and a later assertion expected lowercase
`unknown` where the decoder emitted `Unknown`. Both failed logs are retained.
Output now streams directly into gzip. Completed numerical outputs were reused
for the final report; no new runtime timing claim is made. The retained old ridge
output predates the driver's output-streaming change, which did not alter model
calls or serialization. Its payload equals the new output exactly.
