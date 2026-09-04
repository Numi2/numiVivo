# Surrogate authority boundary

## Principle

A learned or reduced model can accelerate NumiVivo only when its domain, uncertainty, and refresh policy are explicit. A surrogate never becomes the authority for biological constraints, safety monitors, artifact identity, or transactional commit.

## Contract identity

Every surrogate contract binds:

- model fingerprint;
- training-data fingerprint;
- mechanism-pack fingerprint;
- optional host-context fingerprint;
- input feature definitions;
- output feature definitions;
- uncertainty representation;
- acceptance thresholds;
- mandatory authoritative-refresh interval.

A backend whose model fingerprint or tensor width does not match the contract is rejected before prediction.

## Feature contract

Each input and output feature declares:

- identifier;
- unit;
- lower and upper domain bound;
- normalization offset;
- normalization scale.

Normalization parameters are part of the fingerprinted contract. A runtime cannot substitute training-time normalization inferred from model metadata.

## Extrapolation

For each input feature, normalized extrapolation is computed relative to its declared range. The gate takes the maximum over the complete batch.

An input outside the allowed extrapolation radius causes immediate authoritative fallback. The backend is not invoked when the input shape or domain is already invalid.

## Uncertainty

Supported uncertainty forms are:

- none;
- predicted standard deviation;
- ensemble variance;
- conformal radius.

The contract defines a maximum normalized uncertainty. Backend output may contain one uncertainty value per batch row or one per output element. All values must be finite and non-negative.

A contract declaring no uncertainty must use a zero threshold. This prevents an unreported uncertainty estimate from being implied by a positive allowance.

## Output validation

After prediction, the gate verifies:

- returned model fingerprint;
- batch size;
- output width;
- output value count;
- uncertainty value count;
- finite outputs;
- output feature bounds.

An out-of-domain output is rejected even when the reported uncertainty is low.

## Mandatory authoritative refresh

A surrogate can be accepted only for a bounded number of consecutive steps. At the configured interval, the gate requires an authoritative evaluation and resets the acceptance counter.

The authoritative result can be used to:

- continue simulation;
- measure surrogate error;
- update online diagnostics;
- invalidate the contract;
- initiate recalibration through a new artifact lineage.

The current contract does not permit unrecorded online weight updates.

## Backend implementations

### Affine reference backend

The affine backend evaluates a fixed matrix and bias with explicit output uncertainty. It provides a deterministic reference path for reduced models and backend comparisons.

### Core ML backend

The Core ML backend loads a compiled model with an explicit `MLModelConfiguration`, input feature name, output feature name, optional uncertainty feature, and declared tensor widths.

The backend returns predictions only. The authority gate performs all domain and uncertainty acceptance checks. Core ML compute placement does not change semantic authority.

### Future MLX backend

An MLX backend may be added for research workflows and trainable local surrogates. It must implement the same protocol and artifact identity. An MLX model is not permitted to mutate production state directly or bypass the authority gate.

## Runtime integration

The authoritative sequence is:

1. Construct a typed input batch from committed state.
2. Validate shape and input domain.
3. Enforce the refresh interval.
4. Run the backend.
5. Validate prediction identity, uncertainty, and output bounds.
6. Either reject and invoke the authoritative mechanism, or use the prediction as a candidate numerical value.
7. Execute normal runtime constraints and monitors.
8. Commit or reject through the normal transaction.

An accepted surrogate result remains candidate state until the owning transaction commits.

## Failure reasons

The gate reports one reason:

- accepted;
- authoritative interval reached;
- model fingerprint mismatch;
- input shape mismatch;
- input outside declared domain;
- excessive extrapolation;
- excessive uncertainty;
- invalid output;
- output outside bounds;
- backend failure.

The reason is included in execution diagnostics and may be aggregated over a campaign.

## Calibration and replacement

A surrogate contract is immutable. Recalibration creates:

- a new training-data fingerprint;
- a new model fingerprint;
- a new contract fingerprint;
- a new validation report;
- a new deployment decision.

Replacing model bytes while retaining the old fingerprint is prohibited.

## Safety restrictions

A surrogate cannot:

- assert that a monitor is safe without authoritative evaluation;
- change a declared termination response;
- extend its own domain;
- lower its own uncertainty;
- suppress a non-finite authoritative state;
- alter a random-stream namespace;
- discard delayed events;
- change coupling ownership;
- sign or rewrite provenance records;
- promote simulation output into biological validation.

## Current limitations

The source includes a Core ML backend and a deterministic affine backend but has not received integration compilation or numerical validation. Model conversion, feature naming, shape conventions, and compute-unit selection must be checked in the designated Apple environment before production use.
