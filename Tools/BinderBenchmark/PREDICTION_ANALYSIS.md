# Source-bound prediction surface analysis

This increment computes additional structural observations from original predicted
complexes without reading experimental outcomes. It does not yet fit a structural
ranker, run model weights, consume PAE arrays or establish improved binder selection.

## One prediction through the existing CLI

```text
numivivo binder-analyze-prediction PREDICTION.json NEW_BUNDLE
numivivo binder-verify NEW_BUNDLE
```

`VivoBinderPredictionAnalysis.Input` contains `schemaVersion: 1`, an `identity`,
`binderSequence`, an existing `VivoBinderStructureSources.Source`, and `surfacePlan`.
The identity names the candidate, predictor, source sample and target construct.
Use the actual source sample identifier; do not fabricate missing seed or weight
metadata. Original PDB/mmCIF text and its SHA-256 remain inside the source record.
The existing source parser checks declared target-chain sequences, and the analysis
checks the binder sequence. The current input requires exactly one selected binder
chain and one selected target chain, each a complete canonical protein chain.

The strict decoder rejects unknown top-level and structured nested fields, including
experimental labels. Structural calculations have no outcome parameter. Candidate
and model labels remain declarations; a source hash is not an authenticated model
execution receipt or independent evidence of a biological target assignment.

`VivoBinderBundleIO` publishes a `predictionAnalysis` bundle containing the unchanged
`prediction.json`, computed `analysis.json` and the existing receipt. Verification
reparses the original source, rechecks sequences, recomputes geometry and surface,
and compares the complete canonical report. Rehashing a modified feature report
without reconstructing the calculation does not pass. Outputs never overwrite.
The same producing executable is required for strict replay; a rebuilt executable
requires fresh input processing rather than relabelling older evidence.

## Surface and packing definitions

`VivoMolecularInterface.analyzeSurface` extends the existing interface owner and
reuses its occupancy, topology, coordinate, chain and periodicity guards. It uses
explicitly supplied element radii and a probe radius, all in nanometres. The
qualification profile uses C=0.17, N=0.155, O=0.152 and S=0.18 nm with a 0.14 nm probe.
These are declared analysis choices, not native force-field parameter assignment.

The native Shrake-Rupley calculation samples a deterministic golden-angle sphere
around each heavy atom. Spatial cells find possible occluding neighbours. Separate
budgets bound neighbour search, stored neighbour pairs and point tests. It preserves
per-atom isolated and complex solvent-accessible areas, in square nanometres.

For each partner, burial is isolated accessible area minus its area in the complex.
`buriedAreaSumNM2` is the sum of both partner losses. `halfBuriedAreaSumNM2` is exactly
half that sum. Both conventions are named; neither is silently called the other.
The feature dictionary serializes both, while the report also retains partner losses.

The packing observation counts cross-partner pairs for which
`radius_i + radius_j - distance` exceeds the declared overlap threshold. Probe radii
are not added to this overlap calculation. This is a radius-based geometric overlap,
not an interaction energy or a comprehensive chemical clash score.

Computed features include burial, binder buried fraction, overlap count and maximum,
residue contact count, and contacts per summed buried area. A ratio with zero
denominator is unavailable and omitted. Present separated partners retain zero
contacts and zero burial; an absent prediction is not replaced with those values.
Hydrogens, solvent and unselected atoms do not occlude the selected partner surface.
No polar/nonpolar classification, hydrogen-bond typing or unsatisfied polar-group
claim is made by this implementation.

Finite sphere quadrature is not exactly rotation invariant. Analytical sphere tests
bound its error under rotation, and the numerical campaign compares at the declared
resolution. Agreement with another implementation is not a resolution-convergence
study. Higher-resolution and alternative-orientation checks are required before
small feature differences are treated as robust ranking information.

## Multiple samples through the SDK

`VivoBinderPredictionAnalysis.Accumulator` accepts one input at a time. It retains
compact contact/feature records and hashes, rather than every source structure.
It rejects duplicate candidate/predictor/sample/construct identities and changes to
candidate or construct sequence/target identity across predictors. Failed appends
do not change accumulator state.

`finish()` groups by candidate, target construct and predictor. Within each group,
it checks the same sequence and analysis profile, and reports feature minimum,
maximum, mean and population standard deviation with available counts. These are
descriptive sample statistics, not confidence intervals over biological experiments.
Contact-set agreement uses residue offsets in the checked full sequences, avoiding
an arbitrary dependence on different chain names. Pairs with two empty contact sets
are recorded as unavailable agreement, not perfect interface agreement.

Predictors and target constructs remain separate. There is no unanimous-pass filter,
no cross-model pose alignment, and no PAE ingestion yet. The multi-sample path is an
SDK API tested with synthetic fixtures; the public CLI currently analyzes individual
predictions. Persistent batch/series CLI publication remains follow-up work.

Limits include 20 MiB for a prediction request, existing per-source structure limits,
32-4096 surface points per atom, at most 10,000 accumulated predictions, 64 samples
per summarized group and 1,000,000 retained contacts. Exceeding a budget is an explicit
error; the code never silently lowers quadrature or drops candidate inputs.

## Reproduce the published numerical check

```sh
python3 Tools/BinderBenchmark/native_package.py /tmp/numi-prediction-native --test
swift build --package-path /tmp/numi-prediction-native -c release --jobs 2
python3 Tools/BinderBenchmark/fetch_structure_panel.py /tmp/numi-prediction-sources
# The checker needs Biopython 1.86 and NumPy 2.3.5 in its Python environment.
python3 Tools/BinderBenchmark/test_prediction_panel.py -v
python3 Tools/BinderBenchmark/check_prediction_panel.py \
  /tmp/numi-prediction-sources /tmp/numi-prediction-native/.build/release/BinderCLI \
  /tmp/numi-prediction-qualified
```

The checker uses original source bytes, target-construct FASTA and unique sequence
matches. It verifies surface areas with stock Biopython and overlap counts with
independent blocked NumPy pair enumeration. Its reference coordinates and sphere
points have different precision, so atom-level comparisons have explicit tolerances.
It writes the protocol before native analysis, retains each result or failure, and
replays every successfully compared native bundle. It never reads experimental
outcomes from the retained source table. Biopython is a check-time dependency only.

The fixed numerical panel has 18 seed-best Boltz-2 predictions, not the full 270
candidates or the proposed 5,400 multi-predictor samples. See the
[execution record](evidence/2026-09-18-prediction-surface/README.md) for measured scope.

Method reference: Shrake and Rupley, J Mol Biol 79 (1973), 351-371.
Independent implementation: https://biopython.org/docs/latest/api/Bio.PDB.SASA.html
Published source: Anthropic, `claude-protein-binder-design`, CC BY 4.0,
revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`.
