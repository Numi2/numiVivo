# Native structural observations for binder evaluation

This route measures static heavy-atom geometry on `VivoMolecularStructure`. It
reuses `VivoStructureValidator` and the binder evaluator. It does not run a folding
model, repair a protein, assign charges, simulate dynamics or predict binding affinity.

## Owners and execution

- `Structure/VivoMolecularInterface.swift`: bounded nonperiodic interface geometry.
- `Binder/VivoBinderStructuralFeatures.swift`: candidate/sequence checks and derived features.
- `Binder/VivoBinderBundleIO.swift`: immutable input/result bundles and complete replay.

```text
numivivo binder-evaluate-structures IMPORT_BUNDLE PLAN.json STRUCTURES.json NEW_RESULT
numivivo binder-verify NEW_RESULT
```

`IMPORT_BUNDLE` must have been produced with the same executable. A structural
input has `schemaVersion: 1` and an `observations` array. Every observation contains
`candidateID`, `target`, `sourceLabel`, a complete native `structure`, and an
`interfacePlan` with these fields:

```json
{
  "schemaVersion": 1,
  "binderChains": ["B"],
  "targetChains": ["T"],
  "conformerID": "prediction-0",
  "contactDistanceNM": 0.45,
  "shortDistanceNM": 0.2,
  "maximumPairEvaluations": 25000000
}
```

Chain and conformer names must refer to the actual structure; the example names
are not inferred defaults. Native conformer coordinates are in nanometers. Use
existing molecular import code to convert PDB/mmCIF coordinates rather than
inserting Angstrom-valued coordinates into a nanometer container. This command
accepts native structures, not raw PDB/mmCIF files. Retain the original prediction
files separately: this version reconstructs geometry from the embedded structure,
not from an authenticated external model run or original PDB/mmCIF bytes.

## Identity, supported chemistry and exclusions

The candidate must exist in the experimental import and the target string must
match. Exactly one binder chain is required and its canonical residue sequence
must equal the source row's `sequence`. Standard protein residues must contain
all expected heavy-atom names and matching elements; unsupported residues, missing
atoms and extra nonstandard heavy atoms are rejected. `OXT` is permitted. No missing
atoms, protonation states or parameters are invented.

Target chains undergo the same canonical residue checks. Selected target sequences
must agree within one target cohort, but the biological target assignment remains
caller-declared. `sourceLabel` records the caller's model-origin description; it is
not a verified model/weights/execution identity. Sequence agreement is not a proof
of structural correctness or physical preparation.

Selected chains must be disjoint, known and contain heavy atoms. Periodic inputs,
selected heavy-atom alternative locations, partial/zero occupancies, duplicate
residue/atom identities and declared covalent bonds between partners are rejected.
Hydrogens and unselected atoms do not contribute to the measurements. A declared
absence of covalent bonds does not establish that the chemical topology is complete.

All observations use the same contact and short-distance cutoffs. Candidate-specific
cutoffs are rejected. Missing structures are recorded as missing features, never
zero contacts. A present structure with separated partners legitimately has zero
contacts and remains distinguishable from an absent structure.

## Features and units

All names below have the prefix `numi.interface.`:

| Feature | Definition |
| --- | --- |
| `contactAtomPairs` | Selected cross-partner heavy-atom pairs at distance <= contactDistanceNM. |
| `contactResiduePairs` | Unique cross-partner residue pairs with at least one contact. |
| `binderContactFraction` | Binder heavy atoms contacting the target / selected binder heavy atoms. |
| `targetContactFraction` | Target heavy atoms contacting the binder / selected target heavy atoms. |
| `shortDistanceFraction` | Cross-partner heavy-atom pairs at distance < shortDistanceNM / all evaluated cross-partner pairs. |
| `minimumDistanceNM` | Minimum cross-partner heavy-atom distance in nanometers. |
| `centroidDistanceNM` | Distance between unweighted geometric heavy-atom centroids in nanometers. |
| `binderGeometricRgNM` | Root mean squared heavy-atom distance from the binder's unweighted geometric centroid, in nanometers. |

The report also retains exact heavy/contact atom indices, closest residue-pair
distances, evaluated pair count, structure SHA-256, method and plan. Calculations
use Double arithmetic. Cutoffs are declared analysis choices, not calibrated
biological thresholds. Short-distance counts are not force-field clash scores;
centroids and radius of gyration are not mass weighted. No hydrogen-bond, surface
area, electrostatic energy, binding-energy or affinity claim is made.

## Matched comparison

The plan must combine one or more allowed in-silico score features with one or more
of the geometry features above. Keep the baseline an allowed in-silico feature.
For example, the existing three ipSAE features can be compared with those same
features plus `numi.interface.contactAtomPairs` and
`numi.interface.minimumDistanceNM`. Selecting these features is an experiment,
not a claim that they are predictive.

`report.json` contains the combined model. `score-only-report.json` contains a
separately fitted control with the same score features and no geometry. Both use
identical training and held-out candidate IDs, the same target split and ridge
penalty, and training-only scaling. Candidates missing any required combined-model
feature are excluded from both fits/evaluations. Original test rows are retained
for sequence-group purging even when they cannot be scored.

`geometry.json` retains derived observations and unavailable candidate IDs.
`structures.json`, the original CSV/import, plan and reports are bound in the
receipt. `binder-verify` reimports source, rechecks sequences, recomputes geometry
and refits both models. A coordinate change is rejected even when the coordinate
file's recorded hash is updated without rebuilding its dependent reports.

Freeze the analysis and structural-input selection without using held-out outcomes.
The already inspected published panel is retrospective development data. Geometry
must beat both the fixed score and the matched score-only model before an added
predictive-value claim is considered; an independent dataset and prospective
measurements are separate subsequent requirements.

## Bounded implementation

Per structure: at most 100,000 atoms, 400,000 bonds, 100,000 residues, 10,000 chains
and 32 conformers. Selected coordinate magnitudes are bounded by 1,000,000 nm.
A plan admits at most 100,000,000 atom-pair evaluations (default 25,000,000);
contact-residue pairs are limited to 100,000. This implementation uses direct
cross-partner pair enumeration, not an accelerated neighbor search.

Per structural input: at most 10,000 observations, 1,000,000 total atoms and
100,000,000 accumulated pair evaluations. Bundle files are limited to 128 MiB
and structures are resident. These are explicit admission limits, not measured
maximum biological scale. Streaming, repeated conformations, replicas and a full
270-design structural campaign are not qualified by the synthetic checks.

## Executed example and tests

A complete artificial example is runnable without model weights or a GPU:

```sh
python3 Tools/BinderBenchmark/native_package.py /tmp/numi-interface-native --test
python3 Tools/BinderBenchmark/check_structures.py \
  /tmp/numi-interface-native/.build/debug/BinderCLI
```

`check_structures.py` generates its inputs, invokes the real native import and
structural-evaluation commands, verifies outputs and independently checks geometry
and both fitted models in Python. It uses 17 deliberately artificial polyglycine
structures, 12 matched training rows and five matched test rows; one source row
has no structure. Coordinates are software fixtures, not physically prepared proteins.

At commit `cbcf613b5a283f2a45b28ef010006195f5b0e984`, focused CI passed all
50 Swift tests, including 15 interface tests, seven offline Python tests and this
complete synthetic structural campaign. Distances, cutoff boundaries, rigid motion,
sequence/element identity, target consistency, candidate-specific cutoff rejection,
missing structures, matched populations and tamper rejection are exercised.
The [published-data record](evidence/2026-09-18-published/README.md) concerns the
score-only baseline; actual published structures have not yet been evaluated here.
