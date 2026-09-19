# Original PDB/mmCIF sources in binder evaluation

The existing native binder workflow now accepts original structure text rather
than requiring callers to construct `VivoMolecularStructure` JSON. This is a
source-to-geometry/evaluation route, not a structure predictor, molecular
preparation workflow, affinity calculation or biological qualification.

## Run

Prepare a source manifest with one entry per available candidate. Every entry
must name its exact candidate ID, target, caller-declared source label, `pdb`
or `mmcif` format, relative `sourcePath`, expected lowercase `sha256`, exact
`targetChainSequences` mapping, and the existing `interfacePlan`. Paths are
relative to the explicit source root, with no absolute paths, traversal or
symbolic links. No target, model or chain assignment is inferred from a filename.
The existing structural plan still chooses the model features and held-out split.

```sh
python3 Tools/BinderBenchmark/prepare_structure_sources.py \
  source-manifest.json original-structures raw-sources.json

numivivo binder-evaluate-structure-sources \
  IMPORT_BUNDLE PLAN.json raw-sources.json NEW_RESULT
numivivo binder-verify NEW_RESULT
```

The transport helper preserves the exact UTF-8 source bytes, including CRLF.
Native `VivoPDB` and `VivoMMCIF` own parsing and the angstrom-to-nanometer
conversion. The source's chosen conformer must use the native parser's actual
identifier (for example `model-1`). Keep alternate locations, periodic cells,
occupancies, missing atoms and chemistry explicit; this route does not repair,
unwrap, truncate or resolve them. Existing interface/sequence admission rules
continue to apply. The expected sequence is required for every selected target
chain; the binder sequence must match its published experimental candidate.

A source hash binds bytes, not upstream authenticity or model correctness.
Expected target sequences are declared references, not proof of biological
identity. Freeze their origin, model selection and the experimental plan without
using held-out outcomes; already-inspected targets remain development data.

## Complete reconstruction

`VivoBinderBundleIO.evaluateStructureSources` publishes the new
`sourceStructuralEvaluation` kind. It retains original source text and declared
references in `structure-sources.json`, native reconstruction in `structures.json`,
the original experimental CSV, import configuration, evaluation plan, computed
geometry, and both matched reports. Existing bundle kinds remain readable.

Verification reparses the original structure text, repeats sequence admission,
compares the reconstructed native structure, recomputes geometry, and refits both
models. Merely changing coordinates and updating their recorded hashes does not
make stale derived reports valid. Publication uses the existing same-parent
private-directory rename; failures do not publish a completed result.

Missing candidates remain unavailable, not zero observations. The combined and
score-only fits use identical complete-case populations and the existing
training-target exclusion and training-only scaling rules.

## Bounds and execution evidence

This first route is bounded and resident, not a full-panel streaming claim:
16 MiB of UTF-8 text per structure, 64 MiB cumulative raw text, 128 MiB per encoded
input/bundle file, at most 10,000 source entries and 1,000,000 parsed atoms.
The existing per-structure and aggregate interface-work limits are unchanged.
For larger panels, develop the existing owners' streaming interface rather than
silently splitting the evaluation population or enlarging these limits.

Local Linux Swift 6.2.1 qualification: the exact-source binder/parser subset
builds; 58 Swift tests pass, including eight new raw-source tests. Ten Python
transport tests pass. The full synthetic CLI check runs 17 structures through
both embedded and mixed original-PDB/mmCIF inputs, preserves the matched
12-training/5-test population and one missing structure, independently checks
geometry and both fits, and rejects altered source/derived artifacts.

These are software results on deliberately artificial glycine coordinates.
No published structural panel, folding weights, molecular dynamics, Apple/Metal
execution or improvement in experimental candidate selection is established.
