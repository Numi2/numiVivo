# AlphaGenome Atlas + RNA-supported splice evidence

## Status

NumiVivo now has a native, immutable genomic-evidence layer for public-reference
neoantigen research. The native Swift side prepares exact Atlas requests, imports
and reconstructs completed evidence bundles, validates pVACsplice/RegTools evidence,
and binds research reviews to the exact evidence revision.

The external adapters in `ReferenceAdapters/AlphaGenomeAtlas/` close the execution
gap. Native NumiVivo still performs no network request. Predictions remain separate
from measured tumor RNA and are never converted into a treatment recommendation.

Supported scope is synthetic demonstrations and public-reference research. There is
no patient intake, clinical authorization, vaccine manufacture, dosing, or treatment
release path in this feature.

## Atlas workflow

The Atlas bridge uses the official AlphaGenome Python SDK pinned to commit
`aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d`. NumiVivo archives the exact response
rather than inventing a remote model/dataset version; the inspected API does not
provide an immutable service-version selector.

Create an isolated environment:

```sh
python3 -m venv .venv-atlas
. .venv-atlas/bin/activate
python3 -m pip install -r ReferenceAdapters/AlphaGenomeAtlas/requirements.txt
```

Keep the API key only in the process environment:

```sh
export ALPHAGENOME_API_KEY='...'
```

Before retrieval, review the applicable AlphaGenome terms and save the exact bytes
locally as `terms.txt`. Create a separate `usage.json`:

```json
{
  "purpose": "noncommercialResearch",
  "termsURI": "https://deepmind.google.com/science/alphagenome/terms",
  "termsSHA256": "SHA256_OF_THE_REVIEWED_TERMS_FILE",
  "acknowledgedBy": "researcher-identifier",
  "externalSharingApproved": true,
  "snapshotStorageApproved": true,
  "modelTrainingAllowed": false
}
```

This declaration is a software precondition, not legal advice or an entitlement.
Do not set an approval field true unless the applicable terms and institution
permit that activity.

Retrieve the live scorer catalog first rather than guessing scorer identifiers:

```sh
python3 ReferenceAdapters/AlphaGenomeAtlas/atlas_bridge.py catalog \
  --usage usage.json --terms terms.txt --output atlas-catalog
```

Create selectors using exact scorer names from that catalog:

```json
{
  "requestedScorers": ["EXACT_SCORER_NAME"],
  "ontologyTerms": ["CL:0000000"]
}
```

Generate the native request from an already verified parent neoantigen receipt:

```sh
swift run numivivo neoantigen-atlas-request parent-receipt.json \
  --store neoantigen-store --selectors selectors.json --output atlas-request
```

The request supports explicit GRCh38 single-nucleotide substitutions on primary
chr1–22/X/Y only. Unsupported candidates remain explicit exclusions. No liftover,
reverse-complement guessing, indel normalization, or missing-to-zero conversion is
performed.

Retrieve Atlas evidence against the exact case reference FASTA:

```sh
python3 ReferenceAdapters/AlphaGenomeAtlas/atlas_bridge.py fetch \
  --request atlas-request/atlas-request.json \
  --fasta reference.fa --usage usage.json --terms terms.txt \
  --output atlas-bundle
```

The bridge streams and hashes the full uncompressed FASTA and independently checks
every queried reference base. A wrong genome file rejects the batch. Per-variant
reference mismatches are explicit failures. Provider errors are redacted into bounded
categories so arbitrary server text or credentials cannot enter evidence records.

A completed bundle contains `atlas-request.json`, `atlas-capture.json`, and a
`complete.json` written last. Import it natively:

```sh
swift run numivivo neoantigen-evidence-import parent-receipt.json \
  --store neoantigen-store --atlas-bundle atlas-bundle \
  --output evidence-review

swift run numivivo neoantigen-evidence-verify evidence-review/receipt.json \
  --store neoantigen-store
```

`available`, `noData`, and `failed` are distinct. Raw scores and calibrated quantiles
remain distinct; signed quantiles retain direction. Quantiles are not probabilities
of vaccine benefit. Atlas predictions do not replace measured tumor RNA or establish
presentation, immune recognition, safety, or clinical benefit.

## Managed splice workflow

Generate the prepared job:

```sh
swift run numivivo neoantigen-splice-job parent-receipt.json \
  --store neoantigen-store --output splice-job
```

Fill every placeholder in `splice-job/splice-job.json`. Supply the local IEDB
prediction-tool installation separately when executing the runner. The runner computes
and records a complete path-to-file-hash inventory for that directory, rather than
trusting a mutable directory path as provenance.

The runner verifies the six input files, RegTools and pVACsplice executables, the
complete predictor directory inventory, and any additional `resourceFiles` entries.
Only explicitly named local class-I algorithms are accepted by this profile. It
runs RegTools first and then pVACsplice using the same VCF, FASTA, GTF and RNA BAM
context. It does not align RNA, call variants, type HLA, or install predictors.

```sh
python3 ReferenceAdapters/AlphaGenomeAtlas/run_pvacsplice.py \
  --job splice-job/splice-job.json \
  --iedb-install-directory /ABSOLUTE/PATH/TO/IEDB_MHC_TOOLS \
  --output splice-bundle
```

The completion bundle contains:

- `splice-manifest.json`
- `splice-report.tsv` — unaggregated class-I all-epitopes output
- `regtools.tsv` — the exact generated junction evidence
- `transcripts.fa` — the pVACsplice paired protein FASTA
- `splice-execution.json` — bounded executable/input/output provenance
- `complete.json` — written only after the evidence files are complete

Native import checks report rows against RegTools coordinates, read counts, anchors,
transcript identities, supplied ALT/WT proteins, candidate peptides, case HLA and
input hashes. It still does not establish allele-specific splice causality, matched-
normal RNA absence, antigen presentation, immune recognition, or therapeutic effect.

Import splice evidence alone:

```sh
swift run numivivo neoantigen-evidence-import parent-receipt.json \
  --store neoantigen-store --splice-bundle splice-bundle \
  --output splice-review
```

Or generate a new Atlas request that also includes the newly supported splice
candidates, retrieve it, and import both bundles in one evidence revision:

```sh
swift run numivivo neoantigen-atlas-request parent-receipt.json \
  --store neoantigen-store --selectors selectors.json \
  --splice-bundle splice-bundle --output splice-atlas-request

python3 ReferenceAdapters/AlphaGenomeAtlas/atlas_bridge.py fetch \
  --request splice-atlas-request/atlas-request.json --fasta reference.fa \
  --usage usage.json --terms terms.txt --output splice-atlas-bundle

swift run numivivo neoantigen-evidence-import parent-receipt.json \
  --store neoantigen-store --splice-bundle splice-bundle \
  --atlas-bundle splice-atlas-bundle --output combined-review
```

## Verification

Adapter unit tests do not require a live AlphaGenome key or installed pVACtools:

```sh
cd ReferenceAdapters/AlphaGenomeAtlas
python3 -m unittest -v
```

The tests exercise exact request identity, FASTA/reference mismatches, usage gates,
explicit Atlas outcome states, duplicate candidate rejection, prepared-job input
hash rejection, public-reference gating, and a fake RegTools+pVACsplice end-to-end
bundle with completion hashes.

A passing adapter test is not a live-provider, biological, clinical, or regulatory
validation. Live Atlas access, an actual public-reference pVACsplice reproduction,
and the repository's native Apple test suite remain separate qualification gates.

## External references

Contract checked against pVACtools 7.1.x documentation and RegTools' public command
interface. pVACsplice requires an annotated VCF plus RegTools junction TSV and writes
an unaggregated class-I all-epitopes report and `<sample>.transcripts.fa`. RegTools
`cis-splice-effects identify` combines the annotated variants, RNA BAM, reference
FASTA and GTF. Keep exact tool versions and local resource identities in the job.
