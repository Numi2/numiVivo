# Neoantigen research workbench — first implementation

## Purpose and current boundary

Make completed neoantigen-analysis results easier to inspect without asking a
researcher to assemble a review UI. This increment adds a native Swift case
checker and pVACseq import adapter, reproducible artifact storage, and an offline
browser review report. It is **not an end-to-end personalized vaccine service**.

Only `synthetic` and `publicReference` data classes are supported. There is no
patient intake, consent-management system, authenticated clinical review,
treatment authorization, vaccine design/manufacturing, or administration path.
Do not use this release to store identifiable patient information.

## First use on the existing Apple build

```sh
swift run numivivo neoantigen-example --output ./neoantigen-demo
open ./neoantigen-demo/report.html
```

The generated example has three **invented** candidate rows. No HCC1395 data,
biological measurements, real prediction-model execution, or validation are
represented. The example includes missing RNA/expression values, unchanged
wildtype sequence, and matched-normal alternate reads so the interface's
limitations are visible. Synthetic resource digests describe invented fixture
inputs, not real genome files or an actual container image.

The report supports search, evidence-gap filtering, original source-field
inspection, draft dispositions and rationale capture. Browser review is local:
there are no network requests, hosted assets, external model calls, analytics,
or automatic uploads. Drafts are held in memory, not persisted in localStorage;
export before closing or reloading the page.

## Import, archive, reconstruct, review

```sh
swift run numivivo neoantigen-import ./neoantigen-demo/case.json \
  --tsv ./neoantigen-demo/all_epitopes.tsv \
  --store ./neoantigen-store --output ./neoantigen-import
open ./neoantigen-import/report.html

swift run numivivo neoantigen-verify ./neoantigen-import/receipt.json \
  --store ./neoantigen-store

# After exporting a review from the report opened above:
swift run numivivo neoantigen-review ./neoantigen-research-review.json \
  --receipt ./neoantigen-import/receipt.json --store ./neoantigen-store
```

An output directory must be new and outside the artifact store. Fixed export
filenames are `case.json`, `all_epitopes.tsv`, `report.json`, `report.html`, and,
for stored imports, `receipt.json`. Existing results are not overwritten.

`neoantigen-verify` reloads the archived original manifest and TSV, regenerates
the report, and compares canonical bytes. A modified result fails even when it
has a valid new SHA-256 digest. Receipts bind the existing CLI implementation
fingerprint (executable and platform). Use the originating binary for a receipt;
a new binary requires a new import rather than relabeling old evidence.

Reviews reference the exact report digest and candidate identifiers. Each
selected disposition requires a reason. Retaining means `retainForResearch`,
not treatment approval. Reviewer identity is self-declared, not authenticated.
Each changed review is a new immutable artifact. This does not yet provide a
signed, ordered clinical audit log or a multiuser review-conflict protocol.

Exit codes: `0` completed; `1` malformed/unsupported input or I/O error; `2` a
reconstructed/imported report has blocking case findings. Successful integrity
verification does not change clinical qualification, which remains absent.

## Supported input contract

The case is `numivivo.org/neoantigen-case/v1`. Generate the example as a complete
machine-readable template. It declares a pseudonymous case and subject, one
sample record each for tumor DNA, matched-normal DNA and tumor RNA, explicit
GRCh37/GRCh38 identity, reference/annotation digests and release identifiers,
matched-normal HLA source, HLA nomenclature release, external tool version,
execution image digest, predictor versions, annotated-VCF digest, recorded
argument vector, declared input mapping, and exact source-report digest.

These checks establish **metadata consistency**, not sequencing QC. Matching
identifiers cannot establish that samples are genetically matched. Only the
supplied TSV bytes are hashed against a declared source digest by this adapter;
raw read files, references, model resources and external execution must still be
independently verified. The source citation is descriptive, not authorization.

Adapter identifier: `pvacseq.classI.all_epitopes.v1`. It accepts the unaggregated
pVACseq TSV column contract, not pVACbind or the aggregated/pVACview format.
Supported alleles are unique two-field `HLA-A`, `HLA-B` and `HLA-C` names, with
8–15-residue standard-amino-acid peptides. Higher-resolution/suffixed alleles,
class II, custom coordinate schemas and aggregated outputs need separate
qualified adapters; do not silently convert or relabel them. The adapter is
column-contract based; a declared tool version is not a tested-version claim.
The required fields are listed in the production parser, and all source fields
are retained. Start/Stop remain zero-based, half-open, as documented upstream.

Files must be UTF-8, at most 16 MiB, at most 10,000 data rows and 256 unique
headers, with bounded fields. LF and CRLF are supported. Malformed rows, duplicate
complete rows, foreign HLA alleles, invalid coordinates, nonfinite/negative
recognized numerical values, fractional depths and VAF outside [0,1] fail the
entire import rather than silently removing candidates. Numeric IC50 and
percentile fields remain distinct. `NA`/empty remain absent, not zero.

The display preserves source order and exposes missing measurements. Gene
expression does not substitute for mutant RNA evidence. Wildtype-equivalent
peptides and alternate reads in the normal are explicitly flagged. No binding
cutoff or clinical ranking is supplied. An upstream `Evaluation=Accept` field is
retained only as an original source field, never promoted to a reviewed decision
or to observed immune recognition. Empty candidate reports remain explicitly
empty; they do not establish that a tumor lacks neoantigens.

## Implementation

- `Sources/NumiVivoKit/Neoantigen/VivoNeoantigenWorkbench.swift`: typed contracts,
  metadata preflight, bounded TSV parser, evidence gaps, report/review binding.
- `VivoNeoantigenExample.swift`: explicitly synthetic deterministic fixture.
- `VivoNeoantigenHTML.swift`: self-contained browser review renderer.
- `VivoNeoantigenArtifacts.swift`: existing `VivoArtifactStore` integration and
  reconstruction; no parallel store or mutable reference scheme.
- `Sources/NumiVivoCLI/VivoNeoantigenCLICommands.swift`: CLI integration.

HTML carries report bytes as base64 and uses `textContent`, not source-derived
HTML insertion. Browser SHA-256 checks those exact bytes. Export is disabled
when the browser does not provide WebCrypto or the digest fails. A document's
embedded hash does not authenticate the document author or code; native receipt
verification against a trusted originating executable remains distinct.

## Verification performed for this increment

```sh
bash Tools/Neoantigen/check.sh
```

The portable harness compiles the actual workbench, example and HTML renderer
with the repository's existing harness-only OpenSSL/CryptoKit support. **37
checks passed on Linux Swift 6.2.1**, including deterministic reconstruction,
missing-versus-zero values, metadata mismatch, report transplantation, invalid
numeric/coordinate/HLA/peptide data, duplicate headers/rows, size bounds, empty
reports, stale review bindings and script-text containment.

The offline DOM rendering was exercised in Chromium with search, filtering,
draft persistence, mobile-width overflow checks and no page errors. Browser
policy blocked file and loopback navigation in this environment; the available
in-memory page was not a secure context. **Browser WebCrypto/export round-trip
was not verified**. Export correctly stayed disabled without WebCrypto.

The artifact-store and CLI integration sources were syntax-checked. Native
`VivoNeoantigenIntegrationTests` cover reconstruction, a correctly hashed forged
report, review binding and immutable history, but **the Apple package and these
native integration tests were not run here**. No pVACseq execution, clinical
validation or real reference-case benchmark is claimed.

## Next completion gate

Qualify a pinned, public reference-data run against its actual upstream output,
then add a managed upstream execution adapter with raw-input identity/QC evidence,
error/cancellation handling and audited resource provenance. Validate the native
Apple integration and secure-context browser export. Add the guided input
wizard only over those verified contracts. Patient-facing and clinical-team
features require separate permissions, privacy, authentication and institutional
validation work; a research label alone does not supply those capabilities.

## Upstream technical references

Contract checked 2026-09-08 against the pVACtools documentation (served as 7.1.3):

- https://pvactools.readthedocs.io/en/latest/pvacseq/output_files.html
- https://pvactools.readthedocs.io/en/latest/pvacseq/filter_commands.html
- https://pvactools.readthedocs.io/en/latest/pvacview/pvacseq_module/pvacseq_upload.html

These references define upstream columns and distinctions. They do not validate
this implementation or the invented example, and `latest` URLs are not immutable
execution provenance.
