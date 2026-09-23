# PXD004894 phase-2 shared HLA-I peptide protocol

Frozen before any phase-2 PXD004894 RAW search or peptide-result inspection.
The endpoint is **one exact, unmodified, non-cysteine, 8–15 aa peptide** from
one GENCODE TransCODE Phase-2 Primary ORF, observed by direct MS2 in multiple
native melanoma donors and then independently in a melanoma cell-line donor.
A valid negative outcome is `null_no_eligible_shared_peptide`,
`lead_failed_audit_no_replacement`, or
`heldout_confirmation_failed_no_replacement`. No threshold relaxation,
replacement lead, or backup candidate is allowed.

Run `verify_phase2_protocol.py` against the nine real fixed inputs before the
first production RAW and again before adjudication. It checks source hashes,
the 88-run/25-donor split, the 16-run/2-donor held-out split, catalog mappings,
and the TRC BED12/FNA/FAA relationship without reading peptide results. The
first-RAW `run_first_raw_pilot.py` is **exploratory**; its receipts and TXT are
never among the 88 production discovery inputs.

## Frozen inputs and separation

| Input | SHA-256 |
| --- | --- |
| PXD004894 HLA-I RAW inventory, exact TSV order | `686dbc1bcd43bf5046c57afb3627670c3b79bbc3f7b2e311ac346af2d5af4661` |
| PXD004894 donor/HLA summary | `3b258f9740fdc11fe93f7f8ef36aa8f183cdfb66a8b328f068e4bb564a1f3e96` |
| MSV000084787 validation inventory | `172175ffd1411a58d42bb4ddf496e5501745a61790948ccfdbd5fb91b18f4b92` |
| CON/UPR/TRC combined targets FASTA | `a73ea6ab3e8c28262d5effec3a0baf5d98bab25153cf6da5d3b68a435db5193d` |
| Combined catalog protein provenance TSV | `381e001d4f992a7dd42b2ed8f4e8f6db0e897f70ee47bab5921434df579d4e0c` |
| Combined catalog manifest | `444a4b589f1e317a5d97a07c25b64bc639200ed891c5c3cbb6d5f6f48ad8b9b9` |
| TransCODE Phase-2 Primary BED12 | `e05acb0ab51094f0ca96f35ff37864f6c2addfac22e004f09712e93af8634f6f` |
| TransCODE Phase-2 Primary spliced FNA | `558ee88208b60353312f23f6f3459f1ac3352d13de54353b8bcc86b5c34d0d7c` |
| TransCODE Phase-2 Primary FAA source | `521d95a14c49fd625b33c704452430d6fdd89dca2a97a8c66a57897ee3dadbe9` |
| TXT-only Comet parameter template | `810d1563226721f427d55261f46920ceced308b06d2b5f9b9efc5f50af6ac2fd` |
| Independent Sage production JSON template | `73ceb13f07926c2ced63eda02f14464b87c5741f46bc3f25df552aab5a854c29` |

PXD004894 contributes exactly 88 HLA-I RAW files from 25 donors; its 52 HLA-II
files and submitted `Search.zip` are excluded. MSV000084787 contributes exactly
16 prelisted mzML runs from Mel02 and Mel11, held for confirmation **after**
one discovery sequence freezes. Donors are counted by the frozen donor field,
never by injection, fraction, treatment, scan, or file. The PXD donor labels do
not overlap Mel02/Mel11. Labels alone do not establish unrelated patients;
publication provenance must close that claim before final acceptance.

The discovery catalog contains only CON, UPR, and TRC target sequences.
CON and UPR compete but cannot nominate. TRC Primary membership and translation
scores are fixed by GENCODE/TransCODE public Ribo-seq resources; no Mel02/Mel11
TPM threshold, MSV peptide identity, or PXD search result enters selection.
NUO is excluded because its source incorporated MEL2/MEL11 Ribo-seq; SMO is
excluded because its independent provenance was not established. The Primary
FNA/FAA/BED all contain 10,127 matching ORFs; 6,771 TRC proteins of length at
least 8 aa entered the search catalog. The external catalog sources are
[GENCODE](https://www.gencodegenes.org/pages/riboseq_orfs/) and the
[TransCODE Phase-2 publication](https://academic.oup.com/nar/article/54/6/gkag234/8539537).

## Production acquisition and search

1. Process the exact 88-row inventory order, one RAW per transaction. Before
   every transfer, require the live PRIDE file accession, project, RAW category,
   name, byte size, SHA-1, and HTTPS location to equal the row. Bound the `.part`
   download to the pinned byte count; atomically rename only after SHA-1 and
   byte verification, and retain SHA-256. Validate converted mzML XML, nonzero
   MS2, and SHA-256. Keep all failed/unverified files.
2. Before the first production RAW, verify the pinned Sage binary SHA-256
   `7e457a288124a25da216e90ae0b9c301a4ed8774536d031ca72ce964520903ea`
   (asset v0.14.7, runtime v0.14.6), and run a bounded resource preflight
   against the **published catalog**. An earlier Sage probe against a different
   140,319-protein catalog terminated abnormally at 18.18 GB maximum RSS; it
   does not prove this 27,596-protein catalog is executable. No production RAW
   begins if the exact-catalog preflight fails or lacks room for both engines.
   Search every mzML against **the same** combined FASTA with Comet 2026.02
   rev. 2 (`1b93ed1cf690026a75d80e1e0ce3ed57394bcd47ba5c8587441668c006e32f0e`)
   and the pinned TXT-only template, rendering only `database_name` to the
   verified catalog path. The target-only FASTA uses Comet `decoy_search=1`,
   giving concatenated target/decoy competition with `DECOY_`. No PIN, pepXML,
   SQT, or mzIdentML output. No-enzyme 8–15 aa, 600–4000 Da, ±10 ppm precursor,
   HCD b/y fragment settings and the frozen modification set are in the template.
   In **the same RAW transaction**, search the same mzML with Sage using
   `sage_production.json.template`. Render only the catalog, per-RAW output
   directory, and mzML placeholders; hash the rendered config. Its no-enzyme
   chemistry, generated decoys, and all search bounds are fixed in the JSON.
3. Validate exact Comet TXT and Sage result schemas and terminal completion.
   Gzip each peptide-bearing output at level 6; verify compressed and
   decompressed SHA-256. Durably seal **both engines'** complete result sets,
   the transaction receipt, version/command/input/config/output hashes, logs,
   and resource metrics. A production receipt must say
   `analysis_kind=production_88_run`, carry the inventory row identity and
   both engines' output hashes, and be distinct from the exploratory pilot
   schema. Only then remove that transaction's verified public RAW, mzML,
   and losslessly compressed plain outputs. The Comet-only first-RAW pilot
   remains exploratory and cannot substitute for a production receipt.
4. Require at least **7 GiB actual free** before each file; abort if free falls
   below 2 GiB, the combined one-file working set of both engines exceeds
   5,210,336,441 bytes, or
   retained phase-2 output exceeds 4 GiB. Before continuing at a stop gate,
   copy sealed evidence to a separately verified archive and reopen/hash it
   there. The compression estimate is planning only. No partial 88-run set can
   nominate or yield a null result.

The [Comet decoy parameter](https://comet-ms.sourceforge.net/parameters/parameters_201801/decoy_search.php)
defines the 1:1 target/decoy competition used here. No peptide-bearing TXT is
read for nomination until all 88 production receipts and both engines'
compressed outputs pass hash/schema/completion checks. Search logs and
metadata may be checked before then without selecting peptides.

## Locked adjudication

For each of the 88 complete Comet TXT files, group all `num=1` rows by
`(RAW basename, scan, charge)` before source classification. A target PSM
requires one exact peptidoform across the group, target-only protein mappings,
and fewer than five tied rows (five is the output cap). Every other group is a
conservative decoy; preserve every tied row and reject its candidate support.
For targets use the least favorable tied XCorr/E-value; for decoys use the most
favorable. Non-rank-1 rows cannot rescue a group.

Calculate target-decoy q-values once across **all 88 runs**, at tied score
thresholds, using `min(1,(decoys+1)/max(targets,1))` and monotone tail minima.
Do this independently for higher-better XCorr and lower-better E-value at both
rank-1 PSM and exact peptidoform levels. Repeat all four calculations in one
pooled TRC subgroup: pure TRC target groups count as targets; any group with a
TRC mapping and a decoy or CON/UPR ambiguity counts as a subgroup decoy.
Peptidoforms observed as both target and decoy count as decoys. A supporting
PSM requires **all eight** global/TRC PSM/peptidoform q-values ≤0.01, and no
target/decoy identity collision. If the TRC decoy pool is insufficient, the
`+1` correction can yield no passing target; that is a valid null.

An eligible peptide is unmodified, cysteine-free, 8–15 standard amino acids,
and must have direct passing spectra in **at least two distinct PXD donors**.
The exact and I/L-normalized peptide must occur in no CON/UPR sequence and in
exactly one original TRC ORF accession and offset across the complete catalog
provenance. Its BED12 locus/strand and spliced FNA reading frame must resolve
without ambiguity; a peptide requiring the artificial initiating M of a
near-cognate start is excluded from this first endpoint. A tied canonical,
contaminant, variant, PTM, isotope, chimeric, or unresolved I/L explanation
rejects the peptide. Genome coordinates and one ORF are source attribution,
not evidence of HLA restriction.

Rank eligible sequences by (1) PXD donor count descending, (2) worst donor's
best qualifying `delta_cn` descending, (3) TRC ORF accession ascending, then
(4) exact peptide sequence ascending. No expression score from Mel02/Mel11
enters ranking. Independently adjudicate the sealed Sage results from **all
88 PXD runs**, generated from all spectra before any peptide nomination, not
a candidate-selected subset. The frozen Sage template specifies no enzyme,
8–15 aa, 600–4000 Da, precursor/fragment ±10 ppm, charge 1–4, isotope error
zero, b/y ions, fixed C+119.004099, variable M+15.994915, N-terminal
Q−17.026549 and C−61.982635 with at most two variable modifications, and
at least five matched peaks.
The first Comet-ranked lead must also have matching Sage direct spectra in
each counted PXD donor. Recompute Sage spectrum and exact-peptide q-values
from the sealed target/decoy rank-1 hyperscores pooled across all 88 runs,
using the same `+1`, score-tie, and monotone-tail rules; retain target/decoy
identity collisions as decoys. Require both pooled Sage q-values ≤0.01 and
the reported posterior error ≤0.01 for each counted match. Audit open
modifications and alternative canonical, variant,
isotope, and chimeric fragments for every counted scan before freezing one
sequence, source ORF, supporting scans, and all hashes. If the lead fails,
report `lead_failed_audit_no_replacement`; do not try the next ranked peptide.
If there is no eligible Comet sequence, record `null_no_eligible_shared_peptide`.

## Held-out confirmation and claim boundary

Only after the one-candidate receipt is read-only and hashed, search **all 16**
pinned MSV mzML runs against the same full catalog with the same Comet chemistry
and independently with Sage; no candidate-selected run subset, transfer, or
reuse of phase-1 accepted PSMs. Recompute global and pooled TRC target-decoy
Comet q-values and pooled Sage spectrum/exact-peptide q-values within these
16 runs, separate from PXD, using the same methods as discovery. Count only
direct MS2 matches
to the frozen exact unmodified sequence that pass both engines, all eight
Comet q-values ≤0.01, Sage spectrum/peptide q ≤0.01 and posterior error ≤0.01,
and candidate-specific alternative fragment review. Require at least one
Mel02/Mel11 donor, at least **three distinct donors total**, at least **five**
direct supporting spectra total, and observations in separate RAW/mzML files
in at least two donors. Failure is a failed confirmation, with no replacement.
Record `heldout_confirmation_failed_no_replacement`, retaining the negative
spectra and all search evidence.

After confirmation, run the already frozen exact/I-L novelty audit and normal
HLA tissue challenge. A prior exact/I-L report blocks “previously unreported.”
Any qualifying normal-tissue presentation blocks tumor restriction; inadequate
normal coverage remains unknown. HLA allele predictions or motif compatibility
are supporting only. Computational MS2 presentation is not T-cell recognition,
safety, or vaccine efficacy. Existing phase-1 MSV work makes this a
candidate-blind role reversal, not an assertion that these files were never
previously analyzed.

## Mini staging once its owner releases it

The Mini canonical `/Users/n/numiVivo` is stale. Fetch the published commit
into that repository **without checking out canonical main**, create a new
detached worktree with `git worktree add --no-checkout`, set sparse checkout to
`Tools/Neoantigen/MelanomaPhase2`, and only then materialize the exact commit.
The runner verifies a clean Git checkout before network access and binds HEAD
in its receipt. Stage the catalog and small pinned inventory/BED/FNA inputs
outside the checkout; verify every hash above and recheck free space after
staging. Do not stage or execute while another Mini workload owns the machine.
