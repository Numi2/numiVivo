# GSE181897: complete count handoff; prediction admission incomplete

**The released counts and native B-cell aggregation are verified. The proposed
external IFN-beta prediction experiment has not been fitted or scored.** Its
[protocol](PROTOCOL.md) was frozen before count inspection. Original treatment
codes remain untranslated pending primary code-to-intervention evidence.

## Source and verification

[GEO GSE181897](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181897)
reports 64 donors, five immune stimuli and twelve technical pools. Its primary
SOFT metadata describes IFN-beta at 500 IU/mL for nine hours and an unstimulated
control incubated under the same conditions. The original
`GSE181897_concat.4.raw.h5ad.gz` contains 136,142 filtered singlet cells, including
all 64 `exp_id` donors. Each donor ID maps one-to-one to a `free_id` genotype
cluster. This does not requalify demultiplexing or exclude overlap with earlier
participants. Kang and Ye are contributors to both this study and the earlier
Kang experiment; separate collection is not independent-laboratory validation.

All **292,741,570 stored values** are finite, nonnegative integers, with sorted,
unique feature indices in each row. Recomputed RNA and antibody totals equal
`mrna_n_counts` and `adts_n_counts` for **every cell**. The full sparse source
check took 15.64 seconds with about 555 MB peak RSS on the local host; this is a
source-validation observation, not a native-product performance benchmark.

The feature metadata contains an important contradiction: `feature_types` calls
all 20,399 features `Gene Expression`. The independent `genome` field identifies
**20,303 GRCh38 RNA features and 96 BD99AbSeq antibody features**. Those two
selections reproduce the author modality totals exactly. Antibody counts are
preserved in the original source and excluded from the RNA denominator. They
are not silently interpreted as transcript counts.

## Native handoff result

The external extraction retains every author `ct3=B` cell under every original
condition code, including `B_Naive`, `B_Mem` and `PB` sublabels in `ct2`. These
are author-defined populations, not newly learned or prospectively assigned
cell identities. It preserves the selected cells' original identifiers and the
explicit metadata columns listed by `read_metadata.py`; it does not claim that
all original AnnData slots were copied into the RNA transport file.

| Check | Result |
| --- | ---: |
| Selected B-lineage cells / donors | 15,272 / 64 |
| RNA features / source count records | 20,303 / 34,287,682 |
| Observed donor / condition groups | 379 |
| Native nonzero aggregate values | 3,818,645 |
| Every selected source count versus transport | Exact |
| Every native aggregate value and group membership | Exact |
| Native publication / reconstruction / repeated publication | Passed; repeated report bytes identical |
| Donors with both literal `B` and `C` codes | 62; IDs `5` and `23` lack a pair |
| Prior shared panel represented by exact RNA symbols | 11,800 / 11,884; all 84 absent symbols retained |

The native product command `singlecell-h5ad-pseudobulk` streams all 34.3 million
selected records; the source exceeds the common in-memory 32-million-entry
allowance. This exercises the existing streamed aggregation owner, not a new
out-of-core implementation. A separate source pass compares **every selected
row, feature and value** with the externally written H5AD transport. The native
aggregate is compared against independently summed source counts, including
all group membership, donor, condition, sample and pool identities.

The exact native executable is the existing interval-qualified M4 Pro build:
`0f08cc423a10a2910962b25039c90731f2ceff013ec32e9a6df09d547a87ae4b`.
The native HDF5 library SHA256 is
`a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e`.
Native publication, reconstruction and repeated publication took respectively
3.72, 3.75 and 3.69 seconds on the shared physical Mac mini. No controlled speed
comparison, GPU gain, new prediction fit or biological success is claimed.

## Why prediction scoring remains gated

The original `cond` values are `0`, `A`, `B`, `C`, `G`, `P`, `R`. The
[author repository](https://github.com/yelabucsf/clue/tree/6775f8f73096b594a836bfb0b9708e7032971a9e)
states that it is incomplete. Its `create_metadata.ipynb` lists the stimulation
codes and excludes `0` from filtered conditions, but does not define full code
names; the pinned tree lacks the production experiment notebooks. An
[external curator's code](https://github.com/hms-dbmi/dseqr.data/blob/419c822aa034e5c2ffc7b079f8dc2610712cf128/data-raw/azimuth/human_stimulated_pbmc/human_stimulated_pbmc.R)
maps `B` to IFN-beta and `C` to Control, but explicitly comments `C = control?`.
The retained source record pins that curator revision. This is supporting
interpretation, not the missing primary experiment mapping.

Native aggregation therefore uses `source-code:<letter>` without guessing
intervention names. The prospective test requires that mapping to be established
before fitting or scoring. Neither 62 apparent pairs nor canonical interferon
expression patterns resolve experimental identity. Once it is resolved, retain
both incomplete donors and the 84 absent panel genes; use all eligible pairs
under the frozen protocol and report failed predictions and interval coverage.
The broader [biological prediction assessment](../../../../Documentation/BiologicalPrediction.md)
remains unchanged by this input qualification.

## Storage, reproduction and retained failures

The original gzip is **1,011,162,509 bytes**, SHA256
`7fe58432f2f238319e81c9218eb35b5f7fbdae6f10f3d87ce9a6044ee851675b`.
Its decoded stream is 3,063,713,137 bytes, SHA256
`183d7756c750fb0ca57f381512fe784df6249ec5a5478a9caf6a62df55cba56c`.
CRC and both hashes were checked during download; the stored gzip hash was
checked independently. Indexed gzip access avoids retaining a redundant 3 GB
copy. This is Python source preparation; native direct gzip-H5AD input is not
claimed. No source dataset or earlier research evidence was deleted.

The exact workspace used is `/Users/home/numivivo-gse181897-20260911`, with
native execution at `/Users/n/numivivo-gse181897-20260911`. These cohort scripts
are an executable historical recipe: copy them into that study directory, then
run `download.py`, `inspect_metadata.py`, `read_metadata.py`, `validate_counts.py`,
and `prepare_handoff.py` in order. Preparation references the previous frozen
panel at `/Users/home/numivivo-cross-study-ifnb-20260911/inputs/panel.json`.
Do not overwrite completed directories. The Python environment used h5py 3.16.0,
AnnData 0.13.3.post0, NumPy 2.5.3 and SciPy 1.18.1; `indexed_gzip` 1.10.3 was
installed in the study-local `deps` directory. The index can be regenerated from
the original gzip. Transfer `handoff/` to the remote study directory and execute
`run_handoff.py` there; it names the qualified native executable and HDF5 library
explicitly. Copy the report/receipt/plan/log files back, then run
`check_handoff.py` locally. The retained H5AD is an RNA cell transport, not a
pseudobulk-row substitute.

[evidence/2026-09-11](evidence/2026-09-11) retains compact source admission,
plans, checks, logs, identities and the executed qualification scripts. Large
source, transport and complete native reports remain at the manifest-bound
external locations. Original author code is hash-bound at its pinned repository
revision and retained locally, not republished in this archive.

Retained attempts include the rejected decoded-download capacity check, an
inventory serialization failure on HDF5 object references, the initial native
missing-library failure, and a checker failure comparing unsorted SciPy aggregate
indices with canonical native indices. The checker now canonicalizes the
reference's sparse ordering before exact comparison; it changes no count,
identity, tolerance, native input or output. Challenged GEO interactive pages,
the inaccessible preprint route and a non-PDF thesis response were not bypassed.
