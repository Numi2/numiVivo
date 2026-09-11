# Complete Parse IFN-beta count qualification

**PASS for complete selected-count ingestion, native source replay and independent
count reconstruction. No prediction was fitted or scored.** All 725,031 literal
IFN-beta/PBS cells, 12 donors, 24 groups, 40,352 source RNA features and
1,373,870,697 positive records were retained. The matrix contains 3,070,817,047
integer counts. Every cell total, detected-feature count, aggregate coordinate,
group membership and original-to-local row mapping was checked.

The native executable SHA-256 is `20b13e585e1526dba7be566484373a3e4775efebd66f22b3b20dff31c1ee3516`.
Each donor was streamed again from the pinned source, with an exit-zero native
reconstruction and exact agreement of all fetched range hashes and its complete
canonical stream digest. All 3,456 source runs cover the selected rows exactly
once. The twelve stream digests together identify 21,981,931,152 canonical bytes;
these are separate donor streams, not a claimed monolithic or whole-H5AD hash.

## Complete donor results

Peak RSS below measures the native process only, in MiB. It excludes the Python
extractor, metadata cache and other processes; this is not total pipeline memory,
a million-cell metadata bound or a Metal speedup measurement.

| Donor | Cells | Records | Ingest peak RSS, MiB | Replay peak RSS, MiB |
| --- | ---: | ---: | ---: | ---: |
| Donor1 | 110,923 | 241,117,037 | unavailable | 162.25 |
| Donor2 | 77,881 | 142,798,033 | 159.84 | 130.62 |
| Donor3 | 52,416 | 101,972,982 | 119.39 | 93.64 |
| Donor4 | 47,867 | 106,730,964 | 116.33 | 87.31 |
| Donor5 | 66,848 | 109,352,718 | 134.67 | 112.36 |
| Donor6 | 71,569 | 130,790,653 | 141.34 | 116.91 |
| Donor7 | 43,062 | 78,550,276 | 110.72 | 85.05 |
| Donor8 | 49,629 | 97,555,781 | 117.81 | 91.39 |
| Donor9 | 34,279 | 61,860,291 | 102.31 | 84.38 |
| Donor10 | 40,875 | 90,785,124 | 110.03 | 84.88 |
| Donor11 | 50,985 | 88,440,305 | 119.72 | 93.00 |
| Donor12 | 78,697 | 123,916,533 | 161.89 | 130.53 |

Donor1's original producer measurement was interrupted during checker recovery;
its ingestion RSS, exit status and producer digest remain unavailable. Its saved
native bundle and every independent count were recovered and checked, then the
complete fresh source replay passed separately. That replay does not reconstruct
the missing original process measurements.

## Preserved discrepancies and failures

Historical source `tscp_count` exceeds the retained matrix by **32,500** counts;
historical detected-feature totals exceed it by **31,902**. Each historical QC
field differs for **30,634** cells. No cell was removed to make the annotations
agree, and no missing mitochondrial annotation was guessed. The initial
monolithic upstream HTTP 500 failure, native truncation rejection, temporary
memory defect and donor-checker recovery remain in the retained evidence.

The [immutable preparation binding](README.md#immutable-preparation-binding)
rechecks all 7,437 original input files and 37 donor-plan files. Its seven actual
partition substitution controls remain separate from native count replay.

## Retention and restoration

The [count archive manifest](evidence/2026-09-11-counts/manifest.json) retains
10,055 logical files. Its compressed payload is
43,131,777 bytes with SHA-256 `ae94b3d586d2212fc09b1d451f893e7b3badc05ad57eb04550bece351af12b97`. The original preparation
and memory-control archives are required dependencies, preserving all axes,
source maps and the exact executable without another copy in this archive.

A separate restoration verified every archive object and all restored file
hashes, then ran the frozen complete-cohort checker with network connections and
child processes forbidden. It checked **17,556 restored files**,
using 17,555 verified APFS copies and
1 directly extracted file. APFS copies have separate
inodes and are rehashed; no hard links or changes to the original study are used.
A second complete restore without `--reuse-from` also passed, extracting
17,465 files directly and cloning 91
identical objects only within the new destination. The [restoration evidence](evidence/2026-09-11-restoration/manifest.json) records
both offline checks. Neither is another native execution or source download.

```sh
R=Tools/Omics/PerturbationPrediction/ParseIFNB
python "$R/restore_donors.py" --evidence "$R/evidence/2026-09-11-counts" \
  --out /absolute/path/to/new-restored-study
# Optional: --reuse-from /absolute/path/to/existing-study
```

Restoration rejects an existing destination and conflicting dependency contents.
Without reusable verified files it extracts with a fixed 1 MiB buffer and
checks a 200 MB free-space reserve before each extraction. It requires the repository dependency archives and
Python with NumPy. A fresh native replay still requires the pinned remote source.

## Biological scope

These results establish input fidelity and count arithmetic. The unchanged
prediction panel still lacks 409 exact source symbols; the HGNC audit supplies
271 naming candidates, not validated replacement measurements. Exact IFN-beta
dose/reagent identity remains unresolved. No predictive accuracy, independent
duration-model validation, tissue function or clinical outcome follows from this
count qualification. Source data and derivatives: Parse Biosciences, CC BY-NC 4.0.
