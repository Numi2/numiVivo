# SEQC technical bulk-RNA benchmark

This benchmark tests the native negative-binomial differential-expression owner
on a public reference-RNA experiment and compares its effects with edgeR,
limma-voom and DESeq2. It is a technical measurement check, not a biological
donor experiment: all six Illumina RefSeq tables use UHRR/HBRR reference-RNA
materials and their lane columns are summed into 54 technical library
preparations. A technical library is never presented as a cell or a biological
replicate. The full protocol is frozen in [PROTOCOL.md](PROTOCOL.md).

The source is `seqc` 1.46.0 from Bioconductor, archive SHA-256
`f54700f7ef029225455f8f6c83bf0211d106f68dc9656ed9ab9d7f82b5cce7ea`. The six
sites are AGR, BGI, CNL, COH, MAY and NVS, with 25,794 retained RefSeq rows per
site. Paired-end fragment counts remain explicitly `fragmentCount` throughout
the adapter and native request. The qPCR reference has 1,044 rows; 785 pass the
predeclared unique EntrezID, all-eight-detections-`P`, positive-finite-value
rule. Its reference effect is the mean log2 B intensity minus the mean log2 A
intensity. Every excluded row and reason is retained in
`taqman-eligibility.tsv`.

## Reproduce the evidence

The commands below use the prepared evidence root shown in the receipts. Native
inputs and frozen native results are produced by `prepare_counts.R`,
`prepare_native.py` and `run_native.py` after compiling the driver in
`Main.swift`. The reference fits use the exact same reordered counts, design
matrix and native median-ratio size factors:

```text
Rscript extract_taqman.R selected/taqman.rda taqman.tsv
python freeze_references.py --root EVIDENCE --reference-dir reference-closed
python score.py --root EVIDENCE --reference-dir reference-closed
python check_results.py --root EVIDENCE
```

`reference-freeze.json` records all 18 complete edgeR/limma/DESeq2 tables and
the earlier truncated attempts. It explicitly records that qPCR values were
not loaded at freeze time. `native-checks.json` verifies the frozen result
hashes, count/design membership, status inventories and conditional fitted
mean/standard-error arithmetic; its six sites each have 22,207–22,775 tested
features. The native output and reference tables are retained under the
Mac mini evidence root `/Users/n/numivivo-seqc-20260912`.

## Result

All six sites pass the narrow descriptive gate (coverage at least 90%, Spearman
correlation at least 0.90, and at least 90% directional agreement for reference
effects with absolute log2 magnitude at least one) for native MLE, native
shrinkage, edgeR, limma-voom and DESeq2. On the common family available to all
five methods, site-level Spearman correlations span 0.938–0.941 and directional
agreement spans 0.975–0.979. Native shrinkage reduces RMSE relative to native
MLE at every site (1.150–1.193 versus 1.155–1.209 log2 units). The full metrics,
method-specific coverage and common-family comparisons are in
`score-summary.json`; per-site rows are in `scores/*.tsv.gz`.

The gate is intentionally descriptive. AGR edgeR has very large raw effects for
zero-supported rows because this comparison uses `prior.count=0`; those values
are retained in the original-family metrics rather than clipped. The result
does not establish false-discovery control, interval coverage, biological
replication, intervention response, clinical benefit or general outcome
prediction. It adds a reproducible bulk-RNA technical benchmark to the shared
count-model evidence and leaves those broader claims open.
