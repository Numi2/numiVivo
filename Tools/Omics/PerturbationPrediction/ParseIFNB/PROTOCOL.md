# Parse IFN-beta source admission

This protocol admits counts for a future independent expression-prediction test.
It does not fit or score a prediction model. The source roster and conditions
were selected before count-outcome inspection; no cell or donor is omitted on
response magnitude, annotation agreement, QC discrepancy or model performance.

- Use the original Parse 10M PBMC H5AD, not a 100K/1M tutorial subset or publisher
  pseudobulk. Source URL:
  `https://parse-wget.s3.us-west-2.amazonaws.com/10m/Parse_10M_PBMC_cytokines.h5ad`.
- Bind every HTTP range to ETag `"f0644aed6d5db8e2fdf394523f3e08f2-6780"`,
  exact size 227,497,986,816 bytes, status 206, Content-Range and returned ETag.
  Retain each requested range's SHA-256. A multipart ETag is not a content SHA-256.
- Check all 9,697,974 source donor, cytokine, sample and treatment code assignments.
  Preserve every one of the 725,031 cells whose literal cytokine is `PBS` or
  `IFN-beta`, all 12 source donor IDs, and all 40,352 source RNA features.
  Preserve exact original barcodes and source row coordinates. Do not rename
  donors from demographic interpretations or manufacture replicate wells.
- Sum by original donor and condition. All PBS cells contribute. Cell numbers
  do not become independent biological replicate counts. Source cell-type labels
  are not used to select cells or qualify new labels in this admission.
- Read all 1,373,870,697 selected X records in bounded ranges. Require finite,
  positive integral values and valid columns. Sort columns within rows if needed;
  reject duplicates. Emit canonical records to the actual native stream owner.
- Validate each native row cardinality against source CSR pointers. Independently
  sum every retained count and compare every native cell total and every one of
  the 24 donor/condition by 40,352-feature aggregate entries. Retain differences
  from historical source `gene_count` and `tscp_count`; do not change or filter
  cells to make those annotations agree with the retained X matrix.
- Replay the complete stream, bind its SHA-256, recheck every range hash and
  reconstruct native output. Software agreement is separate from biological
  prediction accuracy. Preserve failures and incomplete stages.

[Parse's original methods](https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/)
describe 24-hour stimulation. The [primary author preprint](https://doi.org/10.21203/rs.3.rs-8337240/v1)
places specific cytokine doses in Supplementary Table 2; the IFN-beta dose and
reagent are unresolved in the currently accessible evidence. The pinned
[author repository](https://github.com/theislab/HumanCytokineDict/tree/6f9bc00381227fe8b1aa8dcb4e3f6f168cfc3229)
and [huCIRA cytokine table](https://github.com/theislab/huCIRA/tree/eb0b5b2633312f881de568066672b4b16dbd13ab)
identify the source and IFNB1 name but do not establish that dose.

The unchanged duration panel has 409 symbols absent by exact name after its
`symbol|` namespace wrapper is applied. Resolve identity or explicitly freeze a
supported panel before fitting; do not zero-fill missing genes or silently change
the evaluated endpoint.

Before future scoring, declare the verified intervention, dose, exposure,
culture, source roster, expression units, training inputs, unchanged model,
feature panel, primary baseline comparison and failure criteria. A 24-hour
cross-study result would test combined context transfer at one time; it would
not validate an entire temporal trajectory or tissue/clinical outcomes. The
prior GSE226572 duration model is development evidence and must not be tuned on
these outcomes and then described as independently validated.

Original data are attributed to Parse Biosciences and licensed **CC BY-NC 4.0**.
Retained source extracts and count derivatives keep that attribution and license;
the repository's software license does not relicense them.

## Execution amendment after a source failure

The first monolithic stream received upstream HTTP 500 after 1,166,913,217
records. Native truncation rejection prevented a complete result from being
published. The retry preserves that failed attempt, the original complete cohort,
all features and all biological sample identities. It uses the corrected temporary
memory owner, bounded retries for transient responses, and twelve disjoint complete
donor partitions. Every source row/run must occur exactly once across those
partitions. Native donor outputs are independently checked before a restart receipt
is saved, and all twelve must pass full native replay. Count-stream SHA-256 values
are now per donor; no monolithic digest is fabricated from those hashes. The
execution amendment does not fit a model, filter outcomes or alter scoring gates.
