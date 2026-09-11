# Parse fixed-panel feature identity audit

The original exact-name comparison remains **12,584 present / 409 absent** out
of the frozen duration model's 12,993 response symbols. No count, feature axis,
model panel or prediction has been changed by this audit.

The [HGNC approved set](https://www.genenames.org/download/) supplies current
symbols and explicitly curated previous symbols. Its 2026-09-09 source snapshot
contains 45,054 approved rows. The complete decoded TSV SHA-256 is
`6a1423507780773fbcb373eaef84dba5f9357b36e7219db2aeaae13d8a42d62b`.
HGNC data are [CC0](https://genenames.org/download/custom/). This is a current
nomenclature audit, not a reconstruction of the experiment's original annotation.

Only exact equality against `symbol` or `prev_symbol` is admitted. The audit
never uses `alias_symbol`, case folding, suffix removal, fuzzy matching, count
values, response scores or zero fill. A unique candidate needs one HGNC identity,
one source coordinate and no collision with another exact or candidate panel
coordinate.

| Classification among the 409 exact absences | Features |
| --- | ---: |
| Unique nomenclature candidate | 271 |
| No HGNC approved/previous-symbol match | 106 |
| HGNC identity found, corresponding source symbol absent | 29 |
| Multiple HGNC identities | 1 |
| Multiple possible source coordinates | 1 |
| Conflict with another model-panel coordinate | 1 |

For example, `AARS` is a previous symbol for `AARS1` (HGNC:20), and `AARS1`
exists in the Parse source. Three concrete rejection cases prevent automatic
remapping: `QARS` occurs in two HGNC previous-symbol records; `VARS` has two
source-coordinate candidates; and mapping `TIAF1` to `MYO18A` would collide with
the model's existing `MYO18A` feature. Counts are not merged to hide collisions.

**271 candidates do not establish a compatible full panel.** The other 138
features remain unresolved or ambiguous. Moreover, the released source `var`
contains symbols and `n_cells`, without stable gene identifiers. A nomenclature
relationship alone does not prove equivalent genomic intervals, transcript
aggregation or count semantics between reference versions. A future transfer
protocol must explicitly resolve identity and supported features before fitting
or scoring; this audit does not modify the existing exact-name protocol.
The unresolved IFN-beta dose/reagent gate also remains separate.

## Reproduction and retained evidence

```sh
python Tools/Omics/PerturbationPrediction/ParseIFNB/audit_feature_identity.py \
  --panel /path/to/frozen-duration-panel.json \
  --roster /path/to/parse-source-roster.json \
  --hgnc /path/to/hgnc_complete_set.txt.gz \
  --hgnc-sha256 6a1423507780773fbcb373eaef84dba5f9357b36e7219db2aeaae13d8a42d62b \
  --out /path/to/new-audit.json
```

[The retained artifact](evidence/2026-09-11-feature-identity/manifest.json) includes
the complete compressed HGNC table, retrieval headers, original panel and source
roster, every candidate/classification, exact recipe and independent checker.
The checker scans primary rows directly for all 409 missing features and verifies
every candidate and collision; the complete result and repeated execution agree.
The original Parse metadata retains its Parse Biosciences CC BY-NC 4.0 attribution.
No expression matrix, fitted prediction or held-out score is used for matching.
