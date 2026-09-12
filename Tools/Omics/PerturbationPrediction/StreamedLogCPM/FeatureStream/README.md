# Feature-at-a-time log-CPM bundles

```sh
numivivo-omics singlecell-logcpm-feature-stream --plan plan.json \
  --stream-sha256 SHA256 --output new-bundle < records.bin
```

The plan requires `ordering: "featureMajor"`, featureIDs, groupIDs, rowGroups
and full-axis UInt64 rowTotals. Records must be strictly feature-then-row ordered,
positive and unique. Every declared feature is emitted, including all-zero ones.
Zero-count cells remain in group denominators. The library retains axes, row
assignments/totals and one feature's group sums/corrections, not a resident
condition-by-feature result matrix. Existing admission limits remain unchanged.

The new bundle contains:

- `means.bin`: complete feature-major, group-minor little-endian Float64 means,
  with no record headers; shape is features by groups.
- `summary.json`: method, exact feature/group axes, cell and zero-cell counts.
- `plan.json`: the supplied plan bytes.
- `receipt.json`: format `feature-major-group-f64-le-v1`, dimensions, payload
  size and SHA256 values for input counts, plan, summary and means.

Feature writes are provisional until all row totals and the complete input hash
pass. Private staging is removed on failure; a successful bundle is moved to a
new destination. Existing output is preserved. Keep the executing implementation
identity separately: this bundle receipt binds data, not the executable itself.
The qualification receipt below additionally binds the tested binary and sources.

The actual scoped CLI processes every Norman record: 111,445 cells, 33,694
features, 237 conditions and 361,582,621 positive records. Its 7,985,478 means
are bitwise identical to the standalone feature stream and numerically exact
to the prior resident native result. That prior result passed independent NumPy
verification with maximum absolute error 9.06e-14. All axes and cell counts match.

Peak CLI RSS is 44,105,728 bytes, versus 237,813,760 for chunked JSON with resident
means and 770,605,056 for whole-report encoding. Runtime is 9.35 seconds versus
7.27 and 7.14 seconds respectively: this reduces memory but is slower in the
recorded single-run measurements. These are different output representations,
not an equivalent-format speedup comparison or million-cell qualification.

[verification.json](verification.json), [bundle-receipt.json](bundle-receipt.json)
and [norman.log](norman.log) retain source-bound evidence. [check_cli.py](check_cli.py)
verifies the full bundle, zero features/cells, duplicate/decreasing coordinates,
missing/excess totals, truncation, wrong hash, overwrite protection and an
injected file-write failure. All eight rejection cases leave no final bundle
or staging residue. Cancellation callbacks are present; process-kill and
power-failure durability were not qualified in this run.

Full bundle and original runtime remain under
`/Users/n/numivivo-logcpm-feature-output-20260912/cli-check/norman` and its parent
study. Large payloads are not embedded here; hashes bind them. Reproduction needs
the separately retained Norman count store and a new output directory. No
biological prediction or independent-replication evidence follows from this work.
