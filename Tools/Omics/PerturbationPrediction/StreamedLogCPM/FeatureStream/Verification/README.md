# Native feature-stream reconstruction

```sh
numivivo-omics singlecell-logcpm-feature-verify bundle < records.bin
```

The verifier checks receipt format, bounded metadata, dimensions, payload size,
metadata hashes, and every stored Float64 mean against streamed recomputation.
It also compares complete summary axes, cell counts, zero counts and the input
and output digests. The means file is consumed a feature at a time; the entire
result matrix is not loaded. Truncated and trailing data are rejected.

Complete Norman reconstruction passes: 361,582,621 input records and all
7,985,478 means, with 41,959,424 bytes peak RSS and 9.20 seconds wall time in the
recorded run. [norman.log](norman.log) retains the measurements.
[verification.json](verification.json) binds the binary and test source.
[check.py](check.py) verifies the complete bundle and rejects changed means,
changed denominators, shortened/extended payloads even with recomputed hashes,
symlinked means and malformed input. Tests use separate copies of the small
bundle; the original evidence remains unchanged.

This verifies numerical consistency with the supplied plan and count stream.
It does not authenticate a maliciously replaced plan or establish source
provenance independently: retain a trusted receipt/plan identity separately.
It is not a concurrent-mutation snapshot guarantee or biological validation.
Original execution data remains at
`/Users/n/numivivo-logcpm-feature-verify-20260912` with separately retained Norman
source and feature-stream bundles required for reproduction.
