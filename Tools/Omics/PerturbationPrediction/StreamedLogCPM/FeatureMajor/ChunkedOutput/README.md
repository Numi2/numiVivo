# Bounded numeric JSON output

The existing report schema is written in chunks of at most 4,096 means, avoiding
whole-report JSON encoding. Feature/group metadata and result arrays remain
resident. The CLI writes a private staging file, closes and synchronizes it,
then creates the final hard link without replacing an existing destination.
Failure removes staging; no partial report is published.

The complete Norman CLI result remains numerically exact across all 7,985,478
means. Peak RSS falls from 770,605,056 to 237,813,760 bytes (69.1% lower).
Measured runtime is 7.27 seconds versus 7.14 previously; no speedup is claimed.
These are single-run whole-process measurements on the same dataset and host.
Standalone execution also preserves the report exactly at 230,768,640 bytes RSS.

[receipt.json](receipt.json) binds the new binary, source, input and output.
[native.log](native.log) retains process measurements. [check_cli.py](check_cli.py)
checks full output equality, ordering rejection, existing-output preservation,
and an injected file-size-limited write error. [transaction-check.json](transaction-check.json)
confirms no final report or staging residue remains after that error. It does
not simulate a power failure or every concurrent filesystem race.

Large condition-by-feature arrays and metadata still remain resident; this is
not a fully disk-backed result representation or million-cell qualification.
Complete outputs remain on the execution host under
`/Users/n/numivivo-logcpm-output-20260912`; their hashes are in the receipt.
