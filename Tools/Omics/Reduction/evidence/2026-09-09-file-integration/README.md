# File integration qualification archive

See `../../FILE_INTEGRATION.md` for scope, commands and measured failures.
`qualification.json` separates storage/numerical checks from biological outcomes.
Kang NK-cell preservation and its original global erasure control failed; the
supplementary type-condition erasure control does not replace those outcomes.

`archive-manifest.json` records stored and uncompressed file hashes. For a file
marked `gzip`, decompress it and remove only the final `.gz` suffix to recover
the original artifact bytes. Already-compressed legacy oracles and NPZ files
remain unchanged. Native binaries and real H5AD files are external; exact paths,
hashes, source plans, builds and runtime identities are retained.

`external-inputs.json` maps every omitted native `original.h5ad` to an exact
retained canonical input. Restore those bytes before native replay. Kang's raw
public download hash and exact count/identity preparation checks are in
`prepared/preparation.json`; the complete public source has 24,673 cells.
Hagai uses the prior source-supported paired donor plan, not its older unqualified
graph mapping. Source study SDRF and supplement provenance is retained alongside
this archive.

`harness-history/fixtures` retains the initial synthetic query-overlap rejection.
The final 54-check suite uses distinct query IDs; synthetic fixtures establish
software behavior only. The 22-check original-PCA graph regression is separate.
`real-evidence/commands.json` and `downstream-kang/commands.json` are measured on
Mac mini. Python biological references ran on the laptop; timings are not a
CPU/scverse performance comparison.
