# Per-chunk temporary-buffer lifetime

An explicit autorelease pool around each bounded read/consume step releases
Foundation temporaries during processing. The carry buffer remains outside
the pool, preserving records split across reads. Non-Objective-C platforms
invoke the same body directly.

On the complete 2,711-cell real matrix, alternating three runs of each CLI gives
peak RSS 109,395,968 bytes before and 25,575,424 bytes after (76.6% lower).
Every normalized result value is exact. Wall times are 0.11 seconds before and
0.10–0.11 seconds after; this does not establish a meaningful speedup.
These are warm-cache process-level measurements on one input and host.

[cli-memory.json](cli-memory.json) binds both binaries and the input and retains
all six runs. [receipt.json](receipt.json) records the rebuilt CLI's exact real
matrix check, zero-cell retention, and rejection of truncated, duplicate,
wrong-total, out-of-axis, empty, wrong-hash and overwrite cases.
The previous standalone comparison also retained an initial 0.51-second pooled
run; it is not used as evidence of speed improvement. Full-app and million-cell
memory qualification remain separate.
