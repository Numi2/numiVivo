# Fresh real-Hagai CLI lifecycle

The updated CLI completes fresh PCA from the original H5AD, publishes the MNN
bundle, and reconstructs/verifies it successfully. This uses all 13,863 original
Hagai cells and the retained plans, with the new executable identity throughout.

| Step | Exit | Seconds |
| --- | ---: | ---: |
| H5AD to PCA | 0 | 7.639 |
| PCA to MNN bundle | 0 | 10.111 |
| Full MNN bundle verification | 0 | 10.289 |

Fresh PCA scores, MNN scores, anchors and report match the original artifacts
byte-for-byte. These wall times include CLI bundle handling; they are a single
execution, not a comparative performance benchmark or new biological evaluation.

Three rejection checks return exit 65: an existing destination, a matching-work
budget of one, and an anchor byte changed in a separate copied bundle. Existing
good bundle files retain their hashes. Budget rejection leaves neither its target
nor a staging entry. The altered bundle remains retained as failure evidence.
These are scoped checks, not exhaustive adversarial or concurrent-write coverage.

The first launcher failed Python syntax parsing before any CLI execution. Its
source and error log are retained alongside the corrected driver. The successful
run is recorded in terminal.json; qualification.json records rejection results,
original/fresh hashes and the complete good MNN bundle file hashes. The archive
also includes exact plans, source/CLI identities and command logs. Large source,
binary and bundle files remain in the hash-bound external workspace identified
by manifest.json. No full package or GUI app build is claimed.

Run `python3 verify.py` for retained evidence consistency checks. Fresh execution
requires the recorded source, CLI and HDF5 runtime paths: driver.py constructs new
pca/mnn outputs and rejections.py requires a new tampered destination. Use a fresh
workspace with the recorded protocol and plans rather than overwriting evidence.
