# Wave B archive and preparation qualification

This small executable exercises the actual public trajectory writer/reader,
protocol factory and velocity initializer. It does not substitute a Python
implementation for the Swift code.

On an Apple-silicon development machine after resolving package build diagnostics:

```sh
swift run --package-path Examples/md-archive-qualification MDArchiveQualification
```

The checks cover packed coordinate/velocity round trips, frame and file digests,
partial-read verification semantics, no-clobber publication, frame corruption,
truncated footers, trailing data, zero-frame archives, rejection of lossy Double
input, deterministic thermal preparation, COM and constraint projection, distinct
stage random namespaces, and one-time initialization in the default protocol.

This source harness has not been executed in the development continuation that
added it. Its presence is not a passed-test result. It imports NumiVivoKit and
therefore first requires the Apple package to compile. It does not run an MD
trajectory, instantiate Metal force pipelines, qualify PME accuracy, or establish
NVT/NPT ensemble correctness. Those remain independent qualification work.
