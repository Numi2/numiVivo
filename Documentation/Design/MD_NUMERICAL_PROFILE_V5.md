# Classical MD numerical profile v5

`numivivo.org/md-metal-numerics/v5` adds opt-in compensated positions for fixed-cell
classical NVE/NVT. Omitted `positionPrecision` retains FP32 behavior. The v4 torsion
force repair is retained. Old numerical contracts cannot silently resume.

Each position uses high and low FP32 words. Drift, distance projection, velocity
projection and constraint impulses use compensated differences and additions.
Lattice subtraction preserves the correction before forming the local distance.
Fast math is disabled. Force kernels and velocities remain FP32; this is not a
double-precision force implementation. Five additional float4 buffers require
80 bytes per particle, each checked against Metal allocation limits.

Both position words share the existing arena's transactional A/B owner. Snapshot
coordinates are their Double sum. Checkpoints additionally retain both exact
words, including corrections too small to survive that Double sum. Restore
checks precision, shapes, exact FP32 words and consistency with the coordinate
view. Stage transitions cannot implicitly change precision. Failed preparation,
steps and cancellation preserve the accepted pair and clock.

Initial support excludes moving cells, virtual/Drude sites, external candidate
force providers and minimization. FP32 trajectory archives reject compensated
sampling before execution; direct snapshots and JSON benchmark reports retain
the higher-precision coordinate view. These are explicit unsupported paths,
not conversions or claims of broader qualification.

Required evidence: sub-ULP drifts, exact restart of both words, large-coordinate
and periodic-boundary constraints, rollback, unsupported-path admission, and
the unchanged independent reference panel with explicit precision configuration.
Scientific ensembles, free energies and speed leadership require separate runs.
Native v5 validation is pending at the introduction of this document.
