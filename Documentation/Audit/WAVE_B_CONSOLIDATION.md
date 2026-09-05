# Wave B execution consolidation

This pass continues from e1986144f181321ddd17821dc6a8f9051767c0b7. The source-level findings below are not compiler or numerical qualification results.

## Repaired execution boundaries

1. Snapshot/checkpoint: a reservation now covers readback and snapshot metadata. Checkpoints take clock, accepted step and cell from that immutable snapshot, never by rereading mutable actor fields after an await. A new sample API obtains positions and optional observables from one reserved generation. Observers do not reconstruct virtual sites by writing into accepted state.
2. NPT transaction: deterministic/stochastic integration and pressure proposals are completed on uncommitted candidate buffers. Cell, PME and neighbor plans are local to the candidate. The final buffer flip, phase assignment and clock publication happen without suspension, only after all fallible work and cancellation checks. A rejected pressure proposal restores the dynamics candidate; it does not roll the integration back or publish twice. Allocation, compilation and device failures abort the step rather than being counted as statistical Metropolis rejections. GPU command failure stops this runtime instance.
3. Molecule geometry: the barostat now derives components from bonds, constraints, angle/torsion connectivity and virtual parents. A topology spanning tree unwraps long molecules edge by edge. Non-tree edges detect periodic winding. Scaling preserves molecular internal coordinates, not atom-wise fractional geometry.
4. Constrained drift: position-projection corrections now contribute the corresponding velocity increment using the actual full/half drift interval. This repairs an omitted constraint impulse; it does not by itself establish symplecticity or ensemble correctness of the bounded Jacobi projector.
5. Minimization: the descent direction is projected into the distance-constraint tangent space using an identity metric over massive particles. Its tangent residual is checked. Convergence ignores dependent virtual-site forces while retaining their energy. First reported energy is after explicit initial geometry projection; MD time/RNG step do not advance during minimization.
6. Virtual force race: threads with no parent incidence return without writing. Virtual force slots can therefore be read by physical parents without an overlapping add-zero write from the virtual thread.
7. Observables: fixed-tree GPU reductions return scalar potential and kinetic energy instead of downloading one kinetic scalar per particle. Scratch buffers are retained rather than repeatedly allocated per observation/minimization iteration.
8. Numerical identity: checkpoints record md-metal-numerics/v2 independently of system/configuration hashes. Older checkpoints still decode for inspection but cannot silently claim continued execution under changed numerical semantics. Stage transitions have a distinct explicit reconfiguration record. Thermal initialization is explicit and seed-namespaced, never an implicit checkpoint-open side effect.
9. Package manifest: swiftLanguageModes precedes cxxLanguageStandard according to the PackageDescription initializer signature. All new kernel resources and entry points are registered together.

## Scope corrections

The archived structure format accepts general triclinic cells, but the current MD nearest-image kernels round fractional coordinates independently. That rule does not implement the general closest-lattice-vector problem. Runtime preflight now admits orthogonal periodic cells only and explains why skew-cell execution is blocked. This is an intentional restriction on an overstated execution capability, not a regression in structure interchange.

Neighbor lists are rebuilt at both force positions, preventing sampled displacements from turning into a trajectory-selection rule. Spatial binning remains the periodic builder. The configured rebuild interval is retained for compatibility but is not yet a guarantee of amortized reuse. Grid occupancy/neighbor overflow still abort execution. Independent-distance-constraint count is assumed in reported degrees of freedom; general constraint-rank estimation remains work.

## Still requires development/qualification

The constrained Langevin/integrator scheme requires independent distribution and energy-conservation checks. PME charge deposition uses floating-point atomic updates, so fixed seeds do not imply bitwise-identical cross-device or cross-launch trajectories. Mesh spacing and real-space truncation controls are heuristics, not measured force-error bounds. Force-field exception conventions, reciprocal discretization, pressure sensitivity and long-range dispersion must be compared with independent references. Resource preflight is not yet a single shared budget across every NumiLab runtime. Failed allocations are safe transactionally but not an optimized scheduler.

No Apple package build, shader compilation, MD simulation, statistical test or benchmark was run during this pass. Numerical-contract versioning is a compatibility mechanism, not evidence of physical validity.

## Primary references consulted

- Swift concurrency and actor reentrancy: https://docs.swift.org/swift-book/documentation/the-swift-programming-language/concurrency/
- Apple command resource lifetime: https://developer.apple.com/documentation/metal/setting-up-a-command-structure
- OpenMM force/PME conventions: https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html
- OpenMM integrator reference: https://docs.openmm.org/latest/userguide/theory/04_integrators.html
- PackageDescription Package initializer: https://developer.apple.com/documentation/packagedescription/package
