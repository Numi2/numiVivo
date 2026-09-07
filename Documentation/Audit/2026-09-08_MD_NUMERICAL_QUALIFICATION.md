# Classical MD numerical qualification — 8 September 2026

This campaign tests declared classical molecular-dynamics cases on the physical Apple M4 Pro Mac mini. Independent tests first exposed two defects; the numerical implementation was repaired without widening their limits, changing their seeds or replacing their fixtures. The resulting numerical profile is `numivivo.org/md-metal-numerics/v3`.

## Source and execution identity

The [frozen qualification plan](2026-09-07_MD_QUALIFICATION_PLAN.md) was committed before the first hardware measurements. Initial source `82c2e6c7733354129349b7c658137ee23193fd3a` uses numerical profile v2. It differs from the first test commit only by explicit intermediate coordinate types required by Swift's type checker.

The repaired candidate `be4abc6adac272aab335491263087c9a718b165f` passed a complete release build including all test targets in 220.07 seconds and 42 focused test functions across nine suites in 6.719 seconds. That run covered the new scientific gates, ordinary MD restart, Metal sampler identities, QM/MM archive identities, multipole PME, platform workflows and Apple execution. It used the Metal validation layer and a fresh artifact directory.

Final source `17406c0d59a8f56bba2fbacade6c74424408121d` passed the complete release/all-test-target build in **205.65 seconds**, then **261 executed native test functions across 59 suites**, **57 production workflow CLI checks**, both target-panel tasks and six fresh archive CLI invocations. The five new resource-admission host tests also passed separately before the full run. Swift Testing reports 262 tests including the single explicit opt-in long-archive skip; that campaign's earlier enabled execution remains separately recorded.

Full-suite wall time was **917.222 seconds**, including a recorded **160-second SIGSTOP/SIGCONT pause** during the confirmed CPU-only H3 Gaussian calculation to reserve another Numi task's GPU window. This duration is not a throughput result. The process resumed and completed with exit status zero.

| Final executable | SHA-256 |
| --- | --- |
| Production CLI | `cc154a9b2b0040bf3d51266a95e319421f30259a8e2fd768ac1e13593c3fb22a` |
| Test executable | `34d69e2926957465ef61b066c979a9da40130f16b63fd2d30225d24b492e43f8` |

All 23 shader-resource hashes match the first passing v3 candidate. The final checkout remained clean and both executable hashes were unchanged after the gates. Publication adds only documentation and compact evidence to this tested source.

Hardware: Apple M4 Pro, unified memory, macOS 26.6 (25G72), Xcode 26.6 and Swift 6.3.3. Reported build/test durations describe these executions and are not performance comparisons. Other coordinated work used the same machine outside the reserved native windows.

## Charged PME: independent direct-Ewald comparison

Nine four-particle fixtures are evaluated at 16³, 32³ and 64³ mesh dimensions: neutral, displaced, subcell translated, primary-pair excluded, half-scaled, charged with neutralizing background, charged and half-scaled, whole-lattice translated, and affinely changed cell. Coordinates and charges are exactly represented by the fixture. An independent FP64 direct-Ewald calculation supplies energies and forces; reciprocal truncation and real-space cutoff are separately convergence-checked. Each probe verifies that the accepted MD checkpoint remains unchanged.

Fine-grid limits were fixed at 0.05 kJ/mol absolute energy error, 0.001 normalized RMS force error and 0.003 normalized maximum component error. Force normalization is the reference force-component RMS with a 1 kJ/mol/nm floor. Coarse-to-fine aggregate checks and lattice-translation checks are separate gates.

| Measurement | Initial v2 | Final v3 |
| --- | ---: | ---: |
| Excluded-primary normalized RMS force error | 0.004289691482 | 0.00002717347044 |
| Excluded-primary normalized maximum component error | 0.006705767202 | 0.00005933891940 |
| Excluded-primary signed energy error, kJ/mol | -0.004447019071 | -0.0001211523714 |
| Fine-grid cases meeting all fixed limits | 8/9 | 9/9 |

Excluding the strong primary pair reduced the reference force scale and exposed the remaining cubic-mesh error. This was a mesh-accuracy defect, not evidence that the primary-pair exception convention was wrong. Classical charge spread and energy-derivative gather now use sixth-order cardinal B-splines with the matching discrete nodal deconvolution. The existing multipole basis is shared; its influence arithmetic remains separate and is covered by its own regressions.

The compact PME, thermal and NPT scientific observations are identical between the first passing v3 candidate and final source. This records these two executions, not a general cross-device determinism guarantee.

The 16³-grid aggregate energy RMS error changes from 1.04737 to 1.14606 kJ/mol in the first candidate comparison. Improvement at every mesh resolution is not claimed.

This is static force/energy evidence for the listed cells and configurations. It does not establish interacting NPT accuracy, every mesh setting, a large-system error bound or speed. The sixth-order stencil touches 216 nodes instead of 64. Requested PME tolerance still controls Ewald planning rather than certifying total force error.

## Thermal distributions

The unchanged thermal tests pass on both initial v2 and repaired candidate v3. The finite-time free Langevin test uses 1,024 particles and eight seeds, with 8,192 samples separately evaluated at each of 1, 8 and 32 steps. The expected velocity means and covariances use the analytical finite-time Ornstein–Uhlenbeck distribution.

Fixed-geometry thermalization tests use 128 unequal-mass dimers and 128 connected nonlinear triatomics over 32 seeds, yielding 4,096 molecular samples per geometry. A dense, independent mass-metric projection supplies the constrained Gaussian covariance. Total, center-of-mass and rotational kinetic statistics are checked against the appropriate chi-square distributions; global center-of-mass statistics use the 32 independent seed groups.

The nominal family false-rejection budget under the analytical Gaussian sampling model is 10⁻⁴, the concentration-bound budget is 512 comparisons, and the tests execute 140 statistical comparisons (27 free-particle, 43 dimer and 70 triatomic). The fixed FP32 allowance and geometry/tangent tolerances were not changed. Physical thermal samples and all 13 thermal report files are identical between the initial v2 and first v3 candidate; this retains prior passing evidence rather than demonstrating a thermal-method improvement. Kinetic variances are diagnostic, and the 32-seed global center-of-mass check is only a coarse retained-mode guard. Passing these checks does not establish constrained-trajectory equilibrium or a general random-number quality certificate.

## Molecular NPT state transitions

The conditional transition test uses eight monatomic particles and eight unequal-mass rigid dimers, including molecules crossing the periodic boundary. Charges, Lennard-Jones energies, initial velocities and thermostat friction are zero. Both systems use the same fixed proposal seed, temperature and pressure. Independent SI-unit arithmetic checks pressure work and the molecular `(Nmol + 1) log(V'/V)` acceptance term. Proposed molecular centers and retained geometry are reconstructed independently.

The initial dimer run encountered pressure-proposal geometry rejections, then rejected dynamics at accepted step 11, time 1.75 ps, with constraint status 8. Its raw rejection certificate and accepted before/after checkpoints are retained; the rejected candidate buffers are not exposed by this test. Repeated coordinate reconstruction through the cell and reciprocal cell was not an FP32 identity; even zero-correction projection passes could move stationary atoms and create spurious constraint impulses.

The shared periodic helper now uses reciprocal coordinates to select lattice images and subtracts those vectors from the original physical coordinates. Volume proposals pack the actual proposed cell and use a fused molecular-center displacement. The test's 32 constraint iterations and 10⁻⁶ relative tolerance are unchanged.

Final v3 completes all 128 proposals for each system. Each has 106 accepted and 22 rejected pressure proposals, with exactly matching scalar proposal/acceptance traces. The largest dimer bond-length error is 8.9225×10⁻⁸ nm, below the fixed 6.25×10⁻⁷ nm gate. Maximum dimer speed is 2.3842×10⁻⁷ nm/ps, below 5×10⁻⁵; maximum analytical score error is 6.6614×10⁻¹⁵, below 10⁻⁹.

These zero-force, zero-friction tests qualify proposal scoring and accepted-state boundaries. They do not qualify equilibrium volume distributions, interacting systems, pressure estimates or thermostat/barostat mixing.

## Continuation, planning and workflow contracts

The [v3 numerical contract](../Design/MD_NUMERICAL_PROFILE_V3.md) separates algorithm identity from physical system/configuration identity. Older MD checkpoints remain decodable historical data. Direct restore, stage transfer and protocol continuation reject a missing or old numerical contract. Explicitly importing physical state starts a new numerical run; modifying an old contract field is not migration.

Metal constant-pH and both nuclear potential identity constructors now bind this contract. NCMC, ring-polymer checkpoints and reactive labels inherit the changed endpoint/evaluator identity. QM/MM free-energy requests and archive cursors bind the contract even at completed-window boundaries with no active MD checkpoint. Historical static-trace analysis remains separate. The [sampler identity audit](2026-09-07_MD_SAMPLING_IDENTITY.md) documents the source boundaries and host regression controls.

The final source also forwards explicit PME dimensions to nuclear and workflow reservations and MD command planning, matching the allocating PME engine. Point caps and byte budgets therefore describe the selected mesh. Very fine ignored spacing does not cause a fixed-grid request to derive a huge mesh. Direct implicit planning rejects unrepresentable axes and grid-count overflow without integer conversion traps. Five passing host tests cover both workflow start/continue admission, resource boundaries, actual command construction and direct planner boundaries.

CLI workflow-template help is generated from the actual template registry. A host test constructs and validates every advertised recipe, including `md-electronic-analysis`.

## Reproduction and retained evidence

The final native checkout is `/Users/n/numivivo-completion-final-20260907`; raw evidence is retained separately in `/Users/n/numivivo-md-qualification-evidence-20260907`. Local compact observations and raw JSON mirrors use `/Users/home/numivivo-md-qualification-evidence-20260907`. Source archives, build/test logs, executable hashes and individual shader-resource hashes bind the measurements. Previous executable bundles are retained with their matching shader resources so they cannot silently load a later source resource through SwiftPM's fallback path.

The original `tests-2.log`, `native-2` and `initial-observations.json` retain negative v2 results. `tests-v3.log`, `native-v3`, `observations-v3.json` and `raw-manifest-v3.json` retain the first passing candidate. The final run is in `build-final.log`, `tests-final.log`, `native-final`, `observations-final.json` and `raw-manifest-final.json`, with 457 raw JSON files. Its manifest SHA-256 is `7a4bec6f2b6098d59b5218fe820f18817e5c29e51a9aa44dfe4ec186d2597a28`. `cli-final`, `archive-cli-final` and `target-panel-final.json` retain the CLI checks. The [published compact observations](2026-09-08-md-native-observations.json) bind the original failure, first passing candidate, final native result, CLI checks, pause and file identities. Compact published observations are extraction records; raw observations and a successful native test log are both required to establish the result.

```sh
swift build -c release --build-tests --jobs 2 -Xswiftc -enable-testing
MTL_DEBUG_LAYER=1 NUMIVIVO_TEST_ARTIFACTS=/absolute/fresh/evidence \
  swift test -c release --skip-build --no-parallel
python3 Tools/Platform/check_workflow_cli.py .build/release/numivivo /absolute/fresh/cli
```

The platform workflow requires the six new scientific/identity/resource suites and records shader-resource hashes alongside its source and executable identities. YAML and shell syntax were checked; no hosted CI execution of this publication is claimed here.

The opt-in 100,001-chunk archive campaign remains a separate [measured milestone](2026-09-07_ARCHIVE_STREAMING.md). The current campaign does not regenerate its fixture or reinterpret its original source-bound observations.

## Remaining acceptance boundaries

This milestone does not complete NumiVivo. Interacting and longer constrained ensembles, wider force-provider combinations, useful system scale and performance, a published prepared-system-through-observable scientific route, broader declared runtime semantics and release qualification remain in the [completion roadmap](../COMPLETION_ROADMAP.md). No biological prediction or general chemical-accuracy claim follows from these small numerical fixtures.
