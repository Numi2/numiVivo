# Prepared molecule to conditional reaction observable — 2026-09-08

The published hydrogen example completed preparation → accepted Metal sampling → verified geometry export → fresh H₃ reaction qualification → conditional first-event probabilities through both the native SDK and public CLI. The two executions produced identical checkpoint bytes and identical complete reaction payloads. This qualifies one finite software route under explicitly supplied models, not general chemical or biological accuracy.

## Exact tested scope

- Runtime and test source: `3d387f95309fd138f9dfa3f451c4d569dfdffe14`.
- Public example runner: `f10a0b1aaf038e0e25d7c9235bb9b5de9d528e05`. Its only difference from the runtime source is preserving JSON media metadata when importing existing campaign artifacts; Sources, Tests and Package.swift are identical.
- Apple M4 Pro Mac mini, 12 CPU cores, 24 GB memory; macOS 26.6 (25G72), Xcode 26.6 (17F113), Swift 6.3.3.
- Complete release product/test-target build passed. The final incremental build finished in 42.97 s after the preceding library/CLI compilation; this is not a clean-build timing benchmark.
- Metal API validation was enabled for the full native fixture and public example. GPU ownership was coordinated with Numi Human; reaction qualification and kinetic checks used CPU FP64.

| Gate | Result |
| --- | --- |
| Host connection, archive/export, mapping and output regressions | 57 executed test functions across 7 suites, 2.366 s. One expensive full-route function explicitly skipped in this host gate. |
| Opt-in native full-route gate | 2 executed functions, 2202.453 s; includes preparation rejection and the complete retained campaign. The preparation function also appears in the host gate. |
| Public preparation-to-observable example | 42 checks passed; fresh reaction and encounter tasks, explicit failure cases and retained outputs. |
| Sampling/export CLI regressions | 87 checks passed using the earlier retained native prefix in a separate working store. |
| Existing general workflow CLI regressions | 61 checks passed. |
| Evidence audit | Both 725-file Git source archives verified; 51 binary/debug/resource entries verified; all 72 native and 75 full-route CLI content objects verified; immutable task/receipt/output bindings checked. |

These are selected gates, not a rerun of every repository test. CI configuration now includes the new host suites and full public example; no hosted CI execution is claimed here.

## Model and observed result

[Published inputs](../../Examples/prepared-reaction/README.md) supply the complete synthetic harmonic H₂ preparation library: masses 1 Da, bond length 0.125 nm, force constant 1000 kJ mol⁻¹ nm⁻², zero charges/Lennard-Jones terms, and sampling at 300 K. The existing preparation owner compiles the system; the fixture does not directly manufacture a classical system.

Two seeds, 17 and 29, execute four eight-step production blocks after four equilibration steps at 1/1024 ps. Replica 0 is selected at accepted step 36 and clock 0.16015625 ps. Its eight production frames in four chunks are verified. The finite prefix remains `budgetExhausted` with unmet sampling criteria. It is an initial geometry source only.

Both H₂ endpoint guesses receive those mapped coordinates. The separately supplied H atom and H₃ saddle guess, full-CI/STO-3G model, destination masses 1.008 Da, temperature 298.15 K, standard state, connectivity settings and assumed transmission remain explicit. Fresh minimum/saddle qualification and four connected-path trials with twelve passing comparisons precede the rate. The connection record reports 231800 descent electronic evaluations for its declared reconstruction; it is not an aggregate count of all verification replays.

- Molecular TST rate: `1.811334723956876e-11 Pa^-1 s^-1`.
- Explicit maintained H₂ partial pressure: 101325 Pa.
- Conditional hazard: `1.835334909049304e-06 s^-1`.
- First-event probability at 0.001 s: `1.835334995092808e-09`.
- Transmission coefficient: 1, explicitly assumed classical no-recrossing TST.

The existing one-state kinetic solver agrees with `1-exp(-hazard*time)` within the declared absolute checks. At the shortest positive observation time, the very small probability rounds to zero; this is not an assertion of a physically zero rate or relative rare-event accuracy. Local sensitivities and Poisson truncation estimates do not supply physical/model uncertainty.

The public CLI reaction stage took 1111.64 s and its encounter stage 738.05 s. These include full source reconstruction and concurrent CPU work. They establish no throughput, scaling or speedup claim.

## Failure and provenance checks

The host and native gates reject missing force-field parameters, incomplete/duplicate/out-of-range or element-incompatible atom mappings, substituted replica evidence, unsupported periodic/environment transfers, mismatched contexts, corrupted source/output/receipt records and cancelled publication. The full native campaign rejects an altered numeric rate, changed temperature and changed qualified-component identity. Public commands additionally reject missing reservoirs, incompatible units, output collisions and artifact-store aliases.

Two throwing-expression build issues and an example import media-type mismatch were corrected before the passing runs. Their failed attempts remain in the raw evidence. No numerical qualification tolerance was relaxed.

The saved native campaign is `c05fef93f18d332eda51cfff03afb8795b7eb9b7a0966d111b9838e1d4a3dd04`. The selected checkpoint is `a33a48188cf07c27065cd580b3d300060939962e381bf9d1daf83928fcbf6cce`. The identical native/CLI reaction payload is `fa9e3020cac40b2f5939032d973ac9cf3668ec1a25b3ffdb9e2bd08e4c9c8379`.

The CLI executable SHA256 is `ba39e032dc0c04c860d142354055add29f6b1ccd04d6881308714666031f4c8a`; the loaded test image SHA256 is `b544c4d654b94c3b4b9fe6055b1921ee0abe6e8b2ae2061a8a93617886dbd50c`. SwiftPM launches tests through a separate helper; the campaign's `executableSHA256` identifies that launcher. The audit separately verifies the test image, its sampled Mach-O UUID, and shader resources against the frozen build.

Compact machine-readable observations are [retained here](2026-09-08-prepared-reaction-observations.json). Raw evidence is mirrored at `/Users/{home,n}/numivivo-prepared-reaction-evidence-20260908`; frozen binaries remain on the Mac mini. `raw-manifest.json` covers 1050 files, excluding binary directories, Git bundles, itself and the audit console log; SHA256 `3a19ad6c19943d34ce29c836dde4d6c19d3fff3f5248d37c8a5e3c40c5335391`. Its entries were also rehashed in the local mirror.

The [interface contract](../Design/PREPARED_REACTION_CAMPAIGN.md) preserves distinguished mapped reactant roles, source/destination model differences and maintained-reservoir assumptions. The result is not an ensemble-averaged rate, experimental hydrogen rate, finite-pool kinetic model, reverse-reaction network, protein barrier, target occupancy or biological efficacy claim.
