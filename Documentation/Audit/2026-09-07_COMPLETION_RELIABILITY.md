# Completion reliability validation — 7 September 2026

The reliability changes passed native verification: the baseline full suite passed 224 tests in 52 suites, the final source passed all 22 targeted regressions, and the production workflow CLI passed 57 checks. These are scoped execution results, not product-wide scientific qualification.

## Source and environment

- Final tested Swift product/test source: `d4c5216696618e6bf2faaeed76faeea93cfcef63`.
- System-Python CLI-checker compatibility repair: `63e7aa40675c4dd3d011727c2f71b7f5b49cb1ef`; `Sources`, `Tests` and `Package.swift` are unchanged from the tested Swift source.
- Full-suite baseline: `c8fcddc6429684c3e9a0769d7bf76152e3daff33`.
- Physical host: Apple M4 Pro Mac mini, reached through `ssh macmini`.
- Toolchain: Xcode 26.6, Swift 6.3.3; macOS 26.6, build 25G72.
- The two sources have separate clean checkouts, source archives and binaries. No other checkout was used to build these results.
- Raw logs, native JSON, CLI stores, source archives and binary hashes are retained outside Git in `numivivo-completion-evidence-20260907` on the Mac mini and coordinating MacBook.

The final native source adds the failed-stage-outcome publication repair and its native regression to the full-suite baseline. It does not change electronic, nuclear, stochastic or MBAR numerical methods. The CLI-checker repair uses `os.link` for the same hard-link alias check on system Python 3.9. The following documentation commit does not change executable source.

## Implemented behavior

### Shared workflow cancellation

Each caller owns a separate wait on a shared task generation. Cancelling one caller releases that caller promptly while other callers keep the calculation. Cancelling the final caller signals the numerical task and retains the generation until it drains. A replacement waits for that retirement before starting; it neither overlaps abandoned work nor inherits its cancellation. A recipe executor retains this ownership across successive runs.

Six regressions exercise cancellation, concurrent callers, draining and replacement, pre-cancelled invocation, cache behavior and recipe report behavior. Gated operations retain bounded completion observers and explicit cleanup. Cancellation is cooperative; it does not promise preemption of an already submitted GPU command or rollback of an atomic write that has begun.

### MD recovery and stage gates

Resume validates the plan, initial/entry/current states, completed report sequence, recomputed stage transfers and output boundaries before constructing the runtime. Required minimization convergence cannot be bypassed by rewriting a report and recomputing its hash. Clock validation allows bounded FP64 accumulation error without replaying billions of additions.

An accepted checkpoint can recover a missing sample or observation at its current boundary. Recovery publishes it before advancing and does not duplicate an existing frame. Older missing outputs cannot be reconstructed from that checkpoint and remain errors.

A returned terminal outcome owns the runner until its checkpoint, report and terminal cursor are published. A failed report write cannot grant an unsuccessful minimizer a fresh iteration budget using its unpublished advanced geometry. A repeated write failure retains the prior durable cursor and an explicit persistence diagnostic.

Nine regressions cover rehashed corruption, stage transfer, clock rounding, legitimate blocked minimizer arithmetic, native multi-stage execution, uninterrupted-versus-resumed coordinates and velocities, and actual observation/report write failures. The write failures use exact test-owned empty directories at otherwise absent content addresses; no production failure injection hook is added.

### Dedicated hybrid execution

Seven regressions exercise mixed exact SSA, tau-leap and RK2 authority, conservation, all UInt32 publication bits, exact-event dispatch chunk invariance, checkpoint continuation, rejected restore atomicity, prepared reservations, count overflow and exact-work exhaustion. Stochastic birth/death cases use 8,192 lanes with fixed seeds and six analytical standard-error bounds declared before execution.

These checks run the dedicated Metal backend. Missing suitable hardware is a failure, not a successful skip. They do not establish arbitrary tau-leap accuracy, rare-event tails, every propensity law, spatial execution, live authority migration or cross-device bitwise identity.

### Complete build

A clean full build exposed a pre-existing Swift type-checking timeout in the nested MBAR covariance expression. Explicit typed loops retain the same symmetric covariance calculation while allowing the complete release and test targets to compile.

## Commands and results

| Source | Command | Result |
| --- | --- | --- |
| `d4c5216` | `swift build -c release --build-tests --jobs 2 -Xswiftc -enable-testing` | Passed; 236.15 s. |
| `c8fcddc` | `MTL_DEBUG_LAYER=1 swift test -c release --skip-build --no-parallel` | 224 tests / 52 suites passed; 762.463 s. |
| `d4c5216` | `MTL_DEBUG_LAYER=1 swift test -c release --skip-build --no-parallel --filter 'WorkflowCancellationTests\|MDProtocolResumeTests\|HybridRuntimeTests'` | 22 tests / 3 suites passed; 0.849 s. |
| `63e7aa4` checker / `d4c5216` binary | `python3 Tools/Platform/check_workflow_cli.py .build/release/numivivo OUT` | 57 real CLI, planning, cache, export and failure-preservation checks passed on system Python 3.9.6. |
| `d4c5216` | `numivivo workflow-plan` and `workflow-run Examples/workflows/target-panel.json` with an external store/output | Both nodes succeeded; `allTasksSucceeded: true`. |

`NUMIVIVO_TEST_ARTIFACTS` pointed to external evidence directories for each native run. An initial system-Python 3.9 invocation stopped at its unavailable `Path.hardlink_to` API. That aborted invocation is retained separately and is not counted as a pass. Python 3.13 passed the original checker; the compatibility repair then passed all 57 checks again on the default system Python in another fresh directory.

The final CLI binary SHA-256 is `8dbef380e14ffd7e5d4c9d87c0b93231f40808adecb5405fce0a6ce219334984`. The complete 224-test run belongs to `c8fcddc`; it is not relabeled as a full-suite run of the later finalization change. That change received a clean complete build and the 22 affected/native regressions on `d4c5216`.

Build and test elapsed times are observations from this validation session, not throughput benchmarks. The full-suite calculation includes an expensive independently reconstructed mapped H3 reaction path. Sampling confirmed active Gaussian integral evaluation while buffered console output was unchanged; the calculation was not cancelled or weakened.

## Stochastic observations

| Process | Expected mean | Observed mean | Expected variance | Observed variance |
| --- | ---: | ---: | ---: | ---: |
| Exact SSA Poisson birth, mean 4 | 4 | 4.009766 | 4 | 3.950826 |
| Tau Poisson birth, mean 4 | 4 | 3.997070 | 4 | 3.907207 |
| Tau Poisson birth, mean 64 | 64 | 63.861450 | 64 | 63.512605 |
| Exact first-order survivors, 40 initial | 24.261226 | 24.221436 | 9.546049 | 9.482886 |

All four cases met their predeclared bounds. Full-precision observations and tolerances are retained in [Poisson moments](2026-09-07-hybrid-poisson-moments.json) and [survivor moments](2026-09-07-hybrid-exact-death-moments.json). These statistical checks support the stated fixtures only; they are not performance comparisons.

## Remaining qualification

The [completion roadmap](../COMPLETION_ROADMAP.md) remains active. In particular, this record does not establish the complete MD force/ensemble/PME matrix, biological-scale accuracy, useful maximum system size, performance against independent tools, or a fully qualified molecular-to-biological research campaign. Nuclear and surrogate qualification remains a separate workstream with its own evidence.

Trajectory restart checks validate the link sequence and final coordinate payload; they do not exhaustively reread older coordinate and observation payloads. The current materialized index also limits resume to 100,000 chunks although the writer permits more. A separately isolated streaming-validation repair is in development and is not part of the tested source named here.

Later on 7 September, the separate [streaming archive milestone](2026-09-07_ARCHIVE_STREAMING.md) removed that implicit resume ceiling and verified a real 100,001-chunk continuation. It does not change the historical scope of this reliability record.
