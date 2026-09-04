# Executable mixed-authority reaction example

This is an abstract, sequence-free computational example. The rates and counts are synthetic, not fitted biological data. It contains three independent birth/death components and one constant species. Four independent lanes are stored species-major.

`slow_pool` uses exact SSA, `leap_pool` uses tau-leaping, `bulk_pool` uses deterministic RK2, and `constant_pool` must remain unchanged. An authority override on one reaction applies to its whole connected component, including all propensity-only dependencies. Conflicting overrides are rejected.

From the repository root, on the user's Apple-silicon machine:

```sh
swift run numivivo plan-hybrid Examples/hybrid-reaction-runtime/model.json \
  --authorities Examples/hybrid-reaction-runtime/authorities.json \
  --maximum-step 0.01 --output /tmp/numivivo-hybrid-plan.json

swift run numivivo run-hybrid Examples/hybrid-reaction-runtime/model.json \
  --plan /tmp/numivivo-hybrid-plan.json \
  --initial-counts Examples/hybrid-reaction-runtime/initial-counts.json \
  --publications Examples/hybrid-reaction-runtime/publications.json \
  --lanes 4 --steps 100 --dt 0.01 \
  --output /tmp/numivivo-hybrid-result.json \
  --checkpoint /tmp/numivivo-hybrid-checkpoint.json

swift run numivivo run-hybrid Examples/hybrid-reaction-runtime/model.json \
  --plan /tmp/numivivo-hybrid-plan.json \
  --restore /tmp/numivivo-hybrid-checkpoint.json \
  --lanes 4 --steps 100 --dt 0.01 \
  --output /tmp/numivivo-hybrid-resumed.json
```

The command rejects existing output files unless `--force` is provided. A rejected numerical candidate or exact-work-budget exhaustion returns exit 75 and leaves that step uncommitted. A result report includes the final accepted state, executable plan identity, last certificate, and certificate-chain head. Checkpoints retain both UInt32 and FP32 state; inactive representations are canonical zeros.

Source implementation only: these commands, shaders, numerical methods and checkpoint round trips have not been run during this development pass.
