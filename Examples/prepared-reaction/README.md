# Prepared molecule to conditional reaction probability

This example connects the existing engines in one reproducible route:

1. Prepare an explicitly supplied hydrogen molecule and compile its supplied force-field library.
2. Run two seeded Metal sampling replicas and retain their accepted coordinates.
3. Export one exact accepted state and map its coordinates into both hydrogen-molecule endpoint guesses for a fresh H + H₂ exchange calculation.
4. Qualify the minima and saddle, verify both reaction directions at four path refinements, and calculate the molecular transition-state-theory rate.
5. Calculate the probability that a tagged hydrogen atom has reacted by each requested time in an externally maintained H₂ bath.

On an Apple host with Swift 6 and Metal, from the repository root:

```sh
swift build -c release --jobs 3
python3 Tools/Platform/run_prepared_reaction_example.py \
  --binary .build/release/numivivo --out /tmp/numivivo-prepared-reaction-example
```

Use a new empty output directory. The runner uses public CLI commands, retains every request, result, immutable workflow receipt and command log, and writes `checks.json`. The final time series is in `encounter-result.json`. Full reaction validation can take several minutes because the source chemistry is reconstructed at downstream verification boundaries.

The inputs are [preparation.json](preparation.json) and [campaign.json](campaign.json). The harmonic sampling force field is explicitly synthetic: masses 1 Da, bond length 0.125 nm, force constant 1000 kJ mol⁻¹ nm⁻², and zero charge/Lennard-Jones interactions. It runs at 300 K. Its short sampling prefix deliberately does not satisfy the declared convergence criteria; it supplies initial geometry only.

The destination is the existing full-CI/STO-3G H₃ exchange model with masses 1.008 Da and thermochemistry at 298.15 K. The isolated atom, saddle guess, atom identities, electron counts, connectivity settings, standard pressure and transmission coefficient remain explicit supplied inputs. Geometry transfer does not establish equality between the classical and electronic models.

The rate has units Pa⁻¹ s⁻¹. The example declares a constant H₂ partial pressure of 101325 Pa. Multiplying the rate by that pressure gives a first-event hazard in s⁻¹. The existing kinetic solver calculates survival, reaction probability, flux and local sensitivities. This is a conditional model prediction with assumed transmission; it does not establish an experimental rate, finite-bath depletion, equilibrium populations, biological efficacy or uncertainty coverage.

The general interfaces are `reaction-seed`, `reaction-seed-verify`, `reaction-encounter-request`, and the `vivo.native.reaction-qualification` / `vivo.native.conditional-encounter` workflow operations. Run `numivivo reaction-seed-help` for discovery. The Python runner shows the exact request and recipe formats without adding a second numerical implementation.

An additional retained native integration fixture is opt-in:

```sh
swift build -c release --build-tests --jobs 3 -Xswiftc -enable-testing
NUMIVIVO_FULL_ROUTE_TESTS=1 \
NUMIVIVO_TEST_SOURCE_COMMIT="$(git rev-parse HEAD)" \
NUMIVIVO_TEST_ARTIFACTS=/tmp/numivivo-prepared-reaction-native \
swift test -c release --skip-build --no-parallel --filter PreparedReactionFullRouteTests
```

The ordinary host suite includes missing-parameter rejection. The full native fixture additionally checks changed temperature, component identity and numerical rate rejection. Its successful receipt is written only after all required checks pass.
