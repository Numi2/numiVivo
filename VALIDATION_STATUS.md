# Validation status

## 7 September 2026 reliability verification

[The current reliability audit](Documentation/Audit/2026-09-07_COMPLETION_RELIABILITY.md)
records a clean complete build, **224/224** baseline tests across 52 suites,
**22/22** final-source cancellation/MD-recovery/hybrid regressions, and **57/57**
production workflow CLI checks on the physical M4 Pro Mac mini. The audit binds
each result to its exact source revision and retains the remaining scientific
and product-completion boundaries.

The results below remain historical evidence for their named revisions; they
are not automatically requalified by the later reliability work.

## Earlier verified revision

The verified product implementation is the `24e8030` tree plus the focused
acceptance-boundary regression in commit
`b9870f705a80296adf13adeff19b782e3f2277c8`. Documentation and audit files in
this commit do not change product numerical behavior.

The validation was run on the physical Apple M4 Pro Mac mini, not only on
Apple's hosted virtual device. The checkout was clean after validation and
generated artifact stores were kept outside the repository.

## PASS

- `swift build -c release --jobs 3`
- `swift build -c release --build-tests --jobs 3 -Xswiftc -enable-testing`
- `swift test -c release --jobs 3 -Xswiftc -enable-testing`: **23/23** tests in
  four suites, including real Metal/MD execution and mapped connectivity.
- `swift test -c debug --jobs 3 --sanitize=address --filter ScientificClosureTests`:
  **6/6** host-side scientific-closure tests.
- `Tools/check_reaction_qualification.sh`: **130/130** reported checks:
  36 native, 36 independent PySCF, 21 production CLI/store and 37 equilibrium
  solvent checks.
- Production `h3-residual-barrier`, `h2-global-embedding` and committed polar
  LiH global-embedding workflows executed and stored results successfully.
- Production Metal `engagement-run` completed 12,800 commits on the Apple M4
  Pro with maximum fraction-mass error `2.0734297468e-5` against the declared
  `1e-4` tolerance. Its synthetic assumed-parameter and uncertainty flags
  remain true.
- Paper preflight workflows returned explicit negative readiness for both
  acrylamide-methanethiolate and BTK; no missing input was fabricated.

The corrected sanitizer command is SwiftPM's supported `--sanitize=address`
form. The earlier linker-only invocation is not used as evidence.

## Residual-enriched H3/6-31G campaign

The legacy orbital-only campaign remains a separate negative assessment:
`reducedAccuracyNotEstablished`. It has not been relabeled.

The residual-enriched campaign uses nine explicitly mapped scan geometries, a
90-determinant sector, all six spatial orbital modes and one nested common Fock
space. Its dimensions are 9, 18, 27, 36, 45, 54 and 63. The unchanged targets
are 0.001 Hartree for barrier error, relative-profile error and successive-level
change, with two consecutive genuinely reduced levels required.

Measured levels:

| Dimension | Barrier error | Profile error | Successive change | Reference accuracy |
| ---: | ---: | ---: | ---: | :--- |
| 9 | 3.5685895224e-2 | 3.5685895224e-2 | — | no |
| 18 | 5.5993135732e-3 | 5.5993135727e-3 | 3.0086581651e-2 | no |
| 27 | 3.3185396865e-4 | 3.3185396861e-4 | 5.2674596045e-3 | yes |
| 36 | 4.5776274848e-5 | 4.5776274848e-5 | 2.8607769381e-4 | yes |
| 45 | 3.7354548068e-6 | 3.7354547939e-6 | 4.2040820054e-5 | yes |
| 54 | 1.7522174556e-7 | 1.7522174556e-7 | 3.5602330710e-6 | yes |
| 63 | 6.9341230535e-9 | 6.9341230535e-9 | 1.6828763405e-7 | yes |

The accepted window is dimensions **54 and 63**, with `acceptedRound: 6` and
`reducedAccuracyEstablished: true`. Dimension 27 meets the reference-error
criteria but its successive change from dimension 18 is above tolerance; the
focused regression test preserves this distinction. This is a reduced
projected eigenproblem, not fewer orbital modes, fewer determinant-vector
storage requirements, a demonstrated speedup or a protein-scale result.

## Coherent global density and solvent feedback

The global embedding forms one normalized CI state in the orthonormal union of
overlapping fragment spaces, removes duplicate directions, includes cross-term
density contributions and drives the smooth C-PCM feedback loop. Electron
count, occupation bounds, energy reconstruction, density/potential/energy
convergence and returned-field stationarity are checked.

The committed polar LiH example converged in eight solvent iterations with a
14-dimensional union in a 225-determinant sector:

- final density residual: `8.6602955681e-9`;
- returned-field projected residual: `1.8736287756e-10` Hartree;
- external residual: `0.0353106565` Hartree.

The nonzero external residual is retained evidence that selected-space
stationarity is not full-sector convergence. This formulation is distinct from
the original democratic ECC-DMET energy formula. ECC imports require explicit
inactive occupied columns; fractional occupations are not rounded into invented
determinants.

## Production and CI evidence

The previous M4 Pro Metal target/MD/mapped-connectivity repairs remain covered
by the 23-test suite. The target simulation remains synthetic and retains its
assumed-parameter and uncertainty flags.

Successful workflows for the published acceptance-guard revision `b9870f7` were:

- Barrier convergence conformance: run `34045438757`;
- Scientific closure conformance, including Apple and ASAN jobs: run `34045438791`;
- Complete native product contracts: run `34045438774`.

The unchanged native-chemistry workflow also passed on the product implementation
at `24e8030` in run `34017723597`.

The local measured audit is
`Documentation/Audit/scientific-closure-observations.json`.

## BLOCKED ON SCIENTIFIC INPUT

- Acrylamide-methanethiolate: missing mapped precomplex, transition-state,
  product, quantum protocol and reference-energy inputs.
- BTK: missing prepared snapshot, topology, environment and associated mapped
  author-model inputs.

The preflight checker establishes declared content integrity and readiness only;
it cannot authenticate authorship or perform a paper reproduction.

## NOT SCIENTIFICALLY QUALIFIED

- No exact acrylamide-methanethiolate or BTK author-model calculation was run.
- No kinetic rate, transmission coefficient or rate export is qualified.
- The H3 scan is not an optimized saddle, IRC, Gibbs barrier or endpoint study.
- The global CI/C-PCM density is not a completed global ECC-DMET/PCM functional.
- No large-system scaling or production speedup has been established.
