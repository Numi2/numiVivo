# Property-directed electronic-space refinement

This example compares matched fixed-geometry electronic profiles, expands whole orbital groups according to measured profile sensitivity, and confirms the selected space on held-out points and the full declared candidate union. It is a three-orbital algebraic fixture, not a molecular reaction, transition-state characterization, or chemical-accuracy benchmark.

Run using a built NumiVivo executable on the supported Apple toolchain:

```sh
numivivo chemistry-template algebraic-space-refinement --output refinement.json
numivivo chemistry-refine refinement.json --output refinement.result.json
numivivo chemistry-export-space refinement.result.json --point barrier-point --output anchor.json
numivivo chemistry-template solver-selected-ci --output selected-solver.json
numivivo chemistry-solve anchor.json --solver selected-solver.json --output anchor.result.json
```

The fixture has an orbital with a larger almost-common absolute energy correction and another that changes the relative profile. The refinement should choose active orbitals `[0,2]`, then compare against `[0,1,2]` at all five points. This expectation is an algebraic regression assertion, not an assumed real-chemistry result.

Use `algebraic-selected-space-refinement` to execute the same refinement with selected CI rather than direct full-sector CI. The selected solver reuses the native Davidson/action implementation, diagnoses the complete connected omitted residual, and never calls its PT2 correction a certified error bound.

`chemistry-refine` returns exit code 2 when it retains a report with non-established sensitivity, including failed holdout or bounded in-run exhaustion. Its diagnostic result and receipt are still retained. Invalid input, preflight rejection and other hard failures return exit code 1. Do not use the mere existence of a receipt as a success gate. The default accepted selected-CI solver refuses unconverged results; `solver-selected-ci-probe` explicitly permits a diagnostic result, with exit code 2 when unconverged.

## Orbital information

```sh
numivivo chemistry-template orbital-information-selection --output pairs.json
numivivo chemistry-correlations state.json --selection pairs.json --output information.json
```

`state.json` must contain an explicit `VivoCIState` (the `state` object in a selected-CI result), not the whole solver result. Library callers may pass that state directly. Without `--selection`, the command requests all single-orbital marginals but no pairs. Missing pairs remain unexamined, not zero. Entropies use natural logarithms; mutual information is `S(i)+S(j)-S(i,j)` with no factor of one half.

## Verification commands

Portable numerical subset, including negative tests and recorded source/toolchain identities:

```sh
bash Tools/check_property_refinement.sh /tmp/numivivo-property-refinement
```

Actual Apple CLI integration (not the complete biological/Metal product):

```sh
bash Tools/NativeChemistry/build_scoped_cli.sh /tmp/numivivo-chemistry-cli
python3 Tools/NativeChemistry/check_refinement_cli.py /tmp/numivivo-chemistry-cli/numivivo-chemistry /tmp/numivivo-refinement-cli
```

Python is a test driver only. The CLI and numerical implementation do not invoke Python. The full-package test suite additionally contains `PropertyRefinementWorkflowTests` for schema dispatch, source/state binding, cache replay and export provenance.

The [design contract](../../Documentation/Design/PROPERTY_DIRECTED_REFINEMENT.md) specifies current guarantees and remaining integration boundaries.

## One typed workflow recipe

The same numerical and validation stages are registered in the general native
workflow system. The recipe runs refinement, anchor materialization, selected
CI, Hamiltonian-validated state extraction and orbital-information analysis.
It uses the same algebraic fixture, not an additional chemistry implementation.

```bash
numivivo workflow-template property-refinement --output refinement-recipe.json
numivivo workflow-plan refinement-recipe.json --output refinement-plan.json
numivivo workflow-run refinement-recipe.json --store .numivivo/refinement-artifacts --output refinement-report.json
```

A bounded exploratory report may be retained by the refinement node. An
unconfirmed report fails at anchor export and blocks the downstream solver.
Only a residual-converged state passes the state-extraction node. Repeating a
recipe/store validates existing artifacts before reuse; receipt existence does
not weaken any of these requirements. Existing workflow exports require an
explicit `--force` to replace previously written files.

## Prepared molecular and multistate requests

```sh
numivivo chemistry-template molecular-h2-space-preparation --output molecule.json
numivivo chemistry-prepare-space molecule.json --output prepared-refinement.json
numivivo chemistry-refine prepared-refinement.json --output molecular-refinement.json
numivivo chemistry-export-space molecular-refinement.json --point h2-2 --output molecular-anchor.json
numivivo chemistry-template algebraic-multistate-refinement --output states.json
numivivo chemistry-refine states.json --output states.result.json
numivivo workflow-template molecular-property-refinement --output molecular-recipe.json
numivivo workflow-run molecular-recipe.json --store .numivivo/molecular-refinement --output molecular-run.json
```

The H2 example executes physical Gaussian integrals and adjacent cross-geometry overlaps. Its declared interior point is not a transition state and the small basis is not a barrier-accuracy recommendation. The two-state algebraic example requires expansion even though the two state-specific barrier shifts cancel in their average. V2 state-averaged-CASSCF requests use the existing native optimizer and retain the optimized orbital frame in each exported Hamiltonian.
