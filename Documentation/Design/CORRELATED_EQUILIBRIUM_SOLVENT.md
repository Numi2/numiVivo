# Correlated-density equilibrium solvent

`VivoCorrelatedSolvation` solves a normalized FCI or fixed-orbital CASCI state in
an iteratively updated smooth C-PCM field. The reaction workflow exposes this as
`correlatedSolvent`. It is different from `solvatedPath`, which intentionally
retains the existing RHF-reference-frozen field convention. Neither method is
silently substituted for the other.

For an input correlated density D, the smooth cavity produces V(D). The next CI
calculation solves Hgas + V(D), including its occupied-core contribution when a
CAS partition is declared. Its output density D' produces the new solvent field.
The reported energy is <Hgas> + Gpol(D'), not the effective CI eigenvalue plus
polarization: the electron-field contraction is subtracted before adding Gpol.
Convergence requires density, reaction-potential and energy-change tolerances,
in addition to the true CI residual and the solvent equation residual. An
iteration limit is a failure, not an automatically accepted result.

Validation reconstructs the full fixed-spin CI residual, the AO density, the
surface field, the energy terms and an additional density-feedback step. A
self-reported convergence flag or a small change in energy alone is insufficient.
The initial density used by the final iteration is retained for that check.

The density is a normalized CI probability state. Coupled-cluster response RDMs
and sums of overlapping impurity densities are not admissible substitutes. This
implementation does not claim correlated-density ECC-DMET/PCM self-consistency
for an overlapping partition. That requires a declared physically representable
global density and its consistent energy functional.

## Native commands

From a complete macOS checkout:

```bash
swift build -c release
.build/release/numivivo reaction-template h2-equilibrium-cpcm --output solvent.json
.build/release/numivivo reaction-run solvent.json --output solvent.result.json
.build/release/numivivo reaction-template h2-equilibrium-minimum --output minimum.json
.build/release/numivivo reaction-run minimum.json --output minimum.result.json
```

The existing artifact store checks output acceptance on new execution and reuse.
No Python dependency is introduced into the production commands. Coordinates
are Bohr and energy is Hartree. Explicit solvent configuration and basis data
are part of the request identity.

## Nuclear differentiation

`equilibriumFullCI` is an explicit nuclear solver selection. At every displaced
geometry it equilibrates the full CI density and solvent, then differentiates
the resulting total equilibrium energy with the existing step-halved numerical
gradient and Hessian machinery. The optional `correlatedSolventConfiguration`
controls the inner solve. Supplying it to a different method is rejected.
Legacy `fullCI` and `anchoredECC` solvent calculations retain their original
reference-frozen meaning when decoding older requests.

Truncated CAS nuclear differentiation is not enabled by assigning orbital
indices in an independently re-orthogonalized basis. Such a calculation needs a
specified transported frame or orbital-response method. The nuclear equilibrium
route therefore requires the full orbital space; fixed-geometry CASCI remains
available through `correlatedSolvent`.

This is bounded finite-basis electronic structure, not large-molecule scalable
FCI. Harmonic thermochemistry on this energy surface remains an ideal RRHO model;
non-electrostatic solvation, conformational entropy and a validated kinetic rate
are not supplied by electrostatic field convergence.

## Reference methodology

Independent conformance uses pinned PySCF 2.8.0 integral, CI and PCM operations
with explicitly matched orbital frames, basis coefficients, cavity radii and
quadrature. H2 alone is insufficient to demonstrate nontrivial polarization;
the reference set also includes polar LiH full CI and occupied-core CAS(2,2).
The production implementation does not import this reference code.

Primary reference interfaces:
- https://pyscf.org/user/solvent.html
- https://pyscf.org/_modules/pyscf/solvent/pcm.html

No acrylamide-methanethiolate/BTK barrier or reaction connectivity is implied by
these finite-system conformance problems.
