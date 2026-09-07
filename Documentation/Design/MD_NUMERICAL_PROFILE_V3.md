# Classical MD numerical profile v3

`VivoMDExecutionIdentity.current` identifies the executable numerical algorithm separately from physical system and configuration fingerprints. Profile v3 changes periodic coordinate arithmetic and classical charge PME interpolation. Numerical qualification belongs to an exact source and device; the profile name alone is not evidence of accuracy or reproducibility.

## Periodic geometry

The supported native nearest-image profile remains orthogonal periodic cells. A shared shader helper uses reciprocal coordinates only to select integer lattice images, then subtracts those lattice vectors from the original physical vector with fused arithmetic. A zero selected image leaves the physical vector unchanged. Rounding near a cell face can select an adjacent equivalent image; the helper does not promise an exact fractional-coordinate classifier at that boundary.

Force, constraint, constraint-impulse, minimization, virtual-site, neighbor and PME correction paths share this operation. Reconstructing every coordinate through cell and reciprocal matrices is not an FP32 identity: repeated zero-correction projection passes can otherwise move stationary atoms and create false constraint impulses.

Molecular-center volume proposals pack the actual proposed cell and reciprocal cell used by subsequent validation. Center displacement uses the rounded Double value of `scale - 1`, with a fused shift. Internal molecular positions are wrapped by subtracting lattice offsets in that proposed cell. Constraint tolerance and iteration count retain their explicit configuration values.

## Classical charge PME

Classical charges use sixth-order cardinal B-spline assignment. Gather differentiates the same assignment used by spread, so forces derive from that mesh energy. The influence deconvolves the matching discrete nodal modulus, whose one-dimensional centered form is `(66 + 52*cos(theta) + 2*cos(2*theta))/120`. Self, background and primary-pair exception conventions remain explicit physical terms.

The common sixth-order basis is also used by the existing multipole path; its higher derivatives and multipole influence arithmetic retain their previous implementation. Sharing a basis does not imply that charge and multipole paths have identical numerical evidence.

`VivoPMEPlan.make` declares order six and rejects an explicit order-four request. The existing minimum four-point grid dimension remains accepted: six support nodes wrap and sum periodic aliases on that axis. Explicit grids are not silently enlarged. Coarse-grid acceptance is a representational contract, not a force-accuracy guarantee.

The requested PME tolerance controls the Ewald splitting plan. Mesh resolution, assignment, aliasing, reciprocal truncation and floating-point accumulation have separate errors. The six-by-six-by-six stencil touches 216 contributions instead of the former 64; this is a source-level cost tradeoff, not a measured performance comparison.

## Restart and derived state

An MD checkpoint with a missing or older numerical contract can still be decoded as historical data. Direct restore, stage transfer, protocol plans and protocol cursors require the current contract before continuation. Do not change a historical checkpoint's contract field to make it pass. Explicitly importing physical state starts a new numerical run with retained provenance; it is not continuation of the old trajectory.

Metal-backed constant-pH and nuclear potential identities must also bind this contract, because their generic stored state can be converted into an MD checkpoint. The same requirement reaches NCMC, replica exchange, ring-polymer sampling and reactive labels through their potential or Hamiltonian identity. Generic CPU potential schemas and explicit analysis of historical static traces remain separate.

QM/MM free-energy archive execution must bind the contract in both request and cursor, including completed-window boundaries where no active MD checkpoint remains. This prevents mixing old numerical traces with newly executed windows under an unchanged physical configuration.

## Qualification

The [frozen MD qualification plan](../Audit/2026-09-07_MD_QUALIFICATION_PLAN.md) defines independent charged-PME, finite-time Langevin, constrained thermalization and ideal NPT transition gates. The initial v2 run passed the thermal gates and exposed the PME and NPT failures motivating v3. Its raw observations remain retained. Source implementation, compilation and successful numerical acceptance are separate stages; consult the subsequent native evidence record before claiming a v3 result.
