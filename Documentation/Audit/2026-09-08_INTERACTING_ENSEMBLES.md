# Longer conservation and interacting ensemble qualification

This milestone extends the [v6 reference panel](2026-09-08_FRONTIER_MD_REFERENCE_PANEL.md)
with a 50-times-longer conservation campaign and an independent interacting-NVT
qualification harness. It uses the same owning Metal runtime and smooth-LJ input
models. Python performs reference checks and statistics; it supplies no production
forces.

## Execution controls and retained failures

The original 128-iteration campaign passed the short gate. A separately named
32-iteration variant, using the runtime default, failed preparation or thermalization
on alanine, DNA, membrane, salt water and orthogonal water with candidate status 8.
DHFR and the vacuum complex completed their short runs. All 21 static comparisons
passed; static force agreement did not rescue the five failed dynamics cases.

The subsequent 64-iteration variant passed all seven prepared systems' static
comparisons, 100-step execution, independent checkpoint constraints and original
0.1 ps three-timestep conservation/refinement gate. The constraint tolerance remains
1e-6. Increasing available solver iterations changes work, not the accepted residual.
The 32-iteration failures remain failed and were not retried under that name.

The original water source still fails independent preparation because its cell
cannot support the declared cutoff. Orthogonal water is a separate fixture. The
retained full-panel result therefore stays nonzero even when prepared cases pass.

## Longer conservation

The prospective policy extends 0.1 to 5 ps at 1, 0.5 and 0.25 fs on all seven
prepared systems. It retains the maximum energy deviation of 0.01 kJ/mol per
particle, minimum RMS refinement ratio 2 and numerical floor 1e-5 kJ/mol per
particle. All variants must share the exact prepared position and velocity words.

All seven prepared systems passed all three time steps: 21 completed trajectories.
All variants retained identical prepared position/velocity words. The largest
energy deviation stayed below the unchanged 0.01 limit, and every RMS refinement
ratio exceeded the required factor 2.

| System | Largest energy deviation, kJ/mol per particle | RMS reduction, 1 to 0.5 fs | RMS reduction, 0.5 to 0.25 fs |
|---|---:|---:|---:|
| Alanine | 0.000680452 | 3.680 | 6.537 |
| DHFR | 0.000198314 | 3.785 | 2.383 |
| Vacuum complex | 0.008051145 | 4.364 | 4.231 |
| DNA | 0.002505599 | 4.028 | 4.046 |
| Membrane | 0.001008762 | 4.761 | 4.296 |
| Salt water | 0.000633983 | 3.786 | 7.217 |
| Orthogonal water | 0.000319156 | 4.050 | 2.266 |

![Measured longer conservation refinement](frontier-long-nve-refinement.svg)

The [compact observations](frontier-long-nve-observations.json) bind the raw
qualification, source identities, report hashes and independent endpoint audits.
Five ps remains far shorter than a production molecular research trajectory.

## Interacting NVT

The [ensemble policy](../../Tools/Benchmarks/ensemble_policy.json) fixes orthogonal
rigid water, 300/305 K, three seeds per temperature, 1/0.5 fs, friction 10/ps and
50 ps per run. It excludes the first 10 ps and retains 800 production observations
per replica. The higher-temperature seeds add 1,000,000 to the base seeds, providing
separate temperature streams while matching initial states across time steps.
These seed and 64-iteration execution controls were finalized before NVT data.
Earlier 32-iteration policy snapshots remain visible and were not measured as NVT.

Independent analysis counts six kinetic degrees of freedom for each nondegenerate
rigid triangle. It refuses other constraint graphs. It checks kinetic mean, width
and prescribed gamma-distribution shape, plus the expected two-temperature
potential-energy density-ratio slope. Fixed one-ps block bootstraps and quality
diagnostics retain correlation, drift, replica disagreement and insufficient
precision as inconclusive. Failed observations are retained even when precision
is insufficient to classify the overall result.

Measurement status: admitted and running. Implementation tests and the completed
NVE gate do not establish a measured NVT pass.

## Software and evidence identity

- Native source: `c7ded9f9784d238a329e2ca0ae34b01676165ee2`, numerical contract v6.
- Frozen executable SHA256: `eed5b517985dbda7ed545bd82407d99113ee65b9dc64f972dc5a0eefcfcae60a`.
- Conservation orchestration: `ce39bd62daf9fae2581532c7186acb5e552c731a`.
- Ensemble implementation/policy: `13e041ff39d43d49503f9074818f4551e1c3b2ba`.
- All 21 Python regression checks passed in the isolated NumPy 2.5.3 / SciPy 1.18.1
  environment. Controls reject wrong temperatures, suppressed fluctuations,
  matching-moment wrong shapes, unchanged two-temperature distributions, missing
  observations, altered inputs and incomplete matrices.
- A separate exact-marginal correlated-gamma control passed. Its block uncertainty
  was 1.804 times the inappropriate independent-sample estimate. This is synthetic
  method validation, not a molecular trajectory.
- The first statistical implementation failed one of ten controls because a
  bootstrap logistic fit stalled at floating-point objective resolution. Its log
  remains failed. A curvature-based numerical termination fix passed all prescribed
  fits; no bootstrap samples or acceptance thresholds were replaced.
- Native Sources, Tests and Package.swift match the previously built and fully
  tested v6 runtime. No new native build or full native regression is claimed.

## Independent endpoint energy

The additional offline endpoint auditor uses the original serialized OpenMM
System and the saved final coordinates. It independently sums start/end kinetic
energy from masses and exact velocity words. Potential energy retains the original
0.002 kJ/mol-per-particle limit; kinetic agreement uses the v6 FP32 pairwise-reduction
forward-roundoff bound. Five endpoint controls passed, including invalid images
and plausible but incorrect reported energy.

The initial auditor incorrectly sent atom-wrapped coordinates directly to OpenMM.
All periodic potential comparisons failed, while the vacuum case and all kinetic
comparisons passed. These failures remain in `endpoint-nve-0.json`. OpenMM requires
whole molecules for its default bonded forces and exceptions. The corrected
auditor applies only integer cell translations, checks graph-cycle consistency
and every exception image, and records coordinate/image hashes. It refuses
ambiguous or winding molecular graphs. No parameter or acceptance limit changed.

With the corrected conversion, all 21 NVE trajectory endpoints passed. The largest
potential error was 2.02e-5 kJ/mol per particle, versus the unchanged 0.002 limit.
All 42 start/end kinetic checks passed their independently computed roundoff bound.
The original checkpoint files and native trajectories were reused unchanged.
An analytic boundary-bond control independently recovered 0.03125 kJ/mol; supplying
the split molecule without conversion gave 30.03125 kJ/mol. NVT endpoint audits
remain to be completed.

Remote raw evidence: `/Users/n/numivivo-ensemble-evidence-20260908`.
Local mirror: `/Users/home/numivivo-ensemble-evidence-20260908`.
The previous sealed native/reference evidence remains in the corresponding
`numivivo-frontier-evidence-20260908` roots. New evidence retains commands, policies,
source snapshots, binary/shader identities, requests, XML, reports and failed logs.

## Scientific boundary

Conservation and global kinetic/configurational consistency are necessary tests.
They do not establish complete equilibration, all kinetic-energy partitions,
experimental water properties, NPT, transport, free energies, triclinic execution,
reaction accuracy or matched-accuracy speed leadership. On-step native BAOAB
kinetic energies also cannot be equated directly with OpenMM LFMiddle's half-step
velocity estimator. The remaining frontier roadmap gates stay open.

References: [Merz and Shirts, physical validation](https://doi.org/10.1371/journal.pone.0202764),
[OpenMM velocity convention](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LangevinMiddleIntegrator.html).
[OpenMM periodic-coordinate requirements](https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#periodic).
