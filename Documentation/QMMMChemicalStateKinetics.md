# QM/MM chemical-state thermodynamics and pathway kinetics

NumiVivo keeps conditional QM/MM reaction rates separate from the thermodynamic populations and interconversion kinetics of bound chemical states. A protonation state, tautomer, conformer or other chemically distinct bound state is therefore not hidden inside one averaged activation barrier.

Each elementary pathway entering the kinetic layer must already be represented by a converged `VivoQMMMReplicatedFreeEnergyRateResult`. Parallel pathways within one state add as first-order rates. Population generation and state-exchange kinetics remain separate evidence objects.

## pH-dependent state populations

`VivoQMMMChemicalStateThermodynamics` reweights an **explicitly enumerated** state set in the semigrand canonical proton ensemble. It does not enumerate protonation states or predict pKa values.

For each state the caller supplies an evidence-backed relative semigrand Gibbs free energy at a declared reference pH and a relative bound-proton count. Moving from reference pH to target pH changes the state's semigrand free energy by

```text
Delta G_i(pH) = Delta G_i(reference)
              + n_i R T ln(10) [pH - referencePH]
```

and populations are normalized Boltzmann weights of those shifted free energies. Only differences in free energy and proton count matter: a common free-energy offset or a common proton-count offset cancels exactly.

A one-unit pH increase therefore changes the relative weight of states differing by one bound proton by exactly a factor of ten, before any state-specific free-energy difference is applied. States with identical proton stoichiometry retain the same conformational population ratio across pH under this reweighting contract.

```sh
numivivo qmmm-chemical-state-populations thermodynamics.json \
  --output state-populations.json
```

Workflow operation: `vivo.platform.qmmm-chemical-state-thermodynamics`.

The result includes an immutable population-evidence fingerprint that can be attached to the state populations used by the downstream kinetic network. This operation does not calculate the input state free energies, prove that the supplied state set is complete, or perform constant-pH dynamics.

## Rapid pre-equilibrium

Use `VivoQMMMChemicalStateNetwork` only when state exchange is explicitly justified as fast relative to chemical conversion. For state population `p_s` and conditional pathway rate `k_sj`, the effective rate is

```text
k_eff = sum_s p_s sum_j k_sj
```

The implementation requires the populations to sum to one, verifies state/environment identity, requires origin/evidence for every population, and refuses requests that do not assert rapid pre-equilibrium. Populations may come from the semigrand calculation above or from other measured/fitted/calculated evidence.

```sh
numivivo qmmm-chemical-state-network state-network.json \
  --output effective-rate.json
```

Workflow operation: `vivo.platform.qmmm-chemical-state-network`.

## Slow or intermediate exchange

When interconversion competes with reaction, a population-weighted scalar rate is not generally valid. Use `VivoQMMMChemicalExchangeNetwork` instead.

The request contains an explicit finite set of transient chemical states, evidence-bound initial populations, independently replicated state-specific QM/MM pathways, directed interconversion rates in `1/s`, and observation times. Chemical reaction is represented as absorbing loss from each transient state.

For the row-vector transient probability `p(t)`, the model is

```text
dp/dt = p Q
Q_ij = k_i->j                         i != j
Q_ii = -(sum_j k_i->j + k_chemical,i)
```

where `k_chemical,i` is the sum of the independently qualified parallel QM/MM pathway rates in state `i`.

The runtime evaluates the finite-state CTMC with bounded uniformization. It returns, at every requested time:

- unreacted probability in every state;
- total survival probability;
- reacted probability;
- instantaneous chemical hazard;
- the time-dependent apparent first-order rate `-ln S(t)/t` when defined.

It intentionally does not reduce a multi-exponential survival curve to one global rate constant. If every state's chemical loss rate is identical, arbitrary state exchange correctly leaves the total survival equal to `exp(-k t)`; this invariant is covered by regression tests.

```sh
numivivo qmmm-chemical-exchange-network exchange-network.json \
  --output transient-kinetics.json
```

Workflow operation: `vivo.platform.qmmm-chemical-exchange-network`.

## Transient external validation

A slow/intermediate exchange model must not be validated by comparing one fitted first-order constant when its predicted survival is not single exponential. `VivoQMMMChemicalExchangeValidation` therefore compares the immutable exchange-network result directly with measured survival fractions.

Each target carries an exact observation time, survival fraction, optional probability standard deviation, complete environmental context, and immutable measured-source fingerprint. A target is directly comparable only when its compound, target, variant, site, host, temperature, pH and ionic strength match the simulated ensemble and the target time is one of the exchange run's explicit observation times.

Acceptance can gate both absolute probability error and, when measurement uncertainty is supplied, the standardized residual. Off-grid times and condition-mismatched measurements remain inspectable but are not silently interpolated or counted as validation evidence.

```sh
numivivo qmmm-chemical-exchange-validate validation.json \
  --output transient-validation.json
```

Workflow operation: `vivo.platform.qmmm-chemical-exchange-validation`.

## Evidence and scope

State free energies, equilibrium/initial populations and exchange-rate parameters are all evidence-bound. Measured, fitted and calculated inputs require immutable fingerprints while assumptions remain explicit. State-specific reaction pathways retain their complete replicated PMF/transmission provenance. Measured transient validation targets likewise require immutable evidence snapshots.

This layer still does not enumerate missing protonation states/tautomers/conformers, calculate the supplied microstate free energies, generate exchange barriers/rates, add quantum tunnelling, or discover alternative reaction mechanisms. Those quantities must come from qualified calculations or external evidence. Joint uncertainty propagation across state free energies, populations, exchange rates and reaction pathways also remains a separate posterior/sensitivity problem rather than being silently collapsed.

The kinetic-regime distinction is strict:

```text
fast state exchange       -> rapid-pre-equilibrium population weighting
comparable/slow exchange  -> explicit transient exchange network
unknown exchange regime   -> not qualified for either simplification without evidence
```
