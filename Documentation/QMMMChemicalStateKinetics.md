# QM/MM chemical-state and pathway kinetics

NumiVivo keeps conditional QM/MM reaction rates separate from the population and interconversion kinetics of bound chemical states. A protonation state, tautomer, conformer or other chemically distinct bound state is therefore not hidden inside one averaged activation barrier.

Each elementary pathway entering this layer must already be represented by a converged `VivoQMMMReplicatedFreeEnergyRateResult`. Parallel pathways within one state add as first-order rates. What happens across states depends on the exchange timescale.

## Rapid pre-equilibrium

Use `VivoQMMMChemicalStateNetwork` only when state exchange is explicitly justified as fast relative to chemical conversion. For state population `p_s` and conditional pathway rate `k_sj`, the effective rate is

```text
k_eff = sum_s p_s sum_j k_sj
```

The implementation requires the populations to sum to one, verifies state/environment identity, and refuses requests that do not assert rapid pre-equilibrium. It does not infer proton chemical potentials or equilibrium populations.

CLI:

```sh
numivivo qmmm-chemical-state-network state-network.json \
  --output effective-rate.json
```

Workflow operation: `vivo.platform.qmmm-chemical-state-network`.

## Slow or intermediate exchange

When interconversion competes with reaction, a population-weighted scalar rate is not generally valid. Use `VivoQMMMChemicalExchangeNetwork` instead.

The request contains an explicit finite set of transient chemical states, their initial populations, independently replicated state-specific QM/MM pathways, directed interconversion rates in `1/s`, and observation times. Chemical reaction is represented as absorbing loss from each transient state.

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

CLI:

```sh
numivivo qmmm-chemical-exchange-network exchange-network.json \
  --output transient-kinetics.json
```

Workflow operation: `vivo.platform.qmmm-chemical-exchange-network`.

## Evidence and scope

Exchange-rate parameters use the shared `VivoKineticParameter` contract, so measured, fitted and calculated values require immutable evidence fingerprints while assumptions remain explicit. State-specific reaction pathways retain their complete replicated PMF/transmission provenance.

This layer does not calculate protonation free energies, conformational populations, exchange rates, tunnelling corrections or missing pathways. Those quantities must come from qualified calculations or external evidence. It also does not propagate their joint uncertainty automatically; sensitivity or posterior analysis must treat those inputs explicitly.

The distinction is strict:

```text
fast state exchange       -> rapid-pre-equilibrium population weighting
comparable/slow exchange  -> explicit transient exchange network
unknown exchange regime   -> not qualified for either simplification without evidence
```
