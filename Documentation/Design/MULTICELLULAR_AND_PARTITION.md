# Multicellular and physiological distribution

## Scope

NumiVivo contains two distinct higher-scale model classes:

- a spatial phenotype-population runtime for cell-state density fields;
- a physiological compartment runtime for analyte distribution.

They are intentionally separate. A cell population is not a physiological concentration, and a tissue partition coefficient is not a phenotype transition rate.

## Spatial phenotype populations

The population state is a density field:

```text
population[phenotype, voxel]
```

Every phenotype declares:

- baseline birth rate;
- baseline death rate;
- carrying capacity per voxel;
- diffusion coefficient;
- minimum and maximum density;
- whether the field is externally maintained.

Externally maintained fields are read from a coupling authority and cannot be mutated by local population dynamics.

### Local demographic dynamics

The local demographic term contains bounded logistic growth and death:

```text
birth = b_i n_i (1 - total_population / carrying_capacity_i)
death = d_i n_i
```

The implementation rejects non-finite coefficients and invalid carrying capacities. Candidate values outside declared bounds result in adaptive reduction or rejection rather than silent normalization.

### Phenotype transitions

A transition moves density from one phenotype to another. It contains:

- source and destination phenotype;
- constitutive, activated, or repressed mode;
- optional regulator field;
- baseline rate;
- maximum regulated rate;
- threshold;
- Hill coefficient.

Source and destination changes are equal and opposite at a voxel unless the model explicitly includes an independent birth or death term. Transition rates are bounded before integration.

### Inter-population effects

Interactions are directed coefficients from an actor phenotype to a target phenotype. Positive coefficients promote target growth; negative coefficients suppress it. The coefficient is declared per actor-density per second.

Dense interaction tables are compiled before dispatch. Models that exceed the bounded interaction-table capacity are rejected rather than partially loaded.

### Spatial movement

Phenotype movement supports:

- diffusion on a regular three-dimensional grid;
- upwind advection from a declared velocity field;
- no-flux boundaries;
- periodic boundaries;
- absorbing boundaries.

Grid dimensions, spacing, voxel count, and boundary mode are immutable runtime contract fields. The scheduler computes explicit diffusion and advection limits before dispatch.

### Transaction sequence

A population step follows:

1. Read committed phenotype fields and regulator fields.
2. Compute local birth, death, transitions, and interactions.
3. Apply spatial diffusion and advection.
4. Validate finite values and phenotype bounds.
5. Evaluate coupled monitors.
6. Commit all phenotype fields together or retain the previous state.

## Physiological compartment distribution

Physiological state is represented as:

```text
concentration[analyte, compartment]
```

Every compartment has an explicit volume. Every analyte has a unit and concentration bounds.

### Reversible partition edge

A partition edge links one analyte across two compartments. It declares:

- source and target compartment;
- target/source equilibrium partition coefficient;
- exchange clearance in volume per time;
- source and target unbound fractions;
- evidence class and optional evidence reference.

The signed amount flux is:

```text
J = clearance × (
      source_unbound × source_concentration
      - target_unbound × target_concentration / partition_coefficient
    )
```

Concentration derivatives are:

```text
dC_source/dt = -J / source_volume
dC_target/dt =  J / target_volume
```

The same amount flux is applied with opposite sign. Therefore the edge conserves represented analyte amount up to integration and floating-point error.

### Model restrictions

The compiler rejects:

- zero or negative volumes;
- non-positive partition coefficients;
- negative clearance;
- unbound fractions outside zero through one;
- duplicate physical edges for the same analyte-compartment pair;
- unknown analyte or compartment references;
- invalid evidence classes;
- non-finite values.

A pair is represented once, because adding an opposite duplicate edge would count the same physical exchange twice.

### Integration and stability

The Metal runtime evaluates the partition system with a two-stage RK2 method. A conservative rate bound is computed from exchange clearance, unbound fractions, partition coefficient, and inverse compartment volumes.

The requested time step is intersected with this bound. Negative candidates beyond numerical tolerance, non-finite candidates, or upper-bound violations trigger reduction or rejection.

### Amount certification

A caller may request total analyte amount before and after a committed step:

```text
total_amount[a] = sum_c concentration[a,c] × volume[c]
```

The totals are included in the step certificate. This is an execution check for the declared compartment model; it is not evidence that unrepresented elimination, synthesis, binding, or transport processes are absent in biology.

## Directed physiological flow

The existing physiology runtime also supports directed compartment flow, analyte-specific masks, first-order clearance, flow-over-volume clearance, saturable clearance, dosing events, sparse incidence tables, and candidate publications.

Partition edges and directed flow represent different mechanisms:

- partition edges approach an equilibrium concentration ratio;
- directed flow transports material according to bulk compartment flow;
- clearance removes represented material;
- secretion or uptake crosses the NumiVivo coupling boundary.

A combined model must preserve these distinctions in its artifact and certificate.

## Molecular–physiology coupling

The coordinator maps physiological exposure to molecular inputs and molecular secretion or uptake back to physiological candidates.

Each mapping declares:

- source and destination indices;
- unit transform;
- update mode such as replacement, addition, rate, minimum, or maximum;
- bounds;
- relaxation when fixed-point iteration is enabled.

The coordinator can iterate candidates until its residual is below tolerance. Neither subsystem commits before the joint residual, numerical checks, and safety monitors pass.

## NumiTissue coupling

The population runtime exchanges:

- phenotype distributions;
- proliferation and death drives;
- differentiation requests;
- migration drives;
- extracellular source and sink terms;
- tissue injury and repair state.

NumiTissue remains authoritative for cell identity, lineage, development, and tissue organization. NumiVivo remains authoritative for the molecular-control program and declared intracellular state. Ownership is explicit per channel.

## NumanX coupling

NumanX provides mechanical and transport context, including:

- tissue deformation;
- pressure and perfusion;
- oxygen and nutrient fields;
- temperature;
- velocity fields;
- strain and mechanotransduction inputs;
- damage and injury fields.

NumiVivo may publish molecular source terms, permeability changes, matrix-remodelling drives, and other declared coupling outputs. It does not directly mutate NumanX state.

## NumiBrain coupling

NumiBrain may provide endocrine, autonomic, neural, pain, arousal, and homeostatic signals. NumiVivo may publish inflammatory, metabolic, injury, and molecular-state summaries.

Neural and molecular runtimes advance at different rates. The coupling scheduler defines hold, interpolation, accumulation, or event semantics for every channel.

## Fidelity promotion

Promotion from well-mixed molecular execution to multicellular or physiological execution requires:

- explicit geometry or compartments;
- mapping of molecular species to fields or compartments;
- transport and partition parameters;
- cell-population ownership;
- coupling channels;
- new step bounds;
- updated uncertainty and evidence records.

Promotion is not a change of one runtime flag. It creates a new fingerprinted artifact lineage.

## Current limitations

The population and partition runtimes are source-complete foundational implementations but have not received package-wide compiler, shader, numerical, or performance verification. They currently model declared density and compartment systems; they do not automatically infer tissue anatomy, vascular structure, immune populations, or patient-specific parameters.
