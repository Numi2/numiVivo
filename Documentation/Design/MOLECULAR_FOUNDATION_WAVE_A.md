# Molecular Foundation — Wave A

Status: source implementation on `main`

Wave A establishes the single molecular identity and classical-system boundary used by later NumiVivo MD, QM/MM, embedding and reaction-path code. It replaces ad-hoc package-specific coordinate/topology representations with validated, fingerprinted NumiVivo objects. It does not claim that the package or Apple GPU backend has been built or numerically qualified in this development pass.

## 1. Canonical molecular representation

`VivoMolecularStructure` owns physical chemical atoms, bonds, residues, chains, conformers and an optional triclinic periodic cell. Atom/residue/chain indices are dense and equal to array order. Coordinates are stored in nanometres. File adapters perform explicit boundary conversion.

`VivoStructureValidator` checks index density, bond uniqueness and endpoints, residue/chain bidirectional membership, conformer shape and finite coordinates, occupancy/B-factor constraints and periodic-cell volume. `VivoStructureIndex` derives adjacency and grouping caches after validation rather than serializing duplicate authority.

`VivoMolecularStructureDocument` binds a canonical structure fingerprint to the payload and can also retain the raw source fingerprint and source format. A changed coordinate, atom ordering or topology therefore changes the executable structural identity.

## 2. Trajectories, mapping and geometry

`VivoTrajectory` / `VivoTrajectoryFrame` define ordered frames with nm positions, optional nm/ps velocities, kJ mol^-1 nm^-1 forces, thermodynamic values and per-frame cells. Triclinic minimum-image geometry is implemented in `VivoPeriodicCell`.

`VivoAtomMapper` supports strict semantic and source-serial bijections. Ambiguous mappings fail. `VivoStructureSlicer` rewrites every dependent atom/residue/chain/bond/conformer index and returns explicit old-to-new and new-to-old maps. `VivoAlternateLocationResolver` deterministically resolves PDB/mmCIF alternate sites by occupancy or preferred label.

## 3. Unified atom selections

One selection AST is shared by CLI preprocessing and future MD/QM-region construction. Supported selectors include explicit indices, elements, atom/residue/chain names, residue ranges, hetero/hydrogen state, bonded-neighbor expansion, distance selection with optional periodic minimum image, and Boolean composition.

The deterministic query surface includes forms such as:

```text
element(C,N,O,S)
and(protein,resseq(450,500),not(hydrogen(true)))
within(0.35,resname(LIG),false)
```

Distances in selection queries are nanometres.

## 4. Native structure format boundaries

Implemented source adapters:

- PDB: multi-model identity checking, ATOM/HETATM, CRYST1, CONECT, explicit Å -> nm conversion and PDB output.
- mmCIF: CIF tokenizer, `_atom_site` loops, multiple models and crystallographic cell import.
- SDF/MOL V2000: concrete atoms/bonds, formal charges and `M  ISO`; unsupported query/radical/stereo semantics fail instead of being silently flattened. V3000 currently fails explicitly.
- MOL2: molecule/atom/bond/substructure structural import and output.
- SMILES: graph identity, branches, ring closures, disconnected components, aromatic bonds, isotope/formal charge for the supported subset. Directional double-bond stereo, tetrahedral stereo, bracket hydrogen counts and atom-map labels currently fail explicitly because the authoritative stereo/implicit-H IR is not complete.

`VivoStructureCodec` centralizes format detection, import/export and canonical fingerprints. No adapter automatically substitutes a lower-fidelity output format.

## 5. Topology reconstruction

`VivoBiopolymerTopologyBuilder` adds standard amino-acid backbone/side-chain connectivity, peptide links from chain order, optional disulfide detection and conservative hydrogen attachment. Unknown residues are reported rather than guessed.

`VivoCovalentBondPerception` provides an opt-in geometry route for ligands/unknown residues using conservative covalent radii and simple valence caps. Inter-residue and metal perception are disabled by default. This is structural perception, not a general valence/aromaticity engine.

## 6. Unit-explicit force-field boundary

`VivoClassicalSystem` is the Wave B execution input. Canonical units are:

- distance: nm
- time: ps
- particle mass: dalton
- electric charge: elementary charge
- energy: kJ/mol
- bond force constant: kJ mol^-1 nm^-2 with `1/2 k (r-r0)^2`
- angle force constant: kJ mol^-1 rad^-2 with `1/2 k (theta-theta0)^2`
- proper/improper periodic torsion coefficient: kJ/mol

It contains particle-to-structure ownership, masses, charges, LJ parameters, bonds, angles, torsions, constraints, nonbonded exceptions and optional authoritative type-pair LJ coefficients. `explicitPairTable` exists because a precombined force-field table must not be reverse-engineered into an assumed mixing rule.

`VivoForceFieldLibrary` is the native parameter-library representation. `VivoResidueTemplateAssigner` deterministically types represented residues. `VivoForceFieldCompiler` creates the concrete interaction graph, including 1-2/1-3 exclusions and 1-4 scaling. Missing bonded parameters reject by default.

## 7. AMBER bridge for the CovAngelo preparation boundary

`VivoAmberPrmtop` is a fixed-width `%FLAG` / `%FORMAT` parser for the AMBER7+ sectioned topology format and accepts optional `%COMMENT` records. It does not assume section order or whitespace delimiters for A4 names.

`VivoAmberImporter` imports a standard AMBER parm7/prmtop directly into `VivoMolecularStructureDocument` + `VivoClassicalSystem`. It preserves:

- physical atom identity and residue membership;
- partial charges using the AMBER 18.2223 electrostatic storage scale;
- masses and AMBER atom-type labels;
- massless extra points as classical `virtualSite` particles rather than fake chemical elements;
- explicit bonded and angle terms with conversion to NumiVivo units and potential conventions;
- proper and AMBER-style improper torsion records;
- per-dihedral SCEE/SCNB 1-4 scaling;
- the AMBER excluded-atoms table;
- the complete used Lennard-Jones type-pair A/B table as direct C12/C6 coefficients rather than an assumed Lorentz-Berthelot reconstruction.

`VivoAmberRestart` imports positions and triclinic/orthogonal six-value box data. Coordinates are converted Å -> nm. Presence of restart velocities is detected and the raw values are preserved as a named Wave B boundary; this implementation does not pretend that they already have NumiVivo nm/ps semantics.

This bridge is intentionally aimed at normal AMBER parm7 produced by the AmberTools-style preparation used in the CovAngelo workflow. CHAMBER/AMOEBA extensions are not part of the Wave A execution contract and must not be treated as qualified by the existence of the generic section parser.

## 8. CLI surface

```bash
# Structure import and canonical identity
numivivo structure-import input.pdb --output structure.json
numivivo structure-inspect structure.json
numivivo structure-select structure.json \
  --query 'within(0.35,resname(LIG),false)'

# Structure preparation
numivivo structure-resolve-altloc structure.json --output resolved.json --mapping altloc-map.json
numivivo structure-build-topology resolved.json --output topology.json
numivivo structure-slice topology.json --query 'protein' --output protein.json --mapping protein-map.json

# Native parameter-library path
numivivo forcefield-validate library.json
numivivo forcefield-assign topology.json --library library.json --output assignment.json
numivivo forcefield-compile topology.json --library library.json --output system.json --report compile.json

# Existing AMBER preparation -> NumiVivo execution objects
numivivo amber-import complex.prmtop \
  --restart complex.rst7 \
  --structure complex.structure.json \
  --system complex.classical.json \
  --state complex.initial-state.json \
  --mapping complex.particle-map.json
```

## 9. Wave A boundaries that remain explicit

The following are not silently emulated:

- full stereochemical perception/canonical SMILES generation;
- SDF V3000 and full query chemistry;
- mmCIF writing;
- general small-molecule force-field atom typing / AM1-BCC generation;
- protonation-state prediction;
- arbitrary metal coordination chemistry;
- CHAMBER, AMOEBA and other extended AMBER execution semantics;
- restart-velocity conversion/execution;
- complete virtual-site construction rules for native force-field libraries;
- CMAP and other force-field terms not yet represented by `VivoClassicalSystem`.

The AMBER bridge is the immediate route for ff19SB/ff14SB, GAFF2, OPC and AM1-BCC-prepared systems while NumiVivo replaces those preparation functions incrementally. This prevents Wave B from depending on GROMACS/ACPYPE topology semantics.

## 10. Wave B handoff

Wave B should consume only validated `VivoClassicalSystem` + `VivoClassicalInitialState` objects. Its first execution core should implement GPU-resident positions/velocities/forces, bonded interactions, explicit or mixed LJ interactions, electrostatics, neighbor construction, periodic geometry and integration. The AMBER pair table and virtual-site model added in Wave A are part of that contract and must not be flattened during Metal packing.

The first scientific target remains the preparation/MD portion of the CovAngelo workflow: minimization, NVT, NPT and production sampling of an explicitly solvated protein-ligand system. The paper describes AmberTools-to-GROMACS preparation and a 2 fs periodic MD workflow; NumiVivo now has the source representations needed to begin replacing that execution layer natively.

## References used for interoperability decisions

- CovAngelo preprint, arXiv:2604.10487, software environment and MD workflow.
- AMBER parameter/topology format specification: https://ambermd.org/FileFormats.php
- ParmEd AMBER format implementation/documentation: https://parmed.github.io/ParmEd/
- MDAnalysis AMBER topology parser documentation: https://docs.mdanalysis.org/
