# Native chemistry example

This example executes native Gaussian integrals -> RHF -> embedded Hamiltonian -> MP2/FCI -> exact CI orbital information. The default fixture is H2/STO-3G at 1.4 Bohr. It is not the paper's reaction, a drug-discovery result or a chemical-accuracy claim.

Run the full-library example on the supported Apple toolchain:

```sh
swift run --package-path Examples/native-chemistry native-chemistry --h2 h2-report.json
```

The example also accepts a request JSON in place of `--h2`. Its shape is the `request` object in the generated report: `system`, `basis`, `scf`, and `budget`. Internal positions are Bohr and energies Hartree. The example requires RHF because it sends a common spatial-orbital Hamiltonian to restricted MP2 and CI. Native UHF is available separately through `VivoHartreeFock`.

The generated report includes explicit inputs, converged HF, the integral-derived embedded Hamiltonian, MP2, FCI state, and exact orbital-information matrices.
