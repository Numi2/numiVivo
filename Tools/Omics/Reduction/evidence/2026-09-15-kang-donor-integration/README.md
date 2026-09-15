# Current Kang donor-integration receipt

This receipt records the current native fitted-PCA donor correction on the
public Kang H5AD at source commit `476ad07914bd6edc1a92b075ca2bbc9057afa035`.
The source H5AD SHA-256 is
`1d48c1ff10bcfad7c1bc75863ce702a491c6c71b08e921a5e9754a986a9d19f3`; the
current binary SHA-256 is
`e5fd58fe1ddb38f982dc7145aa4bc04cb87cf4fb1cb80ad3b0c6fb909c6858e4`.
The input has 24,673 cells, 15,706 genes, eight donors and 20 PCA components.
The eight source types are B cells, CD14+ Monocytes, CD4 T cells, CD8 T cells,
Dendritic cells, FCGR3A+ Monocytes, Megakaryocytes and NK cells.

The applied plan uses 100 correction clusters, donor as the covariate, diversity
2, temperature 0.1, ridge 1, seed 7, ten maximum iterations and the declared
work budget of 1,000,000,000. Native integration and
`singlecell-pca-integrate-verify` completed. Membership row-sum error is
`1.9984014443252818e-15`, independent correction error is
`8.881784197001252e-14`, and the maximum ridge residual is
`4.5799824110277385e-15`. The solver stopped on relative objective tolerance;
its eight objective values and seven relative improvements are retained in
`report.json` and `summary.json`.

This is a transductive numerical correction receipt. Original counts and PCA
remain preserved, but there is no held-out donor, prospective mapping, outcome
score, multi-cell-type preservation gate or biological-outcome qualification.
