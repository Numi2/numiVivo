# Predeclared ridge-4 integration candidate

This is a single a priori global-correction candidate on the complete HIRISA
PCA bundle. It changes only the fixed donor-effect ridge penalty from 1 to 4;
the seed, 100 soft clusters, diversity, temperature, iteration limit and
relative tolerance are unchanged. The candidate was frozen before its
preservation outcome was read.

The native Apple-silicon publish and reconstruction both completed on all
**1,612,594 cells**. The output uses the current `f1b882f1` implementation and
the pinned HDF5 library. The independent SVD oracle then passed all 114
donor-held-out folds with maximum coefficient error `5.689893001203927e-16`
and zero confusion error.

## Preservation result

The absolute comparison uses the original-PCA baseline, the frozen author
`celltype.l2` ledger and the existing 0.05 mean/0.10 maximum donor-loss
margins. Author labels are recoverability controls, not biological ground truth.

| Candidate | Supported, sensitive comparisons | Failures | Rare failures | Mean sensitive recall loss | Complete gate |
| --- | ---: | ---: | ---: | ---: | --- |
| native global ridge 1 | 146 | 37 | 10 | -0.00380961 | fail |
| ridge-4 global candidate | 146 | 37 | 10 | -0.00382354 | fail |

Ridge 4 leaves the native failure count unchanged and slightly worsens the
maximum fold loss (`0.6156462585` versus `0.6153370439`). It is therefore
**rejected for promotion**. The small aggregate improvement does not establish
marker, rare-cell, response or phenotype preservation.

The full compressed comparison and oracle, plus source, plan, executable and
output hashes, are in the [compact evidence archive](evidence/2026-09-14-ridge4/manifest.json).
The native output remains on the Mac mini under
`/Users/n/numivivo-hirisa-20260910/integration-ridge4-20260914/native` for
replay; it is not embedded in this repository.

## Reproduce

With the qualified binary and HDF5 dependency, run:

```sh
numivivo singlecell-pca-integrate \
  /Users/n/numivivo-hirisa-20260910/integration-ridge4-20260914/input \
  --plan /Users/n/numivivo-hirisa-20260910/integration-ridge4-20260914/plan.json \
  --output /new/ridge4/native
numivivo singlecell-pca-integrate-verify /new/ridge4/native
```

Then run the frozen absolute retention scorer against the original-PCA
baseline. Keep the complete result and independent oracle even when the gate
fails; this experiment is development evidence, not a biological qualification.
