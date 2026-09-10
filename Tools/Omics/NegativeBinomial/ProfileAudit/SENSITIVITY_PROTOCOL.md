# Objective-checked reference sensitivity

Declared after fixed-mean grid checks found improved gene-wise objectives in
all twenty null splits, before evaluating revised downstream inference.

For every split, first reproduce the original DESeq2 result through its explicit
gene-wise/trend/MAP/Wald stages. Compare all original final dispersions, effects,
standard errors and probabilities before admitting a candidate result.

The sole candidate change is to replace a gene-wise estimate with the already
audited grid value when its fixed-mean adjusted objective improves by more than
1e-4 and the value is within [1e-8,max(10,sampleCount)]. Preserve the original
estimate otherwise, including cases where the grid is worse or outside bounds.
Then rerun the same parametric trend, prior, MAP, outlier and Wald stages on the
complete original eligible family, with unchanged counts/design/normalization.
Record fit warnings, failures, changed genes, trend/prior parameters and all
gene results. Do not change filtering, p-value thresholds or choose splits.

This is a named optimization-sensitivity reference, not unmodified DESeq2 or a
replacement for its frozen baseline. It does not certify global optimization,
general false-discovery control or alternative-model power, and it does not
change native inference.
