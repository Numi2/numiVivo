# Reference precision and stopping sensitivity

The first complete independent pass retains two failed Hagai sham checks
(the same feature under the two policies). NumPy longdouble is only 64 bits
in this Apple environment. Subtracting large log terms causes a 2.46e-7
probability discrepancy. An 80-digit calculation at the exact saved means
finds native probability error of 3.93e-13 for that case.

The reference refinement recomputes any case whose initial LR discrepancy
exceeds 2e-7 or probability discrepancy exceeds 2e-8 using 80-digit mpmath.
Those triggers are one tenth of the original acceptance tolerances. Every
initial difference and refined difference is retained; the original 50 checks
and their checker source remain. No native code or tolerance changes for this
reference arithmetic correction. Check all 50 original cases again.

The first pinned edgeR comparison retains 35 probability disagreements across
Kang default (5), Kang active donor (20) and Hagai (10). The seven Crowell
comparisons pass. Some disagreeing edgeR full fits have scaled scores around
3e-4, whereas the native fits satisfy 1e-7. Default edgeR stopping therefore
needs a separate sensitivity check before attributing those differences to
native inference.

Freeze this sensitivity before running it: rerun all ten reference treatment
analyses with full-model glmFit tol=1e-10 and maxit=1000. Keep every count,
dispersion, offset, design, zero prior count, contrast and original LR/p-value
tolerance. glmLRT's own null fitting remains unchanged. Retain full/null score
diagnostics and compare both original and tightened BH decisions. Original
default-reference failures remain failures; the sensitivity is supplementary
numerical evidence, not a calibration or method-promotion gate.
