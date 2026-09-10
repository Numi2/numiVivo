# Native global QL scale and count-model refit

Freeze before new native/reference results. Base 1de403d566e631dab61ed0152da2afc9acf41538.
The adaptive moment owner is qualified for the complete original 58-arm family.
This stage estimates a global QL scale and refits counts; it does not fit an
unequal-DF empirical prior or perform a QL cohort hypothesis test.

Implement native robust local linear smoothing with tricube proximity weights,
three Tukey-bisquare reweightings (cutoff six times median absolute residual),
span 0.5 and interpolation between anchors separated by delta=0.01*range(x).
Use stable weighted centered regression, explicit constant-design handling,
stable input-index tie ordering and a bounded neighborhood-visit budget. Retain
empty-weight and constant-fit diagnostics; never silently certify an unidentified
weighted regression. Check the polynomial/constant/outlier/tie/ordering boundary
cases against stats::lowess, and every real update's smoother against the same
reference applied to the actual native quarter-root inputs. Numerical comparison
tolerance is 2e-7 relative to max(1,abs(reference)). Retain any disagreement.

For an explicit pseudobulk family, supplied abundance covariates and supplied
positive trend dispersions, fit each gene with the existing native NB owner.
At scale=1, compute native adjusted residuals. Robustly smooth the fourth root
of their quasi-dispersions over abundance; set the new global scale to at least
one, using the type-7 90th percentile of the smooth, raised to power four.
Repeat this update once more at the new scale while holding initial means fixed.
Then refit the count GLMs once at trend dispersion/global scale, and calculate
final adjusted residuals with the same native moment owner. Preserve all fit,
scale-iteration, eligibility, summation and failure diagnostics. A failed gene
must not silently disappear to make a global family pass. Omit residuals with
unavailable quasi-dispersion from scale learning and report those indices;
reference tiny-DF zero substitution is not imported as observed information.

Qualify both trend arms of all 29 prior full-support Kang, Hagai and Crowell
cases at unchanged counts, feature IDs, donor design, offsets and contrast.
Extract abundance coordinates and the two intermediate reference updates from
pinned edgeR 4.10.5, checking that the reconstructed reference scale/refit agrees
with the already retained original QL fits. These abundance values are explicit
comparison inputs, not a newly qualified native abundance estimator.

Compare scale and residual differences from edgeR descriptively because its
moment approximations differ from native summation. Arbitrate smoothing with
R lowess on exactly the native inputs. Check final native fits using the native
score threshold 1e-7 and independent tighter-reference fits at the exact native
scale; compare identified fitted means and contrast effects within 2e-5 relative
to max(1,abs(reference)). Retain original reference stopping discrepancies and
all failures. All 58 arms are required for complete coverage. Measure a pilot
before scheduling the full family, checkpoint terminal receipts by hashes, and
check workload ownership before remote native runs.

Native controlled tests cover scale floor, exact two-update order, unchanged
initial means between updates, explicit family failures, eligible support,
geometry, representable scaled dispersion and work limits. Build/test success
and conditional numerical agreement do not establish robust QL moderation,
FDR/interval calibration, biological truth, general speed superiority or
million-cell execution. Reference GPL package source remains outside Git;
implement from the mathematical method without copying package code.
