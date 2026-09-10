# Native profile checks and reference optimization sensitivity

The [stage diagnosis](../DispersionAudit/README.md) found different native and
DESeq2 dispersion curves. This follow-up establishes that **the pinned default
DESeq2 gene-wise estimates miss better in-bounds adjusted likelihoods on these
null data**. Matching that reference is therefore not a sufficient reason to
change the native estimator. Native inference remains experimental; its twenty
original Hagai sham calls are not removed or explained away by this finding.

## Reference objective checks

The [protocol](PROTOCOL.md) checks every eligible gene in all twenty original
splits: 195,058 gene/split evaluations, with overlapping genes and cells across
splits. Original DESeq2 gene-wise estimates reproduce exactly. The package's
own grid routine is evaluated with its actual previously fitted means, unchanged
counts, paired design and Cox–Reid correction. Python independently computes
each signed objective difference using an integer-count rising-factorial identity
and QR, avoiding near-Poisson log-gamma cancellation.

| Per-split diagnostic | Kang | Hagai |
|---|---:|---:|
| Better grid objectives within the package's bounds | 3,639–3,765 | 3,616–3,764 |
| Of those, original estimates below 1e-6 | 3,608–3,727 | 3,493–3,648 |
| Largest in-bounds log-objective improvement | 10.808 | 5.391 |

An improvement must exceed the declared 1e-4 threshold. Nearly all improved
boundary estimates were initialized at 1e-8. The complete tables retain initial
values, iteration counts, means, original/grid estimates and signed objective
differences. This is direct evidence of missed better objectives in these pinned
runs, not a claim that every default DESeq2 analysis has this behavior.

The raw grid can also return worse values or values outside the wrapper's stated
bounds. Those outcomes are retained, and are not admitted to the sensitivity
candidate. A grid search is not treated as a global-optimality certificate.

## Native numerical checks

The exact current FP64 QR/NB sources were compiled on the physical M4 Pro using
Swift 6.3.3. Selection was frozen before evaluating profiles: sixteen largest
native-interior/reference-boundary disagreements from seed 1 of each study,
plus the twenty original native Hagai sham calls, giving 52 gene/split cases.
Checks cover the original/native/grid-selected dispersions and a 41-point
logarithmic grid over [1e-8,100], plus each independently evaluated native optimum.

All **2,305 valid point checks** pass unchanged numerical tolerances. Maximum
relative mean error is 8.13e-7, effect error 7.46e-7, standard-error error 1.32e-7,
stable log-likelihood error 6.12e-10 and adjusted-profile error 7.99e-8. The largest
independent grid objective above the native selected profile is 4.11e-8, below
the declared 2e-5 acceptance tolerance. This provides bounded numerical evidence,
not a proof of global optimization for arbitrary data or calibrated Wald tests.

The initial statsmodels-only checker retained twelve unqualified high-dispersion
comparisons. Refining its coefficient scores with a separate SciPy root solve
resolves them without changing tolerances or native code. One DESeq2 grid value,
4.741331e-9, lies below the native domain; native correctly rejects it. The initial
failed checker output/source and that rejection remain in the evidence.

The native diagnostic executable took 0.64 seconds and peaked at 33,259,520 bytes
RSS on this run. This is a small CPU diagnostic, not an end-to-end single-cell
performance or GPU benchmark.

## Downstream sensitivity and retained failures

The separate [sensitivity protocol](SENSITIVITY_PROTOCOL.md) admits a grid estimate
only when it improves the checked objective by more than 1e-4 and stays inside
[1e-8,max(10,sampleCount)]. It then reruns the original trend/MAP/Wald stages on
the unchanged eligible family. The staged baseline reproduces original final
results to at most 4.98e-14 absolute error. BH is checked independently.

Hagai completes all ten requested parametric sensitivities. Its original DESeq2
baseline has one BH<0.05 sham call in one split; the objective-checked sensitivity
has **25 calls across nine splits**. Native retains its previously published
twenty calls across eight splits. Totals count gene/split calls, not unique genes.
Better optimization therefore does not imply better null calibration, and the
original reference's lower call count is not a sufficient calibration standard.

**Nine Kang parametric sensitivities failed their declared method constraint:**
seeds 1, 2, 3, 5, 6, 7, 8, 9 and 10 automatically fell back to a local trend inside
DESeq2. Only seed 4 completed the requested parametric candidate. All resulting
outputs (zero calls, including the fallback outputs) and messages are preserved,
but the ten-split parametric Kang sensitivity is not qualified. Every completed
baseline/candidate coefficient fit reports convergence. The final summary status
is `completed-with-parametric-sensitivity-failures`.

The first sensitivity adapter omitted the numeric design argument from the
staged Wald call. All twenty failed attempts are retained before its repair to
pass `modelMatrix=x`. A later diagnostic-only rerun collected coefficient
convergence and reproduced every previously completed inference table exactly.
Disk-pressure failures and hash-verified duplicate consolidation are also logged.

## Reproduction and next work

Use the original null benchmark and stage-audit roots plus their pinned Python/R
environments. Run native compilation/evaluation on the owning Mac mini:

```sh
python reference.py --root /runs/null --stages /runs/stages --out /new/profiles \
  --r-library /qualified/R-library
bash build.sh /new/native-build
/new/native-build/profile-audit /new/profiles/native-input.json /new/profiles/native.json
python check_native.py --root /new/profiles --refine-reference
python sensitivity.py --root /new/profiles --null-root /runs/null \
  --r-library /qualified/R-library
python summarize.py --root /new/profiles
```

The published `summarize.py` also checks retained diagnostic rerun tables when
`initial-complete-sensitivity` exists, or when `--prior-sensitivity` is supplied;
otherwise that supplemental rerun check is explicitly unavailable.
The [evidence manifest](evidence/2026-09-10/manifest.json) binds every reference
table, native input/output, numerical check, failure and sensitivity result.

No native solver repair is supported by these selected numerical checks. The
next inference work must address statistical calibration and robustness using
objective-checked references and independent calibration/power evidence, rather
than tuning native outputs to an inadequately optimized reference. The broader
single-cell and cross-scale goal remains open.
