# AlphaGenome Atlas integration assessment

Assessed 2026-09-09 and sources rechecked 2026-09-11: the [DeepMind announcement, 2026-09-08](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/)
and the [official API repository](https://github.com/google-deepmind/alphagenome).
The single-cell integration below remains proposed. The repository now also has
a [native Atlas evidence layer and external retrieval adapter](Design/ALPHAGENOME_ATLAS.md)
for bounded public-reference genomic research. That implemented path is separate
from RNA/ATAC integration and biological qualification.

## Fit to the roadmap

Atlas precomputes molecular effects for approximately nine billion human
single-nucleotide variants. Its AVI score combines AlphaGenome and AlphaMissense
predictions; feature attributions and sequence motifs support interpretation.
The API exposes Atlas alongside model predictions for expression, splicing,
chromatin features and contact maps. This is directly relevant to the planned
variant → regulation → RNA/cell-state connection and future RNA/ATAC comparisons.

For the current [biological prediction decision](BiologicalPrediction.md#available-information-does-not-imply-a-validated-outcome),
Atlas adds predicted molecular effects for identified DNA variants. It does not
supply missing cytokine dose/reagent metadata, establish equivalence between
ambiguous RNA feature names, or provide measured cellular and tissue outcomes.
The Parse IFN-beta admission and independent RNA-response validation therefore
retain their existing input and outcome requirements.

The proposed role is a versioned external prediction source. AVI can prioritize
variants; modality-specific effects can form testable hypotheses. An AVI score
must not be interpreted as an expression fold change, protein concentration,
reaction rate or tissue outcome. Downstream mechanistic links need separate
models and evidence. These predictions do not validate native differential
expression, donor integration, or held-out perturbation prediction.

The completed [Replogle experiment](../Tools/Omics/PerturbationPrediction/Replogle2020/RESULTS.md)
uses nominal guide-to-gene identities resolved by its separate
[identity audit](../Tools/Omics/PerturbationPrediction/Replogle2020/IDENTITIES.md),
without Atlas input. Full guide sequences and genome-wide specificity remain
unverified. CRISPR interference and a single-base DNA substitution are different
interventions; treating an Atlas variant score as a measured guide response
would not validate the frozen target predictor. The
useful next Atlas experiment is a separate, explicitly identified variant panel
with measured regulatory outcomes. This is an integration proposal, not an
additional NumiVivo prediction result.

## Small first integration

1. Query a bounded, public variant panel through the official Atlas API rather
   than attempting to mirror the announced petabyte-scale resource.
2. Persist exact variant identity, genome assembly, coordinate convention,
   reference/alternate alleles, model/dataset version where provided, retrieval
   time, tissue/cell ontology, target feature IDs, score definitions, request,
   source response and content hashes. Missing source metadata must remain
   explicit; do not invent a version or assume assembly compatibility.
3. Store predictions separately from observed count/assay matrices. Resolve
   gene IDs and biological context explicitly when joining RNA or ATAC evidence.
   Retain AVI and its attributions alongside modality-specific scores.
4. Evaluate a predeclared independent experimental panel with negative controls,
   coverage, direction agreement and appropriate ranking/calibration metrics.
   Audit overlap with model training/evaluation data before claiming held-out
   generalization. Preserve unsupported variants and contexts in coverage reports.

Atlas's exhaustive single-letter scope does not establish exhaustive coverage of
indels, structural variation, haplotypes or combinatorial perturbations. Consumers must distinguish an Atlas lookup from a fresh model prediction.

## Access boundary

The official API README states that access is free for non-commercial use,
requires an API key, and is rate-limited. It also states that outputs and Atlas
information generally may not train other machine-learning models, with exceptions
for specifically designated permissive downloadable artifacts under the terms.
Therefore do not turn retrieved predictions into perturbation-model training
targets by default. Record the applicable artifact/access terms before choosing
that use. The announcement says commercial Atlas availability on Google Cloud
is forthcoming; base-model commercial availability is a separate offering.

No credentials were requested, service calls made, package installed, or Atlas
predictions imported during this assessment. The existing adapter pins the SDK,
checks GRCh38 reference bases, archives responses, and distinguishes available,
missing and failed results. Live service access and independent biological
qualification remain unverified by this assessment. The next single-cell link
is a bounded experimental RNA/ATAC comparison using that evidence foundation;
the native single-cell evidence work remains the immediate development priority.
