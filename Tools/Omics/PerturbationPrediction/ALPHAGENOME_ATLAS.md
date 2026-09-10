# AlphaGenome Atlas and the NumiVivo development sequence

Reviewed 2026-09-10 after the user supplied the
[DeepMind announcement](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/).
The announced resource contains precomputed single-nucleotide variant effects,
the combined AlphaGenome/AlphaMissense AVI score, and feature attributions.
Its natural role in the existing goal is the variant-to-regulation stage, with
later validation of connections to RNA and cellular state.

The [official API repository at the inspected revision](https://github.com/google-deepmind/alphagenome/tree/aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d)
documents access to Atlas and base-model predictions. The README distinguishes
the Apache-licensed client code from API/data terms. It describes non-commercial
API use and restrictions on training other models, with exceptions for artifacts
specifically covered by permissive terms. The live terms page returned only its
sign-in shell during this review, so the exact rights for a proposed artifact
remain to be checked. No API credentials, downloaded Atlas scores, or trained
NumiVivo integration are claimed here.

## Concrete fit

Use an explicit external-prediction record keyed by reference assembly,
chromosome, coordinate convention, reference and alternate alleles, molecular
readout, tissue/cell ontology, track identity, model/data release, and query/raw
response hashes. Keep the molecular readouts and their units alongside AVI;
a single ranking score cannot supply a quantitative cellular response.
Resolve reference allele and assembly before joining to Ensembl features.
Missing tissue coverage or ambiguous variant-to-gene assignments must remain
missing. Fetch/cache only the selected variants needed by a declared study.

The first independent qualification should compare these external predictions
against measured variant effects under a frozen locus/study split, retaining
coverage and simple baselines. Any learned connection to NumiVivo requires
appropriate data rights and its own held-out biological validation. It cannot
inherit qualification from the upstream predictor's published benchmarks.

## Current boundary

Adamson is a CRISPRi suppression screen, and Norman is a CRISPRa screen. Neither
intervention is equivalent to changing one DNA base. Atlas therefore does not
replace their measured perturbation outcomes or supply a validated mapping from
guide targets to expression changes. The frozen Adamson GO-kernel protocol stays
unchanged; no Atlas information has entered its preparation or scoring.
Finish independent-study validation and the existing single-cell foundation
before extending to regulatory variant prediction.
