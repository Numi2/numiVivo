# AlphaGenome Atlas and NumiVivo's cross-scale direction

Reviewed 2026-09-10 in response to the proposed resource. This is an integration
design note; no Atlas adapter, access credential or biological validation is
claimed.

DeepMind's [2026-09-08 announcement](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/)
describes precomputed molecular-effect predictions for roughly nine billion human
single-nucleotide variants. Its AVI score combines AlphaGenome and AlphaMissense,
with feature attributions and regulatory motifs for interpretation. This fits
the variant → regulation portion of NumiVivo's long-term objective.

The useful first connection would be a bounded external variant-evidence import:
retain reference assembly, chromosome, allele identity, coordinate convention,
model/dataset version, requested tissue or cell-type ontology, molecular output,
score definition, attribution units, source URL and response hash. Exact reference
alleles and identifiers must be verified before joining to gene/feature spaces.
These are proposed adapter requirements; actual API fields must be checked
against a pinned client revision before implementation.

Keep regulatory predictions distinct from measured RNA, ATAC, protein abundance
and perturbation outcomes. AVI prioritization alone does not supply a signed,
tissue-specific RNA effect or a validated kinetic parameter. Each downstream
link needs an explicit hypothesis and independent measurement. For perturbation
work, an appropriate later experiment would compare a baseline against a
permitted Atlas-informed method on held-out variants and biological contexts,
with gene/locus leakage checks and unchanged outcome criteria.

The official [API repository](https://github.com/google-deepmind/alphagenome)
provides Atlas access as well as model inference. Its current documentation
requires an API key and describes separate output-use terms, including limits on
non-commercial outputs being used to train other machine-learning models. Any
future adapter must record the applicable artifact terms and permitted use;
the client code's Apache license is not the prediction-data license. A bulk
download or training-data import is not part of this change.

Near-term work remains the verified single-cell pipeline and unseen-perturbation
prediction. Atlas is a candidate external evidence source for the later genomics
and cross-scale stages; it does not replace missing donor strata, cell labels,
integration validation or experimental perturbation targets.
