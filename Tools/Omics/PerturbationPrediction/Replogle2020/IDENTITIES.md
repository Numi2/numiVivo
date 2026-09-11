# Replogle UPR: intended-gene identity evidence

All 30 original non-control guide labels now have a supported nominal gene
identity for the frozen RNA-response benchmark. This resolves the descriptor
join; it does not verify full guide sequences, genome-wide specificity,
knockdown efficacy or causal attribution of an observed RNA response.

## Evidence chain

The original GEO assignments retain `sg`-prefixed guide labels. A
[published FBA reproduction of this UPR experiment](https://fba.readthedocs.io/en/latest/tutorials/crispr_screening/PRJNA609688/tutorial.html)
lists the same 32 labels, including the two controls, alongside 19-base guide
prefixes. FBA explicitly attributes its sequence table to Replogle Supplementary
Table 2 and explains that sequences were truncated to equal length. Its source
is pinned to [FBA commit df9e87a](https://github.com/jlduan/fba/blob/df9e87ac4ebadc91b216814350e2c4d1d35f048c/docs/tutorials/crispr_screening/PRJNA609688/tutorial.rst).
This is an independent reproduction, **not a newly retrieved copy of the
original Replogle supplement**. The original supplement remains unavailable in
this run. FBA's MIT license is retained with the captured source.

All 32 source labels join one-to-one by the explicit `sg` prefix convention;
no guide is dropped or merged. The two control roles retain their original
Replogle-paper authority. For 29 targets, the reproduced gene symbol identifies
one original RNA Ensembl feature. The remaining label, `SRPR`, is not present as
a symbol in the RNA feature table. Its single
[HGNC previous-symbol record](https://rest.genenames.org/fetch/prev_symbol/SRPR)
explicitly identifies `SRPRA`, `HGNC:11307`, `ENSG00000182934`; that exact feature
is present once in the deposited RNA table. This is a recorded nomenclature
history join, not a guessed alias or an alias query to the GO service.

Current Ensembl GRCh38 gene records and reference sequence windows were captured
for all 30 exact Ensembl IDs. Each window covers the annotated gene span plus
2 kb at each end. Every reproduced prefix's final 18 bases matches exactly once
among all 30 captured windows, on one strand of its nominated gene locus. The
two control prefixes have no such match in these windows. No expression count,
RNA effect, score or GO similarity selects the mapping.

## What the sequence audit does and does not show

Fourteen full 19-base prefixes match reference DNA literally. Sixteen have an
exact 18-base core but differ at the initial G. An initiating G is consistent
with the [U6 expression design described by Addgene](https://blog.addgene.org/crispr-cas9-faqs-answered),
but that explanation is an inference for these particular constructs. The
truncated table does not establish their full spacer lengths.

The initial literal-prefix audit is retained, including its sixteen failures
and the missing exact `SRPR` symbol. Subsequent fixed-position PAM audits also
failed: `SRPR` has no NGG at position 20, and `TMED10` has none at positions 20
or 21, despite both matching all 19 prefix bases. PAMs following positions 19,
20 and 21 are therefore recorded as diagnostics; they are not used to assert a
reconstructed full guide sequence. Those failed assumptions remain in the archive.

The admitted claim is the **nominal intended gene**, supported jointly by the
published guide/gene table, exact stable gene identity and reference core match.
It does not depend on claiming that all 19 bases match or that a full protospacer
has been reconstructed. The reference search is restricted to these candidate
loci, not the whole genome. This is neither a uniqueness guarantee across the
genome nor evidence that the guide perturbs only its nominal target. Current
GRCh38 reference captures are not a byte-identical reconstruction of the original
Cell Ranger GRCh38-1.2.0 reference package.

These metadata audits followed cohort preparation but preceded model fitting
and scoring. They do not change the RNA cohort, predictor, regularization,
annotation inclusion rule or primary endpoint in the [frozen protocol](PROTOCOL.md).

## GO capture

All 30 exact Ensembl identities have usable direct GO terms under the unchanged
BP/MF/CC rule: human taxon 9606 and exact returned Ensembl identity, all evidence
classes, exclude NOT and root terms, no ancestry expansion or alias guessing.
The MyGene build is `20260906`. Twenty-eight responses reuse the earlier exact,
hashed captures under the same build; `TIMM23` and `IER3IP1` are new requests.
Reusing gene annotations does not resolve the separate Adamson experimental roles.

No retained direct annotation cites PMID 27984733 or 32231336. That checks only
the captured citation fields; current annotation knowledge may still incorporate
these studies indirectly and is not temporally independent. The descriptor
receipt binds every response, inclusion/exclusion record, feature ID and source
identity. All 30 guides remain in all five gemgroups.

FBA also associates GSM4367980 with barcode suffix `-4`, so experiment number
cannot be assumed to equal gemgroup number. Results retain the five original
gemgroup IDs; named platform mappings have not been generalized from that one
example. They are technical conditions, not independent biological replicates.

The identity and descriptor captures, reference sequences, licenses, initial
failures and final audit scripts are retained in the
[prediction evidence archive](evidence/2026-09-11-target-kernel/manifest.json).
