# Native VCF import and Atlas projection (2026-09-14)

The source revision `2115751232a22ddb28c3119df3e0daf55bd7752f` was built and run on the physical Apple M4 Pro Mac mini. `swift test --filter VivoVCFTests` passed all four tests. The CLI then ran `genomic-vcf-plan` over the preserved fixture and emitted `variant-import.json` plus `atlas-projection.json`: four imported alternate records, two GRCh38 primary-contig SNVs eligible for the current Atlas request identity, and two machine-readable exclusions (an indel and a symbolic allele).

The reader hashes exact source bytes before CRLF normalization, preserves INFO/FORMAT/genotype text and multiallelic identity, requires an explicit GRCh37/GRCh38 assembly and reference SHA-256, and rejects unbounded or malformed input. The projection admits only GRCh38 primary chr1-22/X/Y SNVs. Unsupported assembly, contigs and allele classes remain explicit exclusions; no liftover, variant calling, phasing, genotype interpretation, Atlas network call or biological prediction is performed.

The fixture is synthetic and this receipt is software/provenance evidence. It does not qualify a real cohort, variant effects, regulatory activity, protein or cellular response, phenotype, tissue outcome or clinical utility. `execution.json` binds the command, machine, source hashes, output counts and scope; `manifest.json` binds every evidence file.
