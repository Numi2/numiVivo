# AlphaGenome Atlas roadmap assessment

The maintained assessment is [AlphaGenomeAtlas.md](AlphaGenomeAtlas.md), which
links the [implemented native evidence layer and retrieval adapter](Design/ALPHAGENOME_ATLAS.md).
The 2026-09-10 follow-up note initially overlooked that existing implementation;
this pointer consolidates the duplicate without creating another adapter path.

The single-cell connection remains proposed: preserve variant identity, genomic
assembly, tissue context, source scores and provenance; compare regulatory
predictions with independent RNA/ATAC measurements. Existing Atlas evidence does
not qualify held-out perturbation prediction, and its applicable output-use
terms must be retained before considering a training use.
