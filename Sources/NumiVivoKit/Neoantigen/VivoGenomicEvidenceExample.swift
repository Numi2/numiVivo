import Foundation

public enum VivoGenomicEvidenceExample {
    /// All sequence, count and model values here are invented for software use.
    /// This is not the separate HCC1395 interoperability test slice.
    public static func inputs(binding: VivoGenomicCaseBinding, parentInputs: [VivoAtlasVariantInput],
                              implementationSHA256: String) throws -> [String: Data] {
        try VivoAtlasEvidence.require(binding.dataClass == "synthetic" && binding.hlaAlleles.contains("HLA-A*02:01"), "Example requires the synthetic case and HLA-A*02:01.")
        let header = ["Chromosome", "Start", "Stop", "Reference", "Variant", "Junction", "Junction Start", "Junction Stop", "Junction Score", "Junction Anchor", "Transcript", "Gene Name", "Protein Position", "HLA Allele", "Peptide Length", "Epitope Seq", "WT Protein Length", "ALT Protein Length", "Index", "Median IC50 Score", "Tumor RNA Depth", "Tumor RNA VAF", "Normal VAF", "Gene Expression"]
        let index = "0.DEMO_SPLICE.ENST00000000001.DEMO_JUNC.chr4:400-401.D.inframe_splice_site"
        let row = ["chr4", "400", "401", "T", "C", "DEMO_JUNC", "390", "450", "12", "D", "ENST00000000001.1", "DEMO_SPLICE", "2", "HLA-A*02:01", "9", "ACDEFGHIK", "13", "13", index, "225", "20", "0.3", "NA", "NA"]
        let report = Data((header.joined(separator: "\t") + "\n" + row.joined(separator: "\t") + "\n").utf8)
        let regtools = Data("chrom\tstart\tend\tname\tscore\tstrand\tanchor\ttranscripts\tvariant_info\nchr4\t390\t450\tDEMO_JUNC\t12\t+\tD\tENST00000000001\tchr4:400-401\n".utf8)
        let proteins = Data((">WT." + index + "\nMACDEYGHIKLMN\n>ALT." + index + "\nMACDEFGHIKLMN\n").utf8)
        let manifest = VivoSpliceManifest(schema: "numivivo.org/splice-manifest/v1", binding: binding, provenanceMode: "synthetic",
            sourceCitation: "Invented software fixture: no biological measurements, real transcript reconstruction or model execution.",
            sourceRevision: "synthetic-v1", reportStage: "allEpitopes", reportSHA256: try VivoAtlasEvidence.digest(report),
            regtoolsSHA256: try VivoAtlasEvidence.digest(regtools), transcriptsSHA256: try VivoAtlasEvidence.digest(proteins), executionReceiptSHA256: nil)
        var files = ["splice-manifest.json": try VivoCanonicalJSON.encode(manifest), "splice-report.tsv": report, "regtools.tsv": regtools, "transcripts.fa": proteins]
        let splice = try VivoGenomicEvidenceRevisionBuilder.build(binding: binding, parentInputs: parentInputs, files: files, implementationSHA256: implementationSHA256)
        let request = try VivoAtlasEvidence.plan(binding: binding, inputs: parentInputs + splice.spliceReport!.candidates.map(\.variant), requestedScorers: ["SYNTHETIC_RNA", "SYNTHETIC_SPLICE"], ontologyTerms: ["CL:0000000"])
        let requestData = try VivoCanonicalJSON.encode(request)
        var outcomes: [VivoAtlasQueryOutcome] = [], tables: [VivoAtlasTable] = []
        for (i, query) in request.variants.enumerated() {
            let status: VivoAtlasQueryOutcome.Status = i == 1 ? .noData : (i == 2 ? .failed : .available)
            outcomes.append(.init(variantKey: query.key, status: status, diagnostic: status == .failed ? "syntheticDeadlineExample" : nil))
            if status == .available {
                tables.append(.init(scorer: query.chromosome == "chr4" ? "SYNTHETIC_SPLICE" : "SYNTHETIC_RNA", isSigned: true,
                    observations: [.init(variantKey: query.key, metadata: ["gene_name": .string("DEMONSTRATION_ONLY")])],
                    tracks: [["ontology_curie": .string("CL:0000000"), "biosample_name": .string("Invented demonstration tissue")]],
                    rawScores: [[i == 0 ? -0.4 : 0.6]], quantiles: [[i == 0 ? -0.8 : 0.7]]))
            }
        }
        let capture = VivoAtlasCapture(schema: "numivivo.org/atlas-capture/v1", requestSHA256: try VivoAtlasEvidence.digest(requestData),
            mode: "fixture", provider: "google-deepmind/alphagenome-atlas", clientCommit: VivoAtlasEvidence.clientCommit, clientVersion: "0.0.0-synthetic",
            adapterSHA256: implementationSHA256, retrievedAt: "2026-09-08T00:00:00Z", serviceVersion: nil,
            serviceVersionStatus: "not-exposed-by-pinned-api", termsURI: "https://deepmind.google.com/science/alphagenome/terms",
            permittedUse: "noncommercialResearch", trainingPermitted: false, referenceCheck: "synthetic-not-biological",
            referenceSHA256: binding.referenceSHA256, outcomes: outcomes, tables: tables)
        files["atlas-request.json"] = requestData; files["atlas-capture.json"] = try VivoCanonicalJSON.encode(capture)
        return files
    }
}
