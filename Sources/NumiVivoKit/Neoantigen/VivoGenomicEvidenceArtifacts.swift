import Foundation

public struct VivoGenomicEvidenceReceipt: Codable, Sendable, Equatable {
    public var schema: String
    public var parentReceipt: VivoNeoantigenReceipt
    public var inputs: [String: VivoFingerprint]
    public var revision: VivoFingerprint
    public var implementation: VivoFingerprint
}
public struct VivoAtlasSelectors: Codable, Sendable {
    public var requestedScorers: [String]
    public var ontologyTerms: [String]
}

public enum VivoGenomicEvidenceArtifacts {
    public static func binding(for parent: VivoNeoantigenReport) throws -> VivoGenomicCaseBinding {
        try VivoAtlasEvidence.require(parent.state != .blocked, "Resolve blocking parent-case findings before adding evidence.")
        guard let rna = parent.manifest.samples.first(where: { $0.role == .tumorRNA }) else { throw VivoGenomicEvidenceError.invalid("Parent case has no RNA sample binding.") }
        return .init(parentReportSHA256: try VivoNeoantigenWorkbench.fingerprint(parent), caseFingerprint: parent.caseFingerprint,
            dataClass: parent.manifest.dataClass.rawValue, assembly: parent.manifest.assembly,
            referenceSHA256: parent.manifest.referenceSHA256, hlaAlleles: parent.manifest.hlaAlleles, tumorRNASampleID: rna.sampleID)
    }
    public static func parentInputs(_ parent: VivoNeoantigenReport) -> [VivoAtlasVariantInput] {
        parent.candidates.map { .init(candidateID: $0.id, chromosome: $0.chromosome, start0: $0.start, stop0: $0.stop,
            reference: $0.sourceFields["Reference"] ?? "NA", alternate: $0.sourceFields["Variant"] ?? "NA") }
    }
    public static func request(parent: VivoNeoantigenReport, selectors: VivoAtlasSelectors,
                               spliceFiles: [String: Data] = [:], implementation: VivoFingerprint) throws -> VivoAtlasRequest {
        let binding = try binding(for: parent)
        var inputs = parentInputs(parent)
        if !spliceFiles.isEmpty {
            let revision = try VivoGenomicEvidenceRevisionBuilder.build(binding: binding, parentInputs: inputs,
                files: spliceFiles, implementationSHA256: implementation.hex)
            inputs += revision.spliceReport?.candidates.map(\.variant) ?? []
        }
        return try VivoAtlasEvidence.plan(binding: binding, inputs: inputs, requestedScorers: selectors.requestedScorers, ontologyTerms: selectors.ontologyTerms)
    }
    public static func publish(parentReceipt: VivoNeoantigenReceipt, files: [String: Data],
                               implementation: VivoFingerprint, store: VivoArtifactStore) async throws -> (VivoGenomicEvidenceReceipt, VivoGenomicEvidenceRevision, VivoNeoantigenReport) {
        let parent = try await VivoNeoantigenArtifacts.verify(parentReceipt, implementation: implementation, store: store)
        let revision = try VivoGenomicEvidenceRevisionBuilder.build(binding: binding(for: parent), parentInputs: parentInputs(parent), files: files, implementationSHA256: implementation.hex)
        var inputs: [String: VivoFingerprint] = [:]
        for name in files.keys.sorted() {
            let media = name.hasSuffix(".json") ? "application/json" : (name.hasSuffix(".tsv") ? "text/tab-separated-values" : "text/plain")
            inputs[name] = try await store.put(data: files[name]!, kind: "neoantigen.genomic-source", mediaType: media).fingerprint
        }
        let stored = try await store.put(data: VivoCanonicalJSON.encode(revision), kind: "neoantigen.evidence-revision", mediaType: "application/json")
        let receipt = VivoGenomicEvidenceReceipt(schema: "numivivo.org/neoantigen-evidence-receipt/v1", parentReceipt: parentReceipt,
            inputs: inputs, revision: stored.fingerprint, implementation: implementation)
        _ = try await store.put(data: VivoCanonicalJSON.encode(receipt), kind: "neoantigen.evidence-receipt", mediaType: "application/json")
        return (receipt, revision, parent)
    }
    public static func verify(_ receipt: VivoGenomicEvidenceReceipt, implementation: VivoFingerprint,
                              store: VivoArtifactStore) async throws -> (VivoGenomicEvidenceRevision, VivoNeoantigenReport) {
        try VivoAtlasEvidence.require(receipt.schema == "numivivo.org/neoantigen-evidence-receipt/v1" && receipt.implementation == implementation && !receipt.inputs.isEmpty && Set(receipt.inputs.keys).isSubset(of: VivoGenomicEvidenceRevisionBuilder.inputNames), "Evidence receipt schema, implementation or input contract differs.")
        let parent = try await VivoNeoantigenArtifacts.verify(receipt.parentReceipt, implementation: implementation, store: store)
        var files: [String: Data] = [:]
        for name in receipt.inputs.keys.sorted() { files[name] = try await store.data(for: receipt.inputs[name]!, maximumBytes: VivoAtlasEvidence.maximumDocumentBytes) }
        let regenerated = try VivoGenomicEvidenceRevisionBuilder.build(binding: binding(for: parent), parentInputs: parentInputs(parent), files: files, implementationSHA256: implementation.hex)
        let stored = try await store.data(for: receipt.revision, maximumBytes: 128 * 1024 * 1024)
        try VivoAtlasEvidence.require(try VivoCanonicalJSON.encode(regenerated) == stored, "Evidence revision does not reconstruct from archived parent, source bytes and implementation.")
        return (regenerated, parent)
    }
    public static func validate(_ review: VivoGenomicEvidenceReview, against revision: VivoGenomicEvidenceRevision) throws {
        try VivoGenomicEvidenceRevisionBuilder.validate(review, against: revision)
    }
    public static func record(_ review: VivoGenomicEvidenceReview, receipt: VivoGenomicEvidenceReceipt,
                              implementation: VivoFingerprint, store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        let (revision, _) = try await verify(receipt, implementation: implementation, store: store)
        try validate(review, against: revision)
        return try await store.put(data: VivoCanonicalJSON.encode(review), kind: "neoantigen.evidence-review", mediaType: "application/json")
    }
}

extension VivoGenomicEvidenceArtifacts {
    public static func spliceJobTemplate(parent: VivoNeoantigenReport) throws -> Data {
        let binding = try binding(for: parent)
        let bindingJSON = try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(binding))
        func placeholder(_ name: String) -> [String: String] { ["path": "/SET/" + name, "sha256": "SET_SHA256"] }
        let inputs = ["annotatedVCF": placeholder("annotated.vcf"), "rnaBAM": placeholder("tumor.bam"), "rnaBAI": placeholder("tumor.bam.bai"),
                      "referenceFASTA": ["path": "/SET/reference.fa", "sha256": binding.referenceSHA256], "referenceFAI": placeholder("reference.fa.fai"), "annotationGTF": placeholder("annotation.gtf")]
        let job: [String: Any] = ["schema": "numivivo.org/splice-job/v1", "binding": bindingJSON, "inputs": inputs,
            "tools": ["regtools": placeholder("regtools"), "pvacsplice": placeholder("pvacsplice")], "resourceFiles": [String: String](),
            "sampleName": "SET_TUMOR_VCF_SAMPLE", "normalSampleName": "SET_NORMAL_VCF_SAMPLE", "rnaStrand": "SET_RF_FR_OR_XS",
            "predictors": ["SET_EXPLICIT_PREDICTOR"], "epitopeLengths": [8, 9, 10, 11], "threads": 1,
            "junctionMinimumReads": 10, "variantDistance": 100, "expressionMinimum": 1.0,
            "resourceVersions": ["regtools": "SET_VERSION", "pvacsplice": "SET_VERSION", "SET_EXPLICIT_PREDICTOR": "SET_VERSION"],
            "sourceCitation": parent.manifest.sourceCitation]
        return try JSONSerialization.data(withJSONObject: job, options: [.sortedKeys, .withoutEscapingSlashes])
    }
}
