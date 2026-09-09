import Foundation

/// A new revision never modifies the originating neoantigen report or its reviews.
/// Raw input hashes refer to separately archived bytes in the existing store.
public struct VivoGenomicEvidenceRevision: Codable, Sendable, Equatable {
    public var schema: String
    public var binding: VivoGenomicCaseBinding
    public var implementationSHA256: String
    public var sourceSHA256: [String: String]
    public var parentCandidateIDs: [String]
    public var atlasRequest: VivoAtlasRequest?
    public var atlasCapture: VivoAtlasCapture?
    public var spliceReport: VivoSpliceReport?
    public var limitations: [String]
    public var candidateIDs: [String] { parentCandidateIDs + (spliceReport?.candidates.map(\.id) ?? []) }
}

public enum VivoGenomicEvidenceRevisionBuilder {
    public static let inputNames: Set<String> = ["atlas-request.json", "atlas-capture.json", "splice-manifest.json", "splice-report.tsv", "regtools.tsv", "transcripts.fa", "splice-execution.json"]
    public static func build(binding: VivoGenomicCaseBinding, parentInputs: [VivoAtlasVariantInput],
                             files: [String: Data], implementationSHA256: String) throws -> VivoGenomicEvidenceRevision {
        try VivoAtlasEvidence.validate(binding)
        try VivoAtlasEvidence.require(VivoAtlasEvidence.isDigest(implementationSHA256) && !files.isEmpty && Set(files.keys).isSubset(of: inputNames), "Unexpected evidence input or missing implementation identity.")
        try VivoAtlasEvidence.require(files.values.allSatisfy { $0.count <= VivoAtlasEvidence.maximumDocumentBytes }, "Evidence document exceeds size limit.")
        try VivoAtlasEvidence.require(parentInputs.count <= 10_000 && Set(parentInputs.map(\.candidateID)).count == parentInputs.count && parentInputs.allSatisfy { VivoAtlasEvidence.isDigest($0.candidateID) }, "Invalid parent candidate bindings.")
        var spliceReport: VivoSpliceReport?
        let spliceNames = ["splice-manifest.json", "splice-report.tsv", "regtools.tsv", "transcripts.fa"]
        if spliceNames.contains(where: { files[$0] != nil }) {
            try VivoAtlasEvidence.require(spliceNames.allSatisfy { files[$0] != nil }, "Splice import needs a manifest, unaggregated report, RegTools TSV and paired protein FASTA.")
            let manifest = try VivoGenomicDocuments.decode(VivoSpliceManifest.self, from: files["splice-manifest.json"]!)
            try VivoAtlasEvidence.require(manifest.binding == binding, "Splice evidence belongs to a different case/report.")
            if let digest = manifest.executionReceiptSHA256 {
                guard let bytes = files["splice-execution.json"] else { throw VivoGenomicEvidenceError.invalid("Referenced splice execution receipt is missing.") }
                try VivoAtlasEvidence.require(try VivoAtlasEvidence.digest(bytes) == digest, "Splice execution receipt hash differs.")
                // The receipt is an externally recorded assertion, not a signature
                // or proof that tool execution or sample authentication occurred.
                let receipt = try VivoGenomicDocuments.decode(VivoSpliceExecutionReceipt.self, from: bytes)
                try receipt.validate(against: manifest)
            } else {
                try VivoAtlasEvidence.require(files["splice-execution.json"] == nil, "Unreferenced execution receipt.")
            }
            spliceReport = try VivoSpliceEvidence.analyze(manifest: manifest, manifestData: files["splice-manifest.json"]!, reportTSV: files["splice-report.tsv"]!, regtoolsTSV: files["regtools.tsv"]!, transcriptsFASTA: files["transcripts.fa"]!)
        } else {
            try VivoAtlasEvidence.require(files["splice-execution.json"] == nil, "Execution receipt has no splice import.")
        }
        let inputs = parentInputs + (spliceReport?.candidates.map(\.variant) ?? [])
        try VivoAtlasEvidence.require(Set(inputs.map(\.candidateID)).count == inputs.count, "Parent/splice candidate identifier collision.")
        var request: VivoAtlasRequest?, capture: VivoAtlasCapture?
        if let bytes = files["atlas-request.json"] {
            request = try VivoGenomicDocuments.decode(VivoAtlasRequest.self, from: bytes)
            let expected = try VivoAtlasEvidence.plan(binding: binding, inputs: inputs, requestedScorers: request!.requestedScorers, ontologyTerms: request!.ontologyTerms)
            try VivoAtlasEvidence.require(request == expected, "Atlas request does not reconstruct from this exact parent and splice revision. Regenerate it after adding splice candidates.")
            if let data = files["atlas-capture.json"] {
                capture = try VivoGenomicDocuments.decode(VivoAtlasCapture.self, from: data)
                try VivoAtlasEvidence.validate(capture: capture!, request: request!, requestData: bytes)
            }
        } else {
            try VivoAtlasEvidence.require(files["atlas-capture.json"] == nil, "Atlas capture has no exact originating request.")
        }
        try VivoAtlasEvidence.require(request != nil || spliceReport != nil, "No supported evidence revision inputs.")
        let hashes = try files.mapValues(VivoAtlasEvidence.digest)
        return .init(schema: "numivivo.org/neoantigen-evidence-revision/v1", binding: binding,
            implementationSHA256: implementationSHA256, sourceSHA256: hashes,
            parentCandidateIDs: parentInputs.map(\.candidateID), atlasRequest: request, atlasCapture: capture,
            spliceReport: spliceReport, limitations: VivoAtlasEvidence.limitations + VivoSpliceEvidence.limitations)
    }
}

/// Recorded by the managed external runner. Its hashes are verified during
/// import, but an unsigned receipt does not establish trusted execution.
public struct VivoSpliceExecutionReceipt: Codable, Sendable, Equatable {
    public struct Step: Codable, Sendable, Equatable {
        public var tool: String
        public var executableSHA256: String
        public var arguments: [String]
        public var exitCode: Int
        public var stdoutSHA256: String
        public var stderrSHA256: String
    }
    public var schema: String
    public var binding: VivoGenomicCaseBinding
    public var inputSHA256: [String: String]
    public var outputSHA256: [String: String]
    public var steps: [Step]
    public var resourceVersions: [String: String]
    public var startedAt: String
    public var finishedAt: String
    public func validate(against manifest: VivoSpliceManifest) throws {
        try VivoAtlasEvidence.require(schema == "numivivo.org/splice-execution/v1" && binding == manifest.binding, "Splice execution belongs to another case.")
        try VivoAtlasEvidence.require(outputSHA256["splice-report.tsv"] == manifest.reportSHA256 && outputSHA256["regtools.tsv"] == manifest.regtoolsSHA256 && outputSHA256["transcripts.fa"] == manifest.transcriptsSHA256 && inputSHA256["referenceFASTA"] == binding.referenceSHA256, "Splice execution inputs/outputs do not bind the supplied files.")
        try VivoAtlasEvidence.require(["annotatedVCF", "rnaBAM", "rnaBAI", "referenceFASTA", "referenceFAI", "annotationGTF"].allSatisfy { VivoAtlasEvidence.isDigest(inputSHA256[$0] ?? "") }, "Incomplete splice execution input provenance.")
        try VivoAtlasEvidence.require(steps.map(\.tool) == ["regtools", "pvacsplice"] && steps.allSatisfy { $0.exitCode == 0 && VivoAtlasEvidence.isDigest($0.executableSHA256) && VivoAtlasEvidence.isDigest($0.stdoutSHA256) && VivoAtlasEvidence.isDigest($0.stderrSHA256) && !$0.arguments.isEmpty && $0.arguments.count <= 256 && $0.arguments.allSatisfy { !$0.isEmpty && $0.utf8.count <= 4096 && !$0.contains("\u{0}") } }, "Incomplete or failed RegTools/pVACsplice execution.")
        try VivoAtlasEvidence.require(!resourceVersions.isEmpty && resourceVersions.count <= 128 && resourceVersions.allSatisfy { VivoAtlasEvidence.validToken($0.key) && VivoAtlasEvidence.validToken($0.value) && $0.value.lowercased() != "latest" }, "Exact external resource versions are required.")
        guard let start = ISO8601DateFormatter().date(from: startedAt), let finish = ISO8601DateFormatter().date(from: finishedAt), finish >= start else { throw VivoGenomicEvidenceError.invalid("Invalid splice execution timestamps.") }
    }
}

/// Deliberately separate from parent-only reviews: these decisions bind all of
/// the new evidence, and cannot be applied to a differently reconstructed revision.
public struct VivoGenomicEvidenceReview: Codable, Sendable, Equatable {
    public struct Decision: Codable, Sendable, Equatable {
        public enum Disposition: String, Codable, Sendable { case retainForResearch, exclude, deferReview }
        public var candidateID: String
        public var disposition: Disposition
        public var rationale: String
    }
    public var schema: String
    public var revisionSHA256: String
    public var reviewerID: String
    public var reviewedAt: Date
    public var decisions: [Decision]
}
extension VivoGenomicEvidenceRevisionBuilder {
    public static func validate(_ review: VivoGenomicEvidenceReview, against revision: VivoGenomicEvidenceRevision) throws {
        try VivoAtlasEvidence.require(review.schema == "numivivo.org/neoantigen-evidence-review/v1" && review.revisionSHA256 == (try VivoAtlasEvidence.digest(VivoCanonicalJSON.encode(revision))), "Review targets another evidence revision; old parent-only reviews do not approve new evidence.")
        let ids = Set(revision.candidateIDs)
        try VivoAtlasEvidence.require(VivoAtlasEvidence.validToken(review.reviewerID) && review.reviewedAt.timeIntervalSince1970.isFinite && !review.decisions.isEmpty && review.decisions.count <= ids.count, "Invalid reviewer, timestamp or decision count.")
        var seen = Set<String>()
        for decision in review.decisions {
            try VivoAtlasEvidence.require(ids.contains(decision.candidateID) && seen.insert(decision.candidateID).inserted && !decision.rationale.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty && decision.rationale.utf8.count <= 4096 && !decision.rationale.contains("\u{0}"), "Research decisions require unique exact candidate IDs and bounded rationales.")
        }
    }
}
