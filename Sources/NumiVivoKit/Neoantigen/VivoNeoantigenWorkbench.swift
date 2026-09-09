import Foundation

public enum VivoNeoantigenError: Error, Sendable, CustomStringConvertible {
    case invalid(String)
    public var description: String { switch self { case .invalid(let message): return message } }
}

/// Reference-data research only. These records are declarations, not genetic identity
/// verification, consent management, or evidence that an external program executed.
public struct VivoNeoantigenCase: Codable, Sendable, Equatable {
    public enum DataClass: String, Codable, Sendable { case synthetic, publicReference }
    public struct Sample: Codable, Sendable, Equatable {
        public enum Role: String, Codable, Sendable, CaseIterable { case tumorDNA, matchedNormalDNA, tumorRNA }
        public var role: Role
        public var sampleID: String
        public var subjectID: String
        public var assembly: String
        public var referenceSHA256: String
        public var sourceSHA256: String
    }
    public struct ExternalRun: Codable, Sendable, Equatable {
        public var tool: String
        public var version: String
        public var format: String
        public var executionImageSHA256: String
        public var arguments: [String]
        public var modelVersions: [String: String]
        public var sampleInputs: [String: String]
        public var annotatedVCFSHA256: String
        public var reportSHA256: String
    }
    public var schema: String
    public var caseID: String
    public var subjectID: String
    public var dataClass: DataClass
    public var sourceCitation: String
    public var assembly: String
    public var referenceSHA256: String
    public var annotationRelease: String
    public var annotationSHA256: String
    public var hlaNomenclatureRelease: String
    public var hlaSampleID: String
    public var hlaAlleles: [String]
    public var samples: [Sample]
    public var externalRun: ExternalRun
}

public struct VivoNeoantigenFinding: Codable, Sendable, Equatable {
    public enum Severity: String, Codable, Sendable { case blocking, warning }
    public var code: String
    public var severity: Severity
    public var message: String
    public var nextAction: String
    public var owner: String
}

public struct VivoNeoantigenCandidate: Codable, Sendable, Equatable {
    public var id: String
    public var sourceLine: Int
    public var gene: String
    public var transcript: String
    public var hlaAllele: String
    public var mutantPeptide: String
    public var wildtypePeptide: String?
    /// No coordinate conversion: pVACseq Start/Stop are zero-based, half-open.
    public var chromosome: String
    public var start: Int
    public var stop: Int
    public var predictedMedianBindingNM: Double?
    public var tumorDNAFraction: Double?
    public var tumorRNAFraction: Double?
    public var normalDNAFraction: Double?
    public var geneExpression: Double?
    /// Includes unrecognized upstream columns without interpreting them as decisions.
    public var sourceFields: [String: String]
    public var evidenceGaps: [String]
}

public struct VivoNeoantigenReport: Codable, Sendable, Equatable {
    public enum State: String, Codable, Sendable { case blocked, noCandidates, researchReview }
    public var schema: String
    public var implementationSHA256: String
    public var caseFingerprint: String
    public var sourceSHA256: String
    public var manifest: VivoNeoantigenCase
    public var state: State
    public var findings: [VivoNeoantigenFinding]
    public var candidates: [VivoNeoantigenCandidate]
    public var limitations: [String]
}

public struct VivoNeoantigenReview: Codable, Sendable, Equatable {
    public struct Decision: Codable, Sendable, Equatable {
        public enum Disposition: String, Codable, Sendable { case retainForResearch, exclude, deferReview }
        public var candidateID: String
        public var disposition: Disposition
        public var rationale: String
    }
    public var schema: String
    public var reportSHA256: String
    public var reviewerID: String
    public var reviewedAt: Date
    public var decisions: [Decision]
}

public enum VivoNeoantigenWorkbench {
    public static let maximumTSVBytes = 16 * 1_024 * 1_024
    public static let maximumRows = 10_000
    public static let adapterFormat = "pvacseq.classI.all_epitopes.v1"
    public static let limitations = [
        "Research-use reference-data workbench; not for patient treatment or manufacturing.",
        "External-run metadata and sample identities are declarations, not independently verified execution or genetic matching.",
        "Raw DNA/RNA reads, variant calls, expression quantification and HLA typing are not recomputed here.",
        "Imported binding predictions do not demonstrate antigen presentation, immune recognition, safety or clinical benefit.",
        "No candidate ranking, clinical approval, dosing or vaccine manufacture is performed.",
        "An empty result does not prove absence of tumor neoantigens; examine upstream exclusions and failures.",
        "Hashes establish byte integrity, not scientific correctness, authorship or authorization."
    ]
    static func matches(_ text: String, _ pattern: String) -> Bool {
        text.range(of: pattern, options: .regularExpression) != nil
    }
    static func validDigest(_ value: String) -> Bool { matches(value, "^[a-f0-9]{64}$") }
    static func token(_ value: String) -> Bool { matches(value, "^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$") }
    static func missing(_ value: String) -> Bool { value.isEmpty || value == "NA" }
    static func digest(_ bytes: Data) throws -> String { try VivoCanonicalJSON.fingerprint(bytes).hex }
    public static func fingerprint(_ report: VivoNeoantigenReport) throws -> String {
        try digest(VivoCanonicalJSON.encode(report))
    }

    /// A metadata check, deliberately not described as sequencing QC or sample authentication.
    public static func preflight(_ manifest: VivoNeoantigenCase) -> [VivoNeoantigenFinding] {
        var findings: [VivoNeoantigenFinding] = []
        func block(_ code: String, _ message: String, _ action: String) {
            findings.append(.init(code: code, severity: .blocking, message: message, nextAction: action, owner: "data provider"))
        }
        if manifest.schema != "numivivo.org/neoantigen-case/v1" {
            block("schema", "Unsupported case schema.", "Export a neoantigen-case/v1 manifest.")
        }
        if !token(manifest.caseID) || !token(manifest.subjectID) {
            block("identifiers", "Case and subject identifiers must be bounded pseudonymous tokens.", "Use reference-data identifiers, not names or medical record numbers.")
        }
        if manifest.sourceCitation.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty || manifest.sourceCitation.utf8.count > 4096 {
            block("source", "A bounded source citation or synthetic-fixture declaration is required.", "Identify the public source or state that the data are synthetic.")
        }
        if !["GRCh37", "GRCh38"].contains(manifest.assembly) || !validDigest(manifest.referenceSHA256) || !validDigest(manifest.annotationSHA256) || !token(manifest.annotationRelease) || !token(manifest.hlaNomenclatureRelease) {
            block("reference", "Reference assembly, annotation and HLA release identifiers must be explicit and pinned.", "Supply GRCh37/GRCh38 and exact resource versions and SHA-256 digests; do not substitute genome builds.")
        }
        if manifest.samples.count != 3 || Set(manifest.samples.map(\.role)).count != 3 {
            block("sample_roles", "Exactly one tumor-DNA, matched-normal-DNA and tumor-RNA record is required by this profile.", "Ask the data provider for the missing sample records; remove ambiguous duplicates.")
        }
        if Set(manifest.samples.map(\.sampleID)).count != manifest.samples.count {
            block("sample_ids", "Sample identifiers are duplicated.", "Assign a distinct identifier to each source sample.")
        }
        for sample in manifest.samples.prefix(100) {
            if !token(sample.sampleID) || sample.subjectID != manifest.subjectID || !validDigest(sample.sourceSHA256) {
                block("sample_identity", "A sample identifier, declared subject or source digest is inconsistent.", "Reconcile the declared sample-to-subject mapping with the data provider.")
            }
            if sample.assembly != manifest.assembly || sample.referenceSHA256 != manifest.referenceSHA256 {
                block("reference_mismatch", "Sample and case reference identities differ.", "Regenerate a compatible export using the same reference; do not relabel coordinates.")
            }
        }
        if !manifest.samples.contains(where: { $0.role == .matchedNormalDNA && $0.sampleID == manifest.hlaSampleID }) {
            block("hla_sample", "This profile requires HLA typing assigned to the matched-normal sample.", "Provide an explicitly mapped matched-normal HLA typing result.")
        }
        if manifest.hlaAlleles.isEmpty || manifest.hlaAlleles.count > 6 || Set(manifest.hlaAlleles).count != manifest.hlaAlleles.count || !manifest.hlaAlleles.allSatisfy({ matches($0, "^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}$") }) {
            block("hla_profile", "This adapter accepts unique, two-field HLA-A/B/C allele names only.", "Supply canonical class-I alleles; class-II, suffix and higher-resolution conversions require a separately qualified adapter.")
        }
        let run = manifest.externalRun
        if run.tool != "pvacseq" || run.format != adapterFormat || !matches(run.version, "^[0-9]+\\.[0-9]+\\.[0-9]+([+-][A-Za-z0-9.-]+)?$") || !validDigest(run.executionImageSHA256) || !validDigest(run.annotatedVCFSHA256) || !validDigest(run.reportSHA256) {
            block("run_provenance", "Exact external tool, format, version, image, annotated-VCF and report identities are required.", "Export the pVACseq run provenance; 'latest' and unpinned resources are not accepted.")
        }
        if run.arguments.isEmpty || run.arguments.count > 256 || !run.arguments.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 4096 && !$0.contains("\u{0}") }) || run.modelVersions.isEmpty || run.modelVersions.count > 64 || !run.modelVersions.allSatisfy({ token($0.key) && token($0.value) && $0.value.lowercased() != "latest" }) {
            block("model_provenance", "Recorded arguments and pinned prediction-model versions are required.", "Capture the executed argument vector and exact predictor/resource versions; no commands are executed by this importer.")
        }
        if run.sampleInputs.count != manifest.samples.count || !manifest.samples.allSatisfy({ run.sampleInputs[$0.sampleID] == $0.sourceSHA256 }) {
            block("run_inputs", "The external run's declared input mapping does not match the case samples.", "Reconcile the run manifest with the case before importing results.")
        }
        findings.append(.init(code: "declared_metadata", severity: .warning,
            message: "Source digests and identity mappings are declared metadata; only the supplied TSV bytes are hashed here.",
            nextAction: "Review upstream QC, genetic identity, HLA typing and execution records independently.", owner: "research reviewer"))
        return findings
    }

    public static func analyze(manifest: VivoNeoantigenCase, tsv: Data, implementationSHA256: String) throws -> VivoNeoantigenReport {
        guard validDigest(implementationSHA256), tsv.count <= maximumTSVBytes else {
            throw VivoNeoantigenError.invalid("A pinned implementation digest and TSV at most 16 MiB are required.")
        }
        let source = try digest(tsv)
        var findings = preflight(manifest)
        if source != manifest.externalRun.reportSHA256 {
            findings.append(.init(code: "report_digest", severity: .blocking, message: "The TSV does not match the report digest declared in the case.", nextAction: "Obtain the exact report paired with this case; do not transplant another case's report.", owner: "data provider"))
        }
        let caseFingerprint = try digest(VivoCanonicalJSON.encode(manifest))
        let blocked = findings.contains { $0.severity == .blocking }
        let candidates = blocked ? [] : try parse(tsv, manifest: manifest, caseFingerprint: caseFingerprint)
        return .init(schema: "numivivo.org/neoantigen-report/v1", implementationSHA256: implementationSHA256,
            caseFingerprint: caseFingerprint, sourceSHA256: source, manifest: manifest,
            state: blocked ? .blocked : (candidates.isEmpty ? .noCandidates : .researchReview),
            findings: findings, candidates: candidates, limitations: limitations)
    }

    private static func parse(_ data: Data, manifest: VivoNeoantigenCase, caseFingerprint: String) throws -> [VivoNeoantigenCandidate] {
        guard var text = String(data: data, encoding: .utf8), !text.contains("\u{0}") else {
            throw VivoNeoantigenError.invalid("The report must be UTF-8 TSV without NUL characters.")
        }
        if text.hasPrefix("\u{FEFF}") { text.removeFirst() }
        text = text.replacingOccurrences(of: "\r\n", with: "\n")
        guard !text.contains("\r") else { throw VivoNeoantigenError.invalid("Bare carriage returns are not supported.") }
        var lines = text.components(separatedBy: "\n")
        if lines.last == "" { lines.removeLast() }
        guard !lines.isEmpty, lines.count <= maximumRows + 1 else { throw VivoNeoantigenError.invalid("A header and at most 10,000 rows are supported.") }
        let header = lines[0].components(separatedBy: "\t")
        let required = ["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "Gene Name", "HLA Allele", "Peptide Length", "MT Epitope Seq", "WT Epitope Seq"]
        guard header.count <= 256, Set(header).count == header.count, header.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 256 }), required.allSatisfy(header.contains) else {
            throw VivoNeoantigenError.invalid("Expected unaggregated pVACseq all_epitopes columns with unique, bounded headers. Aggregated/pVACbind reports are not interchangeable.")
        }
        var output: [VivoNeoantigenCandidate] = [], seen = Set<String>()
        for (offset, line) in lines.dropFirst().enumerated() {
            let number = offset + 2
            let cells = line.components(separatedBy: "\t")
            guard cells.count == header.count, cells.allSatisfy({ $0.utf8.count <= 4096 }) else { throw VivoNeoantigenError.invalid("TSV row \(number): wrong field count or an oversized field.") }
            let fields = Dictionary(uniqueKeysWithValues: zip(header, cells))
            func field(_ name: String) -> String { fields[name] ?? "NA" }
            func numberValue(_ name: String, upper: Double? = nil, integral: Bool = false, positive: Bool = false) throws -> Double? {
                let raw = field(name)
                if missing(raw) { return nil }
                guard let value = Double(raw), value.isFinite, value >= 0, (!positive || value > 0), upper.map({ value <= $0 }) ?? true, (!integral || value.rounded() == value) else {
                    throw VivoNeoantigenError.invalid("TSV row \(number): invalid numeric field '\(name)'; missing values must remain NA, not zero.")
                }
                return value
            }
            let peptide = field("MT Epitope Seq"), wildtype = field("WT Epitope Seq"), hla = field("HLA Allele")
            guard let length = Int(field("Peptide Length")), (8...15).contains(length), peptide.count == length,
                  matches(peptide, "^[ACDEFGHIKLMNPQRSTVWY]+$"),
                  missing(wildtype) || (wildtype.count == length && matches(wildtype, "^[ACDEFGHIKLMNPQRSTVWY]+$")),
                  manifest.hlaAlleles.contains(hla), token(field("Chromosome")), token(field("Transcript")), token(field("Gene Name")),
                  matches(field("Reference"), "^[ACGTN-]+$"), matches(field("Variant"), "^[ACGTN-]+$"),
                  let start = Int(field("Start")), let stop = Int(field("Stop")), start >= 0, stop >= start, stop <= 1_000_000_000 else {
                throw VivoNeoantigenError.invalid("TSV row \(number): invalid class-I peptide, case HLA allele, variant identity or zero-based coordinates.")
            }
            let binding = try numberValue("Median MT IC50 Score", positive: true)
            let dna = try numberValue("Tumor DNA VAF", upper: 1)
            let rna = try numberValue("Tumor RNA VAF", upper: 1)
            let normal = try numberValue("Normal VAF", upper: 1)
            let expression = try numberValue("Gene Expression")
            let rnaDepth = try numberValue("Tumor RNA Depth", integral: true)
            let dnaDepth = try numberValue("Tumor DNA Depth", integral: true)
            let normalDepth = try numberValue("Normal Depth", integral: true)
            _ = try numberValue("Transcript Expression")
            // Check recognized per-model and aggregate IC50/percentile columns,
            // retaining their exact names and values; never combine score classes.
            for name in header where name.contains("IC50 Score") && !name.contains("Method") {
                _ = try numberValue(name, positive: true)
            }
            for name in header where name.contains("Percentile") && !name.contains("Method") {
                _ = try numberValue(name, upper: 100)
            }
            for (depth, fraction) in [(dnaDepth, dna), (rnaDepth, rna), (normalDepth, normal)] {
                if depth == 0, let fraction, fraction > 0 { throw VivoNeoantigenError.invalid("TSV row \(number): nonzero VAF with zero read depth.") }
            }
            let rowHash = try digest(VivoCanonicalJSON.encode(fields))
            guard seen.insert(rowHash).inserted else { throw VivoNeoantigenError.invalid("TSV row \(number): duplicate complete row; resolve it upstream.") }
            let id = try digest(Data((caseFingerprint + ":" + String(number) + ":" + rowHash).utf8))
            var gaps = ["Antigen presentation and immune recognition have not been experimentally established by this import."]
            if binding == nil { gaps.append("Median binding prediction is missing; no replacement score was inferred.") }
            if dna == nil || dnaDepth == nil || dnaDepth == 0 { gaps.append("Tumor DNA support is missing or has zero coverage.") }
            if dna == 0 { gaps.append("Reported tumor DNA variant fraction is zero.") }
            if expression == nil { gaps.append("Gene-expression evidence is missing.") }
            if expression == 0 { gaps.append("Reported gene expression is zero.") }
            if rna == nil || rnaDepth == nil || rnaDepth == 0 { gaps.append("Mutant RNA support is missing or has zero coverage; gene expression alone is insufficient.") }
            if rna == 0 { gaps.append("Reported tumor RNA variant fraction is zero.") }
            if normal == nil || normalDepth == nil || normalDepth == 0 { gaps.append("Matched-normal evidence is missing or has zero coverage.") }
            if let normal, normal > 0 { gaps.append("Alternate reads are reported in the matched normal; review somatic origin and technical error.") }
            if missing(wildtype) { gaps.append("No wildtype comparison is available.") }
            if wildtype == peptide { gaps.append("The reported peptide is unchanged from wildtype; this row does not establish a mutation-specific peptide.") }
            output.append(.init(id: id, sourceLine: number, gene: field("Gene Name"), transcript: field("Transcript"),
                hlaAllele: hla, mutantPeptide: peptide, wildtypePeptide: missing(wildtype) ? nil : wildtype,
                chromosome: field("Chromosome"), start: start, stop: stop, predictedMedianBindingNM: binding,
                tumorDNAFraction: dna, tumorRNAFraction: rna, normalDNAFraction: normal, geneExpression: expression,
                sourceFields: fields, evidenceGaps: gaps))
        }
        return output
    }

    /// Retaining an item means retaining it for research review, not approving a vaccine.
    /// Reviewer IDs here are self-declared; authentication/signature enforcement is separate.
    public static func validate(_ review: VivoNeoantigenReview, against report: VivoNeoantigenReport) throws {
        guard review.schema == "numivivo.org/neoantigen-review/v1", review.reportSHA256 == (try fingerprint(report)),
              report.state != .blocked, token(review.reviewerID), review.reviewedAt.timeIntervalSince1970.isFinite,
              !review.decisions.isEmpty, review.decisions.count <= report.candidates.count else {
            throw VivoNeoantigenError.invalid("Review must identify a nonblocked exact report, a reviewer and at least one bounded decision.")
        }
        let ids = Set(report.candidates.map(\.id))
        var seen = Set<String>()
        for decision in review.decisions {
            guard ids.contains(decision.candidateID), seen.insert(decision.candidateID).inserted,
                  !decision.rationale.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
                  decision.rationale.utf8.count <= 4096, !decision.rationale.contains("\u{0}") else {
                throw VivoNeoantigenError.invalid("Each decision must reference a unique candidate in this report and include a bounded rationale.")
            }
        }
    }
}
