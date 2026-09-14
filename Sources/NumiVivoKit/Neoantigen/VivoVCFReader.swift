import Foundation

/// The allele classes accepted by the bounded VCF reader.  `snv` is the only
/// class that can be projected into the current AlphaGenome Atlas request
/// contract; the other classes remain present so unsupported biology is not
/// silently discarded.
public enum VivoVCFAlleleClass: String, Codable, Sendable, Equatable {
    case snv
    case indel
    case symbolic
    case breakend
}

public struct VivoVCFInfoField: Codable, Sendable, Equatable {
    public var key: String
    /// `nil` represents a VCF flag; an empty value is rejected by the parser.
    public var value: String?
}

public struct VivoVCFCall: Codable, Sendable, Equatable {
    public var sampleID: String
    /// Values remain textual, including `.` and genotype allele indexes.  No
    /// ploidy, phasing or somatic interpretation is inferred by this reader.
    public var fields: [String: String]
}

public struct VivoVCFHeader: Codable, Sendable, Equatable {
    public var fileFormat: String
    public var metadataLines: [String]
    public var columns: [String]
    public var sampleIDs: [String]
}

public struct VivoVCFVariant: Codable, Sendable, Equatable {
    /// Stable identity derived from the exact source digest, source line and
    /// alternate index.  It is suitable for an explicit downstream binding.
    public var candidateID: String
    public var sourceLine: Int
    public var sourceRecordSHA256: String
    public var alternateIndex: Int
    /// Primary contigs are represented canonically as `chr1` ... `chrX`; other
    /// valid VCF contigs are retained verbatim for explicit exclusion later.
    public var chromosome: String
    /// VCF coordinates are one-based and inclusive at this boundary.
    public var position1: Int
    public var identifier: String?
    public var reference: String
    public var alternate: String
    public var alleleClass: VivoVCFAlleleClass
    public var quality: Double?
    /// `nil` is the VCF `.` value; otherwise the original semicolon-delimited
    /// FILTER field is represented as individual tokens.
    public var filters: [String]?
    public var info: [VivoVCFInfoField]
    public var formatKeys: [String]
    public var calls: [VivoVCFCall]
}

public struct VivoVCFImport: Codable, Sendable, Equatable {
    public var schema: String
    public var sourceSHA256: String
    public var sourceBytes: Int
    public var assembly: String
    public var referenceSHA256: String
    public var header: VivoVCFHeader
    public var variants: [VivoVCFVariant]
}

public struct VivoVCFAtlasExclusion: Codable, Sendable, Equatable {
    public var candidateID: String
    public var sourceLine: Int
    public var reason: String
}

public struct VivoVCFAtlasProjection: Codable, Sendable, Equatable {
    public var schema: String
    public var sourceSHA256: String
    public var assembly: String
    public var referenceSHA256: String
    public var eligibleVariants: [VivoAtlasVariantInput]
    public var exclusions: [VivoVCFAtlasExclusion]
}

/// Strict, source-preserving VCF ingestion for public-reference research.
///
/// This is deliberately a parser and provenance boundary, not a variant caller,
/// liftover tool, genotype interpreter or phenotype predictor.  The source
/// bytes are hashed before any line-ending normalization.  Only an explicit
/// GRCh38 projection of primary-contig SNVs can enter the current Atlas request
/// planner; all other valid records remain represented as exclusions.
public enum VivoVCFReader {
    public static let maximumDocumentBytes = 64 * 1_024 * 1_024
    public static let maximumRecords = 250_000
    public static let maximumHeaderLines = 10_000
    public static let maximumLineBytes = 1 * 1_024 * 1_024

    private static let requiredColumns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    public static func parse(data: Data, assembly: String, referenceSHA256: String) throws -> VivoVCFImport {
        try VivoAtlasEvidence.require(data.count <= maximumDocumentBytes, "VCF exceeds the 64 MiB bounded import limit.")
        try VivoAtlasEvidence.require(["GRCh37", "GRCh38"].contains(assembly), "VCF import requires an explicit GRCh37 or GRCh38 assembly.")
        try VivoAtlasEvidence.require(VivoAtlasEvidence.isDigest(referenceSHA256), "VCF import requires the SHA-256 of the exact reference FASTA.")
        if data.count >= 2, data[0] == 0x1f, data[1] == 0x8b {
            throw VivoGenomicEvidenceError.invalid("Gzip-compressed VCF is not accepted at this boundary; decompress it explicitly and retain the compressed source hash separately.")
        }
        try VivoAtlasEvidence.require(!data.contains(where: { $0 == 0 }), "VCF contains a NUL byte.")
        guard let originalText = String(data: data, encoding: .utf8), !originalText.isEmpty else {
            throw VivoGenomicEvidenceError.invalid("VCF must be nonempty UTF-8 text.")
        }
        var text = originalText
        if text.contains("\r") {
            text = text.replacingOccurrences(of: "\r\n", with: "\n")
            try VivoAtlasEvidence.require(!text.contains("\r"), "Bare carriage returns are not supported in VCF.")
        }
        var lines = text.split(separator: "\n", omittingEmptySubsequences: false).map(String.init)
        if lines.last == "" { lines.removeLast() }
        try VivoAtlasEvidence.require(!lines.isEmpty, "VCF has no header or records.")

        var metadata: [String] = []
        var fileFormat: String?
        var columns: [String] = []
        var columnLine = 0
        var records: [VivoVCFVariant] = []
        let sourceSHA256 = try VivoAtlasEvidence.digest(data)

        for (offset, line) in lines.enumerated() {
            let lineNumber = offset + 1
            try VivoAtlasEvidence.require(line.utf8.count <= maximumLineBytes, "VCF line \(lineNumber) exceeds the 1 MiB line limit.")
            if columns.isEmpty {
                if line.hasPrefix("##") {
                    try VivoAtlasEvidence.require(!line.contains("\t") && line.utf8.count >= 3, "Malformed VCF metadata line \(lineNumber).")
                    metadata.append(line)
                    if line.hasPrefix("##fileformat=") {
                        try VivoAtlasEvidence.require(fileFormat == nil, "VCF declares fileformat more than once.")
                        let value = String(line.dropFirst("##fileformat=".count))
                        try VivoAtlasEvidence.require(matches(value, "^VCFv4\\.[0-9]+$"), "Unsupported VCF fileformat: \(value).")
                        fileFormat = value
                    }
                    try VivoAtlasEvidence.require(metadata.count <= maximumHeaderLines, "VCF header exceeds the bounded metadata-line limit.")
                    continue
                }
                let parsedColumns = line.split(separator: "\t", omittingEmptySubsequences: false).map(String.init)
                try validateColumns(parsedColumns, lineNumber: lineNumber)
                columns = parsedColumns
                columnLine = lineNumber
                continue
            }

            try VivoAtlasEvidence.require(!line.isEmpty && !line.hasPrefix("#"), "Unexpected header or empty line after VCF column header at line \(lineNumber).")
            try VivoAtlasEvidence.require(records.count < maximumRecords, "VCF exceeds the \(maximumRecords)-record limit.")
            let fields = line.split(separator: "\t", omittingEmptySubsequences: false).map(String.init)
            try VivoAtlasEvidence.require(fields.count == columns.count, "VCF line \(lineNumber) has \(fields.count) columns; expected \(columns.count).")
            records += try parseRecord(fields, columns: columns, sourceLine: lineNumber, sourceRecordSHA256: try VivoAtlasEvidence.digest(Data(line.utf8)), sourceSHA256: sourceSHA256)
        }

        guard let fileFormat else { throw VivoGenomicEvidenceError.invalid("VCF is missing a ##fileformat declaration.") }
        try VivoAtlasEvidence.require(columnLine > 0, "VCF is missing its #CHROM column header.")
        let sampleIDs = columns.count > 9 ? Array(columns.dropFirst(9)) : []
        let header = VivoVCFHeader(fileFormat: fileFormat, metadataLines: metadata, columns: columns, sampleIDs: sampleIDs)
        return VivoVCFImport(schema: "numivivo.org/vcf-import/v1", sourceSHA256: sourceSHA256,
            sourceBytes: data.count, assembly: assembly, referenceSHA256: referenceSHA256,
            header: header, variants: records)
    }

    /// Converts only explicit GRCh38 primary-contig SNVs into the current
    /// Atlas input identity.  Unsupported alleles and contigs are retained as
    /// bounded, machine-readable exclusions rather than being dropped.
    public static func atlasProjection(_ document: VivoVCFImport) throws -> VivoVCFAtlasProjection {
        try VivoAtlasEvidence.require(document.schema == "numivivo.org/vcf-import/v1" && ["GRCh37", "GRCh38"].contains(document.assembly) && VivoAtlasEvidence.isDigest(document.sourceSHA256) && VivoAtlasEvidence.isDigest(document.referenceSHA256), "Invalid VCF import identity.")
        var eligible: [VivoAtlasVariantInput] = [], exclusions: [VivoVCFAtlasExclusion] = []
        for variant in document.variants {
            if document.assembly != "GRCh38" {
                exclusions.append(.init(candidateID: variant.candidateID, sourceLine: variant.sourceLine, reason: "unsupportedAssembly: GRCh38 is required; no automatic liftover"))
                continue
            }
            guard let chromosome = VivoAtlasEvidence.chromosome(variant.chromosome) else {
                exclusions.append(.init(candidateID: variant.candidateID, sourceLine: variant.sourceLine, reason: "unsupportedContig: primary chr1-22/X/Y required"))
                continue
            }
            guard variant.alleleClass == .snv else {
                exclusions.append(.init(candidateID: variant.candidateID, sourceLine: variant.sourceLine, reason: "unsupportedAllele: exact single-base substitution required"))
                continue
            }
            eligible.append(.init(candidateID: variant.candidateID, chromosome: chromosome,
                start0: variant.position1 - 1, stop0: variant.position1,
                reference: variant.reference, alternate: variant.alternate))
        }
        return VivoVCFAtlasProjection(schema: "numivivo.org/vcf-atlas-projection/v1",
            sourceSHA256: document.sourceSHA256, assembly: document.assembly,
            referenceSHA256: document.referenceSHA256, eligibleVariants: eligible, exclusions: exclusions)
    }

    private static func validateColumns(_ columns: [String], lineNumber: Int) throws {
        try VivoAtlasEvidence.require(columns.count >= requiredColumns.count && Array(columns.prefix(requiredColumns.count)) == requiredColumns, "VCF line \(lineNumber) must begin with #CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO.")
        if columns.count > 8 { try VivoAtlasEvidence.require(columns[8] == "FORMAT", "VCF line \(lineNumber) requires FORMAT as column nine when sample columns are present.") }
        let samples = columns.count > 9 ? Array(columns.dropFirst(9)) : []
        try VivoAtlasEvidence.require(samples.allSatisfy { validToken($0) }, "VCF sample identifiers must be bounded tokens.")
        try VivoAtlasEvidence.require(Set(samples).count == samples.count, "VCF sample identifiers must be unique.")
    }

    private static func parseRecord(_ fields: [String], columns: [String], sourceLine: Int,
                                    sourceRecordSHA256: String, sourceSHA256: String) throws -> [VivoVCFVariant] {
        let rawChromosome = fields[0]
        try VivoAtlasEvidence.require(validToken(rawChromosome), "VCF line \(sourceLine) has an invalid contig identifier.")
        guard let position1 = Int(fields[1]), position1 > 0, position1 <= 300_000_000 else {
            throw VivoGenomicEvidenceError.invalid("VCF line \(sourceLine) has an invalid one-based position.")
        }
        let identifier: String? = fields[2] == "." ? nil : try checkedValue(fields[2], name: "ID", line: sourceLine, maximum: 4096)
        let reference = fields[3]
        try VivoAtlasEvidence.require(matches(reference, "^[ACGTN]+$"), "VCF line \(sourceLine) has an unsupported REF allele; uppercase A/C/G/T/N is required.")
        let alternatives = fields[4].split(separator: ",", omittingEmptySubsequences: false).map(String.init)
        try VivoAtlasEvidence.require(!alternatives.isEmpty && alternatives.allSatisfy { !$0.isEmpty && $0 != "." } && Set(alternatives).count == alternatives.count, "VCF line \(sourceLine) has empty, missing or duplicate ALT alleles.")
        let quality: Double? = fields[5] == "." ? nil : try checkedNumber(fields[5], name: "QUAL", line: sourceLine, upper: nil)
        let filters = try parseFilters(fields[6], line: sourceLine)
        let info = try parseInfo(fields[7], line: sourceLine)
        let formatKeys: [String]
        let calls: [VivoVCFCall]
        if columns.count > 8 {
            (formatKeys, calls) = try parseCalls(format: fields[8], sampleValues: Array(fields.dropFirst(9)), sampleIDs: Array(columns.dropFirst(9)), line: sourceLine)
        } else {
            formatKeys = []; calls = []
        }
        var result: [VivoVCFVariant] = []
        for (index, alternate) in alternatives.enumerated() {
            let alternateIndex = index + 1
            let alleleClass = try classify(reference: reference, alternate: alternate, line: sourceLine)
            let identity = "\(sourceSHA256):\(sourceLine):\(alternateIndex):\(alternate)"
            let candidateID = try VivoAtlasEvidence.digest(Data(identity.utf8))
            result.append(.init(candidateID: candidateID, sourceLine: sourceLine,
                sourceRecordSHA256: sourceRecordSHA256, alternateIndex: alternateIndex,
                chromosome: VivoAtlasEvidence.chromosome(rawChromosome) ?? rawChromosome,
                position1: position1, identifier: identifier, reference: reference,
                alternate: alternate, alleleClass: alleleClass, quality: quality,
                filters: filters, info: info, formatKeys: formatKeys, calls: calls))
        }
        return result
    }

    private static func classify(reference: String, alternate: String, line: Int) throws -> VivoVCFAlleleClass {
        if matches(alternate, "^[ACGTN]+$") {
            return reference.count == 1 && alternate.count == 1 ? .snv : .indel
        }
        if alternate == "*" || matches(alternate, "^<[A-Za-z][A-Za-z0-9_.-]*>$") { return .symbolic }
        if alternate.contains("[") || alternate.contains("]") { return .breakend }
        throw VivoGenomicEvidenceError.invalid("VCF line \(line) has an unsupported ALT allele encoding.")
    }

    private static func parseFilters(_ raw: String, line: Int) throws -> [String]? {
        guard raw != "." else { return nil }
        let values = raw.split(separator: ";", omittingEmptySubsequences: false).map(String.init)
        try VivoAtlasEvidence.require(!values.isEmpty && values.allSatisfy { validToken($0) } && Set(values).count == values.count, "VCF line \(line) has invalid or duplicate FILTER tokens.")
        try VivoAtlasEvidence.require(!(values.contains("PASS") && values.count > 1), "VCF line \(line) combines PASS with another FILTER token.")
        return values
    }

    private static func parseInfo(_ raw: String, line: Int) throws -> [VivoVCFInfoField] {
        guard raw != "." else { return [] }
        var output: [VivoVCFInfoField] = [], seen = Set<String>()
        for part in raw.split(separator: ";", omittingEmptySubsequences: false).map(String.init) {
            let pieces = part.split(separator: "=", maxSplits: 1, omittingEmptySubsequences: false).map(String.init)
            let key = pieces[0]
            try VivoAtlasEvidence.require(matches(key, "^[A-Za-z][A-Za-z0-9_.]*$") && seen.insert(key).inserted, "VCF line \(line) has an invalid or duplicate INFO key.")
            let value = pieces.count == 1 ? nil : try checkedValue(pieces[1], name: "INFO.\(key)", line: line, maximum: 4096)
            output.append(.init(key: key, value: value))
        }
        return output
    }

    private static func parseCalls(format: String, sampleValues: [String], sampleIDs: [String], line: Int) throws -> ([String], [VivoVCFCall]) {
        if format == "." {
            try VivoAtlasEvidence.require(sampleValues.allSatisfy { $0 == "." }, "VCF line \(line) uses FORMAT='.' but contains sample values.")
            return ([], sampleIDs.map { .init(sampleID: $0, fields: [:]) })
        }
        let keys = format.split(separator: ":", omittingEmptySubsequences: false).map(String.init)
        try VivoAtlasEvidence.require(!keys.isEmpty && keys.allSatisfy { validToken($0) } && Set(keys).count == keys.count, "VCF line \(line) has invalid or duplicate FORMAT keys.")
        var calls: [VivoVCFCall] = []
        for (sampleID, raw) in zip(sampleIDs, sampleValues) {
            let values = raw.split(separator: ":", omittingEmptySubsequences: false).map(String.init)
            try VivoAtlasEvidence.require(values.count == keys.count, "VCF line \(line) sample \(sampleID) has a FORMAT/value count mismatch.")
            var fields: [String: String] = [:]
            for (key, value) in zip(keys, values) { fields[key] = try checkedValue(value, name: "FORMAT.\(key)", line: line, maximum: 4096) }
            calls.append(.init(sampleID: sampleID, fields: fields))
        }
        return (keys, calls)
    }

    private static func checkedValue(_ value: String, name: String, line: Int, maximum: Int) throws -> String {
        try VivoAtlasEvidence.require(value.utf8.count <= maximum && !value.contains("\u{0}") && !value.contains("\t") && !value.contains("\n") && !value.contains("\r"), "VCF line \(line) has an invalid or oversized \(name) value.")
        return value
    }

    private static func checkedNumber(_ value: String, name: String, line: Int, upper: Double?) throws -> Double {
        guard let number = Double(value), number.isFinite, number >= 0, upper.map({ number <= $0 }) ?? true else {
            throw VivoGenomicEvidenceError.invalid("VCF line \(line) has an invalid \(name) value.")
        }
        return number
    }

    private static func matches(_ text: String, _ expression: String) -> Bool {
        text.range(of: expression, options: .regularExpression) != nil
    }

    private static func validToken(_ text: String) -> Bool {
        matches(text, "^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$")
    }
}
