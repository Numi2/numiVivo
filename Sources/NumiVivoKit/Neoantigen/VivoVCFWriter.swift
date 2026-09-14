import Foundation

/// Canonical, provenance-checked VCF serialization for an imported document.
///
/// The writer preserves the parsed record semantics and multiallelic grouping,
/// but it does not promise byte identity with the source file.  It is an export
/// boundary, not a variant caller, liftover tool, genotype interpreter or
/// phenotype predictor.
public enum VivoVCFWriter {
    public static func encode(_ document: VivoVCFImport) throws -> Data {
        let records = try validate(document)
        var metadata = document.header.metadataLines
        if !metadata.contains(where: { $0.hasPrefix("##fileformat=") }) {
            metadata.insert("##fileformat=\(document.header.fileFormat)", at: 0)
        }
        var lines = metadata
        lines.append(document.header.columns.joined(separator: "\t"))
        lines.append(contentsOf: records.map(render(_:)))
        let data = Data((lines.joined(separator: "\n") + "\n").utf8)
        guard data.count <= VivoVCFReader.maximumDocumentBytes else {
            throw VivoGenomicEvidenceError.invalid("canonical VCF export exceeds the 64 MiB bound")
        }
        return data
    }

    public static func write(_ document: VivoVCFImport, to url: URL) throws {
        guard url.isFileURL else { throw VivoGenomicEvidenceError.invalid("VCF export requires a local file") }
        try encode(document).write(to: url, options: .withoutOverwriting)
    }

    private struct Record {
        let variants: [VivoVCFVariant]
        let hasFormat: Bool
    }

    private static func validate(_ document: VivoVCFImport) throws -> [Record] {
        try require(document.schema == "numivivo.org/vcf-import/v1", "unsupported VCF import schema")
        try require(VivoAtlasEvidence.isDigest(document.sourceSHA256), "VCF import has an invalid source digest")
        try require(["GRCh37", "GRCh38"].contains(document.assembly), "VCF import has an invalid assembly")
        try require(VivoAtlasEvidence.isDigest(document.referenceSHA256), "VCF import has an invalid reference digest")
        try require(document.sourceBytes >= 0 && document.sourceBytes <= VivoVCFReader.maximumDocumentBytes,
                    "VCF import source size is outside the bounded export domain")
        try require(document.header.fileFormat.range(of: "^VCFv4\\.[0-9]+$", options: .regularExpression) != nil,
                    "unsupported VCF fileformat")
        try require(document.header.columns.count >= 8,
                    "VCF header must contain the eight required columns")
        let required = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        try require(Array(document.header.columns.prefix(8)) == required,
                    "VCF header does not begin with the required columns")
        if document.header.columns.count > 8 {
            try require(document.header.columns[8] == "FORMAT", "VCF FORMAT column is missing")
        }
        let samples = document.header.columns.count > 9 ? Array(document.header.columns.dropFirst(9)) : []
        try require(samples == document.header.sampleIDs,
                    "VCF header sample identifiers do not match its columns")
        try require(samples.allSatisfy(validToken) && Set(samples).count == samples.count,
                    "VCF sample identifiers are invalid or duplicated")
        try require(document.variants.count <= VivoVCFReader.maximumRecords,
                    "VCF export exceeds the record bound")

        var seenCandidates = Set<String>()
        var groupIndices: [Int: Int] = [:]
        var groups: [Record] = []
        for variant in document.variants {
            try validate(variant, sourceSHA256: document.sourceSHA256, sampleIDs: samples,
                         hasFormat: document.header.columns.count > 8)
            try require(seenCandidates.insert(variant.candidateID).inserted,
                        "VCF candidate identifiers must be unique")
            if let index = groupIndices[variant.sourceLine] {
                let first = groups[index].variants[0]
                try require(first.sourceRecordSHA256 == variant.sourceRecordSHA256 &&
                            sameRecordMetadata(first, variant),
                            "VCF variants sharing a source line must reconstruct one record")
                groups[index] = Record(variants: groups[index].variants + [variant], hasFormat: groups[index].hasFormat)
            } else {
                groupIndices[variant.sourceLine] = groups.count
                groups.append(Record(variants: [variant], hasFormat: document.header.columns.count > 8))
            }
        }
        for group in groups {
            let ordered = group.variants.sorted { $0.alternateIndex < $1.alternateIndex }
            let expected = Array(1...ordered.count)
            try require(ordered.map(\.alternateIndex) == expected,
                        "VCF alternate indexes must be contiguous within each record")
        }
        try require(document.header.metadataLines.count <= VivoVCFReader.maximumHeaderLines,
                    "VCF metadata exceeds the header bound")
        var formatCount = 0
        for line in document.header.metadataLines {
            try require(line.hasPrefix("##") && line.utf8.count >= 3 &&
                        !line.contains("\0") && !line.contains("\t") &&
                        !line.contains("\n") && !line.contains("\r"),
                        "VCF metadata contains an invalid line")
            if line.hasPrefix("##fileformat=") {
                formatCount += 1
                try require(line == "##fileformat=\(document.header.fileFormat)",
                            "VCF metadata fileformat disagrees with the parsed header")
            }
        }
        try require(formatCount <= 1, "VCF metadata declares fileformat more than once")
        return groups
    }

    private static func validate(_ variant: VivoVCFVariant, sourceSHA256: String,
                                 sampleIDs: [String], hasFormat: Bool) throws {
        try require(VivoAtlasEvidence.isDigest(variant.candidateID) &&
                    VivoAtlasEvidence.isDigest(variant.sourceRecordSHA256) &&
                    variant.sourceLine > 0 && variant.alternateIndex > 0,
                    "VCF variant provenance identity is invalid")
        try require(validToken(variant.chromosome) && (1...300_000_000).contains(variant.position1),
                    "VCF variant coordinate is invalid")
        if let identifier = variant.identifier {
            try checked(identifier, maximum: 4096, separators: [])
        }
        try require(matches(variant.reference, "^[ACGTN]+$"), "VCF REF allele is invalid")
        try checked(variant.alternate, maximum: 1_048_576, separators: [",", ";"])
        try require(variant.alternate != ".", "VCF ALT allele cannot be missing")
        let identity = "\(sourceSHA256):\(variant.sourceLine):\(variant.alternateIndex):\(variant.alternate)"
        try require(variant.candidateID == (try VivoAtlasEvidence.digest(Data(identity.utf8))),
                    "VCF candidate identifier is not bound to source, line, alternate index and allele")
        guard let alleleClass = classify(variant.reference, variant.alternate) else {
            throw VivoGenomicEvidenceError.invalid("VCF ALT allele encoding is unsupported")
        }
        try require(alleleClass == variant.alleleClass,
                    "VCF allele class does not match REF/ALT")
        if let quality = variant.quality {
            try require(quality.isFinite && quality >= 0, "VCF QUAL must be finite and nonnegative")
        }
        if let filters = variant.filters {
            try require(!filters.isEmpty && filters.allSatisfy(validToken) &&
                        Set(filters).count == filters.count &&
                        !(filters.contains("PASS") && filters.count > 1),
                        "VCF FILTER values are invalid")
        }
        var infoKeys = Set<String>()
        for field in variant.info {
            try require(matches(field.key, "^[A-Za-z][A-Za-z0-9_.]*$") && infoKeys.insert(field.key).inserted,
                        "VCF INFO key is invalid or duplicated")
            if let value = field.value { try checked(value, maximum: 4096, separators: [";"]) }
        }
        if !hasFormat {
            try require(variant.formatKeys.isEmpty && variant.calls.isEmpty,
                        "VCF FORMAT data is present without a FORMAT header column")
            return
        }
        try require(variant.formatKeys.allSatisfy(validToken) &&
                    Set(variant.formatKeys).count == variant.formatKeys.count,
                    "VCF FORMAT keys are invalid or duplicated")
        try require(variant.calls.count == sampleIDs.count,
                    "VCF call count does not match the header samples")
        for (call, sampleID) in zip(variant.calls, sampleIDs) {
            try require(call.sampleID == sampleID && Set(call.fields.keys) == Set(variant.formatKeys),
                        "VCF call identity or FORMAT keys do not match")
            for key in variant.formatKeys {
                guard let value = call.fields[key] else { throw VivoGenomicEvidenceError.invalid("VCF call is missing a FORMAT value") }
                try checked(value, maximum: 4096, separators: [":"])
            }
        }
        if variant.formatKeys.isEmpty {
            try require(variant.calls.allSatisfy { $0.fields.isEmpty },
                        "VCF FORMAT='.' calls must not carry values")
        }
    }

    private static func sameRecordMetadata(_ lhs: VivoVCFVariant, _ rhs: VivoVCFVariant) -> Bool {
        lhs.chromosome == rhs.chromosome && lhs.position1 == rhs.position1 &&
        lhs.identifier == rhs.identifier && lhs.reference == rhs.reference &&
        lhs.quality == rhs.quality && lhs.filters == rhs.filters && lhs.info == rhs.info &&
        lhs.formatKeys == rhs.formatKeys && lhs.calls == rhs.calls
    }

    private static func render(_ record: Record) -> String {
        let variants = record.variants.sorted { $0.alternateIndex < $1.alternateIndex }
        let first = variants[0]
        let identifier = first.identifier ?? "."
        let alternate = variants.map(\.alternate).joined(separator: ",")
        let quality = first.quality.map { String($0) } ?? "."
        let filters = first.filters?.joined(separator: ";") ?? "."
        let info = first.info.isEmpty ? "." : first.info.map { field in
            field.value.map { "\(field.key)=\($0)" } ?? field.key
        }.joined(separator: ";")
        var fields = [first.chromosome, String(first.position1), identifier, first.reference,
                      alternate, quality, filters, info]
        if record.hasFormat {
            fields.append(first.formatKeys.isEmpty ? "." : first.formatKeys.joined(separator: ":"))
            for call in first.calls {
                fields.append(first.formatKeys.isEmpty ? "." : first.formatKeys.map { call.fields[$0]! }.joined(separator: ":"))
            }
        }
        return fields.joined(separator: "\t")
    }

    private static func classify(_ reference: String, _ alternate: String) -> VivoVCFAlleleClass? {
        if matches(alternate, "^[ACGTN]+$") {
            return reference.count == 1 && alternate.count == 1 ? .snv : .indel
        }
        if alternate == "*" || matches(alternate, "^<[A-Za-z][A-Za-z0-9_.-]*>$") { return .symbolic }
        if alternate.contains("[") || alternate.contains("]") { return .breakend }
        return nil
    }

    private static func require(_ condition: Bool, _ message: String) throws {
        guard condition else { throw VivoGenomicEvidenceError.invalid(message) }
    }

    private static func checked(_ value: String, maximum: Int, separators: [Character]) throws {
        try require(value.utf8.count <= maximum && !value.contains("\0") &&
                    !value.contains("\t") && !value.contains("\n") && !value.contains("\r") &&
                    separators.allSatisfy { !value.contains($0) },
                    "VCF field contains an invalid or oversized value")
    }

    private static func matches(_ text: String, _ expression: String) -> Bool {
        text.range(of: expression, options: .regularExpression) != nil
    }

    private static func validToken(_ text: String) -> Bool {
        matches(text, "^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$")
    }
}
