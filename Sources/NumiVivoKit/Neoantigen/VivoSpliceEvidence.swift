import Foundation

public struct VivoSpliceManifest: Codable, Sendable, Equatable {
    public var schema: String
    public var binding: VivoGenomicCaseBinding
    public var provenanceMode: String
    public var sourceCitation: String
    public var sourceRevision: String
    public var reportStage: String
    public var reportSHA256: String
    public var regtoolsSHA256: String
    public var transcriptsSHA256: String
    public var executionReceiptSHA256: String?
}
public struct VivoSpliceCandidate: Codable, Sendable, Equatable {
    public var id: String
    public var sourceLine: Int
    public var variant: VivoAtlasVariantInput
    public var gene: String
    public var transcript: String
    public var peptide: String
    public var hlaAllele: String
    public var junctionID: String
    public var junctionStart0: Int
    public var junctionStop0: Int
    public var strand: String
    public var junctionReadCount: Int
    public var peptidePosition1: Int
    public var suppliedWildtypeContainsPeptide: Bool
    public var medianBindingNM: Double?
    public var evidenceState: String
    public var evidenceGaps: [String]
    public var sourceFields: [String: String]
    public var regtoolsFields: [String: String]
}
public struct VivoSpliceReport: Codable, Sendable, Equatable {
    public var schema: String
    public var manifest: VivoSpliceManifest
    public var manifestSHA256: String
    public var candidates: [VivoSpliceCandidate]
    public var limitations: [String]
}

/// Shared bounded TSV reader. It never skips malformed rows or fills missing
/// measurements with zero. Source fields remain available for expert review.
public enum VivoGenomicTSV {
    public static func rows(_ data: Data, required: [String], maximumRows: Int = 10_000) throws -> [[String: String]] {
        try VivoAtlasEvidence.require(data.count <= 16 * 1024 * 1024, "TSV exceeds 16 MiB.")
        guard var text = String(data: data, encoding: .utf8), !text.contains("\u{0}") else { throw VivoGenomicEvidenceError.invalid("Expected UTF-8 TSV without NUL characters.") }
        if text.hasPrefix("\u{FEFF}") { text.removeFirst() }
        text = text.replacingOccurrences(of: "\r\n", with: "\n")
        try VivoAtlasEvidence.require(!text.contains("\r"), "Bare carriage returns are unsupported.")
        var lines = text.components(separatedBy: "\n")
        if lines.last == "" { lines.removeLast() }
        try VivoAtlasEvidence.require(!lines.isEmpty && lines.count <= maximumRows + 1, "TSV row limit exceeded.")
        let header = lines[0].components(separatedBy: "\t")
        try VivoAtlasEvidence.require(header.count <= 256 && Set(header).count == header.count && header.allSatisfy { !$0.isEmpty && $0.utf8.count <= 256 } && required.allSatisfy(header.contains), "Missing, duplicate or oversized TSV headers.")
        return try lines.dropFirst().enumerated().map { offset, line in
            let values = line.components(separatedBy: "\t")
            try VivoAtlasEvidence.require(values.count == header.count && values.allSatisfy { $0.utf8.count <= 4096 }, "Invalid field count or size at TSV line \(offset + 2).")
            return Dictionary(uniqueKeysWithValues: zip(header, values))
        }
    }
    static func number(_ fields: [String: String], _ name: String, upper: Double = 1e12, integer: Bool = false) throws -> Double? {
        guard let text = fields[name], text != "NA", !text.isEmpty else { return nil }
        guard let value = Double(text), value.isFinite, value >= 0, value <= upper, !integer || value.rounded() == value else {
            throw VivoGenomicEvidenceError.invalid("Invalid numeric field: \(name).")
        }
        return value
    }
    static func count(_ fields: [String: String], _ name: String) throws -> Int {
        guard let value = try number(fields, name, upper: 1_000_000_000, integer: true) else { throw VivoGenomicEvidenceError.invalid("Required integer is missing: \(name).") }
        return Int(value)
    }
}

public enum VivoSpliceEvidence {
    public static let limitations = [
        "RNA junction counts are imported alignment-derived observations, not an independent sequencing-QC or genetic-identity check.",
        "RegTools links a genomic interval to a junction; this does not establish allele-specific causality or phasing.",
        "The peptide is checked against the supplied pVACsplice ALT protein; translation, transcript reconstruction and the upstream prediction models are not rerun by this importer.",
        "Absence from the paired supplied WT protein is not absence from the full normal proteome or all normal isoforms.",
        "Matched-normal RNA junction coverage and tumor specificity are not established by this input profile.",
        "No Atlas prediction substitutes for RNA support. Presentation, immune recognition, safety and clinical benefit remain unestablished."
    ]
    static func proteins(_ data: Data) throws -> [String: String] {
        try VivoAtlasEvidence.require(data.count <= 16 * 1024 * 1024, "Transcript FASTA exceeds 16 MiB.")
        guard let text = String(data: data, encoding: .utf8), !text.contains("\u{0}") else { throw VivoGenomicEvidenceError.invalid("Invalid protein FASTA encoding.") }
        var result: [String: String] = [:], current: String?
        for rawLine in text.components(separatedBy: "\n") {
            let line = rawLine.hasSuffix("\r") ? String(rawLine.dropLast()) : rawLine
            if line.isEmpty { continue }
            if line.hasPrefix(">") {
                let key = String(line.dropFirst())
                try VivoAtlasEvidence.require(!key.isEmpty && key.utf8.count <= 1024 && result[key] == nil && result.count < 20_000, "Duplicate, invalid or excessive FASTA header.")
                current = key; result[key] = ""
            } else {
                guard let key = current else { throw VivoGenomicEvidenceError.invalid("FASTA sequence precedes its header.") }
                try VivoAtlasEvidence.require(VivoAtlasEvidence.matches(line, "^[ACDEFGHIKLMNPQRSTVWYX*]+$") && result[key]!.count + line.count <= 100_000, "Invalid or oversized protein sequence.")
                result[key]! += line
            }
        }
        try VivoAtlasEvidence.require(result.values.allSatisfy { !$0.isEmpty && !$0.dropLast().contains("*") }, "Empty protein FASTA record or internal stop codon.")
        return result
    }
    public static func analyze(manifest: VivoSpliceManifest, manifestData: Data, reportTSV: Data,
                               regtoolsTSV: Data, transcriptsFASTA: Data) throws -> VivoSpliceReport {
        try VivoAtlasEvidence.validate(manifest.binding)
        try VivoAtlasEvidence.require(manifestData.count <= 128 * 1024 && (try VivoGenomicDocuments.decode(VivoSpliceManifest.self, from: manifestData)) == manifest, "Splice manifest bytes differ from their declaration.")
        try VivoAtlasEvidence.require(manifest.schema == "numivivo.org/splice-manifest/v1" && ["allEpitopes", "filtered"].contains(manifest.reportStage), "Unsupported splice schema or report stage.")
        try VivoAtlasEvidence.require(["synthetic", "upstreamReference", "externalRecorded"].contains(manifest.provenanceMode) && (manifest.provenanceMode == "synthetic") == (manifest.binding.dataClass == "synthetic"), "Splice provenance class mismatch.")
        try VivoAtlasEvidence.require(!manifest.sourceCitation.isEmpty && manifest.sourceCitation.utf8.count <= 4096 && !manifest.sourceRevision.isEmpty && manifest.sourceRevision.utf8.count <= 256, "Source citation and revision are required.")
        if manifest.provenanceMode == "externalRecorded" {
            try VivoAtlasEvidence.require(VivoAtlasEvidence.isDigest(manifest.executionReceiptSHA256 ?? ""), "External execution requires a recorded receipt digest.")
        }
        try VivoAtlasEvidence.require(try VivoAtlasEvidence.digest(reportTSV) == manifest.reportSHA256 && VivoAtlasEvidence.digest(regtoolsTSV) == manifest.regtoolsSHA256 && VivoAtlasEvidence.digest(transcriptsFASTA) == manifest.transcriptsSHA256, "Splice input bytes do not match their declared hashes.")
        let junctions = try VivoGenomicTSV.rows(regtoolsTSV, required: ["chrom", "start", "end", "name", "score", "strand", "anchor", "transcripts", "variant_info"])
        var byJunction: [String: [String: String]] = [:]
        for junction in junctions {
            let name = junction["name"]!
            try VivoAtlasEvidence.require(VivoAtlasEvidence.validToken(name) && byJunction[name] == nil && ["+", "-"].contains(junction["strand"]!), "Invalid, unstranded or duplicate RegTools junction.")
            let start = try VivoGenomicTSV.count(junction, "start"), stop = try VivoGenomicTSV.count(junction, "end")
            try VivoAtlasEvidence.require(start < stop && VivoAtlasEvidence.chromosome(junction["chrom"]!) != nil, "Invalid RegTools junction coordinates.")
            _ = try VivoGenomicTSV.count(junction, "score")
            byJunction[name] = junction
        }
        let fields = try VivoGenomicTSV.rows(reportTSV, required: ["Chromosome", "Start", "Stop", "Reference", "Variant", "Junction", "Junction Start", "Junction Stop", "Junction Score", "Junction Anchor", "Transcript", "Gene Name", "Protein Position", "HLA Allele", "Peptide Length", "Epitope Seq", "WT Protein Length", "ALT Protein Length", "Index"])
        let fasta = try proteins(transcriptsFASTA), manifestSHA = try VivoAtlasEvidence.digest(manifestData)
        var candidates: [VivoSpliceCandidate] = [], seen = Set<String>()
        for (i, row) in fields.enumerated() {
            let rowDigest = try VivoAtlasEvidence.digest(VivoCanonicalJSON.encode(row))
            try VivoAtlasEvidence.require(seen.insert(rowDigest).inserted, "Duplicate complete pVACsplice row.")
            let id = try VivoAtlasEvidence.digest(Data((manifestSHA + ":" + String(i + 2) + ":" + rowDigest).utf8))
            guard let chromosome = VivoAtlasEvidence.chromosome(row["Chromosome"]!), let junction = byJunction[row["Junction"]!] else { throw VivoGenomicEvidenceError.invalid("Report chromosome or RNA junction cannot be resolved.") }
            let start = try VivoGenomicTSV.count(row, "Start"), stop = try VivoGenomicTSV.count(row, "Stop")
            let jStart = try VivoGenomicTSV.count(row, "Junction Start"), jStop = try VivoGenomicTSV.count(row, "Junction Stop")
            let count = try VivoGenomicTSV.count(row, "Junction Score")
            try VivoAtlasEvidence.require(start <= stop && jStart < jStop && chromosome == VivoAtlasEvidence.chromosome(junction["chrom"]!) && jStart == (try VivoGenomicTSV.count(junction, "start")) && jStop == (try VivoGenomicTSV.count(junction, "end")) && count == (try VivoGenomicTSV.count(junction, "score")) && row["Junction Anchor"] == junction["anchor"], "Junction coordinates, counts or anchor differ from RegTools evidence.")
            let interval = "\(chromosome):\(start)-\(stop)"
            try VivoAtlasEvidence.require(junction["variant_info"]!.components(separatedBy: ",").contains(interval), "RegTools does not link the reported variant interval to this junction.")
            let transcript = row["Transcript"]!, baseTranscript = transcript.components(separatedBy: ".")[0]
            try VivoAtlasEvidence.require(VivoAtlasEvidence.matches(transcript, "^ENST[0-9]+\\.[0-9]+$") && junction["transcripts"]!.components(separatedBy: ",").contains(baseTranscript), "Versioned report transcript does not match a RegTools base transcript ID.")
            let peptide = row["Epitope Seq"]!, length = try VivoGenomicTSV.count(row, "Peptide Length")
            try VivoAtlasEvidence.require((8...15).contains(length) && peptide.count == length && VivoAtlasEvidence.matches(peptide, "^[ACDEFGHIKLMNPQRSTVWY]+$") && manifest.binding.hlaAlleles.contains(row["HLA Allele"]!), "Invalid class-I peptide or case HLA allele.")
            let indexParts = row["Index"]!.components(separatedBy: ".")
            try VivoAtlasEvidence.require(indexParts.count == 7 && Int(indexParts[0]).map { $0 >= 0 } == true &&
                indexParts[1] == row["Gene Name"]! && indexParts[2] == baseTranscript && indexParts[3] == row["Junction"]! &&
                indexParts[4] == interval && indexParts[5] == row["Junction Anchor"]! && VivoAtlasEvidence.validToken(indexParts[6]),
                "Protein FASTA index does not identify the same gene, transcript, junction and variant as the report row.")
            guard let alt = fasta["ALT." + row["Index"]!], let wt = fasta["WT." + row["Index"]!] else { throw VivoGenomicEvidenceError.invalid("Paired ALT/WT transcript proteins are missing for the report index.") }
            let position = try VivoGenomicTSV.count(row, "Protein Position")
            try VivoAtlasEvidence.require(position > 0 && position - 1 <= alt.count && length <= alt.count - (position - 1), "Peptide position is outside the supplied ALT protein.")
            let aa = Array(alt.utf8), expected = String(decoding: aa[(position - 1)..<(position - 1 + length)], as: UTF8.self)
            try VivoAtlasEvidence.require(expected == peptide && alt.count == (try VivoGenomicTSV.count(row, "ALT Protein Length")) && wt.count == (try VivoGenomicTSV.count(row, "WT Protein Length")), "Peptide sequence/position or declared protein length does not match the supplied FASTA.")
            try VivoAtlasEvidence.require(VivoAtlasEvidence.matches(row["Reference"]!, "^[ACGTN-]+$") && VivoAtlasEvidence.matches(row["Variant"]!, "^[ACGTN-]+$"), "Invalid genomic alleles.")
            let binding = try VivoGenomicTSV.number(row, "Median IC50 Score")
            if let binding { try VivoAtlasEvidence.require(binding > 0, "IC50 prediction must be positive or missing.") }
            for key in row.keys where key.contains("IC50 Score") && !key.contains("Method") {
                if let value = try VivoGenomicTSV.number(row, key) { try VivoAtlasEvidence.require(value > 0, "Invalid IC50 prediction.") }
            }
            for key in row.keys where key.contains("Percentile") && !key.contains("Method") { _ = try VivoGenomicTSV.number(row, key, upper: 100) }
            for key in ["Tumor DNA VAF", "Tumor RNA VAF", "Normal VAF"] { _ = try VivoGenomicTSV.number(row, key, upper: 1) }
            for key in ["Tumor DNA Depth", "Tumor RNA Depth", "Normal Depth"] { _ = try VivoGenomicTSV.number(row, key, upper: 1e9, integer: true) }
            for key in ["Gene Expression", "Transcript Expression"] { _ = try VivoGenomicTSV.number(row, key) }
            for (depthKey, vafKey) in [("Tumor DNA Depth", "Tumor DNA VAF"), ("Tumor RNA Depth", "Tumor RNA VAF"), ("Normal Depth", "Normal VAF")] {
                let depth = try VivoGenomicTSV.number(row, depthKey, upper: 1e9, integer: true)
                let fraction = try VivoGenomicTSV.number(row, vafKey, upper: 1)
                if depth == 0, let fraction { try VivoAtlasEvidence.require(fraction == 0, "Nonzero variant fraction with zero read depth.") }
            }
            let wtContains = wt.contains(peptide)
            var gaps = ["Matched-normal RNA junction evidence is not supplied.", "RegTools transcript matching is by unversioned base ID; transcript-version equivalence is not proven.", "Presentation and immune recognition have not been experimentally established."]
            if count == 0 { gaps.append("No supporting RNA junction reads are reported.") }
            if wt.contains("X") { gaps.append("The supplied WT protein contains ambiguous residues; literal sequence absence does not establish biological specificity.") }
            if wtContains { gaps.append("This peptide is present in the paired supplied WT protein.") }
            if binding == nil { gaps.append("Median binding prediction is missing.") }
            if manifest.reportStage == "filtered" { gaps.append("This is a prefiltered upstream subset; excluded candidates are not represented.") }
            if let normal = try VivoGenomicTSV.number(row, "Normal VAF", upper: 1), normal > 0 { gaps.append("Alternate reads are present in matched-normal DNA; somatic origin requires review.") }
            candidates.append(.init(id: id, sourceLine: i + 2,
                variant: .init(candidateID: id, chromosome: chromosome, start0: start, stop0: stop, reference: row["Reference"]!, alternate: row["Variant"]!),
                gene: row["Gene Name"]!, transcript: transcript, peptide: peptide, hlaAllele: row["HLA Allele"]!, junctionID: row["Junction"]!,
                junctionStart0: jStart, junctionStop0: jStop, strand: junction["strand"]!, junctionReadCount: count, peptidePosition1: position,
                suppliedWildtypeContainsPeptide: wtContains, medianBindingNM: binding,
                evidenceState: count == 0 ? "noRNAJunctionSupport" : (wtContains ? "presentInSuppliedWildtype" : "RNAJunctionAndALTSequenceConsistent"),
                evidenceGaps: gaps, sourceFields: row, regtoolsFields: junction))
        }
        return .init(schema: "numivivo.org/splice-report/v1", manifest: manifest, manifestSHA256: manifestSHA, candidates: candidates, limitations: limitations)
    }
}
