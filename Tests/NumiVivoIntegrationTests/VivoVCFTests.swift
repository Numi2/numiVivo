import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct VivoVCFTests {
    private let referenceSHA256 = String(repeating: "a", count: 64)

    private func fixture() -> Data {
        Data("""
        ##fileformat=VCFv4.3
        ##reference=file:///references/GRCh38.fa
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\tnormal
        1\t10\trs1\tA\tC,G\t50\tPASS\tDP=9;SOMATIC\tGT:DP\t0/1:6\t0/0:3
        chrX\t20\t.\tT\tTA\t10\tq10\tDP=4\tGT\t1/1\t0/1
        chrM\t30\t.\tC\t<DEL>\t.\t.\t.\tGT\t.\t./.
        """.utf8)
    }

    @Test func parsesMultiallelicCallsAndProjectsOnlyPrimarySNVs() throws {
        let data = fixture()
        let document = try VivoVCFReader.parse(data: data, assembly: "GRCh38", referenceSHA256: referenceSHA256)
        #expect(document.schema == "numivivo.org/vcf-import/v1")
        #expect(document.sourceBytes == data.count)
        #expect(document.sourceSHA256 == (try VivoAtlasEvidence.digest(data)))
        #expect(document.header.fileFormat == "VCFv4.3")
        #expect(document.header.sampleIDs == ["tumor", "normal"])
        #expect(document.variants.count == 4)
        #expect(document.variants[0].chromosome == "chr1")
        #expect(document.variants[0].alternate == "C")
        #expect(document.variants[0].alleleClass == .snv)
        #expect(document.variants[1].alternate == "G")
        #expect(document.variants[1].alternateIndex == 2)
        #expect(document.variants[2].alleleClass == .indel)
        #expect(document.variants[3].alleleClass == .symbolic)
        #expect(document.variants[0].filters == ["PASS"])
        #expect(document.variants[0].info.map(\.key) == ["DP", "SOMATIC"])
        #expect(document.variants[0].info.map(\.value) == ["9", nil])
        #expect(document.variants[0].calls[0].fields["GT"] == "0/1")
        #expect(document.variants[0].calls[1].fields["DP"] == "3")

        let projection = try VivoVCFReader.atlasProjection(document)
        #expect(projection.schema == "numivivo.org/vcf-atlas-projection/v1")
        #expect(projection.eligibleVariants.count == 2)
        #expect(projection.eligibleVariants.map(\.chromosome) == ["chr1", "chr1"])
        #expect(projection.eligibleVariants.map(\.alternate) == ["C", "G"])
        #expect(projection.exclusions.count == 2)
        #expect(projection.exclusions.contains { $0.reason.hasPrefix("unsupportedAllele:") })
    }

    @Test func assemblyAndUnsupportedAllelesRemainExplicit() throws {
        let document = try VivoVCFReader.parse(data: fixture(), assembly: "GRCh37", referenceSHA256: referenceSHA256)
        let projection = try VivoVCFReader.atlasProjection(document)
        #expect(projection.eligibleVariants.isEmpty)
        #expect(projection.exclusions.count == document.variants.count)
        #expect(projection.exclusions.allSatisfy { $0.reason.hasPrefix("unsupportedAssembly:") })
    }

    @Test func rejectsDuplicateInfoAndMismatchedSampleFields() throws {
        let duplicateInfo = Data("""
        ##fileformat=VCFv4.3
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        1\t1\t.\tA\tC\t.\t.\tDP=1;DP=2
        """.utf8)
        #expect(throws: (any Error).self) {
            try VivoVCFReader.parse(data: duplicateInfo, assembly: "GRCh38", referenceSHA256: referenceSHA256)
        }
        let mismatchedCall = Data("""
        ##fileformat=VCFv4.3
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
        1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1
        """.utf8)
        #expect(throws: (any Error).self) {
            try VivoVCFReader.parse(data: mismatchedCall, assembly: "GRCh38", referenceSHA256: referenceSHA256)
        }
    }

    @Test func sourceBytesAreBoundBeforeLineNormalization() throws {
        let lf = fixture()
        let crlf = Data(String(decoding: lf, as: UTF8.self).replacingOccurrences(of: "\n", with: "\r\n").utf8)
        let a = try VivoVCFReader.parse(data: lf, assembly: "GRCh38", referenceSHA256: referenceSHA256)
        let b = try VivoVCFReader.parse(data: crlf, assembly: "GRCh38", referenceSHA256: referenceSHA256)
        #expect(a.sourceSHA256 != b.sourceSHA256)
        #expect(a.variants.map(\.candidateID) != b.variants.map(\.candidateID))
        #expect(a.variants.map(\.alternate) == b.variants.map(\.alternate))
    }

    @Test func canonicalExportRoundTripsMultiallelicRecords() throws {
        let source = fixture()
        let document = try VivoVCFReader.parse(data: source, assembly: "GRCh38", referenceSHA256: referenceSHA256)
        let encoded = try VivoVCFWriter.encode(document)
        let roundTrip = try VivoVCFReader.parse(data: encoded, assembly: "GRCh38", referenceSHA256: referenceSHA256)
        #expect(roundTrip.header.fileFormat == document.header.fileFormat)
        #expect(roundTrip.header.metadataLines == document.header.metadataLines)
        #expect(roundTrip.header.columns == document.header.columns)
        #expect(roundTrip.variants.count == document.variants.count)
        for (left, right) in zip(document.variants, roundTrip.variants) {
            #expect(right.chromosome == left.chromosome)
            #expect(right.position1 == left.position1)
            #expect(right.identifier == left.identifier)
            #expect(right.reference == left.reference)
            #expect(right.alternate == left.alternate)
            #expect(right.alleleClass == left.alleleClass)
            #expect(right.quality == left.quality)
            #expect(right.filters == left.filters)
            #expect(right.info == left.info)
            #expect(right.formatKeys == left.formatKeys)
            #expect(right.calls == left.calls)
        }
        #expect(roundTrip.sourceSHA256 != document.sourceSHA256)
    }

    @Test func exportRejectsInconsistentProvenanceAndDelimiterInjection() throws {
        let document = try VivoVCFReader.parse(data: fixture(), assembly: "GRCh38", referenceSHA256: referenceSHA256)
        var inconsistent = document
        inconsistent.variants[1].sourceRecordSHA256 = String(repeating: "b", count: 64)
        #expect(throws: (any Error).self) { try VivoVCFWriter.encode(inconsistent) }

        var injected = document
        injected.variants[0].info[0].value = "unsafe;value"
        #expect(throws: (any Error).self) { try VivoVCFWriter.encode(injected) }
    }
}
