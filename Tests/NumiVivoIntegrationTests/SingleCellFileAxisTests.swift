import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellFileAxisTests {
    static let implementation = SingleCellCountStreamTests.implementation
    func workspace(_ body: (URL) throws -> Void) throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("file-axis-test-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        try body(root)
    }
    func header(_ data: VivoSingleCellDataset, totals: Bool = true) throws -> VivoCellAxisHeader {
        let m = data.metadata
        return .init(metadata: .init(id: m.id, evidence: m.evidence, sourceDescription: m.sourceDescription,
            countUnit: m.countUnit, samples: m.samples, features: m.features, cells: []), cellCount: m.cells.count,
            nonzeros: data.matrix.counts.count, hasRowTotals: totals, annotations: Set(m.cells.compactMap(\.group)).sorted(), sourceDeclaration: "Numerical control")
    }
    func rows(_ data: VivoSingleCellDataset, totals: Bool = true) throws -> Data {
        let qc = try VivoSingleCellAnalysis.quality(data)
        return try data.cells.indices.reduce(into: Data()) { output, i in
            output.append(try VivoCanonicalJSON.encode(VivoCellAxisInputRow(cell: data.cells[i], nonzeros: qc[i].detectedFeatures, totalCounts: totals ? qc[i].totalCounts : nil)))
            output.append(10)
        }
    }
    func axis(_ root: URL, _ data: VivoSingleCellDataset = SingleCellCountStreamTests.fixture(), totals: Bool = true) throws -> URL {
        let input = root.appendingPathComponent("cells.jsonl"), output = root.appendingPathComponent("axis")
        try rows(data, totals: totals).write(to: input)
        let file = try FileHandle(forReadingFrom: input); defer { try? file.close() }
        _ = try VivoFileCellAxis.publish(header: header(data, totals: totals), input: file, implementation: Self.implementation, to: output)
        return output
    }
    func publish(_ root: URL, _ axis: URL, _ data: Data) throws -> URL {
        let input = root.appendingPathComponent("counts.bin"), output = root.appendingPathComponent("counts")
        try data.write(to: input); let file = try FileHandle(forReadingFrom: input); defer { try? file.close() }
        _ = try VivoFileCountStream.publish(axis: axis, input: file, implementation: Self.implementation, to: output)
        return output
    }
    @Test func diskQualityAndMembershipMatchIndependentResidentAnalysis() throws {
        try workspace { root in
            let fixture = SingleCellCountStreamTests.fixture(), source = try axis(root)
            let directory = try publish(root, source, SingleCellCountStreamTests.bytes())
            let output = try VivoFileCountSnapshot.open(directory, implementation: Self.implementation)
            let expected = try VivoSingleCellAnalysis.pseudobulk(fixture), qc = try VivoSingleCellAnalysis.quality(fixture)
            #expect(output.report.matrix == expected.matrix)
            #expect(output.report.featureIDs == expected.featureIDs)
            var members = [[Int]](repeating: [], count: output.report.groups.count)
            for i in fixture.cells.indices {
                let actual = try output.quality(i); #expect(actual.quality == qc[i]); members[actual.group].append(i)
                #expect(try output.axis.row(i).cell == fixture.cells[i])
            }
            for i in expected.groups.indices {
                #expect(members[i] == expected.groups[i].sourceCellIndices)
                #expect(output.report.groups[i].sourceCellCount == expected.groups[i].sourceCellIndices.count)
            }
            let file = try FileHandle(forReadingFrom: root.appendingPathComponent("counts.bin")); defer { try? file.close() }
            #expect(try VivoFileCountStream.verify(directory, input: file, implementation: Self.implementation) == output.receipt)
            #expect(throws: (any Error).self) { try output.quality(-1) }
            #expect(throws: (any Error).self) { try output.quality(fixture.cells.count) }
        }
    }
    @Test func fragmentedCanonicalStreamAndAbsentIndependentTotals() throws {
        try workspace { root in
            let source = try axis(root, totals: false), axis = try VivoFileCellAxis.open(source, implementation: Self.implementation)
            let bytes = SingleCellCountStreamTests.bytes(), expected = try VivoSingleCellAnalysis.pseudobulk(SingleCellCountStreamTests.fixture())
            for width in [1, 7, 16, 17, 48] {
                var position = 0
                let result = try VivoFileCountStream.evaluate(axis: axis, qualityURL: root.appendingPathComponent("quality-\(width).bin")) {
                    let stop = min(position + width, bytes.count); defer { position = stop }; return Data(bytes[position..<stop])
                }
                #expect(result.0.matrix == expected.matrix)
                #expect(result.1 == (try VivoCanonicalJSON.fingerprint(bytes)))
            }
        }
    }
    @Test func immutableSnapshotSurvivesSourceMutationAndRehashedDuplicateIsRejected() throws {
        try workspace { root in
            let source = try axis(root), opened = try VivoFileCellAxis.open(source, implementation: Self.implementation)
            let path = source.appendingPathComponent("strings.bin")
            let records = try Data(contentsOf: source.appendingPathComponent("rows.bin"))
            let start = Int(records.vivoLE(UInt64.self, at: 2 * VivoFileCellAxis.recordBytes))
            var bytes = try Data(contentsOf: path); bytes.replaceSubrange(start..<(start + 2), with: Data("c0".utf8)); try bytes.write(to: path)
            #expect(try opened.row(2).cell.barcode == "c2")
            let old = opened.receipt
            let changed = VivoCellAxisReceipt(schemaVersion: old.schemaVersion, encoding: old.encoding, header: old.header,
                rows: old.rows, strings: try VivoCanonicalJSON.fingerprint(bytes), stringBytes: old.stringBytes,
                inputJSONL: old.inputJSONL, implementation: old.implementation)
            try VivoCanonicalJSON.encode(changed).write(to: source.appendingPathComponent("receipt.json"))
            #expect(throws: (any Error).self) { try VivoFileCellAxis.open(source, implementation: Self.implementation) }
        }
    }
    @Test func unicodeEquivalentBarcodesAreDuplicateWithinSample() throws {
        try workspace { root in
            let original = SingleCellCountStreamTests.fixture()
            let h = VivoCellAxisHeader(metadata: try header(original).metadata, cellCount: 2, nonzeros: 0, hasRowTotals: false, annotations: [], sourceDeclaration: "Unicode identity control")
            var input = Data()
            for barcode in ["\u{e9}", "e\u{301}"] {
                input.append(try VivoCanonicalJSON.encode(VivoCellAxisInputRow(cell: .init(barcode: barcode, sampleID: "s"), nonzeros: 0, totalCounts: nil))); input.append(10)
            }
            let path = root.appendingPathComponent("unicode.jsonl"); try input.write(to: path)
            let file = try FileHandle(forReadingFrom: path); defer { try? file.close() }
            let output = root.appendingPathComponent("rejected")
            #expect(throws: (any Error).self) { try VivoFileCellAxis.publish(header: h, input: file, implementation: Self.implementation, to: output) }
            #expect(!FileManager.default.fileExists(atPath: output.path))
            #expect(!(try FileManager.default.contentsOfDirectory(atPath: root.path)).contains { $0.hasPrefix(".numivivo-axis-") })
        }
    }
    @Test func malformedCountsRejectWithoutPublishing() throws {
        try workspace { root in
            let source = try axis(root, totals: false)
            let bad = [Data(SingleCellCountStreamTests.bytes().dropLast()), SingleCellCountStreamTests.bytes() + Data([0]),
                SingleCellCountStreamTests.bytes([(0,0,4),(0,0,2),(2,1,7)]),
                SingleCellCountStreamTests.bytes([(0,1,2),(0,0,4),(2,1,7)]),
                SingleCellCountStreamTests.bytes([(0,0,4),(3,1,2),(2,1,7)]),
                SingleCellCountStreamTests.bytes([(0,0,4),(0,2,2),(2,1,7)]),
                SingleCellCountStreamTests.bytes([(0,0,4),(0,1,0),(2,1,7)]),
                SingleCellCountStreamTests.bytes([(0,0,UInt64.max),(0,1,1),(2,1,7)]),
                SingleCellCountStreamTests.bytes([(0,0,1),(0,1,UInt64.max-1),(2,1,7)])]
            for bytes in bad {
                #expect(throws: (any Error).self) { try publish(root, source, bytes) }
                #expect(!FileManager.default.fileExists(atPath: root.appendingPathComponent("counts").path))
            }
            #expect(!(try FileManager.default.contentsOfDirectory(atPath: root.path)).contains { $0.hasPrefix(".numivivo-file-count-") })
        }
    }
    @Test func rehashedChangedAggregateRequiresNumericalReconstruction() throws {
        try workspace { root in
            let source = try axis(root), directory = try publish(root, source, SingleCellCountStreamTests.bytes())
            let saved = try VivoFileCountSnapshot.open(directory, implementation: Self.implementation), old = saved.receipt
            let report = saved.report, matrix = report.matrix
            var counts = matrix.counts; counts[0] += 1
            let altered = VivoFileCountStreamReport(schemaVersion: report.schemaVersion, method: report.method, cellCount: report.cellCount,
                canonicalNonzeros: report.canonicalNonzeros, countUnit: report.countUnit, featureIDs: report.featureIDs, groups: report.groups,
                matrix: .init(cellCount: matrix.cellCount, featureCount: matrix.featureCount, rowOffsets: matrix.rowOffsets, featureIndices: matrix.featureIndices, counts: counts))
            let bytes = try VivoCanonicalJSON.encode(altered); try bytes.write(to: directory.appendingPathComponent("report.json"))
            let receipt = VivoFileCountStreamReceipt(schemaVersion: old.schemaVersion, encoding: old.encoding, qualityEncoding: old.qualityEncoding,
                streamBytes: old.streamBytes, stream: old.stream, axis: old.axis, quality: old.quality,
                report: try VivoCanonicalJSON.fingerprint(bytes), implementation: old.implementation)
            try VivoCanonicalJSON.encode(receipt).write(to: directory.appendingPathComponent("receipt.json"))
            let file = try FileHandle(forReadingFrom: root.appendingPathComponent("counts.bin")); defer { try? file.close() }
            #expect(throws: (any Error).self) { try VivoFileCountStream.verify(directory, input: file, implementation: Self.implementation) }
        }
    }
    @Test func maximumIntegerTotalIsNotUsedAsMissingValue() throws {
        try workspace { root in
            let original = SingleCellCountStreamTests.fixture()
            let data = VivoSingleCellDataset(id: "max", evidence: .synthetic, sourceDescription: "Exact UInt64 control", countUnit: .umiCount,
                samples: original.samples, features: [.init(id: "a", name: "a")], cells: [.init(barcode: "c", sampleID: "s")],
                matrix: .init(cellCount: 1, featureCount: 1, rowOffsets: [0,1], featureIndices: [0], counts: [UInt64.max]))
            let source = try axis(root, data), directory = try publish(root, source, SingleCellCountStreamTests.bytes([(0,0,UInt64.max)]))
            let opened = try VivoFileCountSnapshot.open(directory, implementation: Self.implementation)
            #expect(try opened.axis.row(0).totalCounts == UInt64.max)
            #expect(try opened.quality(0).quality.totalCounts == UInt64.max)
            #expect(try opened.quality(0).quality.mitochondrialFraction == nil)
        }
    }
    @Test func malformedJSONLRejectsUnknownFieldsMissingNewlineAndWrongTotals() throws {
        try workspace { root in
            let original = SingleCellCountStreamTests.fixture(), h = try header(original), good = try rows(original)
            var changed = try #require(JSONSerialization.jsonObject(with: good.split(separator: 10)[0]) as? [String: Any])
            changed["inferredLabel"] = "not-admitted"
            var unknown = try JSONSerialization.data(withJSONObject: changed); unknown.append(10)
            unknown.append(good.suffix(from: good.firstIndex(of: 10)! + 1))
            let bad = [Data(good.dropLast()), unknown, good + Data("\n".utf8)]
            for (i, bytes) in bad.enumerated() {
                let path = root.appendingPathComponent("bad-\(i).jsonl"); try bytes.write(to: path)
                let file = try FileHandle(forReadingFrom: path); defer { try? file.close() }
                #expect(throws: (any Error).self) { try VivoFileCellAxis.publish(header: h, input: file, implementation: Self.implementation, to: root.appendingPathComponent("rejected-\(i)")) }
            }
            let wrong = VivoCellAxisHeader(metadata: h.metadata, cellCount: h.cellCount, nonzeros: h.nonzeros + 1, hasRowTotals: h.hasRowTotals, annotations: h.annotations, sourceDeclaration: h.sourceDeclaration)
            let path = root.appendingPathComponent("valid.jsonl"); try good.write(to: path)
            let file = try FileHandle(forReadingFrom: path); defer { try? file.close() }
            #expect(throws: (any Error).self) { try VivoFileCellAxis.publish(header: wrong, input: file, implementation: Self.implementation, to: root.appendingPathComponent("wrong-total")) }
        }
    }

    @Test func technicalSamplesAnnotationsAndEmptyGroupsPreserveMembership() throws {
        try workspace { root in
            let samples = [
                VivoOmicsSample(id: "s1", biologicalReplicateID: "d", donorID: "d", condition: "PBS", batchID: "b1", organism: "human"),
                VivoOmicsSample(id: "s2", biologicalReplicateID: "d", donorID: "d", condition: "PBS", batchID: "b2", organism: "human"),
                VivoOmicsSample(id: "s3", biologicalReplicateID: "d", donorID: "d", condition: "IFNB", batchID: "b3", organism: "human")]
            let data = VivoSingleCellDataset(id: "groups", evidence: .synthetic, sourceDescription: "Technical replicate and annotation control", countUnit: .umiCount,
                samples: samples, features: SingleCellCountStreamTests.fixture().features,
                cells: [.init(barcode: "x", sampleID: "s1", group: "B"), .init(barcode: "x", sampleID: "s2", group: "B"),
                        .init(barcode: "z", sampleID: "s1"), .init(barcode: "x", sampleID: "s3", group: "B"), .init(barcode: "empty", sampleID: "s2")],
                matrix: .init(cellCount: 5, featureCount: 2, rowOffsets: [0,1,3,3,4,4], featureIndices: [0,0,1,1], counts: [3,2,4,9]))
            let source = try axis(root, data), directory = try publish(root, source, SingleCellCountStreamTests.bytes([(0,0,3),(1,0,2),(1,1,4),(3,1,9)]))
            let actual = try VivoFileCountSnapshot.open(directory, implementation: Self.implementation)
            let expected = try VivoSingleCellAnalysis.pseudobulk(data), quality = try VivoSingleCellAnalysis.quality(data)
            #expect(actual.report.matrix == expected.matrix)
            var members = [[Int]](repeating: [], count: expected.groups.count)
            for row in data.cells.indices { let value = try actual.quality(row); members[value.group].append(row); #expect(value.quality == quality[row]) }
            for i in expected.groups.indices {
                let a = actual.report.groups[i], b = expected.groups[i]
                #expect(a.biologicalReplicateID == b.biologicalReplicateID && a.donorID == b.donorID && a.organism == b.organism)
                #expect(a.condition == b.condition && a.cellGroup == b.cellGroup && a.sampleIDs == b.sampleIDs && a.batchIDs == b.batchIDs)
                #expect(members[i] == b.sourceCellIndices && a.sourceCellCount == b.sourceCellIndices.count)
            }
            var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(header(data))) as? [String: Any])
            object["silentSelection"] = true
            #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoCellAxisHeader.self, from: JSONSerialization.data(withJSONObject: object)) }
        }
    }

    @Test func dictionaryReferencesRetainOriginalUnicodeBytes() throws {
        try workspace { root in
            let composed = "\u{e9}", decomposed = "e\u{301}"
            let data = VivoSingleCellDataset(id: "unicode-references", evidence: .synthetic, sourceDescription: "Lossless reference spelling", countUnit: .umiCount,
                samples: [.init(id: composed, biologicalReplicateID: "d", condition: "PBS", batchID: "b", organism: "human")],
                features: [.init(id: "a", name: "a")], cells: [.init(barcode: decomposed, sampleID: decomposed, group: decomposed)],
                matrix: .init(cellCount: 1, featureCount: 1, rowOffsets: [0,1], featureIndices: [0], counts: [1]))
            let h = try header(data)
            let alternateDictionary = VivoCellAxisHeader(metadata: h.metadata, cellCount: 1, nonzeros: 1, hasRowTotals: true, annotations: [composed], sourceDeclaration: h.sourceDeclaration)
            let path = root.appendingPathComponent("unicode.jsonl"); try rows(data).write(to: path)
            let file = try FileHandle(forReadingFrom: path); defer { try? file.close() }
            let source = root.appendingPathComponent("axis")
            _ = try VivoFileCellAxis.publish(header: alternateDictionary, input: file, implementation: Self.implementation, to: source)
            let opened = try VivoFileCellAxis.open(source, implementation: Self.implementation), cell = try opened.row(0).cell
            #expect(Array(cell.barcode.utf8) == Array(decomposed.utf8))
            #expect(Array(cell.sampleID.utf8) == Array(decomposed.utf8))
            #expect(Array(cell.group!.utf8) == Array(decomposed.utf8))
        }
    }
}
