import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellMappedReductionTests {
    static func bytes(_ entries: [(UInt32,UInt32,Double)]) -> Data {
        var data=Data()
        for (row,column,value) in entries {
            var r=row.littleEndian,c=column.littleEndian,v=value.bitPattern.littleEndian
            withUnsafeBytes(of: &r) { data.append(contentsOf: $0) }
            withUnsafeBytes(of: &c) { data.append(contentsOf: $0) }
            withUnsafeBytes(of: &v) { data.append(contentsOf: $0) }
        }
        return data
    }
    @Test func mappedProductsMatchExplicitMatrixAndAdjoint() throws {
        let url=FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: url) }
        try Self.bytes([(0,0,2),(0,2,3),(1,1,4),(3,0,5)]).write(to: url)
        let map=try VivoMappedReductionEntries(url: url,expectedEntries: 4,rows: 4,columns: 3)
        #expect(try map.project([2,3,4],shift: 1,rows: 4) == [15,11,-1,9])
        #expect(try map.transpose([1,2,3,4],initial: [0.5,-1,2]) == [22.5,7,5])
        let x=[1.2,-0.3,2.1],y=[-0.2,0.7,3.2,-1.1]
        let left=VivoSingleCellReduction.dot(try map.project(x,shift: 0,rows: 4),y)
        let right=VivoSingleCellReduction.dot(x,try map.transpose(y,initial: [0,0,0]))
        #expect(abs(left-right)<1e-12)
        #expect(map.visits==16 && map.byteCount==64)
    }
    @Test func repeatedProductsCrossMappingWindows() throws {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: url) }
        let n = VivoWindowedCountRecords.windowBytes / 16 + 3
        let writer = try VivoCountRecordWriter(url)
        var sums = [Double](repeating: 0, count: 6)
        for i in 0..<n {
            let row = i % 3, column = (i / 3) % 2, value = Double(i % 5 + 1)
            sums[row * 2 + column] += value
            try writer.append(row: row, feature: column, bits: value.bitPattern)
        }
        _ = try writer.finish()
        let map = try VivoMappedReductionEntries(url: url, expectedEntries: n, rows: 3, columns: 2)
        var expected = [Double]()
        for row in 0..<3 { expected.append(-2.0 + sums[row * 2] * 3.0 + sums[row * 2 + 1] * 7.0) }
        let first = 1.0 + sums[0] * 2.0 + sums[2] * 3.0 + sums[4] * 5.0
        let second = 4.0 + sums[1] * 2.0 + sums[3] * 3.0 + sums[5] * 5.0
        for _ in 0..<3 {
            #expect(try map.project([3,7], shift: 2, rows: 3) == expected)
            #expect(try map.transpose([2,3,5], initial: [1,4]) == [first, second])
        }
        #expect(map.visits == n * 6)
        #expect(throws: (any Error).self) { try map.project([1], shift: 0, rows: 3) }
        #expect(throws: (any Error).self) { try map.transpose([1], initial: [0,0]) }
    }
    @Test func malformedCachesAndBudgetsAreRejected() throws {
        for data in [Data(repeating: 0,count: 15),Self.bytes([(4,0,1)]),Self.bytes([(0,3,1)]),Self.bytes([(0,0,.nan)]),Self.bytes([(0,0,0)])] {
            let url=FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
            defer { try? FileManager.default.removeItem(at: url) }
            try data.write(to: url)
            #expect(throws: (any Error).self) { try VivoMappedReductionEntries(url: url,expectedEntries: 1,rows: 4,columns: 3) }
        }
        var options=VivoH5ADReductionOptions();options.maximumCacheBytes=15
        #expect(throws: (any Error).self) { try options.validate() }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoH5ADReductionOptions.self,from: Data("{\"maximumCacheByte\":16}".utf8)) }
    }
}
