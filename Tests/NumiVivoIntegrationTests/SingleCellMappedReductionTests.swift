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
