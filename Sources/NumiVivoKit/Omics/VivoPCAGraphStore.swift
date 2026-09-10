import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

public struct VivoPCAGraphStoreReport: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let method: String
    public let cells: Int
    public let dimensions: Int
    public let options: VivoSingleCellNeighborOptions
    public let neighborEntries: Int
    public let connectivityEntries: Int
    public let connectedComponents: Int
    public let isolatedCells: Int
    public let distancePairs: Int
    public let neighbors: VivoFingerprint
    public let bandwidths: VivoFingerprint
    public let offsets: VivoFingerprint
    public let edges: VivoFingerprint
    public let qualification: String
}

/// Binary files use the existing complete 16-byte little-endian record contract.
/// Temporary transpose storage and sequential merge avoid resident edge arrays.
/// Cell-scale degree/offset/component arrays remain resident; no cells-by-genes
/// or cells-by-cells matrix is created.
enum VivoPCAGraphStore {
    static let files = ["neighbors.bin", "bandwidths.bin", "offsets.bin", "edges.bin"]
    static func build(root: URL, rows n: Int, dimensions: Int, options: VivoSingleCellNeighborOptions,
                      distancePairs: Int, approximate: Bool, neighborsHash: VivoFingerprint) throws -> VivoPCAGraphStoreReport {
        try options.validate()
        let k = options.neighbors
        guard n >= k, n <= VivoPCAStorageLimits.maximumRows, (1...64).contains(dimensions) else { throw VivoOmicsError.limit("binary graph axes") }
        let records = try VivoWindowedCountRecords(root.appendingPathComponent("neighbors.bin"), entries: n*k)
        func row(_ i: Int) throws -> (indices: [Int], distances: [Double]) {
            var indices: [Int] = [], distances: [Double] = []
            for slot in 0..<k {
                let r = try records.record(i*k+slot), d = Double(bitPattern: r.bits)
                guard r.row == i, r.feature < n, d.isFinite, d >= 0 else { throw VivoOmicsError.invalid("binary neighbor record") }
                if slot == 0 {
                    guard r.feature == i, d == 0 else { throw VivoOmicsError.invalid("binary neighbor self") }
                } else {
                    guard r.feature != i, !indices.contains(r.feature) else { throw VivoOmicsError.invalid("binary neighbor duplicate") }
                    if slot > 1 {
                        guard d > distances[slot-1] || (d == distances[slot-1] && r.feature > indices[slot-1]) else {
                            throw VivoOmicsError.invalid("binary neighbor order")
                        }
                    }
                }
                indices.append(r.feature); distances.append(d)
            }
            return (indices, distances)
        }
        var sum = 0.0
        for i in 0..<n {
            try Task.checkCancellation()
            for d in try row(i).distances { sum += d }
        }
        guard sum.isFinite else { throw VivoOmicsError.invalid("binary graph distance sum") }
        let globalMean = sum/Double(n*k)
        let scratch = root.appendingPathComponent(".graph-transpose-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: scratch, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: scratch) }
        let directedURL = scratch.appendingPathComponent("directed.bin"), incomingURL = scratch.appendingPathComponent("incoming.bin")
        let directed = try VivoCountRecordWriter(directedURL), bandwidths = try VivoCountRecordWriter(root.appendingPathComponent("bandwidths.bin"))
        var degrees = [Int](repeating: 0, count: n)
        for i in 0..<n {
            try Task.checkCancellation()
            let r = try row(i), b = try VivoSingleCellNeighbors.bandwidth(row: r.distances, local: options.localConnectivity, globalMean: globalMean)
            for (column, value) in [b.rho, b.sigma, b.residual].enumerated() { try bandwidths.append(row: i, feature: column, bits: value.bitPattern) }
            for slot in 1..<k {
                let delta = r.distances[slot]-b.rho, weight = delta <= 0 ? 1 : exp(-delta/b.sigma)
                if weight > 0 {
                    try directed.append(row: i, feature: r.indices[slot], bits: weight.bitPattern)
                    degrees[r.indices[slot]] += 1
                }
            }
        }
        let bandwidthHash = try bandwidths.finish(); _ = try directed.finish()
        var starts = [Int](repeating: 0, count: n+1)
        for i in 0..<n { starts[i+1] = starts[i]+degrees[i] }
        guard starts[n] == directed.entries else { throw VivoOmicsError.invalid("transpose degree sum") }
        // Reuse degrees as scatter cursors. Source-major visitation guarantees
        // each incoming row is ordered by source ID without an in-memory sort.
        for i in 0..<n { degrees[i] = starts[i] }
        let fd = open(incomingURL.path, O_RDWR | O_CREAT | O_EXCL | O_NOFOLLOW, 0o600)
        guard fd >= 0 else { throw VivoOmicsError.invalid("graph transpose creation") }
        defer { _ = close(fd) }
        guard ftruncate(fd, off_t(directed.entries*16)) == 0 else { throw VivoOmicsError.limit("graph transpose allocation") }
        let input = try VivoWindowedCountRecords(directedURL, entries: directed.entries)
        for index in 0..<directed.entries {
            if index % 4096 == 0 { try Task.checkCancellation() }
            let r = try input.record(index)
            var encoded = ((UInt64(r.feature) | UInt64(r.row)<<32).littleEndian, r.bits.littleEndian)
            try withUnsafeBytes(of: &encoded) { bytes in
                var done = 0
                while done < 16 {
                    let written = pwrite(fd, bytes.baseAddress!.advanced(by: done), 16-done, off_t(degrees[r.feature]*16+done))
                    if written < 0 && errno == EINTR { continue }
                    guard written > 0 else { throw VivoOmicsError.invalid("graph transpose write") }
                    done += written
                }
            }
            degrees[r.feature] += 1
        }
        for i in 0..<n where degrees[i] != starts[i+1] { throw VivoOmicsError.invalid("graph transpose cursor") }
        let incoming = try VivoWindowedCountRecords(incomingURL, entries: directed.entries)
        let edges = try VivoCountRecordWriter(root.appendingPathComponent("edges.bin")), offsets = try VivoCountRecordWriter(root.appendingPathComponent("offsets.bin"))
        var parent = Array(0..<n), components = n, isolated = 0
        func component(_ v: Int) -> Int {
            var v = v
            while parent[v] != v { parent[v] = parent[parent[v]]; v = parent[v] }
            return v
        }
        var directedCursor = 0
        for i in 0..<n {
            try Task.checkCancellation()
            let before = edges.entries
            try offsets.append(row: i, feature: 0, bits: UInt64(before))
            var forward: [(index: Int, weight: Double)] = []
            while directedCursor < directed.entries {
                let r = try input.record(directedCursor)
                if r.row != i { break }
                forward.append((r.feature, Double(bitPattern: r.bits))); directedCursor += 1
            }
            forward.sort { $0.index < $1.index }
            var a = 0, b = starts[i]
            while a < forward.count || b < starts[i+1] {
                let reverse = b < starts[i+1] ? try incoming.record(b) : nil
                if let reverse {
                    guard reverse.row == i, reverse.feature < n, reverse.feature != i else { throw VivoOmicsError.invalid("transpose coordinates") }
                }
                let j = min(a < forward.count ? forward[a].index : n, reverse?.feature ?? n)
                let f = a < forward.count && forward[a].index == j ? forward[a].weight : 0
                let r = reverse?.feature == j ? Double(bitPattern: reverse!.bits) : 0
                let weight = min(1, max(0, f+r-f*r))
                guard j < n, weight.isFinite, weight > 0 else { throw VivoOmicsError.invalid("binary fuzzy edge") }
                try edges.append(row: i, feature: j, bits: weight.bitPattern)
                let x = component(i), y = component(j)
                if x != y { parent[max(x,y)] = min(x,y); components -= 1 }
                if a < forward.count && forward[a].index == j { a += 1 }
                if reverse?.feature == j { b += 1 }
            }
            if edges.entries == before { isolated += 1 }
        }
        guard directedCursor == directed.entries else { throw VivoOmicsError.invalid("binary directed traversal") }
        try offsets.append(row: n, feature: 0, bits: UInt64(edges.entries))
        let offsetsHash = try offsets.finish(), edgesHash = try edges.finish()
        return .init(schemaVersion: 1, format: "PCA-knn-fuzzy-CSR-records-u32-u32-u64-le/v1",
            method: approximate ? "approximate-HNSW-euclidean-PCA-knn-umap-fuzzy-union-v1" : "exact-euclidean-PCA-knn-umap-fuzzy-union-v1",
            cells: n, dimensions: dimensions, options: options, neighborEntries: n*k, connectivityEntries: edges.entries,
            connectedComponents: components, isolatedCells: isolated, distancePairs: distancePairs,
            neighbors: neighborsHash, bandwidths: bandwidthHash, offsets: offsetsHash, edges: edgesHash,
            qualification: "Fixed-width neighbors, bandwidths and CSR edges are binary files; cell identities are bound by the input receipt. Disk transpose and row merge use bounded edge buffers and 16 MiB file windows. Degree, offset, component arrays, input PCA state and any HNSW index remain resident. Approximate distancePairs counts repeated metric calls; exact counts unique pairs. No million-cell, downstream biological or Metal qualification.")
    }
}
