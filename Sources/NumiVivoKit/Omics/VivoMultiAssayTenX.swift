import Foundation

public enum VivoMultiAssayTenX {
    public static func validate(_ plan: VivoTenXMultiAssayPlan) throws {
        try plan.sample.validate()
        guard plan.schemaVersion == 1, vivoOmicsID(plan.id), !plan.sourceDescription.isEmpty,
              plan.sourceDescription.utf8.count <= 16_384, !plan.assays.isEmpty, plan.assays.count <= 16,
              Set(plan.assays.map(\.id)).count == plan.assays.count,
              Set(plan.assays.map(\.featureType)).count == plan.assays.count else { throw VivoOmicsError.invalid("10x multi-assay plan") }
        let kinds: [String: VivoAssayKind] = ["Gene Expression": .rna, "Antibody Capture": .antibodyCapture,
            "Peaks": .chromatinAccessibility, "CRISPR Guide Capture": .guideCapture]
        for a in plan.assays {
            guard vivoOmicsID(a.id), !a.id.contains("/"), a.id != ".", a.id != "..",
                  vivoOmicsID(a.featureNamespace), kinds[a.featureType] == a.kind,
                  a.genomeAssembly.map(vivoOmicsID) ?? true else { throw VivoOmicsError.invalid("10x feature-type mapping") }
            if a.kind == .chromatinAccessibility {
                guard a.genomeAssembly != nil, a.countUnit != .umiCount else { throw VivoOmicsError.invalid("peak assembly/count units") }
            } else if (a.countUnit == .fragmentCount || a.countUnit == .cutSiteCount) { throw VivoOmicsError.invalid("non-peak fragment/cut-site units") }
        }
    }
    /// Caller supplies an immutable snapshot; only one barcode's sparse entries
    /// and bounded HDF5 slices are transient. Final per-assay CSR arrays are resident.
    public static func readSnapshot(_ source: URL, plan: VivoTenXMultiAssayPlan) throws -> VivoMultiAssayDataset {
        try validate(plan)
        let limits = VivoMultiAssayDataset.limits
        return try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(source.path)
            defer { h.close(file, "H5Fclose") }
            func strings(_ path: String, maximum: Int) throws -> [String] {
                let d = try h.dataset(file, "matrix/" + path); defer { h.close(d, "H5Dclose") }
                guard try h.shape(d, attribute: false).count == 1 else { throw VivoOmicsError.invalid("10x vector rank") }
                return try h.strings(d, maximum: maximum)
            }
            func integers(_ path: String, maximum: Int) throws -> [UInt64] {
                let d = try h.dataset(file, "matrix/" + path); defer { h.close(d, "H5Dclose") }
                guard try h.shape(d, attribute: false).count == 1 else { throw VivoOmicsError.invalid("10x vector rank") }
                return try h.integers(d, maximum: maximum)
            }
            let shape = try integers("shape", maximum: 2)
            guard shape.count == 2, shape[0] > 0, shape[0] <= UInt64(limits.maximumFeatures),
                  shape[1] <= UInt64(limits.maximumCells) else { throw VivoOmicsError.limit("10x shape") }
            let barcodes = try strings("barcodes", maximum: limits.maximumCells)
            let ids = try strings("features/id", maximum: limits.maximumFeatures)
            let names = try strings("features/name", maximum: limits.maximumFeatures)
            let types = try strings("features/feature_type", maximum: limits.maximumFeatures)
            let genomes = try strings("features/genome", maximum: limits.maximumFeatures)
            guard barcodes.count == Int(shape[1]), ids.count == Int(shape[0]), names.count == ids.count,
                  types.count == ids.count, genomes.count == ids.count,
                  Set(types) == Set(plan.assays.map(\.featureType)) else { throw VivoOmicsError.invalid("10x axes or incomplete feature-type mapping") }
            let byType = Dictionary(uniqueKeysWithValues: plan.assays.enumerated().map { ($0.element.featureType, $0.offset) })
            var features = [[VivoAssayFeature]](repeating: [], count: plan.assays.count)
            var localIndex = [Int](repeating: 0, count: ids.count), spaceIndex = localIndex
            for i in ids.indices {
                let s = byType[types[i]]!, spec = plan.assays[s]
                if !genomes[i].isEmpty, genomes[i] != spec.genomeAssembly { throw VivoOmicsError.invalid("source genome differs from declared assembly") }
                var interval: VivoGenomicInterval?
                if spec.kind == .chromatinAccessibility {
                    let parts = ids[i].split(separator: ":", omittingEmptySubsequences: false)
                    guard parts.count == 2 else { throw VivoOmicsError.invalid("10x peak ID requires contig:start-end") }
                    let bounds = parts[1].split(separator: "-", omittingEmptySubsequences: false)
                    guard bounds.count == 2, let start = UInt64(bounds[0]), let end = UInt64(bounds[1]), start < end else { throw VivoOmicsError.invalid("10x peak interval") }
                    interval = .init(contig: String(parts[0]), start: start, end: end)
                }
                spaceIndex[i] = s; localIndex[i] = features[s].count
                features[s].append(.init(id: ids[i], name: names[i], interval: interval))
            }
            let offsets = try integers("indptr", maximum: barcodes.count + 1)
            let values = try h.dataset(file, "matrix/data"), indices = try h.dataset(file, "matrix/indices")
            defer { h.close(values, "H5Dclose"); h.close(indices, "H5Dclose") }
            let ds = try h.shape(values, attribute: false)
            guard ds.count == 1, ds[0] <= UInt64(limits.maximumNonzeros), try h.shape(indices, attribute: false) == ds,
                  offsets.count == barcodes.count + 1, offsets.first == 0, offsets.last == ds[0] else { throw VivoOmicsError.limit("10x sparse source bounds") }
            var rowOffsets = [[Int]](repeating: [0], count: plan.assays.count)
            var columns = [[Int]](repeating: [], count: plan.assays.count)
            var counts = [[UInt64]](repeating: [], count: plan.assays.count)
            for row in barcodes.indices {
                try Task.checkCancellation()
                guard offsets[row] <= offsets[row + 1], offsets[row + 1] <= ds[0] else { throw VivoOmicsError.invalid("10x offsets") }
                var entries: [Int: UInt64] = [:], start = Int(offsets[row])
                while start < Int(offsets[row + 1]) {
                    let end = min(start + 65_536, Int(offsets[row + 1]))
                    let ix = try h.integers(indices, maximum: 65_536, range: start..<end)
                    let vv = try h.integers(values, maximum: 65_536, range: start..<end)
                    for k in ix.indices {
                        guard ix[k] < shape[0] else { throw VivoOmicsError.invalid("10x feature index") }
                        if vv[k] > 0 { entries[Int(ix[k])] = try vivoOmicsSum(entries[Int(ix[k])] ?? 0, vv[k]) }
                    }
                    start = end
                }
                for i in entries.keys.sorted() {
                    let s = spaceIndex[i]
                    columns[s].append(localIndex[i]); counts[s].append(entries[i]!)
                }
                for s in plan.assays.indices { rowOffsets[s].append(counts[s].count) }
            }
            let assays = plan.assays.enumerated().map { s, a in
                VivoAssaySpace(id: a.id, kind: a.kind, featureNamespace: a.featureNamespace, countUnit: a.countUnit,
                    genomeAssembly: a.genomeAssembly, sourceDescription: plan.sourceDescription, features: features[s],
                    observationIndices: Array(barcodes.indices), matrix: .init(cellCount: barcodes.count, featureCount: features[s].count,
                        rowOffsets: rowOffsets[s], featureIndices: columns[s], counts: counts[s]))
            }
            let result = VivoMultiAssayDataset(schemaVersion: 1, id: plan.id, evidence: plan.evidence, sourceDescription: plan.sourceDescription,
                samples: [plan.sample], observations: barcodes.map { .init(identity: .init(barcode: $0, sampleID: plan.sample.id), kind: .cell, position: nil) },
                spatialFrames: [], assays: assays)
            try result.validate(); return result
        }
    }
}
