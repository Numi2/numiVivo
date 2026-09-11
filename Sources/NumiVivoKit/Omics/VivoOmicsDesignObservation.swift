import Foundation

/// Resolves membership through a content-bound count bundle's QC group ordinals.
/// The count bundle retains the cell axis; this reference does not invent cell indices.
public struct VivoOmicsFileMembership: Codable, Sendable, Equatable {
    public let countBundleReceipt: VivoFingerprint
    public let groupIndex: Int
    public let sourceCellCount: Int
    init(countBundleReceipt: VivoFingerprint, groupIndex: Int, sourceCellCount: Int) {
        self.countBundleReceipt = countBundleReceipt; self.groupIndex = groupIndex; self.sourceCellCount = sourceCellCount
    }
    private enum CodingKeys: String, CodingKey { case countBundleReceipt, groupIndex, sourceCellCount }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["countBundleReceipt", "groupIndex", "sourceCellCount"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        countBundleReceipt = try c.decode(VivoFingerprint.self, forKey: .countBundleReceipt)
        groupIndex = try c.decode(Int.self, forKey: .groupIndex); sourceCellCount = try c.decode(Int.self, forKey: .sourceCellCount)
        try validate()
    }
    func validate() throws {
        guard countBundleReceipt.bytes.count == 32, (0..<100_000).contains(groupIndex),
              (1...20_000_000).contains(sourceCellCount) else {
            throw VivoOmicsError.invalid("file-backed design membership")
        }
    }
}

public enum VivoOmicsObservationMembership: Sendable, Equatable {
    case indices([Int])
    case file(VivoOmicsFileMembership)
}

/// Shared biological observation identity for aggregate inference. Legacy JSON
/// retains its original sourceCellIndices; file-backed observations encode an
/// explicit membership reference and never materialize or fabricate that array.
public struct VivoOmicsDesignObservation: Codable, Sendable, Equatable {
    public let biologicalReplicateID: String
    public let donorID: String?
    public let condition: String
    public let organism: String
    public let cellGroup: String?
    public let sampleIDs: [String]
    public let batchIDs: [String]
    public let membership: VivoOmicsObservationMembership
    public var sourceCellCount: Int {
        switch membership { case .indices(let rows): return rows.count; case .file(let reference): return reference.sourceCellCount }
    }
    public var sourceCellIndices: [Int]? {
        if case .indices(let rows) = membership { return rows }; return nil
    }
    init(_ group: VivoPseudobulkGroup) {
        biologicalReplicateID = group.biologicalReplicateID; donorID = group.donorID
        condition = group.condition; organism = group.organism; cellGroup = group.cellGroup
        sampleIDs = group.sampleIDs; batchIDs = group.batchIDs; membership = .indices(group.sourceCellIndices)
    }
    init(_ group: VivoFilePseudobulkGroup, receipt: VivoFingerprint, index: Int) {
        biologicalReplicateID = group.biologicalReplicateID; donorID = group.donorID
        condition = group.condition; organism = group.organism; cellGroup = group.cellGroup
        sampleIDs = group.sampleIDs; batchIDs = group.batchIDs
        membership = .file(.init(countBundleReceipt: receipt, groupIndex: index, sourceCellCount: group.sourceCellCount))
    }
    func validate() throws {
        guard [biologicalReplicateID, condition, organism].allSatisfy(vivoOmicsID),
              donorID.map(vivoOmicsID) ?? true, cellGroup.map(vivoOmicsID) ?? true,
              !sampleIDs.isEmpty, sampleIDs.allSatisfy(vivoOmicsID), batchIDs.allSatisfy(vivoOmicsID),
              Set(sampleIDs).count == sampleIDs.count, Set(batchIDs).count == batchIDs.count else {
            throw VivoOmicsError.invalid("design observation identity")
        }
        switch membership {
        case .indices(let rows):
            var previous = -1
            for row in rows { guard row > previous else { throw VivoOmicsError.invalid("design source-cell ordering") }; previous = row }
        case .file(let reference): try reference.validate()
        }
    }
    private enum CodingKeys: String, CodingKey {
        case biologicalReplicateID, donorID, condition, organism, cellGroup, sampleIDs, batchIDs, sourceCellIndices, fileMembership
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["biologicalReplicateID", "donorID", "condition", "organism", "cellGroup", "sampleIDs", "batchIDs", "sourceCellIndices", "fileMembership"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        biologicalReplicateID = try c.decode(String.self, forKey: .biologicalReplicateID); donorID = try c.decodeIfPresent(String.self, forKey: .donorID)
        condition = try c.decode(String.self, forKey: .condition); organism = try c.decode(String.self, forKey: .organism)
        cellGroup = try c.decodeIfPresent(String.self, forKey: .cellGroup); sampleIDs = try c.decode([String].self, forKey: .sampleIDs); batchIDs = try c.decode([String].self, forKey: .batchIDs)
        guard c.contains(.sourceCellIndices) != c.contains(.fileMembership) else { throw VivoOmicsError.invalid("exactly one observation membership representation is required") }
        if c.contains(.sourceCellIndices) { membership = .indices(try c.decode([Int].self, forKey: .sourceCellIndices)) }
        else { membership = .file(try c.decode(VivoOmicsFileMembership.self, forKey: .fileMembership)) }
        try validate()
    }
    public func encode(to encoder: Encoder) throws {
        var c = encoder.container(keyedBy: CodingKeys.self)
        try c.encode(biologicalReplicateID, forKey: .biologicalReplicateID); try c.encodeIfPresent(donorID, forKey: .donorID)
        try c.encode(condition, forKey: .condition); try c.encode(organism, forKey: .organism); try c.encodeIfPresent(cellGroup, forKey: .cellGroup)
        try c.encode(sampleIDs, forKey: .sampleIDs); try c.encode(batchIDs, forKey: .batchIDs)
        switch membership {
        case .indices(let rows): try c.encode(rows, forKey: .sourceCellIndices)
        case .file(let reference): try c.encode(reference, forKey: .fileMembership)
        }
    }
}
