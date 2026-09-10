import Foundation

public struct VivoSingleCellProgramMember: Codable, Sendable, Equatable {
    public let featureID: String
    public let weight: Double
    public init(featureID: String, weight: Double = 1) { self.featureID=featureID;self.weight=weight }
    private enum CodingKeys: String,CodingKey { case featureID,weight }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["featureID","weight"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        featureID=try c.decode(String.self,forKey: .featureID);weight=try c.decode(Double.self,forKey: .weight)
    }
}

/// A supplied expression signature, not an inferred label or a pathway activity assay.
/// The namespace is a declaration; matching uses exact dataset feature IDs only.
public struct VivoSingleCellProgramDefinition: Codable, Sendable, Equatable {
    public let id: String
    public let organism: String
    public let featureNamespace: String
    public let sourceURI: String
    public let sourceVersion: String
    public let sourceDescription: String
    public let members: [VivoSingleCellProgramMember]
    public var minimumWeightCoverage: Double = 1
    public init(id: String,organism: String,featureNamespace: String,sourceURI: String,sourceVersion: String,
                sourceDescription: String,members: [VivoSingleCellProgramMember],minimumWeightCoverage: Double = 1) {
        self.id=id;self.organism=organism;self.featureNamespace=featureNamespace;self.sourceURI=sourceURI
        self.sourceVersion=sourceVersion;self.sourceDescription=sourceDescription;self.members=members
        self.minimumWeightCoverage=minimumWeightCoverage
    }
    private enum CodingKeys: String,CodingKey { case id,organism,featureNamespace,sourceURI,sourceVersion,sourceDescription,members,minimumWeightCoverage }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["id","organism","featureNamespace","sourceURI","sourceVersion","sourceDescription","members","minimumWeightCoverage"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        id=try c.decode(String.self,forKey: .id);organism=try c.decode(String.self,forKey: .organism)
        featureNamespace=try c.decode(String.self,forKey: .featureNamespace);sourceURI=try c.decode(String.self,forKey: .sourceURI)
        sourceVersion=try c.decode(String.self,forKey: .sourceVersion);sourceDescription=try c.decode(String.self,forKey: .sourceDescription)
        members=try c.decode([VivoSingleCellProgramMember].self,forKey: .members)
        minimumWeightCoverage=try c.decodeIfPresent(Double.self,forKey: .minimumWeightCoverage) ?? 1
    }
    public func validate() throws {
        guard [id,organism,featureNamespace,sourceURI,sourceVersion,sourceDescription].allSatisfy(vivoOmicsID),
              organism != "unreported", URL(string: sourceURI)?.scheme != nil,
              !members.isEmpty,members.count<=10_000,Set(members.map(\.featureID)).count==members.count,
              members.allSatisfy({ vivoOmicsID($0.featureID) && $0.weight.isFinite && $0.weight != 0 }),
              minimumWeightCoverage.isFinite,minimumWeightCoverage>0,minimumWeightCoverage<=1 else {
            throw VivoOmicsError.invalid("program definition, provenance, weights or coverage")
        }
    }
}

public struct VivoSingleCellProgramOptions: Codable, Sendable, Equatable {
    public let definitions: [VivoSingleCellProgramDefinition]
    public var maximumScoreValues: Int = 5_000_000
    public var maximumUpdates: Int = 100_000_000
    public init(definitions: [VivoSingleCellProgramDefinition]) { self.definitions=definitions }
    private enum CodingKeys: String,CodingKey { case definitions,maximumScoreValues,maximumUpdates }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["definitions","maximumScoreValues","maximumUpdates"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        definitions=try c.decode([VivoSingleCellProgramDefinition].self,forKey: .definitions)
        maximumScoreValues=try c.decodeIfPresent(Int.self,forKey: .maximumScoreValues) ?? 5_000_000
        maximumUpdates=try c.decodeIfPresent(Int.self,forKey: .maximumUpdates) ?? 100_000_000
    }
    public func validate() throws {
        guard !definitions.isEmpty,definitions.count<=256,Set(definitions.map(\.id)).count==definitions.count,
              (1...20_000_000).contains(maximumScoreValues),(1...1_000_000_000).contains(maximumUpdates),
              definitions.reduce(0, { $0+$1.members.count })<=100_000 else { throw VivoOmicsError.limit("program definitions or budgets") }
        for d in definitions { try d.validate() }
    }
}

public struct VivoSingleCellResolvedProgram: Codable, Sendable, Equatable {
    public let definition: VivoSingleCellProgramDefinition
    public let definitionFingerprint: VivoFingerprint
    /// Source feature order, paired with effectiveWeights. Missing members never become zero-expression observations.
    public let featureIndices: [Int]
    public let effectiveWeights: [Double]
    public let missingFeatureIDs: [String]
    public let weightCoverage: Double
}
public struct VivoSingleCellProgramResult: Codable, Sendable, Equatable {
    public let method: String
    public let normalizationTarget: Double
    public let programs: [VivoSingleCellResolvedProgram]
    public let cells: [VivoOmicsCellIdentity]
    /// Cells by programs. Nil for a zero-count library; measured zero expression otherwise gives zero.
    public let scores: [[Double?]]
    public let detectedMembers: [[Int]]
    public let updates: Int
    public let qualification: String
}

/// Sparse accumulator shared by resident CSR and canonical streamed CSR/CSC.
enum VivoProgramFeatureMatch: String, Codable, Sendable { case featureID, featureName }

final class VivoSingleCellProgramAccumulator {
    let options: VivoSingleCellProgramOptions
    let programs: [VivoSingleCellResolvedProgram]
    let cells: [VivoOmicsCellIdentity]
    let normalizationTarget: Double
    let hasLibrary: [Bool]
    let memberships: [[(program: Int,weight: Double)]]
    var values: [Double]
    var detected: [Int]
    var updates=0
    init(features: [VivoOmicsFeature],cells: [VivoOmicsCell],samples: [VivoOmicsSample],quality: [VivoCellQuality],
         normalizationTarget: Double,options: VivoSingleCellProgramOptions,featureMatch: VivoProgramFeatureMatch = .featureID) throws {
        try options.validate()
        guard normalizationTarget.isFinite,normalizationTarget>0,quality.count==cells.count,
              zip(cells,quality).allSatisfy({ $0.sampleID==$1.sampleID && $0.barcode==$1.barcode }),
              Set(features.map(\.id)).count==features.count else { throw VivoOmicsError.invalid("program source axes or normalization") }
        let size=cells.count.multipliedReportingOverflow(by: options.definitions.count)
        guard !size.overflow,size.partialValue<=options.maximumScoreValues else { throw VivoOmicsError.limit("program score-value budget") }
        let organisms=Set(samples.map(\.organism))
        let keys=features.map { featureMatch == .featureID ? $0.id : $0.name }
        let positions=Dictionary(grouping: features.indices,by: { keys[$0] })
        var programs: [VivoSingleCellResolvedProgram]=[]
        var memberships=[[(program: Int,weight: Double)]](repeating: [],count: features.count)
        for (p,definition) in options.definitions.enumerated() {
            guard organisms==[definition.organism] else { throw VivoOmicsError.invalid("program organism does not match source samples") }
            // Scaling before summation makes the L1 denominator safe for large finite weights.
            let scale=definition.members.map { abs($0.weight) }.max()!
            let requested=Dictionary(uniqueKeysWithValues: definition.members.map { ($0.featureID,$0.weight/scale) })
            guard requested.values.allSatisfy({ $0 != 0 }) else { throw VivoOmicsError.invalid("program weight dynamic range underflow") }
            let total=definition.members.reduce(0.0) { $0+abs(requested[$1.featureID]!) }
            guard requested.keys.allSatisfy({ (positions[$0]?.count ?? 0)<=1 }) else {
                throw VivoOmicsError.invalid("ambiguous requested program feature name")
            }
            let selected=features.indices.filter { requested[keys[$0]] != nil }
            let matched=selected.reduce(0.0) { $0+abs(requested[keys[$1]]!) }
            let ids=Set(keys),missing=definition.members.map(\.featureID).filter { !ids.contains($0) }.sorted()
            let coverage=missing.isEmpty ? 1 : matched/total
            guard matched>0,coverage>=definition.minimumWeightCoverage,
                  definition.minimumWeightCoverage<1 || missing.isEmpty else { throw VivoOmicsError.invalid("program \(definition.id) missing-feature weight coverage") }
            let weights=selected.map { requested[keys[$0]]!/matched }
            guard weights.allSatisfy({ $0.isFinite && $0 != 0 }) else { throw VivoOmicsError.invalid("program effective weight underflow") }
            let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys,.withoutEscapingSlashes]
            programs.append(.init(definition: definition,definitionFingerprint: try VivoCanonicalJSON.fingerprint(encoder.encode(definition)),
                featureIndices: selected,effectiveWeights: weights,missingFeatureIDs: missing,weightCoverage: coverage))
            for (index,weight) in zip(selected,weights) { memberships[index].append((p,weight)) }
        }
        self.options=options;self.programs=programs;self.memberships=memberships;self.normalizationTarget=normalizationTarget
        self.cells=cells.map { .init(sampleID: $0.sampleID,barcode: $0.barcode) };hasLibrary=quality.map { $0.totalCounts>0 }
        values=Array(repeating: 0,count: size.partialValue);detected=Array(repeating: 0,count: size.partialValue)
    }
    func add(row: Int,feature: Int,logValue: Double) throws {
        guard row>=0,row<cells.count,feature>=0,feature<memberships.count,hasLibrary[row],logValue.isFinite,logValue>0 else {
            throw VivoOmicsError.invalid("program sparse entry or empty library")
        }
        let targets=memberships[feature]
        guard targets.count<=options.maximumUpdates-updates else { throw VivoOmicsError.limit("program sparse-update budget") }
        updates+=targets.count
        for target in targets {
            let index=row*programs.count+target.program
            values[index]+=logValue*target.weight;detected[index]+=1
        }
    }
    func finish() throws -> VivoSingleCellProgramResult {
        guard values.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("nonfinite program score") }
        return .init(method: "signed-l1-mean-log-normalized-expression-v1",normalizationTarget: normalizationTarget,programs: programs,cells: cells,
            scores: cells.indices.map { row in programs.indices.map { hasLibrary[row] ? values[row*programs.count+$0] : nil } },
            detectedMembers: cells.indices.map { row in Array(detected[(row*programs.count)..<((row+1)*programs.count)]) },updates: updates,
            qualification: "Supplied expression programs with exact-ID matching and declared provenance; descriptive scores, not fitted labels, calibrated pathway activity or biological ground truth.")
    }
}

enum VivoSingleCellPrograms {
    static func run(_ processed: VivoSingleCellProcessed,options: VivoSingleCellProgramOptions) throws -> VivoSingleCellProgramResult {
        let data=processed.dataset,normal=processed.normalized
        let accumulator=try VivoSingleCellProgramAccumulator(features: data.features,cells: data.cells,samples: data.samples,
            quality: processed.sourceCellIndices.map { processed.decisions[$0].quality },normalizationTarget: normal.targetSum,options: options)
        for row in data.cells.indices {
            try Task.checkCancellation()
            for k in normal.rowOffsets[row]..<normal.rowOffsets[row+1] { try accumulator.add(row: row,feature: normal.featureIndices[k],logValue: normal.values[k]) }
        }
        return try accumulator.finish()
    }
    static func run(snapshot: URL,mapping: VivoH5ADImportPlan,metadata: VivoSingleCellCountMetadata,quality: [VivoCellQuality],
                    normalizationTarget: Double,options: VivoSingleCellProgramOptions) throws -> VivoSingleCellProgramResult {
        return try accumulate(snapshot: snapshot,mapping: mapping,metadata: metadata,quality: quality,
                              normalizationTarget: normalizationTarget,options: options).finish()
    }
    static func accumulate(snapshot: URL,mapping: VivoH5ADImportPlan,metadata: VivoSingleCellCountMetadata,quality: [VivoCellQuality],
                           normalizationTarget: Double,options: VivoSingleCellProgramOptions,featureMatch: VivoProgramFeatureMatch = .featureID) throws -> VivoSingleCellProgramAccumulator {
        let accumulator=try VivoSingleCellProgramAccumulator(features: metadata.features,cells: metadata.cells,samples: metadata.samples,
            quality: quality,normalizationTarget: normalizationTarget,options: options,featureMatch: featureMatch)
        _ = try VivoSingleCellH5AD.scanSnapshot(snapshot,plan: mapping,limits: VivoH5ADPseudobulk.sourceLimits,onMetadata: {
            guard $0==metadata else { throw VivoOmicsError.invalid("program source metadata changed") }
        },onEntry: { row,feature,count in
            let value=log1p((Double(count)/Double(quality[row].totalCounts))*normalizationTarget)
            try accumulator.add(row: row,feature: feature,logValue: value)
        })
        return accumulator
    }
}
