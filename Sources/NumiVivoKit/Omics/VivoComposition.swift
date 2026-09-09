import Foundation

public struct VivoCompositionContext: Codable, Sendable, Equatable {
    public let id: String
    public let organism: String
    public let featureNamespace: String
    public let perturbationNamespace: String
    public let countUnit: VivoOmicsCountUnit
    public init(id: String,organism: String,featureNamespace: String,perturbationNamespace: String,countUnit: VivoOmicsCountUnit) {
        self.id=id;self.organism=organism;self.featureNamespace=featureNamespace;self.perturbationNamespace=perturbationNamespace;self.countUnit=countUnit
    }
    private enum CodingKeys: String,CodingKey,CaseIterable { case id,organism,featureNamespace,perturbationNamespace,countUnit }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: Set(CodingKeys.allCases.map(\.rawValue)))
        let c=try decoder.container(keyedBy: CodingKeys.self)
        id=try c.decode(String.self,forKey: .id);organism=try c.decode(String.self,forKey: .organism)
        featureNamespace=try c.decode(String.self,forKey: .featureNamespace);perturbationNamespace=try c.decode(String.self,forKey: .perturbationNamespace)
        countUnit=try c.decode(VivoOmicsCountUnit.self,forKey: .countUnit)
    }
    func validate() throws {
        guard [id,organism,featureNamespace,perturbationNamespace].allSatisfy(vivoOmicsID) else { throw VivoOmicsError.invalid("composition context identities") }
    }
}
public struct VivoCompositionTarget: Codable, Sendable, Equatable {
    public let id: String
    public let condition: String
    public init(id: String,condition: String) { self.id=id;self.condition=condition }
    private enum CodingKeys: String,CodingKey { case id,condition }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["id","condition"])
        let c=try decoder.container(keyedBy: CodingKeys.self);id=try c.decode(String.self,forKey: .id);condition=try c.decode(String.self,forKey: .condition)
    }
}
public struct VivoCompositionPreparation: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let context: VivoCompositionContext
    public let controlCondition: String
    public let targets: [VivoCompositionTarget]
    public let provenance: String
    public init(context: VivoCompositionContext,controlCondition: String,targets: [VivoCompositionTarget],provenance: String) {
        schemaVersion=1;self.context=context;self.controlCondition=controlCondition;self.targets=targets;self.provenance=provenance
    }
    private enum CodingKeys: String,CodingKey,CaseIterable { case schemaVersion,context,controlCondition,targets,provenance }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: Set(CodingKeys.allCases.map(\.rawValue)))
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);context=try c.decode(VivoCompositionContext.self,forKey: .context)
        controlCondition=try c.decode(String.self,forKey: .controlCondition);targets=try c.decode([VivoCompositionTarget].self,forKey: .targets)
        provenance=try c.decode(String.self,forKey: .provenance)
    }
    func validate() throws {
        try context.validate()
        guard schemaVersion==1,vivoOmicsID(controlCondition),(2...256).contains(targets.count),
              targets.allSatisfy({ vivoOmicsID($0.id) && vivoOmicsID($0.condition) && $0.condition != controlCondition }),
              Set(targets.map(\.id)).count==targets.count,Set(targets.map(\.condition)).count==targets.count,
              !provenance.isEmpty,provenance.utf8.count<=16_384 else { throw VivoOmicsError.invalid("composition training selection") }
    }
}
/// Rows are condition aggregates: control first, followed by the declared
/// single-target conditions. The sparse storage's cellCount is a row dimension,
/// exactly as in VivoPseudobulkCounts; these are not individual cells.
public struct VivoCompositionTraining: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let selection: VivoCompositionPreparation
    public let source: VivoFingerprint
    public let sourceReport: VivoFingerprint
    public let evidence: VivoOmicsEvidence
    public let featureIDs: [String]
    public let matrix: VivoSparseCounts
    public init(selection: VivoCompositionPreparation,source: VivoFingerprint,sourceReport: VivoFingerprint,evidence: VivoOmicsEvidence,featureIDs: [String],matrix: VivoSparseCounts) {
        schemaVersion=1;self.selection=selection;self.source=source;self.sourceReport=sourceReport;self.evidence=evidence;self.featureIDs=featureIDs;self.matrix=matrix
    }
    private enum CodingKeys: String,CodingKey,CaseIterable { case schemaVersion,selection,source,sourceReport,evidence,featureIDs,matrix }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: Set(CodingKeys.allCases.map(\.rawValue)))
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);selection=try c.decode(VivoCompositionPreparation.self,forKey: .selection)
        source=try c.decode(VivoFingerprint.self,forKey: .source);sourceReport=try c.decode(VivoFingerprint.self,forKey: .sourceReport)
        evidence=try c.decode(VivoOmicsEvidence.self,forKey: .evidence);featureIDs=try c.decode([String].self,forKey: .featureIDs)
        matrix=try c.decode(VivoSparseCounts.self,forKey: .matrix)
    }
    func validate() throws {
        try selection.validate()
        guard schemaVersion==1,!featureIDs.isEmpty,featureIDs.count<=100_000,
              featureIDs.allSatisfy(vivoOmicsID),Set(featureIDs).count==featureIDs.count,
              matrix.cellCount==selection.targets.count+1,matrix.featureCount==featureIDs.count,
              matrix.cellCount<=10_000_000/featureIDs.count else { throw VivoOmicsError.limit("composition condition-feature dimensions") }
        var limits=VivoOmicsLimits();limits.maximumCells=257;limits.maximumFeatures=100_000
        try matrix.validate(limits: limits)
    }
}
public struct VivoCompositionQuery: Codable, Sendable, Equatable {
    public let id: String
    public let targets: [String]
    public init(id: String,targets: [String]) { self.id=id;self.targets=targets }
    private enum CodingKeys: String,CodingKey { case id,targets }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["id","targets"])
        let c=try decoder.container(keyedBy: CodingKeys.self);id=try c.decode(String.self,forKey: .id);targets=try c.decode([String].self,forKey: .targets)
    }
}
public struct VivoCompositionQueryPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let context: VivoCompositionContext
    public let queries: [VivoCompositionQuery]
    public init(context: VivoCompositionContext,queries: [VivoCompositionQuery]) { schemaVersion=1;self.context=context;self.queries=queries }
    private enum CodingKeys: String,CodingKey { case schemaVersion,context,queries }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","context","queries"])
        let c=try decoder.container(keyedBy: CodingKeys.self);schemaVersion=try c.decode(Int.self,forKey: .schemaVersion)
        context=try c.decode(VivoCompositionContext.self,forKey: .context);queries=try c.decode([VivoCompositionQuery].self,forKey: .queries)
    }
}
public enum VivoCompositionMethod: String,Codable,Sendable,CaseIterable {
    case noChange,meanSingleResponse,additiveLogResponse,meanConstituentLogResponse,additiveCPMResponse,shuffledAdditiveLogResponse
}
public struct VivoCompositionModel: Codable,Sendable,Equatable {
    public let schemaVersion: Int
    public let method: String
    public let context: VivoCompositionContext
    public let featureIDs: [String]
    public let targetIDs: [String]
    public let controlCPM: [Double]
    public let targetCPM: [[Double]]
}
public struct VivoCompositionDiagnostic: Codable,Sendable,Equatable {
    public let clippedFeatures: Int
    public let impliedCPMSum: Double
    public let renormalized: Bool
}
public struct VivoCompositionPredictions: Codable,Sendable,Equatable {
    public let method: VivoCompositionMethod
    public let queryRows: [Int]
    public let expression: [[Double]]
    public let diagnostics: [VivoCompositionDiagnostic]
}
public struct VivoCompositionReport: Codable,Sendable,Equatable {
    public let schemaVersion: Int
    public let context: VivoCompositionContext
    public let featureIDs: [String]
    public let queryIDs: [String]
    public let methods: [VivoCompositionPredictions]
}

public enum VivoComposition {
    public static func select(_ report: VivoH5ADPseudobulkReport,plan: VivoCompositionPreparation,source: VivoFingerprint,sourceReport: VivoFingerprint) throws -> VivoCompositionTraining {
        try plan.validate()
        let bulk=report.pseudobulk
        guard bulk.countUnit==plan.context.countUnit else { throw VivoOmicsError.invalid("composition count unit mismatch") }
        let targets=plan.targets.sorted { $0.id<$1.id }
        let conditions=[plan.controlCondition]+targets.map(\.condition)
        let rows=try conditions.map { condition -> Int in
            let matching=bulk.groups.indices.filter { bulk.groups[$0].condition==condition }
            guard matching.count==1 else { throw VivoOmicsError.invalid("composition requires exactly one declared aggregate per condition") }
            return matching[0]
        }
        let first=bulk.groups[rows[0]]
        guard first.organism==plan.context.organism,rows.allSatisfy({
            let g=bulk.groups[$0]
            return g.organism==first.organism && g.cellGroup==first.cellGroup && g.donorID==first.donorID && g.biologicalReplicateID==first.biologicalReplicateID
        }) else { throw VivoOmicsError.invalid("composition conditions cross biological contexts") }
        var limits=VivoOmicsLimits();limits.maximumCells=1_000_000
        try bulk.matrix.validate(limits: limits)
        guard bulk.matrix.cellCount==bulk.groups.count,bulk.matrix.featureCount==bulk.featureIDs.count else { throw VivoOmicsError.invalid("composition aggregate dimensions") }
        var offsets=[0],indices: [Int]=[],values: [UInt64]=[]
        for row in rows {
            let range=bulk.matrix.rowOffsets[row]..<bulk.matrix.rowOffsets[row+1]
            indices.append(contentsOf: bulk.matrix.featureIndices[range]);values.append(contentsOf: bulk.matrix.counts[range]);offsets.append(values.count)
        }
        let selection=VivoCompositionPreparation(context: plan.context,controlCondition: plan.controlCondition,targets: targets,provenance: plan.provenance)
        let training=VivoCompositionTraining(selection: selection,source: source,sourceReport: sourceReport,evidence: report.metadata.evidence,featureIDs: bulk.featureIDs,
            matrix: .init(cellCount: rows.count,featureCount: bulk.featureIDs.count,rowOffsets: offsets,featureIndices: indices,counts: values))
        try training.validate();return training
    }
    public static func fit(_ training: VivoCompositionTraining) throws -> VivoCompositionModel {
        try training.validate()
        let m=training.matrix,p=training.featureIDs.count
        func normalized(_ row: Int) throws -> [Double] {
            try Task.checkCancellation()
            let range=m.rowOffsets[row]..<m.rowOffsets[row+1]
            var total: UInt64=0
            for k in range { total=try vivoOmicsSum(total,m.counts[k]) }
            guard total>0,total<=9_007_199_254_740_992 else { throw VivoOmicsError.invalid("composition empty or inexact library total") }
            var result=[Double](repeating: 0,count: p)
            for k in range { result[m.featureIndices[k]]=Double(m.counts[k])/Double(total)*1_000_000 }
            return result
        }
        let order=training.selection.targets.indices.sorted { training.selection.targets[$0].id<training.selection.targets[$1].id }
        return try .init(schemaVersion: 1,method: "native-condition-CPM-composition-v1",context: training.selection.context,
            featureIDs: training.featureIDs,targetIDs: order.map { training.selection.targets[$0].id },
            controlCPM: normalized(0),targetCPM: order.map { try normalized($0+1) })
    }
    public static func predict(_ model: VivoCompositionModel,plan: VivoCompositionQueryPlan) throws -> VivoCompositionReport {
        try model.context.validate();try plan.context.validate()
        let n=model.targetIDs.count,p=model.featureIDs.count,q=plan.queries.count
        guard model.schemaVersion==1,model.method=="native-condition-CPM-composition-v1",plan.schemaVersion==1,plan.context==model.context,
              (2...256).contains(n),(1...100_000).contains(p),(1...256).contains(q),
              n+1<=10_000_000/p,2+4*q<=20_000_000/p,
              model.targetIDs==model.targetIDs.sorted(),Set(model.targetIDs).count==n,model.targetIDs.allSatisfy(vivoOmicsID),
              Set(model.featureIDs).count==p,model.featureIDs.allSatisfy(vivoOmicsID),
              model.controlCPM.count==p,model.targetCPM.count==n,model.targetCPM.allSatisfy({ $0.count==p }),
              Set(plan.queries.map(\.id)).count==q,plan.queries.allSatisfy({ vivoOmicsID($0.id) }) else {
            throw VivoOmicsError.invalid("composition model, query context or resource bounds")
        }
        for row in [model.controlCPM]+model.targetCPM {
            let tolerance=max(1e-8,Double(p)*Double.ulpOfOne*1_000_000*4)
            guard row.allSatisfy({ $0.isFinite && $0>=0 }),abs(row.reduce(0,+)-1_000_000)<=tolerance else { throw VivoOmicsError.invalid("composition normalized library") }
        }
        let byID=Dictionary(uniqueKeysWithValues: model.targetIDs.enumerated().map { ($0.element,$0.offset) })
        let pairs=try plan.queries.map { query -> (Int,Int) in
            guard query.targets.count==2,query.targets[0] != query.targets[1],let a=byID[query.targets[0]],let b=byID[query.targets[1]] else {
                throw VivoOmicsError.invalid("composition query needs two distinct observed targets")
            };return (a,b)
        }
        let baseline=model.controlCPM.map(log1p)
        let delta=model.targetCPM.map { row in row.indices.map { log1p(row[$0])-baseline[$0] } }
        var mean=[Double](repeating: 0,count: p)
        for row in delta { for j in 0..<p { mean[j]+=row[j] } }
        for j in 0..<p { mean[j]/=Double(n) }
        var methods: [VivoCompositionPredictions]=[]
        for method in VivoCompositionMethod.allCases {
            let shared=method == .noChange || method == .meanSingleResponse
            var expression: [[Double]]=[],diagnostics: [VivoCompositionDiagnostic]=[]
            for i in 0..<(shared ? 1 : q) {
                try Task.checkCancellation()
                let (a,b)=pairs[i]
                var values=[Double](repeating: 0,count: p),clipped=0,sum=0.0
                for j in 0..<p {
                    let raw: Double
                    switch method {
                    case .noChange:raw=baseline[j]
                    case .meanSingleResponse:raw=baseline[j]+mean[j]
                    case .additiveLogResponse:raw=baseline[j]+(delta[a][j]+delta[b][j])
                    case .meanConstituentLogResponse:raw=baseline[j]+(delta[a][j]+delta[b][j])/2
                    case .shuffledAdditiveLogResponse:raw=baseline[j]+(delta[(a+1)%n][j]+delta[(b+1)%n][j])
                    case .additiveCPMResponse:
                        let cpm=model.controlCPM[j]+((model.targetCPM[a][j]-model.controlCPM[j])+(model.targetCPM[b][j]-model.controlCPM[j]))
                        if cpm<0 { clipped+=1 };raw=log1p(max(0,cpm))
                    }
                    if raw<0 { clipped+=1 };values[j]=max(0,raw);sum+=expm1(values[j])
                }
                guard sum.isFinite,values.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("composition nonfinite prediction") }
                expression.append(values);diagnostics.append(.init(clippedFeatures: clipped,impliedCPMSum: sum,renormalized: false))
            }
            methods.append(.init(method: method,queryRows: shared ? [Int](repeating: 0,count: q) : Array(0..<q),expression: expression,
                diagnostics: shared ? Array(repeating: diagnostics[0],count: q) : diagnostics))
        }
        return .init(schemaVersion: 1,context: model.context,featureIDs: model.featureIDs,queryIDs: plan.queries.map(\.id),methods: methods)
    }
}
