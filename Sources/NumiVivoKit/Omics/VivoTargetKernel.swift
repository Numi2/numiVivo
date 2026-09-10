import Foundation

public struct VivoTargetDescriptor: Codable, Sendable, Equatable {
    public enum Status: String, Codable, Sendable { case available, noData, unresolvedIdentity }
    public let targetID: String
    public let featureID: String?
    public let status: Status
    public let terms: [String]
    public init(targetID: String, featureID: String?, status: Status, terms: [String]) {
        self.targetID=targetID; self.featureID=featureID; self.status=status; self.terms=terms
    }
    private enum CodingKeys: String, CodingKey { case targetID, featureID, status, terms }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["targetID","featureID","status","terms"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        targetID=try c.decode(String.self,forKey: .targetID);featureID=try c.decodeIfPresent(String.self,forKey: .featureID)
        status=try c.decode(Status.self,forKey: .status);terms=try c.decode([String].self,forKey: .terms)
    }
    func validate() throws {
        guard vivoOmicsID(targetID),featureID.map(vivoOmicsID) ?? true,terms.count<=4096,
              terms==terms.sorted(),Set(terms).count==terms.count,terms.allSatisfy(vivoOmicsID) else { throw VivoOmicsError.invalid("target descriptor identities or terms") }
        switch status {
        case .available: guard featureID != nil,!terms.isEmpty else { throw VivoOmicsError.invalid("available target descriptor needs identity and terms") }
        case .noData: guard featureID != nil,terms.isEmpty else { throw VivoOmicsError.invalid("no-data descriptor") }
        case .unresolvedIdentity: guard featureID==nil,terms.isEmpty else { throw VivoOmicsError.invalid("unresolved descriptor") }
        }
    }
}
public struct VivoTargetKernelPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let context: VivoCompositionContext
    public let descriptorNamespace: String
    public let source: VivoFingerprint
    public let provenance: String
    public let targets: [VivoTargetDescriptor]
    public let regularization: Double
    public let maximumWork: Int
    public init(context: VivoCompositionContext, descriptorNamespace: String, source: VivoFingerprint,
                provenance: String, targets: [VivoTargetDescriptor], regularization: Double = 1, maximumWork: Int = 200_000_000) {
        schemaVersion=1;self.context=context;self.descriptorNamespace=descriptorNamespace;self.source=source
        self.provenance=provenance;self.targets=targets;self.regularization=regularization;self.maximumWork=maximumWork
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,context,descriptorNamespace,source,provenance,targets,regularization,maximumWork }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","context","descriptorNamespace","source","provenance","targets","regularization","maximumWork"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);context=try c.decode(VivoCompositionContext.self,forKey: .context)
        descriptorNamespace=try c.decode(String.self,forKey: .descriptorNamespace);source=try c.decode(VivoFingerprint.self,forKey: .source)
        provenance=try c.decode(String.self,forKey: .provenance);targets=try c.decode([VivoTargetDescriptor].self,forKey: .targets)
        regularization=try c.decodeIfPresent(Double.self,forKey: .regularization) ?? 1
        maximumWork=try c.decodeIfPresent(Int.self,forKey: .maximumWork) ?? 200_000_000
    }
    func validate() throws {
        try context.validate()
        guard schemaVersion==1,vivoOmicsID(descriptorNamespace),!provenance.isEmpty,provenance.utf8.count<=16_384,
              (2...256).contains(targets.count),Set(targets.map(\.targetID)).count==targets.count,
              targets.reduce(0,{ $0+$1.terms.count })<=100_000,
              regularization.isFinite,(0.001...1000).contains(regularization),(1...100_000_000_000).contains(maximumWork) else { throw VivoOmicsError.invalid("target-kernel plan") }
        for target in targets { try target.validate() }
        let features=targets.compactMap(\.featureID)
        guard Set(features).count==features.count else { throw VivoOmicsError.invalid("target aliases share a feature identity") }
    }
}
public struct VivoTargetKernelQuery: Codable, Sendable, Equatable {
    public let id: String
    public let descriptor: VivoTargetDescriptor
    public init(id: String,descriptor: VivoTargetDescriptor) { self.id=id;self.descriptor=descriptor }
    private enum CodingKeys: String,CodingKey { case id,descriptor }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["id","descriptor"])
        let c=try decoder.container(keyedBy: CodingKeys.self);id=try c.decode(String.self,forKey: .id);descriptor=try c.decode(VivoTargetDescriptor.self,forKey: .descriptor)
    }
}
public struct VivoTargetKernelQueryPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let context: VivoCompositionContext
    public let descriptorNamespace: String
    public let source: VivoFingerprint
    public let queries: [VivoTargetKernelQuery]
    public let maximumWork: Int
    public init(context: VivoCompositionContext,descriptorNamespace: String,source: VivoFingerprint,queries: [VivoTargetKernelQuery],maximumWork: Int = 200_000_000) {
        schemaVersion=1;self.context=context;self.descriptorNamespace=descriptorNamespace;self.source=source;self.queries=queries;self.maximumWork=maximumWork
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,context,descriptorNamespace,source,queries,maximumWork }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","context","descriptorNamespace","source","queries","maximumWork"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);context=try c.decode(VivoCompositionContext.self,forKey: .context)
        descriptorNamespace=try c.decode(String.self,forKey: .descriptorNamespace);source=try c.decode(VivoFingerprint.self,forKey: .source)
        queries=try c.decode([VivoTargetKernelQuery].self,forKey: .queries);maximumWork=try c.decodeIfPresent(Int.self,forKey: .maximumWork) ?? 200_000_000
    }
}
public struct VivoTargetKernelModel: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let base: VivoCompositionModel
    public let descriptorNamespace: String
    public let descriptors: [VivoTargetDescriptor]
    public let supportedRows: [Int]
    public let regularization: Double
    public let inverse: [Double]
    public let maximumInverseResidual: Double
}
public struct VivoTargetKernelBaseline: Codable, Sendable, Equatable {
    public let method: String
    public let expression: [Double]
    public let diagnostic: VivoCompositionDiagnostic
}
public struct VivoTargetKernelQueryResult: Codable, Sendable, Equatable {
    public let id: String
    public let targetID: String
    public let status: String
    public let maximumSimilarity: Double
    public let weights: [Double]?
    public let expression: [Double]?
    public let shuffledExpression: [Double]?
    public let diagnostic: VivoCompositionDiagnostic?
    public let shuffledDiagnostic: VivoCompositionDiagnostic?
}
public struct VivoTargetKernelReport: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let context: VivoCompositionContext
    public let featureIDs: [String]
    public let supportedTrainingTargets: [String]
    public let baselines: [VivoTargetKernelBaseline]
    public let queries: [VivoTargetKernelQueryResult]
    public let qualification: String
}

public enum VivoTargetKernel {
    static let method="native-term-jaccard-ridge-intercept-Double-v1"
    static func similarity(_ a: [String],_ b: [String]) -> Double {
        var i=0,j=0,intersection=0
        while i<a.count && j<b.count {
            if a[i]==b[j] { intersection+=1;i+=1;j+=1 }
            else if a[i]<b[j] { i+=1 } else { j+=1 }
        }
        let union=a.count+b.count-intersection
        return union>0 ? Double(intersection)/Double(union) : 0
    }
    static func kernel(_ descriptors: [VivoTargetDescriptor],regularization: Double) -> [Double] {
        let n=descriptors.count
        var result=[Double](repeating: 0,count: n*n)
        for i in 0..<n { for j in 0...i {
            let value=similarity(descriptors[i].terms,descriptors[j].terms)+(i==j ? regularization : 0)
            result[i*n+j]=value;result[j*n+i]=value
        } }
        return result
    }
    static func inverse(_ matrix: [Double],n: Int) throws -> [Double] {
        var l=[Double](repeating: 0,count: n*n)
        for i in 0..<n { for j in 0...i {
            var value=matrix[i*n+j]
            for k in 0..<j { value-=l[i*n+k]*l[j*n+k] }
            if i==j {
                guard value.isFinite,value>0 else { throw VivoOmicsError.invalid("target-kernel positive-definite solve") }
                l[i*n+j]=sqrt(value)
            } else { l[i*n+j]=value/l[j*n+j] }
        } }
        var result=[Double](repeating: 0,count: n*n)
        for column in 0..<n {
            try Task.checkCancellation()
            var y=[Double](repeating: 0,count: n),x=y
            for i in 0..<n {
                var value=i==column ? 1.0 : 0.0
                for j in 0..<i { value-=l[i*n+j]*y[j] }
                y[i]=value/l[i*n+i]
            }
            for i in (0..<n).reversed() {
                var value=y[i]
                for j in (i+1)..<n { value-=l[j*n+i]*x[j] }
                x[i]=value/l[i*n+i]
            }
            for i in 0..<n { result[i*n+column]=x[i] }
        }
        return result
    }
    static func residual(_ matrix: [Double],inverse: [Double],n: Int) throws -> Double {
        var error=0.0
        for i in 0..<n { for j in 0..<n {
            var value=0.0
            for k in 0..<n { value+=matrix[i*n+k]*inverse[k*n+j] }
            error=max(error,abs(value-(i==j ? 1 : 0)))
        } }
        guard error.isFinite,error<1e-8,inverse.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("target-kernel inverse residual") }
        return error
    }
    public static func fit(_ training: VivoCompositionTraining,plan: VivoTargetKernelPlan) throws -> VivoTargetKernelModel {
        try training.validate();try plan.validate()
        guard plan.context==training.selection.context,Set(plan.targets.map(\.targetID))==Set(training.selection.targets.map(\.id)) else { throw VivoOmicsError.invalid("target-kernel training context or descriptor coverage") }
        let n=plan.targets.count,p=training.featureIDs.count
        guard 3*n*n*n+2*n*p<=plan.maximumWork else { throw VivoOmicsError.limit("target-kernel fit work budget") }
        let base=try VivoComposition.fit(training)
        let byID=Dictionary(uniqueKeysWithValues: plan.targets.map { ($0.targetID,$0) })
        let descriptors=base.targetIDs.map { byID[$0]! }
        let supported=descriptors.indices.filter { descriptors[$0].status == .available }
        guard supported.count>=2 else { throw VivoOmicsError.invalid("target-kernel requires two supported training targets") }
        let matrix=kernel(supported.map { descriptors[$0] },regularization: plan.regularization)
        let inverse=try Self.inverse(matrix,n: supported.count)
        return try .init(schemaVersion: 1,method: method,base: base,descriptorNamespace: plan.descriptorNamespace,descriptors: descriptors,
            supportedRows: supported,regularization: plan.regularization,inverse: inverse,maximumInverseResidual: residual(matrix,inverse: inverse,n: supported.count))
    }
    static func validate(_ model: VivoTargetKernelModel) throws {
        let base=model.base,n=base.targetIDs.count,p=base.featureIDs.count,s=model.supportedRows.count
        try base.context.validate()
        guard model.schemaVersion==1,model.method==method,base.schemaVersion==1,base.method=="native-condition-CPM-composition-v1",
              (2...256).contains(n),(1...100_000).contains(p),n+1<=10_000_000/p,
              base.targetIDs==base.targetIDs.sorted(),Set(base.targetIDs).count==n,base.targetIDs.allSatisfy(vivoOmicsID),
              Set(base.featureIDs).count==p,base.featureIDs.allSatisfy(vivoOmicsID),base.controlCPM.count==p,
              base.targetCPM.count==n,base.targetCPM.allSatisfy({ $0.count==p }),model.descriptors.count==n,
              model.descriptors.map(\.targetID)==base.targetIDs,vivoOmicsID(model.descriptorNamespace),
              model.regularization.isFinite,(0.001...1000).contains(model.regularization),
              (2...n).contains(s),model.supportedRows==model.descriptors.indices.filter({ model.descriptors[$0].status == .available }),
              model.inverse.count==s*s,model.maximumInverseResidual.isFinite,model.maximumInverseResidual>=0,model.maximumInverseResidual<1e-8 else { throw VivoOmicsError.invalid("target-kernel model dimensions or identities") }
        for descriptor in model.descriptors { try descriptor.validate() }
        let features=model.descriptors.compactMap(\.featureID)
        guard Set(features).count==features.count else { throw VivoOmicsError.invalid("target-kernel duplicate feature identity") }
        for row in [base.controlCPM]+base.targetCPM {
            guard row.allSatisfy({ $0.isFinite && $0>=0 }),abs(row.reduce(0,+)-1_000_000)<=max(1e-8,Double(p)*Double.ulpOfOne*4_000_000) else { throw VivoOmicsError.invalid("target-kernel normalized library") }
        }
        let actualResidual=try residual(kernel(model.supportedRows.map { model.descriptors[$0] },regularization: model.regularization),inverse: model.inverse,n: s)
        guard actualResidual==model.maximumInverseResidual else { throw VivoOmicsError.invalid("target-kernel residual diagnostic") }
    }
    static func clipped(_ values: [Double]) throws -> ([Double],VivoCompositionDiagnostic) {
        let clipped=values.map { max(0,$0) },sum=clipped.reduce(0) { $0+expm1($1) }
        guard values.allSatisfy(\.isFinite),sum.isFinite else { throw VivoOmicsError.invalid("target-kernel prediction nonfinite") }
        return (clipped,.init(clippedFeatures: values.filter { $0<0 }.count,impliedCPMSum: sum,renormalized: false))
    }
    public static func predict(_ model: VivoTargetKernelModel,plan: VivoTargetKernelQueryPlan) throws -> VivoTargetKernelReport {
        let n=model.base.targetIDs.count,p=model.base.featureIDs.count,s=model.supportedRows.count,q=plan.queries.count
        try plan.context.validate()
        guard plan.schemaVersion==1,plan.context==model.base.context,plan.descriptorNamespace==model.descriptorNamespace,
              (1...64).contains(q),(1...100_000_000_000).contains(plan.maximumWork),
              (1...100_000).contains(p),(2...256).contains(n),(2...256).contains(s),
              3+2*q<=20_000_000/p,s*s*s+2*q*s*p<=plan.maximumWork,
              Set(plan.queries.map(\.id)).count==q,plan.queries.allSatisfy({ vivoOmicsID($0.id) }),
              plan.queries.reduce(0,{ $0+$1.descriptor.terms.count })<=100_000 else { throw VivoOmicsError.invalid("target-kernel query context or resource budget") }
        try validate(model)
        let descriptors=model.supportedRows.map { model.descriptors[$0] }
        let knownFeatures=Set(model.descriptors.compactMap(\.featureID))
        for query in plan.queries {
            try query.descriptor.validate()
            guard !model.base.targetIDs.contains(query.descriptor.targetID),query.descriptor.featureID.map({ !knownFeatures.contains($0) }) ?? true else { throw VivoOmicsError.invalid("target-kernel query was observed in training") }
        }
        let baseline=model.base.controlCPM.map(log1p)
        let delta=model.base.targetCPM.map { row in row.indices.map { log1p(row[$0])-baseline[$0] } }
        var allMean=[Double](repeating: 0,count: p),supportedMean=allMean
        for row in delta { for j in 0..<p { allMean[j]+=row[j] } }
        for i in model.supportedRows { for j in 0..<p { supportedMean[j]+=delta[i][j] } }
        for j in 0..<p { allMean[j]=baseline[j]+allMean[j]/Double(n);supportedMean[j]=baseline[j]+supportedMean[j]/Double(s) }
        var baselines: [VivoTargetKernelBaseline]=[]
        for (name,values) in [("noChange",baseline),("meanSingleResponse",allMean),("meanSupportedResponse",supportedMean)] {
            let (expression,diagnostic)=try clipped(values);baselines.append(.init(method: name,expression: expression,diagnostic: diagnostic))
        }
        var u=[Double](repeating: 0,count: s)
        for i in 0..<s { for j in 0..<s { u[i]+=model.inverse[i*s+j] } }
        let total=u.reduce(0,+)
        guard total.isFinite,total>0 else { throw VivoOmicsError.invalid("target-kernel intercept solve") }
        var results: [VivoTargetKernelQueryResult]=[]
        for query in plan.queries {
            try Task.checkCancellation()
            let similarities=descriptors.map { similarity(query.descriptor.terms,$0.terms) },maximum=similarities.max() ?? 0
            guard query.descriptor.status == .available,maximum>0 else {
                results.append(.init(id: query.id,targetID: query.descriptor.targetID,status: query.descriptor.status == .available ? "noSharedTerms" : query.descriptor.status.rawValue,
                    maximumSimilarity: maximum,weights: nil,expression: nil,shuffledExpression: nil,diagnostic: nil,shuffledDiagnostic: nil));continue
            }
            var weights=[Double](repeating: 0,count: s)
            for j in 0..<s { for i in 0..<s { weights[j]+=similarities[i]*model.inverse[i*s+j] } }
            let adjustment=(1-weights.reduce(0,+))/total
            for j in 0..<s { weights[j]+=adjustment*u[j] }
            guard weights.allSatisfy(\.isFinite),abs(weights.reduce(0,+)-1)<1e-8 else { throw VivoOmicsError.invalid("target-kernel query weights") }
            var ordinary=[Double](repeating: 0,count: p),shuffled=ordinary
            for i in 0..<s { for j in 0..<p {
                ordinary[j]+=weights[i]*delta[model.supportedRows[i]][j]
                shuffled[j]+=weights[i]*delta[model.supportedRows[(i+1)%s]][j]
            } }
            for j in 0..<p { ordinary[j]+=baseline[j];shuffled[j]+=baseline[j] }
            let (expression,diagnostic)=try clipped(ordinary),(shuffledExpression,shuffledDiagnostic)=try clipped(shuffled)
            results.append(.init(id: query.id,targetID: query.descriptor.targetID,status: "predicted",maximumSimilarity: maximum,
                weights: weights,expression: expression,shuffledExpression: shuffledExpression,diagnostic: diagnostic,shuffledDiagnostic: shuffledDiagnostic))
        }
        return .init(schemaVersion: 1,context: model.base.context,featureIDs: model.base.featureIDs,supportedTrainingTargets: descriptors.map(\.targetID),baselines: baselines,queries: results,
            qualification: "Fixed term-Jaccard kernel ridge with unpenalized intercept, native Double arithmetic and independent query descriptors. Condition means and original feature denominators; no dense cell-by-gene matrix. Seen target IDs and known feature aliases reject. Missing or disjoint descriptors have no target-specific prediction. Clipping is reported without CPM reclosure. Generic baselines and a matched fixed shuffle are retained. Descriptor provenance is bound, not externally authenticated. Resident condition matrices; no unseen-study, context-transfer, uncertainty, single-cell distribution, causal, kinetic or Metal qualification.")
    }
}
