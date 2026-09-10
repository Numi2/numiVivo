import Foundation

public struct VivoPerturbationPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let featureNamespace: String
    public let perturbationID: String
    public let controlCondition: String
    public let treatmentCondition: String
    public let provenance: String
    public init(mapping: VivoH5ADImportPlan,featureNamespace: String,perturbationID: String,controlCondition: String,treatmentCondition: String,provenance: String) {
        schemaVersion=1;self.mapping=mapping;self.featureNamespace=featureNamespace;self.perturbationID=perturbationID
        self.controlCondition=controlCondition;self.treatmentCondition=treatmentCondition;self.provenance=provenance
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,mapping,featureNamespace,perturbationID,controlCondition,treatmentCondition,provenance }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","mapping","featureNamespace","perturbationID","controlCondition","treatmentCondition","provenance"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);mapping=try c.decode(VivoH5ADImportPlan.self,forKey: .mapping)
        featureNamespace=try c.decode(String.self,forKey: .featureNamespace);perturbationID=try c.decode(String.self,forKey: .perturbationID)
        controlCondition=try c.decode(String.self,forKey: .controlCondition);treatmentCondition=try c.decode(String.self,forKey: .treatmentCondition)
        provenance=try c.decode(String.self,forKey: .provenance)
    }
    func validate() throws {
        guard schemaVersion==1,[featureNamespace,perturbationID,controlCondition,treatmentCondition].allSatisfy(vivoOmicsID),
              controlCondition != treatmentCondition,!provenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,provenance.utf8.count<=16_384,
              mapping.samples.allSatisfy({ $0.condition==controlCondition || $0.condition==treatmentCondition }) else {
            throw VivoOmicsError.invalid("perturbation plan schema, identity, provenance or conditions")
        }
    }
    var trainingPlan: VivoH5ADPseudobulkPlan { .init(mapping: mapping) }
}
public struct VivoPerturbationQueryPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let featureNamespace: String
    public let perturbationID: String
    public init(mapping: VivoH5ADImportPlan,featureNamespace: String,perturbationID: String) {
        schemaVersion=1;self.mapping=mapping;self.featureNamespace=featureNamespace;self.perturbationID=perturbationID
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,mapping,featureNamespace,perturbationID }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","mapping","featureNamespace","perturbationID"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);mapping=try c.decode(VivoH5ADImportPlan.self,forKey: .mapping)
        featureNamespace=try c.decode(String.self,forKey: .featureNamespace);perturbationID=try c.decode(String.self,forKey: .perturbationID)
    }
    func validate() throws {
        guard schemaVersion==1,vivoOmicsID(featureNamespace),vivoOmicsID(perturbationID) else { throw VivoOmicsError.invalid("perturbation query schema or identity") }
    }
}
public struct VivoPerturbationModel: Codable, Sendable, Equatable {
    public let method: String
    public let plan: VivoPerturbationPlan
    public let source: VivoFingerprint
    public let featureIDs: [String]
    public let organism: String
    public let cellGroup: String?
    public let trainingDonors: [String]
    public let selectedFeatureIndices: [Int]
    public let contextCenters: [Double]
    public let contextScales: [Double]
    public let contexts: [[Double]]
    public let meanResponse: [Double]
    public let medianResponse: [Double]
    public let dualCoefficients: [[Double]]
    public let maximumSolveResidual: Double
    public let qualification: String
}
public struct VivoPerturbationEstimate: Codable, Sendable, Equatable {
    public let baseline: String
    public let unclippedResponse: [Double]
    public let predictedTreated: [Double]
    public let predictedResponse: [Double]
    public let impliedCPMSum: Double
}
public struct VivoPerturbationPrediction: Codable, Sendable, Equatable {
    public let group: VivoPseudobulkGroup
    public let libraryCounts: UInt64
    public let control: [Double]
    public let estimates: [VivoPerturbationEstimate]
}
public struct VivoPerturbationReport: Codable, Sendable, Equatable {
    public let method: String
    public let referenceSource: VivoFingerprint
    public let featureIDs: [String]
    public let predictions: [VivoPerturbationPrediction]
    public let qualification: String
}

public enum VivoPerturbation {
    /// Only donor-level pseudobulks are materialized. Source cells stay sparse.
    static func dense(_ bulk: VivoPseudobulkCounts) throws -> ([[UInt64]],[[Double]],[UInt64]) {
        let n=bulk.groups.count,m=bulk.featureIDs.count
        guard n>0,n<=128,m>0,n<=4_000_000/m else { throw VivoOmicsError.limit("perturbation pseudobulk storage budget") }
        var counts=Array(repeating: [UInt64](repeating: 0,count: m),count: n),logs=Array(repeating: [Double](repeating: 0,count: m),count: n),totals=[UInt64](repeating: 0,count: n)
        for i in 0..<n {
            try Task.checkCancellation()
            for k in bulk.matrix.rowOffsets[i]..<bulk.matrix.rowOffsets[i+1] {
                let j=bulk.matrix.featureIndices[k],v=bulk.matrix.counts[k];counts[i][j]=v;totals[i]=try vivoOmicsSum(totals[i],v)
            }
            guard totals[i]>0 else { throw VivoOmicsError.invalid("perturbation pseudobulk has an empty library") }
            for j in 0..<m { logs[i][j]=log1p(Double(counts[i][j])/Double(totals[i])*1_000_000) }
        }
        return (counts,logs,totals)
    }
    static func model(_ report: VivoH5ADPseudobulkReport,plan: VivoPerturbationPlan,source: VivoFingerprint) throws -> VivoPerturbationModel {
        try model(report.pseudobulk,plan: plan,source: source)
    }
    /// The caller validates and isolates aggregate membership before fitting.
    static func model(_ bulk: VivoPseudobulkCounts,plan: VivoPerturbationPlan,source: VivoFingerprint) throws -> VivoPerturbationModel {
        try plan.validate()
        let m=bulk.featureIDs.count,groups=bulk.groups
        let organisms=Set(groups.map(\.organism)),cellGroups=Set(groups.map(\.cellGroup))
        guard organisms.count==1,cellGroups.count==1,groups.allSatisfy({ $0.donorID != nil }),
              groups.allSatisfy({ $0.condition==plan.controlCondition || $0.condition==plan.treatmentCondition }) else {
            throw VivoOmicsError.invalid("perturbation training requires one organism/cell group and explicit paired donors")
        }
        let donors=Set(groups.compactMap(\.donorID)).sorted(),n=donors.count
        guard (2...64).contains(n),n<=2_000_000/m,n<=100_000_000/m/n else { throw VivoOmicsError.limit("perturbation donor or fit work budget") }
        var controls: [Int]=[],treated: [Int]=[]
        for donor in donors {
            let c=groups.indices.filter { groups[$0].donorID==donor && groups[$0].condition==plan.controlCondition }
            let t=groups.indices.filter { groups[$0].donorID==donor && groups[$0].condition==plan.treatmentCondition }
            guard c.count==1,t.count==1,groups[c[0]].biologicalReplicateID==groups[t[0]].biologicalReplicateID else {
                throw VivoOmicsError.invalid("perturbation training needs exactly one matched control/treatment pair per donor")
            }
            controls.append(c[0]);treated.append(t[0])
        }
        let (counts,logs,_)=try dense(bulk)
        var selected: [Int]=[],mean=[Double](repeating: 0,count: m),median=mean
        var response=Array(repeating: mean,count: n)
        for j in 0..<m {
            if j%256==0 { try Task.checkCancellation() }
            var sum: UInt64=0,expressing=0
            for i in 0..<n {
                let pair=try vivoOmicsSum(counts[controls[i]][j],counts[treated[i]][j]);sum=try vivoOmicsSum(sum,pair)
                if pair>0 { expressing+=1 }
                response[i][j]=logs[treated[i]][j]-logs[controls[i]][j];mean[j]+=response[i][j]
            }
            mean[j]/=Double(n)
            if sum>=10,expressing>=2 { selected.append(j) }
            let ordered=response.map { $0[j] }.sorted()
            median[j]=n%2==0 ? (ordered[n/2-1]+ordered[n/2])/2 : ordered[n/2]
        }
        guard !selected.isEmpty else { throw VivoOmicsError.invalid("perturbation context has no training-expressed features") }
        let width=selected.count,divisor=sqrt(Double(width))
        var centers=[Double](repeating: 0,count: width),scales=centers,contexts=Array(repeating: centers,count: n)
        for (k,j) in selected.enumerated() {
            for i in controls { centers[k]+=logs[i][j] };centers[k]/=Double(n)
            for i in controls { let delta=logs[i][j]-centers[k];scales[k]+=delta*delta }
            scales[k]=sqrt(scales[k]/Double(n))
            if controls.allSatisfy({ logs[$0][j]==logs[controls[0]][j] }) { centers[k]=logs[controls[0]][j];scales[k]=1 }
            else if scales[k]==0 { scales[k]=1 }
            for i in 0..<n { contexts[i][k]=(logs[controls[i]][j]-centers[k])/scales[k]/divisor }
        }
        var system=[Double](repeating: 0,count: n*n)
        for i in 0..<n { for j in 0..<n { system[i*n+j]=VivoSingleCellReduction.dot(contexts[i],contexts[j])+(i==j ? 1:0) } }
        let eigen=try VivoSingleCellReduction.symmetricEigen(system,n: n)
        guard eigen.values.allSatisfy({ $0.isFinite && $0>0 }) else { throw VivoOmicsError.invalid("perturbation ridge system is not positive definite") }
        var dual=Array(repeating: [Double](repeating: 0,count: m),count: n),maximumResidual=0.0
        for gene in 0..<m {
            if gene%256==0 { try Task.checkCancellation() }
            for component in 0..<n {
                var projection=0.0
                for i in 0..<n { projection+=eigen.vectors[i*n+component]*(response[i][gene]-mean[gene]) }
                projection/=eigen.values[component]
                for i in 0..<n { dual[i][gene]+=eigen.vectors[i*n+component]*projection }
            }
            for i in 0..<n {
                let target=response[i][gene]-mean[gene]
                var actual=0.0;for j in 0..<n { actual+=system[i*n+j]*dual[j][gene] }
                maximumResidual=max(maximumResidual,abs(actual-target)/max(1,abs(target)))
            }
        }
        guard maximumResidual.isFinite,maximumResidual<=1e-9 else { throw VivoOmicsError.invalid("perturbation ridge solve residual") }
        return .init(method: "paired-donor-log1p-CPM-response-baselines-alpha1-v1",plan: plan,source: source,featureIDs: bulk.featureIDs,
            organism: organisms.first!,cellGroup: groups.first!.cellGroup,trainingDonors: donors,selectedFeatureIndices: selected,
            contextCenters: centers,contextScales: scales,contexts: contexts,meanResponse: mean,medianResponse: median,dualCoefficients: dual,
            maximumSolveResidual: maximumResidual,qualification: "Known perturbation donor-response baselines, equal donor weighting; no Bayesian uncertainty, unseen-perturbation or mechanistic qualification.")
    }
    static func evaluate(snapshot: URL,plan: VivoPerturbationQueryPlan,model: VivoPerturbationModel) throws -> VivoPerturbationReport {
        try validateQuery(plan,model: model)
        let report=try VivoH5ADPseudobulk.evaluateSnapshot(snapshot,plan: .init(mapping: plan.mapping))
        return try evaluate(report.pseudobulk,plan: plan,model: model)
    }
    private static func validateQuery(_ plan: VivoPerturbationQueryPlan,model: VivoPerturbationModel) throws {
        try plan.validate()
        guard plan.featureNamespace==model.plan.featureNamespace,plan.perturbationID==model.plan.perturbationID,
              plan.mapping.countUnit==model.plan.mapping.countUnit else { throw VivoOmicsError.invalid("perturbation query identity, namespace or unit mismatch") }
        guard plan.mapping.samples.allSatisfy({ $0.condition==model.plan.controlCondition }) else {
            throw VivoOmicsError.invalid("perturbation prediction accepts control-only query samples")
        }
    }
    static func evaluate(_ bulk: VivoPseudobulkCounts,plan: VivoPerturbationQueryPlan,model: VivoPerturbationModel) throws -> VivoPerturbationReport {
        try validateQuery(plan,model: model)
        guard Set(bulk.featureIDs)==Set(model.featureIDs),bulk.featureIDs.count==model.featureIDs.count,
              bulk.groups.allSatisfy({ $0.organism==model.organism && $0.cellGroup==model.cellGroup && $0.donorID != nil && $0.condition==model.plan.controlCondition }),
              Set(bulk.groups.compactMap(\.donorID)).count==bulk.groups.count else {
            throw VivoOmicsError.invalid("perturbation query feature universe, organism, cell group or donor multiplicity mismatch")
        }
        guard Set(bulk.groups.compactMap(\.donorID)).isDisjoint(with: model.trainingDonors) else { throw VivoOmicsError.invalid("perturbation query donor overlaps training") }
        let (_,logs,totals)=try dense(bulk),m=model.featureIDs.count,n=model.trainingDonors.count
        guard bulk.groups.count<=20_000_000/m else { throw VivoOmicsError.limit("perturbation prediction output budget") }
        let lookup=Dictionary(uniqueKeysWithValues: bulk.featureIDs.enumerated().map { ($0.element,$0.offset) })
        let order=model.featureIDs.map { lookup[$0]! },selected=model.selectedFeatureIndices
        var predictions: [VivoPerturbationPrediction]=[]
        for row in bulk.groups.indices {
            try Task.checkCancellation()
            let control=order.map { logs[row][$0] }
            let context=selected.indices.map { (control[selected[$0]]-model.contextCenters[$0])/model.contextScales[$0]/sqrt(Double(selected.count)) }
            let similarities=model.contexts.map { VivoSingleCellReduction.dot(context,$0) }
            var ridge=model.meanResponse
            for i in 0..<n { for gene in 0..<m { ridge[gene]+=similarities[i]*model.dualCoefficients[i][gene] } }
            var estimates: [VivoPerturbationEstimate]=[]
            for (name,response) in [("noChange",[Double](repeating: 0,count: m)),("meanResponse",model.meanResponse),("medianResponse",model.medianResponse),("contextRidge",ridge)] {
                let predicted=(0..<m).map { max(0,control[$0]+response[$0]) },applied=(0..<m).map { predicted[$0]-control[$0] }
                let implied=predicted.reduce(0.0) { $0+expm1($1) }
                guard response.allSatisfy(\.isFinite),predicted.allSatisfy(\.isFinite),implied.isFinite else { throw VivoOmicsError.invalid("perturbation prediction is nonfinite") }
                estimates.append(.init(baseline: name,unclippedResponse: response,predictedTreated: predicted,predictedResponse: applied,impliedCPMSum: implied))
            }
            predictions.append(.init(group: bulk.groups[row],libraryCounts: totals[row],control: control,estimates: estimates))
        }
        return .init(method: model.method,referenceSource: model.source,featureIDs: model.featureIDs,predictions: predictions,
            qualification: "Frozen known-perturbation response estimates for supplied held-out control donors. Log1p-CPM point estimates are not reclosed compositions or raw count libraries; implied CPM sums are reported. No treated query outcomes used. No calibrated intervals, single-cell distributions, unseen perturbation identity or causal/mechanistic claims.")
    }
}
