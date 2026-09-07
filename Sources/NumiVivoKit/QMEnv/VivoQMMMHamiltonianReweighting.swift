import Foundation

public struct VivoQMMMHamiltonianEnergyObservation:Codable,Sendable,Equatable {
    public let sourceConfigurationFingerprint:VivoFingerprint
    public let referenceEnergyKJPerMol:Double
    public let candidateEnergyKJPerMol:Double
    public init(sourceConfigurationFingerprint:VivoFingerprint,referenceEnergyKJPerMol:Double,candidateEnergyKJPerMol:Double) {
        self.sourceConfigurationFingerprint=sourceConfigurationFingerprint
        self.referenceEnergyKJPerMol=referenceEnergyKJPerMol;self.candidateEnergyKJPerMol=candidateEnergyKJPerMol
    }
}

public struct VivoQMMMHamiltonianReweightingRequest:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-hamiltonian-reweighting/v1"
    public var schema:String
    public var referenceHamiltonianFingerprint:VivoFingerprint
    public var candidateHamiltonianFingerprint:VivoFingerprint
    public var samplingExecutionFingerprint:VivoFingerprint
    public var temperatureK:Double
    public var observations:[VivoQMMMHamiltonianEnergyObservation]
    public var minimumEffectiveSamples:Double
    public var maximumNormalizedWeight:Double
    public init(referenceHamiltonianFingerprint:VivoFingerprint,candidateHamiltonianFingerprint:VivoFingerprint,
                samplingExecutionFingerprint:VivoFingerprint,temperatureK:Double,
                observations:[VivoQMMMHamiltonianEnergyObservation],minimumEffectiveSamples:Double=100,
                maximumNormalizedWeight:Double=0.05) {
        schema=Self.schema;self.referenceHamiltonianFingerprint=referenceHamiltonianFingerprint
        self.candidateHamiltonianFingerprint=candidateHamiltonianFingerprint;self.samplingExecutionFingerprint=samplingExecutionFingerprint
        self.temperatureK=temperatureK;self.observations=observations;self.minimumEffectiveSamples=minimumEffectiveSamples
        self.maximumNormalizedWeight=maximumNormalizedWeight
    }
}

public struct VivoQMMMHamiltonianReweightingResult:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-hamiltonian-reweighting-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let freeEnergyDifferenceKJPerMol:Double
    public let effectiveSamples:Double
    public let maximumNormalizedWeight:Double
    public let meanEnergyDifferenceKJPerMol:Double
    public let energyDifferenceStandardDeviationKJPerMol:Double
    public let converged:Bool
    public let issues:[String]
    public let interpretation:String
    public let evidenceFingerprint:VivoFingerprint
}

/// One-direction Zwanzig reweighting of a reference equilibrium ensemble into a
/// candidate Hamiltonian. The result is accepted only when normalized importance
/// weights retain enough effective samples and no single configuration dominates.
/// Low overlap is a mandatory fresh-sampling signal, not permission to extrapolate.
public enum VivoQMMMHamiltonianReweighting {
    private static let gasConstantKJ=0.00831446261815324
    public static let interpretation="One-direction equilibrium Hamiltonian free-energy perturbation over explicitly retained configurations. Effective-sample and maximum-weight gates determine whether the reference ensemble supports the candidate Hamiltonian; failure requires new sampling and is not converted into a correction."

    public static func calculate(_ request:VivoQMMMHamiltonianReweightingRequest)throws->VivoQMMMHamiltonianReweightingResult {
        guard request.schema==VivoQMMMHamiltonianReweightingRequest.schema,
              request.referenceHamiltonianFingerprint != request.candidateHamiltonianFingerprint,
              request.temperatureK.isFinite,request.temperatureK>0,
              request.observations.count>=2,request.observations.count<=10_000_000,
              Set(request.observations.map(\.sourceConfigurationFingerprint)).count==request.observations.count,
              request.minimumEffectiveSamples.isFinite,request.minimumEffectiveSamples>=2,
              request.maximumNormalizedWeight.isFinite,request.maximumNormalizedWeight>0,request.maximumNormalizedWeight<1,
              request.observations.allSatisfy({$0.referenceEnergyKJPerMol.isFinite && $0.candidateEnergyKJPerMol.isFinite}) else {
            throw VivoChemistryError.invalid("QM/MM Hamiltonian reweighting identity, observations or acceptance")
        }
        let beta=1/(gasConstantKJ*request.temperatureK)
        let delta=request.observations.map{$0.candidateEnergyKJPerMol-$0.referenceEnergyKJPerMol}
        let logWeights=delta.map{-beta*$0}
        guard let maximum=logWeights.max(),maximum.isFinite else { throw VivoChemistryError.convergence("Hamiltonian reweighting log weights") }
        let shifted=logWeights.map{exp($0-maximum)},sum=shifted.reduce(0,+),square=shifted.reduce(0){$0+$1*$1}
        guard sum.isFinite,sum>0,square.isFinite,square>0 else { throw VivoChemistryError.convergence("Hamiltonian reweighting normalization") }
        let normalized=shifted.map{$0/sum},effective=sum*sum/square,largest=normalized.max() ?? 1
        let logMean=maximum+log(sum)-log(Double(delta.count))
        let freeEnergy = -logMean/beta
        let mean=delta.reduce(0,+)/Double(delta.count)
        let variance=delta.count>1 ? delta.reduce(0.0){$0+pow($1-mean,2)}/Double(delta.count-1):0
        let sd=sqrt(max(0,variance))
        guard freeEnergy.isFinite,effective.isFinite,largest.isFinite,mean.isFinite,sd.isFinite else {
            throw VivoChemistryError.convergence("Hamiltonian reweighting estimator")
        }
        var issues:[String]=[]
        if effective<request.minimumEffectiveSamples { issues.append("candidate Hamiltonian effective reweighting sample count is below threshold") }
        if largest>request.maximumNormalizedWeight { issues.append("candidate Hamiltonian is dominated by an individual reference configuration") }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence:Codable { let schema:String;let request:VivoQMMMHamiltonianReweightingRequest;let delta:[Double];let freeEnergy:Double;let effective:Double;let largest:Double;let issues:[String] }
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema:"numivivo.org/qmmm-hamiltonian-reweighting-evidence/v1",request:request,delta:delta,
            freeEnergy:freeEnergy,effective:effective,largest:largest,issues:issues)))
        return .init(schema:VivoQMMMHamiltonianReweightingResult.schema,requestFingerprint:requestID,
            freeEnergyDifferenceKJPerMol:freeEnergy,effectiveSamples:effective,maximumNormalizedWeight:largest,
            meanEnergyDifferenceKJPerMol:mean,energyDifferenceStandardDeviationKJPerMol:sd,
            converged:issues.isEmpty,issues:issues,interpretation:interpretation,evidenceFingerprint:evidenceID)
    }

    public static func validate(_ result:VivoQMMMHamiltonianReweightingResult,
                                request:VivoQMMMHamiltonianReweightingRequest)throws {
        guard result==(try calculate(request)) else { throw VivoChemistryError.invalid("Hamiltonian reweighting does not reconstruct") }
    }
}
