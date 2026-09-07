import Foundation
import NumiVivoKit

struct VivoQMMMFreeEnergyCLICommands {
    static func handles(_ name:String?)->Bool {
        ["qmmm-free-energy-analyze","qmmm-free-energy-rate","qmmm-free-energy-replicated-rate",
         "qmmm-transmission-analyze","qmmm-transmission-apply","qmmm-chemical-qualify",
         "qmmm-chemical-state-network","qmmm-chemical-exchange-network","qmmm-chemical-exchange-validate",
         "qmmm-free-energy-help"].contains(name ?? "")
    }
    private func load<T:Decodable & Sendable>(_ type:T.Type,_ path:String)throws->T {
        try VivoValidatedArtifactLoader.decode(type,at:URL(fileURLWithPath:path),limits:.init(
            maximumBytes:512*1024*1024,maximumNodes:12_000_000,maximumArrayElements:4_000_000)).value
    }
    private func write<T:Encodable>(_ value:T,_ path:String?)throws {
        let data=try VivoCanonicalJSON.encode(value)
        guard let path else { FileHandle.standardOutput.write(data);FileHandle.standardOutput.write(Data("\n".utf8));return }
        try VivoKineticsDocumentIO.write(data,to:URL(fileURLWithPath:path),overwrite:false)
    }
    func run(arguments:[String]) async -> Int32 {
        do {
            guard let command=arguments.first,Self.handles(command) else { throw VivoChemistryError.invalid("QM/MM free-energy command") }
            if command=="qmmm-free-energy-help" {
                guard arguments.count==1 else { throw VivoChemistryError.invalid("qmmm-free-energy-help has no options") }
                FileHandle.standardOutput.write(Data(Self.help.utf8));return 0
            }
            guard arguments.count>=2,!arguments[1].hasPrefix("--") else { throw VivoChemistryError.invalid("QM/MM free-energy input JSON is required") }
            var output:String?,kinetics:String?,i=2
            while i<arguments.count {
                guard i+1<arguments.count else { throw VivoChemistryError.invalid("QM/MM free-energy option requires a value") }
                switch arguments[i] {
                case "--output":guard output==nil else {throw VivoChemistryError.invalid("duplicate --output")};output=arguments[i+1]
                case "--kinetics":guard kinetics==nil else {throw VivoChemistryError.invalid("duplicate --kinetics")};kinetics=arguments[i+1]
                default:throw VivoChemistryError.invalid("unknown QM/MM free-energy option")
                }
                i+=2
            }
            switch command {
            case "qmmm-free-energy-analyze":
                guard kinetics==nil else { throw VivoChemistryError.invalid("--kinetics applies only to QM/MM rate commands") }
                let request=try load(VivoQMMMFreeEnergyAnalysisRequest.self,arguments[1]),result=try request.calculate()
                try write(result,output);return result.converged ? 0:75
            case "qmmm-transmission-analyze":
                guard kinetics==nil else { throw VivoChemistryError.invalid("--kinetics does not apply to transmission evidence analysis") }
                let request=try load(VivoQMMMDynamicalTransmissionAnalysisRequest.self,arguments[1]),result=try request.calculate()
                try write(result,output);return result.converged ? 0:75
            case "qmmm-transmission-apply":
                guard kinetics==nil else { throw VivoChemistryError.invalid("--kinetics applies after flux-normalized rate calculation, not transmission application") }
                let request=try load(VivoQMMMComputedTransmissionApplicationRequest.self,arguments[1])
                try write(request.calculate(),output);return 0
            case "qmmm-free-energy-rate":
                let request=try load(VivoQMMMFreeEnergyRateRequest.self,arguments[1]),result=try VivoQMMMFreeEnergyRate.calculate(request)
                try write(result,output)
                if let kinetics {
                    guard let output else { throw VivoChemistryError.invalid("--kinetics requires --output so the updated pack has a distinct path") }
                    let pack=try load(VivoCovalentKineticPack.self,kinetics)
                    try write(VivoQMMMFreeEnergyRate.applying(result,request:request,to:pack),output+".kinetics.json")
                }
                return 0
            case "qmmm-free-energy-replicated-rate":
                let request=try load(VivoQMMMReplicatedFreeEnergyRateRequest.self,arguments[1])
                let result=try VivoQMMMReplicatedFreeEnergyRate.calculate(request)
                try write(result,output)
                if let kinetics {
                    guard let output else { throw VivoChemistryError.invalid("--kinetics requires --output so the updated pack has a distinct path") }
                    let pack=try load(VivoCovalentKineticPack.self,kinetics)
                    try write(VivoQMMMReplicatedFreeEnergyRate.applying(result,request:request,to:pack),output+".kinetics.json")
                }
                return result.converged ? 0:75
            case "qmmm-chemical-qualify":
                guard kinetics==nil else { throw VivoChemistryError.invalid("--kinetics is not an external-validation operation") }
                let request=try load(VivoQMMMChemicalQualificationRequest.self,arguments[1])
                let result=try VivoQMMMChemicalQualification.calculate(request)
                try write(result,output);return result.converged ? 0:75
            case "qmmm-chemical-state-network":
                guard kinetics==nil else { throw VivoChemistryError.invalid("--kinetics requires an explicit downstream kinetic model adapter") }
                let request=try load(VivoQMMMChemicalStateNetworkRequest.self,arguments[1])
                try write(VivoQMMMChemicalStateNetwork.calculate(request),output);return 0
            case "qmmm-chemical-exchange-network":
                guard kinetics==nil else { throw VivoChemistryError.invalid("--kinetics requires an explicit downstream kinetic model adapter") }
                let request=try load(VivoQMMMChemicalExchangeNetworkRequest.self,arguments[1])
                try write(VivoQMMMChemicalExchangeNetwork.calculate(request),output);return 0
            case "qmmm-chemical-exchange-validate":
                guard kinetics==nil else { throw VivoChemistryError.invalid("--kinetics is not a transient external-validation operation") }
                let request=try load(VivoQMMMChemicalExchangeValidationRequest.self,arguments[1])
                let result=try VivoQMMMChemicalExchangeValidation.calculate(request)
                try write(result,output);return result.converged ? 0:75
            default:throw VivoChemistryError.invalid("unresolved QM/MM free-energy command")
            }
        } catch {
            FileHandle.standardError.write(Data("QM/MM free-energy workflow failed: \(error)\n".utf8));return 65
        }
    }
    private static let help="""
    QM/MM activation free-energy, transmission and chemical qualification workflows
      numivivo qmmm-free-energy-analyze analysis-request.json --output activation-pmf.json
      numivivo qmmm-transmission-analyze transmission-evidence.json --output transmission.json
      numivivo qmmm-transmission-apply transmission-application.json --output calculated-rate-request.json
      numivivo qmmm-free-energy-rate rate-request.json --output rate.json [--kinetics kinetic-pack.json]
      numivivo qmmm-free-energy-replicated-rate replicas.json --output replicated-rate.json [--kinetics kinetic-pack.json]
      numivivo qmmm-chemical-qualify qualification.json --output qualification-result.json
      numivivo qmmm-chemical-state-network state-network.json --output effective-rate.json
      numivivo qmmm-chemical-exchange-network exchange-network.json --output transient-kinetics.json
      numivivo qmmm-chemical-exchange-validate validation.json --output transient-validation.json

    qmmm-free-energy-analyze reconstructs retained umbrella data with native unbinned
    MBAR and exits 75 if sampling/overlap gates fail. qmmm-transmission-analyze
    reconstructs paired +v/-v signed reactive-flux histories and requires a stable
    late-time plateau; it does not serialize a Born-Oppenheimer provider.
    qmmm-transmission-apply binds a converged transmission coefficient and sampled
    dividing-surface velocity normalization to the exact PMF. Replicated-rate
    qualification requires disjoint stochastic executions and preserves unresolved
    within-replica uncertainty as unknown. qmmm-chemical-qualify evaluates
    predeclared protocol/electronic/QM-region sensitivity and only directly compares
    condition-matched first-order chemical-rate measurements. qmmm-chemical-state-network
    computes sum_s p_s sum_path k_s,path only when rapid pre-equilibrium is explicitly
    asserted. qmmm-chemical-exchange-network propagates explicit transient state exchange
    when interconversion competes with chemical conversion. qmmm-chemical-exchange-validate
    compares that potentially non-single-exponential survival curve only to measured,
    condition-matched survival fractions at explicit simulated times. Outputs are
    canonical no-clobber JSON.
    """
}
