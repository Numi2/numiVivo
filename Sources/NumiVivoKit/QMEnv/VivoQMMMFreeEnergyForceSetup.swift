import Foundation
@preconcurrency import Metal

/// One assembly surface for free-energy sampling. Every mode returns the exact
/// retained classical system integrated by VivoMDMetalRuntime plus an additive
/// Born-Oppenheimer provider bound to that system. The umbrella layer therefore
/// remains independent of whether the electronic region is fixed or adaptive.
public enum VivoQMMMFreeEnergyElectronicMode: Sendable {
    case fixed(plan:VivoQMMMHamiltonianPlan,configuration:VivoQMMMDynamicsElectronicConfiguration)
    case adaptive(configuration:VivoAdaptiveQMMMConfiguration,basis:VivoAdaptiveBasisSpecification,
                  electronic:VivoQMMMDynamicsElectronicConfiguration)
}

public struct VivoQMMMFreeEnergyForceSetup: Sendable {
    public let retainedSystem:VivoClassicalSystem
    public let provider:VivoMDCandidateForceProvider
    public let modeDescription:String
    public init(retainedSystem:VivoClassicalSystem,provider:VivoMDCandidateForceProvider,modeDescription:String) {
        self.retainedSystem=retainedSystem;self.provider=provider;self.modeDescription=modeDescription
    }
}

public enum VivoQMMMFreeEnergyForceFactory {
    public static func make(document:VivoMolecularStructureDocument,
                            sourceSystem:VivoClassicalSystem,
                            initialGeometry:VivoMDCandidateGeometry,
                            dynamics:VivoMDConfiguration,
                            mode:VivoQMMMFreeEnergyElectronicMode,
                            device requested:MTLDevice? = nil) async throws -> VivoQMMMFreeEnergyForceSetup {
        try dynamics.validate();try VivoClassicalSystemValidator.validate(sourceSystem)
        switch mode {
        case .fixed(let plan,let electronic):
            try plan.validate(document:document,source:sourceSystem,budget:electronic.budget)
            let provider:VivoMDCandidateForceProvider
            if electronic.periodic.reciprocalMesh == nil {
                provider=try .hartreeFock(document:document,sourceSystem:sourceSystem,plan:plan,configuration:electronic)
            } else {
                provider=try await .hartreeFockPME(document:document,sourceSystem:sourceSystem,plan:plan,
                    configuration:electronic,device:requested)
            }
            let retained=plan.retainedSystem
            let retainedID=try retained.fingerprint()
            guard provider.retainedSystemFingerprint==retainedID else {
                throw VivoChemistryError.invalid("fixed QM/MM free-energy provider does not own the plan retained system")
            }
            switch plan.configuration.boundary {
            case .finiteCluster:
                guard initialGeometry.periodicCell==nil else {
                    throw VivoChemistryError.invalid("finite-cluster QM/MM free-energy setup received a periodic geometry")
                }
            case .periodicElectrostatic:
                guard initialGeometry.periodicCell != nil,dynamics.electrostatics == .pme else {
                    throw VivoChemistryError.invalid("periodic QM/MM free-energy setup requires a periodic geometry and retained PME")
                }
            }
            try provider.validate(system:retained,configuration:dynamics,cell:initialGeometry.periodicCell)
            let method=electronic.periodic.reciprocalMesh == nil ? "fixed-region analytic HF QM/MM" : "fixed-region analytic HF QM/MM with reciprocal Metal PME"
            return .init(retainedSystem:retained,provider:provider,modeDescription:method)

        case .adaptive(let adaptive,let basis,let electronic):
            guard initialGeometry.periodicCell != nil else {
                throw VivoChemistryError.invalid("adaptive QM/MM free-energy setup requires a periodic geometry")
            }
            let setup=try await VivoAdaptiveQMMMDynamics.make(document:document,system:sourceSystem,
                initialGeometry:initialGeometry,dynamics:dynamics,adaptive:adaptive,basis:basis,
                electronic:electronic,device:requested)
            let retainedID=try setup.retainedBaselineSystem.fingerprint()
            guard setup.provider.retainedSystemFingerprint==retainedID else {
                throw VivoChemistryError.invalid("adaptive QM/MM free-energy provider does not own its retained baseline system")
            }
            try setup.provider.validate(system:setup.retainedBaselineSystem,configuration:dynamics,cell:initialGeometry.periodicCell)
            return .init(retainedSystem:setup.retainedBaselineSystem,provider:setup.provider,
                modeDescription:"adaptive complete-partition analytic HF QM/MM with C2 transition forces")
        }
    }
}
