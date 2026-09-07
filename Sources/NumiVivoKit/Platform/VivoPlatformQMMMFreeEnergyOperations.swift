import Foundation

/// Deterministic artifact-DAG adapters for retained PMF analysis and rate
/// qualification. Executable force providers remain Swift runtime resources and
/// are intentionally not serialized through workflow artifacts.
public enum VivoPlatformQMMMFreeEnergyOperations {
    public static func definitions(implementationFingerprint id:VivoFingerprint)->[VivoWorkflowDefinition] {
        [VivoPlatformOperations.pure(identifier:"vivo.platform.qmmm-free-energy-analyze",id:id,
            inputs:["request":"vivo.qmmm-free-energy-analysis-request"],
            outputs:[.init(name:"result",kind:"vivo.qmmm-activation-free-energy-result")],
            summary:"Reconstruct retained umbrella traces with native unbinned MBAR and numerical acceptance diagnostics.",
            configure:VivoPlatformOperations.empty,calculate:{ _,inputs,_ in
                let request=try VivoPlatformOperations.input(VivoQMMMFreeEnergyAnalysisRequest.self,"request",inputs)
                return ["result":try VivoCanonicalJSON.encode(request.calculate())]
            }),
         VivoPlatformOperations.pure(identifier:"vivo.platform.qmmm-free-energy-rate",id:id,
            inputs:["request":"vivo.qmmm-free-energy-rate-request"],
            outputs:[.init(name:"result",kind:"vivo.qmmm-free-energy-rate-result")],
            summary:"Flux-normalized conventional TST rate from one context-qualified converged QM/MM PMF.",
            configure:VivoPlatformOperations.empty,calculate:{ _,inputs,_ in
                let request=try VivoPlatformOperations.input(VivoQMMMFreeEnergyRateRequest.self,"request",inputs)
                return ["result":try VivoCanonicalJSON.encode(VivoQMMMFreeEnergyRate.calculate(request))]
            }),
         VivoPlatformOperations.pure(identifier:"vivo.platform.qmmm-transmission-analyze",id:id,
            inputs:["request":"vivo.qmmm-dynamical-transmission-analysis-request"],
            outputs:[.init(name:"result",kind:"vivo.qmmm-dynamical-transmission-result")],
            summary:"Reconstruct flux-weighted dividing-surface recrossing evidence already generated under the bound BO Hamiltonian.",
            configure:VivoPlatformOperations.empty,calculate:{ _,inputs,_ in
                let request=try VivoPlatformOperations.input(VivoQMMMDynamicalTransmissionAnalysisRequest.self,"request",inputs)
                return ["result":try VivoCanonicalJSON.encode(request.calculate())]
            }),
         VivoPlatformOperations.pure(identifier:"vivo.platform.qmmm-apply-transmission",id:id,
            inputs:["request":"vivo.qmmm-computed-transmission-application-request"],
            outputs:[.init(name:"rateRequest",kind:"vivo.qmmm-free-energy-rate-request")],
            summary:"Apply only a converged computed recrossing coefficient to its exact PMF rate request.",
            configure:VivoPlatformOperations.empty,calculate:{ _,inputs,_ in
                let request=try VivoPlatformOperations.input(VivoQMMMComputedTransmissionApplicationRequest.self,"request",inputs)
                return ["rateRequest":try VivoCanonicalJSON.encode(request.calculate())]
            }),
         VivoPlatformOperations.pure(identifier:"vivo.platform.qmmm-replicated-free-energy-rate",id:id,
            inputs:["request":"vivo.qmmm-replicated-free-energy-rate-request"],
            outputs:[.init(name:"result",kind:"vivo.qmmm-replicated-free-energy-rate-result")],
            summary:"Independent disjoint-seed PMF replica agreement and replicated kinetic evidence.",
            configure:VivoPlatformOperations.empty,calculate:{ _,inputs,_ in
                let request=try VivoPlatformOperations.input(VivoQMMMReplicatedFreeEnergyRateRequest.self,"request",inputs)
                return ["result":try VivoCanonicalJSON.encode(VivoQMMMReplicatedFreeEnergyRate.calculate(request))]
            }),
         VivoPlatformOperations.pure(identifier:"vivo.platform.qmmm-apply-replicated-rate",id:id,
            inputs:["request":"vivo.qmmm-replicated-free-energy-rate-request",
                    "result":"vivo.qmmm-replicated-free-energy-rate-result",
                    "kinetics":"vivo.covalent-kinetic-pack"],
            outputs:[.init(name:"kinetics",kind:"vivo.covalent-kinetic-pack")],
            summary:"Replace only bound-complex inactivation after replicated QM/MM rate agreement passes.",
            configure:VivoPlatformOperations.empty,calculate:{ _,inputs,_ in
                let request=try VivoPlatformOperations.input(VivoQMMMReplicatedFreeEnergyRateRequest.self,"request",inputs)
                let result=try VivoPlatformOperations.input(VivoQMMMReplicatedFreeEnergyRateResult.self,"result",inputs)
                let kinetics=try VivoPlatformOperations.input(VivoCovalentKineticPack.self,"kinetics",inputs)
                return ["kinetics":try VivoCanonicalJSON.encode(VivoQMMMReplicatedFreeEnergyRate.applying(result,request:request,to:kinetics))]
            })]
    }
}
