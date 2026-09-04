import CryptoKit
import Foundation
import Metal
import NumiVivoKit

private struct ExampleSummary: Encodable {
    let partitionModelFingerprint: String
    let partitionStep: VivoPartitionStepCertificate
    let finalConcentrations: [Float]
    let checkpointFingerprint: String
    let checkpointBytes: Int
    let campaignFingerprint: String
    let candidateCount: Int
    let jobCount: Int
    let hybridPlanFingerprint: String
    let hybridModes: [String]
    let surrogateAccepted: Bool
    let surrogateReason: String
    let couplingChannelCount: Int
}

@main
struct DigitalTissueHomeostasis {
    static func main() async throws {
        guard let device = MTLCreateSystemDefaultDevice() else {
            throw ExampleError.metalUnavailable
        }

        // One abstract mediator is distributed between three computational
        // physiological compartments. These values are model inputs, not a
        // physical dosing or experimental protocol.
        let partitionModel = try VivoPhysiologicalPartitionModel(
            name: "abstract-mediator-distribution",
            compartments: [
                .init(id: "central", volumeCubicMetres: 3.0e-3),
                .init(id: "interstitial", volumeCubicMetres: 9.0e-3),
                .init(id: "local-tissue", volumeCubicMetres: 2.0e-4)
            ],
            analytes: [
                .init(
                    id: "repair-mediator",
                    unit: "normalized",
                    minimumConcentration: 0,
                    maximumConcentration: 10
                )
            ],
            edges: [
                .init(
                    id: "central-to-interstitial",
                    analyte: "repair-mediator",
                    sourceCompartment: "central",
                    targetCompartment: "interstitial",
                    partitionCoefficient: 0.8,
                    clearanceCubicMetresPerSecond: 1.0e-5,
                    evidenceClass: "hypothetical"
                ),
                .init(
                    id: "interstitial-to-local",
                    analyte: "repair-mediator",
                    sourceCompartment: "interstitial",
                    targetCompartment: "local-tissue",
                    partitionCoefficient: 1.4,
                    clearanceCubicMetresPerSecond: 2.0e-6,
                    evidenceClass: "hypothetical"
                )
            ]
        )
        let partitionRuntime = try await VivoPhysiologicalPartitionRuntime.make(
            model: partitionModel,
            configuration: .init(
                timeStep: 0.25,
                minimumTimeStep: 1e-6,
                maximumTimeStep: 1,
                maximumAttempts: 16
            ),
            initialConcentrations: [1, 0, 0],
            device: device
        )
        let partitionStep = try await partitionRuntime.step(
            certifyAmountConservation: true
        )
        let partitionSnapshot = try await partitionRuntime.snapshot()

        let stateSection = VivoCheckpointSectionPayload.float32(
            id: "partition.concentrations",
            values: partitionSnapshot.concentrations
        )
        let modelSection = try VivoCheckpointSectionPayload.canonicalJSON(
            id: "partition.model",
            value: partitionModel
        )
        let checkpoint = try VivoCheckpointCodec.encode(.init(
            runtime: .physiologicalPartition,
            artifactFingerprint: partitionModel.fingerprint,
            stepIndex: UInt64(partitionSnapshot.stepIndex),
            logicalTime: Double(partitionSnapshot.absoluteTime),
            sections: [modelSection, stateSection],
            metadata: [
                "purpose": "computational-reference",
                "device": device.name
            ]
        ))
        let decodedCheckpoint = try VivoCheckpointCodec.decode(checkpoint)

        let artifacts = VivoCampaignArtifactReferences(
            programPackFingerprint: fingerprint("program-pack"),
            mechanismPackFingerprint: fingerprint("mechanism-pack"),
            hostContextFingerprint: fingerprint("host-context"),
            couplingPackFingerprint: fingerprint("coupling-pack")
        )
        let campaign = try VivoCampaignCompiler().compile(.init(
            name: "homeostasis-parameter-screen",
            artifacts: artifacts,
            sampling: .halton,
            axes: [
                .init(
                    parameterID: "response-gain",
                    unit: "1",
                    lowerBound: 0.25,
                    upperBound: 2.0,
                    sampleCount: 16
                ),
                .init(
                    parameterID: "clearance-rate",
                    unit: "1/s",
                    transform: .logarithmic10,
                    lowerBound: 1e-5,
                    upperBound: 1e-2,
                    sampleCount: 16
                )
            ],
            candidateCount: 16,
            replicateCount: 4,
            baseSeed: 0x4e_56_49_56_4f,
            maximumJobs: 10_000,
            tags: ["scope": "sequence-free-computational-reference"]
        ))

        let hybridPlanner = try VivoHybridStochasticPlanner()
        let hybridPlan = try hybridPlanner.plan(
            speciesCount: 4,
            reactions: [
                .init(
                    reactionIndex: 0,
                    speciesIndices: [0, 1],
                    critical: true,
                    minimumReactantCount: 12,
                    medianReactantCount: 28,
                    expectedFiringsPerStep: 4,
                    maximumPropensityPerSecond: 20
                ),
                .init(
                    reactionIndex: 1,
                    speciesIndices: [1, 2],
                    critical: false,
                    minimumReactantCount: 20_000,
                    medianReactantCount: 30_000,
                    expectedFiringsPerStep: 4_000,
                    maximumPropensityPerSecond: 200
                ),
                .init(
                    reactionIndex: 2,
                    speciesIndices: [3],
                    critical: false,
                    spatial: true,
                    minimumReactantCount: 2_000,
                    medianReactantCount: 5_000,
                    expectedFiringsPerStep: 500,
                    maximumPropensityPerSecond: 40
                )
            ]
        )

        let surrogateModelFingerprint = fingerprint("affine-surrogate")
        let surrogateContract = try VivoSurrogateContract(
            id: "bounded-homeostasis-surrogate",
            modelFingerprint: surrogateModelFingerprint,
            trainingDataFingerprint: fingerprint("training-data"),
            mechanismPackFingerprint: artifacts.mechanismPackFingerprint,
            hostContextFingerprint: artifacts.hostContextFingerprint,
            inputs: [
                .init(id: "signal", unit: "1", minimum: 0, maximum: 1),
                .init(id: "stress", unit: "1", minimum: 0, maximum: 1)
            ],
            outputs: [
                .init(id: "predicted-response", unit: "1", minimum: 0, maximum: 1)
            ],
            uncertaintyKind: .conformalRadius,
            maximumNormalizedUncertainty: 0.05,
            maximumNormalizedExtrapolation: 0,
            maximumConsecutiveAcceptedSteps: 16,
            mandatoryAuthorityInterval: 8
        )
        let surrogateBackend = try VivoAffineSurrogateBackend(
            modelFingerprint: surrogateModelFingerprint,
            inputWidth: 2,
            outputWidth: 1,
            weights: [0.65, -0.25],
            bias: [0.2],
            outputUncertainty: [0.02]
        )
        let surrogateGate = try VivoSurrogateAuthorityGate(
            contract: surrogateContract,
            backend: surrogateBackend
        )
        let surrogateDecision = await surrogateGate.evaluate(try .init(
            batchSize: 2,
            inputWidth: 2,
            values: [0.5, 0.1, 0.8, 0.3]
        ))

        let couplingChannels = try VivoNumiLabStandardChannels.make(shape: .init(
            spatialSampleCount: 64,
            phenotypeCount: 3,
            endocrineChannelCount: 4,
            neuralChannelCount: 8
        ))

        let summary = ExampleSummary(
            partitionModelFingerprint: partitionModel.fingerprint,
            partitionStep: partitionStep,
            finalConcentrations: partitionSnapshot.concentrations,
            checkpointFingerprint: decodedCheckpoint.packageFingerprint,
            checkpointBytes: checkpoint.count,
            campaignFingerprint: campaign.fingerprint,
            candidateCount: campaign.candidates.count,
            jobCount: campaign.jobs.count,
            hybridPlanFingerprint: hybridPlan.fingerprint,
            hybridModes: hybridPlan.cohorts.map { $0.mode.rawValue },
            surrogateAccepted: surrogateDecision.accepted,
            surrogateReason: surrogateDecision.reason.rawValue,
            couplingChannelCount: couplingChannels.count
        )
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        FileHandle.standardOutput.write(try encoder.encode(summary))
        FileHandle.standardOutput.write(Data("\n".utf8))
    }

    private static func fingerprint(_ text: String) -> String {
        SHA256.hash(data: Data(text.utf8)).map { String(format: "%02x", $0) }.joined()
    }
}

private enum ExampleError: Error {
    case metalUnavailable
}
