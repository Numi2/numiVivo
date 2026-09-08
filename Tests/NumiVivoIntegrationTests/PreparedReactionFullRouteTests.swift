import Foundation
import CryptoKit
import Darwin
import Testing
@testable import NumiVivoKit

/// Published synthetic input -> real preparation/Metal sampling -> exact export
/// -> fresh electronic/path qualification -> conditional first-event observable.
/// Partial sampling supplies geometry guesses only; the two models stay distinct.
@Suite(.serialized) struct PreparedReactionFullRouteTests: Sendable {
    private struct Campaign: Codable, Sendable {
        struct Sampling: Codable, Sendable {
            let replicaSeeds: [UInt64]
            let sourceTimesPS: [Double]
            let timeStepPS: Double
            let temperatureK: Double
            let frictionPerPS: Double
            let cutoffNM: Double
            let equilibrationSteps: UInt64
            let stepsPerBlock: UInt64
            let sampleEvery: UInt64
            let maximumBlocks: Int
            let requiredConsecutivePasses: Int
            let minimumRetainedFramesPerReplica: Int
            let minimumReplicas: Int
            let maximumRHat: Double
            let minimumEffectiveSamplesPerReplica: Double
            let maximumAutocorrelationLag: Int
            let maximumMeanStandardError: Double
            let trajectoryChunkBytes: Int
        }
        struct Encounter: Codable, Sendable {
            let reservoirPressurePa: Double
            let observationTimesSeconds: [Double]
        }
        let schema: String
        let sampling: Sampling
        let reactionTemplate: String
        let modelTransferStatement: String
        let encounter: Encounter
    }
    private struct Stage: Encodable {
        let task: VivoChemistryTask
        let taskFingerprint: VivoFingerprint
        let workflowReceiptFingerprint: VivoFingerprint
        let outputs: [VivoChemistryOutputReceipt]
        let reused: Bool
        init(_ task: VivoChemistryTask, _ result: VivoChemistryTaskResult) {
            self.task = task; taskFingerprint = result.taskFingerprint
            workflowReceiptFingerprint = result.receiptFingerprint; outputs = result.outputs; reused = result.reused
        }
    }
    private struct Implementation: Codable {
        let schema = "numivivo.org/test-evidence/prepared-reaction-implementation/v1"
        let sourceCommit: String
        let numericalContract: String
        let executableSHA256: String
        let operatingSystem: String
    }
    private struct Evidence: Encodable {
        let schema = "numivivo.org/test-evidence/prepared-reaction-campaign/v1"
        let status = "success"
        let implementation: Implementation
        let implementationFingerprint: VivoFingerprint
        let store: String
        let publishedPreparationBytes: VivoFingerprint
        let publishedCampaignBytes: VivoFingerprint
        let preparation: Stage
        let preparedStructureFingerprint: VivoFingerprint
        let compiledSystemFingerprint: VivoFingerprint
        let samplingRequestFingerprint: VivoFingerprint
        let samplingCheckpointFingerprint: VivoFingerprint
        let samplingReceiptFingerprint: VivoFingerprint
        let samplingTermination: String
        let declaredObservableCriteriaSatisfied: Bool
        let selectedReplica: Int
        let selectedSeed: UInt64
        let selectedStep: UInt64
        let selectedTimePS: Double
        let originalMDCheckpointFingerprint: VivoFingerprint
        let samplingExportReceiptFingerprint: VivoFingerprint
        let samplingExportVerification: VivoMolecularSamplingExportVerification
        let reactionSeedReceiptFingerprint: VivoFingerprint
        let reactionSeedVerification: VivoSamplingReactionSeedVerification
        let reaction: Stage
        let connectedPathTrials: Int
        let connectedPathComparisons: Int
        let rateConstant: Double
        let rateUnits: String
        let temperatureK: Double
        let molecularity: Int
        let assumedTransmissionCoefficient: Double
        let encounter: Stage
        let sourceReactionPayloadFingerprint: VivoFingerprint
        let conditionalHazardPerSecond: Double
        let observations: [VivoLinearKineticObservation]
        let negativeControls: [String: Bool]
        let interpretation = "Real preparation and accepted partial synthetic harmonic MD feed both H2 geometry guesses for fresh full-CI/STO-3G H3 exchange qualification. The isolated H components and saddle remain explicit independent destination inputs. Conditional probabilities describe a tagged H in an externally maintained ideal H2 reservoir under the declared TST model, standard state and assumed transmission. This is a complete finite software route, not a converged sampling ensemble, physical H2 force-field validation, experimental rate, diffusion simulation, finite-bath kinetics, efficacy or uncertainty qualification."
    }

    private var examples: URL {
        URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
            .deletingLastPathComponent().appendingPathComponent("Examples/prepared-reaction")
    }
    private func inputs() throws -> (Data, VivoMolecularPreparationRequest, Data, Campaign) {
        let preparation = try Data(contentsOf: examples.appendingPathComponent("preparation.json"))
        let campaign = try Data(contentsOf: examples.appendingPathComponent("campaign.json"))
        try #require(preparation.count < 64 * 1024 && campaign.count < 64 * 1024)
        return (preparation, try VivoCanonicalJSON.decode(VivoMolecularPreparationRequest.self, from: preparation),
                campaign, try VivoCanonicalJSON.decode(Campaign.self, from: campaign))
    }
    private func write<T: Encodable>(_ value: T, _ name: String, root: URL) throws {
        try VivoKineticsDocumentIO.write(VivoCanonicalJSON.encode(value), to: root.appendingPathComponent(name))
    }
    private func task(_ operation: VivoChemistryOperation, inputs: [VivoChemistryTaskInput],
                      budget: VivoChemistryBudget = .init()) -> VivoChemistryTask {
        .init(operation: operation.identifier, version: operation.version,
            implementationFingerprint: operation.implementationFingerprint, inputs: inputs,
            configuration: .object([:]), outputs: operation.outputs,
            resources: .init(budget: budget, maximumInputBytes: budget.maximumBytes, maximumOutputBytes: budget.maximumBytes))
    }
    private func output(_ name: String, _ result: VivoChemistryTaskResult) throws -> VivoChemistryOutputReceipt {
        try #require(result.outputs.first { $0.name == name })
    }
    private func rejects(_ body: () throws -> Void) throws -> Bool {
        var rejected = false
        do { try body() } catch is VivoChemistryError { rejected = true }
        try #require(rejected)
        return rejected
    }

    @Test func publishedPreparationRequiresItsExplicitBondParameter() throws {
        let (bytes, request, _, _) = try inputs()
        let valid = try VivoMolecularPreparation.prepare(request)
        let compiled = try #require(valid.compiledForceField)
        try #require(compiled.system.particles.count == 2 && compiled.system.bonds.count == 1)
        try #require(compiled.system.bonds[0].forceConstant == 1000 && compiled.system.bonds[0].lengthNM == 0.125)
        var missing = request
        missing.forceField = try #require(request.forceField)
        missing.forceField!.bondParameters = []
        var rejected = false
        do { _ = try VivoMolecularPreparation.prepare(missing) }
        catch VivoArtifactValidationError.unresolved(let message) {
            rejected = message.contains("missing bond parameter")
        }
        try #require(rejected, "a missing published parameter must not become an implicit interaction")
        #expect(try VivoCanonicalJSON.decode(VivoMolecularPreparationRequest.self, from: bytes) == request)
        #expect(request.forceField?.bondParameters.count == 1)
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_FULL_ROUTE_TESTS"] == "1"))
    func preparedSamplingSeedsFreshReactionAndConditionalEncounter() async throws {
        let environment = ProcessInfo.processInfo.environment
        let commit = try #require(environment["NUMIVIVO_TEST_SOURCE_COMMIT"])
        try #require(commit.utf8.count == 40 && commit.utf8.allSatisfy { (48...57).contains($0) || (97...102).contains($0) })
        let retained = environment["NUMIVIVO_TEST_ARTIFACTS"]
        let parent: URL
        if let retained {
            try #require(!retained.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty)
            parent = URL(fileURLWithPath: retained, isDirectory: true).standardizedFileURL
            try FileManager.default.createDirectory(at: parent, withIntermediateDirectories: true)
        } else { parent = FileManager.default.temporaryDirectory }
        let root = parent.appendingPathComponent("prepared-reaction-campaign-\(UUID().uuidString)")
        guard Darwin.mkdir(root.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        defer { if retained == nil { try? FileManager.default.removeItem(at: root) } }
        print("NUMIVIVO_PREPARED_REACTION_ROOT=\(root.path)")
        let storeRoot = root.appendingPathComponent("store")
        guard Darwin.mkdir(storeRoot.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        let store = try VivoArtifactStore(rootURL: storeRoot, createIfNeeded: false)
        let workflow = VivoChemistryWorkflow(store: store)
        let binary = try FileHandle(forReadingFrom: URL(fileURLWithPath: CommandLine.arguments[0]))
        defer { try? binary.close() }
        var digest = SHA256()
        while let bytes = try binary.read(upToCount: 4 * 1024 * 1024), !bytes.isEmpty { digest.update(data: bytes) }
        let implementation = Implementation(sourceCommit: commit, numericalContract: VivoMDExecutionIdentity.current,
            executableSHA256: digest.finalize().map { String(format: "%02x", $0) }.joined(),
            operatingSystem: ProcessInfo.processInfo.operatingSystemVersionString)
        let implementationID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(implementation))
        let (preparationBytes, preparationRequest, campaignBytes, campaign) = try inputs()
        try #require(campaign.schema == "numivivo.org/examples/prepared-reaction/v1")
        try #require(campaign.reactionTemplate == "h3-connected-rate")
        try VivoKineticsDocumentIO.write(preparationBytes, to: root.appendingPathComponent("preparation-input.json"))
        try VivoKineticsDocumentIO.write(campaignBytes, to: root.appendingPathComponent("campaign-input.json"))
        let preparationInput = try await store.put(data: preparationBytes, kind: "prepared-reaction-example-preparation", mediaType: "application/json")
        let campaignInput = try await store.put(data: campaignBytes, kind: "prepared-reaction-example-campaign", mediaType: "application/json")

        // The classical system is compiled by the existing preparation owner;
        // no direct synthetic VivoClassicalSystem construction occurs here.
        let calculation = VivoPreparedMolecularCalculation.prepare(request: preparationRequest)
        let preparedInput = try await store.put(data: VivoCanonicalJSON.encode(calculation),
            kind: "vivo.prepared-molecular-request", mediaType: "application/json")
        let preparationOperation = VivoPreparedMolecularWorkflow.operation(implementationFingerprint: implementationID)
        let preparationTask = task(preparationOperation,
            inputs: [.init(name: "request", artifact: preparedInput.fingerprint, kind: preparedInput.kind)])
        let preparationRun = try await workflow.run(preparationTask, using: preparationOperation)
        try #require(!preparationRun.reused)
        let preparationPort = try output("result", preparationRun)
        let preparedWrapper = try VivoCanonicalJSON.decode(VivoPreparedMolecularCalculationResult.self,
            from: await workflow.payload(artifact: preparationPort.artifact, expectedKind: preparationPort.kind))
        guard case .prepared(let prepared) = preparedWrapper else { throw VivoChemistryError.invalid("preparation omitted its prepared result") }
        let compiled = try #require(prepared.compiledForceField), system = compiled.system
        try #require(!prepared.requiresCoordinateRelaxation && prepared.sourceToPrepared == [0, 1])
        try #require(system.particles.count == 2 && system.particles.allSatisfy {
            $0.massDa == 1 && $0.chargeE == 0 && $0.sigmaNM == 0 && $0.epsilonKJPerMol == 0
        })
        try #require(system.bonds.count == 1 && system.bonds[0].lengthNM == 0.125 && system.bonds[0].forceConstant == 1000)
        try write(prepared, "preparation-result.json", root: root)
        let systemID = try system.fingerprint(), policy = campaign.sampling
        try #require(policy.replicaSeeds == [17, 29] && policy.sourceTimesPS == [0.125, 0.375])
        let initial = policy.sourceTimesPS.map {
            VivoClassicalInitialState(systemFingerprint: systemID,
                positionsNM: prepared.structure.conformers[0].positionsNM, sourceTimePS: $0)
        }
        let dynamics = VivoMDConfiguration(timeStepPS: policy.timeStepPS, cutoffNM: policy.cutoffNM, neighborSkinNM: 0,
            electrostatics: .cutoff, ensemble: .nvt, thermostat: .langevinMiddle,
            targetTemperatureK: policy.temperatureK, frictionPerPS: policy.frictionPerPS, neighborListEnabled: false)
        let samplingRequest = VivoMolecularSamplingRunRequest(structure: prepared.structure, system: system,
            initialStates: initial, replicaSeeds: policy.replicaSeeds, md: dynamics,
            contextIdentifier: "published prepared harmonic H2 partial-prefix geometry-seed campaign",
            observables: [.init(identifier: "distance", kind: .distance(atomA: 0, atomB: 1),
                maximumMeanStandardError: policy.maximumMeanStandardError)],
            convergence: .init(minimumRetainedFramesPerReplica: policy.minimumRetainedFramesPerReplica,
                minimumReplicas: policy.minimumReplicas, maximumRHat: policy.maximumRHat,
                minimumEffectiveSamplesPerReplica: policy.minimumEffectiveSamplesPerReplica,
                maximumAutocorrelationLag: policy.maximumAutocorrelationLag),
            equilibrationSteps: policy.equilibrationSteps, stepsPerBlock: policy.stepsPerBlock,
            sampleEvery: policy.sampleEvery, maximumBlocks: policy.maximumBlocks,
            requiredConsecutivePasses: policy.requiredConsecutivePasses, trajectoryChunkBytes: policy.trajectoryChunkBytes)
        try write(samplingRequest, "sampling-request.json", root: root)
        try Task.checkCancellation()
        let sampling = try await VivoMolecularSamplingRunner.run(samplingRequest, store: store)
        try write(sampling, "sampling-receipt.json", root: root)
        try #require(sampling.status == .budgetExhausted, "native sampling failed: \(sampling.diagnostic ?? "unknown")")
        try #require(sampling.completedBlocks == 4 && sampling.consecutivePasses == 0 && sampling.sampling?.converged == false)
        let samplingReceipt = try await store.put(data: VivoCanonicalJSON.encode(sampling), kind: "molecular-sampling-receipt", mediaType: "application/json")
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: store, checkpoint: sampling.checkpoint)
        let selected = try await reader.selectReplica(index: 0, validation: .allPayloads, receipt: samplingReceipt.fingerprint)
        try #require(selected.replicaSeed == 17 && selected.mdCheckpoint.acceptedStep == 36)
        try #require(selected.mdCheckpoint.timePS == 0.125 + 36.0 / 1024)
        try #require(selected.trajectoryValidation.verifiedPayloads == 4 && !selected.inspection.declaredObservableCriteriaSatisfied)
        let exported = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: implementationID)
        try write(exported.receipt, "sampling-export.json", root: root)
        let verified = try await VivoMolecularSamplingExporter.verifiedExport(exported.receipt, store: store,
            implementationFingerprint: implementationID, minimumValidation: .allPayloads)
        try #require(verified.selection.mdCheckpointData == selected.mdCheckpointData)
        try #require(verified.verification.sourceTermination == "budgetExhausted" && !verified.verification.declaredObservableCriteriaSatisfied)
        try write(verified.verification, "sampling-export-verification.json", root: root)

        let destination = try VivoReactionQualificationWorkflow.template(campaign.reactionTemplate)
        let seedRequest = VivoSamplingReactionSeedRequest(sourceExport: exported.receipt, destination: destination,
            assignments: [.init(target: .connectedEndpoint(identifier: "H0-H1_plus_H2", componentIndex: 0), sourceAtomByNucleus: [0, 1]),
                          .init(target: .connectedEndpoint(identifier: "H0_plus_H1-H2", componentIndex: 1), sourceAtomByNucleus: [0, 1])],
            modelTransferStatement: campaign.modelTransferStatement)
        try write(seedRequest, "reaction-seed-request.json", root: root)
        let seeded = try await VivoSamplingReactionSeed.publish(seedRequest, source: verified, implementationFingerprint: implementationID)
        try write(seeded.receipt, "reaction-seed-receipt.json", root: root)
        let seedVerification = try await VivoSamplingReactionSeed.verify(seeded.receipt, store: store,
            implementationFingerprint: implementationID, minimumValidation: .allPayloads)
        try write(seedVerification, "reaction-seed-verification.json", root: root)
        let reactionPort = try #require(seeded.receipt.outputs.first { $0.name == "request" && $0.kind == "vivo.reaction-calculation-request" })
        let reactionBytes = try await workflow.payload(artifact: reactionPort.artifact, expectedKind: reactionPort.kind)
        let reactionRequest = try VivoCanonicalJSON.decode(VivoReactionCalculationRequest.self, from: reactionBytes)
        try VivoKineticsDocumentIO.write(reactionBytes, to: root.appendingPathComponent("reaction-request.json"))
        guard case .connectedReaction(let connected) = reactionRequest.calculation,
              case .connectedReaction(let original) = destination.calculation else { throw VivoChemistryError.invalid("connected seed destination changed kind") }
        try #require(connected.saddle == original.saddle && connected.connectivity == original.connectivity && connected.transmission == original.transmission)
        try #require(connected.endpoints[0].components[1] == original.endpoints[0].components[1])
        try #require(connected.endpoints[1].components[0] == original.endpoints[1].components[0])
        for seed in [connected.endpoints[0].components[0], connected.endpoints[1].components[1]] {
            try #require(seed.qualification.massesDa == [1.008, 1.008] && seed.qualification.thermochemistry.temperatureK == 298.15)
            for index in 0..<2 {
                let p = selected.mdCheckpoint.positionsNM[index], q = seed.qualification.model.system.nuclei[index].positionBohr
                let scale = 1e-9 / VivoNuclearUnits.bohrM
                try #require(abs(q.x - p.x * scale) < 1e-12 && abs(q.y - p.y * scale) < 1e-12 && abs(q.z - p.z * scale) < 1e-12)
            }
        }
        let reactionOperation = VivoReactionQualificationWorkflow.operation(implementationFingerprint: implementationID)
        let reactionTask = task(reactionOperation,
            inputs: [.init(name: "request", artifact: reactionPort.artifact, kind: reactionPort.kind)], budget: reactionRequest.budget)
        try Task.checkCancellation()
        let reactionRun = try await workflow.run(reactionTask, using: reactionOperation)
        try #require(!reactionRun.reused, "this new store must freshly qualify the sampled geometry seeds")
        let qualifiedPort = try output("result", reactionRun)
        let qualifiedBytes = try await workflow.payload(artifact: qualifiedPort.artifact, expectedKind: qualifiedPort.kind)
        try VivoKineticsDocumentIO.write(qualifiedBytes, to: root.appendingPathComponent("reaction-result.json"))
        let result = try VivoCanonicalJSON.decode(VivoReactionCalculationResult.self, from: qualifiedBytes)
        guard case .connectedReaction(let rate) = result else { throw VivoChemistryError.invalid("fresh connected qualification omitted its TST result") }
        try #require(rate.request.connectivity.converged && rate.request.connectivity.trials.count == 4)
        try #require(rate.request.connectivity.comparisons.count == 12 && rate.request.connectivity.comparisons.allSatisfy(\.passed))
        try #require(rate.molecularity == 2 && rate.rateUnits == "Pa^-1 s^-1" && rate.temperatureK == 298.15)
        try #require(rate.rateConstant.isFinite && rate.rateConstant > 0 && rate.standardStateValue == 101325)
        try #require(!rate.hasIndependentTransmissionEvidence && rate.transmissionCoefficient == 1)

        let components = try VivoConditionalEncounter.componentBindings(in: rate)
        try #require(components.count == 2 && components[0].atomIdentifiers == ["H0", "H1"] && components[1].atomIdentifiers == ["H2"])
        let encounterRequest = VivoConditionalEncounterRequest(identifier: "tagged-H-in-maintained-H2-bath",
            reactantEndpointIdentifier: connected.reactantEndpointIdentifier, temperatureK: rate.temperatureK,
            taggedComponent: components[1], reservoirs: [.init(component: components[0], value: campaign.encounter.reservoirPressurePa, unit: .pascal)],
            observationTimesSeconds: campaign.encounter.observationTimesSeconds)
        try write(encounterRequest, "encounter-request.json", root: root)
        let encounterInput = try await store.put(data: VivoCanonicalJSON.encode(encounterRequest),
            kind: VivoConditionalEncounterWorkflow.requestKind, mediaType: "application/json")
        let encounterOperation = VivoConditionalEncounterWorkflow.operation(implementationFingerprint: implementationID)
        let encounterTask = task(encounterOperation, inputs: [
            .init(name: "request", artifact: encounterInput.fingerprint, kind: encounterInput.kind),
            .init(name: "reaction", artifact: qualifiedPort.artifact, kind: qualifiedPort.kind)], budget: reactionRequest.budget)
        let encounterRun = try await workflow.run(encounterTask, using: encounterOperation)
        try #require(!encounterRun.reused)
        let encounterPort = try output("result", encounterRun)
        let encounterBytes = try await workflow.payload(artifact: encounterPort.artifact, expectedKind: encounterPort.kind)
        try VivoKineticsDocumentIO.write(encounterBytes, to: root.appendingPathComponent("encounter-result.json"))
        let encounterOutput = try VivoCanonicalJSON.decode(VivoConditionalEncounterWorkflowResult.self, from: encounterBytes)
        let encounter = encounterOutput.encounter, hazard = rate.rateConstant * campaign.encounter.reservoirPressurePa
        let qualifiedPayloadFingerprint = try VivoCanonicalJSON.fingerprint(qualifiedBytes)
        let rateFingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(rate))
        try #require(encounterOutput.sourceReactionPayloadFingerprint == qualifiedPayloadFingerprint)
        try #require(encounter.sourceTransitionStateFingerprint == rateFingerprint)
        try #require(encounter.transmission == rate.request.transmission && encounter.rateUnits == rate.rateUnits && encounter.molecularity == 2)
        try #require(abs(encounter.conditionalHazardPerSecond - hazard) <= 1e-12 * max(1, hazard))
        try #require(encounter.kinetics.observations.map(\.timeSeconds) == campaign.encounter.observationTimesSeconds)
        for observation in encounter.kinetics.observations {
            let exponent = hazard * observation.timeSeconds, survival = exp(-exponent), reacted = -expm1(-exponent)
            try #require(abs(observation.survivalProbability - survival) <= 1e-10)
            try #require(abs(observation.reactedProbability - reacted) <= 1e-10)
            try #require(abs(observation.survivalProbability + observation.reactedProbability - 1) <= 1e-10)
            try #require(abs(observation.reactiveFluxPerSecond - hazard * survival) <= 1e-9 * max(1, hazard * survival))
            let derivative = try #require(observation.derivatives.first { $0.parameterIdentifier == "natural-log-rate" })
            try #require(abs(derivative.reactedProbability - exponent * survival) <= 1e-9)
        }
        var negatives: [String: Bool] = [:]
        let wrongTemperature = VivoConditionalEncounterRequest(identifier: encounterRequest.identifier,
            reactantEndpointIdentifier: encounterRequest.reactantEndpointIdentifier, temperatureK: rate.temperatureK + 1,
            taggedComponent: components[1], reservoirs: encounterRequest.reservoirs,
            observationTimesSeconds: encounterRequest.observationTimesSeconds)
        negatives["changed-temperature-rejected"] = try rejects { _ = try VivoConditionalEncounter.calculate(wrongTemperature, source: rate) }
        let wrongTag = VivoConditionalEncounterComponent(index: components[1].index, atomIdentifiers: components[1].atomIdentifiers,
            qualifiedPointFingerprint: preparationInput.fingerprint)
        let wrongComponent = VivoConditionalEncounterRequest(identifier: encounterRequest.identifier,
            reactantEndpointIdentifier: encounterRequest.reactantEndpointIdentifier, temperatureK: rate.temperatureK,
            taggedComponent: wrongTag, reservoirs: encounterRequest.reservoirs, observationTimesSeconds: encounterRequest.observationTimesSeconds)
        negatives["changed-component-identity-rejected"] = try rejects { _ = try VivoConditionalEncounter.calculate(wrongComponent, source: rate) }
        var rateJSON = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(rate)) as? [String: Any])
        rateJSON["rateConstant"] = rate.rateConstant * 2
        let forgedRate = try VivoCanonicalJSON.decode(VivoTransitionStateTheoryResult.self, from: JSONSerialization.data(withJSONObject: rateJSON))
        negatives["changed-numeric-rate-rejected"] = try rejects { _ = try VivoConditionalEncounter.calculate(encounterRequest, source: forgedRate) }

        let evidence = Evidence(implementation: implementation, implementationFingerprint: implementationID,
            store: storeRoot.path, publishedPreparationBytes: preparationInput.fingerprint, publishedCampaignBytes: campaignInput.fingerprint,
            preparation: Stage(preparationTask, preparationRun), preparedStructureFingerprint: try VivoStructureCodec.fingerprint(prepared.structure),
            compiledSystemFingerprint: systemID, samplingRequestFingerprint: sampling.requestFingerprint,
            samplingCheckpointFingerprint: sampling.checkpoint, samplingReceiptFingerprint: samplingReceipt.fingerprint,
            samplingTermination: sampling.status.rawValue, declaredObservableCriteriaSatisfied: selected.inspection.declaredObservableCriteriaSatisfied,
            selectedReplica: selected.replicaIndex, selectedSeed: selected.replicaSeed, selectedStep: selected.mdCheckpoint.acceptedStep,
            selectedTimePS: selected.mdCheckpoint.timePS, originalMDCheckpointFingerprint: selected.mdCheckpointFingerprint,
            samplingExportReceiptFingerprint: exported.artifact.fingerprint, samplingExportVerification: verified.verification,
            reactionSeedReceiptFingerprint: seeded.artifact.fingerprint, reactionSeedVerification: seedVerification,
            reaction: Stage(reactionTask, reactionRun), connectedPathTrials: rate.request.connectivity.trials.count,
            connectedPathComparisons: rate.request.connectivity.comparisons.count, rateConstant: rate.rateConstant,
            rateUnits: rate.rateUnits, temperatureK: rate.temperatureK, molecularity: rate.molecularity,
            assumedTransmissionCoefficient: rate.transmissionCoefficient, encounter: Stage(encounterTask, encounterRun),
            sourceReactionPayloadFingerprint: encounterOutput.sourceReactionPayloadFingerprint,
            conditionalHazardPerSecond: encounter.conditionalHazardPerSecond, observations: encounter.kinetics.observations,
            negativeControls: negatives)
        try Task.checkCancellation()
        let evidenceBytes = try VivoCanonicalJSON.encode(evidence)
        let receipt = try await store.put(data: evidenceBytes, kind: "prepared-reaction-campaign-receipt", mediaType: "application/json")
        try VivoKineticsDocumentIO.write(evidenceBytes, to: root.appendingPathComponent("prepared-reaction-campaign-receipt.json"))
        print("NUMIVIVO_PREPARED_REACTION_RECEIPT=\(root.appendingPathComponent("prepared-reaction-campaign-receipt.json").path) SHA256=\(receipt.fingerprint.hex)")
    }
}
